nextflow.preview.output = true

// BBDUK filters the reads for those containing k-mers of the target sequence
// Input: Paired-end reads (zstd compressed and interleaved), target sequence/k-mers of target seqeunce, k used for BBDuk
//    First input argument is a tuple (list of sample_divs, list of paths)
// Output: Filtered forward and reverse reads that match the reference
process BBDUK {
  label "medium"
  label "BBTools"
  input:
    tuple val(sample_divs), path(reads_files)
    path(ref_fasta_path)
    val(k)
  output:
    path("*_1.fastq"), emit: fwd_reads
    path("*_2.fastq"), emit: rev_reads
  script:
    """
    # Convert space-separated lists to arrays
    IFS=' ' read -ra samples <<< "${sample_divs.join(' ')}"
    IFS=' ' read -ra files <<< "${reads_files.join(' ')}"
    
    # Process each file with its corresponding sample_div
    for i in "\${!samples[@]}"; do
      sample_div="\${samples[\$i]}"
      reads_file="\${files[\$i]}"
      
      # Define input/output
      io="in=stdin.fastq ref=${ref_fasta_path} outm=\${sample_div}_1.fastq outm2=\${sample_div}_2.fastq"
      # Define parameters
      par="ordered=t rcomp=t minkmerhits=1 mm=f interleaved=t k=${k} t=${task.cpus} -Xmx${task.memory.toGiga()}g"
      # Execute
      zstdcat "\${reads_file}" | bbduk.sh \${io} \${par}
    done
    """
}

// CONCAT_READS concatenates several forward and reverse reads into a single forward and reverse read file
// Read headers are annotated with the basename of the file they came from
process CONCAT_READS{
  label "small"
  label "coreutils"
  input:
    path(fwd_reads)
    path(rev_reads)
  output:
    path("reads_1.fastq"), emit: final_fwd_read
    path("reads_2.fastq"), emit: final_rev_read
  script:
    """
    # Process forward reads
    for fwd in ${fwd_reads}; do
      filename=\$(basename "\$fwd" _1.fastq)
      awk -v fname="\$filename" 'NR%4==1 {print \$0 " " fname; next} {print}' "\$fwd"
    done > "reads_1.fastq"

    # Process reverse reads
    for rev in ${rev_reads}; do
      filename=\$(basename "\$rev" _2.fastq)
      awk -v fname="\$filename" 'NR%4==1 {print \$0 " " fname; next} {print}' "\$rev"
    done > "reads_2.fastq"
    """
}

workflow {
  main:
    // Load input files and create channel
    reads = Channel.fromPath(params.s3_files)
      .splitCsv(header: false, sep: "\t")
      .map { row -> tuple(row[0], row[1]) }
    
    // Batch inputs
    batch_size = params.batch_size ?: 10
    batched_reads = reads
      .collate(batch_size)
      .map { batch ->
        def sample_divs = batch.collect { it[0] }
        def files = batch.collect { file(it[1]) }
        tuple(sample_divs, files)
      }
    
    // Process batches
    bbduk_results = BBDUK(batched_reads, params.ref_fasta_path, params.kmer)
    
    // Collect and filter outputs, removing empty (no-hit) files
    fwd_reads = bbduk_results.fwd_reads
      .flatten()
      .filter { it.size() > 0 }
      .toSortedList { a, b -> 
        a.name.tokenize('_')[0] <=> b.name.tokenize('_')[0]
      }
    rev_reads = bbduk_results.rev_reads
      .flatten()
      .filter { it.size() > 0 }
      .toSortedList { a, b ->
        a.name.tokenize('_')[0] <=> b.name.tokenize('_')[0]
      }
    assert fwd_reads.size() == rev_reads.size() : "Forward and reverse reads must have equal length (fwd: ${fwd_reads.size()}, rev: ${rev_reads.size()})"
    CONCAT_READS(fwd_reads, rev_reads)

  publish:
    fwd_reads = CONCAT_READS.out.final_fwd_read
    rev_reads = CONCAT_READS.out.final_rev_read
}

output {
  fwd_reads {
    path "results"
  }
  rev_reads {
    path "results"
  }
}
