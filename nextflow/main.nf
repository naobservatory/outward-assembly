nextflow.preview.output = true

// BBDUK filters the reads for those containing the target sequence (or k-mers of the target sequence)
// Input: Paired-end reads (zstd compressed and interleaved), target sequence/k-mers of target seqeunce, k used for BBDuk
// Output: Filtered forward and reverse reads that match the reference
// Parameters:
//   - Uses exact kmer matching (mm=f)
//   - Maintains read order (ordered=t)
//   - Processes both forward and reverse complement (rcomp=t)
//   - Expects interleaved input (interleaved=t)
//   - Requires at least 1 kmer hit (minkmerhits=1)
process BBDUK {
    label "medium"
    label "BBTools"
    input:
        tuple val(sample_div), path(reads)
        path(ref_fasta_path)
        val(k)
    output:
        tuple val(sample_div), path("${sample_div}_{1,2}.fastq")
    shell:
        '''
        # Define input/output
        io="in=stdin.fastq ref=!{ref_fasta_path} outm=!{sample_div}_1.fastq outm2=!{sample_div}_2.fastq "
        # Define parameters
        par="ordered=t rcomp=t minkmerhits=1 mm=f interleaved=t k=!{k} t=!{task.cpus} -Xmx!{task.memory.toGiga()}g"
        # Execute
        zstdcat !{reads} | bbduk.sh ${io} ${par}
        '''
}

// CONCAT_READS concatenates several forward and reverse reads into a single forward and reverse read file, adding the filename to the read id
// Input: Several forward and reverse read files (from BBDUK process)
// Output: Two consolidated FASTQ files (one forward, one reverse)
process CONCAT_READS{
  label "small"
  label "coreutils"
  input:
    path(fwd_reads)
    path(rev_reads)
  output:
    path("reads_1.fastq"), emit: final_fwd_read
    path("reads_2.fastq"), emit: final_rev_read
  shell:
    '''
    # Process forward reads
    for fwd in !{fwd_reads}; do
      filename=$(basename "$fwd" _1.fastq)
      awk -v fname="$filename" 'NR%4==1 {print $0 " " fname; next} {print}' "$fwd"
    done > "reads_1.fastq"

    # Process reverse reads
    for rev in !{rev_reads}; do
      filename=$(basename "$rev" _2.fastq)
      awk -v fname="$filename" 'NR%4==1 {print $0 " " fname; next} {print}' "$rev"
    done > "reads_2.fastq"
    '''
}

workflow {
  main:
    reads = Channel.fromPath(params.s3_files).splitCsv(header: false, sep: "\t"). map {
      row -> tuple(row[0], row[1])
    }
    bbduk_reads = BBDUK(reads, params.ref_fasta_path, params.kmer)

    // Filter out empty reads
    filtered_bbduk_reads = bbduk_reads
      .filter { _sample_div, fastq_files -> 
          fastq_files[0].size() > 0 && fastq_files[1].size() > 0
      }

    fwd_reads = filtered_bbduk_reads
        .toSortedList { a, b -> a[0] <=> b[0] }  // Sort by sample_div
        .flatMap()  // Flatten the sorted list back to a channel
        .map { _sample_div, fastq_files -> 
            fastq_files[0]  // This selects the first FASTQ file (ending with _1.fastq)
        }
        .collect()

    rev_reads = filtered_bbduk_reads
        .toSortedList { a, b -> a[0] <=> b[0] }  // Sort by sample_div
        .flatMap()  // Flatten the sorted list back to a channel
        .map { _sample_div, fastq_files -> 
            fastq_files[1]  // This selects the second FASTQ file (ending with _2.fastq)
        }
        .collect()

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
