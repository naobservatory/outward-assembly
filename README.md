# Outward assembly
Welcome to the outward assembly repository! This repo provides Python tooling to assemble the genomic context around a provided _seed sequence_ from a collection of read pairs. The seed sequence could be a flagged chimeric junction, a kmer identified by reference-free growth detection, etc.

When the collection of reads is very large -- e.g. billions of reads from one delivery, or tens of billions of reads across deliveries -- doing a full metagenomic joint assembly is slow and computationally expensive. Instead of jointly assembling all reads, outward assembly attempts to iteratively grow contigs outward from the provided seed. The basic algorithm is simple:
1. Start with the seed as your initial contig.
2. Iteratively:
    * Find all read pairs that share a kmer with your contigs.
    * Assemble these read pairs.
    * Filter contigs to those that contain the seed.
3. Continue until either maximum iterations are reached or the assembly algorithm does not make progress from one iteration to another ("convergence").

Although the basic algorithm is simple, in practice, getting good assembly results is a bit more complex (see [algorithm details](docs/algorithm_details.md) for how we handle some of these complexities) and might require running outward assembly multiple times with different parameters. To handle this, we've introduced a way to automate the labor-intensive iterative process of:

1. Running outward assembly;
2. Evaluating outputs;
3. Deciding whether to accept the outputs as final or modify parameters and re-run.

Whether you're running outward assembly once or iteratively, you can find more information in the [usage docs](docs/usage.md).

## Quick Start

1. Install dependencies: `uv sync --extra dev`
2. Create tools environment: `mamba env create -n oa-tools -f oa_tools_env.yml --channel-priority flexible`
3. Activate tools environment: `mamba activate oa-tools`
4. Run with: `uv run your_script.py`

See [installation docs](docs/installation.md) for detailed setup instructions.

## Documentation

- [Installation](docs/installation.md)
- [Usage](docs/usage.md)
- [Algorithm details](docs/algorithm_details.md)  
- [Changelog](CHANGELOG.md)
