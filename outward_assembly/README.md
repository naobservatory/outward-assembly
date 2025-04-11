# Developer notes

This is a collection of notes for developers working on the outward assembly pipeline.

## Nomenclature

The following nomenclature is used in the outward assembly pipeline:
* outer iteration: An outer iteration is a single run of the outward assembly pipeline and will only be referenced when we are running outward assembly in *automated* mode.
* inner iteration: An inner iteration refers to the amount of iterations it takes for one run of the outward assembly pipeline to complete.

Hence, *normal* mode will only have one outer iteration and a variable number of inner iterations, while *automated* mode will potentially have multiple outer iterations and different amounts of inner iterations for each outer iteration.

## Code notes
* The naive algorithm used to construct the contig overlap graph requires checking all pairs of contigs for overlap and thus has runtime `O(n^2)` in `n` the number of contigs. Typically this step is super fast because we have few contigs, but if you're finding outward assembly stalling with a single Python process utilizing ~100% of a CPU core, this is a plausible candidate.
    * One way to accidentally get lots of contigs: use a too-small `read_subset_k` argument when calling `outward_assembly`, which can result in BBDuk selecting lots of reads unrelated to the seed sequence.
* For code formatting, we use [Black](https://github.com/psf/black) and [isort](https://github.com/PyCQA/isort).

## Testing

We use the [pytest framework](https://docs.pytest.org/en/stable/) to write tests.

You can run tests with `pytest` in the base repo dir, or `pytest -s` to see the real time logs from the `outward_assembly` end to end test.

We currently have all tests enabled to run in both local and batch mode, hence you'll have to have gone through the [installation process](../docs/installation.md#optional-batch-profile) for batch mode.

To run the tests, copy the below`.env.example` file to `.env` and update the variables.

```.env.example
BATCH_QUEUE='my-batch-queue' # Update this to your aws batch queue
BATCH_WORKDIR='s3://my-bucket/my-outward-workdir' # Update this to your s3 directory for batch
```

Then, run the following command in the base repo dir:

```
pytest -s
```

to run all tests.

## Improvements
* [Google doc](https://docs.google.com/document/d/1AiQUWMNUhbwYZBqLleZ1K-XXnND84z0tURidb2OD8sw/edit?tab=t.0) with some ideas for performance and sensitivity enhancements.
* Richer test coverage; existing test coverage is pretty modest.
* Set up GitHub actions testing
* Rearranging code for greater clarity and shorter files
* Better logging:
    * More comprehensive
    * Principled usage of log levels
    * Optionally not discarding all logs from parallel processes (currently written to `/dev/null`)
* Refactor use of working directory so that every iteration of outward assembly has its own subdirectory and we're not updating files like `current_contigs.fasta` in a stateful way.
