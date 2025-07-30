import pytest

from outward_assembly.io_helpers import (
    process_s3_paths,
    _count_lines,
    concat_and_tag_fastq,
)


@pytest.mark.fast
@pytest.mark.unit
def test_process_s3_paths():
    """Processing s3 paths should produce unique valid filenames."""
    s3_paths = [
        # same suffix, we should disambiguate
        "s3://a/b/c.fastq.zst",
        "s3://a/x/c.fastq.zst",
    ]
    records = process_s3_paths(s3_paths)
    assert len(set([rec.filename for rec in records])) == len(s3_paths)
    assert sorted(s3_paths) == sorted([rec.s3_path for rec in records])
    for rec in records:
        # filename should be a valid file name; could check more conditions
        assert "/" not in rec.filename
        assert len(rec.filename) <= 245  # need extra space for other prefixes/suffixes


@pytest.mark.fast
@pytest.mark.unit
def test_process_s3_paths_invalid():
    with pytest.raises(ValueError):
        process_s3_paths(["s3://not/valid/extension.jpg"])


@pytest.mark.fast
@pytest.mark.unit
def test_count_lines(temp_empty_file):
    assert _count_lines(temp_empty_file) == 0
    with open(temp_empty_file, "a") as f:
        f.write("hi")  # no newline
    assert _count_lines(temp_empty_file) == 1
    with open(temp_empty_file, "a") as f:
        f.write("\n")  # finish the same line
    assert _count_lines(temp_empty_file) == 1


@pytest.mark.fast
@pytest.mark.unit
def test_concat_and_tag_fastq_forward_reads(temp_workdir):
    sample1_fwd = temp_workdir / "sample1_1.fastq"
    sample2_fwd = temp_workdir / "sample2_1.fastq"
    output_fwd = temp_workdir / "combined_1.fastq"

    # Create test forward read files
    with open(sample1_fwd, "w") as f:
        f.write("@read1\nACGT\n+\nIIII\n@read2\nTGCA\n+\nIIII\n")

    with open(sample2_fwd, "w") as f:
        f.write("@read3\nGGGG\n+\nIIII\n")

    concat_and_tag_fastq([sample1_fwd, sample2_fwd], output_fwd)

    with open(output_fwd, "r") as f:
        content = f.read()

    expected = "@read1 sample1\nACGT\n+\nIIII\n@read2 sample1\nTGCA\n+\nIIII\n@read3 sample2\nGGGG\n+\nIIII\n"
    assert content == expected


@pytest.mark.fast
@pytest.mark.unit
def test_concat_and_tag_fastq_empty_file(temp_workdir):
    empty_fwd = temp_workdir / "empty_1.fastq"
    sample_fwd = temp_workdir / "sample_1.fastq"
    output_fwd = temp_workdir / "combined_1.fastq"

    # Create empty and non-empty files
    with open(empty_fwd, "w") as f:
        pass  # Empty file

    with open(sample_fwd, "w") as f:
        f.write("@read1\nACGT\n+\nIIII\n")

    concat_and_tag_fastq([empty_fwd, sample_fwd], output_fwd)

    with open(output_fwd, "r") as f:
        content = f.read()

    expected = "@read1 sample\nACGT\n+\nIIII\n"
    assert content == expected


@pytest.mark.fast
@pytest.mark.unit
def test_concat_and_tag_fastq_file_not_found(temp_workdir):
    nonexistent_fwd = temp_workdir / "nonexistent_1.fastq"
    output_fwd = temp_workdir / "combined_1.fastq"

    with pytest.raises(RuntimeError, match="Error processing files"):
        concat_and_tag_fastq([nonexistent_fwd], output_fwd)
