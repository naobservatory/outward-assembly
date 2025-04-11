import pytest

from outward_assembly.io_helpers import process_s3_paths, _count_lines


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


def test_process_s3_paths_invalid():
    with pytest.raises(ValueError):
        process_s3_paths(["s3://not/valid/extension.jpg"])


def test_count_lines(temp_empty_file):
    assert _count_lines(temp_empty_file) == 0
    with open(temp_empty_file, "a") as f:
        f.write("hi")  # no newline
    assert _count_lines(temp_empty_file) == 1
    with open(temp_empty_file, "a") as f:
        f.write("\n")  # finish the same line
    assert _count_lines(temp_empty_file) == 1
