import pytest
from Bio.Seq import Seq

from outward_assembly.basic_seq_operations import is_subseq


@pytest.mark.fast
@pytest.mark.unit
def test_is_subseq_empty():
    """Test behavior with empty strings"""
    assert is_subseq("", "ACGT")
    assert not is_subseq("ACGT", "")
    assert is_subseq("", "")


@pytest.mark.fast
@pytest.mark.unit
def test_is_subseq_exact_match():
    """Test when haystack is exactly the needle"""
    needle = Seq("AAGCT")
    haystack = Seq("AAGCT")
    assert is_subseq(needle, haystack)
    assert is_subseq(needle, haystack, check_rc=False)


@pytest.mark.fast
@pytest.mark.unit
def test_is_subseq_contained():
    """Test when needle is contained within larger haystack"""
    needle = Seq("AAGCT")
    haystack = Seq("TTAAGCTCC")
    assert is_subseq(needle, haystack)
    assert is_subseq(needle, haystack, check_rc=False)


@pytest.mark.fast
@pytest.mark.unit
def test_is_subseq_rc_only():
    """Test when only reverse complement is present"""
    needle = Seq("AAGT")  # RC ACTT
    haystack = Seq("CCACTTGG")
    assert is_subseq(needle, haystack, check_rc=True)
    assert not is_subseq(needle, haystack, check_rc=False)


@pytest.mark.fast
@pytest.mark.unit
def test_is_subseq_str_inputs():
    """Verify function works with string inputs as well as Seq objects"""
    assert is_subseq("AAGT", "CCACTTGG", check_rc=True)
    assert not is_subseq("AAGT", "CCACTTGG", check_rc=False)
