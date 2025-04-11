from pathlib import Path
import pytest
import tempfile


@pytest.fixture
def temp_workdir():
    """Create and clean up a temporary working directory."""
    with tempfile.TemporaryDirectory() as tmpdir:
        yield Path(tmpdir)


@pytest.fixture
def temp_empty_file(temp_workdir):
    """Create and return path to a temporary empty file."""
    path = temp_workdir / "testfile.txt"
    path.touch()
    yield path
    path.unlink(missing_ok=True)
