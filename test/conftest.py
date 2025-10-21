import pytest
import yaml
from pathlib import Path

@pytest.fixture(scope="session")
def test_config():
    """Load test config for Nextflow pipeline tests."""
    config_path = Path(__file__).parent / "test_config.yml"
    with open(config_path) as f:
        return yaml.safe_load(f)

@pytest.fixture
def workdir(tmp_path): 
    """Temporary work directory for Nextflow runs. Will be cleaned up after tests."""
    return tmp_path

@pytest.fixture
def FlowCelldir(tmp_path): 
    """Temporary work directory for Nextflow runs. Will be cleaned up after tests."""
    tmp_fcpath = tmp_path / "testFC"
    tmp_fcpath.mkdir()
    return tmp_fcpath