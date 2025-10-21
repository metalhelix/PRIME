import subprocess
from pathlib import Path

def test_pipeline_run_bulk(test_config, workdir):
    # Run Nextflow with test config
    cmd = [
        "nextflow", "run", "/n/ngs/tools/PRIME/RUN/main.nf",
        "--fcid", test_config["params"]["fcid_bulk"],
        "--fcpath", FlowCelldir,
        "--fastq_dir", test_config["input"]["bulk_fastq_dir"],
        "--work-dir", workdir 
    ]
    
    print(f"Testing bulk-RNAseq, command: {' '.join(cmd)}", )  # For debugging
    
    result = subprocess.run(cmd, cwd=Path(__file__).parents[1], capture_output=True)

    # Check run was successful
    assert result.returncode == 0, f"Nextflow bulk-RNAseq test failed:\n{result.stderr.decode()}"

def test_pipeline_run_scRNA(test_config, workdir):
    cmd = [
        "nextflow", "run", "/n/ngs/tools/PRIME/RUN/main.nf",
        "--fcid", test_config["params"]["fcid_single_cell"],
        "--fcpath", FlowCelldir,
        "--fastq_dir", test_config["input"]["single_cell_fastq_dir"],
        "--work-dir", workdir,
        "--run_type", "10X"   
    ]

    print(f"Testing single-cell RNAseq, command: {' '.join(cmd)}")  # For debugging
    
    result = subprocess.run(cmd, cwd=Path(__file__).parents[1], capture_output=True)
    
    assert result.returncode == 0, f"Nextflow single-cell RNAseq failed:\n{result.stderr.decode()}"