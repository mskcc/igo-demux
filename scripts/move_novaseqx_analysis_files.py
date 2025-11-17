import os
import re
import subprocess
from pathlib import Path

SEQUENCERS = ["bono", "fauci2"]

SEQUENCER_BASE = "/igo/sequencers"
STATS_BASE = "/igo/staging/stats"
FASTQ_BASE = "/igo/staging/FASTQ"

def run_cmd(cmd):
    print(f"Running: {cmd}")
    subprocess.run(cmd, shell=True, check=True)

def discover_runs():
    """
    Find all run folders under each sequencer.
    A run folder looks like: yymmdd_SEQUENCERNAME_RUNID
    """
    runs = []
    pattern = re.compile(r"^\d{6}_(?P<seq>[A-Za-z0-9]+)_(?P<runid>.+)$")

    for seq in SEQUENCERS:
        seq_path = Path(SEQUENCER_BASE) / seq
        if not seq_path.exists():
            continue

        for entry in seq_path.iterdir():
            if entry.is_dir():
                m = pattern.match(entry.name)
                if m:
                    sequencer = m.group("seq")
                    run_id = m.group("runid")
                    runs.append((seq, sequencer, run_id, entry))
    return runs


def copy_analysis_files(source_run_path, dest_stats_dir):
    """
    Copy ALL files under:
        Analysis/1/Data/Project_*/DragenGermline/*_IGO_*/
    into DRAGEN folder (flat).
    """
    dragen_dir = Path(dest_stats_dir) / "DRAGEN"
    dragen_dir.mkdir(parents=True, exist_ok=True)

    search_root = Path(source_run_path) / "Analysis" / "1" / "Data"

    for project_dir in search_root.glob("Project_*"):
        dragen_germline_root = project_dir / "DragenGermline"

        for sample_dir in dragen_germline_root.glob("*_IGO_*"):
            # Copy files inside sample folder
            for root, dirs, files in os.walk(sample_dir):
                for f in files:
                    src_file = Path(root) / f
                    cmd = f'rsync -avh "{src_file}" "{dragen_dir}/"'
                    run_cmd(cmd)


def copy_fastqs(source_run_path, dest_fastq_dir):
    """
    Copy all FASTQ files under:
        Analysis/1/Data/Project_*/DragenGermline/fastq/*.fastq.gz
    """
    dest_fastq_dir.mkdir(parents=True, exist_ok=True)

    search_root = Path(source_run_path) / "Analysis" / "1" / "Data"

    for project_dir in search_root.glob("Project_*"):
        fastq_dir = project_dir / "DragenGermline" / "fastq"
        if fastq_dir.exists():
            for fastq in fastq_dir.glob("*.fastq.gz"):
                cmd = f'rsync -avh "{fastq}" "{dest_fastq_dir}/"'
                run_cmd(cmd)


def main():
    runs = discover_runs()
    if not runs:
        print("No runs found.")
        return

    for seq, seq_name, run_id, run_path in runs:
        seq_upper = seq_name.upper()
        dest_stats_dir = Path(STATS_BASE) / f"{seq_upper}_{run_id}"
        dest_fastq_dir = Path(FASTQ_BASE) / f"{seq_upper}_{run_id}"

        print(f"\n=== Processing {run_path} ===")
        print(f"Stats dir:   {dest_stats_dir}")
        print(f"FASTQ dir:   {dest_fastq_dir}")

        dest_stats_dir.mkdir(parents=True, exist_ok=True)

        # 1. Copy analysis → DRAGEN folder
        copy_analysis_files(run_path, dest_stats_dir)

        # 2. Copy fastqs → FASTQ folder
        copy_fastqs(run_path, dest_fastq_dir)

        print(f"✔ Completed copying for {seq_upper}_{run_id}\n")


if __name__ == "__main__":
    main()
