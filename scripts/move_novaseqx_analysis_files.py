import os
import re
import logging
import subprocess
from pathlib import Path

# -----------------------------
# CONFIG
# -----------------------------
SEQUENCERS = ["bono", "fauci2"]

SEQUENCER_BASE = Path("/igo/sequencers")
STATS_BASE      = Path("/igo/staging/stats")
FASTQ_BASE      = Path("/igo/staging/FASTQ")

COPY_COMPLETE_FILENAME = "CopyComplete.txt"

# Demux files to copy
DEMUX_FILES = [
    "Demultiplex_Stats.csv",
    "Top_Unknown_Barcodes.csv"
]

# -----------------------------
# LOGGING SETUP
# -----------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
)

logger = logging.getLogger("novaseqx_copy")


# -----------------------------
# UTILITY
# -----------------------------
def run_cmd(cmd):
    logger.info(f"Running: {cmd}")
    subprocess.run(cmd, shell=True, check=True)


# -----------------------------
# RUN DISCOVERY
# -----------------------------
def discover_runs():
    """
    Discover run folders named:
        yymmdd_SEQUENCERNAME_runId
    Returns list of tuples:
        (sequencer_folder_name, extracted_seq_name, run_id, run_path)
    """
    runs = []
    pattern = re.compile(r"^\d{6}_(?P<seq>[A-Za-z0-9]+)_(?P<runid>.+)$")

    for sequencer in SEQUENCERS:
        seq_root = SEQUENCER_BASE / sequencer
        if not seq_root.exists():
            continue

        for entry in seq_root.iterdir():
            if not entry.is_dir():
                continue

            m = pattern.match(entry.name)
            if m:
                extracted_seq = m.group("seq")
                run_id = m.group("runid")
                runs.append((sequencer, extracted_seq, run_id, entry))

    return runs


# -----------------------------
# CHECK RUN COMPLETION
# -----------------------------
def run_has_copy_complete(run_path: Path) -> bool:
    """
    A run is considered ready only if:
        {run}/Analysis/1/CopyComplete.txt exists
    """
    marker = run_path / "Analysis" / "1" / COPY_COMPLETE_FILENAME
    exists = marker.exists()

    if exists:
        logger.info(f"✔ Found CopyComplete.txt for run: {run_path}")
    else:
        logger.info(f"✖ No CopyComplete.txt — skipping run: {run_path}")

    return exists


# -----------------------------
# ANALYSIS FILE COPY
# -----------------------------
def copy_analysis_files(source_run_path: Path, dest_stats_dir: Path):
    """
    Copy all files under:
        Analysis/1/Data/Project_*/DragenGermline/*_IGO_*/
    into:
        stats/{SEQ}_{RUN}/DRAGEN/
    """

    search_root = source_run_path / "Analysis" / "1" / "Data"
    dragen_dir = dest_stats_dir / "DRAGEN"
    dragen_dir.mkdir(parents=True, exist_ok=True)

    logger.info(f"Copying ANALYSIS files → {dragen_dir}")

    for project_dir in search_root.glob("Project_*"):
        germline_root = project_dir / "DragenGermline"

        for sample_dir in germline_root.glob("*_IGO_*"):
            for root, dirs, files in os.walk(sample_dir):
                for f in files:
                    src = Path(root) / f
                    cmd = f'rsync -avh "{src}" "{dragen_dir}/"'
                    run_cmd(cmd)


# -----------------------------
# COPY DEMUX FILES
# -----------------------------
def copy_demux_files(source_run_path: Path, dest_fastq_root: Path):
    """
    Copies:
        Demultiplex_Stats.csv
        Top_Unknown_Barcodes.csv
    from:
        Analysis/1/Data/Demux/
    into:
        FASTQ/<SEQ>_<RUNID>/
    """
    demux_dir = source_run_path / "Analysis" / "1" / "Data" / "Demux"
    if not demux_dir.exists():
        logger.info(f"No Demux directory found at {demux_dir}")
        return

    logger.info(f"Copying Demux files → {dest_fastq_root}")

    for filename in DEMUX_FILES:
        src = demux_dir / filename
        if src.exists():
            cmd = f'rsync -avh "{src}" "{dest_fastq_root}/"'
            run_cmd(cmd)
        else:
            logger.info(f"Demux file not found: {src}")


# -----------------------------
# FASTQ FILE COPY
# -----------------------------
def copy_fastqs(source_run_path: Path, dest_fastq_root: Path):
    """
    Organize FASTQ destination as:
        FASTQ/{SEQ}_{RUN}/Project_<name>/Sample_<sampleFolderName>/
    """

    search_root = source_run_path / "Analysis" / "1" / "Data"

    logger.info(f"Copying FASTQ files → {dest_fastq_root}")

    for project_dir in search_root.glob("Project_*"):
        project_name = project_dir.name
        fastq_root = project_dir / "DragenGermline" / "fastq"

        if not fastq_root.exists():
            continue

        # /FASTQ/SEQ_RUN/Project_12345_A/
        project_dest = dest_fastq_root / project_name
        project_dest.mkdir(parents=True, exist_ok=True)

        # fastq/*_IGO_*
        for sample_dir in fastq_root.glob("*_IGO_*"):
            dest_sample = project_dest / sample_dir.name
            dest_sample.mkdir(parents=True, exist_ok=True)

            for fastq in sample_dir.glob("*.fastq.gz"):
                cmd = f'rsync -avh "{fastq}" "{dest_sample}/"'
                run_cmd(cmd)


# -----------------------------
# PROCESS EACH RUN
# -----------------------------
def process_run(sequencer_folder, seq_name, run_id, run_path: Path):
    seq_upper = seq_name.upper()

    # Check CopyComplete.txt
    if not run_has_copy_complete(run_path):
        return  # Skip entirely

    logger.info(f"=== Processing run {seq_upper}_{run_id} ({run_path}) ===")

    dest_stats_dir = STATS_BASE / f"{seq_upper}_{run_id}"
    dest_fastq_dir = FASTQ_BASE / f"{seq_upper}_{run_id}"

    dest_stats_dir.mkdir(parents=True, exist_ok=True)
    dest_fastq_dir.mkdir(parents=True, exist_ok=True)

    logger.info(f"→ Stats destination:  {dest_stats_dir}")
    logger.info(f"→ FASTQ destination:  {dest_fastq_dir}")

    # Copy analysis (DRAGEN)
    copy_analysis_files(run_path, dest_stats_dir)

    # Copy FASTQs
    copy_fastqs(run_path, dest_fastq_dir)

    # Copy Demux stats + unknown barcodes
    copy_demux_files(run_path, dest_fastq_dir)

    logger.info(f"✔ Finished copying for {seq_upper}_{run_id}\n")


# -----------------------------
# MAIN
# -----------------------------
def main():
    logger.info("=== Starting NovaSeqX analysis + FASTQ copy script ===")

    runs = discover_runs()
    if not runs:
        logger.info("No run directories found.")
        return

    for seq, seq_name, run_id, run_path in runs:
        process_run(seq, seq_name, run_id, run_path)

    logger.info("=== Completed run processing ===")


if __name__ == "__main__":
    main()
