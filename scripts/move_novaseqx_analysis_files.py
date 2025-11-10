import subprocess
import os
import sys
import glob

def detect_sequencer(run_id):
    """Infer sequencer name (bono or fauci2) from the run ID."""
    run_id_lower = run_id.lower()
    if "bono" in run_id_lower:
        return "bono"
    elif "fauci2" in run_id_lower:
        return "fauci2"
    else:
        print("Could not determine sequencer type from Run ID.")
        print("Run ID must contain either 'BONO' or 'FAUCI2'.")
        sys.exit(1)

def copy_run_stats(run_id):
    sequencer = detect_sequencer(run_id)
    print(f"Detected sequencer: {sequencer.upper()}")

    # Create destination folder
    dest_dir = f"/igo/staging/stats/{run_id}"
    os.makedirs(dest_dir, exist_ok=True)

    # Build source path pattern dynamically
    source_pattern = f"/igo/sequencers/{sequencer}/*_{run_id}/Analysis/1/Data/Project_*/DragenGermline/*_IGO_*/germline_seq/"
    source_dirs = glob.glob(source_pattern)

    if not source_dirs:
        print(f"No matching source directories found for pattern:\n  {source_pattern}")
        sys.exit(1)

    for src in source_dirs:
        print(f"Copying from: {src}")

        rsync_cmd = [
            "rsync",
            "-avh",          # archive, verbose, human-readable
            "--progress",    # show progress
            f"{src}",        # include trailing slash to copy contents
            dest_dir
        ]

        try:
            subprocess.run(rsync_cmd, check=True)
            print(f"✅ Successfully copied from {src} to {dest_dir}\n")
        except subprocess.CalledProcessError as e:
            print(f"❌ rsync failed for {src}: {e}")
            sys.exit(1)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python3 move_novaseqx_analysis_files.py <RunID>")
        sys.exit(1)

    run_id = sys.argv[1].strip()
    copy_run_stats(run_id)
