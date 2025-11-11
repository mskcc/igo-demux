import os
import sys
import glob
import subprocess

def move_novaseqx_analysis_files(run_id: str):
    """
    Move sequencing analysis results from:
      /igo/sequencers/<sequencer>/<RunID>/Analysis/1/Data/Project_*/DragenGermline/*_IGO_*/germline_seq/
    to:
      /igo/staging/stats/<sequencer>_<RunID>/Project_*/<_IGO_*>

    The project and sample directory names are preserved exactly as they appear in the source.
    """

    # Normalize RunID and detect sequencer
    run_id_upper = run_id.upper()
    if "BONO" in run_id_upper:
        sequencer = "bono"
    elif "FAUCI2" in run_id_upper:
        sequencer = "fauci2"
    else:
        print(f"❌ ERROR: Unknown sequencer in RunID: {run_id}")
        sys.exit(1)

    # Define base source and destination directories
    base_src = f"/igo/sequencers/{sequencer}/{run_id}/Analysis/1/Data"
    base_dest = f"/igo/staging/stats/{sequencer}_{run_id}"

    # Validate source existence
    if not os.path.exists(base_src):
        print(f"❌ Source directory not found: {base_src}")
        sys.exit(1)

    os.makedirs(base_dest, exist_ok=True)
    print(f"✅ Created or verified destination root: {base_dest}")

    # Find all project-level directories
    project_dirs = sorted(glob.glob(os.path.join(base_src, "Project_*")))
    if not project_dirs:
        print(f"⚠️ No Project_* directories found under {base_src}")
        sys.exit(0)

    for project_dir in project_dirs:
        project_name = os.path.basename(project_dir)
        project_dest = os.path.join(base_dest, project_name)
        os.makedirs(project_dest, exist_ok=True)

        print(f"\n Processing project: {project_name}")

        # Find all sample directories (matching *_IGO_*)
        sample_dirs = sorted(glob.glob(os.path.join(project_dir, "DragenGermline", "*_IGO_*")))
        if not sample_dirs:
            print(f"  ⚠️ No *_IGO_* directories found under {project_dir}/DragenGermline/")
            continue

        for sample_dir in sample_dirs:
            sample_name = os.path.basename(sample_dir)
            sample_dest = os.path.join(project_dest, sample_name)
            os.makedirs(sample_dest, exist_ok=True)

            src_germline_seq = os.path.join(sample_dir, "germline_seq")
            if not os.path.exists(src_germline_seq):
                print(f"  ⚠️ germline_seq directory missing in {sample_dir}")
                continue

            # rsync copy while preserving all files inside germline_seq
            rsync_cmd = [
                "rsync",
                "-avh",
                "--progress",
                f"{src_germline_seq}/",
                f"{sample_dest}/"
            ]

            print(f" Copying sample {sample_name} ...")
            try:
                subprocess.run(rsync_cmd, check=True)
                print(f"  ✅ Copied {sample_name} → {sample_dest}")
            except subprocess.CalledProcessError as e:
                print(f"  ❌ rsync failed for {sample_name}: {e}")

    print("\n All eligible project/sample directories processed successfully.")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python3 scripts/move_novaseqx_analysis_files.py <RunID>")
        sys.exit(1)

    run_id = sys.argv[1].strip()
    move_novaseqx_analysis_files(run_id)
