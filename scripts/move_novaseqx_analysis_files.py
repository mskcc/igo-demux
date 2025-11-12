import os
import sys
import glob
import subprocess

def move_novaseqx_analysis_files(run_id: str):
    """
    Move all germline_seq files from:
      /igo/sequencers/<sequencer>/<RunID>/Analysis/1/Data/Project_*/DragenGermline/*_IGO_*/germline_seq/
    into a single flat destination folder:
      /igo/staging/stats/<sequencer>_<RunID>/DRAGEN/
    """

    # Detect sequencer automatically from RunID
    run_id_upper = run_id.upper()
    if "BONO" in run_id_upper:
        sequencer = "bono"
    elif "FAUCI2" in run_id_upper:
        sequencer = "fauci2"
    else:
        print(f"❌ ERROR: Unknown sequencer in RunID: {run_id}")
        sys.exit(1)

    base_src = f"/igo/sequencers/{sequencer}/{run_id}/Analysis/1/Data"
    base_dest = f"/igo/staging/stats/{sequencer}_{run_id}/DRAGEN"

    if not os.path.exists(base_src):
        print(f"❌ Source directory not found: {base_src}")
        sys.exit(1)

    os.makedirs(base_dest, exist_ok=True)
    print(f"✅ Destination directory ready: {base_dest}")

    # Find all germline_seq folders under Project_*/DragenGermline/*_IGO_*/
    search_pattern = os.path.join(base_src, "Project_*", "DragenGermline", "*_IGO_*", "germline_seq")
    germline_dirs = sorted(glob.glob(search_pattern))

    if not germline_dirs:
        print(f"⚠️ No germline_seq directories found under {base_src}")
        sys.exit(0)

    print(f"🔍 Found {len(germline_dirs)} germline_seq directories to copy.")

    for src_dir in germline_dirs:
        project_name = src_dir.split("/")[-4]   # Project_*
        sample_name = src_dir.split("/")[-2]    # *_IGO_*
        print(f"  Copying {project_name}/{sample_name}")

        rsync_cmd = [
            "rsync",
            "-avh",
            "--progress",
            f"{src_dir}/",
            f"{base_dest}/"
        ]

        try:
            subprocess.run(rsync_cmd, check=True)
            print(f"  ✅ Done copying {project_name}/{sample_name}")
        except subprocess.CalledProcessError as e:
            print(f"  ❌ rsync failed for {project_name}/{sample_name}: {e}")

    print("\n All files successfully copied into DRAGEN folder.")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python3 scripts/move_novaseqx_analysis_files.py <RunID>")
        sys.exit(1)

    run_id = sys.argv[1].strip()
    move_novaseqx_analysis_files(run_id)
