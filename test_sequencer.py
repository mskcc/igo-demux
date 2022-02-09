import sequencer
import os
import tempfile


def test_create_sequencer_done_file():

    tmpdir = tempfile.mkdtemp()
    # create a CopyComplete.txt that should be picked up in
    sequencer_path = os.path.join(tmpdir, "sequencer_a")
    os.mkdir(sequencer_path)
    run_dir_a = os.mkdir(os.path.join(sequencer_path, "run1"))
    run_dir_b = os.mkdir(os.path.join(sequencer_path, "run2"))
    run_dir_b_complete = os.path.join(sequencer_path, "run2/CopyComplete.txt")
    open(run_dir_b_complete, "x")
    print("Created {}".format(run_dir_b_complete))
    sequencers = {
        "sequencers": [
            {
                "name": "ayyan",
                "path": sequencer_path,
                "last_file": "CopyComplete.txt"
            }
        ]
    }
    print(sequencers)
    ready_to_demux_list = sequencer.find_completed_runs(sequencers, 10)

    assert(len(ready_to_demux_list) == 1)
