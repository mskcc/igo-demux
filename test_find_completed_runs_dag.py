import find_completed_runs_dag
import json
import os

def test_read_config():
    config = find_completed_runs_dag.read_config("dag_run.conf")
    assert 7 == len(config["sequencers"])

def test_find_completed_runs():
    find_completed_runs_dag.find_completed_runs(find_completed_runs_dag.read_config("dag_run.conf"), 30)
    assert 0 == 0

def test_create_sequencer_done_file(tmpdir):
    sequencer_json = '{"sequencers":[{"name": "sequencer_a","path":\"'+str(tmpdir.realpath())+'/sequencer_a\","last_file": "RTAComplete.txt"}]}'
    # create a RTAComplete.txt that should be picked up 
    p = tmpdir.mkdir("sequencer_a")
    run_dir_a = p.join("run1").mkdir()
    run_dir_b = p.join("run2").mkdir()
    run_dir_b_complete = run_dir_b.join("RTAComplete.txt")
    run_dir_b_complete.write("COMPLETE")
    sequencers_json = json.loads(sequencer_json)
    runs_to_demux = find_completed_runs_dag.find_completed_runs(sequencers_json,30)
    assert(os.path.isfile(tmpdir.join("sequencer_a").join("run2").join("RTAComplete.txt")))
    assert 1 == len(runs_to_demux)