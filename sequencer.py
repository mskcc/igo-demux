import json
import os
import datetime
import pathlib

def read_config(config_file):
    print("Reading sequencer config file:", config_file)
    with open(config_file, "r") as json_file:
      data = json.load(json_file)
      print("Sequencers: ")
      for sequencer in data["sequencers"]:
        print(" ", sequencer["name"])
    return data

def find_completed_runs(sequencers, mins_ago):
    time_minutes_ago = datetime.datetime.now() - datetime.timedelta(minutes=int(mins_ago))
    print("{} - Searching for runs completed since {}".format(datetime.datetime.now(), time_minutes_ago))
    
    ready_to_demux = list()
    for sequencer in sequencers["sequencers"]:
        print(" ", sequencer["name"])
        sequencer_path = sequencer["path"]
        last_file = sequencer["last_file"]
        if os.path.isdir(sequencer_path):
            runs = os.listdir(sequencer_path)
            for run_folder in runs:
                full_run_path = os.path.join(sequencer_path,run_folder,last_file)
                print("Checking if run {} just completed.".format(full_run_path))
                if os.path.isfile(full_run_path):
                    fname = pathlib.Path(full_run_path)
                    filetime = datetime.datetime.fromtimestamp(fname.stat().st_ctime)
                    if filetime >= time_minutes_ago:
                        print("Found run ready to demux:", full_run_path)
                        ready_to_demux.append(full_run_path)                
        else:
            print("Path does not exist:", sequencer_path)
    return ready_to_demux