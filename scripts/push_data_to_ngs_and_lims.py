#!/usr/bin/env python3

import sys
import time
from subprocess import call


def main(sequencer, run):
	
	# lets push txt data files to NGS, and then to LIMS
	delphi_endpoint = "curl http://delphi.mskcc.org:8080/ngs-stats/picardstats/updaterun/" + sequencer + "/" + run
	print(delphi_endpoint)
	call(delphi_endpoint, shell = True)

	# give time for NGS to post data
	time.sleep(100)
	lims_endpoint = "curl -k https://igolims.mskcc.org:8443/LimsRest/updateLimsSampleLevelSequencingQc?runId=" + run
	print(lims_endpoint)
	call(lims_endpoint, shell = True)
	
	

if __name__ == "__main__":
	
	sequencer = sys.argv[1]
	run = sys.argv[2]
	
	main(sequencer, run)
	