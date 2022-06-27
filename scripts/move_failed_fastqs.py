"""

"""

STATS_DIR = "/igo/staging/stats"

def move_failed_fastqs(igo_id, run):
    if not igo_id or not run:
        return "igo_id and run are required arguments."

    return "Completed move_failed_fastqs"



#optionally invoke directly, for example:
#python move_failed_fastqs.py igo_id run
if __name__ == '__main__':
    igo_id = sys.argv[1]
    run = sys.argv[2]
    move_failed_fastqs(igo_id, run)