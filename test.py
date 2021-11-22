from airflow.decorators import dag, task
from airflow.models import DagBag
from demux import *

def test_no_import_errors():
    dag_bag = DagBag()
    dag = dag_bag.get_dag(dag_id='demultiplex')
    print(dag)
    assert len(dag_bag.import_errors) == 0, "No Import Failures"
    assert dag is not None

test_no_import_errors()