from airflow.decorators import dag, task
from airflow.models import DagBag
from find_completed_runs_dag import *

def test_no_import_errors():
    dag_bag = DagBag()
    dag = dag_bag.get_dag(dag_id='find_completed_runs')
    print(dag)
    assert len(dag_bag.import_errors) == 0, "No Import Failures"
    assert dag is not None