# restart the scheduler if it is not running
pidof 'airflow scheduler' > /dev/null && echo "Service is running" || bash /igo/work/igo/igo-delivery/scripts/scheduler_start.sh