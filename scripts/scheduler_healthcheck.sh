# restart the scheduler if it is not running
pidof 'airflow scheduler' > /dev/null && echo "Service is running" || bash /igo/work/igo/igo-demux/scripts/scheduler_start.sh