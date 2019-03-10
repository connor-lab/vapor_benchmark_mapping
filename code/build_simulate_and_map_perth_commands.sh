PIDS="98 96 94 92 90 88 86 84"
for PID in $PIDS; do
        for rep in {1..1000}; do
                echo "bash code/simulate_and_map_perth.sh res/perth1609.fa res/perth1609.fa 1 $PID results_perth/HA_Perth1609vsPerth1609_mutall.$PID.$rep.csv $PID$rep"
        done;
done;

