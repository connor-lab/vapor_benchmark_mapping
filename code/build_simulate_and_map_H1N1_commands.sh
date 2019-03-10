N_HUMAN=16679
N_SWINE=4054
N_AVIAN=552
for i in {1..16679}; do
        echo "bash code/simulate_and_map_H1N1.sh res/cali0709.fa res/HA_HuH1N1_nf.fa $i results_H1N1/HA_HuH1N1vsCali0709.$i.csv HuH1N1$i"
done

for i in {1..4054}; do
        echo "bash code/simulate_and_map_H1N1.sh res/cali0709.fa res/HA_SwH1N1_nf.fa $i results_H1N1/HA_SwH1N1vsCali0709.$i.csv SwH1N1$i"
done

for i in {1..552}; do
        echo "bash code/simulate_and_map_H1N1.sh res/cali0709.fa res/HA_AvH1N1_nf.fa $i results_H1N1/HA_AvH1N1vsCali0709.$i.csv AvH1N1$i"
done
