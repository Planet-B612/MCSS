# seed_num=200
# MC_round=10000
for del_num in 10000 20000
do
    ./run -del_num $del_num -seed_num 200 -MC_round 10000 -eps_for_verification 0.1 -eta_for_verification 10000
done