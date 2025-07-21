dataset=("facebook" "dblp" "flickr" "nethept" "epinions" "youtube" "pokec" "orkut" "livejournal" "friendster" "DBLP_sym" "Youtube_sym" "twitter" "citeseer" "Flickr_sym" "wikitalk" "wikitalkar")
seed_num=500
MC_round=100000
eta=20000
eps=0.1

dataset_No=8
./run -dataset_No $dataset_No -seed_num $seed_num -MC_round $MC_round -eps_for_verification $eps -eta_for_verification $eta |tee -a ./log/${dataset[$dataset_No]}.log

dataset_No=4
./run -dataset_No $dataset_No -seed_num $seed_num -MC_round $MC_round -eps_for_verification $eps -eta_for_verification $eta|tee -a ./log/${dataset[$dataset_No]}.log

dataset_No=10
./run -dataset_No $dataset_No -seed_num $seed_num -MC_round $MC_round -eps_for_verification $eps -eta_for_verification $eta|tee -a ./log/${dataset[$dataset_No]}.log

dataset_No=11
./run -dataset_No $dataset_No -seed_num $seed_num -MC_round $MC_round -eps_for_verification $eps -eta_for_verification $eta|tee -a ./log/${dataset[$dataset_No]}.log