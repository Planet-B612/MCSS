dataset=("facebook" "dblp" "flickr" "nethept" "epinions" "youtube" "pokec" "orkut" "livejournal" "friendster" "DBLP_sym" "Youtube_sym" "twitter" "citeseer" "Flickr_sym" "wikitalk" "wikitalkar" "sample" "Twitter")
#18: twitter

# if [[ "$1" -lt 0 ]] || [[ "$1" -ge ${#dataset[@]} ]]; then
#     echo "Error: Invalid dataset index"
#     exit 1
# fi
# if [[ "$1" -lt 0 ]] || [[ "$1" -ge ${#dataset[@]} ]]; then
#     echo "Error: Invalid dataset index"
#     exit 1
# fi

# if ! [[ "$2" =~ ^[0-9]+$ ]]; then
#     echo "Error: eta must be a number"
#     exit 1
# fi
# if ! [[ "$2" =~ ^[0-9]+$ ]]; then
#     echo "Error: eta must be a number"
#     exit 1
# fi
cur_date=$( date +"%m-%d")

# batch=2
eps=0.9
eta=$1
Rand_cost=0  #bool Rnd_cost = false;
batch=4
k=$((eta/2000))
set="tw2"
if [[ "$1" -eq 18 ]]; then
    batch=4
    eps=0.9
fi
start_time=3
end_time=5
model="LT"
if [[ "$model" == "LT" ]]; then
    echo "Should be using LT model."
fi
delta_amp=1
# OUTPUT="/home/cfeng/mRR_Regen/code/log_mine/mine_${dataset[$1]}_${eta}_b${batch}_eps${eps}_${start_time}_${end_time}_${cur_date}.log"
OUTPUT="/home/cfeng/mRR_Regen/code/log_mine/mine_${dataset[18]}_${eta}_b${batch}_eps${eps}_${model}_${cur_date}.log"
sudo cset proc -s "$set" -e -- ./run -dataset_No "18" -eta_0 "$eta" -model "$model" -batch "$batch" -eps "$eps" -times "$start_time" -run_times "$end_time" -Rand_cost "$Rand_cost" -regen "0.6" | tee -a "$OUTPUT"
# ./run -dataset_No "18" -eta_0 "$eta" -model "$model" -batch "$batch" -eps "$eps" -times "$start_time" -run_times "$end_time" -Rand_cost "$Rand_cost" -regen "0.6" | tee -a "$OUTPUT"