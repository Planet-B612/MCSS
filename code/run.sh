dataset=("facebook" "dblp" "flickr" "nethept" "epinions" "youtube" "pokec" "orkut" "livejournal" "friendster" "DBLP_sym" "Youtube_sym" "twitter" "citeseer" "Flickr_sym" "wikitalk" "wikitalkar" "sample" "Twitter")
#18: twitter

# if [[ "$1" -lt 0 ]] || [[ "$1" -ge ${#dataset[@]} ]]; then
#     echo "Error: Invalid dataset index"
#     exit 1
# fi

# if ! [[ "$2" =~ ^[0-9]+$ ]]; then
#     echo "Error: eta must be a number"
#     exit 1
# fi
cur_date=$( date +"%m-%d")

# batch=2
mRR_time_test=0
dataset_no=4
eps=0.7
eta=$1
batch=2
# set="youtube5"
start_time=0
end_time=10
set="tw"
Rand_cost=0  #bool Rnd_cost = false;
delta_amp=1
model="IC"
adapt_IM=1
if [[ "$model" == "IC" ]]; then
    echo "Should be using IC model."
fi
# OUTPUT="/home/cfeng/mRR_Regen/code/log_mine/mine_${dataset[$1]}_${eta}_b${batch}_eps${eps}_${start_time}_${end_time}_${cur_date}.log"
OUTPUT="/home/cfeng/mRR_Regen/code/log_mine/mine_${dataset[$dataset_no]}_${eta}_b${batch}_eps${eps}_amp${delta_amp}_${model}_${cur_date}.log"
# ./run -dataset_No "$dataset_no" -eta_0 "$eta" -batch "$batch" -eps "$eps" -times "$start_time" -run_times "$end_time" -Rnd_cost "$Rand_cost" -regen "0.6" -delta_amp "$delta_amp" -adapt_IM "$adapt_IM" | tee -a "$OUTPUT"
./run -dataset_No "$dataset_no" -eta_0 "$eta" -batch "$batch" -eps "$eps" -times "$start_time" -run_times "$end_time" -Rnd_cost "$Rand_cost" -regen "0.6" -delta_amp "$delta_amp" -mRR_time_test "$mRR_time_test" | tee -a "$OUTPUT"



# sudo cset proc -s "$set" -e -- ./run -dataset_No "8" -eta_0 "20000" -batch "$batch" -eps "$eps" -times "$start_time" -run_times "$end_time" -Rand_cost "0" -regen "0.6" -delta_amp "$delta_amp" | tee -a "$OUTPUT"
# ./runlog -dataset_No "18" -eta_0 "$eta" -model "$model" -batch "$batch" -eps "$eps" -times "$start_time" -run_times "$end_time" -Rnd_cost "$Rand_cost" -regen "0.6" | tee -a "$OUTPUT"