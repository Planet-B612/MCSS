dataset=("facebook" "dblp" "flickr" "nethept" "epinions" "youtube" "pokec" "orkut" "livejournal" "friendster" "DBLP_sym" "Youtube_sym" "twitter" "citeseer" "Flickr_sym" "wikitalk" "wikitalkar" "sample" "Twitter")
#18: twitter

if [[ "$1" -lt 0 ]] || [[ "$1" -ge ${#dataset[@]} ]]; then
    echo "Error: Invalid dataset index"
    exit 1
fi

if ! [[ "$2" =~ ^[0-9]+$ ]]; then
    echo "Error: eta must be a number"
    exit 1
fi
cur_date=$( date +"%m-%d")

# batch=2
eps=0.7
eta=$2
batch=2
k=$((eta/2000))
set="lj9"
if [[ "$1" -eq 4 ]]; then
    batch=2
    set="epinions${k}"
fi
if [[ "$1" -eq 10 ]]; then
    batch=2
    set="dblp${k}"
fi
if [[ "$1" -eq 11 ]]; then
    batch=8
    set="mcss${k}"
fi
if [[ "$1" -eq 8 ]]; then
    k=$((eta/10000))
    batch=8
    set="lj${k}"
fi
if [[ "$1" -eq 18 ]]; then
    batch=4
    eps=0.7
fi
# set="youtube5"
delta_amp=10
start_time=0
end_time=10
model="IC"
if [[ "$model" == "IC" ]]; then
    echo "Should be using IC model."
fi
# OUTPUT="/home/cfeng/mRR_Regen/code/log_mine/mine_${dataset[$1]}_${eta}_b${batch}_eps${eps}_${start_time}_${end_time}_${cur_date}.log"
OUTPUT="/home/cfeng/mRR_Regen/code/log_mine/mine_${dataset[$1]}_${eta}_b${batch}_eps${eps}_amp${delta_amp}_${model}_${cur_date}.log"
sudo cset proc -s "$set" -e -- ./run -dataset_No "$1" -eta_0 "$2" -batch "$batch" -eps "$eps" -times "$start_time" -run_times "$end_time" -Rand_cost "0" -regen "0.6" -delta_amp "$delta_amp" | tee -a $OUTPUT