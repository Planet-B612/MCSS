 
# for j in 1 2 3 4 5 6
# do
# nohup ./run -eta_0 0.08 > stdout.txt 2>stderr.txt &
# nohup ./run -eta_0 0.10 > stdout.txt 2>stderr.txt &
# nohup ./run -eta_0 0.12 > stdout.txt 2>stderr.txt &
# nohup ./run -eta_0 0.14 > stdout.txt 2>stderr.txt &
# nohup ./run -eta_0 0.16 > stdout.txt 2>stderr.txt &
# nohup ./run -eta_0 0.18 > stdout.txt 2>stderr.txt &
#local eta
# eta=$((0.08+$j*0.02))
# ./run -eta_0 0.12
# ./run -eta_0 0.12
# ./run -eta_0 0.12
# ./run -eta_0 0.12
# ./run -eta_0 0.08 -E_PCG 1 
# ./run -eta_0 0.10 -E_PCG 1 
# ./run -eta_0 0.12 -E_PCG 1
# ./run -eta_0 0.14 -E_PCG 1 
# ./run -eta_0 0.16 -E_PCG 1 
# ./run -eta_0 0.18 -E_PCG 1 
# ./run -eta_0 0.08 -E_PCG 0
# ./run -eta_0 0.10 -E_PCG 0 
# ./run -eta_0 0.12 -E_PCG 0 
# ./run -eta_0 0.14 -E_PCG 0 
# ./run -eta_0 0.16 -E_PCG 0 
# ./run -eta_0 0.18 -E_PCG 0 


# ./run -eta_0 0.08
# ./run -eta_0 0.10
# ./run -eta_0 0.12
#./home/cfeng/MCSS/a.out -eta_0 0.08
#./home/cfeng/MCSS/a.out -eta_0 0.10
#./home/cfeng/MCSS/a.out -eta_0 0.12
# ./run -eta_0 0.14
# ./run -eta_0 0.16
# ./run -eta_0 0.18

# done

if [ "$1" -lt 0 ] || [ "$1" -ge ${#dataset[@]} ]; then
    echo "Error: Invalid dataset index"
    exit 1
fi

if ! [[ "$2" =~ ^[0-9]+$ ]]; then
    echo "Error: eta must be a number"
    exit 1
fi

batch=2
eps=0.7
cur_date=$( date +"%m-%d")
eta=$2
k=$((((eta-3000)/500)+1))
set="${dataset[$1]}${k}"
if [ "$1" -eq 10 ]; then
    set="dblp${k}"
fi
if [ "$1" -eq 11 ]; then
    set="youtube${k}"
fi
batch=$3
OUTPUT="log_mine/mine_${dataset[$1]}_${eta}_b${batch}_eps${eps}_${cur_date}.log"
sudo cset proc -s "$set" -e -- ./asm -dataset_No "$1" -eta "$2" -batch "$batch" -epsilon "$eps" -start_time "0" -time "10" -Rand_cost "0"| tee -a $OUTPUT
