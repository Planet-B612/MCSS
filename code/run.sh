 
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

dataset=("facebook" "dblp" "flickr" "nethept" "epinions" "youtube" "pokec" "orkut" "livejournal" "friendster" "DBLP_sym" "Youtube_sym" "twitter" "citeseer" "Flickr_sym" "wikitalk" "wikitalkar")
data_No=4 # 10: DBLP_sym, 11: Youtube_sym, 4: epinions, 8: livejournal
Times=(0 1 2 3 4 5 6 7 8 9)
set=set1
ratio=1e-4
batch=4
eta=10000
OUTPUT=../log_Q/our_${dataset[$data_No]}_${eta}_${ratio}_b${batch}
{
    for times in ${Times[@]}
    do
        sudo cset proc -s $set -e -- ./run -dataset_No $data_No -eta_0 $eta -batch $batch -times $times -q_ratio $ratio
    done
}|tee -a $OUTPUT
