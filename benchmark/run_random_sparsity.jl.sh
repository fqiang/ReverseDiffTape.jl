#run_random_sparsity.jl.run

#generate file
# julia -e 'include("data_gen.jl"); n=4000; for i=1:6; genRandData(n,2^i); end'
# julia -e 'include("data_gen.jl"); k=32; n=500; for i=1:6; genRandData(n,k);n=n*2; end'

opt=$1  #1 - jump, 2 - ep
n=4000
for k in 2 4 8 16 32 64; do 
    if (( $opt == 1 )); then
        echo "JuMP - random_sparsity.jl - $n $k"
    elif (( $opt == 2 )); then
        echo "Tape - random_sparsity.jl - $n $k"
    fi
    julia random_sparsity.jl $n $k $opt
    echo "====================================== "
done 



#k=32
#for n in 500 1000 2000 4000 8000 16000; do
#    if (( $opt == 1 )); then
#        echo "JuMP - random_sparsity.jl - $n $k"
#    elif (( $opt == 2 )); then
#        echo "Tape - random_sparsity.jl - $n $k"
#    fi
#    julia random_sparsity.jl $n $k $opt
#    echo "====================================== "
#done 
