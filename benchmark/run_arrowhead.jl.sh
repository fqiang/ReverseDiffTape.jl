#run_arrowhead.jl.sh

opt=$1  #1 - jump, 2 - ep
k=16
for n in 2000 4000 8000 16000 32000 64000; do # 8000 16000 32000 64000; do
    if (( $opt == 1 )); then
        echo "JuMP - arrowhead.jl - $n $k"
    elif (( $opt == 2 )); then
        echo "Tape - arrowhead.jl - $n $k"
    fi
        #statements
    julia arrowhead.jl $n $k $opt
    echo "====================================== "
done 

n=32000
for k in 2 4 8 16 32 64; do
    if (( $opt == 1 )); then
        echo "JuMP - arrowhead.jl - $n $k"
    elif (( $opt == 2 )); then
        echo "Tape - arrowhead.jl - $n $k"
    fi
        #statements
    julia arrowhead.jl $n $k $opt
    echo "====================================== "
done 
