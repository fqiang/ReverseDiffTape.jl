#run_opf_me.jl.sh

opt=$1 #1 - jump, 2 - ep
for i in 662 6620 66200; do
	if (( $opt == 1 )); then
	   echo "JuMP - opf_me.jl - $i"
       julia opf_me.jl $i $opt 
    elif (( $opt == 2 )); then
       echo "Tape - opf_me.jl - $i"
       julia opf_me.jl $i $opt 
    fi	
    echo "====================================== "
done
