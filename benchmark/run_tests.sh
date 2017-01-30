opt=$1 #1 - jump, 2 - ep
julia opf_me.jl 662 $opt
julia arrowhead.jl 2000 8 $opt
julia -e 'include("data_gen.jl"); genRandData(40,2);'
julia random_sparsity.jl 40 2 $opt
rm -f randset40_2*
julia -e 'include("data_gen.jl"); genLogData(2,200);'
julia logmod.jl 2 200 $opt
rm -f log2_200*