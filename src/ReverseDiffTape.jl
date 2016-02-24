module ReverseDiffTape

importall Base

# package code goes here
function __init__()
    # Logging.configure(level=DEBUG)
    @show "loading ReverseDiffTape.jl"
end

export
#constant
    OP, S_TO_OC,
#types
    TYPE_V, TYPE_P, TYPE_O,
#Objects
    AD, AD_O, AD_P, AD_V, EdgeSet, Tape,
#Functions
    append_array,
    tapeBuilder, #building tape from Julia expression 
    feval, 
    grad_reverse,  grad_structure,
    hess_structure_lower, hess_reverse, clean_hess_eset,
    hess_structure2, hess_reverse2, reset_hess2, prepare_reeval_hess2,
    report_tape_mem


@inline function append_array{I,V}(dest::Vector{V},d_offset::I,src::Vector{V},s_offset::I, n::I)
    for i=1:n
        @inbounds dest[i+d_offset] = src[i+s_offset]
    end
end


include("./types.jl")
include("./operator.jl")
include("./func_eval.jl")
include("./reverse_grad.jl")
include("./reverse_hess_ep.jl")
include("./reverse_hess_ep2.jl")

end # module
