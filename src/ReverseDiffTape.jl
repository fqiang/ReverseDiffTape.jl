module ReverseDiffTape

importall Base

# package code goes here
function __init__()
    # Logging.configure(level=DEBUG)
    println("loading ReverseDiffTape.jl")
end

export
#constant
    OP, S_TO_OC,
#types
    TYPE_V, TYPE_P, TYPE_O,
#Objects
    AD, AD_O, AD_P, AD_V, EdgeSet, Tape,
#Functions
    tapeBuilder, #building tape from Julia expression 
    feval, grad_reverse, hess_reverse, grad_structure, hess_structure_lower, clean_hess_eset,append_array


include("./types.jl")
include("./operator.jl")
include("./func_eval.jl")
include("./reverse_grad.jl")
include("./reverse_hess_ep.jl")

end # module
