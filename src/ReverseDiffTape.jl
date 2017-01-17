module ReverseDiffTape

importall Base

# package code goes here
function __init__()
    # Logging.configure(level=DEBUG)
    @show "loading ReverseDiffTape.jl"
end

export @timing, @asserting, 
#constant
    OP, S_TO_OC,
#types
    TYPE_V, TYPE_P, TYPE_O,
#Objects
    AD, AD_O, AD_P, AD_V, Tape, HessStorage,
#Functions
    append_array,
    tapeBuilder, tapeBuilderNoHess,appendMultParam,appendTapeMultParam,buildSumTape, #building tape from Julia expression 
    mergeTapes,getMaxWorkingSize,resizeHessStorage,
    feval, 
    grad_reverse, grad_reverse_dense, grad_structure,
    hess_structure, hess_reverse, reset_hess,
    tape_report_mem


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


end # module
