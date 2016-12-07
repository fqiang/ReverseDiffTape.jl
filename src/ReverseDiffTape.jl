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
    AD, AD_O, AD_P, AD_V, Tape,
#Functions
    append_array,
    tapeBuilder, #building tape from Julia expression 
    tapeBuilderSimple, 
    mergeTapes,
    feval, 
    grad_reverse,  grad_structure,
    hess_structure, hess_reverse, reset_hess, prepare_reeval_hess,
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

# warming up 
# @printf "warming up ... \n"
# imm = Vector{Float64}();
# p = Vector{Float64}();
# x = [1.1, 2.2, 3.3];
# tape1 = Tape{Int,Float64}(imm=imm,with_timing=false, bh_type=1);
# tape2 = Tape{Int,Float64}(imm=imm,with_timing=false, bh_type=2);
# expr=parse("exp(x[1]+x[2]+x[3])*exp(x[1]*x[2]*x[3])")
# tapeBuilder(expr,tape1,p)
# tapeBuilder(expr,tape2,p)

# @printf " feval \n"
# # @time feval(tape::Tape{Int,Float64},x::Vector{Float64},p::Vector{Float64});
# feval(tape1::Tape{Int,Float64},x::Vector{Float64},p::Vector{Float64});

# @printf " grad_structure \n"
# # @time grad_structure(tape::Tape{Int,Float64}); 
# tape1.nzg = -one(Int)
# grad_structure(tape1::Tape{Int,Float64}); 

# @printf " grad_reverse \n"
# # @time grad_reverse(tape::Tape{Int,Float64},x::Vector{Float64},p::Vector{Float64});
# grad_reverse(tape1::Tape{Int,Float64},x::Vector{Float64},p::Vector{Float64});

# @printf " hess_structure \n"
# # @time hess_structure2(tape::Tape{Int,Float64});
# reset_hess2(tape1::Tape{Int,Float64});
# reset_hess3(tape2::Tape{Int,Float64})
# hess_structure2(tape1::Tape{Int,Float64});
# hess_structure3(tape2::Tape{Int,Float64})

# @printf " hess_reverse \n"
# prepare_reeval_hess2(tape1::Tape{Int,Float64});
# prepare_reeval_hess3(tape2::Tape{Int,Float64})
# hess_reverse2(tape1::Tape{Int,Float64},x::Vector{Float64},p::Vector{Float64});
# hess_reverse3(tape2::Tape{Int,Float64},x::Vector{Float64},p::Vector{Float64});
# # hess_reverse2(tape::Tape{Int,Float64},x::Vector{Float64},p::Vector{Float64});

# @printf "end warming up ... \n"
#


end # module
