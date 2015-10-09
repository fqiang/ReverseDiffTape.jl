module ADOLJ

using DataStructures
using Logging

importall Base

# package code goes here
function __init__()
	Logging.configure(level=DEBUG)
	debug("loading ADOLJ.jl")
end

export
#types
	TT_TYPE, TV_TYPE,
#Objects
	AD_P, AD_V, Edge,
#Functions
	feval, grad_reverse, reverse_hess_ep, nzg, nzh


include("./types.jl")
include("./operator.jl")
include("./func_eval.jl")
include("./reverse_grad.jl")
include("./reverse_hess_ep.jl")


end # module
