#random_sparsity.jl

# include("interface_tape.jl")
# using TapeInterface

using JuMP
using Ipopt

N = parse(Int,ARGS[1])
K = parse(Int,ARGS[2])
RandSet= round(Int,readdlm(string("randset",N,"_",K,".txt")))


if parse(Int, ARGS[3]) == 1
    m = Model(solver=IpoptSolver())
elseif parse(Int, ARGS[3]) == 2
    include("interface_tape.jl")
    using TapeInterface
    m = Model(solver=TapeSolver(IpoptSolver()))
elseif parse(Int, ARGS[3]) == 3
    include("interface_tape_dict.jl")
    using TapeInterface
    m = Model(solver=TapeSolver(IpoptSolver()))
else
    @assert false
end

# @show RandSet

@variable(m, -100<=x[1:N]<=100)

@NLobjective(m, Min, sum{ (x[i]-1)^2 + prod{x[j], j=RandSet[i,:]}, i=1:N} )

solve(m)

