#arrowhead.jl


using JuMP
using Ipopt

N = parse(Int, ARGS[1])
K = parse(Int, ARGS[2])

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

@variable(m, -100<=x[1:N+K]<=100)

@NLobjective(m, Min, sum{ cos(sum{ x[i+j], j=1:K}) + sum{ (x[i] + x[j])^2, j=1:K} , i=1:N} )

solve(m)

