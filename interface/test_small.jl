
include("interface.jl")

using TapeInterface
using JuMP
using Ipopt
N     = 1000
alpha = 350

m = Model(solver=TapeSolver(IpoptSolver()))

@defVar(m, -100 <= x0 <= 100)
@defVar(m, -100 <= x[1:N] <= 100)
@defVar(m, -100 <= y[1:N] <= 100)
@defVar(m, -100 <= z[1:N] <= 100)


@setNLObjective(m, Min, x0*(sum{ sin(x[i])*cos(y[i]) + sin(x[i]+z[i]), i = 1:N}))
@addNLConstraint(m, x0*(sum{ x[i] , i=1:N})  <= 100)
@addNLConstraint(m, x0*(sum{ y[i] , i=1:N})  <= 100)
@addNLConstraint(m, x0*(sum{ z[i] , i=1:N})  <= 100)

solve(m)


include("comparison.jl")
