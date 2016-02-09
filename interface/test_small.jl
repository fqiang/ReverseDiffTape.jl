
include("interface.jl")

using TapeInterface
using JuMP
using Ipopt

m = Model(solver=TapeSolver(IpoptSolver()))
# m = Model(solver=IpoptSolver())
@defVar(m, 1<=x<=10)
@defVar(m, 1<=y<=10)
@defVar(m, 1<=z<=10)
@defVar(m, 1<=v<=10)
@setNLObjective(m, Min, -0.5*x*y-x*y+x*z*v)
@addNLConstraint(m,sin(x*y)-v+cos(x*z) <= 100)
@addNLConstraint(m,x^2*y^2+y*z-v <= 100)
solve(m)


include("comparison.jl")
