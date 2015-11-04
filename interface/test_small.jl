
using TapeInterface
using JuMP
using Ipopt
m = Model(solver=TapeSolver(IpoptSolver()))
# m = Model(solver=IpoptSolver())
@defVar(m, 1<=x<=10)
@defVar(m, 1<=y<=10)
@setNLObjective(m, Min, 0.5*x*y-x/y)
@addNLConstraint(m,sin(x*y) <= 100)
@addNLConstraint(m,x^2*y^2 <= 100)
solve(m)
