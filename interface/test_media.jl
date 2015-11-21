
# from JuMP example/optcontrol.jl

include("interface.jl")

using TapeInterface
using JuMP
using Ipopt
N     = 5
ni    = N
h     = 1/ni
alpha = 350

m = Model(solver=TapeSolver(IpoptSolver()))

@defVar(m, -1 <= t[1:(ni+1)] <= 1)
@defVar(m, -0.05 <= x[1:(ni+1)] <= 0.05)
@defVar(m, u[1:(ni+1)])

@setNLObjective(m, Min, sum{ 0.5*h*(u[i+1]^2 + u[i]^2) + 0.5*alpha*h*(cos(t[i+1]) + cos(t[i])), i = 1:ni})

# cons1
for i in 1:ni
    @addNLConstraint(m, x[i+1] - x[i] - (0.5h)*(sin(t[i+1])+sin(t[i])) == 0)
end
# cons2
for i in 1:ni
    @addConstraint(m, t[i+1] - t[i] - (0.5h)*u[i+1] - (0.5h)*u[i] == 0)
end

solve(m)

# include("comparion.jl")