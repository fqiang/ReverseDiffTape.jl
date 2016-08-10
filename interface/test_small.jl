
# include("interface.jl")

# using TapeInterface
# using JuMP
# using Ipopt
# N     = 300
# M = 100
# alpha = 350
# S = 100

# m = Model(solver=TapeSolver(IpoptSolver()))

# @defVar(m, -100 <= x0[1:N] <= 100)
# @defVar(m, -100 <= x0[1:N] <= 100)
# @defVar(m, -100 <= x[1:M] <= 100)
# @defVar(m, -100 <= y[1:M] <= 100)
# @defVar(m, -100 <= z[1:M] <= 100)


# # @setNLObjective(m, Min, x0*(sum{ sin(x[i])*cos(y[i]) + sin(x[i]+z[i]), i = 1:N}))
# # @setNLObjective(m, Min, x0 + sum{ sin(x[i])*cos(y[i]) + sin(x[i]+z[i]), i = 1:M})

# for i=1:N
#     @addNLConstraint(m, x0*(x[i])  <= 100)
#     @addNLConstraint(m, x0*(y[i])  <= 100)
#     @addNLConstraint(m, x0*(z[i])  <= 100)

#     # @addNLConstraint(m, x0*sum{x[j] ,j=1:M} <= 100)
#     # @addNLConstraint(m, x0*sum{y[j] ,j=1:M} <= 100)
#     # @addNLConstraint(m, x0*sum{z[j] ,j=1:M} <= 100)

#     # @addNLConstraint(m, x0[i]*sum{x[j] ,j=1:M} <= 100)
#     # @addNLConstraint(m, x0[i]*sum{y[j] ,j=1:M} <= 100)
#     # @addNLConstraint(m, x0[i]*sum{z[j] ,j=1:M} <= 100)
# end

# solve(m)


# include("comparison.jl")


include("interface.jl")

using TapeInterface
using JuMP
using Ipopt
N     = 800
alpha = 350

m = Model(solver=TapeSolver(IpoptSolver()))

@defVar(m, -100 <= x0 <= 100)
@defVar(m, -100 <= x[1:N] <= 100)
@defVar(m, -100 <= y[1:N] <= 100)
@defVar(m, -100 <= z[1:N] <= 100)


# @setNLObjective(m, Min, x0*(sum{ sin(x[i])*cos(y[i]) + sin(x[i]+z[i]), i = 1:N}) + x0*(sum{ x[i] , i=1:N}) + x0*(sum{ y[i] , i=1:N}) + x0*(sum{ z[i] , i=1:N}) )
@setNLObjective(m, Min, x0*(sum{ sin(x[i])*cos(y[i]) + sin(x[i]+z[i]), i = 1:N}) + x0*(sum{ x[i] + y[i] +z[i] , i=1:N}) )

# @addNLConstraint(m, x0*(sum{ x[i] , i=1:N})  <= 100)
# @addNLConstraint(m, x0*(sum{ y[i] , i=1:N})  <= 100)
# @addNLConstraint(m, x0*(sum{ z[i] , i=1:N})  <= 100)


solve(m)


include("comparison.jl")