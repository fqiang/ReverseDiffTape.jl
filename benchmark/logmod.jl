#logmod.jl



using JuMP
using Ipopt

M = parse(Int,ARGS[1])
N = parse(Int,ARGS[2])

include(string("log",M,"_",N,"_dat.jl"))

m = Model()

simplify = true
@show simplify
if parse(Int, ARGS[3]) == 1
    m = Model(solver=IpoptSolver(), simplify_nonlinear_expressions=simplify)
    #use jump
elseif parse(Int, ARGS[3]) == 2
    include("interface_tape.jl")
    using TapeInterface
    m = Model(solver=TapeSolver(IpoptSolver()),simplify_nonlinear_expressions=simplify)
elseif parse(Int, ARGS[3]) == 3
    include("interface_tape_dict.jl")
    using TapeInterface
    m = Model(solver=TapeSolver(IpoptSolver()),simplify_nonlinear_expressions=simplify)
else
    @assert false
end

@variable(m, theta[1:N], start=-0.0001)
# @variable(m, theta2[1:N], start=-0.0001)

# @variable(m, sum1[1:M], start=1.0)
# for i = 1:M
#     @constraint(m, y[i] * sum{ theta[j] * x[i,j], j = 1:N } == sum1[i])
# end

# testing expression
# @NLobjective(m, Min, sum{ log(1+1/exp(sum{theta[j], j=1:N})) ,i=1:M})
#@NLobjective(m, Min, sum{ log(1+1/exp(sum{theta[j]+theta2[j], j=1:N})) ,i=1:M})
# @NLobjective(m, Min, sum{ log(1+1/exp(sum{theta[j] + theta2[j], j=1:N})) ,i=1:M})

#original logrithmic regression (no parameter)
#@NLobjective(m, Min, sum{ theta[i]^2, i = 1:N} + sum{ log(1+1/exp(sum{theta[j], j=1:N})), i =1:M})
#original logrithmic regression 
@NLobjective(m, Min, lambda*(sum{ theta[i]^2, i = 1:N}) + sum{ log(1+1/exp(y[i] * sum{ theta[j] * x[i,j]  ,j = 1:N } )), i =1:M})


#remvoing parameter y
#@NLobjective(m, Min, lambda*(sum{ theta[i]^2, i = 1:N}) + sum{ log(1+1/exp(sum{ theta[j] * x[i,j]  ,j = 1:N } )), i =1:M})


#log x  -> x^2
#@NLobjective(m, Min, lambda*(sum{ theta[i]^2, i = 1:N}) + sum{ (1+1/exp(y[i] * sum{ theta[j] * x[i,j]  ,j = 1:N } ))^2, i =1:M})

#exp x - > x^2
#@NLobjective(m, Min, lambda*(sum{ theta[i]^2, i = 1:N}) + sum{ log(1+1/(y[i] * sum{ theta[j] * x[i,j]  ,j = 1:N } )^2), i =1:M})

# 1/exp x -> x^2
# @NLobjective(m, Min, lambda*(sum{ theta[i]^2, i = 1:N}) + sum{ log(1+(y[i] * sum{ theta[j] * x[i,j]  ,j = 1:N } )^2 ), i =1:M})

solve(m)
