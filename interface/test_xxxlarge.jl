
include("interface.jl")

#
#  An NLP model from JuMP
#

using TapeInterface
using JuMP
using Ipopt

ns = Vector();
ks = Vector();
jump = Vector();
tape = Vector();
nvar = Vector();
nnzs = Vector();
nnzs_nor = Vector();
mem = Vector();

N=10000
K=8
for t=1:7
    N     = 2*N
    # N = 1280000
    # K = 2*K
    # ni    = N
    # h     = 1/ni
    # alpha = 350

    m = Model(solver=TapeSolver(IpoptSolver(max_iter = 2)))

    # m = Model(solver=IpoptSolver(max_iter=2))

    # @defVar(m, -1 <= t[1:(ni+1)] <= 1)
    # @defVar(m, -0.05 <= x[1:(ni+1)] <= 0.05)
    # @defVar(m, u[1:(ni+1)])

    # @setNLObjective(m, Min, sum{ 0.5*h*(u[i+1]^2 + u[i]^2) + 0.5*alpha*h*(cos(t[i+1]) + cos(t[i])), i = 1:ni} + sum{x[i+1] - x[i] - (0.5h)*(sin(t[i+1])+sin(t[i])) , i=1:ni} )

    @defVar(m, -100<=x[1:N+K+1]<=100)

    # @setNLObjective(m, Min, sum{ cos(x[i] - x[i-1])^4 , i = 2:N})
    # @setNLObjective(m, Min, sum{ sin( sum{  x[i+j]   ,j=1:k}   )^3   , i=1:N-k})

    # @setNLObjective(m,Min, sum{  sin(  sum{  x[1] - x[i]^2   ,j=1:k} )  , i=1:N })

    @setNLObjective(m, Min, sum{sum{ x[i] * x[j+1], j=i:i+K},  i=1:N})

    # cons1
    # for i in 1:ni
        # @addNLConstraint(m, x[i+1] - x[i] - (0.5h)*(sin(t[i+1])+sin(t[i])) }== 0)
    # end
    # cons2
    # for i in 1:ni
    #     @addConstraint(m, t[i+1] - t[i] - (0.5h)*u[i+1] - (0.5h)*u[i] == 0)
    # end
    @show "do gc now"
    gc()
    solve(m)

    d = m.internalModel.evaluator;
    one_H = ones(Float64,length(d.laghess_I))
    csc = sparse(d.laghess_I,d.laghess_J,one_H,d.numVar,d.numVar)
    assert(d.numVar == m.numCols)

    @show N
    @show K
    @show m.numCols
    @show d.laghess_nnz
    @show length(csc.nzval)
    @show d.eval_hesslag_timer
    @show d.jeval_hesslag_timer
    d_mb = report_mem(d) / 1024 /1024 #mb
    @show d_mb


    push!(ns,N)
    push!(ks,K)
    push!(nvar, m.numCols)
    push!(nnzs, d.laghess_nnz)
    push!(nnzs_nor,length(csc.nzval))
    push!(jump,d.jeval_hesslag_timer)
    push!(tape,d.eval_hesslag_timer)
    push!(mem, d_mb)

    nothing
end

hcat(ns,ks,jump,tape,nvar,nnzs,nnzs_nor,mem)

# include("comparison.jl")
