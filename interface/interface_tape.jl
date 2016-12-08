#tape interface

module TapeInterface

using ReverseDiffTape
import MathProgBase

const EPS=1e-5


export TapeSolver

###############################################################
#
#  TapeSolver
#
################################################################
type TapeSolver <: MathProgBase.AbstractMathProgSolver
    s::MathProgBase.AbstractMathProgSolver
end

type TapeData <: MathProgBase.AbstractMathProgModel
    evaluator   ##cannot be typed ? why
    m::MathProgBase.AbstractMathProgModel
end


function MathProgBase.NonlinearModel(solver::TapeSolver)
    # @show "TapeSolver model"
    return TapeData(nothing,MathProgBase.NonlinearModel(solver.s))
end 

function MathProgBase.setwarmstart!(m::TapeData, x) 
    return MathProgBase.setwarmstart!(m.m,x)
end

function MathProgBase.status(m::TapeData)
    e = m.evaluator
    @printf "hesslag_structure_time= %s \n" e.hesslag_structure_time
    @printf "eval_hesslag_time= %s \n" e.eval_hesslag_time
    return MathProgBase.status(m.m)
end

function MathProgBase.getobjval(m::TapeData)
    return MathProgBase.getobjval(m.m)
end

function MathProgBase.getsolution(m::TapeData)
    return MathProgBase.getsolution(m.m)
end

function MathProgBase.getreducedcosts(m::TapeData)
    return MathProgBase.getreducedcosts(m.m)
end

function MathProgBase.getconstrduals(m::TapeData)
    return MathProgBase.getconstrduals(m.m)
end

function MathProgBase.loadproblem!(m::TapeData, numVar, numConstr, l, u, lb, ub, sense, jd::MathProgBase.AbstractNLPEvaluator)
    # @show "TapeData - loadproblem"
    tape_evaluator = TapeNLPEvaluator(jd,numVar,numConstr)
    m.evaluator = tape_evaluator
    MathProgBase.loadproblem!(m.m,numVar,numConstr,l,u,lb,ub,sense,tape_evaluator)
end

function MathProgBase.optimize!(m::TapeData)
    MathProgBase.optimize!(m.m)
end

###############################################################
#
#  TapeEvaluator
#
################################################################
type TapeNLPEvaluator <: MathProgBase.AbstractNLPEvaluator
    jd::MathProgBase.AbstractNLPEvaluator
    numVar::Int
    numConstr::Int
    lag_start::Int
    pvals::Vector{Float64}
    lag_tt::Tape{Int,Float64}
    ttstarts::Vector{Int}
    ttends::Vector{Int}
    trends::Vector{Int}
    hs::HessStorage{Int,Float64}
    
    nl_idxes::Vector{Int}
    init::Int
    enable_timing_stats::Bool
    
    jac_I::Vector{Int}
    jac_J::Vector{Int}
    jac_nnz::Int

    laghess_I::Vector{Int}
    laghess_J::Vector{Int}
    laghess_nnz::Int

    #working array
    stk::Vector{Float64}
    vals::Vector{Float64}
    imm::Vector{Float64}

    # timers
    tape_build::Float64
    eval_f_time::Float64
    eval_g_time::Float64
    eval_grad_f_time::Float64
    jac_structure_time::Float64
    eval_jac_g_time::Float64
    hesslag_structure_time::Float64
    eval_hesslag_time::Float64

    #jump nlp evaluator timer
    jd_init::Float64
end

function TapeNLPEvaluator(nlpe::MathProgBase.AbstractNLPEvaluator,numVar,numConstr;with_timing=false)
    e = TapeNLPEvaluator(nlpe, 
        numVar,numConstr,0,
        Vector{Float64}(), #parameter values
        
        Tape{Int,Float64}(), #lag_tt
        Vector{Int}(),
        Vector{Int}(),
        Vector{Int}(),
        HessStorage{Int,Float64}(0,-1),

        Vector{Int}(), #nl_idxes
        -1, #init
        with_timing,
        
        Vector{Int}(), Vector{Int}(),-1, #Jacobian
        Vector{Int}(), Vector{Int}(),-1, #Hessian

        Vector{Float64}(),Vector{Float64}(),Vector{Float64}(), #working array


        0.0, 0.0, 0.0, 0.0,    #timer
        0.0, 0.0, 0.0, 0.0     
        ,0.0                    
        ) 
    finalizer(e, evaluator_cleanup)
    return e  
end

function evaluator_cleanup(e::TapeNLPEvaluator)
    @timing e.enable_timing_stats begin
        f = open("tape_build.txt","w")
        writedlm(f,e.tape_build)
        close(f)

        f = open("eval_f_time.txt","w")
        writedlm(f, e.eval_f_time)
        close(f)

        f = open("eval_g_time.txt","w")
        writedlm(f, e.eval_g_time)
        close(f)
        
        f = open("eval_grad_f_time.txt","w")
        writedlm(f,e.eval_grad_f_time)
        close(f)

        f = open("jac_structure_time.txt","w")
        writedlm(f, e.jac_structure_time)
        close(f)

        f = open("eval_jac_g_time.txt","w")
        writedlm(f, e.eval_jac_g_time)
        close(f)

        f = open("hesslag_structure_time.txt","w")
        writedlm(f, e.hesslag_structure_time)
        close(f)

        f = open("eval_hesslag_time.txt","w")
        writedlm(f, e.eval_hesslag_time)
        close(f)
    end
end


function MathProgBase.initialize(d::TapeNLPEvaluator, requested_features::Vector{Symbol})
    @assert d.init == -1
    # @assert length(d.pvals) == length(d.ttstarts) == length(d.ttends) == length(d.trends) == 0  ## assume no values in pvals
    jd = d.jd
    gnvar = MathProgBase.numvar(jd.m)

    @timing d.enable_timing_stats tic()
    MathProgBase.initialize(jd,[:ExprGraph])
    @timing d.enable_timing_stats d.jd_init += toq()
    
    objexpr = MathProgBase.obj_expr(jd)

#indexing approach
    pholder = Vector{Int}()
    bigT = tapeBuilderNoHess(objexpr,d.pvals,gnvar)
    push!(d.ttstarts, 1)
    push!(d.ttends,length(bigT.tt))
    push!(d.trends,length(bigT.tr))
    appendMultParam(bigT,pholder)   

    ntts = 1
    obj_nvnode = bigT.nvnode

    for i = 1:d.numConstr
        conexpr = MathProgBase.constr_expr(jd,i)
        if !MathProgBase.isconstrlinear(jd,i) 
            push!(d.nl_idxes,i)
        end
        j = length(conexpr.args) == 3?1:3
        @timing d.enable_timing_stats tic()
        tt = tapeBuilderNoHess(conexpr.args[j],d.pvals, gnvar)
        @timing d.enable_timing_stats d.tape_build += toq()
        appendTapeMultParam(bigT, tt, d.pvals,d.ttstarts,d.ttends,d.trends,pholder)
        ntts +=1
    end
    @assert length(pholder) == 1+d.numConstr

    j = length(d.pvals)
    d.lag_start = j + 2
    for i=1:length(pholder)
        @assert bigT.tt[pholder[i]] == -2
        bigT.tt[pholder[i]] = j + i
    end
    append!(d.pvals,ones(Float64,ntts))
    @assert length(d.pvals) == (d.lag_start-1+d.numConstr)
    d.lag_tt = buildSumTape(bigT.tt, bigT.depth,ntts,gnvar,d.hs)
    @assert length(d.ttstarts) == length(d.ttends) == length(d.trends) == (1 + d.numConstr)
##### otherwise
    
    # @timing d.enable_timing_stats tic()
    # d.obj_tt = tapeBuilderNoHess(objexpr,d.pvals, gnvar)
    # @timing d.enable_timing_stats d.tape_build += toq()
    # # @show d.pvals

    # for i =1:d.numConstr
    #     conexpr = MathProgBase.constr_expr(jd,i)
    #     if !MathProgBase.isconstrlinear(jd,i) 
    #         push!(d.nl_idxes,i)
    #     end
    #     j = length(conexpr.args)==3?1:3
    #     @timing d.enable_timing_stats tic()
    #     tt = tapeBuilderNoHess(conexpr.args[j],d.pvals, gnvar)
    #     @timing d.enable_timing_stats d.tape_build += toq()
    #     push!(d.constr_tt,tt)
    # end

    # #build single lag
    # d.lag_start = length(d.pvals) + 2
    # d.lag_tt = mergeTapes(d.obj_tt, d.constr_tt,d.pvals, gnvar, d.hs)
    @show d.lag_tt.depth
    @assert (d.lag_start + d.numConstr - 1) == length(d.pvals)  "$(length(d.pvals)) $(d.numConstr) $(d.lag_start)"
    #end building single lag

    resize!(d.stk, d.lag_tt.depth + d.lag_tt.maxoperands)
    resize!(d.vals, d.lag_tt.nnode)
    resize!(d.imm, d.lag_tt.immlen)


    # mx_vallen, mx_depth, mx_ops, mx_immlen = getMaxWorkingSize(d.obj_tt,d.lag_tt,d.constr_tt)
    # @assert mx_vallen == d.lag_tt.nnode && mx_depth == d.lag_tt.depth && mx_immlen == d.lag_tt.immlen && mx_ops == d.lag_tt.maxoperands


    # resize!(d.jac_I,(d.lag_tt.nvnode - d.obj_tt.nvnode))
    # resize!(d.jac_J,(d.lag_tt.nvnode - d.obj_tt.nvnode))

    resize!(d.jac_I,(d.lag_tt.nvnode - obj_nvnode))
    resize!(d.jac_J,(d.lag_tt.nvnode - obj_nvnode))    

    # @show "TapeNLPEvaluator initialize ",d.tape_build
    d.init = 1  #set initialized 
    nothing
end

MathProgBase.features_available(d::TapeNLPEvaluator) = [:Grad, :Jac, :Hess]

#evaluate the objective function given on iterate x
function MathProgBase.eval_f(d::TapeNLPEvaluator, x)
    @timing d.enable_timing_stats tic()
    tt = d.lag_tt.tt
    @inbounds ts = d.ttstarts[1]
    @inbounds te = d.ttends[1]    
    v = feval(ts,te,tt,x,d.pvals,d.stk)
    @timing d.enable_timing_stats d.eval_f_time += toq()
    
    return v
end

#evaluate the objective gradient on given iterate x. Results is set to g
function MathProgBase.eval_grad_f(d::TapeNLPEvaluator, g, x)
    # @timing d.enable_timing_stats tic()
    tt = d.lag_tt.tt
    @inbounds ts = d.ttstarts[1]
    @inbounds te = d.ttends[1]
    grad_reverse_dense(ts,te,tt,x,d.pvals,d.vals,d.stk, d.imm, g)
    # @timing d.enable_timing_stats d.eval_grad_f_time += toq()
    
    return nothing
end

#constraint evaluation
function MathProgBase.eval_g(d::TapeNLPEvaluator, g, x)
    @timing d.enable_timing_stats tic()
    tt = d.lag_tt.tt
    @inbounds for i=1:d.numConstr
        # @inbounds tt = d.constr_tt[i].tt
        @inbounds ts = d.ttstarts[i+1]
        @inbounds te = d.ttends[i+1]
        @inbounds g[i]=feval(ts,te,tt,x,d.pvals,d.stk)
    end
    @timing d.enable_timing_stats d.eval_g_time += toq()
    
    return nothing 
end


function MathProgBase.jac_structure(d::TapeNLPEvaluator)
    @assert d.jac_nnz == -1
    @assert length(d.jac_I) == length(d.jac_J)
    # @timing d.enable_timing_stats tic()
    tt = d.lag_tt.tt
    start = 1
    for i=1:d.numConstr
        # @inbounds tt = d.constr_tt[i].tt
        @inbounds ts = d.ttstarts[i+1]
        @inbounds te = d.ttends[i+1]
        nz = grad_structure(ts, te, tt, start, d.jac_J)
        @inbounds fill!(sub(d.jac_I, start:(start+nz-1)),i) 
        start += nz
    end
    d.jac_nnz = length(d.jac_I)
    # @timing d.enable_timing_stats d.jac_structure_time += toq()    
    
    @assert length(d.jac_I) == length(d.jac_J) == (start - 1) "$(tt_i.nzg) $(tt_i.nvnode) $(length(d.jac_J)) $(length(d.jac_I))"
    return d.jac_I, d.jac_J
end

function MathProgBase.eval_jac_g(d::TapeNLPEvaluator, J, x)
    @assert d.jac_nnz != -1  #structure already computed
    @assert length(J) == d.jac_nnz

    tt = d.lag_tt.tt
    j = 1
    for i = 1:d.numConstr
        # @inbounds tt = d.constr_tt[i].tt
        @inbounds ts = d.ttstarts[i+1]
        @inbounds te = d.ttends[i+1]
        # @timing d.enable_timing_stats tic()
        nz = grad_reverse(ts,te,tt,x,d.pvals,d.vals,d.stk,d.imm,j,J)
        j += nz
        # @timing d.enable_timing_stats d.eval_jac_g_time += toq()
    end
    @assert j-1 == d.jac_nnz "$j $(d.jac_nnz)"
    return
end

function MathProgBase.eval_hesslag_prod(
    d::TapeNLPEvaluator,
    h::Vector{Float64}, # output vector
    x::Vector{Float64}, # current solution
    v::Vector{Float64}, # rhs vector
    σ::Float64,         # multiplier for objective
    μ::Vector{Float64}) # multipliers for each constraint

    #@assert false 
    MathProgBase.eval_hesslag_prod(d,h,x,v,σ,μ)
end

function MathProgBase.hesslag_structure(d::TapeNLPEvaluator)
    @assert d.laghess_nnz == -1
    @assert length(d.laghess_I) == length(d.laghess_J) == 0
    
    @timing true tic()
      
    #with single constraints block
    tt = d.lag_tt.tt
    tr = d.lag_tt.tr
    nz = hess_structure(1, length(tt), tt, length(tr), tr, d.hs, d.laghess_I, d.laghess_J)
    #end single constraint block
    
    d.laghess_nnz = length(d.laghess_I)
    @timing true d.hesslag_structure_time += toq()
   
    @assert d.laghess_nnz == length(d.laghess_J) == length(d.laghess_I) == nz
    return  d.laghess_I, d.laghess_J
end

function MathProgBase.eval_hesslag(
    d::TapeNLPEvaluator,
    H::Vector{Float64},         # Sparse hessian entry vector
    x::Vector{Float64},         # Current solution
    obj_factor::Float64,        # Lagrangian multiplier for objective
    lambda::Vector{Float64})    # Multipliers for each constraint
   
    @assert (length(d.pvals)-d.lag_start+1) == length(lambda)
    @assert length(H) == d.laghess_nnz    
    
    @timing true tic()
    #with single constraint block
    @inbounds d.pvals[d.lag_start-1] = obj_factor 

    j=1
    @simd for i=d.lag_start:length(d.pvals)
        @inbounds d.pvals[i] = lambda[j]
        j+=1
    end
    @assert (j-1) == d.numConstr 
    tt = d.lag_tt.tt
    tr = d.lag_tt.tr
    nz = hess_reverse(1, length(tt), tt, length(tr), tr ,x,d.pvals,d.stk, d.vals, d.imm, d.hs, 1,H)
    #end with single constraint block

    @timing true d.eval_hesslag_time += toq()
    return
end

end #end module
#########

