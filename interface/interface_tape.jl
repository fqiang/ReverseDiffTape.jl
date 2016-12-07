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
    obj_tt::Tape{Int,Float64}
    lag_tt::Tape{Int,Float64}
    constr_tt::Vector{Tape{Int,Float64}}
    nl_idxes::Vector{Int}
    init::Int
    enable_timing_stats::Bool
    
    jac_I::Vector{Int}
    jac_J::Vector{Int}
    jac_nnz::Int

    laghess_I::Vector{Int}
    laghess_J::Vector{Int}
    laghess_nnz::Int

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
    stk = Vector{Float64}()
    vals = Vector{Float64}()
    imm = Vector{Float64}()
    e = TapeNLPEvaluator(nlpe, 
        numVar,numConstr,0,
        Vector{Float64}(), #parameter values
        
        Tape{Int,Float64}(stk,vals,imm), #objective tape
        Tape{Int,Float64}(stk,vals,imm), #lag_tt
        Vector{Tape{Int,Float64}}(),  #constraints tape

        Vector{Int}(), #nl_idxes
        -1, #init
        with_timing,
        
        Vector{Int}(), Vector{Int}(),-1, #Jacobian
        
        Vector{Int}(), Vector{Int}(),-1, #Hessian

        0.0, 0.0, 0.0, 0.0,    #timer
        0.0, 0.0, 0.0, 0.0     #
        ,0.0                    #jump timer
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
    # @show "TapeNLPEvaluator - initialize"
    jd = d.jd
    gnvar = MathProgBase.numvar(jd.m)

    @timing d.enable_timing_stats tic()
    MathProgBase.initialize(jd,[:ExprGraph])
    @timing d.enable_timing_stats d.jd_init += toq()
    # @show "JuMP evaluator initialize ", d.jd_init

    #let's building up the tape
    objexpr = MathProgBase.obj_expr(jd)
    # @show objexpr
    
    #@assert length(d.pvals) == 0  ## assume no values in pvals
    # @show "build objective"
    @timing d.enable_timing_stats tic()
    tapeBuilder(objexpr,d.obj_tt,d.pvals, gnvar)
    @timing d.enable_timing_stats d.tape_build += toq()
    # @show d.pvals

    #compute gradient structure for f
    # grad_structure(d.obj_tt)

    for i =1:d.numConstr
        conexpr = MathProgBase.constr_expr(jd,i)
        if !MathProgBase.isconstrlinear(jd,i) 
            push!(d.nl_idxes,i)
        end
        j = length(conexpr.args)==3?1:3
        tt = Tape{Int,Float64}(d.obj_tt.stk,d.obj_tt.vals,d.obj_tt.imm)
        @timing d.enable_timing_stats tic()
        tapeBuilderSimple(conexpr.args[j],tt,d.pvals, gnvar)
        @timing d.enable_timing_stats d.tape_build += toq()
        push!(d.constr_tt,tt)
    end

    #build single lag
    d.lag_start = length(d.pvals) + 2
    mergeTapes(d.lag_tt,d.obj_tt, d.constr_tt,d.pvals, gnvar)
    @show d.lag_tt.depth
    @assert (d.lag_start + length(d.constr_tt) - 1) == length(d.pvals)
    #end building single lag

    resizeWorkingMemory(d.obj_tt,d.lag_tt,d.constr_tt)
    

    # @show "TapeNLPEvaluator initialize ",d.tape_build
    d.init = 1  #set initialized 
    nothing
end

MathProgBase.features_available(d::TapeNLPEvaluator) = [:Grad, :Jac, :Hess]

#evaluate the objective function given on iterate x
function MathProgBase.eval_f(d::TapeNLPEvaluator, x)
    @timing d.enable_timing_stats tic()
    v = feval(d.obj_tt,x,d.pvals)
    @timing d.enable_timing_stats d.eval_f_time += toq()
    
    return v
end

#evaluate the objective gradient on given iterate x. Results is set to g
function MathProgBase.eval_grad_f(d::TapeNLPEvaluator, g, x)
    fill!(g,0.0)
        
    @timing d.enable_timing_stats tic()
    grad_reverse_dense(d.obj_tt,x,d.pvals,g)
    @timing d.enable_timing_stats d.eval_grad_f_time += toq()
    
    return nothing
end

#constraint evaluation
function MathProgBase.eval_g(d::TapeNLPEvaluator, g, x)
    @timing d.enable_timing_stats tic()
    @inbounds for i=1:d.numConstr
        @inbounds g[i]=feval(d.constr_tt[i],x,d.pvals)
    end
    @timing d.enable_timing_stats d.eval_g_time += toq()
    
    return nothing 
end


function MathProgBase.jac_structure(d::TapeNLPEvaluator)
    @assert d.jac_nnz == -1
    for i=1:d.numConstr
        @inbounds tt_i = d.constr_tt[i]
        @timing d.enable_timing_stats tic()
        grad_structure(tt_i, d.jac_J)
        @timing d.enable_timing_stats d.jac_structure_time += toq()
        
        v = Vector{Int}(tt_i.nzg)
        fill!(v,i)
        append!(d.jac_I,v) 
    end

    d.jac_nnz = length(d.jac_I)
    @assert length(d.jac_I) == length(d.jac_J) "$(tt_i.nzg) $(tt_i.nvnode) $(length(d.jac_J)) $(length(d.jac_I))"
    return d.jac_I, d.jac_J
end

function MathProgBase.eval_jac_g(d::TapeNLPEvaluator, J, x)
    @assert d.jac_nnz != -1  #structure already computed
    @assert length(J) == d.jac_nnz

    j = 1
    for i = 1:d.numConstr
        @inbounds tt_i = d.constr_tt[i]
        @timing d.enable_timing_stats tic()
        grad_reverse(tt_i,x,d.pvals,j,J)
        j += tt_i.nzg
        @timing d.enable_timing_stats d.eval_jac_g_time += toq()
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
    hess_structure(d.lag_tt,d.laghess_I, d.laghess_J)
    #end single constraint block
    
    d.laghess_nnz = length(d.laghess_I)
    @timing true d.hesslag_structure_time += toq()
   
    @assert d.laghess_nnz == length(d.laghess_J) == length(d.laghess_I) == d.lag_tt.nzh
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
    prepare_reeval_hess(d.lag_tt)
    @inbounds d.pvals[d.lag_start-1] = obj_factor 

    j=1
    @simd for i=d.lag_start:length(d.pvals)
        @inbounds d.pvals[i] = lambda[j]
        j+=1
    end
    @assert (j-1) == d.numConstr 
    hess_reverse(d.lag_tt,x,d.pvals,1,H)
    #end with single constraint block

    @timing true d.eval_hesslag_time += toq()
    return
end

end #end module
#########

