
module TapeInterface

using ReverseDiffTape
import MathProgBase

const EPS=1e-10

###############################################################
#
#  TapeSolver
#
################################################################
type TapeSolver <: MathProgBase.AbstractMathProgSolver
	s::MathProgBase.AbstractMathProgSolver
end
export TapeSolver

type TapeData <: MathProgBase.AbstractMathProgModel
    evaluator   ##cannot be typed ? why
    m::MathProgBase.AbstractMathProgModel
end


function MathProgBase.model(solver::TapeSolver)
	println("TapeSolver model")
	return TapeData(nothing,MathProgBase.model(solver.s))
end 

function MathProgBase.setwarmstart!(m::TapeData, x) 
    return MathProgBase.setwarmstart!(m.m,x)
end

function MathProgBase.status(m::TapeData)
    # @show m.evaluator
    # @show m.evaluator.eval_f_timer
    # @show m.evaluator.eval_g_timer
    # @show m.evaluator.eval_grad_f_timer
    # @show m.evaluator.eval_jac_g_timer
    # @show m.evaluator.eval_hesslag_timer
    
    # @show m.evaluator.jd.eval_f_timer
    # @show m.evaluator.jd.eval_g_timer
    # @show m.evaluator.jd.eval_grad_f_timer
    # @show m.evaluator.jd.eval_jac_g_timer
    # @show m.evaluator.jd.eval_hesslag_timer
    
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

function MathProgBase.loadnonlinearproblem!(m::TapeData, numVar, numConstr, l, u, lb, ub, sense, jd::MathProgBase.AbstractNLPEvaluator)
    println("TapeData - loadnonlinearproblem")
    # if(jd.eval_f_timer>0.0)
    #     @show typeof(m.evaluator)
    #     @show m.evaluator
    #     tape_evaluator::MathProgBase.SolverInterface.AbstractNLPEvaluator = m.evaluator
    #     println("tape_evaluator arleady set")
    # else
        tape_evaluator = TapeNLPEvaluator(jd,numVar,numConstr)
        m.evaluator = tape_evaluator
        println("tape_evaluator is set")
    # end
    MathProgBase.loadnonlinearproblem!(m.m,numVar,numConstr,l,u,lb,ub,sense,tape_evaluator)
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
    pvals::Vector{Float64}
    obj_tt::Tape{Int,Float64}
    constr_tt::Vector{Tape{Int,Float64}}
    init::Int
    # laghess_tt::Tape{Int,Float64}  #full expr
    # p::Vector{Float64}             #obj factor + lambda

    jac_I::Vector{Int}
    jac_J::Vector{Int}
    jac_nnz::Int

    laghess_I::Vector{Int}
    laghess_J::Vector{Int}
    laghess_nnz::Int

    # timers
    tape_build::Float64
    eval_f_timer::Float64
    eval_g_timer::Float64
    eval_grad_f_timer::Float64
    eval_jac_g_timer::Float64
    eval_hesslag_timer::Float64

    jeval_f_timer::Float64
    jeval_g_timer::Float64
    jeval_grad_f_timer::Float64
    jeval_jac_g_timer::Float64
    jeval_hesslag_timer::Float64


end

function TapeNLPEvaluator(nlpe::MathProgBase.AbstractNLPEvaluator,numVar,numConstr)
 	return TapeNLPEvaluator(nlpe, 
        numVar,numConstr,
        Vector{Float64}(), #parameter values
        
        Tape{Int,Float64}(), #objective tape
        Vector{Tape{Int,Float64}}(),  #constraints tape

        -1, #init

        # Tape{Int,Float64}(),  #full expr
        # Vector{Float64}(),    #obj factor + lambda

        Vector{Int}(), Vector{Int}(),-1, #Jacobian
        
        Vector{Int}(), Vector{Int}(),-1, #Hessian
        0.0,                        #tape build time
        0.0, 0.0, 0.0, 0.0, 0.0,   #my timer
        0.0, 0.0, 0.0, 0.0, 0.0)   #jump timer
end

function MathProgBase.initialize(d::TapeNLPEvaluator, requested_features::Vector{Symbol})
	println("TapeNLPEvaluator - initialize")
	# @show d.eval_f_timer
 #    @show d.eval_grad_f_timer
 #    @show d.eval_g_timer
 #    @show d.eval_jac_g_timer 
 #    @show d.eval_hesslag_timer

    jd = d.jd
    MathProgBase.initialize(jd,[:Grad,:Jac,:ExprGraph,:Hess])
    # @show jd.eval_f_timer
    # @show jd.eval_grad_f_timer
    # @show jd.eval_g_timer
    # @show jd.eval_jac_g_timer 
    # @show jd.eval_hesslag_timer

    # if(d.eval_f_timer > 0.0)
    #     println("return - already init")
    #     return  #already initialized
    # end
    # @show requested_features

	
    # @show jd.eval_f_timer
    # @show jd.eval_grad_f_timer
    # @show jd.eval_g_timer
    # @show jd.eval_jac_g_timer 
    # @show jd.eval_hesslag_timer

	#let's building up the tape
    if(d.init == -1)
        tic()
    	objexpr = MathProgBase.obj_expr(jd)
        # @show objexpr
       
    	assert(length(d.pvals) == 0)  ## assume no values in pvals
    	tapeBuilder(objexpr,d.obj_tt,d.pvals)
        
        # hesslag_expr = AD(d.obj_tt.tt)*AD_P(d.p,1.0)

        for i =1:1:d.numConstr
           conexpr = MathProgBase.constr_expr(jd,i)
           # @show conexpr.args[1]
           # @show dump(conexpr)
           j = length(conexpr.args)==3?1:3
           tt = Tape{Int,Float64}()
           tapeBuilder(conexpr.args[j],tt,d.pvals)
           push!(d.constr_tt,tt)
           # hesslag_expr = hesslag_expr + AD(tt.tt)*AD_P(d.p,1.0)
           # @show d.constr_tt[i]
        end
        # d.laghess_tt = tapeBuilder(hesslag_expr.data)
        @assert length(d.constr_tt) == d.numConstr
        # @assert length(d.p) == 1+d.numConstr

        tprep = toq()
        println("TapeNLPEvaluator - initialize takes ",tprep)
        d.tape_build = tprep

        d.init = 1  #set initialized tapes
    end

    # reset timers
    d.eval_f_timer = 0.0
    d.eval_grad_f_timer = 0.0
    d.eval_g_timer = 0.0
    d.eval_jac_g_timer = 0.0
    d.eval_hesslag_timer = 0.0

    d.jeval_f_timer = 0.0
    d.jeval_grad_f_timer = 0.0
    d.jeval_g_timer = 0.0
    d.jeval_jac_g_timer = 0.0
    d.jeval_hesslag_timer = 0.0

    nothing
end

MathProgBase.features_available(d::TapeNLPEvaluator) = [:Grad, :Jac, :Hess, :HessVec]

#evaluate the objective function given on iterate x
function MathProgBase.eval_f(d::TapeNLPEvaluator, x)
    tic()
    v = feval(d.obj_tt,x,d.pvals)
    d.eval_f_timer += toq()
    
    tic()
    jv = MathProgBase.eval_f(d.jd,x)
    d.jeval_f_timer += toq()
    assert(jv == v)
    return v
end

#evaluate the objective gradient on given iterate x. Results is set to g
function MathProgBase.eval_grad_f(d::TapeNLPEvaluator, g, x)
    # @show x
    assert(length(g) == length(x))
    tape = d.obj_tt
    if tape.nzg==-1
        grad_structure(tape)
    end
    assert(tape.nzg!=-1)

    tic()
    grad_reverse(tape,x,d.pvals)
    d.eval_grad_f_timer += toq()

    #converting to dense
    sparse_g = sparsevec(tape.g_I,tape.g,length(x))
    for i = 1:length(sparse_g)
        @inbounds g[i] = sparse_g[i]
    end
    ####################

    jg = Array{Float64,1}(length(x))
    tic()
    MathProgBase.eval_grad_f(d.jd,jg,x)
    d.jeval_grad_f_timer += toq()
    # @show jg
    
    # temp = Array{Float64,1}(length(g))
    # fill!(temp,EPS)
    assertArrayEqualEps(g,jg)
    return
end

#constraint evaluation
function MathProgBase.eval_g(d::TapeNLPEvaluator, g, x)
    assert(length(g) == d.numConstr)
    tic()
    @inbounds for i=1:d.numConstr
        @inbounds g[i]=feval(d.constr_tt[i],x,d.pvals)
    end
    d.eval_g_timer += toq()
    # @show g

    tic()
    jg = Array{Float64,1}(length(g))
    MathProgBase.eval_g(d.jd,jg,x)
    d.jeval_g_timer += toq()
    # @show jg

    # temp = Array{Float64,1}(length(g))
    # fill!(temp,EPS)
    assertArrayEqualEps(g,jg)
    return
end


function MathProgBase.jac_structure(d::TapeNLPEvaluator)
    # @show "MathProgBase.jac_structure"
    if(d.jac_nnz != -1) 
        return d.jac_I, d.jac_J
    end

    for i=1:d.numConstr
        tape_i = d.constr_tt[i]
        grad_structure(tape_i)
        v = Vector{Int}(length(tape_i.g_I))
        fill!(v,i)
        append!(d.jac_I,v)
        append!(d.jac_J,tape_i.g_I)
    end
    # @show I
    # @show J
    
    assert(length(d.jac_I) == length(d.jac_J))
    V = ones(Float64,length(d.jac_I))
    csc = sparse(d.jac_I,d.jac_J,V)

    (jI, jJ) = MathProgBase.jac_structure(d.jd)
    assert(length(jI) == length(jJ))
    jV = ones(Float64,length(jI))
    jcsc = sparse(jI,jJ,jV)
    # @show jI
    # @show jJ
    
    # @show csc
    # @show jcsc
    # @show csc.colptr
    # @show jcsc.colptr
    assert(csc.colptr == jcsc.colptr)
    assert(csc.rowval == jcsc.rowval)
    assert(csc.m == jcsc.m)
    assert(csc.n == jcsc.n)

    d.jac_nnz = length(d.jac_I)
    return d.jac_I, d.jac_J
end

function MathProgBase.eval_jac_g(d::TapeNLPEvaluator, J, x)
    # @show x
    assert(d.jac_nnz != -1)  #structure already computed
    assert(length(J) == d.jac_nnz)
    
    tic()
    J_len = 0
    @inbounds for i = 1:d.numConstr
        tape_i = d.constr_tt[i]
        grad_reverse(tape_i,x,d.pvals)
        append_array(J,J_len,tape_i.g,0,length(tape_i.g))
        J_len += length(tape_i.g)
    end
    d.eval_jac_g_timer += toq()
    
    # @show J
    csc = sparse(d.jac_I,d.jac_J,J)

    jJ = zeros(Float64,length(d.jd.jac_I))
    tic()
    MathProgBase.eval_jac_g(d.jd,jJ,x)
    d.jeval_jac_g_timer += toq()

    # @show jJ
    jcsc = sparse(d.jd.jac_I,d.jd.jac_J,jJ)

    # @show csc.nzval
    # @show jcsc.nzval
    assert(csc.m == jcsc.m)
    assert(csc.n == jcsc.n)

    # assert(csc.colptr == jcsc.colptr)
    # assert(csc.rowval == jcsc.rowval)
   
    assertArrayEqualEps(csc.nzval,jcsc.nzval)
   
    # @show d.eval_jac_g_timer
    # @show d.jd.eval_jac_g_timer
    return
end

function MathProgBase.eval_hesslag_prod(
    d::TapeNLPEvaluator,
    h::Vector{Float64}, # output vector
    x::Vector{Float64}, # current solution
    v::Vector{Float64}, # rhs vector
    σ::Float64,         # multiplier for objective
    μ::Vector{Float64}) # multipliers for each constraint

    assert(false)    
    MathProgBase.eval_hesslag_prod(d,h,x,v,σ,μ)
end

function MathProgBase.hesslag_structure(d::TapeNLPEvaluator)
    # @show "MathProgBase.hesslag_structure" 
    if(d.laghess_nnz != -1)
        return d.laghess_I, d.laghess_J        
    end
   
    hess_structure2(d.obj_tt)
    I = d.obj_tt.h_I
    J = d.obj_tt.h_J

    for i=1:d.numConstr
        hess_structure2(d.constr_tt[i])
        append!(I,d.constr_tt[i].h_I)
        append!(J,d.constr_tt[i].h_J)
    end
    
    #
    # hess_structure2(d.laghess_tt)
    # I = d.laghess_tt.h_I
    # J = d.laghess_tt.h_J

    assert(length(I) == length(J))
    V = ones(Float64,length(I))
    csc = sparse(I,J,V,d.numVar,d.numVar)
    # @show I
    # @show J


    (jI, jJ) = MathProgBase.hesslag_structure(d.jd)
    assert(length(jI) == length(jJ))
    jV = ones(Float64,length(jI))
    jcsc = sparse(jI,jJ,jV,d.numVar,d.numVar)
    # @show jI
    # @show jJ
   

    # @show csc
    # @show jcsc
    # @show csc.colptr
    # @show jcsc.colptr
    assert(csc.m == jcsc.m)
    assert(csc.n == jcsc.n)
    # assert(csc.colptr == jcsc.colptr)
    # assert(csc.rowval == jcsc.rowval)

    
    d.laghess_nnz = length(I)
    d.laghess_I = I
    d.laghess_J = J
    return  d.laghess_I, d.laghess_J
end

function MathProgBase.eval_hesslag(
    d::TapeNLPEvaluator,
    H::Vector{Float64},         # Sparse hessian entry vector
    x::Vector{Float64},         # Current solution
    obj_factor::Float64,        # Lagrangian multiplier for objective
    lambda::Vector{Float64})    # Multipliers for each constraint
    

    @assert length(lambda) == d.numConstr
    @assert length(H) == length(d.laghess_J)
    @assert length(H) == length(d.laghess_I)
    # @assert length(d.p) == d.numConstr + 1
    
    prepare_reeval_hess2(d.obj_tt)
    for i=1:d.numConstr
        prepare_reeval_hess2(d.constr_tt[i])
    end

    # prepare_reeval_hess2(d.laghess_tt)
    ## clean up done

    # d.p[1] = obj_factor
    # for i = 2:length(d.p)
    #     @inbounds d.p[i] = lambda[i-1]
    # end

    tic()
    hess_reverse2(d.obj_tt,x,d.pvals,obj_factor)
    m=1
    for i=1:length(d.obj_tt.hess)
        @inbounds H[m] = d.obj_tt.hess[i]
        m += 1
    end
    
    for i=1:d.numConstr
        @inbounds tt = d.constr_tt[i]
        @inbounds hess_reverse2(tt,x,d.pvals,lambda[i])
        @inbounds for j=1:length(tt.hess) 
            @inbounds H[m] = tt.hess[j]
            m += 1
        end
    end

    # hess_reverse2(d.laghess_tt,x,d.p)
    d.eval_hesslag_timer += toq()
    # H = d.laghess_tt.hess

    csc = sparse(d.laghess_I,d.laghess_J,H)
    # @show csc

    tic()
    jH = Array{Float64,1}(length(d.jd.hess_I))
    MathProgBase.eval_hesslag(d.jd,jH,x,obj_factor,lambda)
    d.jeval_hesslag_timer += toq()
    
    # @show d.jd.hess_I
    # @show d.jd.hess_J
    # @show jH
    jcsc = sparse(d.jd.hess_I,d.jd.hess_J,jH)
    # @show jcsc
    # @show csc.colptr, jcsc.colptr
    # @show csc.rowval,jcsc.rowval
    # assert(csc.colptr == jcsc.colptr)
    # assert(csc.rowval == jcsc.rowval)
    # @show csc.m,jcsc.m
    # @show csc.n,jcsc.n
    # assert(csc.m == jcsc.m)
    # assert(csc.n == jcsc.n)
    
    # @show d.eval_hesslag_timer
    # @show d.jd.eval_hesslag_timer
    return
end

function assertArrayEqualEps(a1,a2)
    assert(length(a1)==length(a2))
    for i=1:length(a1)
        v1 = a1[i]
        v2 = a2[i]
        assert(abs(v1-v2)<EPS)
    end
end

end

#########

