
module TapeInterface

using ReverseDiffTape
import MathProgBase

const EPS=1e-5


export  TapeSolver, 
        reset_timer, report_mem

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

type TapeNLPEvaluator <: MathProgBase.AbstractNLPEvaluator
    jd::MathProgBase.AbstractNLPEvaluator
    numVar::Int
    numConstr::Int
    pvals::Vector{Float64}
    obj_tt::Tape{Int,Float64}
    constr_tt::Vector{Tape{Int,Float64}}
    nl_idxes::Vector{Int}
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
    init_timer::Float64
    eval_f_timer::Float64
    eval_g_timer::Float64
    eval_grad_f_timer::Float64
    eval_jac_g_timer::Float64
    eval_hesslag_timer::Float64

    jinit_timer::Float64
    jeval_f_timer::Float64
    jeval_g_timer::Float64
    jeval_grad_f_timer::Float64
    jeval_jac_g_timer::Float64
    jeval_hesslag_timer::Float64

    f_times::Int
    g_times::Int
    grad_f_times::Int
    jac_g_times::Int
    hesslag_times::Int

    #jump's temp storage
    jfval::Float64
    jgrad_f::Vector{Float64}
    jg::Vector{Float64}
    jjac_I::Vector{Int}
    jjac_J::Vector{Int}
    jjac_g::Vector{Float64}
    jlaghess_I::Vector{Int}
    jlaghess_J::Vector{Int}
    jlaghess_h::Vector{Float64}

end


function TapeNLPEvaluator(nlpe::MathProgBase.AbstractNLPEvaluator,numVar,numConstr)
    return TapeNLPEvaluator(nlpe, 
        numVar,numConstr,
        Vector{Float64}(), #parameter values
        
        Tape{Int,Float64}(), #objective tape
        Vector{Tape{Int,Float64}}(),  #constraints tape

        Vector{Int}(), #nl_idxes
        -1, #init

        # Tape{Int,Float64}(),  #full expr
        # Vector{Float64}(),    #obj factor + lambda

        Vector{Int}(), Vector{Int}(),-1, #Jacobian
        
        Vector{Int}(), Vector{Int}(),-1, #Hessian
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0,   #my timer
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0,   #jump timer
        0,0,0,0,0,                      #counter
        0.0, 
        Vector{Float64}(), 
        Vector{Float64}(),
        Vector{Int}(),
        Vector{Int}(),
        Vector{Float64}(),
        Vector{Int}(),
        Vector{Int}(),
        Vector{Float64}()
        )
end


function MathProgBase.NonlinearModel(solver::TapeSolver)
	@show "TapeSolver model"
	return TapeData(nothing,MathProgBase.NonlinearModel(solver.s))
end 

function MathProgBase.setwarmstart!(m::TapeData, x) 
    return MathProgBase.setwarmstart!(m.m,x)
end

function MathProgBase.status(m::TapeData)
    e = m.evaluator
    @show e.init_timer
    @show e.eval_f_timer
    @show e.eval_g_timer
    @show e.eval_grad_f_timer
    @show e.eval_jac_g_timer
    @show e.eval_hesslag_timer
    
    @show e.jinit_timer
    @show e.jeval_f_timer
    @show e.jeval_g_timer
    @show e.jeval_grad_f_timer
    @show e.jeval_jac_g_timer
    @show e.jeval_hesslag_timer

    @show e.f_times
    @show e.g_times
    @show e.grad_f_times
    @show e.jac_g_times
    @show e.hesslag_times

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
    @show "TapeData - loadproblem"
    tape_evaluator = TapeNLPEvaluator(jd,numVar,numConstr)
    m.evaluator = tape_evaluator
    MathProgBase.loadproblem!(m.m,numVar,numConstr,l,u,lb,ub,sense,tape_evaluator)
end

function MathProgBase.optimize!(m::TapeData)
    MathProgBase.optimize!(m.m)
end

function MathProgBase.initialize(d::TapeNLPEvaluator, requested_features::Vector{Symbol})
	@show "TapeNLPEvaluator - initialize"

    tic()
    jd = d.jd
    MathProgBase.initialize(jd,[:Grad,:Jac,:ExprGraph,:Hess])
    jd_init = toq()
    @show "JuMP evaluator initialize ", jd_init

	#let's building up the tape
    if(d.init == -1)
        tic()
    	objexpr = MathProgBase.obj_expr(jd)
        # @show objexpr
       
    	assert(length(d.pvals) == 0)  ## assume no values in pvals
    	tapeBuilder(objexpr,d.obj_tt,d.pvals)
        
        # hesslag_expr = AD(d.obj_tt.tt)*AD_P(d.p,1.0)

        for i =1:d.numConstr
            conexpr = MathProgBase.constr_expr(jd,i)
            if !MathProgBase.isconstrlinear(jd,i) 
                push!(d.nl_idxes,i)
            end
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

        d.init_timer = toq()
        @show "TapeNLPEvaluator - tape build ",d.init_timer

        #initilize JuMP temp storage
        resize!(d.jgrad_f,d.numVar)
        resize!(d.jg,d.numConstr)
        (d.jjac_I,d.jjac_J) = MathProgBase.jac_structure(d.jd)
        resize!(d.jjac_g,length(d.jjac_I))

        (d.jlaghess_I, d.jlaghess_J) = MathProgBase.hesslag_structure(d.jd)
        resize!(d.jlaghess_h,length(d.jlaghess_J))

        d.init = 1  #set initialized tapes
    end

    # reset timers
    reset_timer(d)

    nothing
end

MathProgBase.features_available(d::TapeNLPEvaluator) = [:Grad, :Jac, :Hess, :HessVec]

#evaluate the objective function given on iterate x
function MathProgBase.eval_f(d::TapeNLPEvaluator, x)
    #jump
    tic()
    d.jfval = MathProgBase.eval_f(d.jd,x)
    d.jeval_f_timer += toq()
    
    #tape
    tic()
    v = feval(d.obj_tt,x,d.pvals)
    d.eval_f_timer += toq()
    
    # assert(d.jfval == v)
    d.f_times += 1
    return v
end

#evaluate the objective gradient on given iterate x. Results is set to g
function MathProgBase.eval_grad_f(d::TapeNLPEvaluator, g, x)
    assert(length(g) == length(x))
    fill!(g,0.0)
	assert(sum(g) == 0.0)

    #JuMP    
    # jg = Array{Float64,1}(length(x))
    tic()
    MathProgBase.eval_grad_f(d.jd,d.jgrad_f,x)
    d.jeval_grad_f_timer += toq()
    
    #tape
    tape = d.obj_tt
    if tape.nzg==-1
        grad_structure(tape)
    end

    tic()
    grad_reverse(tape,x,d.pvals)
    d.eval_grad_f_timer += toq()

    for i=1:length(tape.g_I)
        @inbounds g[tape.g_I[i]] = tape.g[i]
    end
   
    assertArrayEqualEps(g,d.jgrad_f)
    d.grad_f_times += 1
    return
end

#constraint evaluation
function MathProgBase.eval_g(d::TapeNLPEvaluator, g, x)
    assert(length(g) == d.numConstr)

    #jump
    # jg = Array{Float64,1}(length(g))
    tic()
    MathProgBase.eval_g(d.jd,d.jg,x)
    d.jeval_g_timer += toq()

    #tape
    tic()
    @inbounds for i=1:d.numConstr
        @inbounds g[i]=feval(d.constr_tt[i],x,d.pvals)
    end
    d.eval_g_timer += toq()
    
    # @show g,jg
    assertArrayEqualEps(g,jg)
    d.g_times += 1
    return
end


function MathProgBase.jac_structure(d::TapeNLPEvaluator)
    #jump
    # jV = ones(Float64,length(d.jjac_I))
    # jcsc = sparse(d.jjac_I,d.jjac_J,jV)


    #tape
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
    
    assert(length(d.jac_I) == length(d.jac_J))
    # V = ones(Float64,length(d.jac_I))
    # csc = sparse(d.jac_I,d.jac_J,V)

    d.jac_nnz = length(d.jac_I)
    return d.jac_I, d.jac_J
end

function MathProgBase.eval_jac_g(d::TapeNLPEvaluator, J, x)
    assert(d.jac_nnz != -1)  #structure already computed
    assert(length(J) == d.jac_nnz)
    
    #jump
    tic()
    MathProgBase.eval_jac_g(d.jd,d.jjac_g,x)
    d.jeval_jac_g_timer += toq()

    #tape
    tic()
    J_len = 0
    @inbounds for i = 1:d.numConstr
        tape_i = d.constr_tt[i]
        grad_reverse(tape_i,x,d.pvals)
        append_array(J,J_len,tape_i.g,0,length(tape_i.g))
        J_len += length(tape_i.g)
    end
    d.eval_jac_g_timer += toq()
    
    # csc = sparse(d.jac_I,d.jac_J,J)
    # jcsc = sparse(d.jd.jac_I,d.jd.jac_J,jJ)
    # matrixEqaul(csc,jcsc)

    # @show d.eval_jac_g_timer
    # @show d.jd.eval_jac_g_timer
    d.jac_g_times += 1
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
    if(d.laghess_nnz != -1)
        return d.laghess_I, d.laghess_J        
    end
    
    tic()
    hess_structure2(d.obj_tt)
    I = d.obj_tt.h_I
    J = d.obj_tt.h_J
    
    for i=1:length(d.nl_idxes)
        @inbounds tt = d.constr_tt[d.nl_idxes[i]]
        hess_structure2(tt)
        append!(I,tt.h_I)
        append!(J,tt.h_J)
    end
    hesslag_time = toq()
    @show hesslag_time
    
    # hess_structure2(d.laghess_tt)
    # I = d.laghess_tt.h_I
    # J = d.laghess_tt.h_J

    # assert(length(I) == length(J))
    # V = ones(Float64,length(I))
    # csc = sparse(I,J,V,d.numVar,d.numVar)
    

    # jV = ones(Float64,length(d.jlaghess_I))
    # jcsc = sparse(d.jlaghess_I,d.jlaghess_J,jV,d.numVar,d.numVar)
    
    
    d.laghess_nnz = length(I)
    d.laghess_I = I
    d.laghess_J = J

    # @show "start drawing"
    # writedlm("hess_structure.txt", csc)
    # @show "end drawing"
    @show "Hessian nnz list", length(d.laghess_I)

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
    @assert length(d.jlaghess_h) == length(d.jlaghess_I)
    @assert length(d.jlaghess_h) == length(d.jlaghess_J)

    #jump
    tic()    
    MathProgBase.eval_hesslag(d.jd,d.jlaghess_h,x,obj_factor,lambda)
    d.jeval_hesslag_timer += toq()
    
    #tape
    # @assert length(d.p) == d.numConstr + 1
    prepare_reeval_hess2(d.obj_tt)
    for i=1:length(d.nl_idxes)
        @inbounds prepare_reeval_hess2(d.constr_tt[d.nl_idxes[i]])
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
    
    for i=1:length(d.nl_idxes)
        @inbounds tt = d.constr_tt[d.nl_idxes[i]]
        @inbounds hess_reverse2(tt,x,d.pvals,lambda[i])
        for j=1:length(tt.hess) 
            @inbounds H[m] = tt.hess[j]
            m += 1
        end
    end
    d.eval_hesslag_timer += toq()
    # hess_reverse2(d.laghess_tt,x,d.p)
    # H = d.laghess_tt.hess
    # d.eval_hesslag_timer += toq()

    # csc = sparse(d.laghess_I,d.laghess_J,H)
    # jcsc = sparse(d.jd.hess_I,d.jd.hess_J,jH)
   
    # matrixEqaul(csc,jcsc)
    # @show d.eval_hesslag_timer
    # @show d.jd.eval_hesslag_timer
    d.hesslag_times += 1
    return
end

function assertArrayEqualEps(a1,a2)
    @assert length(a1)==length(a2)
    for i=1:length(a1)
        v1 = a1[i]
        v2 = a2[i]
        if (isnan(v1) && isnan(v2)) || (isinf(v1) && isinf(v2))
            return 
        else
            assert(abs(v1-v2)<EPS)
        end
    end
end

function matrixEqaul(mat1,mat2)
    @show "checking matrix equal"
    (n1,m1) = size(mat1)
    (n2,m2) = size(mat2)
    assert(n1 == n2)
    assert(m2 == m2)
    for i=1:n1
        for j=1:m1
            @assert abs(mat1[i,j]-mat2[i,j])<EPS  "$(i), $(j)"
        end
    end
end

function reset_timer(d)
    d.init_timer = 0.0
    d.eval_f_timer = 0.0
    d.eval_grad_f_timer = 0.0
    d.eval_g_timer = 0.0
    d.eval_jac_g_timer = 0.0
    d.eval_hesslag_timer = 0.0

    d.jinit_timer = 0.0
    d.jeval_f_timer = 0.0
    d.jeval_grad_f_timer = 0.0
    d.jeval_g_timer = 0.0
    d.jeval_jac_g_timer = 0.0
    d.jeval_hesslag_timer = 0.0

    d.f_times = 0
    d.grad_f_times = 0
    d.g_times = 0
    d.jac_g_times = 0
    d.hesslag_times = 0
end


function report_mem(e::TapeNLPEvaluator)
    nb = 0::Int
    nb += sizeof(e)

    nb += sizeof(e.pvals)
    nb += report_tape_mem(e.obj_tt)
    for t in e.constr_tt
        nb += report_tape_mem(t)
    end
    nb += sizeof(e.nl_idxes)
    nb += sizeof(e.jac_I)
    nb += sizeof(e.jac_J)
    nb += sizeof(e.laghess_I)
    nb += sizeof(e.laghess_J)
    @show "only tape - ",nb/1024/1024

    # jfval::Float64
    nb += sizeof(e.jgrad_f)
    nb += sizeof(e.jg)
    nb += sizeof(e.jjac_I)
    nb += sizeof(e.jjac_J)
    nb += sizeof(e.jjac_g)
    nb += sizeof(e.jlaghess_I)
    nb += sizeof(e.jlaghess_J)
    nb += sizeof(e.jlaghess_h)
    @show "with jump temp storage", nb/1024/1024

    return nb
end


end #end module
#########

