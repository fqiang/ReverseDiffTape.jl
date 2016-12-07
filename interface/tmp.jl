#temp2.jl
include("test_small.jl") 

using ReverseDiffTape

jd = JuMP.NLPEvaluator(m)
MathProgBase.initialize(jd,[:ExprGraph])
ex=MathProgBase.obj_expr(jd);

tt = Tape{Int,Float64}();
p = Vector{Float64}();
tapeBuilder(ex,tt,p);

x=Vector{Float64}()
for i=1:tt.nvar
	push!(x,1.0*i)
end
@time feval(tt,x,p)

@time grad_structure(tt)

@time  grad_reverse(tt,x,p)

@time hess_structure(tt)

@time hess_reverse(tt,x,p)

ReverseDiffTape.prepare_reeval_hess(tt)
@time hess_reverse(tt,x,p)


ReverseDiffTape.reset_hess(tt)
Profile.clear_malloc_data()
@time hess_structure(tt)
@time hess_reverse(tt,x,p)




