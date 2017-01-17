#temp2.jl
include("test_small.jl") 

using ReverseDiffTape

jd = JuMP.NLPEvaluator(m)
MathProgBase.initialize(jd,[:ExprGraph])
ex=MathProgBase.obj_expr(jd);

stk = Vector{Float64}();
vals = Vector{Float64}();
imm = Vector{Float64}();

tt = Tape{Int,Float64}(stk,vals,imm);
p = Vector{Float64}();
tapeBuilder(ex,tt,p,MathProgBase.numvar(jd.m));
resizeWorkingMemory(tt,Tape{Int,Float64}(stk,vals,imm),Vector{Tape{Int,Float64}}())

x=Vector{Float64}()
for i=1:tt.nvar
	push!(x,1.0*i)
end

@time feval(tt,x,p)

gI = Vector{Int}(tt.nvnode);
@time grad_structure(tt,1,gI)
g = zeros(length(gI));
@time  grad_reverse(tt,x,p,1,g)

g_dense = zeros(MathProgBase.numvar(jd.m));
@time grad_reverse_dense(tt,x,p,g_dense)



@assert sum(abs(sparsevec(gI, g)-  g_dense)) == 0.0



h_I = Vector{Int}();
h_J = Vector{Int}();
@time hess_structure(tt,h_I,h_J)
@assert length(h_I) == length(h_J)

h = zeros(length(h_I));
@time hess_reverse(tt,x,p,1,h)


ReverseDiffTape.reset_hess(tt)
resize!(h_I,0);
resize!(h_J,0);
fill!(h,0.0);
Profile.clear_malloc_data()
@time hess_structure(tt,h_I,h_J)
@time hess_reverse(tt,x,p,1,h)




