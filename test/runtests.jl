using ReverseDiffTape 
using Base.Test

# write your own tests here
using FactCheck

## Test for forward function evaluation
facts("Function evaluataion ") do
	tt = TT_TYPE()
	pvals = TV_TYPE()
	vvals = TV_TYPE()
	x1 = AD_V(tt, vvals, 1.1)
	x2 = AD_V(tt, vvals, 2.2)
	x3 = AD_V(tt, vvals, 3.3)
	p1 = AD_P(tt, pvals, 1)
	p2 = AD_P(tt, pvals, 2)
	c = sin(x1)+cos(x2^p2) * p1 - x3*p2
	val = feval(tt,vvals,pvals)
	@fact sin(1.1)+cos(2.2^2)*1-3.3*2 --> val
end

facts("Reverse gradient sin(a)*cos(b) ") do
	tt = TT_TYPE()
	pvals = TV_TYPE()
	vvals = TV_TYPE()
	a = AD_V(tt,vvals,1.1)
	b = AD_V(tt,vvals,2.2)
	c=sin(a)*cos(b)
	grad = Array{VV_TYPE,1}(length(vvals))
	grad_reverse(tt,vvals,pvals,grad)
	I = Array{IDX_TYPE,1}()
	grad_structure(tt,I)
	@fact length(I) --> 2
	@fact length(grad) --> length(vvals)
	@fact grad[1] --> cos(2.2)*cos(1.1)
	@fact grad[2] --> sin(1.1)*(-sin(2.2))
end

facts("Reverse gradient a^3") do
	tt = TT_TYPE()
	pvals = TV_TYPE()
	vvals = TV_TYPE()
	a = AD_V(tt,vvals,1.1)
	b = AD_P(tt,pvals,3)
	c = a^b
	grad = Array{VV_TYPE,1}(length(vvals))
	grad_reverse(tt,vvals,pvals,grad)
	@fact length(grad) --> length(vvals)
	@fact grad[1] --> 3*1.1^2
end

facts("Reverse gradient x1*x2*x3*x4") do
	tt = TT_TYPE()
	pvals = TV_TYPE()
	vvals = TV_TYPE()
	x1 = AD_V(tt,vvals,1.1)
	x2 = AD_V(tt,vvals,2.2)
	x3 = AD_V(tt,vvals,3.3)
	x4 = AD_V(tt,vvals,4.4)
	c = x1*x2*x3*x4
	grad = Array{VV_TYPE,1}(length(vvals))
	grad_reverse(tt,vvals,pvals,grad)
	@fact length(grad) --> length(vvals)
	@fact grad[1] --> 2.2*3.3*4.4
	@fact grad[2] --> 1.1*3.3*4.4
	@fact grad[3] --> 1.1*2.2*4.4
	@fact grad[4] --> 1.1*2.2*3.3	
end

facts("Reverse gradient sin(x1*x2*x3*x4)") do
	tt = TT_TYPE()
	pvals = TV_TYPE()
	vvals = TV_TYPE()
	x1 = AD_V(tt,vvals,1.1)
	x2 = AD_V(tt,vvals,2.2)
	x3 = AD_V(tt,vvals,3.3)
	x4 = AD_V(tt,vvals,4.4)
	c = sin(x1*x2*x3*x4)
	grad = Array{VV_TYPE,1}(length(vvals))
	grad_reverse(tt,vvals,pvals,grad)
	@fact length(grad) --> length(vvals)
	@fact grad[1] --> cos(1.1*2.2*3.3*4.4)*2.2*3.3*4.4
	@fact grad[2] --> cos(1.1*2.2*3.3*4.4)*1.1*3.3*4.4
	@fact grad[3] --> cos(1.1*2.2*3.3*4.4)*1.1*2.2*4.4
	@fact grad[4] --> cos(1.1*2.2*3.3*4.4)*1.1*2.2*3.3	
end


facts("Reverse gradient cos(x1*x2*x3*x4)") do
	tt = TT_TYPE()
	pvals = TV_TYPE()
	vvals = TV_TYPE()
	x1 = AD_V(tt,vvals,1.1)
	x2 = AD_V(tt,vvals,2.2)
	x3 = AD_V(tt,vvals,3.3)
	x4 = AD_V(tt,vvals,4.4)
	c = cos(x1*x2*x3*x4)
	grad = Array{VV_TYPE,1}(length(vvals))
	grad_reverse(tt,vvals,pvals,grad)
	@fact length(grad) --> length(vvals)
	@fact grad[1] --> -sin(1.1*2.2*3.3*4.4)*2.2*3.3*4.4
	@fact grad[2] --> -sin(1.1*2.2*3.3*4.4)*1.1*3.3*4.4
	@fact grad[3] --> -sin(1.1*2.2*3.3*4.4)*1.1*2.2*4.4
	@fact grad[4] --> -sin(1.1*2.2*3.3*4.4)*1.1*2.2*3.3	
end



facts("Hessian EP algorithm x1^2") do
	tt = TT_TYPE()
	pvals = TV_TYPE()
	vvals = TV_TYPE()
	x1 = AD_V(tt,vvals,1.1)
	p2 = AD_P(tt,pvals,2)
	c = x1^p2
	eset = EdgeSet()
	reverse_hess_ep(tt,vvals,pvals,eset)
	h11 = Edge(tt,1,1)
	@fact length(eset) --> 1
	@fact haskey(eset,h11) --> true
	@fact eset[h11] --> 2
end


facts("Hessian EP algorithm x1^2*x2^2") do
	tt = TT_TYPE()
	pvals = TV_TYPE()
	vvals = TV_TYPE()
	x1 = AD_V(tt,vvals,1.1)
	x2 = AD_V(tt,vvals,2.2)
	p2 = AD_P(tt,pvals,2)
	c = x1^p2*x2^p2
	eset = EdgeSet()
	reverse_hess_ep(tt,vvals,pvals,eset)
	h11 = Edge(tt,1,1)
	h12 = Edge(tt,1,2)
	h21 = Edge(tt,2,1)
	h22 = Edge(tt,2,2)
	@fact length(eset) --> 3
	@fact haskey(eset,h11) --> true
	@fact haskey(eset,h12) --> true
	@fact haskey(eset,h21) --> true
	@fact haskey(eset,h22) --> true
	@fact eset[h11] --> 2.2*2.2*2
	@fact eset[h12] --> 2*1.1*2*2.2
	@fact eset[h21] --> 2*1.1*2*2.2
	@fact eset[h22] --> 1.1*1.1*2 
end

facts("Hessian EP algorithm sin(x1)") do
	tt = TT_TYPE()
	pvals = TV_TYPE()
	vvals = TV_TYPE()
	x1 = AD_V(tt,vvals,1.1)
	c = sin(x1)
	eset = EdgeSet()
	reverse_hess_ep(tt,vvals,pvals,eset)
	h11 = Edge(tt,1,1)
	@fact length(eset) --> 1
	@fact haskey(eset,h11) --> true
	@fact eset[h11] --> -sin(1.1) 
end


facts("Hessian EP algorithm cos(x1)") do
	tt = TT_TYPE()
	pvals = TV_TYPE()
	vvals = TV_TYPE()
	x1 = AD_V(tt,vvals,1.1)
	c = cos(x1)
	eset = EdgeSet()
	reverse_hess_ep(tt,vvals,pvals,eset)
	h11 = Edge(tt,1,1)
	@fact length(eset) --> 1
	@fact haskey(eset,h11) --> true
	@fact eset[h11] --> -cos(1.1) 
end

facts("Hessian EP algorithm x1*x2") do
	tt = TT_TYPE()
	pvals = TV_TYPE()
	vvals = TV_TYPE()
	x1 = AD_V(tt,vvals,1.1)
	x2 = AD_V(tt,vvals,2.2)
	c = x1*x2
	eset = EdgeSet()
	reverse_hess_ep(tt,vvals,pvals,eset)
	h12 = Edge(tt,1,2)
	@fact length(eset) --> 1
	@fact eset[h12] --> 1
end


facts("Hessian EP algorithm sin(x1*x2)") do
	tt = TT_TYPE()
	pvals = TV_TYPE()
	vvals = TV_TYPE()
	x1 = AD_V(tt,vvals,1.1)
	x2 = AD_V(tt,vvals,2.2)
	c = sin(x1*x2)
	eset = EdgeSet()
	reverse_hess_ep(tt,vvals,pvals,eset)
	h11 = Edge(tt,1,1)
	h12 = Edge(tt,1,2)
	h21 = Edge(tt,2,1)
	h22 = Edge(tt,2,2)
	@fact length(eset) --> 3
	@fact eset[h12] --> cos(1.1*2.2) - 1.1*2.2*sin(1.1*2.2)
	@fact eset[h11] --> -2.2^2*sin(1.1*2.2)
	@fact eset[h21] --> cos(1.1*2.2) - 1.1*2.2*sin(1.1*2.2)
	@fact eset[h22] --> -1.1^2*sin(1.1*2.2)
end

facts("Hessian EP algorithm cos(sin(x1))") do
	tt = TT_TYPE()
	pvals = TV_TYPE()
	vvals = TV_TYPE()
	x1 = AD_V(tt,vvals,1.1)
	c = cos(sin(x1))
	eset = EdgeSet()
	reverse_hess_ep(tt,vvals,pvals,eset)
	h11 = Edge(tt,1,1)
	@fact length(eset) --> 1
	@fact eset[h11] --> sin(sin(1.1))*sin(1.1) - cos(sin(1.1))*(cos(1.1)^2)
end


facts("Hessian EP algorithm cos(sin(x1))") do
	tt = TT_TYPE()
	pvals = TV_TYPE()
	vvals = TV_TYPE()
	x1 = AD_V(tt,vvals,1.1)
	x2 = AD_V(tt,vvals,2.2)
	c = cos(sin(x1*x2))
	eset = EdgeSet()
	reverse_hess_ep(tt,vvals,pvals,eset)
	h11 = Edge(tt,1,1)
	@fact length(eset) --> 3
	@fact eset[h11] --> roughly(2.2^2*sin(sin(1.1*2.2))*sin(1.1*2.2)-2.2^2*cos(sin(1.1*2.2))*cos(1.1*2.2)^2)
end

facts("Hessian EP algorithm x1*x1") do
	tt = TT_TYPE()
	pvals = TV_TYPE()
	vvals = TV_TYPE()
	x1 = AD_V(tt,vvals,1.1)
	c = x1*x1
	eset = EdgeSet()
	reverse_hess_ep(tt,vvals,pvals,eset)
	h11 = Edge(tt,1,1)
	@fact length(eset) --> 1
	@fact eset[h11] --> 2
end

facts("Hessian EP algorithm cos(x1*x2)") do
	tt = TT_TYPE()
	pvals = TV_TYPE()
	vvals = TV_TYPE()
	x1 = AD_V(tt,vvals,1.1)
	x2 = AD_V(tt,vvals,2.2)
	c = cos(x1*x2)
	eset = EdgeSet()
	reverse_hess_ep(tt,vvals,pvals,eset)
	h11 = Edge(tt,1,1)
	h12 = Edge(tt,1,2)
	h21 = Edge(tt,2,1)
	h22 = Edge(tt,2,2)
	@fact length(eset) --> 3
	@fact haskey(eset,h11) --> true
	@fact haskey(eset,h12) --> true
	@fact haskey(eset,h21) --> true
	@fact haskey(eset,h22) --> true
	@fact eset[h11] --> -cos(1.1*2.2)*2.2*2.2
	@fact eset[h12] --> -(sin(1.1*2.2)+cos(1.1*2.2)*1.1*2.2)
	@fact eset[h21] --> -(sin(1.1*2.2)+cos(1.1*2.2)*1.1*2.2)
	@fact eset[h22] --> -cos(1.1*2.2)*1.1*1.1 
end


facts("Hessian EP algorithm x1*x1*x2*x2") do
	tt = TT_TYPE()
	pvals = TV_TYPE()
	vvals = TV_TYPE()
	x1 = AD_V(tt,vvals,1.1)
	x2 = AD_V(tt,vvals,2.2)
	p2 = AD_P(tt,pvals,2)
	c = x1*x1*x2*x2
	eset = EdgeSet()
	reverse_hess_ep(tt,vvals,pvals,eset)
	h11 = Edge(tt,1,1)
	h12 = Edge(tt,1,2)
	h21 = Edge(tt,2,1)
	h22 = Edge(tt,2,2)
	@fact length(eset) --> 3
	@fact haskey(eset,h11) --> true
	@fact haskey(eset,h12) --> true
	@fact haskey(eset,h21) --> true
	@fact haskey(eset,h22) --> true
	@fact eset[h11] --> 2.2*2.2*2
	@fact eset[h12] --> 2*1.1*2*2.2
	@fact eset[h21] --> 2*1.1*2*2.2
	@fact eset[h22] --> 1.1*1.1*2 
end

facts("Hessian EP algorithm x1*x1*x1") do
	tt = TT_TYPE()
	pvals = TV_TYPE()
	vvals = TV_TYPE()
	x1 = AD_V(tt,vvals,1.1)
	p2 = AD_P(tt,pvals,2)
	c = x1*x1*x1
	eset = EdgeSet()
	reverse_hess_ep(tt,vvals,pvals,eset)
	h11 = Edge(tt,1,1)
	@fact length(eset) --> 1
	@fact haskey(eset,h11) --> true
	@fact eset[h11] --> 6*1.1
end

facts("Hessian EP algorithm sin(x1)+cos(x2^2)*1-x3*2") do
	tt = TT_TYPE()
	pvals = TV_TYPE()
	vvals = TV_TYPE()
	x1 = AD_V(tt, vvals, 1.1)
	x2 = AD_V(tt, vvals, 2.2)
	x3 = AD_V(tt, vvals, 3.3)
	p1 = AD_P(tt, pvals, 1)
	p2 = AD_P(tt, pvals, 2)
	c = sin(x1)+cos(x2^p2) * p1 - x3*p2
	eset = EdgeSet()
	reverse_hess_ep(tt,vvals,pvals,eset)
	h11 = Edge(tt,1,1)
	h22 = Edge(tt,2,2)
	nzeset = Set{Edge}()
	hess_structure(tt,nzeset)
	@fact length(nzeset) --> 2
	@fact in(h11,nzeset)  --> true
	@fact in(h22,nzeset)  --> true
	@fact length(eset) --> 2
	@fact eset[h11] --> -sin(1.1)
	@fact eset[h22] --> -2*sin(2.2^2)-4*2.2^2*cos(2.2^2)
	# print(eset)
end
FactCheck.exitstatus()
