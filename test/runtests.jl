using ReverseDiffTape 
using Base.Test
using Calculus

# write your own tests here
using FactCheck

global y1 = 1.1
global y2 = 2.2
global y3 = 3.3
global y4 = 4.4

function count(eset)
	c = 0
	for (i,iset) in eset
		for j in iset
			c+=1
		end
	end
	return c
end

## Test for forward function evaluation
facts("Function evaluataion ") do
	p = Vector{Float64}()
	x = Vector{Float64}()
	x1 = AD_V(x, 1.1)
	x2 = AD_V(x, 2.2)
	x3 = AD_V(x, 3.3)
	p1 = AD_P(p, 1)
	p2 = AD_P(p, 2)
	c = sin(x1)+cos(x2^p2) * p1 - x3*p2
	tt = tapeBuilder(c.data)
	val = feval(tt,x,p)
	@fact sin(1.1)+cos(2.2^2)*1-3.3*2 --> val
end

facts("Reverse gradient") do
	x = Vector{Float64}()
	p = Vector{Float64}()
	x1 = AD_V(x,1.1)
	x2 = AD_V(x,2.2)
	x3 = AD_V(x,3.3)
	x4 = AD_V(x,4.4)
	p1 = AD_P(p,1.0)
	p2 = AD_P(p,2.0)
	p3 = AD_P(p,3.0)

	#############################################	
	c=sin(x1)*cos(x2)
	tt = tapeBuilder(c.data)
	grad_structure(tt)
	ag = Vector{Float64}(tt.nvar)
	grad_reverse(tt,x,p,ag)

	@fact tt.nvar --> 2
	@fact length(tt.g_I) --> tt.nvnode
	@fact length(tt.g) --> length(tt.g_I)
	g = Vector{Float64}(tt.nvar)
	for i=1:length(g)
		g[i] = @eval $(differentiate("sin(y1)*cos(y2)",[:y1,:y2])[i])
	end
	for i in 1:tt.nvar
		@fact g[i] --> ag[i]
	end

	#############################################	
	c = x1^p3
	tt = tapeBuilder(c.data)
	grad_structure(tt)
	ag = Vector{Float64}(tt.nvar)
	grad_reverse(tt,x,p,ag)

	@fact tt.nvar --> 1
	@fact length(tt.g_I) --> tt.nvnode
	@fact length(tt.g) --> length(tt.g_I)
	g = Vector{Float64}(tt.nvar)
	for i=1:length(g)
		g[i] = @eval $(differentiate("y1^3",[:y1])[i])
	end
	for i in 1:tt.nvar
		@fact g[i] --> ag[i]
	end


	#############################################	
	c = x1^p2
	tt = tapeBuilder(c.data)
	grad_structure(tt)
	ag = Vector{Float64}(tt.nvar)
	grad_reverse(tt,x,p,ag)

	@fact tt.nvar --> 1
	@fact length(tt.g_I) --> tt.nvnode
	@fact length(tt.g) --> length(tt.g_I)
	g = Vector{Float64}(tt.nvar)
	for i=1:length(g)
		g[i] = @eval $(differentiate("y1^2",[:y1])[i])
	end
	for i in 1:tt.nvar
		@fact g[i] --> ag[i]
	end

	#############################################	
	c = x1*x2*x3*x4
	tt = tapeBuilder(c.data)
	grad_structure(tt)
	ag = Vector{Float64}(tt.nvar)
	grad_reverse(tt,x,p,ag)

	@fact tt.nvar --> 4
	@fact length(tt.g_I) --> tt.nvnode
	@fact length(tt.g) --> length(tt.g_I)
	g = Vector{Float64}(tt.nvar)
	for i=1:length(g)
		g[i] = @eval $(differentiate("y1*y2*y3*y4",[:y1,:y2,:y3,:y4])[i])
	end
	for i in 1:tt.nvar
		@fact g[i] --> ag[i]
	end

	#############################################
	c = sin(x1*x2*x3*x4)
	tt = tapeBuilder(c.data)
	grad_structure(tt)
	ag = Vector{Float64}(tt.nvar)
	grad_reverse(tt,x,p,ag)

	@fact tt.nvar --> 4
	@fact length(tt.g_I) --> tt.nvnode
	@fact length(tt.g) --> length(tt.g_I)
	g = Vector{Float64}(tt.nvar)
	for i=1:length(g)
		g[i] = @eval $(differentiate("sin(y1*y2*y3*y4)",[:y1,:y2,:y3,:y4])[i])
	end
	for i in 1:tt.nvar
		@fact g[i] --> ag[i]
	end
	#############################################
	c = cos(x1*x2*x3*x4)
	tt = tapeBuilder(c.data)
	grad_structure(tt)
	ag = Vector{Float64}(tt.nvar)
	grad_reverse(tt,x,p,ag)

	@fact tt.nvar --> 4
	@fact length(tt.g_I) --> tt.nvnode
	@fact length(tt.g) --> length(tt.g_I)
	g = Vector{Float64}(tt.nvar)
	for i=1:length(g)
		g[i] = @eval $(differentiate("cos(y1*y2*y3*y4)",[:y1,:y2,:y3,:y4])[i])
	end
	for i in 1:tt.nvar
		@fact g[i] --> ag[i]
	end
end

facts("Hessian EP algorithm x1^2*x2^2") do
	p = Vector{Float64}()
	x = Vector{Float64}()
	x1 = AD_V(x,1.1)
	x2 = AD_V(x,2.2)
	p2 = AD_P(p,2)
	c = x1^p2*x2^p2
	tt = tapeBuilder(c.data)

	eset = Dict{Int,Set{Int}}()
	hess_structure_lower(tt,eset)
	h = EdgeSet{Int,Float64}()
	hess_reverse(tt,x,p,h)
	
	@fact count(eset) --> 3
	@fact h[1][1] --> 2.2*2.2*2
	@fact h[2][1] --> 2*1.1*2*2.2
	@fact h[2][2] --> 1.1*1.1*2 
end

# facts("Hessian EP algorithm sin(x1)") do
# 	p = Vector{Float64}()
# 	x = Vector{Float64}()
# 	x1 = AD_V(x,1.1)
# 	c = sin(x1)
# 	tt = tapeBuilder(c.data)
# 	eset = Dict{Int,Set{Int}}()
# 	hess_structure_lower(tt,eset)
# 	h = EdgeSet{Int,Float64}()
# 	hess_reverse(tt,x,p,h)
# 	@fact count(eset) --> 1
# 	@fact h[1][1] --> -sin(1.1) 
# end


# facts("Hessian EP algorithm cos(x1)") do
# 	p = Vector{Float64}()
# 	x = Vector{Float64}()
# 	x1 = AD_V(x,1.1)
# 	c = cos(x1)
# 	tt = tapeBuilder(c.data)
# 	eset = Dict{Int,Set{Int}}()
# 	hess_structure_lower(tt,eset)
# 	h = EdgeSet{Int,Float64}()
# 	hess_reverse(tt,x,p,h)
# 	@fact count(eset) --> 1
# 	@fact h[1][1] --> -cos(1.1) 
# end

# facts("Hessian EP algorithm x1*x2") do
# 	p = Vector{Float64}()
# 	x = Vector{Float64}()
# 	x1 = AD_V(x,1.1)
# 	x2 = AD_V(x,2.2)
# 	c = x1*x2
# 	tt = tapeBuilder(c.data)
# 	eset = Dict{Int,Set{Int}}()
# 	hess_structure_lower(tt,eset)
# 	h = EdgeSet{Int,Float64}()
# 	hess_reverse(tt,x,p,h)
# 	@fact count(eset) --> 1
# 	@fact h[2][1] --> 1.0
# end


# facts("Hessian EP algorithm sin(x1*x2)") do
# 	p = Vector{Float64}()
# 	x = Vector{Float64}()
# 	x1 = AD_V(x,1.1)
# 	x2 = AD_V(x,2.2)
# 	c = sin(x1*x2)
# 	tt = tapeBuilder(c.data)
# 	eset = Dict{Int,Set{Int}}()
# 	hess_structure_lower(tt,eset)
# 	h = EdgeSet{Int,Float64}()
# 	hess_reverse(tt,x,p,h)
# 	@fact count(eset) --> 3
# 	@fact h[2][1] --> cos(1.1*2.2) - 1.1*2.2*sin(1.1*2.2)
# 	@fact h[1][1] --> -2.2^2*sin(1.1*2.2)
# 	@fact h[2][2] --> -1.1^2*sin(1.1*2.2)
# end

# facts("Hessian EP algorithm cos(sin(x1))") do
# 	p = Vector{Float64}()
# 	x = Vector{Float64}()
# 	x1 = AD_V(x,1.1)
# 	c = cos(sin(x1))
# 	tt = tapeBuilder(c.data)
# 	eset = Dict{Int,Set{Int}}()
# 	hess_structure_lower(tt,eset)
# 	h = EdgeSet{Int,Float64}()
# 	hess_reverse(tt,x,p,h)
# 	@fact count(eset) --> 1
# 	@fact h[1][1] --> sin(sin(1.1))*sin(1.1) - cos(sin(1.1))*(cos(1.1)^2)
# end


# facts("Hessian EP algorithm cos(sin(x1*x2))") do
# 	p = Vector{Float64}()
# 	x = Vector{Float64}()
# 	x1 = AD_V(x,1.1)
# 	x2 = AD_V(x,2.2)
# 	c = cos(sin(x1*x2))
# 	tt = tapeBuilder(c.data)
# 	eset = Dict{Int,Set{Int}}()
# 	hess_structure_lower(tt,eset)
# 	h = EdgeSet{Int,Float64}()
# 	hess_reverse(tt,x,p,h)
# 	@fact count(eset) --> 3
# 	@fact h[1][1] --> roughly(2.2^2*sin(sin(1.1*2.2))*sin(1.1*2.2)-2.2^2*cos(sin(1.1*2.2))*cos(1.1*2.2)^2)
# end

# facts("Hessian EP algorithm x1*x1") do
# 	p = Vector{Float64}()
# 	x = Vector{Float64}()
# 	x1 = AD_V(x,1.1)
# 	c = x1*x1
# 	tt = tapeBuilder(c.data)
# 	eset = Dict{Int,Set{Int}}()
# 	hess_structure_lower(tt,eset)
# 	h = EdgeSet{Int,Float64}()
# 	hess_reverse(tt,x,p,h)
# 	@fact count(eset) --> 1
# 	@fact h[1][1] --> 2.0
# end

# facts("Hessian EP algorithm cos(x1*x2)") do
# 	p = Vector{Float64}()
# 	x = Vector{Float64}()
# 	x1 = AD_V(x,1.1)
# 	x2 = AD_V(x,2.2)
# 	c = cos(x1*x2)
# 	tt = tapeBuilder(c.data)
# 	eset = Dict{Int,Set{Int}}()
# 	hess_structure_lower(tt,eset)
# 	h = EdgeSet{Int,Float64}()
# 	hess_reverse(tt,x,p,h)
# 	@fact count(eset) --> 3
# 	@fact h[1][1] --> -cos(1.1*2.2)*2.2*2.2
# 	@fact h[2][1] --> -(sin(1.1*2.2)+cos(1.1*2.2)*1.1*2.2)
# 	@fact h[2][2] -->  -cos(1.1*2.2)*1.1*1.1 
# end


# facts("Hessian EP algorithm x1*x1*x2*x2") do
# 	p = Vector{Float64}()
# 	x = Vector{Float64}()
# 	x1 = AD_V(x,1.1)
# 	x2 = AD_V(x,2.2)
# 	p2 = AD_P(p,2)
# 	c = x1*x1*x2*x2
# 	tt = tapeBuilder(c.data)
# 	eset = Dict{Int,Set{Int}}()
# 	hess_structure_lower(tt,eset)
# 	h = EdgeSet{Int,Float64}()
# 	hess_reverse(tt,x,p,h)
# 	@fact count(eset) --> 3
# 	@fact h[1][1] --> 2.2*2.2*2
# 	@fact h[2][1] --> 2*1.1*2*2.2
# 	@fact h[2][2] --> 1.1*1.1*2 
# end

# facts("Hessian EP algorithm x1*x1*x1") do
# 	p = Vector{Float64}()
# 	x = Vector{Float64}()
# 	x1 = AD_V(x,1.1)
# 	p2 = AD_P(p,2)
# 	c = x1*x1*x1
# 	tt = tapeBuilder(c.data)
# 	eset = Dict{Int,Set{Int}}()
# 	hess_structure_lower(tt,eset)
# 	h = EdgeSet{Int,Float64}()
# 	hess_reverse(tt,x,p,h)

# 	@fact count(eset) --> 1
# 	@fact h[1][1] --> 6*1.1
# end

# facts("Hessian EP algorithm sin(x1)+cos(x2^2)*1-x3*2") do
# 	p = Vector{Float64}()
# 	x = Vector{Float64}()
# 	x1 = AD_V(x, 1.1)
# 	x2 = AD_V(x, 2.2)
# 	x3 = AD_V(x, 3.3)
# 	p1 = AD_P(p, 1.0)
# 	p2 = AD_P(p, 2.0)
# 	c = sin(x1)+cos(x2^p2) * p1 - x3*p2
# 	# c = sin(x1)+cos(x2^p2) - x3*p2
# 	tt = tapeBuilder(c.data)
# 	eset = Dict{Int,Set{Int}}()
# 	hess_structure_lower(tt,eset)
# 	h = EdgeSet{Int,Float64}()
# 	hess_reverse(tt,x,p,h)

# 	@fact count(eset) --> 2
# 	@fact h[1][1] --> -sin(1.1)
# 	@fact h[2][2] --> -2*sin(2.2^2)-4*2.2^2*cos(2.2^2)
# end
FactCheck.exitstatus()
