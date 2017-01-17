using ReverseDiffTape 
using Base.Test
using Calculus

# write your own tests here
using FactCheck
import Base.count

global y1 = 1.1
global y2 = 2.2
global y3 = 3.3
global y4 = 4.4


function to_dense_array(sparse)
	dense = Vector{Float64}(length(sparse))
	fill!(dense,0.0)
	for i=1:length(sparse)
		dense[i] = sparse[i]
	end
	return dense
end

function SparseMatrix.sparse(I,J,V, M, N;keepzeros=false)
    if(!keepzeros)
        return sparse(I,J,V,M,N)
    else
        full = sparse(I,J,ones(Float64,length(I)),M,N)
        actual = sparse(I,J,V,M,N)
        fill!(full.nzval,0.0)

        for c = 1:N
            for i=nzrange(actual,c)
                r = actual.rowval[i]
                v = actual.nzval[i]
                if(v!=0)
                    full[r,c] = v
                end
            end  
            # full.nzval[crange] = actual.nzval[crange] 
        end        
        return full
    end
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
	
	grad_reverse(tt,x,p)
	ag = to_dense_array(sparsevec(tt.g_I, tt.g,tt.nvar))


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
	grad_reverse(tt,x,p)
	ag = to_dense_array(sparsevec(tt.g_I, tt.g,tt.nvar))

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
	grad_reverse(tt,x,p)
	ag = to_dense_array(sparsevec(tt.g_I, tt.g,tt.nvar))

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
	grad_reverse(tt,x,p)
	ag = to_dense_array(sparsevec(tt.g_I, tt.g,tt.nvar))


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
	grad_reverse(tt,x,p)
	ag = to_dense_array(sparsevec(tt.g_I, tt.g,tt.nvar))


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
	grad_reverse(tt,x,p)
	ag = to_dense_array(sparsevec(tt.g_I, tt.g,tt.nvar))


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

include("runtest_ep.jl")
### Version 2 - reimplement Hessian EP
# include("runtest_ep2.jl")


FactCheck.exitstatus()
