
#forward evaluation for a scalar function
function forward_evaluate{V,I}(tt::Array{I,1}, vvals::Array{V,1}, pvals::Array{V,1})
	idx = one(I)
	vals = Array{V,1}()
	v = [zero(V)]

	sizehint!(v,1)
	sizehint!(vals,round(I,length(tt)/10)) #just a random guess
	
	while(idx <= length(tt))
		# @show idx
		@inbounds ntype = tt[idx]
		idx += 1
		v[1] = zero(V)
		if(ntype == TYPE_P)
			@inbounds v[1] = pvals[tt[idx]]
			idx += 1
			@inbounds push!(vals,v[1])
			idx += 1 #skip TYPE_P
		elseif(ntype == TYPE_V)
			@inbounds v[1] = vvals[tt[idx]]
			idx += 1
			@inbounds push!(vals,v[1])
			idx += 1 #skip TYPE_V
		elseif(ntype == TYPE_O)
			@inbounds oc = tt[idx]
			idx += 1
			@inbounds n = tt[idx]
			idx += 1
			idx += 1 #skip TYPE_O
			# @show OP[oc],vals
			@inbounds evaluate(OP[oc],vals,length(vals)-n+1,v)
			@inbounds resize!(vals,length(vals)-n+1)
			# @show v
			@inbounds vals[end] = v[1]
		end
		# @show vals
	end
	return vals[1]
end

## Interface method
function feval{V,I}(tt::Array{I,1}, vvals::Array{V,1}, pvals::Array{V,1})
	val = forward_evaluate(tt,vvals,pvals)
	return val
end
