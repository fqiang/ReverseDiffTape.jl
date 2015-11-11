
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


#forward pass on the scalar function
function forward_pass{V,I}(tt::Array{I,1}, vvals::Array{V,1}, pvals::Array{V,1}, vals::Array{V,1})
	idx = 1::I
	v = [zero(V)]
	empty!(vals)
	sizehint!(vals,round(I,length(tt)/3.5))
	sizehint!(v,1)

	@inbounds while(idx <= length(tt))
		# @show idx
		ntype = tt[idx]
		idx += 1
		if(ntype == TYPE_P)
			# tic()
			v[1] = pvals[tt[idx]]
			idx += 1
			push!(vals,v[1])
			idx += 1 #skip TYPE_P
		elseif(ntype == TYPE_V)
			v[1] = vvals[tt[idx]]
			idx += 1
			push!(vals,v[1])
			idx += 1 #skip TYPE_V
		elseif(ntype == TYPE_O)
			oc = tt[idx]
			idx += 1
			n = tt[idx]
			idx += 1
			idx += 1 #skip TYPE_O
			evaluate(OP[oc],vals,length(vals)-n+1,v)
			push!(vals,v[1])
		end
	end
end

function reverse_pass{V,I}(tt::Array{I,1},v::Array{V,1},g::Array{V,1})
	idx = length(tt)
	while(idx > 0)
		ntype = tt[idx]
		idx += 1
		if(ntype == TYPE_P)

		elseif(ntype == TYPE_V)

		elseif(ntype == TYPE_O)

		end
	end
end


# function evaluate{V,I}(s::Symbol, vals::Array{V,1}, i::I, v::Array{V,1})
# 	# n = length(nvals)
# 	# @show s,vals,i
# 	if(s == :+)
# 		@inbounds v[1] = zero(V)
# 		@simd for j=i:1:length(vals)
# 			@inbounds v[1]+=vals[j]
# 		end
# 	elseif(s == :-)
# 		@inbounds return v[1] = (-)(vals[i],vals[i+1])
# 	elseif(s == :*)
# 		@inbounds v[1] = vals[i]
# 		@simd for j=i+1:1:length(vals)
# 			@inbounds v[1] *= vals[j]
# 		end
# 	elseif(s == :/)
# 		@inbounds v[1]=(/)(vals[i],vals[i+1])
# 	elseif(s == :^)
# 		@inbounds v[1]=(^)(vals[i],vals[i+1])
# 	elseif(s == :sin)
# 		@inbounds v[1]=sin(vals[i])
# 	elseif(s == :cos)
# 		@inbounds v[1]=cos(vals[i])
# 	end
# end


## Interface method
function feval{V,I}(tt::Array{I,1}, vvals::Array{V,1}, pvals::Array{V,1})
	val = forward_evaluate(tt,vvals,pvals)
	return val
end
