
#forward evaluation for a scalar function
function forward_evaluate{V,I}(tape::Tape{I}, vvals::Array{V,1}, pvals::Array{V,1})
	tt = tape.tt
	idx = one(I)
	vals = Array{V,1}()
	v = [zero(V)]

	sizehint!(v,1)
	sizehint!(vals,tape.maxoperands)
	
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
function feval{V,I}(tape::Tape{I}, vvals::Array{V,1}, pvals::Array{V,1})
	val = forward_evaluate(tape,vvals,pvals)
	return val
end
