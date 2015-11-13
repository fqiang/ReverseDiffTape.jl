
#forward evaluation for a scalar function
function forward_pass_0ord{I,V}(tape::Tape{I}, vvals::Array{V,1}, pvals::Array{V,1})
	tt = tape.tt
	idx = one(I)
	
	stk = Array{V,1}()
	# sizehint!(stk, tape.maxoperands+20)
	stk = MyArray{V}(tape.maxoperands+20) #to use MyArray
	v = [zero(V)]
	sizehint!(v,1)
	
	@inbounds while(idx <= length(tt))
		# @show idx
		ntype = tt[idx]
		idx += 1
		v[1] = zero(V)
		# @show ntype
		if(ntype == TYPE_P)
			@inbounds v[1] = pvals[tt[idx]]
			idx += 1
			# @show v[1]
			@inbounds push!(stk,v[1])
			idx += 1 #skip TYPE_P
		elseif(ntype == TYPE_V)
			@inbounds v[1] = vvals[tt[idx]]
			idx += 1
			# @show v[1]
			@inbounds push!(stk,v[1])
			idx += 1 #skip TYPE_V
		elseif(ntype == TYPE_O)
			@inbounds oc = tt[idx]
			idx += 1
			@inbounds n = tt[idx]
			idx += 1
			idx += 1 #skip TYPE_O
			# @show OP[oc],stk
			# @inbounds eval_0ord(OP[oc],stk,length(stk)-n+1,length(stk),v) 
			@inbounds eval_0ord(OP[oc],stk.a,length(stk)-n+1,length(stk),v)  #using MyArray
			@inbounds resize!(stk,length(stk)-n+1)
			# @show v
			@inbounds stk[end] = v[1]
		end
	end
	return stk[1]
end

## Interface method
function feval{I,V}(tape::Tape{I}, vvals::Array{V,1}, pvals::Array{V,1})
	val = forward_pass_0ord(tape,vvals,pvals)
	return val
end
