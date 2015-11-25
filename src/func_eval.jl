
#forward evaluation for a scalar function
function forward_pass_0ord{I,V}(tape::Tape{I,V}, vvals::Array{V,1}, pvals::Array{V,1})
	tt = tape.tt
	idx = one(I)
	
	stk = tape.stk
    stklen = zero(I)
	# sizehint!(stk, tape.maxoperands+20)
	#stk = MyArray{V}(tape.maxoperands+20) #to use MyArray
	
	@inbounds while(idx <= length(tt))
		# @show idx
		ntype = tt[idx]
		idx += 1
		# @show ntype
		if(ntype == TYPE_P)
			@inbounds val = pvals[tt[idx]]
			idx += 1
            stklen += 1
            @inbounds stk[stklen] = val
			idx += 1 #skip TYPE_P
		elseif(ntype == TYPE_V)
			@inbounds val = vvals[tt[idx]]
			idx += 1
            stklen += 1
            @inbounds stk[stklen] = val
			idx += 1 #skip TYPE_V
		elseif(ntype == TYPE_O)
			@inbounds oc = tt[idx]
			idx += 1
			@inbounds n = tt[idx]
			idx += 1
			idx += 1 #skip TYPE_O
            if n == 1 # 1-argument functions
                @inbounds stk[stklen] = eval_0ord(OP[oc],stk[stklen])
            else
                # @show OP[oc],stk
                @inbounds val = eval_0ord(OP[oc],stk,stklen-n+1,stklen)  #using MyArray
                stklen -= n-1
                @inbounds stk[stklen] = val
            end
		end
	end
	return stk[1]
end

## Interface method
function feval{I,V}(tape::Tape{I,V}, vvals::Array{V,1}, pvals::Array{V,1})
	val = forward_pass_0ord(tape,vvals,pvals)
	return val
end
