
## forward evaluation 
function evaluate(tt::TT_TYPE,idx::IDX_TYPE, vvals::TV_TYPE, pvals::TV_TYPE)
	# debug("enter - ",idx)
	assert(idx != 0)
	ret = NaN
	ntype = tt[idx] #type node
	idx -= 1
	if(ntype == TYPE_P)
		ret = pvals[tt[idx]]
		idx -=1
	elseif(ntype == TYPE_V)
		ret = vvals[tt[idx]]
		idx -= 1
	elseif(ntype == TYPE_OB)
		oc = tt[idx]
		assert(B_OP_START<= oc <= B_OP_END)
		idx -= 1
		ridx = tt[idx]
		idx -= 1
		lidx = tt[idx]
		idx -= 1
		# debug("before right - ",ridx)
		(rval,idx) = evaluate(tt,ridx,vvals,pvals)
		# debug("before left - ",lidx)
		(lval,idx) = evaluate(tt,lidx,vvals,pvals)
		ret = evaluate(OP[oc],lval,rval)
	else
		oc = tt[idx]
		assert(U_OP_START <= oc <= U_OP_END)
		idx -= 1
		lidx = tt[idx]
		idx -= 1
		# debug("unary before  - ",lidx)
		(lval,idx) = evaluate(tt,lidx,vvals,pvals)
		ret = evaluate(OP[oc],lval)
	end
	assert(!isnan(ret))
	# debug("exit - ",idx)
	return ret,idx
end


## Interface method
function feval(tt::TT_TYPE, vvals::TV_TYPE, pvals::TV_TYPE)
	(val,idx) = evaluate(tt,convert(IDX_TYPE,length(tt)),vvals, pvals)
	return val
end
