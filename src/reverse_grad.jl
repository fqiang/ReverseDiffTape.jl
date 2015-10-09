
# forward pass on the tape tt, to build ss stack
function reverse_grad_0(tt::TT_TYPE, idx::IDX_TYPE,ss::TV_STACK,vvals::TV_TYPE,pvals::TV_TYPE)
	# debug("enter - ",idx)
	assert(idx != 0)
	val = NaN
	ntype = tt[idx] #type node
	idx -= 1
	if(ntype == TYPE_P)
		val = pvals[tt[idx]]
		idx -=1
		# push!(ss, val)
	elseif(ntype == TYPE_V)
		vidx = tt[idx]
		idx -= 1
		val = vvals[vidx]
		# push!(ss,val)
	elseif(ntype == TYPE_OB)
		oc = tt[idx]
		assert(B_OP_START<= oc <= B_OP_END)
		idx -= 1
		ridx = tt[idx]
		idx -= 1
		lidx = tt[idx]
		idx -= 1
		# debug("before right - ",ridx)
		rval = reverse_grad_0(tt,ridx,ss,vvals, pvals)
		# debug("before left - ",lidx)
		lval = reverse_grad_0(tt,lidx,ss,vvals, pvals)
		(val,ld,rd) = reverse_grad_calc(oc,lval,rval)
		# push!(ss,val)
		push!(ss,rd)
		push!(ss,ld)
	elseif(ntype == TYPE_OU)
		oc = tt[idx]
		assert(U_OP_START <= oc <= U_OP_END)
		idx -= 1
		lidx = tt[idx]
		idx -= 1
		# debug("unary before  - ",lidx)
		lval = reverse_grad_0(tt,lidx,ss,vvals, pvals)
		(val,ld) = reverse_grad_calc(oc,lval)
		# push!(ss,val)
		push!(ss,ld)
	else
		assert(false)
	end
	assert(!isnan(val))
	# debug("exit - ",idx)
	return val
end

function reverse_grad_calc(oc::OP_TYPE,lval::VV_TYPE, rval::VV_TYPE)
	ld = NaN
	rd = NaN
	val = evaluate(OP[oc],lval,rval)
	if(OP[oc]==:+)
		ld = 1
		rd = 1
	elseif(OP[oc]==:-)
		ld = 1
		rd = -1
	elseif(OP[oc]==:*)
		ld = rval
		rd = lval
	elseif(OP[oc]==:/)
		ld = 1/rval
		rd = -lval/rval^2
	elseif(OP[oc]==:^)
		ld = rval*(lval^(rval-1))
		rd = val*(log(lval))
	else
		assert(false)
	end
	assert(ld!=NaN && rd!=NaN)
	return (val,ld,rd)
end

function reverse_grad_calc(oc::OP_TYPE,lval::VV_TYPE)
	ld = NaN
	val = evaluate(OP[oc],lval)
	if(OP[oc]==:sin)
		ld = cos(lval)
	elseif(OP[oc]==:cos)
		ld = -sin(lval)
	else
		assert(false)
	end
	assert(ld != NaN)
	return (val,ld)
end

function reverse_grad_1(tt::TT_TYPE, idx::IDX_TYPE, ss::TV_STACK, adj::VV_TYPE, vvals::TV_TYPE, pvals::TV_TYPE, grad::TV_TYPE)
	assert(length(vvals) == length(grad))
	# debug("enter - ", adj)
	ntype = tt[idx]
	idx -= 1
	if(ntype ==TYPE_P)  
		#do nothing , this is a parameter
	elseif(ntype == TYPE_V)
		vidx = tt[idx]
		idx -= 1
		# debug("vidx ",vidx, "  ",adj)
		grad[vidx] += adj
	elseif(ntype == TYPE_OB)
		dl = pop!(ss)
		dr = pop!(ss)
		idx -= 1 #skip oc
		ridx = tt[idx]
		idx -= 1
		lidx = tt[idx]
		idx -= 1
		l_adj = adj*dl
		r_adj = adj*dr
		# debug("before left ",lidx, " ", adj," ", dl," ", l_adj)
		reverse_grad_1(tt,lidx,ss,l_adj,vvals,pvals,grad)
		# debug("before right ",ridx," ", adj," ", dr," ", r_adj)
		reverse_grad_1(tt,ridx,ss,r_adj,vvals,pvals,grad)
	elseif(ntype == TYPE_OU)
		dl = pop!(ss)
		idx -= 1 #skip oc
		lidx = tt[idx]
		idx -= 1
		l_adj = adj*dl
		reverse_grad_1(tt,lidx,ss,l_adj,vvals,pvals,grad)
	else
		assert(false)
	end
	# debug("exit")
end

function grad_nz(tt::TT_TYPE, idx::IDX_TYPE, nz::UInt)
	assert(length(vvals) == length(grad))
	# debug("enter - ", adj)
	ntype = tt[idx]
	idx -= 1
	if(ntype ==TYPE_P)  
		#do nothing , this is a parameter
	elseif(ntype == TYPE_V)
		vidx = tt[idx]
		nz += 1
	elseif(ntype == TYPE_OB)
		idx -= 1
		ridx = tt[idx]
		idx -= 1
		lidx = tt[idx]
		idx -= 1
		nzg(tt,lidx,nz)
		nzg(tt,ridx,nz)
	elseif(ntype == TYPE_OU)
		idx -= 1 #skip oc
		lidx = tt[idx]
		idx -= 1
		nzg(tt,lidx,nz)
	else
		assert(false)
	end
end

#Interface function
function nzg(tt::TT_TYPE)
	nz = 0
	nzg(tt,convert(IDX_TYPE,length(tt)),nz)
	return nz
end

function grad_reverse(tt::TT_TYPE,vvals::TV_TYPE,pvals::TV_TYPE)
	ss = TV_STACK(VV_TYPE)
	reverse_grad_0(tt,convert(IDX_TYPE,length(tt)),ss,vvals,pvals)
	adj = 1.0
	grad = zeros(Float64,length(vvals))
	reverse_grad_1(tt,convert(IDX_TYPE,length(tt)),ss,adj,vvals,pvals,grad)
	return grad
end
