#edge pusing algorithm for Hessian reverse AD
immutable Edge
	tt::TT_TYPE
    lidx::IDX_TYPE
    ridx::IDX_TYPE
end

function Base.isequal(e1::Edge, e2::Edge)
	# println("isequal, ", e1, ":",e2)
	assert(e1.tt == e2.tt)
	return (e1.lidx == e2.lidx && e1.ridx == e2.ridx) || (e1.lidx == e2.ridx && e1.ridx == e2.lidx)
end

function Base.hash(e::Edge)
	v = ""
	if e.lidx > e.ridx
		v = string(e.lidx,e.ridx)
	else
		v = string(e.ridx,e.lidx)
	end
	return hash(v)
end

typealias EdgeSet Dict{Edge,VV_TYPE}


function isEndPointOnEdge(e::Edge, idx::IDX_TYPE)
	return (e.lidx == idx ) || (e.ridx == idx  )
end

function isSelfEdge(e::Edge)
	return e.lidx == e.ridx 
end

function removeEdges(eset::EdgeSet, idx::IDX_TYPE)
	edges = keys(eset)
	state = start(edges)
	while !done(edges,state)
		(e, state) = next(edges, state)
		if (isEndPointOnEdge(e,idx))
			delete!(eset,e)
		end
	end
end


function removeEdges(eset::Set{Edge}, idx::IDX_TYPE)
	state = start(eset)
	while !done(eset,state)
		(e, state) = next(eset, state)
		if (isEndPointOnEdge(e,idx))
			delete!(eset,e)
		end
	end
end

function isBothVariable(e::Edge)
	lidx = e.lidx
	ridx = e.ridx
	return e.tt[lidx] == TYPE_V && e.tt[ridx] == TYPE_V
end

function Base.show(io::IO, e::Edge)
	print(io,"Edge (",e.lidx, ",", e.ridx,")")
end
function Base.show(io::IO, eset::EdgeSet)
	state = start(eset)
	while !done(eset,state)
		(pair,state) = next(eset,state)
		println(io,pair[1],"=",pair[2])
	end
end

function createEdge(eset::EdgeSet, tt::TT_TYPE, lidx::IDX_TYPE,ridx::IDX_TYPE, w::VV_TYPE=0.0)
	e = Edge(tt,lidx,ridx)
	# println("create edge: ",e)
	if(!haskey(eset,e))
		# println("create edge: ",e)
		eset[e] = w
	else
		# println("already exist: ",e)
	end
	return e
end

function createEdge(eset::Set{Edge}, tt::TT_TYPE, lidx::IDX_TYPE,ridx::IDX_TYPE)
	e = Edge(tt,lidx,ridx)
	# println("create edge: ",e)
	push!(eset,e)
	return e
end

function reverse_hess_ep0(tt::TT_TYPE, idx::IDX_TYPE, ss::TV_STACK, vvals::TV_TYPE, pvals::TV_TYPE)
	assert(idx != 0)
	val = NaN
	ntype = tt[idx]
	idx -= 1
	if(ntype == TYPE_P)
		val = pvals[tt[idx]]
		idx -= 1
	elseif(ntype == TYPE_V)
		val = vvals[tt[idx]]
		idx -= 1
	elseif(ntype == TYPE_OB)
		oc = tt[idx]
		assert(B_OP_START<= oc <= B_OP_END)
		idx -= 1
		ridx = tt[idx]
		idx -= 1
		lidx = tt[idx]
		idx -= 1
		rval = reverse_hess_ep0(tt,ridx,ss,vvals,pvals)
		lval = reverse_hess_ep0(tt,lidx,ss,vvals,pvals)
		(val,hl,hr,hll,hrr,hlr) = reverse_hess_calc(oc, lval, rval)
		# push!(ss,adj)
		push!(ss,hl)
		push!(ss,hr)
		push!(ss,hll)
		push!(ss,hrr)
		push!(ss,hlr)
	elseif(ntype == TYPE_OU)
		oc = tt[idx]
		assert(U_OP_START<= oc <= U_OP_END)
		idx -= 1
		lidx = tt[idx]
		idx -= 1
		lval = reverse_hess_ep0(tt,lidx,ss,vvals,pvals)
		(val,hl,hll) = reverse_hess_calc(oc,lval)
		# push!(ss,adj)
		push!(ss,hl)
		push!(ss,hll)
	else
		assert(false)
	end
	assert(!isnan(val))
	return val
end

function reverse_hess_calc(oc::OP_TYPE, lval::Real, rval::Real)
	val;hl;hll;hrr;hlr = NaN
	val = evaluate(OP[oc],lval,rval)
	if(OP[oc]==:+)
		hl = 1
		hr = 1
		hll = 0
		hrr = 0
		hlr = 0
	elseif(OP[oc]==:-)
		hl = 1
		hr = -1
		hll = 0
		hrr = 0
		hlr = 0
	elseif(OP[oc]==:*)
		hl = rval
		hr = lval
		hll = 0
		hrr = 0
		hlr = 1
	elseif(OP[oc]==:/)
		hl = 1/rval
		hr = -lval/(rval^2)
		hll = 0
		hrr = 2*lval/(rval^3)
		hlr = -1/(rval^2)
	elseif(OP[oc]==:^)
		hl = rval*(lval^(rval-1))
		hr = 0
		hll = rval*(rval-1)*(lval^(rval-2))
		hrr = 0
		hlr = 0
	else
		assert(false)
	end
	assert(val!=NaN && hl != NaN && hll!=NaN && hrr!=NaN &&hlr!=NaN)
	return (val,hl,hr,hll,hrr,hlr)
end


function reverse_hess_calc(oc::OP_TYPE, lval::Real)
	val;hl;hll = NaN
	val = evaluate(OP[oc],lval)
	if(OP[oc]==:sin)
		hl = cos(lval)
		hll = -sin(lval)
	elseif(OP[oc]==:cos)
		hl = -sin(lval)
		hll = -cos(lval)
	else
		assert(false)	
	end
	assert(val!=NaN && hl!=NaN && hll!=NaN)
	return (val,hl,hll)
end

function hess_nz(tt::TT_TYPE, idx::IDX_TYPE, eset::Set{Edge})
	ntype = tt[idx]
	this_idx = idx
	idx -= 1
	if(ntype == TYPE_P)
		removeEdges(eset,this_idx)
	elseif(ntype == TYPE_V)
		#nothing
	elseif(ntype == TYPE_OB)
		oc = tt[idx]
		idx -= 1 #skip oc
		ridx = tt[idx]
		idx -= 1
		lidx = tt[idx]
		idx -= 1

		state = start(eset)
		while !done(eset,state)
			(e,state) = next(eset,state)
			if(isEndPointOnEdge(e,this_idx))
				if(isSelfEdge(e))
					e1 = createEdge(eset,tt,lidx,ridx)
					e2 = createEdge(eset,tt,lidx,lidx)
					e3 = createEdge(eset,tt,ridx,ridx)
				else
					oidx = e.lidx == this_idx?e.ridx:e.lidx
					e1 = createEdge(eset,tt,lidx,oidx)
					e2 = createEdge(eset,tt,ridx,oidx)
				end
				delete!(eset,e)
			end
		end

		#creating
		if(OP[oc] == :+ || OP[oc] == :-)
			#do nothing for linear operator
		elseif(OP[oc] == :*)
			e1 = createEdge(eset,tt,lidx,ridx)
		elseif(OP[oc] == :/)
			e1 = createEdge(eset,tt,lidx,ridx)
			e3 = createEdge(eset,tt,ridx,ridx)
		elseif(OP[oc] == :^)
			e1 = createEdge(eset,tt,lidx,ridx)
			e2 = createEdge(eset,tt,lidx,lidx)
			e3 = createEdge(eset,tt,ridx,ridx)
		else
			assert(false);
		end

		hess_nz(tt,lidx,eset)
		hess_nz(tt,ridx,eset)
	elseif(ntype == TYPE_OU)
		oc = tt[idx]
		idx -= 1 #skip oc
		lidx = tt[idx]

		#pushing
		state = start(eset)
		while !done(eset,state)
			(e,state) = next(eset,state)
			if(isEndPointOnEdge(e,this_idx))
				if(isSelfEdge(e))
					#case 2
					e2 = createEdge(eset,tt,lidx,lidx)
				else
					#case 1, 3
					oidx = e.lidx == this_idx? e.ridx:e.lidx
					e1 = createEdge(eset,tt,lidx,oidx)
				end
				delete!(eset, e)
			end
		end

		#creating
		if(OP[oc]==:sin)
			e1 = createEdge(eset,tt,lidx,lidx)
		elseif (OP[oc] == :cos)
			e1 = createEdge(eset,tt,lidx,lidx)
		else
			assert(false)
		end

		hess_nz(tt,lidx,eset)
	else
		assert(false)
	end
end

function hess_nz1(eset::Set{Edge})
	veset = Set{Edge}()
	state = start(eset)
	while !done(eset,state)
		(e,state) = next(eset,state)
		lidx = e.lidx
		ridx = e.ridx
		assert(e.tt[lidx] == TYPE_V)
		assert(e.tt[ridx] == TYPE_V)
		lidx -= 1
		ridx -= 1
		lvidx = e.tt[lidx]
		rvidx = e.tt[ridx]
		ve = Edge(e.tt,lvidx,rvidx)
		push!(veset,ve)
	end
	return veset
end

function reverse_hess_ep1(tt::TT_TYPE, idx::IDX_TYPE, ss::TV_STACK, adj::VV_TYPE, eset::EdgeSet)
	# debug("enter - ", adj)
	this_idx = idx
	ntype = tt[idx]
	idx -= 1
	if(ntype ==TYPE_P)  
		#deleting all edges that has an endpoint to this node
		removeEdges(eset,this_idx)
	elseif(ntype == TYPE_V)
		#do nothing
	elseif(ntype == TYPE_OB)
		hlr = pop!(ss)
		hrr = pop!(ss)
		hll = pop!(ss)
		hr = pop!(ss)
		hl = pop!(ss)

		oc = tt[idx]
		assert(B_OP_START<= oc <= B_OP_END)
		idx -= 1
		ridx = tt[idx]
		idx -= 1
		# rt = tt[ridx]
		lidx = tt[idx]
		idx -= 1
		# lt = tt[lidx]
		
		#pushing
		state = start(eset)
		# println("binary op before pushing - eset(",length(eset),")")
		while !done(eset,state)
			(pair, state) = next(eset,state)
			e = pair[1]
			w = pair[2]
			if(isEndPointOnEdge(e,this_idx))
				if(isSelfEdge(e))
					#case 2
					e1 = createEdge(eset,tt,lidx,ridx)
					e2 = createEdge(eset,tt,lidx,lidx)
					e3 = createEdge(eset,tt,ridx,ridx)
					eset[e1] += (hl*hr*w);
					eset[e2] += (hl*hl*w);
					eset[e3] += (hr*hr*w);
				else
					#case 1 , 3
					oidx = e.lidx == this_idx? e.ridx: e.lidx
					e1 = createEdge(eset,tt,lidx,oidx)
					e2 = createEdge(eset,tt,ridx,oidx)
					e1.lidx == e1.ridx? eset[e1] += (2*hl*w): eset[e1] += (hl*w)
					e2.lidx == e2.ridx? eset[e2] += (2*hr*w): eset[e2] += (hr*w)
				end
				delete!(eset,e)
			end
		end
		# println("binary op after pushing - eset(",length(eset),")")
		# println("before remove set rmv(",length(rmv),") - eset(",length(eset),")")
		# setdiff!(eset, rmv)
		# println("after remove set rmv(",length(rmv),") - eset(",length(eset),")")

		#creating
		if(OP[oc] == :+ || OP[oc] == :-)
			#do nothing for linear operator
		elseif(OP[oc] == :*)
			e1 = createEdge(eset,tt,lidx,ridx)
			# e2 = createEdge(eset,tt,ridx,lidx)
			eset[e1] += (adj*hlr)
			# eset[e2] += (adj*hlr)
		elseif(OP[oc] == :/)
			e1 = createEdge(eset,tt,lidx,ridx)
			e3 = createEdge(eset,tt,ridx,ridx)
			eset[e1] += (adj*hlr)
			eset[e3] += (adj*hrr)
		elseif(OP[oc] == :^)
			e1 = createEdge(eset,tt,lidx,ridx)
			e2 = createEdge(eset,tt,lidx,lidx)
			e3 = createEdge(eset,tt,ridx,ridx)
			eset[e1] += (adj*hlr)
			eset[e2] += (adj*hll)
			eset[e3] += (adj*hrr)
		else
			assert(false);
		end
		
		#adjoint
		ladj = adj*hl;
		radj = adj*hr;

		#recursive
		reverse_hess_ep1(tt, lidx, ss, ladj, eset)
		reverse_hess_ep1(tt, ridx, ss, radj, eset)
	elseif(ntype == TYPE_OU)
		hll = pop!(ss)
		hl = pop!(ss)
		oc = tt[idx]
		idx -= 1
		assert(U_OP_START<= oc <= U_OP_END)
		lidx = tt[idx]
		idx -= 1

		#pushing
		state = start(eset)
		while !done(eset,state)
			(pair, state) = next(eset,state)
			e = pair[1]
			w = pair[2]
			if(isEndPointOnEdge(e,this_idx))
				if(isSelfEdge(e))
					#case 2
					e2 = createEdge(eset,tt,lidx,lidx)
					eset[e2] += hl*hl*w
				else
					#case 1, 3
					oidx = e.lidx == this_idx? e.ridx:e.lidx
					e1 = createEdge(eset,tt,lidx,oidx)
					e1.lidx == e1.ridx? eset[e1] += 2*hl*w: eset[e1] += hl*w
				end
				delete!(eset, e)
			end
		end

		#creating
		if(OP[oc]==:sin)
			e1 = createEdge(eset,tt,lidx,lidx)
			eset[e1] += (adj*hll)
		elseif (OP[oc] == :cos)
			e1 = createEdge(eset,tt,lidx,lidx)
			eset[e1] += (adj*hll)
		else
			assert(false)
		end

		#adjoint
		ladj = adj*hl;
		
		#recursive
		reverse_hess_ep1(tt,lidx,ss,ladj,eset)
	else
		assert(false)
	end
	# debug("exit")
end


function reverse_hess_ep2(eset::EdgeSet)
	veset = EdgeSet()
	state = start(eset)
	while !done(eset,state)
		(pair, state) = next(eset, state)
		e = pair[1]
		w = pair[2]
		lidx = e.lidx
		ridx = e.ridx
		assert(e.tt[lidx] == TYPE_V)
		assert(e.tt[ridx] == TYPE_V)
		lidx -= 1
		ridx -= 1
		lvidx = e.tt[lidx]
		rvidx = e.tt[ridx]
		if(lidx != ridx && lvidx == rvidx) 
			w += w
		end
		# println("reverse_hess_ep2 - adding (",lvidx,",",rvidx,") = ",w)
		ve = Edge(e.tt,lvidx,rvidx)
		if(!haskey(veset,ve))
			veset[ve] = w
		else
			veset[ve] += w
		end
	end
	return veset
end

#Interface function
function nzh(tt::TT_TYPE)
	eset = Set{Edge}()
	hess_nz(tt, convert(IDX_TYPE,length(tt)), eset)
	veset = hess_nz1(eset)
end


function reverse_hess_ep(tt::TT_TYPE,vvals::TV_TYPE,pvals::TV_TYPE)
	ss = TV_STACK(VV_TYPE)
	reverse_hess_ep0(tt,convert(IDX_TYPE,length(tt)),ss,vvals,pvals)
	adj = 1.0
	eset = EdgeSet()
	reverse_hess_ep1(tt,convert(IDX_TYPE,length(tt)),ss,adj,eset)
	assert(isempty(ss) == true)
	veset = reverse_hess_ep2(eset)
	return veset
end