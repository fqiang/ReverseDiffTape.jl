typealias IDX_TYPE Int
typealias VV_TYPE Float64
typealias OP_TYPE Int
typealias TT_TYPE Array{IDX_TYPE,1}
typealias TV_STACK Stack
typealias TV_TYPE Array{VV_TYPE,1}


#edge pusing algorithm for Hessian reverse AD

typealias EdgeSet{I,V} Dict{I,Dict{I,V}}

function forward_pass_2ord{I,V}(tape::Tape{I}, vvals::Array{V,1}, pvals::Array{V,1}, imm::Array{V,1}, tr::Array{I,1}, eset::EdgeSet{I,V})
	tt = tape.tt
	idx = one(I)
	# empty!(imm)  #used for immediate derrivative
	# empty!(tr)
	stk = Vector{V}(tape.maxoperands+20) #used for value evaluation
	istk = Vector{I}(tape.maxoperands+20)
	stklen = zero(I)
	immlen = zero(I)
	trlen = zero(I)
	
	@inbounds while(idx <= length(tt))
		# @show idx
		# println("++++++++++++++++++++++++++++++++++++")
		ntype = tt[idx]
		# eset[idx] = Dict{I,V}() #initialize edge set
		idx += 1
		if(ntype == TYPE_P)
			stklen += 1
			istk[stklen] = idx-1
			stk[stklen] = pvals[tt[idx]]
			idx += 2 #skip TYPE_P
		elseif(ntype == TYPE_V)
			stklen += 1
			istk[stklen] = idx-1
			stk[stklen] = vvals[tt[idx]]
			idx += 2 #skip TYPE_V
		elseif(ntype == TYPE_O)
			oc = tt[idx]
			idx += 1
			n = tt[idx]
			idx += 2 #skip TYPE_O
			# @show OP[oc], stklen-n+1, n
			# @show stk
			counter = zero(I)
			if(n==1)
				@inbounds (counter,stk[stklen]) = eval_2ord(OP[oc],stk[stklen],imm,immlen+1)
				@inbounds tr[trlen+1] = istk[stklen]
			else
				@inbounds (counter,val) = eval_2ord(OP[oc],stk,stklen-n+1,stklen,imm,immlen+1)
				# @inbounds tr[trlen+1:trlen+n] = @inbounds istk[stklen-n+1:stklen]
				append_array(tr,trlen+1,istk,stklen-n+1,n)
				stklen -= n-1
				@inbounds stk[stklen] = val
			end
			trlen += n 
			immlen += counter
			istk[stklen] = idx - 4
			# @show imm		
			# @show stk[stklen]	
		end
		# @show stklen
		# println("++++++++++++++++++++++++++++++++++++")
	end
	tape.imm2ord = immlen
	assert(tape.nnode-1==trlen)
	# @show stklen
	return stk[1]
end

function make_edge{I,V}(eset::EdgeSet{I,V},i1::I,i2::I,w::V)
	assert(i1>=i2)
	if(haskey(eset[i1],i2))
		eset[i1][i2] += w
	else
		eset[i1][i2] = w
	end
end

@inline function append_array{I,V}(dest::Vector{V},d_offset::I,src::Vector{V},s_offset::I, n::I)
	for i=0:n-1
		@inbounds dest[i+d_offset] = src[i+s_offset]
	end
end


function reverse_pass_2ord{I,V}(tape::Tape{I},imm::Array{V,1},tr::Array{I,1},eset::EdgeSet{I,V})
	assert(tt.nnode-1 == trlen)
	tt = tape.tt
	idx = length(tt)
	trlen = tt.nnode-1
	immlen = tt.imm2ord

	adjs = Array{V,1}()
	sizehint!(adjs, tape.maxoperands+20)  #the acutally size should be (depth_of_tree - current_depth + maxoperands)
	push!(adjs,one(V))  #init value 

	@inbounds while(idx > 0)
		i = idx
		ntype = tt[idx]
		idx -= 1
		adj = pop!(adjs)
		if(ntype == TYPE_P)
			idx -= 2
			delete!(eset,i)
		elseif(ntype == TYPE_V)
			idx -= 2
		elseif(ntype == TYPE_O)
			n = tt[idx]
			idx -= 3 #skip TYPE_O 
			trlen -= n
			#pushing	
			set = eset[i] #Dict{I,V} 
			state = start(set)
			while !done(set,state)
				(p,w) = next(set,state)
				if(i==p)
					for tr_i = trlen-n+1:trlen
						n1 = tr[tr_i]
						h1 = imm[tr_i]
						assert(isempty(eset[n1]))
						make_edge(eset,n1,n1,w*h1*h1)
						for tr_ii=tr_i+1:trlen
							n2 = tr[tr_ii]
							h2 = imm[tr_ii]
							assert(isempty(eset[n2]))
							assert(n2>n1)
							make_edge(eset,n2,n1,w*h1*h2)
						end
					end
				else
					for tr_i = trlen-n+1:trlen
						n1 = tr[tr_i]
						h1 = imm[tr_i]
						if(n1==p) #the circle case
							make_edge(eset,p,p,2*w*h1)
						else
							make_edge(eset,n1,p,w*h1)
						end
					end
				end
			end
			delete!(eset,i)	

			#creating
			if n==1
				make_edge(eset,)
			end
		# 	#creating
		# if(OP[oc] == :+ || OP[oc] == :-)
		# 	#do nothing for linear operator
		# elseif(OP[oc] == :*)
		# 	e1 = createEdge(tt,lidx,ridx)
		# 	w1 = adj*hlr
		# 	aeset[e1] = haskey(aeset,e1)?aeset[e1]+w1:w1
		# 	# eset[e1] += (adj*hlr)
		# elseif(OP[oc] == :/)
		# 	e1 = createEdge(tt,lidx,ridx)
		# 	w1 = adj*hlr
		# 	e3 = createEdge(tt,ridx,ridx)
		# 	w3 = adj*hrr
		# 	aeset[e1] = haskey(aeset,e1)?aeset[e1]+w1:w1
		# 	aeset[e3] = haskey(aeset,e3)?aeset[e3]+w3:w3
		# 	# eset[e1] += (adj*hlr)
		# 	# eset[e3] += (adj*hrr)
		# elseif(OP[oc] == :^)
		# 	e1 = createEdge(tt,lidx,ridx)
		# 	w1 = adj*hlr
		# 	e2 = createEdge(tt,lidx,lidx)
		# 	w2 = adj*hll
		# 	e3 = createEdge(tt,ridx,ridx)
		# 	w3 = adj*hrr
		# 	aeset[e1] = haskey(aeset,e1)?aeset[e1]+w1:w1
		# 	aeset[e2] = haskey(aeset,e2)?aeset[e2]+w2:w2
		# 	aeset[e3] = haskey(aeset,e3)?aeset[e3]+w3:w3
		# 	# eset[e1] += (adj*hlr)
		# 	# eset[e2] += (adj*hll)
		# 	# eset[e3] += (adj*hrr)
		# if(OP[oc]==:sin)
		# 			e1 = createEdge(tt,lidx,lidx)
		# 			w1 = adj*hll
		# 			# eset[e1] += (adj*hll)
		# 			aeset[e1] = haskey(aeset,e1)?aeset[e1]+w1:w1
		# 		elseif (OP[oc] == :cos)
		# 			e1 = createEdge(tt,lidx,lidx)
		# 			w1 = adj*hll
		# 			aeset[e1] = haskey(aeset,e1)?aeset[e1]+w1:w1
		# 			# eset[e1] += (adj*hll)
		# 		else
		# 			assert(false)
		# 		end

			# #adj
			# for i=length(imm)-n+1:length(imm)
			# 	push!(adjs,imm[i]*adj)
			# end
			# resize!(imm,length(imm)-n)
			# resize!(tr,length(tr)-n)
		end
	end
end




immutable Edge{I} 
    lidx::I
    ridx::I

    function Edge(l::I,r::I)
    	return new(l,r)
    end
end

immutable WEdge{I,V}
	e::Edge{I}
	w::V

	function WEdge(e::Edge{I})
		return new(e,zero(V))
	end
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

function Base.isequal(e1::Edge, e2::Edge)
	# assert(e1.tt == e2.tt)  # important ! - can be on different tape if following same context
	return (e1.lidx == e2.lidx && e1.ridx == e2.ridx) || (e1.lidx == e2.ridx && e1.ridx == e2.lidx)
end

function Base.hash(e::Edge)
	if e.lidx > e.ridx
		v = (e.lidx,e.ridx)
	else
		v = (e.ridx,e.lidx)
	end
	return hash(v)
end

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
		
		deset = Set{Edge}()
		aeset = EdgeSet()

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
					e1 = createEdge(tt,lidx,ridx)
					w1 = hl*hr*w
					e2 = createEdge(tt,lidx,lidx)
					w2 = hl*hl*w
					e3 = createEdge(tt,ridx,ridx)
					w3 = hr*hr*w
					aeset[e1] = haskey(aeset,e1)?aeset[e1]+w1:w1
					aeset[e2] = haskey(aeset,e2)?aeset[e2]+w2:w2
					aeset[e3] = haskey(aeset,e3)?aeset[e3]+w3:w3
				else
					#case 1 , 3
					oidx = e.lidx == this_idx? e.ridx: e.lidx
					e1 = createEdge(tt,lidx,oidx)
					e2 = createEdge(tt,ridx,oidx)
					w1 = e1.lidx == e1.ridx? (2*hl*w): (hl*w)
					w2 = e2.lidx == e2.ridx? (2*hr*w): (hr*w)
					aeset[e1] = haskey(aeset,e1)? aeset[e1]+w1:w1
					aeset[e2] = haskey(aeset,e2)? aeset[e2]+w2:w2
				end
				# delete!(eset,e)
				push!(deset,e)
			end
		end

		for k in deset  #delete the pushed edges
			delete!(eset,k)
		end
		# println("binary op after pushing - eset(",length(eset),")")
		# println("before remove set rmv(",length(rmv),") - eset(",length(eset),")")
		# setdiff!(eset, rmv)
		# println("after remove set rmv(",length(rmv),") - eset(",length(eset),")")

		#creating
		if(OP[oc] == :+ || OP[oc] == :-)
			#do nothing for linear operator
		elseif(OP[oc] == :*)
			e1 = createEdge(tt,lidx,ridx)
			w1 = adj*hlr
			aeset[e1] = haskey(aeset,e1)?aeset[e1]+w1:w1
			# eset[e1] += (adj*hlr)
		elseif(OP[oc] == :/)
			e1 = createEdge(tt,lidx,ridx)
			w1 = adj*hlr
			e3 = createEdge(tt,ridx,ridx)
			w3 = adj*hrr
			aeset[e1] = haskey(aeset,e1)?aeset[e1]+w1:w1
			aeset[e3] = haskey(aeset,e3)?aeset[e3]+w3:w3
			# eset[e1] += (adj*hlr)
			# eset[e3] += (adj*hrr)
		elseif(OP[oc] == :^)
			e1 = createEdge(tt,lidx,ridx)
			w1 = adj*hlr
			e2 = createEdge(tt,lidx,lidx)
			w2 = adj*hll
			e3 = createEdge(tt,ridx,ridx)
			w3 = adj*hrr
			aeset[e1] = haskey(aeset,e1)?aeset[e1]+w1:w1
			aeset[e2] = haskey(aeset,e2)?aeset[e2]+w2:w2
			aeset[e3] = haskey(aeset,e3)?aeset[e3]+w3:w3
			# eset[e1] += (adj*hlr)
			# eset[e2] += (adj*hll)
			# eset[e3] += (adj*hrr)
		else
			assert(false);
		end
		
		for (k,v) in aeset  #adding the pushed and created edges
			# println("reverse_hess_ep1 - add, ",k," weight ",v)
			eset[k] = haskey(eset,k)?eset[k]+v:v
			# @show eset
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

		deset = Set{Edge}()
		aeset = EdgeSet()

		#pushing
		state = start(eset)
		while !done(eset,state)
			(pair, state) = next(eset,state)
			e = pair[1]
			w = pair[2]
			if(isEndPointOnEdge(e,this_idx))
				if(isSelfEdge(e))
					#case 2
					e2 = createEdge(tt,lidx,lidx)
					w2 = hl*hl*w
					aeset[e2] = haskey(aeset,e2)?aeset[e2]+w2:w2
					# eset[e2] += hl*hl*w
				else
					#case 1, 3
					oidx = e.lidx == this_idx? e.ridx:e.lidx
					e1 = createEdge(tt,lidx,oidx)
					w1 = e1.lidx == e1.ridx? 2*hl*w: hl*w
					aeset[e1] = haskey(aeset,e1)?aeset[e1]+w1:w1
				end
				# delete!(eset, e)
				push!(deset,e)
			end
		end

		for k in deset  #delete the pushed edges
			delete!(eset,k)
		end

		#creating
		if(OP[oc]==:sin)
			e1 = createEdge(tt,lidx,lidx)
			w1 = adj*hll
			# eset[e1] += (adj*hll)
			aeset[e1] = haskey(aeset,e1)?aeset[e1]+w1:w1
		elseif (OP[oc] == :cos)
			e1 = createEdge(tt,lidx,lidx)
			w1 = adj*hll
			aeset[e1] = haskey(aeset,e1)?aeset[e1]+w1:w1
			# eset[e1] += (adj*hll)
		else
			assert(false)
		end

		for (k,v) in aeset  #adding the pushed and created edges
			eset[k] = haskey(eset,k)?eset[k]+v:v
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


function reverse_hess_ep2(eset::EdgeSet,factor::VV_TYPE, veset::EdgeSet)
	# veset = EdgeSet()
	# @show eset
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
		# println("reverse_hess_ep2 - make edge ",lvidx," -- ",rvidx)

		if(!haskey(veset,ve))
			veset[ve] = w*factor
		else
			veset[ve] += w*factor
		end
	end
	# @show veset
	return veset
end


# function reverse_hess_ep0{I}(tape::Tape{I}, idx::IDX_TYPE, ss::TV_STACK, vvals::TV_TYPE, pvals::TV_TYPE)
# 	assert(idx != 0)
# 	val = NaN
# 	ntype = tt[idx]
# 	idx -= 1
# 	if(ntype == TYPE_P)
# 		val = pvals[tt[idx]]
# 		idx -= 1
# 	elseif(ntype == TYPE_V)
# 		val = vvals[tt[idx]]
# 		idx -= 1
# 	elseif(ntype == TYPE_OB)
# 		oc = tt[idx]
# 		assert(B_OP_START<= oc <= B_OP_END)
# 		idx -= 1
# 		ridx = tt[idx]
# 		idx -= 1
# 		lidx = tt[idx]
# 		idx -= 1
# 		rval = reverse_hess_ep0(tt,ridx,ss,vvals,pvals)
# 		lval = reverse_hess_ep0(tt,lidx,ss,vvals,pvals)
# 		(val,hl,hr,hll,hrr,hlr) = reverse_hess_calc(oc, lval, rval)
# 		# push!(ss,adj)
# 		push!(ss,hl)
# 		push!(ss,hr)
# 		push!(ss,hll)
# 		push!(ss,hrr)
# 		push!(ss,hlr)
# 	elseif(ntype == TYPE_OU)
# 		oc = tt[idx]
# 		assert(U_OP_START<= oc <= U_OP_END)
# 		idx -= 1
# 		lidx = tt[idx]
# 		idx -= 1
# 		lval = reverse_hess_ep0(tt,lidx,ss,vvals,pvals)
# 		(val,hl,hll) = reverse_hess_calc(oc,lval)
# 		# push!(ss,adj)
# 		push!(ss,hl)
# 		push!(ss,hll)
# 	else
# 		assert(false)
# 	end
# 	assert(!isnan(val))
# 	return val
# end

# function reverse_hess_calc(oc::OP_TYPE, lval::Real, rval::Real)
# 	hl = NaN
# 	hll = NaN
# 	hrr = NaN
# 	hlr = NaN
# 	val = evaluate(OP[oc],lval,rval)::VV_TYPE
# 	if(OP[oc]==:+)
# 		hl = 1.0
# 		hr = 1.0
# 		hll = 0.0
# 		hrr = 0.0
# 		hlr = 0.0
# 	elseif(OP[oc]==:-)
# 		hl = 1.0
# 		hr = -1.0
# 		hll = 0.0
# 		hrr = 0.0
# 		hlr = 0.0
# 	elseif(OP[oc]==:*)
# 		hl = rval
# 		hr = lval
# 		hll = 0.0
# 		hrr = 0.0
# 		hlr = 1.0
# 	elseif(OP[oc]==:/)
# 		hl = 1/rval
# 		hr = -lval/(rval^2)
# 		hll = 0.0
# 		hrr = 2*lval/(rval^3)
# 		hlr = -1/(rval^2)
# 	elseif(OP[oc]==:^)
# 		hl = rval*(lval^(rval-1))
# 		hr = 0.0
# 		hll = rval*(rval-1)*(lval^(rval-2))
# 		hrr = 0.0
# 		hlr = 0.0
# 	else
# 		assert(false)
# 	end
# 	assert(val!=NaN && hl != NaN && hll!=NaN && hrr!=NaN &&hlr!=NaN)
# 	return (val,hl,hr,hll,hrr,hlr)
# end


# function reverse_hess_calc(oc::OP_TYPE, lval::Real)
# 	hl  = NaN
# 	hll = NaN
# 	val = evaluate(OP[oc],lval)::Float64
# 	if(OP[oc]==:sin)
# 		hl = cos(lval)
# 		hll = -sin(lval)
# 	elseif(OP[oc]==:cos)
# 		hl = -sin(lval)
# 		hll = -cos(lval)
# 	else
# 		assert(false)	
# 	end
# 	assert(val!=NaN && hl!=NaN && hll!=NaN)
# 	return (val,hl,hll)
# end

function hess_nz(tt::TT_TYPE, idx::IDX_TYPE, eset::Set{Edge})
	# @show eset
	ntype = tt[idx]
	this_idx = idx
	idx -= 1
	if(ntype == TYPE_P)
		removeEdges(eset,this_idx)
	elseif(ntype == TYPE_V)
		#nothing
	elseif(ntype == TYPE_OB)
		# @show "TYPE_OB"
		oc = tt[idx]
		idx -= 1 #skip oc
		ridx = tt[idx]
		idx -= 1
		lidx = tt[idx]
		idx -= 1

		deset = Set{Edge}()
		aeset = Set{Edge}()
		
		#pushing
		state = start(eset)
		# i = 0
		while !done(eset,state)
			# @show state
			# @show eset
			# i += 1
			# println("iter ",i," - length ",length(eset))
			(e,state) = next(eset,state)
			if(isEndPointOnEdge(e,this_idx))
				if(isSelfEdge(e))
					e1 = createEdge(tt,lidx,ridx)
					e2 = createEdge(tt,lidx,lidx)
					e3 = createEdge(tt,ridx,ridx)
					push!(aeset, e1)
					push!(aeset, e2)
					push!(aeset, e3)
				else
					oidx = e.lidx == this_idx?e.ridx:e.lidx
					e1 = createEdge(tt,lidx,oidx)
					e2 = createEdge(tt,ridx,oidx)
					push!(aeset, e1)
					push!(aeset, e2)
				end
				# println("before delete iter",i," - length ",length(eset))
				# eset = delete!(eset,e)
				# println("after delete iter ",i," - length ",length(eset))
				push!(deset,e)
			end
		end
		setdiff!(eset, deset)  #remove deleted edge
		
		#creating
		if(OP[oc] == :+ || OP[oc] == :-)
			#do nothing for linear operator
		elseif(OP[oc] == :*)
			e1 = createEdge(tt,lidx,ridx)
			push!(aeset,e1)
		elseif(OP[oc] == :/)
			e1 = createEdge(tt,lidx,ridx)
			e3 = createEdge(tt,ridx,ridx)
			push!(aeset,e1)
			push!(aeset,e3)
		elseif(OP[oc] == :^)
			e1 = createEdge(tt,lidx,ridx)
			e2 = createEdge(tt,lidx,lidx)
			e3 = createEdge(tt,ridx,ridx)
			push!(aeset,e1)
			push!(aeset,e2)
			push!(aeset,e3)
		else
			assert(false);
		end
		union!(eset, aeset)  #ading pushed edges and created edges

		hess_nz(tt,lidx,eset)
		hess_nz(tt,ridx,eset)
	elseif(ntype == TYPE_OU)
		# @show "TYPE_OU"
		oc = tt[idx]
		idx -= 1 #skip oc
		lidx = tt[idx]

		deset = Set{Edge}()
		aeset = Set{Edge}()
		#pushing
		state = start(eset)
		while !done(eset,state)
			(e,state) = next(eset,state)
			if(isEndPointOnEdge(e,this_idx))
				if(isSelfEdge(e))
					#case 2
					e2 = createEdge(tt,lidx,lidx)
					push!(aeset,e2)
				else
					#case 1, 3
					oidx = e.lidx == this_idx? e.ridx:e.lidx
					e1 = createEdge(tt,lidx,oidx)
					push!(aeset,e1)
				end
				# delete!(eset, e)
				push!(deset,e)
			end
		end
		setdiff!(eset,deset)  #remove deleted edge
		#creating
		if(OP[oc]==:sin)
			e1 = createEdge(tt,lidx,lidx)
			push!(aeset,e1)
		elseif (OP[oc] == :cos)
			e1 = createEdge(tt,lidx,lidx)
			push!(aeset,e1)
		else
			assert(false)
		end
		union!(eset,aeset) #adding pushed edges and created edges

		hess_nz(tt,lidx,eset)
	else
		assert(false)
	end
end


function hess_nz1(eset::Set{Edge},veset::Set{Edge})
	# println("before hess_nz1", eset)
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
		# println("hess_nz1 - make edge ",lvidx," -- ",rvidx)
		push!(veset,ve)
	end
	# println("after hess_nz1", veset)
end

#Interface function
function hess_structure(tt::TT_TYPE, veset::Set{Edge})
	eset = Set{Edge}()
	hess_nz(tt, convert(IDX_TYPE,length(tt)), eset)
	hess_nz1(eset,veset)
end

function reverse_hess_ep(tt::TT_TYPE,vvals::TV_TYPE,pvals::TV_TYPE, veset::EdgeSet)
	reverse_hess_ep(tt::TT_TYPE,vvals::TV_TYPE,pvals::TV_TYPE, 1.0, veset::EdgeSet)
end

function reverse_hess_ep(tt::TT_TYPE,vvals::TV_TYPE,pvals::TV_TYPE, factor::VV_TYPE, veset::EdgeSet)
	ss = TV_STACK(VV_TYPE)
	reverse_hess_ep0(tt,convert(IDX_TYPE,length(tt)),ss,vvals,pvals)
	adj = 1.0
	eset = EdgeSet()
	reverse_hess_ep1(tt,convert(IDX_TYPE,length(tt)),ss,adj,eset)
	assert(isempty(ss) == true)
	reverse_hess_ep2(eset,factor,veset)
end