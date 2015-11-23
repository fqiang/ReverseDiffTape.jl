typealias IDX_TYPE Int

#edge pusing algorithm for Hessian reverse AD

typealias EdgeSet{I,V} Dict{I,Dict{I,V}}

function forward_pass_2ord{I,V}(tape::Tape{I,V}, vvals::Array{V,1}, pvals::Array{V,1}, imm::Array{V,1})
	tt = tape.tt
	idx = one(I)
	# empty!(imm)  #used for immediate derrivative
	# empty!(tr)
	stk = Vector{V}(tape.maxoperands+20) #used for value evaluation
	stklen = zero(I)
	immlen = zero(I)
	
	@inbounds while(idx <= length(tt))
		# @show idx
		# println("++++++++++++++++++++++++++++++++++++")
		ntype = tt[idx]
		# eset[idx] = Dict{I,V}() #initialize edge set
		idx += 1
		if(ntype == TYPE_P)
			stklen += 1
			stk[stklen] = pvals[tt[idx]]
			idx += 2 #skip TYPE_P
		elseif(ntype == TYPE_V)
			stklen += 1
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
			else
				@inbounds (counter,val) = eval_2ord(OP[oc],stk,stklen-n+1,stklen,imm,immlen+1)
				stklen -= n-1
				@inbounds stk[stklen] = val
			end
			immlen += counter
			# @show imm		
			# @show stk[stklen]	
		end
		# @show stklen
		# println("++++++++++++++++++++++++++++++++++++")
	end
	tape.imm2ord = immlen
	# @show stklen
	return stk[1]
end


@inline function append_array{I,V}(dest::Vector{V},d_offset::I,src::Vector{V},s_offset::I, n::I)
	for i=1:n
		@inbounds dest[i+d_offset] = src[i+s_offset]
	end
end

@inline function push_diag{I,V}(eset::Dict{I,Dict{I,V}},i1::I)
	# @show i1
	eset[i1][i1] = 0.0
end

@inline function push_edge{I,V}(eset::Dict{I,Dict{I,V}},i1::I,i2::I,lo::Bool)
	# @show i1,i2,lo
	assert(i1!=i2)
	if(lo) #lower
		i1<i2? eset[i2][i1] = 0.0:eset[i2][i1] = 0.0
	else  #upper
		i1<i2? eset[i1][i2] = 0.0:eset[i1][i2] = 0.0
	end
end

@inline function push_live_var{I}(liveVar::Dict{I,Set{I}},i1::I,i2::I)
	push!(liveVar[i1],i2)
end

function hess_struct{I,V}(tape::Tape{I,V},eset::Dict{I,Set{I}},lo::Bool)
	tt = tape.tt
	tr = tape.tr
	idx = length(tt)
	trlen = length(tr)
	vidx = Set{I}()
	# @show tape.nnode-1, length(tr)
	assert(tape.nnode-1 == length(tr))
	# assert(isempty(eset))
	
	while (idx > 0)
		ntype = tt[idx]
		idx -= 1
		if(ntype == TYPE_P)
			idx -= 2
		elseif(ntype == TYPE_V)
			idx -= 2
			push!(vidx,idx+1)
		elseif(ntype == TYPE_O)
			n = tt[idx]
			idx -= 1
			oc = tt[idx]
			idx -= 2

			#pushing
			i = idx + 1
			lvi = tape.liveVar[i] #live var set at i
			# for j=trlen-n+1:trlen  #construct live var set for each children
			# 	ci = tr[j]  #child of i
			# 	tape.liveVar[ci] = Dict{I,V}() 
			# end

			for p in lvi  #for each 
				if(i==p)
					for j0=trlen-n+1:trlen  #construct live var set for each children
						ci = tr[j0] #child of i
						push_live_var(tape.liveVar,ci,ci)
						push_diag(tape.eset,ci)
						for j1=j0+1:trlen
							cii = tr[j1] #child of i
							assert(ci<cii) #ci always smaller since tr building
							push_live_var(tape.liveVar,ci,cii)
							push_live_var(tape.liveVar,cii,ci)
							push_edge(tape.eset,ci,cii,lo)
						end 
					end
				else
					for j0 = trlen-n+1:trlen	
						ci = tr[j0]
						push_live_var(tape.liveVar,ci,p)
						push_edge(tape.eset,ci,p,lo)
					end
				end
			end

			#creating
			if (U_OP_START<=oc<=U_OP_END)  #for 1-ary operator
				# @show OP[oc],d
				d = S_TO_DIFF_FLAG[OP[oc]]
				if(d[1]==false) #not zero
					ci = tr[trlen]
					push_live_var(tape.liveVar,ci,ci)
					push_diag(tape.eset,ci)
				end
			elseif (OP[oc] == :*) #special case for *
				for j0=trlen-n+1:trlen  #construct live var set for each children
					ci = tr[j0]
					for j1=j0+1:trlen
						cii = tr[j1]
						assert(ci<cii)
						# @show OP[oc],ci,cii
						push_live_var(tape.liveVar,ci,cii)
						push_live_var(tape.liveVar,cii,ci)
						push_edge(tape.eset,ci,cii,lo)
					end
				end
			elseif (OP[oc] ==:+ || OP[oc] ==:-)
				#skil no second order
			else #+ has no edges created
				d = S_TO_DIFF_FLAG[OP[oc]]
				ri = tr[trlen]
				li = tr[trlen-1]
				assert(li<ri)
				# @show OP[oc],d
				if(d[1] == false)
					# push!(tape.liveVar[li],li=>0.0)
					push_live_var(tape.liveVar,li,li)
					push_diag(tape.eset,li) #dxx
				end
				if(d[2] == false)
					push_live_var(tape.liveVar,li,ri)
					push_live_var(tape.liveVar,ri,li)
					push_edge(tape.eset,li,ri,lo) 
				end
				if(d[3] == false)
					push_live_var(tape.liveVar,ri,ri)
					push_diag(tape.eset,ri)
				end
			end
			trlen -= n
		end
	end

	# @show vidx
	for i in vidx
		# @show i
		eseti = tape.eset[i]
		# @show eseti
		for j in keys(eseti)
			# @show j
			ii = tt[i+1]
			if tt[j] == TYPE_V
				jj = tt[j+1]
				if(!haskey(eset,ii))
					eset[ii] = Set{I}()
				end
				# @show ii,jj
				push!(eset[ii],jj)
			end
		end
	end	
end


function reverse_pass_2ord{I,V}(tape::Tape{I,V},imm::Array{V,1},eset::EdgeSet{I,V})
	tr = tape.tr
	tt = tape.tt
	idx = length(tt)
	immlen = tt.imm2ord
	vidx = Set{I}()
	assert(length(imm) == immlen)

	adjs = Vector{V}(tape.maxoperands)
	adjlen = 1
	adjs[1] = one(V)

	@inbounds while(idx > 0)
		ntype = tt[idx]
		idx -= 1
		adj = adjs[adjlen]
		adjlen -= 1

		if(ntype == TYPE_P)
			idx -= 2
		elseif(ntype == TYPE_V)
			idx -= 2
			push!(vidx,tt[idx+1])
		elseif(ntype == TYPE_O)
			n = tt[idx]
			idx -= 1
			oc = tt[idx]
			idx -= 2

			#pushing
			i = idx + 1 #current node idx
			lvi = tape.liveVar[i] #live var set at i

			for p in lvi  #for each 
				if(i==p)

				else
					
				end
			end

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


# function reverse_hess_ep1(tt::TT_TYPE, idx::IDX_TYPE, ss::TV_STACK, adj::VV_TYPE, eset::EdgeSet)
# 	# debug("enter - ", adj)
# 	this_idx = idx
# 	ntype = tt[idx]
# 	idx -= 1
# 	if(ntype ==TYPE_P)  
# 		#deleting all edges that has an endpoint to this node
# 		removeEdges(eset,this_idx)
# 	elseif(ntype == TYPE_V)
# 		#do nothing
# 	elseif(ntype == TYPE_OB)
# 		hlr = pop!(ss)
# 		hrr = pop!(ss)
# 		hll = pop!(ss)
# 		hr = pop!(ss)
# 		hl = pop!(ss)

# 		oc = tt[idx]
# 		assert(B_OP_START<= oc <= B_OP_END)
# 		idx -= 1
# 		ridx = tt[idx]
# 		idx -= 1
# 		# rt = tt[ridx]
# 		lidx = tt[idx]
# 		idx -= 1
# 		# lt = tt[lidx]
		
# 		deset = Set{Edge}()
# 		aeset = EdgeSet()

# 		#pushing
# 		state = start(eset)
# 		# println("binary op before pushing - eset(",length(eset),")")
# 		while !done(eset,state)
# 			(pair, state) = next(eset,state)
# 			e = pair[1]
# 			w = pair[2]
# 			if(isEndPointOnEdge(e,this_idx))
# 				if(isSelfEdge(e))
# 					#case 2
# 					e1 = createEdge(tt,lidx,ridx)
# 					w1 = hl*hr*w
# 					e2 = createEdge(tt,lidx,lidx)
# 					w2 = hl*hl*w
# 					e3 = createEdge(tt,ridx,ridx)
# 					w3 = hr*hr*w
# 					aeset[e1] = haskey(aeset,e1)?aeset[e1]+w1:w1
# 					aeset[e2] = haskey(aeset,e2)?aeset[e2]+w2:w2
# 					aeset[e3] = haskey(aeset,e3)?aeset[e3]+w3:w3
# 				else
# 					#case 1 , 3
# 					oidx = e.lidx == this_idx? e.ridx: e.lidx
# 					e1 = createEdge(tt,lidx,oidx)
# 					e2 = createEdge(tt,ridx,oidx)
# 					w1 = e1.lidx == e1.ridx? (2*hl*w): (hl*w)
# 					w2 = e2.lidx == e2.ridx? (2*hr*w): (hr*w)
# 					aeset[e1] = haskey(aeset,e1)? aeset[e1]+w1:w1
# 					aeset[e2] = haskey(aeset,e2)? aeset[e2]+w2:w2
# 				end
# 				# delete!(eset,e)
# 				push!(deset,e)
# 			end
# 		end

# 		for k in deset  #delete the pushed edges
# 			delete!(eset,k)
# 		end
# 		# println("binary op after pushing - eset(",length(eset),")")
# 		# println("before remove set rmv(",length(rmv),") - eset(",length(eset),")")
# 		# setdiff!(eset, rmv)
# 		# println("after remove set rmv(",length(rmv),") - eset(",length(eset),")")

# 		#creating
# 		if(OP[oc] == :+ || OP[oc] == :-)
# 			#do nothing for linear operator
# 		elseif(OP[oc] == :*)
# 			e1 = createEdge(tt,lidx,ridx)
# 			w1 = adj*hlr
# 			aeset[e1] = haskey(aeset,e1)?aeset[e1]+w1:w1
# 			# eset[e1] += (adj*hlr)
# 		elseif(OP[oc] == :/)
# 			e1 = createEdge(tt,lidx,ridx)
# 			w1 = adj*hlr
# 			e3 = createEdge(tt,ridx,ridx)
# 			w3 = adj*hrr
# 			aeset[e1] = haskey(aeset,e1)?aeset[e1]+w1:w1
# 			aeset[e3] = haskey(aeset,e3)?aeset[e3]+w3:w3
# 			# eset[e1] += (adj*hlr)
# 			# eset[e3] += (adj*hrr)
# 		elseif(OP[oc] == :^)
# 			e1 = createEdge(tt,lidx,ridx)
# 			w1 = adj*hlr
# 			e2 = createEdge(tt,lidx,lidx)
# 			w2 = adj*hll
# 			e3 = createEdge(tt,ridx,ridx)
# 			w3 = adj*hrr
# 			aeset[e1] = haskey(aeset,e1)?aeset[e1]+w1:w1
# 			aeset[e2] = haskey(aeset,e2)?aeset[e2]+w2:w2
# 			aeset[e3] = haskey(aeset,e3)?aeset[e3]+w3:w3
# 			# eset[e1] += (adj*hlr)
# 			# eset[e2] += (adj*hll)
# 			# eset[e3] += (adj*hrr)
# 		else
# 			assert(false);
# 		end
		
# 		for (k,v) in aeset  #adding the pushed and created edges
# 			# println("reverse_hess_ep1 - add, ",k," weight ",v)
# 			eset[k] = haskey(eset,k)?eset[k]+v:v
# 			# @show eset
# 		end

# 		#adjoint
# 		ladj = adj*hl;
# 		radj = adj*hr;

# 		#recursive
# 		reverse_hess_ep1(tt, lidx, ss, ladj, eset)
# 		reverse_hess_ep1(tt, ridx, ss, radj, eset)
# 	elseif(ntype == TYPE_OU)
# 		hll = pop!(ss)
# 		hl = pop!(ss)
# 		oc = tt[idx]
# 		idx -= 1
# 		assert(U_OP_START<= oc <= U_OP_END)
# 		lidx = tt[idx]
# 		idx -= 1

# 		deset = Set{Edge}()
# 		aeset = EdgeSet()

# 		#pushing
# 		state = start(eset)
# 		while !done(eset,state)
# 			(pair, state) = next(eset,state)
# 			e = pair[1]
# 			w = pair[2]
# 			if(isEndPointOnEdge(e,this_idx))
# 				if(isSelfEdge(e))
# 					#case 2
# 					e2 = createEdge(tt,lidx,lidx)
# 					w2 = hl*hl*w
# 					aeset[e2] = haskey(aeset,e2)?aeset[e2]+w2:w2
# 					# eset[e2] += hl*hl*w
# 				else
# 					#case 1, 3
# 					oidx = e.lidx == this_idx? e.ridx:e.lidx
# 					e1 = createEdge(tt,lidx,oidx)
# 					w1 = e1.lidx == e1.ridx? 2*hl*w: hl*w
# 					aeset[e1] = haskey(aeset,e1)?aeset[e1]+w1:w1
# 				end
# 				# delete!(eset, e)
# 				push!(deset,e)
# 			end
# 		end

# 		for k in deset  #delete the pushed edges
# 			delete!(eset,k)
# 		end

# 		#creating
# 		if(OP[oc]==:sin)
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

# 		for (k,v) in aeset  #adding the pushed and created edges
# 			eset[k] = haskey(eset,k)?eset[k]+v:v
# 		end

# 		#adjoint
# 		ladj = adj*hl;
		
# 		#recursive
# 		reverse_hess_ep1(tt,lidx,ss,ladj,eset)
# 	else
# 		assert(false)
# 	end
# 	# debug("exit")
# end


# function reverse_hess_ep2(eset::EdgeSet,factor::VV_TYPE, veset::EdgeSet)
# 	# veset = EdgeSet()
# 	# @show eset
# 	state = start(eset)
# 	while !done(eset,state)
# 		(pair, state) = next(eset, state)
# 		e = pair[1]
# 		w = pair[2]
# 		lidx = e.lidx
# 		ridx = e.ridx
# 		assert(e.tt[lidx] == TYPE_V)
# 		assert(e.tt[ridx] == TYPE_V)
# 		lidx -= 1
# 		ridx -= 1
# 		lvidx = e.tt[lidx]
# 		rvidx = e.tt[ridx]
# 		if(lidx != ridx && lvidx == rvidx) 
# 			w += w
# 		end
# 		# println("reverse_hess_ep2 - adding (",lvidx,",",rvidx,") = ",w)
# 		ve = Edge(e.tt,lvidx,rvidx)
# 		# println("reverse_hess_ep2 - make edge ",lvidx," -- ",rvidx)

# 		if(!haskey(veset,ve))
# 			veset[ve] = w*factor
# 		else
# 			veset[ve] += w*factor
# 		end
# 	end
# 	# @show veset
# 	return veset
# end



#Interface function
function hess_structure_lower{I,V}(tape::Tape{I,V}, eset::Dict{I,Set{I}})
	hess_struct(tape,eset,true)
end

function hess_reverse{I,V}(tape::Tape{I,V},vvals::Vector{V},pvals::Vector{V}, eset::EdgeSet{I,V})
	imm = Vector{V}(tape.imm2ord)
	forward_pass_2ord(tape,vvals,pvals,imm)
	reverse_pass_2ord(tape,imm,eset)
end
