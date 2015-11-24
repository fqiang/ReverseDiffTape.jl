#edge pusing algorithm for Hessian reverse AD


@inline function append_array{I,V}(dest::Vector{V},d_offset::I,src::Vector{V},s_offset::I, n::I)
	for i=1:n
		@inbounds dest[i+d_offset] = src[i+s_offset]
	end
end

@inline function push_diag{I,V}(eset::Dict{I,Dict{I,V}},i1::I)
	# @show i1
	# assert(haskey(eset,i1))
	eset[i1][i1] = 0.0
end

@inline function push_edge{I,V}(eset::Dict{I,Dict{I,V}},i1::I,i2::I)
	# @show i1,i2
	# assert(i1!=i2)
	# @show "push_edge", i1,i2
	if i1<i2
		# assert(haskey(eset,i2))
		eset[i2][i1] = 0.0
	else
		# assert(haskey(eset,i1))
		eset[i1][i2] = 0.0
	end
end

@inline function push_live_var{I}(liveVar::Dict{I,Set{I}},i1::I,i2::I)
	# assert(haskey(liveVar,i1))
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
	
	@inbounds while (idx > 0)
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
			
			for p in lvi  #for each 
				if(i==p)
					for j0=trlen-n+1:trlen  #construct live var set for each children
						ci = tr[j0] #child of i
						push_live_var(tape.liveVar,ci,ci)
						push_diag(tape.eset,ci)
						for j1=j0+1:trlen
							cii = tr[j1] #child of i
							# assert(ci<cii) #ci always smaller since tr building
							push_live_var(tape.liveVar,ci,cii)
							push_live_var(tape.liveVar,cii,ci)
							push_edge(tape.eset,ci,cii)
						end 
					end
				else
					for j0 = trlen-n+1:trlen	
						ci = tr[j0]
						push_live_var(tape.liveVar,ci,p)
						push_live_var(tape.liveVar,p,ci)
						push_edge(tape.eset,ci,p)
					end
				end
			end

			#creating
			if (n==1)  #for 1-ary operator
				# @show OP[oc],d
				ci = tr[trlen]
				push_live_var(tape.liveVar,ci,ci)
				push_diag(tape.eset,ci)
			elseif (OP[oc] ==:+ || OP[oc] ==:-)
				#zeros
			elseif (OP[oc] == :*) #special case for *
				for j0=trlen-n+1:trlen  #construct live var set for each children
					ci = tr[j0]
					for j1=j0+1:trlen
						cii = tr[j1]
						# assert(ci<cii)
						# @show OP[oc],ci,cii
						push_live_var(tape.liveVar,ci,cii)
						push_live_var(tape.liveVar,cii,ci)
						push_edge(tape.eset,ci,cii)
					end
				end
			else #binary
				ri = tr[trlen]
				li = tr[trlen-1]
				# assert(li<ri)
			
				push_live_var(tape.liveVar,li,li)
				push_diag(tape.eset,li) #dxx
				
				push_live_var(tape.liveVar,li,ri)
				push_live_var(tape.liveVar,ri,li)
				push_edge(tape.eset,li,ri) #dxy
				
				push_live_var(tape.liveVar,ri,ri)
				push_diag(tape.eset,ri) #dyy
			end
			trlen -= n
		end
	end
	# @show tape.liveVar
	# @show tape.eset

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
				# @show ii,jj
				if(lo)
					if(!haskey(eset,ii))
						eset[ii] = Set{I}()
					end
					push!(eset[ii],jj)
				else
					if(!haskey(eset,jj))
						eset[jj] = Set{I}()
					end
					push!(eset[jj],ii)
				end
			end
		end
	end	
end

@inline function getw{I,V}(eset::Dict{I,Dict{I,V}},i1::I,i2::I)
	# @show "in w",i1,i2
	if(i1>=i2)
		# assert(haskey(eset,i1))
		# assert(haskey(eset[i1],i2))
		return eset[i1][i2]
	else
		# assert(haskey(eset,i2))
		# assert(haskey(eset[i2],i1))
		return eset[i2][i1]
	end
end

@inline function incr{I,V}(eset::Dict{I,Dict{I,V}},i1::I,i2::I,w::V)
	# @show "incr",i1,i2, w
	if i1>=i2 
		# assert(haskey(eset,i1))
		# assert(haskey(eset[i1],i2))
		eset[i1][i2]+=w 
	else
		# assert(haskey(eset,i2))
		# assert(haskey(eset[i2],i1))
		eset[i2][i1]+=w
	end
end


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
	# @show tape.imm2ord,immlen
	# assert(tape.imm2ord>=immlen)
	tape.imm2ord = immlen
	# @show stklen
	resize!(imm,immlen)
	return stk[1]
end


function reverse_pass_2ord{I,V}(tape::Tape{I,V}, imm::Array{V,1}, factor::V, eset::Dict{I,Dict{I,V}})
	tr = tape.tr
	tt = tape.tt
	idx = length(tt)
	trlen = length(tr)
	immlen = tape.imm2ord
	assert(length(imm) == immlen)

	vidx = Set{I}()

	adjs = Vector{V}(tape.maxoperands+20)
	adjlen = 1
	adjs[1] = one(V)

	@inbounds while(idx > 0)
		# println("++++++++++++++++++++++++++++++++++++")
		# @show idx
		# @show trlen,immlen, adjlen
		ntype = tt[idx]
		idx -= 1
		adj = adjs[adjlen]
		adjlen -= 1
		# @show adj

		if(ntype == TYPE_P)
			idx -= 2
		elseif(ntype == TYPE_V)
			idx -= 2
			push!(vidx, idx+1)
			# @show idx, vidx
		elseif(ntype == TYPE_O)
			n = tt[idx]
			idx -= 1
			oc = tt[idx]
			idx -= 2

			# @show OP[oc],n
			# @show tr

			#pushing
			i = idx + 1 #current node idx
			lvi = tape.liveVar[i] #live var set at i

			for p in lvi  #for each upper live vars
				# @show p, i
				w = getw(tape.eset,i,p)
				if(i==p)
					if(n==1) #1-ary operator
						# @show "pushing 1-ary", OP[oc],tr[trlen]
						incr(tape.eset,tr[trlen],tr[trlen],imm[immlen-1]*imm[immlen-1]*w)
					else  #2 or more
						if(OP[oc]==:+ )
							for k=trlen-n+1:trlen
								incr(tape.eset,tr[k],tr[k],w)
								j0 = j + 1
								for k0=k+1:trlen
									incr(tape.eset,tr[k],tr[k0],w)
									j0 += 1
								end
								j += 1
							end
						elseif(OP[oc] ==:-)
							l = tr[trlen-1]
							r = tr[trlen]
							incr(tape.eset,l,l,w)
							incr(tape.eset,r,r,w)
							incr(tape.eset,l,r,-1.0*w)
						elseif(OP[oc] == :*)
							j = immlen - round(I,n+n*(n-1)/2)+1
							for k=trlen-n+1:trlen
								incr(tape.eset,tr[k],tr[k],imm[j]*imm[j]*w)
								j0 = j + 1
								for k0=k+1:trlen
									incr(tape.eset,tr[k],tr[k0],imm[j]*imm[j0]*w)
									j0 += 1
								end
								j += 1
							end
						else #binary
							l = tr[trlen-1]
							r = tr[trlen]
							dl = imm[immlen-4]
							dr = imm[immlen-3]
							incr(tape.eset,l,l,dl*dl*w)
							incr(tape.eset,r,r,dr*dr*w)
							incr(tape.eset,l,r,dl*dr*w)
						end
					end
				else
					if(n==1)
						# @show trlen,immlen
						# @show imm
						# @show tr
						# @show imm[immlen-1], tr[trlen], w
						incr(tape.eset,tr[trlen],p,imm[immlen-1]*imm[immlen-1]*w)
					else
						if(OP[oc]==:+)
							for k=trlen-n+1:trlen
								incr(tape.eset,tr[k],p,w)
								# assert(p!=tr[k])
							end
						elseif(OP[oc]==:-)
							l = tr[trlen-1]
							r = tr[trlen]
							incr(tape.eset,l,p,w)
							incr(tape.eset,r,p,-1.0*w)
							# assert(p!=l && p!=r)
						elseif(OP[oc] == :*)
							j = immlen - round(I,n+n*(n-1)/2)+1
							for k=trlen -n+1:trlen
								incr(tape.eset,tr[k],p,imm[j]*w)
								# assert(p!=tr[k])
							end
						else #binary
							l = tr[trlen-1]
							r = tr[trlen]
							dl = imm[immlen-4]
							dr = imm[immlen-3]
							incr(tape.eset,l,p,dl*w)
							incr(tape.eset,r,p,dr*w)
							# assert(p!=l && p!=r)
						end
					end
				end
			end #end pushing

			#creating
			if n==1
				# @show "creating 1-ary", OP[oc],tr[trlen]
				incr(tape.eset,tr[trlen],tr[trlen],adj*imm[immlen])
			else
				if(OP[oc] == :+ || OP[oc] ==:-)
					#zero
				elseif(OP[oc]==:*)
					j = immlen - round(I,n*(n-1)/2) + 1
					for k=trlen-n+1:trlen
						for k0=k+1:trlen
							# @show "creating n-ary",n, OP[oc],tr[trlen] 
							incr(tape.eset,tr[k],tr[k0],adj*imm[j])
							j+=1
						end
					end
				else #binary
					l = tr[trlen-1]
					r = tr[trlen]
					dll = imm[immlen-2]
					dlr = imm[immlen-1]
					drr = imm[immlen]
					incr(tape.eset,l,l,adj*dll)
					incr(tape.eset,l,r,adj*dlr)
					incr(tape.eset,r,r,adj*drr)
				end
			end

			# @show adjlen
			#adj
			imm_counter = zero(I)
			if n==1
				adjlen += 1
				adjs[adjlen] = imm[immlen-1]*adj
				imm_counter = 2
			else
				if OP[oc]==:+ 
					for m=1:n
						adjlen += 1
						adjs[adjlen] = adj
					end
				elseif OP[oc] ==:-
					adjlen += 1
					adjs[adjlen] = adj
					adjlen += 1
					adjs[adjlen] = -1.0*adj
				elseif OP[oc] ==:*
					j=immlen-round(I,n+n*(n-1)/2)+1
					# @show immlen,round(I,n+n*(n-1)/2),j
					for m=1:n
						# @show adjlen, j, m
						adjlen += 1
						adjs[adjlen] = imm[j]*adj
						j+=1
					end
					# @show immlen, adjlen
					imm_counter = round(I,n+n*(n-1)/2)
				else
					adjlen += 1
					adjs[adjlen] = imm[immlen-4]*adj
					adjlen += 1
					adjs[adjlen] = imm[immlen-3]*adj
					imm_counter = 5
				end
			end

			# @show OP[oc],n
			# @show trlen
			# @show immlen,imm_counter
			# @show tr
			# @show imm


			#update
			trlen -= n
			immlen -= imm_counter
		end #end TYPE_O
		# println("++++++++++++++++++++++++++++++++++++")
	end  #end while
	assert(immlen == 0 && trlen == 0)
	# @show tape.eset

	# @show vidx
	for i in vidx
		# @show i
		eseti = tape.eset[i]
		# @show eseti
		for (j,w) in eseti
			# @show j
			ii = tt[i+1]
			if tt[j] == TYPE_V
				jj = tt[j+1]
				if(!haskey(eset,ii))
					eset[ii] = Dict{I,V}()
				end
				if(i!=j && ii==jj)
					haskey(eset[ii],jj)?eset[ii][jj] += 2.0*w*factor : eset[ii][jj] = 2.0*w*factor
				else
					haskey(eset[ii],jj)?eset[ii][jj] += w*factor: eset[ii][jj] = w*factor
				end
				# @show ii,jj, eset[ii][jj]
			end
		end
	end	
end



#Interface function
function hess_structure_lower{I,V}(tape::Tape{I,V}, eset::Dict{I,Set{I}})
	hess_struct(tape,eset,true)
end

function hess_reverse{I,V}(tape::Tape{I,V},vvals::Vector{V},pvals::Vector{V},eset::Dict{I,Dict{I,V}})
	hess_reverse(tape,vvals,pvals,1.0,eset)
end
function hess_reverse{I,V}(tape::Tape{I,V},vvals::Vector{V},pvals::Vector{V}, factor::V,eset::Dict{I,Dict{I,V}})
	imm = Vector{V}(tape.imm2ord)
	forward_pass_2ord(tape,vvals,pvals,imm)
	reverse_pass_2ord(tape,imm,factor,eset)
end
