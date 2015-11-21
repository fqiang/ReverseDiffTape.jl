
# forward pass on the tape tt, to build ss stack

#forward pass on the scalar function
function forward_pass_1ord{V,I}(tape::Tape{I}, vvals::Array{V,1}, pvals::Array{V,1}, imm::Array{V,1})
	tt = tape.tt
	idx = one(I)
	stk = Vector{V}(tape.maxoperands+20) #used for value evaluation
	stklen = 0
	immlen = 0

	@inbounds while(idx <= length(tt))
		# @show idx
		# println("++++++++++++++++++++++++++++++++++++")
		ntype = tt[idx]
		idx += 1
		if(ntype == TYPE_P)
			# tic()
			val = pvals[tt[idx]]
			idx += 2 #skip TYPE_P
			stklen += 1
			stk[stklen] = val
		elseif(ntype == TYPE_V)
			val = vvals[tt[idx]]
			idx += 2 #skip TYPE_V
			stklen += 1
			stk[stklen] = val
		elseif(ntype == TYPE_O)
			oc = tt[idx]
			idx += 1
			n = tt[idx]
			idx += 1
			idx += 1 #skip TYPE_O
			# @show OP[oc], stklen-n+1, stklen
			# @show OP[oc], stk[stklen-n+1:stklen]
			if(n==1)
				@inbounds stk[stklen] = eval_1ord(OP[oc],stk[stklen],imm,immlen+1)
			else
				@inbounds val = eval_1ord(OP[oc],stk,stklen-n+1,stklen,imm,immlen+1)
				stklen -= n-1
				@inbounds stk[stklen] = val
			end
			# @show stk[stklen]
			# @show imm[immlen+1:immlen+n]
			immlen += n
			# @show immlen
			# @show imm
		end
		# @show stklen
		# @show stk
		# println("++++++++++++++++++++++++++++++++++++")
	end
	return stk[1]
end

function reverse_pass_1ord{V,I}(tape::Tape{I},imm::Array{V,1},g::Array{Tuple{I,V},1})
	# @show imm
	# assert(length(imm) == tape.nnode -1)
	tt = tape.tt
	idx = length(tt)
	immlen = length(imm)
	adjs = Vector{V}(tape.nnode)
	adjlen = 1
	adjs[adjlen] = one(V)

	@inbounds while(idx > 0)
		ntype = tt[idx]
		idx -= 1
		adj = adjs[adjlen]
		adjlen -= 1
		if(ntype == TYPE_P)
			idx -= 2
		elseif(ntype == TYPE_V)
			# @show tt[idx],adj
			# adj=isnan(adj)?0.0:adj
			push!(g,(tt[idx],adj))
			idx -= 2
		elseif(ntype == TYPE_O)
			n = tt[idx]
			idx -= 3 #skip TYPE_O 
			@simd for i=immlen-n+1:immlen
				adjlen += 1
				adjs[adjlen] = imm[i]*adj
			end
			immlen -= n
		end
	end
end

function grad_struct{I}(tape::Tape{I}, ilist::Array{I,1}) #repeated indexes, in reverse tracing order
	tt = tape.tt
	idx = length(tt)
	@inbounds while(idx > 0)
		ntype = tt[idx]
		idx -= 1
		if(ntype == TYPE_V)
			push!(ilist,tt[idx])
			idx -= 2
		elseif(ntype == TYPE_P)
			idx -= 2
		elseif(ntype == TYPE_O)
			idx -= 3
		end
	end
end

#Interface function
function grad_structure{I}(tape::Tape{I}, iset::Set{I}) #non repeat version
	ilist = Array{I,1}()
	sizehint!(ilist, tape.nvnode)
	grad_struct(tape,ilist)

	empty!(iset)
	@inbounds for i in 1:length(ilist)
		push!(iset,ilist[i])
	end
end

function grad_structure{I}(tape::Tape{I}, ilist::Array{I,1})  #repeat version
	empty!(ilist)
	grad_struct(tape,ilist)
end

function grad_reverse{V,I}(tape::Tape{I},vvals::Array{V,1},pvals::Array{V,1}, g::Array{Tuple{I,V},1}) #sparse version
	imm = Vector{V}(tape.nnode-1)	
	forward_pass_1ord(tape,vvals,pvals,imm)
	reverse_pass_1ord(tape,imm,g)
end

function grad_reverse{V,I}(tape::Tape{I},vvals::Array{V,1},pvals::Array{V,1}, g::Array{V,1})  #dense version
	grad = Array{Tuple{I,V},1}()
	sizehint!(grad,tape.nvnode)
	grad_reverse(tape,vvals,pvals,grad)
	
	# @show grad
	fill!(g,zero(V))
	# @show g
	@inbounds @simd for i = 1:length(grad)
		g[grad[i][1]] += grad[i][2]
	end
	# @show g
end