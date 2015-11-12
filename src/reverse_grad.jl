
# forward pass on the tape tt, to build ss stack

#forward pass on the scalar function
function forward_pass_1ord{V,I}(tape::Tape{I}, vvals::Array{V,1}, pvals::Array{V,1}, imm::Array{V,1})
	tt = tape.tt
	idx = one(I)
	v = [zero(V)]
	stk = Array{V,1}() #used for value evaluation
	empty!(imm)  #used for immediate derrivative
	sizehint!(stk,tape.maxoperands+20)
	sizehint!(v,1)

	@inbounds while(idx <= length(tt))
		# @show idx
		# println("++++++++++++++++++++++++++++++++++++")
		ntype = tt[idx]
		idx += 1
		if(ntype == TYPE_P)
			# tic()
			push!(stk,pvals[tt[idx]])	
			idx += 2 #skip TYPE_P
		elseif(ntype == TYPE_V)
			push!(stk,vvals[tt[idx]])
			idx += 2 #skip TYPE_V
		elseif(ntype == TYPE_O)
			oc = tt[idx]
			idx += 1
			n = tt[idx]
			idx += 1
			idx += 1 #skip TYPE_O
			# @show OP[oc],stk
			# @show OP[oc], length(stk)-n+1, n
			eval_1ord(OP[oc],stk,length(stk)-n+1,imm,v)
			# @show imm
			# @show v
			resize!(stk,length(stk)-n+1)
			stk[end] = v[1]
		end
		# println("++++++++++++++++++++++++++++++++++++")
	end
	# @show stk
end

function reverse_pass_1ord{V,I}(tape::Tape{I},imm::Array{V,1},g::Array{Tuple{I,V},1})
	tt = tape.tt
	idx = length(tt)
	adjs = Array{V,1}()
	sizehint!(adjs, tape.maxoperands+20)  #the acutally size should be (depth_of_tree - current_depth + maxoperands)
	push!(adjs,one(V))  #init value 

	@inbounds while(idx > 0)
		ntype = tt[idx]
		idx -= 1
		adj = pop!(adjs)
		# val = pop!(v)
		if(ntype == TYPE_P)
			idx -= 2
		elseif(ntype == TYPE_V)
			push!(g,(tt[idx],adj))
			idx -= 2
		elseif(ntype == TYPE_O)
			n = tt[idx]
			idx -= 3 #skip TYPE_O 
			@simd for i=length(imm)-n+1:1:length(imm)
				push!(adjs,imm[i]*adj)
			end
			resize!(imm,length(imm)-n)
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
	@inbounds @simd for i in 1:1:length(ilist)
		push!(iset,ilist[i])
	end
end

function grad_structure{I}(tape::Tape{I}, ilist::Array{I,1})  #repeat version
	empty!(ilist)
	grad_struct(tape,ilist)
end

function grad_reverse{V,I}(tape::Tape{I},vvals::Array{V,1},pvals::Array{V,1}, g::Array{Tuple{I,V},1}) #sparse version
	empty!(g)
	imm = Array{V,1}()
	sizehint!(imm,tape.nnode-1)			
	forward_pass_1ord(tape,vvals,pvals,imm)
	reverse_pass_1ord(tape,imm,g)
end

function grad_reverse{V,I}(tape::Tape{I},vvals::Array{V,1},pvals::Array{V,1}, g::Array{V,1})  #dense version
	grad = Array{Tuple{I,V},1}()
	sizehint!(grad,tape.nvnode)
	grad_reverse(tape,vvals,pvals,grad)
	
	fill!(g,zero(V))
	@inbounds @simd for i = 1:1:length(grad)
		g[grad[i][1]] += grad[i][2]
	end
end