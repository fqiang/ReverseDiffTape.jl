
# forward pass on the tape tt, to build ss stack

#forward pass on the scalar function
function forward_pass{V,I}(tape::Tape{I}, vvals::Array{V,1}, pvals::Array{V,1}, imm::Array{V,1})
	tt = tape.tt
	idx = one(I)
	v = [zero(V)]
	stk = Array{V,1}() #used for value evaluation
	empty!(imm)  #used for immediate derrivative
	sizehint!(stk,tape.maxoperands)
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
			eval_idd(OP[oc],stk,length(stk)-n+1,imm,v)
			# @show imm
			# @show v
			resize!(stk,length(stk)-n+1)
			stk[end] = v[1]
		end
		# println("++++++++++++++++++++++++++++++++++++")
	end
	# @show stk
end

#show be auto generated code
function eval_idd{V,I}(s::Symbol,v::Array{V,1},i::I,imm::Array{V,1},r::Array{V,1})
	if(s==:sin)
		ld = cos(v[i])
		push!(imm,ld)
		r[1] = sin(v[i])
	elseif(s==:cos)
		ld = -sin(v[i])
		push!(imm,ld)
		r[1] = cos(v[i])
	elseif(s ==:+)
		ret = zero(V)
		@simd for j=i:1:length(v)
			push!(imm,one(V))
			ret += v[j]
		end
		r[1] = ret
	elseif(s ==:-)
		push!(imm,one(V))
		push!(imm,-one(V))
		r[1] = (-)(v[i],v[i+1])
	elseif(s==:*)
		# for i=length(v)-n+1:1:length(v)
			# ret = one(V)
			# for j=length(v)-n+1:1:i
			# 	ret *= v[j]
			# end
			# for j=i+1:1:length(v)
			# 	ret *= v[j]
			# end
		ret = one(V)
		@simd for j = i:1:length(v)
			ret *= v[j]
		end
		r[1] = ret
		@simd for j = i:1:length(v)
			push!(imm,r[1]/v[j])
		end
	elseif(s==:/)
		r[1] = (/)(v[i],v[i+1])
		push!(imm,one(V)/v[i+1])
		push!(imm,-v[i]/v[i+1]^2)
		
	elseif(s==:^)
		r[1] = (^)(v[i],v[i+1])
		push!(imm,v[i+1]*(v[i]^(v[i+1]-1)))
		push!(imm,r[1]*(log(v[i])))
	end
end


function reverse_pass{V,I}(tape::Tape{I},imm::Array{V,1},g::Array{Tuple{I,V},1})
	tt = tape.tt
	idx = length(tt)
	adjs = Array{V,1}()
	sizehint!(adjs, round(I,length(tt)/10))  #random guessing size, the acutally size should be (depth_of_tree - current_depth + maxoperands)
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

# function grad_struct{I}(tt::Array{I,1}, idx::IDX_TYPE, I::Set{IDX_TYPE}) #non-repeat version
# 	assert(length(vvals) == length(grad))
# 	# debug("enter - ", adj)
# 	ntype = tt[idx]
# 	idx -= 1
# 	if(ntype ==TYPE_P)  
# 		#do nothing , this is a parameter
# 	elseif(ntype == TYPE_V)
# 		vidx = tt[idx]
# 		push!(I,vidx)
# 	elseif(ntype == TYPE_OB)
# 		idx -= 1
# 		ridx = tt[idx]
# 		idx -= 1
# 		lidx = tt[idx]
# 		idx -= 1
# 		nzg(tt,lidx,I)
# 		nzg(tt,ridx,I)
# 	elseif(ntype == TYPE_OU)
# 		idx -= 1 #skip oc
# 		lidx = tt[idx]
# 		idx -= 1
# 		nzg(tt,lidx,I)
# 	else
# 		assert(false)
# 	end
# end

# function grad_struct(tt::TT_TYPE, idx::IDX_TYPE, I::Array{IDX_TYPE,1}) #repeat version
# 	# debug("enter - ", adj)
# 	ntype = tt[idx]
# 	idx -= 1
# 	if(ntype ==TYPE_P)  
# 		#do nothing , this is a parameter
# 	elseif(ntype == TYPE_V)
# 		vidx = tt[idx]
# 		push!(I,vidx) 
# 	elseif(ntype == TYPE_OB)
# 		idx -= 1
# 		ridx = tt[idx]
# 		idx -= 1
# 		lidx = tt[idx]
# 		idx -= 1
# 		grad_struct(tt,lidx,I)
# 		grad_struct(tt,ridx,I)
# 	elseif(ntype == TYPE_OU)
# 		idx -= 1 #skip oc
# 		lidx = tt[idx]
# 		idx -= 1
# 		grad_struct(tt,lidx,I)
# 	else
# 		assert(false)
# 	end
# end


#Interface function
function grad_structure{I}(tape::Tape{I}, iset::Set{I}) #non repeat version
	ilist = Array{I,1}()
	sizehint!(ilist, tape.nvnode)
	grad_struct(tape,ilist)

	empty!(iset)
	@simd for i in ilist
		push!(iset,i)
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
	forward_pass(tape,vvals,pvals,imm)
	reverse_pass(tape,imm,g)
end

function grad_reverse{V,I}(tape::Tape{I},vvals::Array{V,1},pvals::Array{V,1}, g::Array{V,1})  #dense version
	grad = Array{Tuple{I,V},1}()
	sizehint!(grad,tape.nvnode)
	grad_reverse(tape,vvals,pvals,grad)
	
	fill!(g,zero(V))
	for (i,v) in grad
		g[i] += v
	end
end