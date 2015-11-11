
# forward pass on the tape tt, to build ss stack

#forward pass on the scalar function
function forward_pass{V,I}(tt::Array{I,1}, vvals::Array{V,1}, pvals::Array{V,1}, imm::Array{V,1})
	idx = one(I)
	v = [zero(V)]
	stk = Array{V,1}() #used for value evaluation
	empty!(imm)  #used for immediate derrivative
	sizehint!(imm,round(I,length(tt)/3.5))
	sizehint!(stk,round(I,length(tt)/10))
	sizehint!(v,1)

	@inbounds while(idx <= length(tt))
		# @show idx
		println("++++++++++++++++++++++++++++++++++++")
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
			@show OP[oc], length(stk)-n+1, n
			eval_idd(OP[oc],stk,length(stk)-n+1,imm,v)
			@show imm
			@show v
			resize!(stk,length(stk)-n+1)
			stk[end] = v[1]
		end
		println("++++++++++++++++++++++++++++++++++++")
	end
	@show stk
end


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


function reverse_pass{V,I}(tt::Array{I,1},imm::Array{V,1},g::Array{Tuple{I,V},1})
	idx = length(tt)
	adjs = Array{V,1}()
	sizehint!(adjs, round(I,length(tt)/10))  #random guessing size
	push!(adjs,one(V))  #init value 

	while(idx > 0)
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

function grad_nz(tt::TT_TYPE, idx::IDX_TYPE, I::Set{IDX_TYPE}) #non-repeat version
	assert(length(vvals) == length(grad))
	# debug("enter - ", adj)
	ntype = tt[idx]
	idx -= 1
	if(ntype ==TYPE_P)  
		#do nothing , this is a parameter
	elseif(ntype == TYPE_V)
		vidx = tt[idx]
		push!(I,vidx)
	elseif(ntype == TYPE_OB)
		idx -= 1
		ridx = tt[idx]
		idx -= 1
		lidx = tt[idx]
		idx -= 1
		nzg(tt,lidx,I)
		nzg(tt,ridx,I)
	elseif(ntype == TYPE_OU)
		idx -= 1 #skip oc
		lidx = tt[idx]
		idx -= 1
		nzg(tt,lidx,I)
	else
		assert(false)
	end
end

function grad_struct(tt::TT_TYPE, idx::IDX_TYPE, I::Array{IDX_TYPE,1}) #repeat version
	# debug("enter - ", adj)
	ntype = tt[idx]
	idx -= 1
	if(ntype ==TYPE_P)  
		#do nothing , this is a parameter
	elseif(ntype == TYPE_V)
		vidx = tt[idx]
		push!(I,vidx) 
	elseif(ntype == TYPE_OB)
		idx -= 1
		ridx = tt[idx]
		idx -= 1
		lidx = tt[idx]
		idx -= 1
		grad_struct(tt,lidx,I)
		grad_struct(tt,ridx,I)
	elseif(ntype == TYPE_OU)
		idx -= 1 #skip oc
		lidx = tt[idx]
		idx -= 1
		grad_struct(tt,lidx,I)
	else
		assert(false)
	end
end


#Interface function
function grad_structure(tt::TT_TYPE, I::Set{IDX_TYPE}) #non repeat version
	grad_nz(tt,convert(IDX_TYPE,length(tt)),I)
end

function grad_structure(tt::TT_TYPE, I::Array{IDX_TYPE,1})  #repeat version
	grad_struct(tt,convert(IDX_TYPE,length(tt)),I)
end

function grad_reverse(tt::TT_TYPE,vvals::TV_TYPE,pvals::TV_TYPE, grad::Array{Tuple{IDX_TYPE,VV_TYPE},1}) #sparse version
	ss = TV_STACK(VV_TYPE)
	reverse_grad_0(tt,convert(IDX_TYPE,length(tt)),ss,vvals,pvals)
	adj = 1.0
	reverse_grad_1(tt,convert(IDX_TYPE,length(tt)),ss,adj,vvals,pvals,grad)
	# @show grad
end

function grad_reverse(tt::TT_TYPE,vvals::TV_TYPE,pvals::TV_TYPE, grad::TV_TYPE)  #dense version
	fill!(grad,0.0)
	ss = TV_STACK(VV_TYPE)
	reverse_grad_0(tt,convert(IDX_TYPE,length(tt)),ss,vvals,pvals)
	adj = 1.0
	g = Array{Tuple{IDX_TYPE,VV_TYPE},1}()
	reverse_grad_1(tt,convert(IDX_TYPE,length(tt)),ss,adj,vvals,pvals,g)
	for (i,v) in g
		grad[i] = v
	end
end