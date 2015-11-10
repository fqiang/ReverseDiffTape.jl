
## forward evaluation 
function evaluate(tt::TT_TYPE,idx::IDX_TYPE, vvals::TV_TYPE, pvals::TV_TYPE)
	# debug("enter - ",idx)
	# assert(idx != 0)	
	ret::VV_TYPE = 0.0
	ntype = tt[idx] #type node
	idx -= 1
	if(ntype == TYPE_P)
		ret = pvals[tt[idx]]
		idx -=1
	elseif(ntype == TYPE_V)
		ret = vvals[tt[idx]]
		idx -= 1
	# elseif(ntype == TYPE_OU)
	# 	oc = tt[idx]
	# 	assert(U_OP_START <= oc <= U_OP_END)
	# 	idx -= 1
	# 	lidx = tt[idx]
	# 	idx -= 1
	# 	# debug("unary before  - ",lidx)
	# 	(lval,idx) = evaluate(tt,lidx,vvals,pvals)
	# 	ret = evaluate(OP[oc],lval)
	# elseif(ntype == TYPE_OB)
	# 	oc = tt[idx]
	# 	assert(B_OP_START<= oc <= B_OP_END)
	# 	idx -= 1
	# 	ridx = tt[idx]
	# 	idx -= 1
	# 	lidx = tt[idx]
	# 	idx -= 1
	# 	# debug("before right - ",ridx)
	# 	(rval,idx) = evaluate(tt,ridx,vvals,pvals)
	# 	# debug("before left - ",lidx)
	# 	(lval,idx) = evaluate(tt,lidx,vvals,pvals)
	# 	ret = evaluate(OP[oc],lval,rval)
	else
		num = ntype - 2
		oc = tt[idx]
		idx -= 1
		# assert(oc==1 || oc ==3) #:+ or :* symbol
		ovals = Array{VV_TYPE,1}()
		nidx = 0
		for i = 1:1:num
			oidx = tt[idx]
			idx -= 1
			(rval,nidx) = evaluate(tt,oidx,vvals,pvals)
			push!(ovals,rval)
		end
		idx = nidx
		# @show OP[oc]
		ret = evaluate(OP[oc],ovals...)
		# assert(false)
	end

	# assert(!isnan(ret))
	# debug("exit - ",idx)
	return ret,idx
end


# function evaluate(s::Symbol, nvals::Array{VV_TYPE,1}, i)
# 	# n = length(nvals)
# 	if(s == :+)
# 		result = 0.0
# 		@simd for j=i:1:length(nvals)
# 			@inbounds result+=nvals[j]
# 		end
# 		return result
# 		# return Base.sum(nvals[i:1:end])
# 	elseif(s == :-)
# 		@inbounds return (-)(nvals[i],nvals[i+1])
# 	elseif(s == :*)
# 		result = nvals[i]
# 		@simd for j=i+1:1:length(nvals)
# 			@inbounds result *= nvals[j]
# 		end
# 		@inbounds return result
# 		# return Base.prod(nvals[i:1:end])
# 	elseif(s == :/)
# 		@inbounds return (/)(nvals[i],nvals[i+1])
# 	elseif(s == :^)
# 		@inbounds return (^)(nvals[i],nvals[i+1])
# 	elseif(s == :sin)
# 		@inbounds return sin(nvals[i])
# 	elseif(s == :cos)
# 		@inbounds return cos(nvals[i])
# 	end
# end

# function sum_n(x::Array{VV_TYPE,1},n)
# 	result = 0.0
# 	for i in 1:1:n
# 		result += x[i]
# 	end
# 	return result
# end


# function prod_n(x::Array{VV_TYPE,1},n)
# 	result = 0.0
# 	for i in 1:1:n
# 		result *= x[i]
# 	end
# 	return result
# end

# @eval 	function evaluate(oc, args...)
# 			s = OP[oc]
# 			out = Expr(:if)
# 			comp = Expr(:comparison)

# 			if(s == :+)
# 				return (+)(args...)
# 			elseif(s == :-)
# 				return (-)(args...)
# 			elseif(s == :*)
# 				return (*)(args...)
# 			elseif(s == :/)
# 				return (/)(args...)
# 			elseif(s == :^)
# 				return (^)(args...)
# 			elseif(s == :sin)
# 				return sin(args...)
# 			elseif(s == :cos)
# 				return cos(args...)
# 			end
# 		end


# function evaluate(s::Symbol, nvals)
# 	if(s == :+)
# 		return Base.sum(nvals)
# 	elseif(s == :-)
# 		return (-)(nvals[1],nvals[2])
# 	elseif(s == :*)
# 		return Base.prod(nvals)
# 	elseif(s == :/)
# 		return (/)(nvals[1],nvals[2])
# 	elseif(s == :^)
# 		return (^)(nvals[1],nvals[2])
# 	elseif(s == :sin)
# 		return sin(nvals[1])
# 	elseif(s == :cos)
# 		return cos(nvals[1])
# 	end
# end


function forward_evaluate(tt::TT_TYPE, vvals::TV_TYPE, pvals::TV_TYPE)
	# vals = Array{VV_TYPE,1}()
	idx = 1::Int
	vals = Array{Float64,1}()
	v = Array{Float64,1}(1)
	sizehint!(v,1)
	sizehint!(vals,10010)

	while(idx <= length(tt))
		# @show idx
		@inbounds ntype = tt[idx]
		idx += 1
		v[1] = 0.0
		if(ntype == TYPE_P)
			@inbounds v[1] = pvals[tt[idx]]
			idx += 1
			@inbounds push!(vals,v[1])
			idx += 1 #skip TYPE_P
		elseif(ntype == TYPE_V)
			@inbounds v[1] = vvals[tt[idx]]
			idx += 1
			@inbounds push!(vals,v[1])
			idx += 1 #skip TYPE_V
		elseif(ntype == TYPE_O)
			@inbounds oc = tt[idx]
			idx += 1
			@inbounds n = tt[idx]
			idx += 1
			idx += 1 #skip TYPE_O
			# s = OP[oc]
		 	# v = evaluate(OP[oc],vals[length(vals)-n+1:1:length(vals)])::Float64
		 	@inbounds evaluate(OP[oc],vals,length(vals)-n+1,v)

			# @show OP[oc]
			# @show nvals

			# if(OP[oc] == :+)
			# 	v = Base.sum(nvals)
			# elseif(OP[oc] == :-)
			# 	v = (-)(nvals[1],nvals[2])
			# elseif(OP[oc] == :*)
			# 	v = Base.prod(nvals)
			# elseif(OP[oc] == :/)
			# 	v =  (/)(nvals[1],nvals[2])
			# elseif(OP[oc] == :^)
			# 	v =  (^)(nvals[1],nvals[2])
			# elseif(OP[oc] == :sin)
			# 	v =  sin(nvals[1])
			# elseif(OP[oc] == :cos)
			# 	v =  cos(nvals[1])
			# end

			@inbounds resize!(vals,length(vals)-n+1)
			# @show v
			@inbounds vals[end] = v[1]
		end
	end
	# @show vals
	return vals[1]
end



function evaluate(s::Symbol, nvals::Array{VV_TYPE,1}, i, result::Array{VV_TYPE,1})
	# n = length(nvals)
	if(s == :+)
		@inbounds result[1] = 0.0
		@simd for j=i:1:length(nvals)
			@inbounds result[1]+=nvals[j]
		end
	elseif(s == :-)
		@inbounds return (-)(nvals[i],nvals[i+1])
	elseif(s == :*)
		@inbounds result[1] = nvals[i]
		@simd for j=i+1:1:length(nvals)
			@inbounds result[1] *= nvals[j]
		end
	elseif(s == :/)
		@inbounds result[1]=(/)(nvals[i],nvals[i+1])
	elseif(s == :^)
		@inbounds result[1]=(^)(nvals[i],nvals[i+1])
	elseif(s == :sin)
		@inbounds result[1]=sin(nvals[i])
	elseif(s == :cos)
		@inbounds result[1]=cos(nvals[i])
	end
end


## Interface method
function feval(tt::TT_TYPE, vvals::TV_TYPE, pvals::TV_TYPE)
	# (val::Float64,idx) = evaluate(tt,length(tt),vvals, pvals)
	# nvals = Array{Float64,1}(1000000)
	# nvals = Array{Float64,1}(1000000)
	val = forward_evaluate(tt,vvals,pvals)
	return val
end
