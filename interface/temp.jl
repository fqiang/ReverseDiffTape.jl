 
function bar1(tt::Array{Int,1})
	tt_cp = Array{Int,1}()
	for i in tt
       push!(tt_cp,i)
    end
    return tt_cp
end


function bar2(tt::Array{Int,1})
	tt_cp = Array{Int,1}()
	i=1
	while i<=length(tt)
       push!(tt_cp,tt[i])
       i += 1
    end
    return tt_cp
end


function forward_eval(tt::Array{Int,1})
	vals = Array{Float64,1}()
	idx = 1
	while idx<=length(tt)
		# @show idx
		ntype = tt[idx]
		idx += 1
		if(ntype == 1)
			# v = vvals[tt[idx]]
			v = 1.0
			idx += 1
			push!(vals,v)
			idx += 1 #skip TYPE_V
		elseif(ntype == 2)
			# v = pvals[tt[idx]]
			v = 1.0
			idx += 1
			push!(vals,v)
			idx += 1 #skip TYPE_P
		else #(ntype == 3)
			# nvals = Array{Float64,1}()
			oc = tt[idx]
			idx += 1
			n = tt[idx]
			idx += 1
			idx += 1 #skip TYPE_O
			# s = OP[oc]
			# while n>0
			# 	prepend!(nvals,[pop!(vals)])
			# 	n -= 1
			# end
			v = 1
			# # @show OP[oc]
			# # @show nvals

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
			# # @show v
			push!(vals,v)
		# else
		# 	assert(false)
		end
	end
	# @show vals
	return pop!(vals)
end


function forward_eval1(tt::Array{Int,1}, vals::Array{Float64,1}, vvals::Array{Float64,1}, pvals::Array{Float64,1})
	# vals = Array{VV_TYPE,1}()
	idx = 1
	while i<=length(tt)
		# @show idx
		ntype = tt[idx]
		idx += 1
		if(ntype == 2)
			v = pvals[tt[idx]]
			idx += 1
			push!(vals,v)
			idx += 1 #skip TYPE_P
		elseif(ntype == 1)
			v = vvals[tt[idx]]
			idx += 1
			push!(vals,v)
			idx += 1 #skip TYPE_V
		elseif(ntype == 3)
			nvals = Array{Float64,1}()
			oc = tt[idx]
			idx += 1
			n = tt[idx]
			idx += 1
			idx += 1 #skip TYPE_O
			# s = OP[oc]
			while n>0
				prepend!(nvals,[pop!(vals)])
				n -= 1
			end
			
			# @show OP[oc]
			# @show nvals

			if(OP[oc] == :+)
				v = Base.sum(nvals)
			elseif(OP[oc] == :-)
				v = (-)(nvals[1],nvals[2])
			elseif(OP[oc] == :*)
				v = Base.prod(nvals)
			elseif(OP[oc] == :/)
				v =  (/)(nvals[1],nvals[2])
			elseif(OP[oc] == :^)
				v =  (^)(nvals[1],nvals[2])
			elseif(OP[oc] == :sin)
				v =  sin(nvals[1])
			elseif(OP[oc] == :cos)
				v =  cos(nvals[1])
			end
			# @show v
			push!(vals,v)
		else
			assert(false)
		end
	end
	# @show vals
	return pop!(vals)
end