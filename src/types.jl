
#type alias
typealias VV_TYPE Float64
typealias IDX_TYPE Int
typealias OP_TYPE Int
typealias TT_TYPE Array{IDX_TYPE,1}
typealias TV_STACK Stack
typealias TV_TYPE Array{VV_TYPE,1}


#the AD types below
const TYPE_V = 1	#variable node
const TYPE_P = 2	#param node
const TYPE_O = 3
const TYPE_OU = 3	#unary op
const TYPE_OB = 4  	#binary op

abstract Placeholder
immutable AD_V <: Placeholder
	tt::TT_TYPE
	idx::IDX_TYPE  #index on vvals
	t::Uint    #type code

	function AD_V(tt::TT_TYPE, vvals::TV_TYPE,val) #provide variable value
		push!(vvals,val)
		new(tt,length(vvals),TYPE_V)
	end
	function AD_V(tt::TT_TYPE,idx) #without variable value
		new(tt,idx,TYPE_V)
	end
end

immutable AD_P <: Placeholder
	tt::TT_TYPE
	idx::IDX_TYPE  #index on pvals
	t::Uint    #type code

	function AD_P(tt::TT_TYPE, pvals::TV_TYPE,val)
		push!(pvals,val)
		new(tt,length(pvals),TYPE_P)
	end
end

immutable AD_O <: Placeholder
	tt::TT_TYPE
	idx::IDX_TYPE  #index on tape
    t::Uint    #type code	 

  #   function AD_O(tt,oc,lidx,ridx)
  #   	push!(tt,lidx)
		# push!(tt,ridx)
		# push!(tt,oc)
		# push!(tt,TYPE_OB)
  #   	new(tt,length(tt),TYPE_OB)
  #   end
  #   function AD_O(tt,oc,lidx)
  #   	push!(tt,lidx)
		# push!(tt,oc)
		# push!(tt,TYPE_OU)
  #   	new(tt,length(tt),TYPE_OB)
  #   end
    function AD_O(tt,oc,args...)
    	# @show args
    	num = length(args)
    	type_c = num + 2
    	# @show num
    	assert(num>=1)
    	for i in args
    		push!(tt,i)
    	end
    	push!(tt,oc)
    	push!(tt,type_c) #based on TYPE_OU + 2
    	new(tt,length(tt),type_c)
    end
end

const AD_TYPES = [Type{AD_V},Type{AD_P},Type{AD_O}]

# function AD_V(tt::TT_TYPE, vvals::TV_TYPE,val=NaN)
# 	push!(vvals,val)
# 	this = AD_V(tt,length(vvals),TYPE_V)
# 	return this
# end

# function AD_P(tt::TT_TYPE, pvals::TV_TYPE,val=NaN)
# 	push!(pvals,val)
# 	this = AD_P(tt,length(pvals),TYPE_P)
# 	return this
# end

# function AD_O(tt::TT_TYPE, oc, lidx::IDX_TYPE)
# 	push!(tt,lidx)
# 	push!(tt,oc)
# 	push!(tt,TYPE_OU)
# 	this = AD_O(tt,length(tt),TYPE_OU)
# 	return this
# end


# function AD_O(tt::TT_TYPE,oc, lidx::IDX_TYPE, ridx::IDX_TYPE)
# 	push!(tt,lidx)
# 	push!(tt,ridx)
# 	push!(tt,oc)
# 	push!(tt,TYPE_OB)
# 	this = AD_O(tt,length(tt),TYPE_OB)
# 	return this
# end

function Base.show(io::IO,m::AD_V)
	assert(m.t == TYPE_V)
	print(io, "AD_V[",m.idx,"]")
end

function Base.show(io::IO,m::AD_P)
	assert(m.t == TYPE_P)
	print(io, "AD_P[",m.idx,"]")
end

function Base.show(io::IO,m::AD_O)
	assert(m.t >= TYPE_OU )
	print(io, "AD_O[",m.idx,"]")
end
