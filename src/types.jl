
#type alias
typealias VV_TYPE Float64
typealias IDX_TYPE UInt
typealias OP_TYPE UInt
typealias TT_TYPE Array{UInt,1}
typealias TV_STACK Stack
typealias TV_TYPE Array{VV_TYPE,1}


#the AD types below
const TYPE_V = 1	#variable node
const TYPE_P = 2	#param node
const TYPE_OU = 3	#unary op
const TYPE_OB = 4  	#binary op

abstract Placeholder
immutable AD_V <: Placeholder
	tt::TT_TYPE
	idx::IDX_TYPE  #index on vvals
	t::Uint    #type code
end
immutable AD_P <: Placeholder
	tt::TT_TYPE
	idx::IDX_TYPE  #index on pvals
	t::Uint    #type code
end
immutable AD_O <: Placeholder
	tt::TT_TYPE
	idx::IDX_TYPE  #index on tape
    t::Uint    #type code	 
end

const AD_TYPES = [Type{AD_V},Type{AD_P},Type{AD_O}]

function AD_V(tt::TT_TYPE, vvals::TV_TYPE,val=NaN)
	push!(vvals,val)
	this = AD_V(tt,length(vvals),TYPE_V)
	return this
end

function AD_P(tt::TT_TYPE, pvals::TV_TYPE,val=NaN)
	push!(pvals,val)
	this = AD_P(tt,length(pvals),TYPE_P)
	return this
end

function AD_O(tt::TT_TYPE, oc, lidx::UInt)
	push!(tt,lidx)
	push!(tt,oc)
	push!(tt,TYPE_OU)
	this = AD_O(tt,length(tt),TYPE_OU)
	return this
end


function AD_O(tt::TT_TYPE,oc, lidx::UInt, ridx::UInt)
	push!(tt,lidx)
	push!(tt,ridx)
	push!(tt,oc)
	push!(tt,TYPE_OB)
	this = AD_O(tt,length(tt),TYPE_OB)
	return this
end

function Base.show(io::IO,m::AD_V)
	assert(m.t == TYPE_V)
	print(io, "AD_V[",m.idx,"]")
end

function Base.show(io::IO,m::AD_P)
	assert(m.t == TYPE_P)
	print(io, "AD_P[",m.idx,"]")
end

function Base.show(io::IO,m::AD_O)
	assert(m.t == TYPE_OU || m.t == TYPE_OB)
	print(io, "AD_O[",m.idx,"]")
end
