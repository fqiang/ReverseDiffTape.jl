
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

type Tape{I<:Int}
	tt::Array{I,1}
	nvar::I
	nvnode::I
	nnode::I
	maxoperands::I
	
	function Tape()
		return new(Array{I,1}(),zero(I),zero(I),zero(I),zero(I))
	end
end


abstract Placeholder{I,V}
immutable AD_V{I,V} <: Placeholder{I,V}
	tape::Tape{I}
	idx::I  #index on vvals
	t::I    #type code

	function AD_V(tape::Tape{I}, vvals::Array{V,1}, val::V) #provide variable value
		push!(vvals,val)
		push!(tape.tt,TYPE_V)
		push!(tape.tt,length(vvals))
		push!(tape.tt,TYPE_V)
		tape.nvar += 1
		tape.nnode += 1
		return new(tape,length(vvals),TYPE_V)
	end
	function AD_V(tape::Tape{I},idx::I) #without variable value
		push!(tape.tt,TYPE_V)
		push!(tape.tt,idx)
		push!(tape.tt,TYPE_V)
		tape.nar += 1
		tape.nnode += 1
		return new(tape,idx,TYPE_V)
	end
end

immutable AD_P{I,V} <: Placeholder{I,V}
	tape::Tape{I}
	idx::I  #index on pvals
	t::I    #type code

	function AD_P(tape::Tape{I}, pvals::Array{V,1},val::V)
		push!(pvals,val)
		push!(tape.tt,TYPE_P)
		push!(tape.tt,length(pvals))
		push!(tape.tt,TYPE_P)
		tape.nnode += 1
		new(tape,length(pvals),TYPE_P)
	end
end

immutable AD_O{I,V} <: Placeholder{I,V}
	tape::Tape{I}
	idx::I  #index on tape
    t::I    #type code	 

    function AD_O(tape::Tape{I},s::Symbol,args...)
    	# @show args
    	n = length(args)
    	# @show n
    	assert(n>=1)
    	push!(tape.tt,TYPE_O)
    	push!(tape.tt,S_TO_OC[s])
    	push!(tape.tt,n) 
    	push!(tape.tt,TYPE_O)
    	tape.maxoperands < n?tape.maxoperands=n:nothing
    	tape.nnode += 1
    	new(tape,length(tape.tt),TYPE_O)
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
	assert(m.t == TYPE_O )
	print(io, "AD_O[",m.idx,"]")
end
