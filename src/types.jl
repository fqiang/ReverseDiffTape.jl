
#type alias
typealias OP_TYPE Int


#the AD types below
const TYPE_V = 1	#variable node
const TYPE_P = 2	#param node
const TYPE_O = 3


##############################################################################
# 
# MyArray type encapsulate a julia array 1-dim
#
##############################################################################
importall Base

type MyArray{V}
	a::Array{V,1}
	len::Int
	maxlen::Int
end

call{V}(::Type{MyArray{V}},max_sz::Int) = MyArray{V}(Array{V,1}(max_sz),0,max_sz)

function push!{V}(a::MyArray{V},v::V)
	# println("push!, a, v")
	# @show v
	# @show a
	# assert(a.len+1<=a.maxlen)
	a.len += 1
	@inbounds a.a[a.len] = v
	# @show a
end

function resize!{V}(a::MyArray{V},sz::Int)
	a.len = sz
end

function length(a::MyArray)
	return a.len
end

function getindex{V}(a::MyArray{V},i::Int)
	# assert(i<=a.maxlen)
	return a.a[i]
end

function setindex!{V}(a::MyArray{V},v::V,i::Int)
	# assert(i<=a.maxlen)
	@inbounds a.a[i] = v
end

endof(a::MyArray) = length(a)


function Base.show(io::IO,m::MyArray)
	println(io,m.a[1:m.len],",",m.len,",",m.maxlen)
end

##############################################################################

type Tape{I<:Int}
	tt::Array{I,1}
	nvar::I
	nvnode::I
	nnode::I
	maxoperands::I
	imm2ord::I
	
	function Tape()
		return new(Array{I,1}(),zero(I),zero(I),zero(I),zero(I),zero(I))
	end

	function Tape(data::Array{I,1})
		this = new(data,zero(I),zero(I),zero(I),zero(I),zero(I))
		analysize_tape(this)
		return this
	end
end

function analysize_tape{I}(tape::Tape{I})
	tt = tape.tt
	idx = one(I)
	iset = Set{I}()
	@inbounds while(idx <= length(tt))
		# @show idx
		ntype = tt[idx]
		idx += 1
		if(ntype == TYPE_P)
			idx += 2 #skip TYPE_P
		elseif(ntype == TYPE_V)
			push!(iset,tt[idx])
			idx += 2 #skip TYPE_V
			tape.nvnode += 1
		elseif(ntype == TYPE_O)
			idx += 1  #skip oc
			n = tt[idx]
			idx += 2  #skip TYPE_O
			tape.imm2ord = n + round(I,(n+1)*n/2)
			tape.maxoperands<n?tape.maxoperands=n:nothing
		end
		tape.nnode += 1
	end
	tape.nvar += length(iset)
end

immutable AD{I}
	data::Array{I,1}
end

function AD_V{V}(vvals::Array{V,1}, val) #provide variable value
	push!(vvals,val)
	I = typeof(length(vvals))
	this = AD{I}(Array{I,1}())
	push!(this.data,TYPE_V)
	push!(this.data,length(vvals))
	push!(this.data,TYPE_V)
	return this
end
function AD_V{I}(idx::I) #without variable value
	this = AD{I}(Array{I,1}())
	push!(this.data,TYPE_V)
	push!(this.data,idx)
	push!(this.data,TYPE_V)
	return this
end
function AD_P{V}(pvals::Array{V,1},val)
	push!(pvals,val)
	I = typeof(length(pvals))
	this = AD{I}(Array{I,1}())
	push!(this.data,TYPE_P)
	push!(this.data,length(pvals))
	push!(this.data,TYPE_P)
	return this
end

function AD_O{I}(s::Symbol,l::AD{I})
	@show s
	this = AD{I}(Array{I,1}())
	append!(this.data,l.data)
	push!(this.data,TYPE_O)
	push!(this.data,S_TO_OC[s])
	push!(this.data,1) #1 operand simply 
	push!(this.data,TYPE_O)
end
function AD_O{I,N}(s::Symbol,args::NTuple{N,AD{I}})
	# @show s
	# @show args
	# @show N
	# @show n
	assert(N>1)
	this = AD{I}(Array{I,1}())
	@simd for i = 1:1:N
		append!(this.data,args[i].data)
	end
	push!(this.data,TYPE_O)
	push!(this.data,S_TO_OC[s])
	push!(this.data,N) 
	push!(this.data,TYPE_O)
	return this
end

function Base.show(io::IO,m::AD)
	print(io, m.data)
end
