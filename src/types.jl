
#type alias
typealias OP_TYPE Int
typealias IDX_TYPE Int
typealias EdgeSet{I,V} Dict{I,Dict{I,V}}

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

type Tape{I,V}
	tt::Vector{I}
	stk::Vector{V}
	g_I::Vector{I}
	g::Vector{V}
	imm1ord::Vector{V}
	imm1ordlen::I

	tr::Vector{I}
	trlen::I

	eset::Dict{I,Dict{I,V}}
	liveVar::Dict{I,Set{I}}
	imm2ord::Vector{V}
	imm2ordlen::I
	h::Dict{I,Dict{I,V}}

	nvar::I
	nvnode::I
	nnode::I
	maxoperands::I
	fstkmax::I
	nzg::I

	function Tape()
		return new(Vector{I}(),
			Vector{V}(),
			Vector{I}(),
			Vector{V}(),
			Vector{V}(),zero(I),
			Vector{I}(),zero(I),
			Dict{I,Dict{I,V}}(),Dict{I,Set{I}}(),
			Vector{V}(),zero(I),
			Dict{I,Dict{I,V}}(),
			zero(I),zero(I),zero(I),zero(I),zero(I),-one(I))
	end

	function Tape(data::Vector{I})
		this = new(data,
			Vector{V}(),
			Vector{I}(),
			Vector{V}(),
			Vector{V}(),zero(I),
			Vector{I}(),zero(I),
			Dict{I,Dict{I,V}}(),Dict{I,Set{I}}(),
			Vector{V}(),zero(I),
			Dict{I,Dict{I,V}}(),
			zero(I),zero(I),zero(I),zero(I),zero(I),-one(I))
		analysize_tape(this)
		return this
	end
end

function analysize_tape{I,V}(tape::Tape{I,V})
	tt = tape.tt
	idx = one(I)
	istk = Vector{I}()
	iset = Set{I}()
	@inbounds while(idx <= length(tt))
		# @show idx
		ntype = tt[idx]
		idx += 1
		if(ntype == TYPE_P)
			idx += 2 #skip TYPE_P
			push!(istk,idx-3)
		elseif(ntype == TYPE_V)
			push!(iset,tt[idx])
			idx += 2 #skip TYPE_V
			tape.nvnode += 1
			push!(istk,idx-3)
		elseif(ntype == TYPE_O)
			idx += 1  #skip oc
			n = tt[idx]
			idx += 2  #skip TYPE_O
			tape.imm2ordlen += n + round(I,(n+1)*n/2)  #max estimation
			tape.maxoperands = max(n,tape.maxoperands)
			
			tape.fstkmax = max(tape.fstkmax,length(istk))
			t = Vector{I}() #slow but works
			for i = 1:n
				cidx = pop!(istk)
				push!(t,cidx)
			end
			append!(tape.tr, reverse!(t))
			push!(istk,idx-4)
		end
		tape.nnode += 1
	end
	tape.nvar = length(iset)

	# init
	resize!(tape.imm1ord, tape.nnode-1)
	tape.imm1ordlen = length(tape.imm1ord)
	resize!(tape.imm2ord, tape.imm2ordlen)
	resize!(tape.g_I, tape.nvnode)
	resize!(tape.g, tape.nvnode)
	resize!(tape.stk, tape.nnode) 
	tape.trlen = length(tape.tr)
	tape.nzg = -1
	# verification
	assert(length(tape.imm1ord) == tape.nnode-1)
	assert(length(tape.tr) == tape.nnode-1)	
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
	this = AD{I}(Array{I,1}())
	append!(this.data,l.data)
	push!(this.data,TYPE_O)
	push!(this.data,S_TO_OC[s])
	push!(this.data,1) #1 operand simply 
	push!(this.data,TYPE_O)
	return this
end

function AD_O{I}(s::Symbol,l::AD{I},r::AD{I})
	AD_O(s,tuple(l,r))
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


##########################################################################################
#
# tape builder from a Julia Expr type
#	tape memory property is initialized 
#
##########################################################################################

#building tape with Julia expression
function tapeBuilder{I,V}(expr::Expr,tape::Tape{I,V}, pvals::Array{V,1})
	# @show expr
	vset = Set{I}()
	istk = Vector{I}()
	tapeBuilder(expr,tape, pvals, vset, istk)
	assert(length(tape.tr)==tape.nnode-1)
	
	tape.nvar = length(vset)
	tape.imm1ordlen = tape.nnode -1
	resize!(tape.imm1ord, tape.imm1ordlen)
	resize!(tape.imm2ord, tape.imm2ordlen)

	resize!(tape.stk, tape.nnode)
	resize!(tape.g_I, tape.nvnode)
	resize!(tape.g,tape.nvnode)
	tape.nzg = -1  #-1 until call grad_structure
	
	# @show tape
end

function tapeBuilder{I,V}(expr::Expr,tape::Tape{I,V}, pvals::Array{V,1},vset::Set{I},istk::Vector{I})
	tt = tape.tt
	head = expr.head
	if(head == :ref)  #a JuMP variable
		assert(length(expr.args) == 2)
		vidx = expr.args[2]
		push!(tt,TYPE_V)
		push!(tt,vidx)
		push!(tt,TYPE_V)
		
		push!(istk,length(tt)-2)
		push!(vset,vidx)
		tape.nvnode += 1
		tape.nnode += 1
	elseif(head == :call)
		# @show expr.args[2]
		op = expr.args[1]
		n = length(expr.args)-1
		assert(typeof(op)==Symbol)
		if(op==:+ && n<2) 
			#simpliy the expression eliminate 1-ary + node
			tapeBuilder(expr.args[2],tape,pvals,vset,istk)
		else
			for i in 2:length(expr.args)
				tapeBuilder(expr.args[i],tape,pvals,vset,istk)
			end
			push!(tt,TYPE_O)
			push!(tt,S_TO_OC[op])
			push!(tt,n)
			push!(tt,TYPE_O)

			tape.fstkmax = max(tape.fstkmax,length(istk))
			t = Vector{I}()
			for i=1:n
				cidx = pop!(istk)
				push!(t,cidx)
			end
			append!(tape.tr,reverse!(t))
			push!(istk,length(tt)-3)
			tape.nnode += 1
			tape.maxoperands = max(length(expr.args)-1, tape.maxoperands)
			tape.imm2ordlen += n + round(I,n*(n+1)/2)
			# tape.eset[length(tt)-3] = Dict{I,V}()
			# @show length(tt) - 3
			if op==:+ && n<2
				@show expr
				assert(false)
			end
		end
    else
    	println("error !")
    	dump(expr)
    	assert(false)
    end
    nothing
end

function tapeBuilder{I,V}(expr::Real, tape::Tape{I,V}, pvals::Array{V,1},vset::Set{I},istk::Vector{I}) #a JuMP parameter
	tt = tape.tt
	push!(tt,TYPE_P)
	push!(tt,length(pvals)+1)
	push!(tt,TYPE_P)

	push!(istk,length(tt)-2)
	push!(pvals,expr)
	tape.nnode += 1
end

##########################################################################################
#
# tape builder from types
#
##########################################################################################
function tapeBuilder{I}(data::Array{I,1})
	tape = Tape{I,Float64}(data)
	return tape
end

##########################################################################################

