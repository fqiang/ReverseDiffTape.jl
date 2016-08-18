
#type alias
typealias OP_TYPE Int
typealias IDX_TYPE Int
typealias EdgeSet{I,V} Dict{I,Dict{I,V}}

#the AD types below
const TYPE_V = 1    #variable node
const TYPE_P = 2    #param node
const TYPE_O = 3

##############################################################################
# 
# MyArray type encapsulate a julia array 1-dim
#
##############################################################################
# importall Base

# type MyArray{V}
#     a::Array{V,1}
#     len::Int
#     maxlen::Int
# end

# call{V}(::Type{MyArray{V}},max_sz::Int) = MyArray{V}(Array{V,1}(max_sz),0,max_sz)

# @inline function push!{V}(a::MyArray{V},v::V)
#     # println("push!, a, v")
#     # @show v
#     # @show a
#     # assert(a.len+1<=a.maxlen)
#     a.len += 1
#     @inbounds a.a[a.len] = v
#     # @show a
# end

# @inline function resize!{V}(a::MyArray{V},sz::Int)
#     assert(sz<=a.maxlen)
#     a.len = sz
# end

# @inline function length(a::MyArray)
#     return a.len
# end

# @inline function getindex{V}(a::MyArray{V},i::Int)
#     # assert(i<=a.maxlen)
#     @inbounds return a.a[i]
# end

# @inline function setindex!{V}(a::MyArray{V},v::V,i::Int)
#     # assert(i<=a.maxlen)
#     @inbounds a.a[i] = v
# end

# endof(a::MyArray) = length(a)


# function Base.show(io::IO,m::MyArray)
#     println(io,m.a[1:m.len],",",m.len,",",m.maxlen)
# end

##############################################################################

type mPair{I,V}
    i::I
    w::V
    function mPair()
        return new(zero(I),0.0)  #not valid entry
    end
    function mPair(idx,ww)
        return new(idx,ww)
    end
end

type EP4{I,V}
    s::Set{I}
    v::SparseMatrixCSC{V,I}
    function EP4()
        return new(Set{I}(),sparsevec([1],[0.0]))
    end
    function EP4(ep4::EP4)
        set = ep4.s
        idxes = Vector{I}(length(set))
        vals = Vector{V}(length(set))
        i::I = 1
        for idx in set
            idxes[i] = idx
            vals[i] = zero(V)
            i += 1
        end
        if(!isempty(idxes))
            ep4.v = sparsevec(idxes, vals)
        end
        return ep4
    end
end

type Tape{I,V}
    tt::Vector{I}
    stk::Vector{V}

    g_I::Vector{I}
    g::Vector{V}
    nzg::I
    
    h_I::Vector{I}
    h_J::Vector{I}
    hess::Vector{V}
    nzh::I

    #use by ep2
    node_idx_to_number::SparseMatrixCSC{I,I}
    bh::Vector{Vector{mPair{I,V}}}  #big hessian matrix for everybody
    bh_idxes::Vector{I}   #current horizontal indicies
    # bh_length::Vector{I}  

    #used by ep3
    bh3::Vector{Dict{I,V}} 

    #used by ep4
    bh4::Vector{EP4{I,V}}

    h_type::I

    imm::Vector{V}
    imm1ordlen::I
    imm2ordlen::I

    tr::Vector{I}
    trlen::I

    eset::Dict{I,Dict{I,V}}
    liveVar::Dict{I,Set{I}}
    h::Dict{I,Dict{I,V}}


    nvar::I
    nvnode::I
    nnode::I
    maxoperands::I
    fstkmax::I

    t_push_edge::V
    
    function Tape()
        return new(
            Vector{I}(),  #tt
            Vector{V}(),  #stack - used for adjoints in reverse sweep
            
            Vector{I}(),  #grad_I
            Vector{V}(),  #grad value
            -one(I),      #grad indicator
            
            Vector{I}(),  #hess_I
            Vector{I}(),  #hess_J
            Vector{V}(),  #hess value
            -one(I),      #hess indicator
            sparsevec([1],[-1]),  #node_idx_to_number


            Vector{Vector{mPair{I,V}}}(),  #big hessian matrix
            Vector{I}(),    #current horizontal indicies
            # Vector{I}(),    #horizontal lengths

            Vector{Dict{I,V}}(), #bh3

            Vector{EP4{I,V}}(), #bh4

            zero(I), #h_type

            Vector{V}(), #imm , using for both 1st and 2nd order
            zero(I),     #1st order length
            zero(I),     #2nd order length
            
            Vector{I}(), #reverse order level trace, tr vector
            zero(I),     #trlen


            Dict{I,Dict{I,V}}(),  #eset
            Dict{I,Set{I}}(),     #liveVar
            Dict{I,Dict{I,V}}(),  #h

            zero(I),zero(I),zero(I),zero(I),zero(I)
            ,zero(V)
            )
    end

    function Tape(data::Vector{I})
        this = new(
            data,         #tt
            Vector{V}(),  #stack
            
            Vector{I}(),  #grad_I
            Vector{V}(),  #grad value
            -one(I),      #grad indicator
            
            Vector{I}(),  #hess_I
            Vector{I}(),  #hess_J
            Vector{V}(),  #hess value
            -one(I),      #hess indicator
            sparsevec([1],[-1]),  #node_idx_to_number
            Vector{Vector{mPair{I,V}}}(),  #big hessian matrix
            Vector{I}(),    #current horizontal indicies

            Vector{Dict{I,V}}(), #bh3

            Vector{EP4{I,V}}(), #bh4

            zero(I), #h_type

            Vector{V}(), #imm , using for both 1st and 2nd order
            zero(I),     #1st order length
            zero(I),     #2nd order length
            
            Vector{I}(), #reverse order level trace, tr vector
            zero(I),     #trlen


            Dict{I,Dict{I,V}}(),  #eset
            Dict{I,Set{I}}(),     #liveVar
            Dict{I,Dict{I,V}}(),  #h

            zero(I),zero(I),zero(I),zero(I),zero(I)
            )
        analysize_tape(this)
        return this
    end
end

function analysize_tape{I,V}(tape::Tape{I,V})
    tt = tape.tt
    idx = one(I)
    istk = Vector{I}()
    v_idx_max = zero(I)
    immlen_2nd = zero(I)  # for sure >= 2ord 1ord
   	immlen_1st = zero(I)
 
    node_idxes = Vector{I}()
    node_numbers = Vector{I}()
    @inbounds while(idx <= length(tt))
        # @show idx
        @inbounds ntype = tt[idx]
        idx += 1
        if(ntype == TYPE_P)
            idx += 2 #skip TYPE_P
            push!(istk,idx-3)
            node_idx = idx - 3
        elseif(ntype == TYPE_V)
            v_idx_max = max(v_idx_max,tt[idx])
            idx += 2 #skip TYPE_V
            tape.nvnode += 1
            push!(istk,idx-3)
            node_idx = idx - 3
        elseif(ntype == TYPE_O)
			@inbounds op_sym = OP[tt[idx]]
            idx += 1  #skip oc
            @inbounds n = tt[idx]
            idx += 2  #skip TYPE_O
			if n==1
				immlen_2nd += 2 
			else
				if op_sym==:*
            		immlen_2nd += div(n*(n+1),2)  #max estimation
				elseif op_sym == :+ || op_sym ==:-	
					immlen_2nd += 0
				else
					immlen_2nd += 5
				end
			end
			immlen_1st += n			

            tape.maxoperands = max(n,tape.maxoperands)
            
            tape.fstkmax = max(tape.fstkmax,length(istk))
            t = Vector{I}() #slow but works
            for i = 1:n
                cidx = pop!(istk)
                push!(t,cidx)
            end
            append!(tape.tr, reverse!(t))
            push!(istk,idx-4)
            node_idx = idx - 4
        end
        tape.nnode += 1
        push!(node_numbers, tape.nnode) #node id number start with 1
        push!(node_idxes,node_idx)      #on tape index
    end
    tape.nvar = v_idx_max
    node_idxes = node_idxes + tape.nvar #shift nodes up to make room for independent nodes
    node_numbers = node_numbers + tape.nvar #shift  up 
    
    # @show node_idxes,node_numbers
    prepend!(node_numbers,collect(1:tape.nvar))    #top tape.nvar is the independent nodes
    prepend!(node_idxes,collect(1:tape.nvar))
    # @show node_idxes,node_numbers
    tape.node_idx_to_number = sparsevec(node_idxes,node_numbers)  #independent nodes mapping
    
    # used by ep2
    tape.bh = Vector{Vector{mPair{Int,Float64}}}(tape.nnode+tape.nvar)
    # tape.bh_length = round(Int,readdlm("log4000_1_bh_length.txt")[:,1])
    for i=1:tape.nnode+tape.nvar 
        tape.bh[i] = Vector{mPair{Int,Float64}}()
        # for j=1:tape.bh_length[i]
        #     push!(tape.bh[i],mPair{Int,Float64}())
        # end
    end
    # fill!(tape.bh_length,0)
    # tape.bh_length = zeros(Int,tape.nnode+tape.nvar)


    tape.bh_idxes = zeros(Int,tape.nnode+tape.nvar)
    


    # used by ep3
    tape.bh3 = Vector{Dict{I,V}}(tape.nnode + tape.nvar);
    for i=1:tape.nnode + tape.nvar
        @inbounds tape.bh3[i] = Dict{I,V}();
    end

    #used by ep4
    tape.bh4 = Vector{EP4{I,V}}(tape.nnode + tape.nvar);
    for i = 1:tape.nnode + tape.nvar
        @inbounds tape.bh4[i] = EP4{I,V}();
    end

    # init
	# @show max(immlen_1st,immlen_2nd)
    resize!(tape.imm, max(immlen_1st,immlen_2nd))
    resize!(tape.g_I, tape.nvnode)
    resize!(tape.g, tape.nvnode)
    resize!(tape.stk, tape.nnode) 
    tape.trlen = length(tape.tr)
    # tape.nzg = -1
    # verification
    assert(length(tape.tr) == tape.nnode-1)  #root node is not on tr
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
#   tape memory property is initialized 
#
##########################################################################################

#building tape with Julia expression
function tapeBuilder{I,V}(expr::Expr, pvals::Vector{V}, i::I=0)
    tape = Tape{I,V}()
    tapeBuilder(expr,tape,pvals)
    return tape
end

function tapeBuilder{I,V}(expr::Expr,tape::Tape{I,V}, pvals::Vector{V})
    # @show expr
    assert(length(tape.tt)==0)
    istk = Vector{I}()
    tapeBuilder(expr,tape, pvals, istk)
    analysize_tape(tape)
end

function tapeBuilder{I,V}(expr::Expr,tape::Tape{I,V}, pvals::Vector{V},istk::Vector{I})
    tt = tape.tt
    head = expr.head
    if(head == :ref)  #a JuMP variable
        assert(length(expr.args) == 2)
        vidx = expr.args[2]
        push!(tt,TYPE_V)
        push!(tt,vidx)
        push!(tt,TYPE_V)
        
        push!(istk,length(tt)-2)
        # push!(vset,vidx)
        # tape.nvnode += 1
        # tape.nnode += 1
    elseif(head == :call)
        # @show expr.args[2]
        op = expr.args[1]
        n = length(expr.args)-1
        assert(typeof(op)==Symbol)
        if(op==:+ && n==1) 
            #simpliy the expression eliminate 1-ary + node
            tapeBuilder(expr.args[2],tape,pvals,istk)
        elseif (op == :+ && n==0)
            #adding 0.0 to tape
            tapeBuilder(0.0,tape,pvals,istk)
        elseif (op == :- && n==1)
            tapeBuilder(expr.args[2],tape,pvals,istk)

            push!(tt,TYPE_O)
            push!(tt,S_TO_OC[op])
            push!(tt,n)
            push!(tt,TYPE_O)
            
            tape.fstkmax = max(tape.fstkmax,length(istk))
            cidx = pop!(istk)
            push!(istk,length(tt)-3)
        else
            for i in 2:length(expr.args)
                tapeBuilder(expr.args[i],tape,pvals,istk)
            end
            push!(tt,TYPE_O)
            push!(tt,S_TO_OC[op])
            push!(tt,n)
            push!(tt,TYPE_O)

            tape.fstkmax = max(tape.fstkmax,length(istk))
            # t = Vector{I}()
            for i=1:n
                cidx = pop!(istk)
                # push!(t,cidx)
            end
            # append!(tape.tr,reverse!(t))
            push!(istk,length(tt)-3)
            # tape.nnode += 1
            # tape.maxoperands = max(length(expr.args)-1, tape.maxoperands)
            # tape.imm2ordlen += n + round(I,n*(n+1)/2)
            # tape.eset[length(tt)-3] = Dict{I,V}()
            # @show length(tt) - 3
        end
    else
        @show "error !"
        dump(expr)
        assert(false)
    end
    nothing
end

function tapeBuilder{I,V}(expr::Real, tape::Tape{I,V}, pvals::Vector{V},istk::Vector{I}) #a JuMP parameter
    tt = tape.tt
    push!(tt,TYPE_P)
    push!(tt,length(pvals)+1)
    push!(tt,TYPE_P)

    push!(istk,length(tt)-2)
    push!(pvals,expr)
    # tape.nnode += 1
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
function report_tape_mem(tape)
    report_tape_mem(tape,tape.h_type)
end

function report_tape_mem{I,V}(tape::Tape{I,V},ep::I)
    si = sizeof(I)
    sf = sizeof(V)
    i  = 0
    tb = 0

    tb_now = length(tape.tt) * si
    @show "tape ", length(tape.tt), tb_now
    tb += tb_now

    tb_now = length(tape.stk) *si
    @show "stk ", length(tape.stk), tb_now
    tb += tb_now

    tb_now = length(tape.g_I) * si + length(tape.g)*sf + 1*si
    tb += tb_now
    
    tb_now = length(tape.h_I) * si + length(tape.h_J) * si + length(tape.hess)*sf+ 1*si
    tb += tb_now
    @show " hessian nnz", length(tape.h_I)

    tb_now = sizeof(tape.node_idx_to_number)
    tb_now += length(tape.node_idx_to_number.nzval)  * sf
    tb_now += length(tape.node_idx_to_number.colptr)  * si
    tb_now += length(tape.node_idx_to_number.rowval) * si
    tb += tb_now

    if ep==2
        tb_now = sizeof(tape.bh)
        for j = 1:length(tape.bh)
            tb_now += sizeof(tape.bh[j])  #length(bh[j]) * 8
            tb_now += length(tape.bh[j]) * sizeof(mPair{I,V})
        end
        @show "bh - bytes ", tb_now
        tb += tb_now
    elseif ep==3
        tb_now = sizeof(tape.bh3)
        for d in tape.bh3
            tb_now += sizeof(d)
            tb_now += length(d) * (sizeof(I) + sizeof(V))
        end
        @show "bh - bytes ", tb_now
        tb += tb_now    
    end

    tb_now = length(tape.bh_idxes) * si
    @show "bh_idxes ", length(tape.bh_idxes), tb_now
    tb += tb_now

    tb_now = length(tape.imm) * sf
    @show "imm ", length(tape.imm), tb_now,  "needed length " , tape.imm2ordlen
    tb += tb_now

    tb += 2*si #imm1ordlen , imm2ordlen

    tb_now = length(tape.tr) * si
    @show "tr ",length(tape.tr), tb_now
    tb += tb_now
    tb += 1*si  #trlen

    tb += 5* si #nvar, nvnode, nnode, maxoperands, fstkmax

    @show " tape = ", tb ," bytes"

    tb_ep1 = sizeof(tape.eset)
    tb_ep1 += sizeof(tape.liveVar)
    tb_ep1 += sizeof(tape.h)
    @show " with ep1 data ", tb_ep1, " bytes "

    return tb+tb_ep1
end
