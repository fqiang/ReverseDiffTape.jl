
macro timing(cond,code)
    return quote
        if $(esc(cond))
            $(esc(code))
       end
    end
end

#the AD types below
const TYPE_V = 1    #variable node
const TYPE_P = 2    #param node
const TYPE_O = 3

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

    #bh
    bh_type::I # 0 - two vector - 1 dict
    bhi::Vector{Vector{I}}
    bhv::Vector{Vector{V}}
    bh_idxes::Vector{I}   #current horizontal indicies
    bh::Vector{Dict{I,V}}
    
    imm::Vector{V}
    imm1ordlen::I
    imm2ordlen::I

    tr::Vector{I}

    nvar::I
    nvnode::I
    nnode::I
    maxoperands::I
    fstkmax::I

    # timing stat
    enable_timing_stats::Bool
    ep_reverse_times::Vector{V}
    ep_forward_time::V
    ep_recovery_time::V
    ep_structure::V
    ep_structure_recovery::V
    ep_findnz::V
    ep_n::I


    function Tape(;imm=Vector{V}(), tt=Vector{I}(), with_timing=false, bh_type=1)
        tape = new(
            tt,  #tt
            Vector{V}(),  #stack - used for adjoints in reverse sweep
            
            Vector{I}(),  #grad_I
            Vector{V}(),  #grad value
            -one(I),      #grad indicator
            
            Vector{I}(),  #hess_I
            Vector{I}(),  #hess_J
            Vector{V}(),  #hess value
            -one(I),      #hess indicator

            bh_type, #bh_type
            Vector{Vector{I}}(),
            Vector{Vector{V}}(),
            Vector{I}(),    #current horizontal indicies
            Vector{Dict{I,V}}(),
            
            
            imm, #imm , using for both 1st and 2nd order
            zero(I),     #1st order length
            zero(I),     #2nd order length
            
            Vector{I}(), #reverse order level trace, tr vector
           
            zero(I),zero(I),zero(I),zero(I),zero(I)

            ,with_timing,Vector{V}(), zero(V), zero(V), zero(V), zero(V), zero(V), zero(I) #timing stat
            )
        finalizer(tape, tape_cleanup)
        return tape
    end
end

function tape_cleanup{I,V}(tape::Tape{I,V})
    #output timing 
    @timing tape.enable_timing_stats begin
        f = open("ep_reverse_times.txt","w")
        writedlm(f,tape.ep_reverse_times)
        close(f)

        f = open("ep_forward_time.txt","w")
        writedlm(f, tape.ep_forward_time)
        close(f)

        f = open("ep_recovery_time.txt","w")
        writedlm(f, tape.ep_recovery_time)
        close(f)
        
        f = open("ep_n.txt","w")
        writedlm(f,tape.ep_n)
        close(f)

        f = open("ep_structure.txt","w")
        writedlm(f, tape.ep_structure)
        close(f)

        f = open("ep_structure_recovery.txt","w")
        writedlm(f, tape.ep_structure_recovery)
        close(f)

        f = open("tt.txt","w")
        writedlm(f, tape.tt)
        close(f)

        f = open("tr.txt","w")
        writedlm(f, tape.tr)
        close(f)

        f = open("ep_findnz.txt","w")
        writedlm(f, tape.ep_findnz)
        close(f)

        if tape.bh_type == 1
            f = open("ep_bhi.txt","w")
            v = Vector{I}(length(tape.bhi))
            for i =1:length(tape.bhi)
                v[i] = length(tape.bhi[i])
                for j = 1:length(tape.bhi[i])
                    write(f, "$(tape.bhi[i][j]) ")
                end
                write(f,"\n")
            end
            close(f)

            f = open("ep_num_edges.txt","w")
            writedlm(f, v)
            close(f)
        end
    end
end

function tape_analysize{I,V}(tape::Tape{I,V})
    tt = tape.tt
    idx = one(I)
    istk = Vector{I}()
    v_idx_max = zero(I)
    immlen_2nd = zero(I)  # for sure >= 2ord 1ord
   	immlen_1st = zero(I)
    
    pnode_num = 0
    tnode_num = 0
    onode_num = 0
    vnode_num = 0
    @inbounds while idx <= length(tt)
        @inbounds ntype = tt[idx]
        idx += 1
        if ntype == TYPE_P
            pnode_num += 1
            idx += 3
        elseif ntype == TYPE_O
            onode_num += 1
            idx += 4
        elseif ntype == TYPE_V
            vnode_num += 1
            v_idx_max = max(v_idx_max,tt[idx])
            idx += 2
        else
            @assert false
        end
        tnode_num += 1
    end
    @assert idx == length(tt) + 1
    # @show vnode_num, pnode_num, onode_num, tnode_num
    @assert vnode_num + pnode_num + onode_num == tnode_num
    tape.nvnode = vnode_num
    tape.nvar = v_idx_max
    tape.nnode = tnode_num
    bhlen = pnode_num + onode_num + v_idx_max

    #updating ID
    node_id = tape.nvar+1
    idx = one(I)
    @inbounds while(idx <= length(tt))
        # @show idx
        @inbounds ntype = tt[idx]
        idx += 1
        
        if(ntype == TYPE_P)
            @assert tt[idx] == -1
            tt[idx] = node_id 
            node_id += 1
            push!(istk, tt[idx])
            idx += 3 #skip TYPE_P
            # @show "TYPE_P:", pnum, node_idx
        elseif(ntype == TYPE_V)
            push!(istk, tt[idx])
            idx += 2 #skip TYPE_V
        elseif(ntype == TYPE_O)
            @assert tt[idx] == -1
            tt[idx] = node_id 
            node_id += 1
            idx += 1  #skip id
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
                cid = pop!(istk)
                push!(t,cid)
            end
            append!(tape.tr, reverse!(t))
            push!(istk,tt[idx - 4])
        else
            @assert false
        end
    end
    
    # bh
    # @show tape.bh_type
    if tape.bh_type == 1
        tape.bhi = Vector{Vector{Int}}(bhlen)
        tape.bhv = Vector{Vector{Float64}}(bhlen)
        for i=1:bhlen
            @inbounds tape.bhi[i] = Vector{Int}()
            @inbounds tape.bhv[i] = Vector{Float64}()
        end
        tape.bh_idxes = zeros(Int,bhlen)  
    elseif tape.bh_type == 2 
        tape.bh = Vector{Dict{Int,Float64}}(bhlen);
        for i = 1:bhlen
            @inbounds tape.bh[i] = Dict{Int,Float64}();
        end
    else
        @assert false tape.bh_type
    end

    #timing vector for ep reverse
    @timing tape.enable_timing_stats tape.ep_reverse_times = zeros(Float64, bhlen)

    # init
	immlen = max(immlen_1st,immlen_2nd)
    if length(tape.imm) < immlen
        resize!(tape.imm, immlen)
        @show "imm size $immlen"
    end
    resize!(tape.g_I, tape.nvnode)
    resize!(tape.g, tape.nvnode)
    resize!(tape.stk, tape.nnode) 
    
    # verification
    assert(length(tape.tr) == tape.nnode-1)  #root node is not on tr
end


immutable AD{I}
    data::Vector{I}
end

function AD_V{V}(vvals::Vector{V}, val::V) #provide variable value
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
    push!(this.data,-1)
    push!(this.data,idx)
    push!(this.data,TYPE_V)
    return this
end
function AD_P{V}(pvals::Vector{V},val)
    push!(pvals,val)
    I = typeof(length(pvals))
    this = AD{I}(Array{I,1}())
    push!(this.data,TYPE_P)
    push!(this.data,-1)
    push!(this.data,length(pvals))
    push!(this.data,TYPE_P)
    return this
end

function AD_O{I}(s::Symbol,l::AD{I})
    this = AD{I}(Array{I,1}())
    append!(this.data,l.data)
    push!(this.data,TYPE_O)
    push!(this.data,-1)
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
    this = AD{I}(Vector{I}())
    @simd for i = 1:1:N
        append!(this.data,args[i].data)
    end
    push!(this.data,TYPE_O)
    push!(this.data,-1)
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

function tapeBuilder{I,V}(expr::Expr,tape::Tape{I,V}, pvals::Vector{V})
    assert(length(tape.tt)==0)
    istk = Vector{I}()
    tape_builder(expr,tape, pvals, istk)
    @assert length(istk) == 1
    tape_analysize(tape)
end

function tape_builder{I,V}(expr,tape::Tape{I,V}, pvals::Vector{V},istk::Vector{I})
    tt = tape.tt
    if isa(expr, Real)
        push!(tt,TYPE_P)
        push!(tt,-1)
        push!(tt,length(pvals)+1)
        push!(tt,TYPE_P)

        push!(istk,length(tt)-3)
        push!(pvals,expr)
    elseif isa(expr, Expr) && expr.head == :ref  #a JuMP variable
        # @show "TYPE_V", expr, tape.nextid
        assert(length(expr.args) == 2)
        vidx = expr.args[2]
        push!(tt,TYPE_V)
        push!(tt,vidx)
        push!(tt,TYPE_V)
        push!(istk,length(tt)-2)
    elseif isa(expr, Expr) && expr.head == :call
        # @show expr.args[2]
        # @show "TYPE_O", expr, tape.nextid
        op = expr.args[1]
        n = length(expr.args)-1
        assert(typeof(op)==Symbol)
        if(op==:+ && n==1) 
            #simpliy the expression eliminate 1-ary + node
            tape_builder(expr.args[2],tape,pvals,istk)
        elseif (op == :+ && n==0)
            #adding 0.0 to tape
            tape_builder(0.0,tape,pvals,istk)
        elseif (op == :- && n==1)
            tape_builder(expr.args[2],tape,pvals,istk)

            push!(tt,TYPE_O)
            push!(tt,-1)
            push!(tt,S_TO_OC[op])
            push!(tt,n)
            push!(tt,TYPE_O)
            
            tape.fstkmax = max(tape.fstkmax,length(istk))
            cidx = pop!(istk)
            push!(istk,length(tt)-4)
            # @show length(tt)
        else
            for i in 2:length(expr.args)
                tape_builder(expr.args[i],tape,pvals,istk)
            end
            push!(tt,TYPE_O)
            push!(tt,-1)
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
            push!(istk,length(tt)-4)
        end
    else
        @show "error !"
        dump(expr)
        assert(false)
    end
    nothing
end

##########################################################################################
#
# tape builder from types
#
##########################################################################################
function tapeBuilder{I}(data::Array{I,1}, t::I=1)
    tape = Tape{I,Float64}(tt=data, bh_type=t)
    tape_analysize(tape)
    return tape
end

##########################################################################################
function tape_report_mem(tape)
    report_tape_mem(tape,tape.h_type)
end

function tape_report_mem{I,V}(tape::Tape{I,V},ep::I)
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
