
macro timing(cond,code)
    return quote
        if $(esc(cond))
            $(esc(code))
       end
    end
end

macro asserting(cond, code)
    return quote
        if $(esc(cond))
            @assert $(esc(code))
        end
    end
end

#the AD types below
const TYPE_P = 0   #param node
const TYPE_V = 1    #variable node
const TYPE_O = 2
const TYPE_O1 = 3  
const TYPE_O2 = 4  
const TYPE_O3 = 5
const TYPE_O4 = 6
const TYPE_O5 = 7
const TYPE_O6 = 8
const TYPE_O7 = 9

type Tape{I,V}
    tt::Vector{I}
    tr::Vector{I}
    nvar::I
    nvnode::I
    nnode::I
    maxoperands::I
    depth::I
    immlen::I

    nzg::I
    nzh::I

    #bh
    bh_type::I # 0 - two vector - 1 dict
    bhi::Vector{Vector{I}}
    bhv::Vector{Vector{V}}
    bh_idxes::Vector{I}   #current horizontal indicies
    
    stk::Vector{V}
    vals::Vector{V}
    imm::Vector{V}

    # timing stat
    debug::Bool
    enable_timing_stats::Bool
    ep_reverse_times::Vector{V}
    ep_forward_time::V
    ep_recovery_time::V
    ep_structure::V
    ep_structure_recovery::V
    ep_findnz::V
    ep_n::I


    function Tape(stk::Vector{V}, vals::Vector{V}, imm::Vector{V}; tt=Vector{I}(), with_timing=false, bh_type=1, with_debug=true)
        tape = new(
            tt,  #tt
            Vector{I}(), #reverse order level trace, tr vector
            zero(I),zero(I),zero(I),zero(I),zero(I),zero(I),
            
            -one(I),      #grad indicator
            -one(I),      #hess indicator

            bh_type, #bh_type
            Vector{Vector{I}}(),
            Vector{Vector{V}}(),
            Vector{I}(),    #current horizontal indicies
            
            stk, 
            vals,
            imm,
            
            with_debug, with_timing,Vector{V}(), zero(V), zero(V), zero(V), zero(V), zero(V), zero(I) #timing stat
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

function tape_analysize{I,V}(tape::Tape{I,V}, gnvar::I, with_hess_mem::Bool)
    tt = tape.tt
    tr = tape.tr
    idx = one(I)
    mxn = one(I)
    mxn_times = one(I)
    vmax = zero(I)
    vnodes = zero(I)
    nnodes = zero(I)
    n = [0,0,0,0,0,0,0,0,0]

    @inbounds while idx <= length(tt)
        @inbounds ntype = tt[idx]
        # @show ntype
        nnodes += 1
        if ntype == TYPE_P
            @assert false
        elseif ntype == TYPE_V
            n[TYPE_V] += 1
            vmax = max(vmax, tt[idx+1])
            vnodes += 1
            idx += 3
        elseif ntype == TYPE_O
            n[TYPE_O] +=1
            @assert tt[idx+3] != 0
            mxn = max(mxn,tt[idx+3])
            if OP[tt[idx+2]]==:* 
                mxn_times = max(mxn_times,tt[idx+3])
            end
            idx += 5
        elseif ntype == TYPE_O1
            n[TYPE_O1] += 1
            idx += 5
        elseif ntype == TYPE_O2
            n[TYPE_O2] += 1
            idx += 5 
        elseif ntype == TYPE_O3
            n[TYPE_O3] += 1
            vmax = max(vmax,tt[idx+4])
            vnodes += 1
            idx += 6

        elseif ntype == TYPE_O4
            n[TYPE_O4] += 1
            vmax = max(vmax,tt[idx+4])
            vnodes += 1
            idx += 6

        elseif ntype == TYPE_O5
            n[TYPE_O5] += 1
            vmax = max(vmax,tt[idx+3])
            vnodes += 1
            idx += 5            

        elseif ntype == TYPE_O6
            n[TYPE_O6] += 1
            vmax = max(vmax,tt[idx+3])
            vnodes += 1
            idx += 5
        elseif ntype == TYPE_O7
            n[TYPE_O7] += 1
            vmax = max(vmax,tt[idx+3])
            vnodes += 1
            idx += 5
        else
            @assert false
        end
        # @show vmax
    end
    tape.maxoperands = mxn
    tape.nvar = gnvar
    tape.nvnode = vnodes
    tape.nnode = nnodes
    
    @assert gnvar >= vmax
    @assert idx == length(tt) + 1
    # @show n
    
    #updating ID
    nid = gnvar
    idx = one(I)
    istk = Vector{I}()
    @inbounds while(idx <= length(tt))
        @inbounds ntype = tt[idx]
        # @show istk
        # @show ntype
        if ntype == TYPE_P
            @assert false

        elseif ntype == TYPE_V
            push!(istk, tt[idx+1])
            idx += 3

        elseif ntype == TYPE_O
            nid += 1
            tt[idx+1] = nid
            num = tt[idx+3]
            istklen = length(istk)
            append!(tape.tr,sub(istk,istklen-num+1:istklen))
            resize!(istk,istklen-num)
            push!(istk, nid)
            idx += 5

        elseif ntype == TYPE_O1
            nid += 1
            tt[idx+1] = nid
            push!(tape.tr,pop!(istk))
            push!(istk,nid)
            idx += 5

        elseif ntype == TYPE_O2
            nid += 1
            tt[idx+1] = nid
            push!(tape.tr,pop!(istk))
            push!(istk,nid)
            idx += 5

        elseif ntype == TYPE_O3
            nid += 1
            tt[idx+1] = nid
            push!(istk,nid)
            idx += 6

        elseif ntype == TYPE_O4
            nid += 1
            tt[idx+1] = nid
            push!(istk,nid)
            idx += 6

        elseif ntype == TYPE_O5
            nid += 1
            tt[idx+1] = nid
            push!(tape.tr,pop!(istk))
            push!(istk,nid)
            idx += 5            

        elseif ntype == TYPE_O6
            nid += 1
            tt[idx+1] = nid
            push!(tape.tr,pop!(istk))
            push!(istk,nid)
            idx += 5

        elseif ntype == TYPE_O7
            nid += 1
            tt[idx+1] = nid
            push!(istk,nid)
            idx += 5
        else
            @assert false
        end
    end
    tape.immlen = max(5, convert(I,(mxn_times-1)*mxn_times/2))

    if with_hess_mem 
        @assert tape.bh_type == 1
        bhlen = nid
        tape.bhi = Vector{Vector{Int}}(bhlen)
        tape.bhv = Vector{Vector{Float64}}(bhlen)
        for i=1:bhlen
            @inbounds tape.bhi[i] = Vector{Int}()
            @inbounds tape.bhv[i] = Vector{Float64}()
        end
        tape.bh_idxes = zeros(Int,bhlen)  
        #timing vector for ep reverse
        @timing tape.enable_timing_stats tape.ep_reverse_times = zeros(Float64, bhlen)
    end
    
    # verification
    @assert length(tape.tr) == tape.nnode-1  #root node is not on tr
end

function resizeWorkingMemory{I,V}(obj_tt::Tape{I,V},lag_tt::Tape{I,V}, cons_tt::Vector{Tape{I,V}})
    mx_vallen = max(obj_tt.nnode,lag_tt.nnode)
    mx_depth  = max(obj_tt.depth,lag_tt.depth)
    mx_ops    = max(obj_tt.maxoperands,lag_tt.maxoperands)
    mx_immlen = max(obj_tt.immlen,lag_tt.immlen)
    # @show mx_vallen, mx_depth, mx_ops, mx_immlen, length(obj_tt.vals), length(obj_tt.stk), length(obj_tt.imm)

    for i=1:length(cons_tt)
        mx_immlen = max(mx_immlen,cons_tt[i].immlen)
        mx_vallen = max(mx_vallen,cons_tt[i].nnode)
        mx_depth = max(mx_depth,cons_tt[i].depth)
        mx_ops = max(mx_ops,cons_tt[i].maxoperands)
    end
    resize!(obj_tt.stk, mx_depth + mx_ops)
    resize!(obj_tt.vals, mx_vallen)
    resize!(obj_tt.imm, mx_immlen)
    # @show mx_vallen, mx_depth, mx_ops, mx_immlen, length(obj_tt.vals), length(obj_tt.stk), length(obj_tt.imm)

    for i = 1:length(cons_tt)
        @assert length(cons_tt[i].imm) == length(obj_tt.imm) == length(lag_tt.imm)
        @assert length(cons_tt[i].stk) == length(obj_tt.stk) == length(lag_tt.stk) 
        @assert length(cons_tt[i].vals) == length(obj_tt.vals) == length(lag_tt.vals)
    end
end

##########################################################################################
#
# tape builder from a Julia Expr type
#   tape memory property is initialized 
#
##########################################################################################

function tapeBuilderSimple{I,V}(expr::Expr,tape::Tape{I,V}, pvals::Vector{V}, gnvar::I)
    @assert length(tape.tt)==0
    tape_builder(expr,tape, pvals)
    tape_analysize(tape, gnvar, false)
end

function tapeBuilder{I,V}(expr::Expr,tape::Tape{I,V}, pvals::Vector{V}, gnvar::I)
    @assert length(tape.tt)==0
    tape_builder(expr,tape, pvals)
    tape_analysize(tape, gnvar, false)
end

function mergeTapes{I,V}(tape::Tape{I,V}, obj_tape::Tape{I,V}, tapes::Vector{Tape{I,V}}, pvals::Vector{V}, gnvar::I)
    @assert length(tape.tt) == 0
    tt = tape.tt
    mxd = zero(I)
    
    #objective
    mxd = max(obj_tape.depth,mxd)
    append!(tt, obj_tape.tt)
    push!(pvals,1.0)

    push!(tt,TYPE_O1)
    push!(tt,-1)
    push!(tt,S_TO_OC[:*])
    push!(tt,length(pvals))
    push!(tt,TYPE_O1)

    #constraints
    for i=1:length(tapes)
        mxd = max(tapes[i].depth,mxd)
        append!(tt, tapes[i].tt)
        push!(pvals,1.0)

        push!(tt,TYPE_O1)
        push!(tt,-1)
        push!(tt,S_TO_OC[:*])
        push!(tt,length(pvals))
        push!(tt,TYPE_O1)
    end

    #adding all
    if length(tapes) > 0
        push!(tt, TYPE_O)
        push!(tt,-1)
        push!(tt,S_TO_OC[:+])
        push!(tt,length(tapes)+1)
        push!(tt, TYPE_O)
    end

    tape.depth = mxd
    tape_analysize(tape, gnvar, true)
end

type B_NODE{I,V}
    t::I
    vidx::I
    pval::V
end

function tape_builder{I,V}(expr,tape::Tape{I,V}, pvals::Vector{V}, d::I=0)
    # @show expr
    # @show d
    tt = tape.tt
    d += 1
    tape.depth = max(tape.depth,d)

    if isa(expr, Real)
        return B_NODE(0,0,expr)
    elseif isa(expr, Expr) && expr.head == :ref  #a JuMP variable
        @assert length(expr.args) == 2
        vidx = expr.args[2]
        return B_NODE(1,vidx,0.0)
    elseif isa(expr, Expr) && expr.head == :call
        op = expr.args[1]
        n = length(expr.args)-1
        @assert typeof(op)==Symbol
        if (n==1 || n==0) && op==:+ 
            #simpliy the expression eliminate 1-ary + node
            return tape_builder(expr.args[2],tape,pvals,d)
        elseif op == :- && n==1
            node = tape_builder(expr.args[2],tape,pvals,d)

            if node.t == 0
                return B_NODE(0,0,-node.pval)
            elseif node.t == 1
                push!(tt,TYPE_O7)
                push!(tt,-1)
                push!(tt.S_TO_OC[op])
                push!(tt,node.vidx)
                push!(tt,TYPE_O7)
                return B_NODE(9,0,0.0)
            elseif node.t == 2 || node.t == 3 || node.t == 4 || node.t==5 || node.t == 6 || node.t == 7 || node.t == 8 || node.t == 9
                push!(tt, TYPE_O)
                push!(tt,-1)
                push!(tt,S_TO_OC[op])
                push!(tt,n)
                push!(tt,TYPE_O)
                return B_NODE(2,0,0.0)
            else
                @assert false
            end
        elseif n == 1
            node = tape_builder(expr.args[2],tape,pvals,d)
            if node.t == 0
                pval = eval(op)(node.pval)
                return B_NODE(0,0,pval)
            elseif node.t == 1
                push!(tt, TYPE_O7)
                push!(tt, -1)
                push!(tt, S_TO_OC[op])
                push!(tt, node.vidx)
                push!(tt, TYPE_O7)
                return B_NODE(9,0,0.0)
            elseif node.t == 2 || node.t == 3 || node.t == 4 || node.t==5 || node.t == 6 || node.t == 7 || node.t == 8 || node.t == 9
                push!(tt, TYPE_O)
                push!(tt, -1)
                push!(tt, S_TO_OC[op])
                push!(tt, n)
                push!(tt, TYPE_O)
                return B_NODE(2,0,0.0)
            else 
                @assert false
            end
        elseif n == 2
            ln = tape_builder(expr.args[2],tape,pvals,d)
            rn = tape_builder(expr.args[3],tape,pvals,d)
            if ln.t == rn.t == 0
                pval = eval(op)(ln.pval,rn.pval)
                return B_NODE(0,0,pval)
            elseif ln.t == 0 && rn.t != 1
                push!(pvals, ln.pval)

                push!(tt,TYPE_O1)
                push!(tt,-1)
                push!(tt,S_TO_OC[op])
                push!(tt,length(pvals))
                push!(tt,TYPE_O1)
                return B_NODE(3,0,0.0)
            elseif ln.t == 0 && rn.t == 1
                push!(pvals, ln.pval)

                push!(tt,TYPE_O3)
                push!(tt, -1)
                push!(tt, S_TO_OC[op])
                push!(tt, length(pvals))
                push!(tt, rn.vidx)
                push!(tt, TYPE_O3)
                return B_NODE(5,0,0.0)
            elseif rn.t == 0 && ln.t != 1
                push!(pvals, rn.pval)

                push!(tt,TYPE_O2)
                push!(tt,-1)
                push!(tt,S_TO_OC[op])
                push!(tt,length(pvals))
                push!(tt,TYPE_O2)
                return B_NODE(4,0,0.0)
            elseif rn.t == 0 && ln.t == 1
                push!(pvals,rn.pval)

                push!(tt, TYPE_O4)
                push!(tt,-1)
                push!(tt,S_TO_OC[op])
                push!(tt,length(pvals))
                push!(tt,ln.vidx)
                push!(tt,TYPE_O4)
                return B_NODE(6,0,0.0)
            elseif ln.t == 1 && rn.t == 1
                push!(tt,TYPE_V)
                push!(tt,ln.vidx)
                push!(tt,TYPE_V)

                push!(tt,TYPE_V)
                push!(tt,rn.vidx)
                push!(tt,TYPE_V)

                push!(tt,TYPE_O)
                push!(tt,-1)
                push!(tt,S_TO_OC[op])
                push!(tt,n)
                push!(tt,TYPE_O)
                return B_NODE(2,0,0.0)
            elseif ln.t == 1 && rn.t != 1
                push!(tt,TYPE_O5)
                push!(tt,-1)
                push!(tt,S_TO_OC[op])
                push!(tt,ln.vidx)
                push!(tt,TYPE_O5)
                return B_NODE(7,0,0.0)
            elseif rn.t == 1 && ln.t != 1
                push!(tt,TYPE_O6)
                push!(tt, -1)
                push!(tt,S_TO_OC[op])
                push!(tt,rn.vidx)
                push!(tt,TYPE_O6)
                return B_NODE(8,0,0.0)
            else
                push!(tt,TYPE_O)
                push!(tt,-1)
                push!(tt,S_TO_OC[op])
                push!(tt,n)
                push!(tt,TYPE_O)
                return B_NODE(2,0,0.0)
            end
        else #n>2
            @assert (op == :* || op == :+) && n>2
            rstk = Vector{B_NODE}()
            pval = zero(V)
            pnode = zero(I)
            for i in 2:length(expr.args)
                node = tape_builder(expr.args[i],tape,pvals,d)
                if node.t == 0 
                    pnode==zero(I)? pval = node.pval : ( op == :+? pval += node.pval : pval *= node.pval )
                    pnode += 1
                elseif node.t == 1
                    push!(rstk,node)
                else
                    #nothing
                end
            end

            nv = length(rstk)
            on = n - pnode - nv

            if pnode == n
                return B_NODE(0,0,pval)
            elseif nv == 1
                if on == 1 && pnode == 0
                    push!(tt,TYPE_O5)
                    push!(tt,-1)
                    push!(tt,S_TO_OC[op])
                    push!(tt,rstk[1].vidx)
                    push!(tt,TYPE_O5)
                    return B_NODE(7,0,0.0)
                elseif on == 1 && pnode > 0
                    push!(tt,TYPE_O5)
                    push!(tt,-1)
                    push!(tt,S_TO_OC[op])
                    push!(tt,rstk[1].vidx)
                    push!(tt,TYPE_O5)

                    push!(pvals, pval)
                    push!(tt,TYPE_O1)
                    push!(tt,-1)
                    push!(tt,S_TO_OC[op])
                    push!(tt,length(pvals))
                    push!(tt,TYPE_O1)
                    return B_NODE(3,0,0.0)
                elseif on > 1 && pnode == 0
                    push!(tt,TYPE_O)
                    push!(tt,-1)
                    push!(tt,S_TO_OC[op])
                    push!(tt,on)
                    push!(tt,TYPE_O)

                    push!(tt,TYPE_O5)
                    push!(tt,-1)
                    push!(tt,S_TO_OC[op])
                    push!(tt,rskt[1].vidx)
                    push!(tt,TYPE_O5)
                    return B_NODE(7,0,0.0)
                elseif on == 0
                    @assert pnode > 0
                    push!(pvals,pval)
                    push!(tt,TYPE_O3)
                    push!(tt,-1)
                    push!(tt,S_TO_OC[op])
                    push!(tt,length(pvals))
                    push!(tt,rstk[1].vidx)
                    push!(tt,TYPE_O3)
                    return B_NODE(5,0,0.0)
                else
                    # @show expr, pnode, nv , on , n
                    @assert on > 1 && pnode > 1 && nv == 1
                    push!(tt,TYPE_O)
                    push!(tt,-1)
                    push!(tt,S_TO_OC[op])
                    push!(tt,on)
                    push!(tt,TYPE_O)

                    push!(tt,TYPE_O5)
                    push!(tt,-1)
                    push!(tt,S_TO_OC[op])
                    push!(tt,rskt[1].vidx)
                    push!(tt,TYPE_O5)

                    push!(pvals, pval)
                    push!(tt,TYPE_O1)
                    push!(tt,-1)
                    push!(tt,S_TO_OC[op])
                    push!(tt,length(pvals))
                    push!(tt,TYPE_O1)
                    return B_NODE(3,0,0.0)
                end
            elseif nv > 1 
                for item in rstk
                    push!(tt, TYPE_V)
                    push!(tt, item.vidx)
                    push!(tt, TYPE_V)
                end
                @assert n-pnode > 1
                push!(tt,TYPE_O)
                push!(tt,-1)
                push!(tt,S_TO_OC[op])
                push!(tt,n-pnode)
                push!(tt,TYPE_O)

                if pnode > 0
                    push!(pvals, pval)
                    push!(tt, TYPE_O1)
                    push!(tt,-1)
                    push!(tt, S_TO_OC[op])
                    push!(tt, length(pvals))
                    push!(tt, TYPE_O1)
                    return B_NODE(3,0,0)
                else
                    @assert pnode == 0
                    return B_NODE(2,0,0.0)
                end
            elseif nv == 0
                if on == 1 
                    @assert pnode > 1
                    push!(pvals, pval)
                    push!(tt, TYPE_O1)
                    push!(tt,-1)
                    push!(tt, S_TO_OC[op])
                    push!(tt, length(pvals))
                    push!(tt, TYPE_O1)
                    return B_NODE(3,0,0)
                else
                    @assert on > 1
                    @assert n-pnode > 0
                    push!(tt,TYPE_O)
                    push!(tt,-1)
                    push!(tt,S_TO_OC[op])
                    push!(tt,n-pnode)
                    push!(tt,TYPE_O)


                    if pnode > 0
                        push!(pvals, pval)
                        push!(tt, TYPE_O1)
                        push!(tt,-1)
                        push!(tt, S_TO_OC[op])
                        push!(tt, length(pvals))
                        push!(tt, TYPE_O1)
                        return B_NODE(3,0,0)
                    else
                        # @show expr, pnode, nv , on , n
                        @assert pnode == 0 && on == n
                        return B_NODE(2,0,0.0)
                    end
                end
            end
        end
    else
        @show "error !"
        dump(expr)
        assert(false)
    end
    d -= 1
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
#
# operator overload interface
#
##########################################################################################
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
# function tape_report_mem(tape)
#     report_tape_mem(tape,tape.h_type)
# end

# function tape_report_mem{I,V}(tape::Tape{I,V},ep::I)
#     si = sizeof(I)
#     sf = sizeof(V)
#     i  = 0
#     tb = 0

#     tb_now = length(tape.tt) * si
#     @show "tape ", length(tape.tt), tb_now
#     tb += tb_now

#     tb_now = length(tape.stk) *si
#     @show "stk ", length(tape.stk), tb_now
#     tb += tb_now

#     tb_now = length(tape.g_I) * si + length(tape.g)*sf + 1*si
#     tb += tb_now
    
#     tb_now = length(tape.h_I) * si + length(tape.h_J) * si + length(tape.hess)*sf+ 1*si
#     tb += tb_now
#     @show " hessian nnz", length(tape.h_I)

#     tb_now = sizeof(tape.node_idx_to_number)
#     tb_now += length(tape.node_idx_to_number.nzval)  * sf
#     tb_now += length(tape.node_idx_to_number.colptr)  * si
#     tb_now += length(tape.node_idx_to_number.rowval) * si
#     tb += tb_now

#     if ep==2
#         tb_now = sizeof(tape.bh)
#         for j = 1:length(tape.bh)
#             tb_now += sizeof(tape.bh[j])  #length(bh[j]) * 8
#             tb_now += length(tape.bh[j]) * sizeof(mPair{I,V})
#         end
#         @show "bh - bytes ", tb_now
#         tb += tb_now
#     elseif ep==3
#         tb_now = sizeof(tape.bh3)
#         for d in tape.bh3
#             tb_now += sizeof(d)
#             tb_now += length(d) * (sizeof(I) + sizeof(V))
#         end
#         @show "bh - bytes ", tb_now
#         tb += tb_now    
#     end

#     tb_now = length(tape.bh_idxes) * si
#     @show "bh_idxes ", length(tape.bh_idxes), tb_now
#     tb += tb_now

#     tb_now = length(tape.imm) * sf
#     @show "imm ", length(tape.imm), tb_now,  "needed length " , tape.imm2ordlen
#     tb += tb_now

#     tb += 2*si #imm1ordlen , imm2ordlen

#     tb_now = length(tape.tr) * si
#     @show "tr ",length(tape.tr), tb_now
#     tb += tb_now
#     tb += 1*si  #trlen

#     tb += 5* si #nvar, nvnode, nnode, maxoperands, fstkmax

#     @show " tape = ", tb ," bytes"

#     tb_ep1 = sizeof(tape.eset)
#     tb_ep1 += sizeof(tape.liveVar)
#     tb_ep1 += sizeof(tape.h)
#     @show " with ep1 data ", tb_ep1, " bytes "

#     return tb+tb_ep1
# end
