
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

type HessStorage{I,V}
    bh_type::I
    idxmap::Vector{I}
    HH::Vector{V}
    gnvar::I
    bhi::Vector{Vector{I}}
    bhv::Vector{Vector{V}}
    bh_idxes::Vector{I}
end

function call{I,V}(::Type{HessStorage{I,V}}, bhlen::I, gnvar::I;bh_type::I=1) 
    hs = HessStorage{I,V}(bh_type, Vector{V}(), Vector{I}(), gnvar, Vector{Vector{I}}(bhlen), Vector{Vector{V}}(bhlen), Vector{I}(bhlen))
    for i=1:bhlen
        @inbounds hs.bhi[i] = Vector{I}()
        @inbounds hs.bhv[i] = Vector{V}()
    end
    fill!(hs.bh_idxes,0)
    return hs
end

type Tape{I,V}
    tt::Vector{I}
    tr::Vector{I}
    stklen::I
    nvar::I
    nvnode::I
    nnode::I
    immlen::I
    depth::I
    
    # nzg::I

    # timing stat
    debug::Bool
    enable_timing_stats::Bool
    ep_forward_time::V
    ep_recovery_time::V
    ep_structure::V
    ep_structure_recovery::V
    ep_findnz::V
    ep_n::I

    function Tape(tt::Vector{I}, tr::Vector{I}, stklen::I, gnvar::I, vnodes::I, nnodes::I, immlen::I, depth::I; with_timing=false, with_debug=false)
        tape = new(
            tt, tr, stklen, gnvar, vnodes, nnodes, immlen, depth,
            
            # -one(I),      #grad indicator
             
            with_debug, with_timing, zero(V), zero(V), zero(V), zero(V), zero(V), zero(I) #timing stat
            )
        finalizer(tape, tape_cleanup)
        return tape
    end

    function Tape()
        return new(Vector{I}(),Vector{I}(),0,0,0,0,0,0,
            false,false, zero(V), zero(V), zero(V), zero(V), zero(V), zero(I))
    end
end

function tape_analysize{I}(tt::Vector{I}, tr::Vector{I}, gnvar::I)
    @assert length(tr) == 0
    idx = one(I)
    mx_stklen = one(I)
    stklen = zero(I)
    mxn_times = one(I)
    vmax = zero(I)
    vnodes = zero(I)
    nnodes = zero(I)
    immlen = zero(I)
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
            stklen +=1
            idx += 3
        elseif ntype == TYPE_O
            n[TYPE_O] +=1
            @assert tt[idx+3] != 0
            if OP[tt[idx+2]]==:* 
                mxn_times = max(mxn_times,tt[idx+3])
            end
            stklen -= tt[idx+3]-1
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
            stklen += 1
            idx += 6

        elseif ntype == TYPE_O4
            n[TYPE_O4] += 1
            vmax = max(vmax,tt[idx+4])
            vnodes += 1
            stklen += 1
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
            stklen += 1
            idx += 5
        else
            @assert false
        end
        # @show vmax
        mx_stklen = max(mx_stklen, stklen)
    end
    @assert idx == length(tt) + 1
    @assert gnvar >= vmax && sum(n) == nnodes

    immlen = max(5, convert(I,(mxn_times-1)*mxn_times/2)+mxn_times)
    
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
            append!(tr,sub(istk,istklen-num+1:istklen))
            resize!(istk,istklen-num)
            push!(istk, nid)
            idx += 5

        elseif ntype == TYPE_O1
            nid += 1
            tt[idx+1] = nid
            push!(tr,pop!(istk))
            push!(istk,nid)
            idx += 5

        elseif ntype == TYPE_O2
            nid += 1
            tt[idx+1] = nid
            push!(tr,pop!(istk))
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
            push!(tr,pop!(istk))
            push!(istk,nid)
            idx += 5            

        elseif ntype == TYPE_O6
            nid += 1
            tt[idx+1] = nid
            push!(tr,pop!(istk))
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
    @assert length(istk) == 1 && idx == length(tt) + 1
    push!(tr,pop!(istk))
    @assert length(tr) == nnodes
    
    # @show nid, mxn, gnvar, vnodes, nnodes, immlen, stklen, mx_stklen
    return nid, mx_stklen, gnvar, vnodes, nnodes, immlen
end

function getMaxWorkingSize{I,V}(obj_tt::Tape{I,V},lag_tt::Tape{I,V}, cons_tt::Vector{Tape{I,V}})
    mx_vallen = max(obj_tt.nnode,lag_tt.nnode)
    mx_depth  = max(obj_tt.depth,lag_tt.depth)
    mx_stklen    = max(obj_tt.stklen,lag_tt.stklen)
    mx_immlen = max(obj_tt.immlen,lag_tt.immlen)
    # @show mx_vallen, mx_depth, mx_ops, mx_immlen, length(obj_tt.vals), length(obj_tt.stk), length(obj_tt.imm)

    for i=1:length(cons_tt)
        mx_immlen = max(mx_immlen,cons_tt[i].immlen)
        mx_vallen = max(mx_vallen,cons_tt[i].nnode)
        mx_depth = max(mx_depth,cons_tt[i].depth)
        mx_ops = max(mx_ops,cons_tt[i].stklen)
    end
    return mx_vallen, mx_depth, mx_ops, mx_immlen
end

function resizeHessStorage{I,V}(hs::HessStorage{I,V},bhlen::I,gnvar::I)
    @assert hs.bh_type == 1
    resize!(hs.bhi,bhlen)
    resize!(hs.bhv,bhlen)
    hs.gnvar = gnvar
    for i=1:bhlen
        @inbounds hs.bhi[i] = Vector{Int}()
        @inbounds hs.bhv[i] = Vector{Float64}()
    end
    resize!(hs.bh_idxes,bhlen)
    fill!(hs.bh_idxes,0)
end

##########################################################################################
#
# tape builder from a Julia Expr type
#   tape memory property is initialized 
#
##########################################################################################

function tapeBuilderNoHess{I,V}(expr::Expr, pvals::Vector{V}, gnvar::I)
    depth = [zero(I)]
    tt = Vector{I}()
    tr = Vector{I}()
    tape_builder(expr,tt, pvals,depth)
    nid, stklen, gnvar, vnodes, nnodes, immlen = tape_analysize(tt, tr, gnvar)
    tape = Tape{I,V}(tt, tr, stklen, gnvar, vnodes, nnodes, immlen, depth[1])
    return tape
end

function appendMultParam{I,V}(tape::Tape{I,V},pholder::Vector{I})
    tt = tape.tt
    push!(tt, TYPE_O1)
    push!(tt, -1)
    push!(tt, S_TO_OC[:*])
    push!(tt, -2)
    push!(pholder,length(tt))
    push!(tt, TYPE_O1)
end

function appendTapeMultParam{I,V}(to::Tape{I,V}, from::Tape{I,V}, pvals::Vector{V},ttstarts::Vector{I},ttends::Vector{I},trends::Vector{I},pholder::Vector{I})
    to.stklen = max(to.stklen,from.stklen)
    @assert to.nvar == from.nvar
    tt = to.tt
    tr = to.tr
    to.nvnode = max(to.nvnode, from.nvnode)
    to.nnode = max(to.nnode,from.nnode)
    to.immlen = max(to.immlen, from.immlen)
    to.depth = max(to.depth,from.depth)
    
    push!(ttstarts, length(tt)+1)
    append!(to.tt, from.tt)
    push!(ttends,   length(tt))
    append!(to.tr, from.tr)
    push!(trends,   length(tr))

    push!(tt, TYPE_O1)
    push!(tt, -1)
    push!(tt, S_TO_OC[:*])
    push!(tt, -2)
    push!(pholder,length(tt))
    push!(tt, TYPE_O1)
end

function buildSumTape{I,V}(tt::Vector{I}, mx_depth::I, n::I, gnvar::I, hs::HessStorage{I,V})
    tr = Vector{I}()
    if n > 1
        push!(tt, TYPE_O)
        push!(tt, -1)
        push!(tt, S_TO_OC[:+])
        push!(tt, n)
        push!(tt, TYPE_O)
        mx_depth += 1
    end
    nid, stklen, gnvar, vnodes, nnodes, immlen = tape_analysize(tt, tr, gnvar)
    tape = Tape{I,V}(tt, tr, stklen, gnvar, vnodes, nnodes, immlen, mx_depth)
    resizeHessStorage(hs,nid,gnvar)
    return tape
end

function mergeTapes{I,V}(obj_tape::Tape{I,V}, tapes::Vector{Tape{I,V}}, pvals::Vector{V}, gnvar::I,hs::HessStorage{I,V})
    tt = Vector{I}()
    tr = Vector{I}()
    mx_depth = zero(I)
    
    #objective
    mx_depth = max(obj_tape.depth,mx_depth)
    append!(tt, obj_tape.tt)
    push!(pvals,1.0)

    push!(tt,TYPE_O1)
    push!(tt,-1)
    push!(tt,S_TO_OC[:*])
    push!(tt,length(pvals))
    push!(tt,TYPE_O1)

    #constraints
    for i=1:length(tapes)
        mx_depth = max(tapes[i].depth, mx_depth)
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

    nid, stklen, gnvar, vnodes, nnodes, immlen = tape_analysize(tt, tr, gnvar)
    tape = Tape{I,V}(tt, tr, stklen, gnvar, vnodes, nnodes, immlen, mx_depth)
    resizeHessStorage(hs,nid,gnvar)

    return tape
end

type B_NODE{I,V}
    t::I
    vidx::I
    pval::V
end

function tape_builder{I,V}(expr,tt::Vector{I}, pvals::Vector{V}, depth::Vector{I}, d::I=0)
    # @show expr
    # @show d
    @assert length(depth) == 1
    d += 1
    depth[1] = max(depth[1],d)

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
            return tape_builder(expr.args[2],tt,pvals,depth,d)
        elseif op == :- && n==1
            node = tape_builder(expr.args[2],tt,pvals,depth,d)

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
            node = tape_builder(expr.args[2],tt,pvals,depth,d)
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
            ln = tape_builder(expr.args[2],tt,pvals,depth,d)
            rn = tape_builder(expr.args[3],tt,pvals,depth,d)
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
                node = tape_builder(expr.args[i],tt,pvals,depth,d)
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



function tape_cleanup{I,V}(tape::Tape{I,V})
    #output timing 
    @timing tape.enable_timing_stats begin
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
