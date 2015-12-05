
# forward pass on the tape tt, to build ss stack

#forward pass on the scalar function
function forward_pass_1ord{I,V}(tape::Tape{I,V}, vvals::Vector{V}, pvals::Vector{V})
    tt = tape.tt
    idx = one(I)
    stk = tape.stk
    stklen = 0
    imm = tape.imm1ord
    immlen = 0

    @inbounds while(idx <= length(tt))
        # @show idx
        # println("++++++++++++++++++++++++++++++++++++")
        ntype = tt[idx]
        idx += 1
        if(ntype == TYPE_P)
            # tic()
            val = pvals[tt[idx]]
            idx += 2 #skip TYPE_P
            stklen += 1
            stk[stklen] = val
        elseif(ntype == TYPE_V)
            val = vvals[tt[idx]]
            idx += 2 #skip TYPE_V
            stklen += 1
            stk[stklen] = val
        elseif(ntype == TYPE_O)
            oc = tt[idx]
            idx += 1
            n = tt[idx]
            idx += 1
            idx += 1 #skip TYPE_O
            # @show OP[oc], n, stklen-n+1, stklen
            # @show OP[oc], n, stk[stklen-n+1:stklen]
            if(n==1)
                # @show OP[oc]
                @inbounds stk[stklen] = eval_1ord(OP[oc],stk[stklen],imm,immlen+1)
            else
                # @show OP[oc],n, stklen, immlen
                @inbounds val = eval_1ord(OP[oc],stk,stklen-n+1,stklen,imm,immlen+1)
                stklen -= n-1
                @inbounds stk[stklen] = val
            end
            # @show stk[stklen]
            # @show imm[immlen+1:immlen+n]
            immlen += n
            # @show immlen
            # @show imm
        end
        # @show stklen
        # @show stk
        # println("++++++++++++++++++++++++++++++++++++")
    end
    return stk[1]
end

function reverse_pass_1ord{I,V}(tape::Tape{I,V})
    # @show imm
    # assert(length(imm) == tape.nnode -1)
    tt = tape.tt
    idx = length(tt)
    imm = tape.imm1ord
    immlen = length(imm)
    
    adjs = tape.stk
    adjlen = 1
    adjs[adjlen] = one(V)

    nnz = zero(I)

    @inbounds while(idx > 0)
        ntype = tt[idx]
        idx -= 1
        adj = adjs[adjlen]
        adjlen -= 1
        if(ntype == TYPE_P)
            idx -= 2
        elseif(ntype == TYPE_V)
            # @show tt[idx],adj
            # adj=isnan(adj)?0.0:adj
            nnz += 1
            tape.g[nnz]=adj
            idx -= 2
        elseif(ntype == TYPE_O)
            n = tt[idx]
            idx -= 3 #skip TYPE_O 
            @simd for i=immlen-n+1:immlen
                adjlen += 1
                adjs[adjlen] = imm[i]*adj
            end
            immlen -= n
        end
    end
    nothing
end

function grad_struct{I,V}(tape::Tape{I,V}) #repeated indexes, in reverse tracing order
    tt = tape.tt
    idx = length(tt)
    nnz = zero(I)
    @inbounds while(idx > 0)
        ntype = tt[idx]
            idx -= 1
        if(ntype == TYPE_V)
            nnz += 1
            tape.g_I[nnz] = tt[idx]
            idx -= 2
        elseif(ntype == TYPE_P)
            idx -= 2
        elseif(ntype == TYPE_O)
            idx -= 3
        end
    end
    tape.nzg = tape.nvnode
    assert(length(tape.g) == tape.nzg && nnz==tape.nzg)
end

#Interface function
function grad_structure{I,V}(tape::Tape{I,V})  #repeat version
    grad_struct(tape)
    return tape.g_I
end

function grad_reverse{I,V}(tape::Tape{I,V},vvals::Vector{V},pvals::Vector{V}) #sparse version
    forward_pass_1ord(tape,vvals,pvals)
    reverse_pass_1ord(tape)
    return tape.g
end

function grad_reverse{I,V}(tape::Tape{I,V},vvals::Vector{V},pvals::Vector{V}, g::Vector{V})  #dense version
    # assert(length(tape.g_I)==tape.nvnode)
    # assert(length(tape.g)==tape.nvnode)
    grad_reverse(tape,vvals,pvals)
    for i = 1:length(tape.g_I)
        @inbounds t = tape.g[i]
        @inbounds g[tape.g_I[i]] += t
    end
    nothing
end