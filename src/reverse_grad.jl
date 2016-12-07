
# forward pass on the tape tt, to build ss stack

#forward pass on the scalar function
function forward_pass_1ord{I,V}(tape::Tape{I,V}, vvals::Vector{V}, pvals::Vector{V})
    tt = tape.tt
    idx = one(I)
    stk = tape.stk
    vals = tape.vals
    vallen = zero(I)
    stklen = zero(I)

   @inbounds while(idx <= length(tt))
        ntype = tt[idx]
        # @show ntype
        if ntype == TYPE_P
            @assert false
        
        elseif ntype == TYPE_V
            stklen += 1
            @inbounds stk[stklen] = vvals[tt[idx+1]]
            idx += 3

        elseif ntype == TYPE_O
            @inbounds oc = tt[idx+2]
            @inbounds n = tt[idx+3]

            if n == 1 # 1-argument functions
                @inbounds val = eval_0ord(OP[oc],stk[stklen])
                vallen += 1
                @inbounds vals[vallen] = stk[stklen]
                @inbounds stk[stklen] = val
            elseif n == 2
                # @show OP[oc],stk
                @inbounds val = eval_0ord(OP[oc],stk[stklen-1],stk[stklen]) 
                @inbounds vals[vallen+1] = stk[stklen-1]
                @inbounds vals[vallen+2] = stk[stklen]
                vallen += 2
                stklen -= 1
                @inbounds stk[stklen] = val
            else 
                @inbounds val = eval_0ord(OP[oc],stk,stklen-n+1,stklen)
                @simd for i=1:n
                    @inbounds vals[vallen+i] = stk[stklen-n+i]
                end
                vallen += n
                stklen -= n-1
                @inbounds stk[stklen] = val
            end
            idx += 5

        elseif ntype == TYPE_O1
            @inbounds oc = tt[idx+2]
            @inbounds pval = pvals[tt[idx+3]]
            @inbounds val = eval_0ord(OP[oc],pval,stk[stklen])
            vallen += 1
            @inbounds vals[vallen] = stk[stklen]
            @inbounds stk[stklen] = val
            idx += 5
        
        elseif ntype == TYPE_O2
            @inbounds oc = tt[idx+2]
            @inbounds pval = pvals[tt[idx+3]]
            @inbounds val = eval_0ord(OP[oc],stk[stklen],pval)
            vallen += 1
            @inbounds vals[vallen] = stk[stklen]
            @inbounds stk[stklen] = val
            idx += 5
        
        elseif ntype == TYPE_O3
            @inbounds oc = tt[idx+2]
            @inbounds pval = pvals[tt[idx+3]]
            @inbounds vval = vvals[tt[idx+4]]
            stklen += 1
            @inbounds stk[stklen] = eval_0ord(OP[oc],pval,vval)
            idx += 6
        
        elseif ntype == TYPE_O4
            @inbounds oc = tt[idx+2]
            @inbounds pval = pvals[tt[idx+3]]
            @inbounds vval = vvals[tt[idx+4]]
            stklen += 1
            @inbounds stk[stklen] = eval_0ord(OP[oc],vval,pval)
            idx += 6

        elseif ntype == TYPE_O5
            @inbounds oc = tt[idx+2]
            @inbounds vval = vvals[tt[idx+3]]
            @inbounds val = eval_0ord(OP[oc],vval,stk[stklen])
            vallen += 1
            @inbounds vals[vallen] = stk[stklen]
            @inbounds stk[stklen] = val
            idx += 5

        elseif ntype == TYPE_O6
            @inbounds oc = tt[idx+2]
            @inbounds vval = vvals[tt[idx+3]]
            @inbounds val = eval_0ord(OP[oc],stk[stklen],vval)
            vallen += 1
            @inbounds vals[vallen] = stk[stklen]
            @inbounds stk[stklen] = val
            idx += 5

        elseif ntype == TYPE_O7
            @inbounds oc = tt[idx+2]
            @inbounds vval = vvals[tt[idx+3]]
            stklen += 1
            @inbounds stk[stklen] = eval_0ord(OP[oc],vval)
            idx += 5

        # else 
        #     @assert false
        end
    end
    @assert stklen == 1 && vallen + 1 == tape.nnode
    @inbounds ret = vals[vallen + 1] = stk[stklen]
    return ret
end

function reverse_pass_1ord{I,V}(tape::Tape{I,V},vvals::Vector{V}, pvals::Vector{V},start::I, g::Vector{V}) 
    tt = tape.tt
    idx = length(tt)
    vals = tape.vals
    vallen = tape.nnode - 1
    imm = tape.imm
    
    adjs = tape.stk
    adjlen = one(I)
    @inbounds adjs[adjlen] = one(V) #initial value

    i = start
    @inbounds while idx > 0
        ntype = tt[idx]
        # @show i, start, ntype, idx
        # @show ntype, vallen , adjlen
        @inbounds adj = adjs[adjlen]
        
        if ntype == TYPE_P
            @assert false

        elseif ntype == TYPE_V
            @inbounds g[i] = adj
            i += 1
            adjlen -= 1
            idx -= 3

        elseif ntype == TYPE_O
            @inbounds n = tt[idx-1]
            @inbounds oc = tt[idx-2]

            # @show n
            if n == 1
                @inbounds eval_1ord(OP[oc],vals[vallen],imm)
                vallen -= 1
                @inbounds adjs[adjlen] = imm[1]*adj
            elseif n == 2
                @inbounds eval_1ord(OP[oc],vals[vallen-1],vals[vallen],imm)
                vallen -= 2

                @inbounds adjs[adjlen] = imm[1]*adj
                adjlen +=1
                @inbounds adjs[adjlen] = imm[2]*adj
            else
                @inbounds eval_1ord(OP[oc],vals,vallen-n+1,vallen,imm)
                vallen -= n
                @inbounds op_sym = OP[oc]
                if op_sym == :+
                    @simd for j=1:n
                        @inbounds adjs[adjlen+j-1] = adj
                    end
                else
                    @assert op_sym == :*
                    @simd for j=1:n
                        @inbounds adjs[adjlen+j-1] = imm[j]*adj
                    end
                end
                adjlen += n-1
            end
            idx -= 5

        elseif ntype == TYPE_O1
            @inbounds oc = tt[idx-2]
            @inbounds pval = pvals[tt[idx-1]]
            @inbounds eval_1ord(OP[oc],pval,vals[vallen],imm)
            vallen -= 1
            @inbounds adjs[adjlen] = imm[2]*adj
            idx -= 5

        elseif ntype == TYPE_O2
            @inbounds oc = tt[idx-2]
            @inbounds pval = pvals[tt[idx-1]]
            @inbounds eval_1ord(OP[oc],vals[vallen],pval,imm)
            vallen -= 1
            @inbounds adjs[adjlen] = imm[1]*adj
            idx -= 5

        elseif ntype == TYPE_O3
            @inbounds oc = tt[idx-3]
            @inbounds vval = vvals[tt[idx-1]]
            @inbounds pval = pvals[tt[idx-2]]
            @inbounds eval_1ord(OP[oc],pval,vval,imm)
            @inbounds g[i] = imm[2]*adj
            i += 1
            adjlen -= 1
            idx -= 6

        elseif ntype == TYPE_O4
            @inbounds oc = tt[idx-3]
            @inbounds vval = vvals[tt[idx-1]]
            @inbounds pval = pvals[tt[idx-2]]
            @inbounds eval_1ord(OP[oc],vval,pval,imm)
            @inbounds g[i] = imm[1]*adj
            i += 1
            adjlen -= 1
            idx -= 6

        elseif ntype == TYPE_O5
            @inbounds oc = tt[idx-2]
            @inbounds vval = vvals[tt[idx-1]]
            @inbounds eval_1ord(OP[oc],vval,vals[vallen],imm)
            @inbounds g[i] = imm[1]*adj
            i += 1
            @inbounds adjs[adjlen] = imm[2]*adj
            vallen -= 1
            idx -= 5

        elseif ntype == TYPE_O6
            @inbounds oc = tt[idx-2]
            @inbounds vval = vvals[tt[idx-1]]
            @inbounds eval_1ord(OP[oc],vals[vallen],vval,imm)
            @inbounds g[i] = imm[2]*adj
            i += 1
            @inbounds adjs[adjlen] = imm[1]*adj
            vallen -= 1
            idx -= 5

        elseif ntype == TYPE_O7
            @inbounds oc = tt[idx-2]
            @inbounds vval = vvals[tt[idx-1]]
            @inbounds eval_1ord(OP[oc],vval,imm)
            @inbounds g[i] = imm[1]*adj
            i += 1
            adjlen -= 1
            idx -= 5

        else 
            @assert false
        end
    end
    @assert idx == 0 && adjlen == 0 && vallen == 0
    @assert (i-start) == tape.nzg  "$i $start $(tape.nzg) $(tape.nvnode)"
    nothing
end

function reverse_pass_1ord_dense{I,V}(tape::Tape{I,V},vvals::Vector{V}, pvals::Vector{V}, g::Vector{V})  #repeated indicies
    # @show imm
    # assert(length(imm) == tape.nnode -1)
    tt = tape.tt
    idx = length(tt)
    vals = tape.vals
    vallen = tape.nnode - 1
    imm = tape.imm
    
    adjs = tape.stk
    adjlen = one(I)
    @inbounds adjs[adjlen] = one(V) #initial value

    
    @inbounds while idx > 0
        ntype = tt[idx]
        # @show ntype, vallen , adjlen
        @inbounds adj = adjs[adjlen]
        
        if ntype == TYPE_P
            @assert false

        elseif ntype == TYPE_V
            @inbounds g[tt[idx-1]] += adj
            adjlen -= 1
            idx -= 3

        elseif ntype == TYPE_O
            @inbounds n = tt[idx-1]
            @inbounds oc = tt[idx-2]

            # @show n
            if n == 1
                @inbounds eval_1ord(OP[oc],vals[vallen],imm)
                vallen -= 1
                @inbounds adjs[adjlen] = imm[1]*adj
            elseif n == 2
                @inbounds eval_1ord(OP[oc],vals[vallen-1],vals[vallen],imm)
                vallen -= 2

                @inbounds adjs[adjlen] = imm[1]*adj
                adjlen +=1
                @inbounds adjs[adjlen] = imm[2]*adj
            else
                @inbounds eval_1ord(OP[oc],vals,vallen-n+1,vallen,imm)
                vallen -= n
                @inbounds op_sym = OP[oc]
                if op_sym == :+
                    @simd for i=1:n
                        @inbounds adjs[adjlen+i-1] = adj
                    end
                else
                    @assert op_sym == :*
                    @simd for i=1:n
                        @inbounds adjs[adjlen+i-1] = imm[i]*adj
                    end
                end
                adjlen += n-1
            end
            idx -= 5

        elseif ntype == TYPE_O1
            @inbounds oc = tt[idx-2]
            @inbounds pval = pvals[tt[idx-1]]
            @inbounds eval_1ord(OP[oc],pval,vals[vallen],imm)
            vallen -= 1
            @inbounds adjs[adjlen] = imm[2]*adj
            idx -= 5

        elseif ntype == TYPE_O2
            @inbounds oc = tt[idx-2]
            @inbounds pval = pvals[tt[idx-1]]
            @inbounds eval_1ord(OP[oc],vals[vallen],pval,imm)
            vallen -= 1
            @inbounds adjs[adjlen] = imm[1]*adj
            idx -= 5

        elseif ntype == TYPE_O3
            @inbounds oc = tt[idx-3]
            @inbounds vidx = tt[idx-1]
            @inbounds vval = vvals[vidx]
            @inbounds pval = pvals[tt[idx-2]]
            @inbounds eval_1ord(OP[oc],pval,vval,imm)
            @inbounds g[vidx] += imm[2]*adj
            adjlen -= 1
            idx -= 6

        elseif ntype == TYPE_O4
            @inbounds oc = tt[idx-3]
            @inbounds vidx = tt[idx-1]
            @inbounds vval = vvals[vidx]
            @inbounds pval = pvals[tt[idx-2]]
            @inbounds eval_1ord(OP[oc],vval,pval,imm)
            @inbounds g[vidx] += imm[1]*adj
            adjlen -= 1
            idx -= 6

        elseif ntype == TYPE_O5
            @inbounds oc = tt[idx-2]
            @inbounds vidx = tt[idx-1]
            @inbounds vval = vvals[vidx]
            @inbounds eval_1ord(OP[oc],vval,vals[vallen],imm)
            @inbounds g[vidx] += imm[1]*adj
            @inbounds adjs[adjlen] = imm[2]*adj
            vallen -= 1
            idx -= 5

        elseif ntype == TYPE_O6
            @inbounds oc = tt[idx-2]
            @inbounds vidx = tt[idx-1]
            @inbounds vval = vvals[vidx]
            @inbounds eval_1ord(OP[oc],vals[vallen],vval,imm)
            @inbounds g[vidx] += imm[2]*adj
            @inbounds adjs[adjlen] = imm[1]*adj
            vallen -= 1
            idx -= 5

        elseif ntype == TYPE_O7
            @inbounds oc = tt[idx-2]
            @inbounds vidx = tt[idx-1]
            @inbounds vval = vvals[vidx]
            @inbounds eval_1ord(OP[oc],vval,imm)
            @inbounds g[vidx] += imm[1]*adj
            adjlen -= 1
            idx -= 5

        else 
            @assert false
        end
    end
    @assert idx == 0 
    @assert adjlen == 0
    @assert vallen == 0
    
    nothing
end


function grad_struct{I,V}(tape::Tape{I,V}, iJ::Vector) #repeated indexes, in reverse tracing order
    if tape.nzg!=-1
        return tape.nzg
    end
    tt = tape.tt
    idx = length(tt)
    start = length(iJ) + 1
    append!(iJ,zeros(tape.nvnode)) 
    i = start
    @inbounds while(idx > 0)
        ntype = tt[idx]
        # @show i, start, ntype,idx
        if ntype == TYPE_P
            @assert false
        elseif ntype == TYPE_V
            @inbounds iJ[i] = tt[idx-1]
            i += 1
            idx -= 3
        elseif ntype == TYPE_O
            idx -= 5
        elseif ntype == TYPE_O1 || ntype == TYPE_O2
            idx -= 5
        elseif ntype == TYPE_O3 || ntype == TYPE_O4
            @inbounds iJ[i] = tt[idx-1]
            i+=1
            idx -= 6
        elseif ntype == TYPE_O5 || ntype == TYPE_O6 || ntype == TYPE_O7
            @inbounds iJ[i] = tt[idx-1]
            i+=1
            idx -= 5
        # else
        #     @assert false
        end
    end
    tape.nzg = i - start
    @assert tape.nzg == tape.nvnode == length(iJ) - start + 1 "$(tape.nzg) $i $start $(tape.nvnode) $(length(iJ))"
    return tape.nzg
end

#Interface function
function grad_structure{I,V}(tape::Tape{I,V}, iJ::Vector{I})  
    return grad_struct(tape, iJ)
end

function grad_reverse{I,V}(tape::Tape{I,V},vvals::Vector{V},pvals::Vector{V}, start::I, g::Vector{V}) #sparse version
    forward_pass_1ord(tape,vvals,pvals)
    reverse_pass_1ord(tape,vvals,pvals,start,g) 
end

function grad_reverse_dense{I,V}(tape::Tape{I,V},vvals::Vector{V},pvals::Vector{V}, g::Vector{V}) #sparse version
    forward_pass_1ord(tape,vvals,pvals)
    reverse_pass_1ord_dense(tape,vvals,pvals,g) 
end
