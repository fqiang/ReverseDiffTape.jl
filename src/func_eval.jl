
#forward evaluation for a scalar function
function forward_pass_0ord{I,V}(tape::Tape{I,V}, vvals::Vector{V}, pvals::Vector{V})
    tt = tape.tt
    idx = one(I)
    
    stk = tape.stk
    stklen = zero(I)
    
    @inbounds while(idx <= length(tt))
        ntype = tt[idx]
        # @show ntype
        # @show stklen, stk
        
        # if ntype == TYPE_P
        #     @assert false
        
        # else
        if ntype == TYPE_V
            stklen += 1
            @inbounds stk[stklen] = vvals[tt[idx+1]]
            idx += 3

        elseif ntype == TYPE_O
            @inbounds oc = tt[idx+2]
            @inbounds n = tt[idx+3]

            if n == 1 # 1-argument functions
                @inbounds stk[stklen] = eval_0ord(OP[oc],stk[stklen])
            elseif n == 2
                # @show OP[oc],stk
                @inbounds val = eval_0ord(OP[oc],stk[stklen-1],stk[stklen])  
                stklen -= 1
                @inbounds stk[stklen] = val
            else 
                @inbounds val = eval_0ord(OP[oc],stk,stklen-n+1,stklen)
                stklen -= n-1
                @inbounds stk[stklen] = val
            end
            idx += 5

        elseif ntype == TYPE_O1
            @inbounds oc = tt[idx+2]
            @inbounds pval = pvals[tt[idx+3]]
            @inbounds stk[stklen] = eval_0ord(OP[oc],pval,stk[stklen])
            idx += 5
        
        elseif ntype == TYPE_O2
            @inbounds oc = tt[idx+2]
            @inbounds pval = pvals[tt[idx+3]]
            @inbounds stk[stklen] = eval_0ord(OP[oc],stk[stklen],pval)
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
            @inbounds stk[stklen] = eval_0ord(OP[oc],vval,stk[stklen])
            idx += 5

        elseif ntype == TYPE_O6
            @inbounds oc = tt[idx+2]
            @inbounds vval = vvals[tt[idx+3]]
            @inbounds stk[stklen] = eval_0ord(OP[oc],stk[stklen],vval)
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
    return stk[1]
end

## Interface method
function feval{I,V}(tape::Tape{I,V}, vvals::Array{V,1}, pvals::Array{V,1})
    return forward_pass_0ord(tape,vvals,pvals)
end
