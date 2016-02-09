#edge pusing algorithm for Hessian reverse AD

function reset_hess2(tape)
    assert(length(tape.live_vars) == length(tape.bh) == length(tape.bh_idxes) == tape.nnode+tape.nvar)  
    for i=1:tape.nnode+tape.nvar 
        tape.live_vars[i] = Vector{Int}() 
        tape.bh[i] = Vector{mPair{Int,Float64}}()
    end
    fill!(tape.bh_idxes,zero(Int))

    tape.h_I = Vector{I}() #hess_I
    tape.h_J = Vector{I}() #hess_J
    tape.hess = Vector{Float64}() #hess value
    tape.nzh = -one(Int)     #hess indicator
end

function prepare_reeval_hess2(tape)
    fill!(tape.bh_idxes,zero(Int))
end

@inline function push_edge(tape,to,from)
    # @show "edge - ",to," <--- ", from
    @inbounds push!(tape.bh[to],mPair{Int,Float64}(from,0.0))    
    # @show tape.bh
end

function hess_struct2{I,V}(tape::Tape{I,V})
    if(tape.nzh != -1)
        return tape.nzh
    end

    tt = tape.tt
    tr = tape.tr
    idx = length(tt)
    trlen = length(tr)
    # @show tape.nnode-1, length(tr)
    assert(tape.nnode-1 == length(tr))
    
    while (idx > 0)
        @inbounds ntype = tt[idx]
        idx -= 1
        if(ntype == TYPE_P)
            idx -= 2
        elseif(ntype == TYPE_V)
            @inbounds v_idx = tt[idx]
            idx -= 1
            i_idx = idx + tape.nvar
            idx -= 1

            #pushing edges pointing to this v node to the independent variables below
            # @show "pushing ", i_idx
            @inbounds i_num = tape.node_idx_to_number[i_idx]
            assert(i_num!=0 && v_idx <= tape.nvar)
            @inbounds lvi = tape.bh[i_num]
            for j = 1:length(lvi)
                # @inbounds (p_idx,w) = lvi[j]
                @inbounds p_idx = lvi[j].i
                @inbounds p_num = tape.node_idx_to_number[p_idx]
                if p_idx == i_idx
                    push_edge(tape,v_idx,v_idx)
                    # push!(tape.live_vars[v_idx],v_idx)  #v_num == v_idx  - by construction
                    # push!(tape.bh[v_idx],(v_idx,0.0))
                else
                    push_edge(tape,p_num,v_idx)
                    # push!(tape.live_vars[p_num],v_idx)
                    # push!(tape.bh[p_num],(v_idx,0.0))
                end
            end
        elseif(ntype == TYPE_O)
            @inbounds n = tt[idx]
            idx -= 1
            @inbounds oc = tt[idx]
            idx -= 1

            #pushing
            i_idx = idx  + tape.nvar #node i's idx on tape
            # @show "pushing",i_idx
            @inbounds i_num = tape.node_idx_to_number[i_idx]  #node i's number
            assert(i_num!=0)
            idx -= 1
            @inbounds lvi = tape.bh[i_num]
            for j = 1:length(lvi)
                # @inbounds (p_idx,w) = lvi[j] #index of node p
                @inbounds p_idx = lvi[j].i #index of node p
                @inbounds p_num = tape.node_idx_to_number[p_idx]
                assert(p_num!=0)
                if(p_idx == i_idx)
                    for j0=trlen-n+1:trlen 
                        @inbounds ci_idx = tr[j0] + tape.nvar
                        @inbounds ci_num = tape.node_idx_to_number[ci_idx]
                        assert(ci_num!=0)
                        # push!(tape.live_vars[ci_num],ci_idx)
                        push_edge(tape,ci_num,ci_idx)
                        for j1=j0+1:trlen
                            @inbounds cii_idx = tr[j1] + tape.nvar
                            @inbounds cii_num = tape.node_idx_to_number[cii_idx]
                            assert(cii_num!=0)
                            assert(ci_idx < cii_idx)
                            # push!(tape.live_vars[cii_num],ci_idx)    # ci_idx -> cii_idx, will handle cii_idx first
                            push_edge(tape,cii_num,ci_idx)
                        end
                    end
                else
                    for j0 = trlen -n + 1:trlen
                        @inbounds ci_idx = tr[j0] + tape.nvar
                        @inbounds ci_num = tape.node_idx_to_number[ci_idx]
                        assert(p_idx <= ci_idx)             
                        # push!(tape.live_vars[ci_num],p_idx)    #p_idx -> ci_idx
                        push_edge(tape,ci_num,p_idx)
                    end
                end
            end
            # @show "after pushing", tape.live_vars

            #creating 
            # @show "creating ",i_idx
            @inbounds op_sym = OP[oc]
            # @show op_sym
            if (op_sym ==:+ || op_sym == :-)
                #zeros
            elseif (n == 1 && op_sym != :-) #1-ary operator
                @inbounds ci_idx = tr[trlen] + tape.nvar
                @inbounds ci_num = tape.node_idx_to_number[ci_idx]
                assert(ci_num!=0)
                # push!(tape.live_vars[ci_num],ci_idx)
                push_edge(tape,ci_num,ci_idx)
            elseif (op_sym == :*)
                # @show "times ", n
                for j0 = trlen -n + 1:trlen
                    @inbounds ci_idx = tr[j0] + tape.nvar
                    @inbounds ci_num = tape.node_idx_to_number[ci_idx]
                    assert(ci_num!=0)
                    for j1=j0+1:trlen
                        @inbounds cii_idx = tr[j1] + tape.nvar
                        @inbounds cii_num = tape.node_idx_to_number[cii_idx]
                        assert(cii_num!=0)
                        assert(ci_idx < cii_idx)
                        # push!(tape.live_vars[cii_num],ci_idx)
                        push_edge(tape,cii_num,ci_idx)
                        # @show "push ",cii_num, cii_idx, "<--", ci_idx
                    end
                end 
            elseif (op_sym == :/) # binary operator /
                assert(false) #not tested , --- not implemented in operator
                assert(n==2)
                @inbounds ri_idx = tr[trlen] + tape.nvar
                @inbounds ri_num = tape.node_idx_to_number[ri_idx]
                assert(ri_num!=0)
                @inbounds li_idx = tr[trlen-1] + tape.nvar
                @inbounds li_num = tape.node_idx_to_number[li_idx]
                assert(li_num!=0) 
                assert(li_idx < ri_idx)
                # push!(tape.live_vars[ri_num],li_idx)
                push_edge(tape,ri_num,li_idx)
                # push!(tape.live_vars[ri_num],ri_idx)
                push_edge(tape,ri_num,li_idx)
            else # other binary
                assert(n==2)
                @inbounds ri_idx = tr[trlen] + tape.nvar
                @inbounds ri_num = tape.node_idx_to_number[ri_idx]
                assert(ri_num!=0)
                @inbounds li_idx = tr[trlen-1] + tape.nvar
                @inbounds li_num = tape.node_idx_to_number[li_idx]
                assert(li_num!=0)
                assert(li_idx < ri_idx)
                # push!(tape.live_vars[ri_num],ri_idx)
                push_edge(tape,li_num,li_idx)
                # push!(tape.live_vars[li_num],ri_idx)
                push_edge(tape,ri_num,li_idx)   #li (from) --> ri (to)
                # push!(tape.live_vars[li_num],li_idx)
                push_edge(tape,ri_num,ri_idx)
            end
            # @show "after creating", tape.live_vars
            trlen -= n
        end
    end #end while loop

    for i=1:tape.nvar
        @inbounds lvi = tape.bh[i] 
        for j = 1:length(lvi)
            # @inbounds (v_idx,w) = lvi[j]
            @inbounds v_idx = lvi[j].i
            if(v_idx<=tape.nvar)
                push!(tape.h_I,i)
                push!(tape.h_J,v_idx)
            end
        end
    end
    tape.nzh = length(tape.h_I)
    tape.hess = Vector{V}(tape.nzh)
    return tape.nzh
end

function forward_pass2_2ord{I,V}(tape::Tape{I,V}, vvals::Array{V,1}, pvals::Array{V,1})
    tt = tape.tt
    idx = one(I)
    # empty!(imm)  #used for immediate derrivative
    # empty!(tr)
    stk = tape.stk
    stklen = zero(I)
    imm = tape.imm2ord
    immlen = zero(I)
    
    @inbounds while(idx <= length(tt))
        # @show idx
        # println("++++++++++++++++++++++++++++++++++++")
        ntype = tt[idx]
        # eset[idx] = Dict{I,V}() #initialize edge set
        idx += 1
        if(ntype == TYPE_P)
            stklen += 1
            stk[stklen] = pvals[tt[idx]]
            idx += 2 #skip TYPE_P
        elseif(ntype == TYPE_V)
            stklen += 1
            stk[stklen] = vvals[tt[idx]]
            idx += 2 #skip TYPE_V
        elseif(ntype == TYPE_O)
            oc = tt[idx]
            idx += 1
            n = tt[idx]
            idx += 2 #skip TYPE_O
            # @show OP[oc], stklen-n+1, n
            # @show stk
            counter = zero(I)
            if(n==1)
                @inbounds (counter,stk[stklen]) = eval_2ord(OP[oc],stk[stklen],imm,immlen+1)
            else
                @inbounds (counter,val) = eval_2ord(OP[oc],stk,stklen-n+1,stklen,imm,immlen+1)
                stklen -= n-1
                @inbounds stk[stklen] = val
            end
            immlen += counter
            # @show imm     
            # @show stk[stklen] 
        end
        # @show stklen
        # println("++++++++++++++++++++++++++++++++++++")
    end
    # @show tape.imm2ord,immlen
    # assert(tape.imm2ord>=immlen)
    tape.imm2ordlen = immlen
    # @show stklen
    # resize!(imm,immlen) #if not resize memory will be hold in julia , and not claimed by gc
    return stk[1]
end

@inline function update(tape,to,from,w)
    # @show "update - ", to, "<-- ", from, w
    # @show tape.bh_idxes[to]
    # @inbounds 
    tape.bh_idxes[to] += 1  
    # @show tape.bh_idxes[to]  
    # @show tape.bh  
    # @inbounds 
    # assert(tape.bh[to][tape.bh_idxes[to]].i == from)
    tape.bh[to][tape.bh_idxes[to]].w = w
end

@inline function update_diag(tape,to,from,w)
    # @show "update diag "
    update(tape,to,from,w)
end

function reverse_pass2_2ord{I,V}(tape::Tape{I,V}, factor::V)
    assert(tape.nzh != -one(I))
    tr = tape.tr
    tt = tape.tt
    idx = length(tt)
    trlen = length(tr)
    imm = tape.imm2ord
    immlen = tape.imm2ordlen
    # assert(length(imm) == immlen)  

    adjs = tape.stk
    adjlen = 1
    adjs[1] = one(V)  #initialize adjoint = 1

    while(idx > 0)
        # println("++++++++++++++++++++++++++++++++++++")
        # @show idx
        # @show trlen,immlen, adjlen
        @inbounds ntype = tt[idx]
        idx -= 1
        
        #adjoints
        adj = adjs[adjlen]
        adjlen -= 1
        # @show adj


        if(ntype == TYPE_P)
            idx -= 2
        elseif(ntype == TYPE_V)
            @inbounds v_idx = tt[idx]
            idx -= 1
            i_idx = idx + tape.nvar
            idx -= 1
            # @show "pushing ",i_idx
            @inbounds i_num = tape.node_idx_to_number[i_idx]
            @inbounds lvi = tape.bh[i_num]
            for j = 1:length(lvi)
                # @inbounds (p_idx,w) = lvi[j]
                @inbounds p = lvi[j]
                p_idx = p.i
                w = p.w
                @inbounds p_num = tape.node_idx_to_number[p_idx]
                if p_idx == i_idx
                    update_diag(tape,v_idx,v_idx,w)
                elseif p_num == v_idx  #p belong to i's child
                    update(tape,v_idx,v_idx,2.0*w)
                else
                    update(tape,p_num,v_idx,w)
                end
            end
            # @show idx, vidx
        elseif(ntype == TYPE_O)
            @inbounds n = tt[idx]
            idx -= 1
            @inbounds oc = tt[idx]
            @inbounds op_sym = OP[oc]
            idx -= 1

            # @show OP[oc],n
            # @show tr

            i_idx = idx + tape.nvar
            idx -= 1
            @inbounds i_num = tape.node_idx_to_number[i_idx]
            #pushing   
            # @show "pushing ",i_idx,op_sym,n
            @inbounds lvi = tape.bh[i_num]
            for j = 1:length(lvi)
                # @inbounds (p_idx,w) = liv[j]
                @inbounds p = lvi[j]
                p_idx = p.i
                w = p.w
                @inbounds p_num = tape.node_idx_to_number[p_idx]
                if p_idx == i_idx
                    if (n==1 && op_sym != :-)
                        @inbounds ci_idx = tr[trlen] + tape.nvar
                        @inbounds ci_num = tape.node_idx_to_number[ci_idx]
                        @inbounds w_bar = imm[immlen-1]*imm[immlen-1]*w
                        update_diag(tape,ci_num,ci_idx,w_bar)
                    elseif op_sym == :+
                        for j0=trlen - n + 1:trlen
                            @inbounds ci_idx = tr[j0] + tape.nvar
                            @inbounds ci_num = tape.node_idx_to_number[ci_idx]
                            update_diag(tape,ci_num,ci_idx,w)
                            for j1=j0+1:trlen
                                @inbounds cii_idx = tr[j1] + tape.nvar
                                @inbounds cii_num = tape.node_idx_to_number[cii_idx]
                                update(tape,cii_num,ci_idx,w)
                            end  #j1 += 1
                        end  #j0 += 1
                    elseif op_sym == :-
                        @inbounds li_idx = tr[trlen-1] + tape.nvar
                        @inbounds li_num = tape.node_idx_to_number[li_idx]
                        @inbounds ri_idx = tr[trlen] + tape.nvar
                        @inbounds ri_num = tape.node_idx_to_number[ri_idx]
                        update_diag(tape,li_num,li_idx,w)
                        update_diag(tape,ri_num,ri_idx,w)
                        update(tape,ri_num,li_idx,-1.0*w)
                    elseif op_sym == :*
                        k = immlen - round(I,n+n*(n-1)/2)+1
                        for j0=trlen-n+1:trlen
                            @inbounds ci_idx = tr[j0] + tape.nvar
                            @inbounds ci_num = tape.node_idx_to_number[ci_idx]
                            @inbounds w_bar0 = imm[k] * imm[k] * w
                            update_diag(tape,ci_num,ci_idx,w_bar0)
                            k0 = k + 1
                            for j1=j0+1:trlen
                                @inbounds cii_idx = tr[j1] + tape.nvar
                                @inbounds cii_num = tape.node_idx_to_number[cii_idx]
                                @inbounds w_bar1 = imm[k]*imm[k0]*w
                                update(tape,cii_num,ci_idx,w_bar1)
                                k0 += 1
                            end
                            k += 1
                        end
                    else #other binary
                        @inbounds li_idx = tr[trlen-1] + tape.nvar
                        @inbounds li_num = tape.node_idx_to_number[li_idx]
                        @inbounds ri_idx = tr[trlen] + tape.nvar
                        @inbounds ri_num = tape.node_idx_to_number[ri_idx]
                        @inbounds dl = imm[immlen-4]
                        @inbounds dr = imm[immlen-3]
                        update_diag(tape,li_num,li_idx,dl*dl*w)
                        update(tape,ri_num,li_idx,dl*dr*w)
                        update_diag(tape,ri_num,ri_idx,dr*dr*w)
                    end
                else   #when p_idx != i_idx
                    if(n==1)
                        @inbounds ci_idx = tr[trlen] + tape.nvar
                        @inbounds ci_num = tape.node_idx_to_number[ci_idx]
                        @inbounds w_bar = imm[immlen-1]*imm[immlen-1]*w
                        update(tape,ci_num,p_idx,w_bar)
                    else
                        if(op_sym==:+)
                            for k=trlen-n+1:trlen
                                @inbounds ci_idx = tr[k] + tape.nvar
                                @inbounds ci_num = tape.node_idx_to_number[ci_idx]
                                # @show ci_idx , ci_num, p_idx, w
                                update(tape,ci_num,p_idx,w)
                            end
                        elseif(op_sym==:-)
                            assert(n==2)
                            @inbounds li_idx = tr[trlen-1] + tape.nvar
                            @inbounds li_num = tape.node_idx_to_number[li_idx]
                            @inbounds ri_idx = tr[trlen] + tape.nvar
                            @inbounds ri_num = tape.node_idx_to_number[ri_idx]
                            update(tape,li_num,p_idx,w)
                            update(tape,ri_num,p_idx,-1.0*w)
                        elseif(op_sym == :*)
                            j = immlen - round(I,n+n*(n-1)/2)+1
                            for k=trlen -n+1:trlen
                                @inbounds ci_idx = tr[k] + tape.nvar
                                @inbounds ci_num = tape.node_idx_to_number[ci_idx]
                                @inbounds w_bar = imm[j] * w
                                update(tape,ci_num,p_idx,w_bar)
                            end
                        else #other binary
                            @inbounds li_idx = tr[trlen-1] + tape.nvar
                            @inbounds li_num = tape.node_idx_to_number[li_idx]
                            @inbounds ri_idx = tr[trlen] + tape.nvar
                            @inbounds ri_num = tape.node_idx_to_number[ri_idx]
                            @inbounds dl = imm[immlen-4]
                            @inbounds dr = imm[immlen-3]
                            update(tape,li_num,p_idx,dl*w)
                            update(tape,ri_num,p_idx,dr*w)
                        end
                    end
                end #end p_idx == i_idx
            end

            #creating
            # @show "creating ",i_idx,op_sym,n
            if(op_sym == :+ || op_sym ==:-)
                #zeros
            elseif(n==1 && op_sym!=:-) #1-ary op
                @inbounds ci_idx = tr[trlen] + tape.nvar
                @inbounds ci_num = tape.node_idx_to_number[ci_idx]
                @inbounds w = adj*imm[immlen]
                update_diag(tape,ci_num,ci_idx,w)
            elseif(op_sym == :*)
                k = immlen - round(I,n*(n-1)/2)+1
                for j0=trlen-n+1:trlen
                    @inbounds ci_idx = tr[j0]+tape.nvar
                    @inbounds ci_num = tape.node_idx_to_number[ci_idx]
                    for j1=j0+1:trlen
                        @inbounds cii_idx = tr[j1]+tape.nvar
                        @inbounds cii_num = tape.node_idx_to_number[cii_idx]
                        @inbounds w = adj*imm[k]
                        update(tape,cii_num,ci_idx,w)
                        k+=1
                    end
                end
            elseif(op_sym == :/)
                assert(false) #not implemented in operator.jl
            else #other binary
                assert(n == 2)
                @inbounds ri_idx = tr[trlen] + tape.nvar
                @inbounds ri_num = tape.node_idx_to_number[ri_idx]
                @inbounds li_idx = tr[trlen-1] + tape.nvar
                @inbounds li_num = tape.node_idx_to_number[li_idx]
                @inbounds dll = imm[immlen-2]
                @inbounds dlr = imm[immlen-1]
                @inbounds drr = imm[immlen]
                update_diag(tape,li_num,li_idx,adj*dll)
                update(tape,ri_num,li_idx,adj*dlr)
                update_diag(tape,ri_num,ri_idx,adj*drr)
            end

            #update adjoints
            # @show adjlen, adjs
            imm_counter = zero(I)
            if op_sym == :+
                for j=1:n
                    adjlen += 1
                    @inbounds adjs[adjlen] = adj
                end
            elseif op_sym == :-
                adjlen += 1
                @inbounds adjs[adjlen] = adj
                adjlen += 1
                @inbounds adjs[adjlen] = -1.0*adj
            elseif n==1
                adjlen += 1
                @inbounds adjs[adjlen] = imm[immlen-1]*adj
                imm_counter = 2
            elseif op_sym == :*
                j = immlen - round(I,n+n*(n-1)/2)+1
                for m=1:n
                    adjlen += 1
                    # @show j,imm[j],adj
                    @inbounds adjs[adjlen] = imm[j] * adj
                    j+=1
                end
                imm_counter = round(I,n+n*(n-1)/2)
            else #other binary
                assert(n==2)
                adjlen += 1
                @inbounds adjs[adjlen] = imm[immlen-4]*adj
                adjlen += 1
                @inbounds adjs[adjlen] = imm[immlen-3]*adj
                imm_counter = 5
            end #end creating
            # @show adjlen, adjs
          
            # @show OP[oc],n
            # @show trlen
            # @show immlen,imm_counter
            # @show tr
            # @show imm

            #update counters
            trlen -= n
            immlen -= imm_counter
        end #end TYPE_O
        # println("++++++++++++++++++++++++++++++++++++")
    end  #end while
    assert(immlen == 0 && trlen == 0)
    # @show tape.eset

    # @show vidx
    nz = one(I)
    for i = 1:tape.nvar
        @inbounds lvi = tape.bh[i]
        for j=1:length(lvi)
            # @inbounds (v_idx,w) = lvi[j]
            @inbounds p = lvi[j]
            v_idx = p.i
            w = p.w
            if(v_idx<=tape.nvar)
                @inbounds tape.hess[nz] = w*factor
                nz += 1
            end
        end
    end
    return tape.nzh
end



#Interface function
function hess_structure2{I,V}(tape::Tape{I,V})
    return hess_struct2(tape)
end

function hess_reverse2{I,V}(tape::Tape{I,V},vvals::Vector{V},pvals::Vector{V})
    hess_reverse2(tape,vvals,pvals,1.0)
end

function hess_reverse2{I,V}(tape::Tape{I,V},vvals::Vector{V},pvals::Vector{V}, factor::V)
    forward_pass2_2ord(tape,vvals,pvals)
    reverse_pass2_2ord(tape,factor)
end
