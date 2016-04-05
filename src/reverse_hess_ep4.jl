#edge pusing algorithm for Hessian reverse AD

function reset_hess3(tape)
    assert(length(tape.bh3) == tape.nnode+tape.nvar)  
    for i=1:tape.nnode+tape.nvar 
        empty!(tape.bh3[i])
    end
    
    empty!(tape.h_I)  #hess_I
    empty!(tape.h_J)  #hess_J
    empty!(tape.hess) #hess value
    tape.nzh = -one(Int)     #hess indicator
end

function prepare_reeval_hess3(tape)
    for i=1:tape.nnode+tape.nvar
        d = tape.bh3[i]
        for k in keys(d)
            d[k] = 0.0
        end
    end
end

@inline function push_edge3(tape,to,from)
    # @show "push_edge3 - ",to," <--- ", from
    @inbounds tape.bh3[to][from] = 0.0  
    # @show tape.bh
end

function hess_struct3{I,V}(tape::Tape{I,V})
    tape.h_type = 3
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
            @inbounds d = tape.bh3[i_num]
            for p_idx in keys(d)
                # @inbounds (p_idx,w) = lvi[j]
                @inbounds p_num = tape.node_idx_to_number[p_idx]
                if p_idx == i_idx
                    push_edge3(tape,v_idx,v_idx)
                    # push!(tape.live_vars[v_idx],v_idx)  #v_num == v_idx  - by construction
                    # push!(tape.bh[v_idx],(v_idx,0.0))
                else
                    push_edge3(tape,p_num,v_idx)
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
            @inbounds d = tape.bh3[i_num]
            for p_idx in keys(d)
                # @inbounds (p_idx,w) = lvi[j] #index of node p
                @inbounds p_num = tape.node_idx_to_number[p_idx]
                assert(p_num!=0)
                if(p_idx == i_idx)
                    for j0=trlen-n+1:trlen 
                        @inbounds ci_idx = tr[j0] + tape.nvar
                        @inbounds ci_num = tape.node_idx_to_number[ci_idx]
                        assert(ci_num!=0)
                        # push!(tape.live_vars[ci_num],ci_idx)
                        push_edge3(tape,ci_num,ci_idx)
                        for j1=j0+1:trlen
                            @inbounds cii_idx = tr[j1] + tape.nvar
                            @inbounds cii_num = tape.node_idx_to_number[cii_idx]
                            assert(cii_num!=0)
                            assert(ci_idx < cii_idx)
                            # push!(tape.live_vars[cii_num],ci_idx)    # ci_idx -> cii_idx, will handle cii_idx first
                            push_edge3(tape,cii_num,ci_idx)
                        end
                    end
                else  #when i_idx != p_idx
                    for j0 = trlen -n + 1:trlen
                        @inbounds ci_idx = tr[j0] + tape.nvar
                        @inbounds ci_num = tape.node_idx_to_number[ci_idx]
                        assert(p_idx <= ci_idx)             
                        # push!(tape.live_vars[ci_num],p_idx)    #p_idx -> ci_idx
                        push_edge3(tape,ci_num,p_idx)
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
                push_edge3(tape,ci_num,ci_idx)
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
                        push_edge3(tape,cii_num,ci_idx)
                        # @show "push ",cii_num, cii_idx, "<--", ci_idx
                    end
                end 
            elseif (op_sym == :/) # binary operator /
                # assert(false) #not tested , --- not implemented in operator
                assert(n==2)
                @inbounds ri_idx = tr[trlen] + tape.nvar
                @inbounds ri_num = tape.node_idx_to_number[ri_idx]
                assert(ri_num!=0)
                @inbounds li_idx = tr[trlen-1] + tape.nvar
                @inbounds li_num = tape.node_idx_to_number[li_idx]
                assert(li_num!=0) 
                assert(li_idx < ri_idx)
                # push!(tape.live_vars[ri_num],li_idx)
                push_edge3(tape,ri_num,li_idx)
                # push!(tape.live_vars[ri_num],ri_idx)
                push_edge3(tape,ri_num,ri_idx)
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
                push_edge3(tape,li_num,li_idx)
                # push!(tape.live_vars[li_num],ri_idx)
                push_edge3(tape,ri_num,li_idx)   #li (from) --> ri (to)
                # push!(tape.live_vars[li_num],li_idx)
                push_edge3(tape,ri_num,ri_idx)
            end
            # @show "after creating", tape.live_vars
            trlen -= n
        end
    end #end while loop

    for i=1:tape.nvar
        @inbounds d = tape.bh3[i] 
        for v_idx in keys(d)
            if(v_idx <=tape.nvar)
                if(v_idx < i)
                    push!(tape.h_I,i)
                    push!(tape.h_J,v_idx)
                else
                    push!(tape.h_J,i)
                    push!(tape.h_I,v_idx)
                end
            end
        end
    end
    tape.nzh = length(tape.h_I)
    resize!(tape.hess, tape.nzh)
    return tape.nzh
end

function forward_pass3_2ord{I,V}(tape::Tape{I,V}, vvals::Array{V,1}, pvals::Array{V,1})
    tt = tape.tt
    idx = one(I)
    stk = tape.stk
    stklen = zero(I)
    imm = tape.imm
    immlen = zero(I)
    
    while(idx <= length(tt))
        # @show idx
        # println("++++++++++++++++++++++++++++++++++++")
        @inbounds ntype = tt[idx]
        # eset[idx] = Dict{I,V}() #initialize edge set
        idx += 1
        if(ntype == TYPE_P)
            stklen += 1
            @inbounds stk[stklen] = pvals[tt[idx]]
            idx += 2 #skip TYPE_P
        elseif(ntype == TYPE_V)
            stklen += 1
            @inbounds stk[stklen] = vvals[tt[idx]]
            idx += 2 #skip TYPE_V
        elseif(ntype == TYPE_O)
            @inbounds oc = tt[idx]
            idx += 1
            @inbounds n = tt[idx]
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
    resize!(tape.imm,max(tape.imm1ordlen,tape.imm2ordlen)) #if not resize memory will be hold in julia , and not claimed by gc
    return @inbounds stk[1]
end

@inline function update3(tape,to,from,w)
    # @show "update3 - ", to, "<-- ", from, tape.bh3[to][from], w, tape.bh3[to][from]+w
    @inbounds tape.bh3[to][from] += w
end

function reverse_pass3_2ord{I,V}(tape::Tape{I,V}, factor::V)
    assert(tape.nzh != -one(I))
    tr = tape.tr
    tt = tape.tt
    idx = length(tt)
    trlen = length(tr)
    imm = tape.imm
    immlen = tape.imm2ordlen
    # assert(length(imm) == immlen)  

    adjs = tape.stk
    adjlen = 1
    @inbounds adjs[1] = one(V)  #initialize adjoint = 1

    while(idx > 0)
        # println("++++++++++++++++++++++++++++++++++++")
        # @show idx
        # @show trlen,immlen, adjlen
        @inbounds ntype = tt[idx]
        idx -= 1
        
        #adjoints
        @inbounds adj = adjs[adjlen]
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
            @inbounds d = tape.bh3[i_num]
            for p_idx in keys(d)
                # @inbounds (p_idx,w) = lvi[j]
                w = d[p_idx]
                @inbounds p_num = tape.node_idx_to_number[p_idx]
                if p_idx == i_idx
                    update3(tape,v_idx,v_idx,w)
                elseif p_idx == v_idx  #p belong to i's child
                    update3(tape,v_idx,v_idx,2.0*w)
                else
                    update3(tape,p_num,v_idx,w)
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
            @inbounds d = tape.bh3[i_num]
            @inbounds tr0_idx = tr[trlen] + tape.nvar
            @inbounds tr0_num = tape.node_idx_to_number[tr0_idx]
            
            imm_counter = zero(I)

            if n==1
                #pushing
                @inbounds t0 = imm[immlen-1]
                for p_idx in keys(d)
                    w = d[p_idx]
                    # @inbounds p_num = tape.node_idx_to_number[p_idx]
                    if p_idx == i_idx                       
                        w_bar = t0*t0*w
                        update3(tape,tr0_num,tr0_idx,w_bar)     
                    else
                        w_bar = t0*w
                        update3(tape,tr0_num,p_idx,w_bar)
                    end
                end
                #creating
                if(op_sym != :-)
                    @inbounds w_bar = adj*imm[immlen]
                    update3(tape,tr0_num,tr0_idx,w_bar)
                end

                #updating 
                adjlen += 1
                @inbounds adjs[adjlen] = t0*adj
                imm_counter = 2
            elseif n == 2
                @inbounds li_idx = tr[trlen-1] + tape.nvar
                @inbounds li_num = tape.node_idx_to_number[li_idx]

                if op_sym == :/
                    #pushing
                    @inbounds dl = imm[immlen-4]
                    @inbounds dr = imm[immlen-3]
                    
                    for p_idx in keys(d)
                        w = d[p_idx]
                    
                        dlw = dl*w
                        drw = dr*w
                        if p_idx == i_idx
                            update3(tape,li_num,li_idx,dl*dlw)
                            update3(tape,tr0_num,li_idx,dl*drw)
                            update3(tape,tr0_num,tr0_idx,dr*drw)
                        else 
                            update3(tape,li_num,p_idx,dlw)
                            update3(tape,tr0_num,p_idx,drw)
                        end
                    end
                    #creating
                    @inbounds dll = imm[immlen-2]
                    @inbounds dlr = imm[immlen-1]
                    @inbounds drr = imm[immlen]
                    update3(tape,tr0_num,li_idx,adj*dlr)
                    update3(tape,tr0_num,tr0_idx,adj*drr)
                    
                    #updating
                    imm_counter = 5
                    adjlen += 1
                    @inbounds adjs[adjlen] = dl*adj
                    adjlen += 1
                    @inbounds adjs[adjlen] = dr*adj
                elseif op_sym == :*
                    #pushing
                    @inbounds dl = imm[immlen-2]
                    @inbounds dr = imm[immlen-1]
                    
                    for p_idx in keys(d)
                        w = d[p_idx]
                    
                        dlw = dl*w
                        drw = dr*w
                        if p_idx == i_idx
                            update3(tape,li_num,li_idx,dl*dlw)
                            update3(tape,tr0_num,li_idx,dl*drw)
                            update3(tape,tr0_num,tr0_idx,dr*drw)
                        else 
                            update3(tape,li_num,p_idx,dlw)
                            update3(tape,tr0_num,p_idx,drw)
                        end
                    end
                    #creating
                    update3(tape,tr0_num,li_idx,adj)  #adj*1.0
                    
                    #updating
                    imm_counter = 3
                    adjlen += 1
                    @inbounds adjs[adjlen] = dl*adj
                    adjlen += 1
                    @inbounds adjs[adjlen] = dr*adj
                elseif  op_sym == :-
                    #pushing
                    for p_idx in keys(d)
                        w = d[p_idx]
                        if p_idx == i_idx
                            update3(tape,li_num,li_idx,w)
                            update3(tape,tr0_num,li_idx,-1.0*w)
                            update3(tape,tr0_num,tr0_idx,w)
                        else
                            update3(tape,li_num,p_idx,w)
                            update3(tape,tr0_num,p_idx,-1.0*w)
                        end
                    end
                    #creating
                        #zeros

                    #updating
                    adjlen += 1
                    @inbounds adjs[adjlen] = adj
                    adjlen += 1
                    @inbounds adjs[adjlen] = -1.0*adj
                    imm_counter = 0
                elseif op_sym == :+
                    #pushing
                    for p_idx in keys(d)
                        w = d[p_idx]
                        if p_idx == i_idx
                            update3(tape,li_num,li_idx,w)
                            update3(tape,tr0_num,li_idx,w)
                            update3(tape,tr0_num,tr0_idx,w)
                        else
                            update3(tape,li_num,p_idx,w)
                            update3(tape,tr0_num,p_idx,w)
                        end
                    end
                    #creating
                        #zeros

                    #updating
                    adjlen += 1
                    @inbounds adjs[adjlen] = adj
                    adjlen += 1
                    @inbounds adjs[adjlen] = adj
                    imm_counter = 0

                else  #other binary ops
                    #pushing
                    @inbounds dl = imm[immlen-4]
                    @inbounds dr = imm[immlen-3]
                    
                    for p_idx in keys(d)
                        w = d[p_idx]
                    
                        dlw = dl*w
                        drw = dr*w
                        if p_idx == i_idx
                            update3(tape,li_num,li_idx,dl*dlw)
                            update3(tape,tr0_num,li_idx,dl*drw)
                            update3(tape,tr0_num,tr0_idx,dr*drw)
                        else 
                            update3(tape,li_num,p_idx,dlw)
                            update3(tape,tr0_num,p_idx,drw)
                        end
                    end
                    #creating
                    @inbounds dll = imm[immlen-2]
                    @inbounds dlr = imm[immlen-1]
                    @inbounds drr = imm[immlen]
                    #@show li_num, li_idx, adj, dll
                    update3(tape,li_num,li_idx,adj*dll)
                    update3(tape,tr0_num,li_idx,adj*dlr)
                    update3(tape,tr0_num,tr0_idx,adj*drr)

                    #updating
                    imm_counter = 5
                    adjlen += 1
                    @inbounds adjs[adjlen] = dl*adj
                    adjlen += 1
                    @inbounds adjs[adjlen] = dr*adj
                end  #end else binary ops
            else  #other + or * with n - operands
                if op_sym == :+
                    #pushing
                    for p_idx in keys(d)
                        w = d[p_idx]
                        # @inbounds p_num = tape.node_idx_to_number[p_idx]
                        if p_idx == i_idx     
                            for j0=trlen-n+1:trlen
                                @inbounds ci_idx =tr[j0] + tape.nvar
                                @inbounds ci_num = tape.node_idx_to_number[ci_idx]
                                update3(tape,ci_num,ci_idx,w)
                                for j1=j0+1:trlen
                                    @inbounds cii_idx = tr[j1] + tape.nvar
                                    @inbounds cii_num = tape.node_idx_to_number[cii_idx]
                                    update3(tape,cii_num,ci_idx,w)
                                end #j1 +=1
                            end #j0+=1
                        else #p_idx != i_idx
                            for k=trlen-n+1:trlen
                                @inbounds ci_idx = tr[k] + tape.nvar
                                @inbounds ci_num = tape.node_idx_to_number[ci_idx]
                                # @show ci_idx , ci_num, p_idx, w
                                update3(tape,ci_num,p_idx,w)
                            end
                        end
                    end
                    #creating
                        #zero

                    #updating
                    for j=1:n
                        adjlen += 1
                        @inbounds adjs[adjlen] = adj
                    end
                    imm_counter = 0
                elseif op_sym == :*
                    #pushing        
                    r = immlen - round(I,n+n*(n-1)/2)+1        
                    for p_idx in keys(d)
                        w = d[p_idx]
                        
                        k = r
                        if p_idx == i_idx
                            for j0=trlen-n+1:trlen
                                @inbounds ci_idx = tr[j0] + tape.nvar
                                @inbounds ci_num = tape.node_idx_to_number[ci_idx]
                                @inbounds t0 = imm[k]
                                t1 = t0*w
                                @inbounds w_bar0 = t0 * t1
                                update3(tape,ci_num,ci_idx,w_bar0)
                                k0 = k + 1
                                for j1=j0+1:trlen
                                    @inbounds cii_idx = tr[j1] + tape.nvar
                                    @inbounds cii_num = tape.node_idx_to_number[cii_idx]
                                    @inbounds w_bar1 = t1*imm[k0]
                                    update3(tape,cii_num,ci_idx,w_bar1)
                                    k0 += 1
                                end
                                k += 1
                            end
                        else #p_idx != i_idx 
                            for j=trlen -n+1:trlen
                                @inbounds ci_idx = tr[j] + tape.nvar
                                @inbounds ci_num = tape.node_idx_to_number[ci_idx]
                                @inbounds w_bar = imm[k] * w
                                update3(tape,ci_num,p_idx,w_bar)
                                k += 1
                            end
                        end
                    end
                    #creating
                    k = r + n
                    for j0=trlen-n+1:trlen
                        @inbounds ci_idx = tr[j0]+tape.nvar
                        @inbounds ci_num = tape.node_idx_to_number[ci_idx]
                        for j1=j0+1:trlen
                            @inbounds cii_idx = tr[j1]+tape.nvar
                            @inbounds cii_num = tape.node_idx_to_number[cii_idx]
                            @inbounds w = adj*imm[k]
                            update3(tape,cii_num,ci_idx,w)
                            k+=1
                        end
                    end

                    #updating
                    k = r
                    for j=1:n
                        adjlen += 1
                        # @show j,imm[j],adj
                        @inbounds adjs[adjlen] = imm[k] * adj
                        k+=1
                    end
                    imm_counter = immlen - r + 1
                end
            end

            #update counters
            trlen -= n
            immlen -= imm_counter
        end #end TYPE_O
        # println("++++++++++++++++++++++++++++++++++++")
    end  #end while
    # assert(immlen == 0 && trlen == 0)
    # @show tape.eset

    # @show vidx
    nz = one(I)
    for i = 1:tape.nvar
        @inbounds d = tape.bh3[i]
        for v_idx in keys(d)
            # @inbounds (v_idx,w) = lvi[j]
            if(v_idx<=tape.nvar)
                w = d[v_idx]
                @inbounds tape.hess[nz] = w*factor
                nz += 1
            end
        end
    end
    return tape.nzh
end



#Interface function
function hess_structure3{I,V}(tape::Tape{I,V})
    return hess_struct3(tape)
end

function hess_reverse3{I,V}(tape::Tape{I,V},vvals::Vector{V},pvals::Vector{V})
    hess_reverse3(tape,vvals,pvals,1.0)
end

function hess_reverse3{I,V}(tape::Tape{I,V},vvals::Vector{V},pvals::Vector{V}, factor::V)
    # @time forward_pass2_2ord(tape,vvals,pvals)
    # @time reverse_pass2_2ord(tape,factor)
    forward_pass3_2ord(tape,vvals,pvals)
    reverse_pass3_2ord(tape,factor)
end
