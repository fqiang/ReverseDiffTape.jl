#edge pusing algorithm for Hessian reverse AD

function reset_hess2(tape)
    assert(length(tape.bh) == length(tape.bh_idxes) == tape.nnode+tape.nvar)  
    for i=1:tape.nnode+tape.nvar 
        @inbounds tape.bh[i] = Vector{mPair{Int,Float64}}()
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

@inline function push_edge2(tape,to,from)
    # @show "push_edge2 - ",to," <--- ", from
    @inbounds push!(tape.bh[to],mPair{Int,Float64}(from,0.0)) 
end

function hess_struct2{I,V}(tape::Tape{I,V})
    if(tape.nzh != -1)
        return tape.nzh
    end

    tt = tape.tt
    tr = tape.tr
    idx = length(tt)
    trlen = length(tr)
    assert(tape.nnode-1 == length(tr))
    
    while (idx > 0)
        @inbounds ntype = tt[idx]
        idx -= 1
        if(ntype == TYPE_P)
            idx -= 3
        elseif(ntype == TYPE_V)
            @inbounds v_idx = tt[idx]
            idx -= 1
            @inbounds i_id = tt[idx]
            idx -= 1
            # @assert i_id!=0 && v_idx <= tape.nvar && i_id > tape.nvar
            idx -= 1
            @inbounds lvi = tape.bh[i_id]
            for j = 1:length(lvi)
                @inbounds p_id = lvi[j].i
                if p_id == i_id
                    push_edge2(tape,v_idx,v_idx)
                else
                    push_edge2(tape,p_id,v_idx)
                end
            end
        elseif(ntype == TYPE_O)
            @inbounds n = tt[idx]
            idx -= 1
            @inbounds oc = tt[idx]
            idx -= 1

            #pushing
            @inbounds i_id = tt[idx]
            idx -= 1
            # assert(i_id!=0 && i_id > tape.nvar)
            idx -= 1
            @inbounds lvi = tape.bh[i_id]
            for j = 1:length(lvi)
                @inbounds p_id = lvi[j].i
                # @show p_id
                if(p_id == i_id)
                    for j0=trlen-n+1:trlen 
                        @inbounds ci_id = tr[j0]
                        push_edge2(tape,ci_id, ci_id)
                        for j1=j0+1:trlen
                            @inbounds cii_id = tr[j1]
                            push_edge2(tape,cii_id, ci_id)
                        end
                    end
                else  #when i_id != p_id
                    for j0 = trlen -n + 1:trlen
                        @inbounds ci_id = tr[j0]
                        push_edge2(tape,ci_id, p_id)
                    end
                end
            end
           
            #creating 
            @inbounds op_sym = OP[oc]
            # @show op_sym
            if (op_sym ==:+ || op_sym == :-)
                #zeros
            elseif (n == 1 && op_sym != :-) #1-ary operator
                @inbounds ci_id = tr[trlen]
                push_edge2(tape,ci_id, ci_id)
            elseif (op_sym == :*)
                # @show "times ", n
                for j0 = trlen -n + 1:trlen
                    @inbounds ci_id = tr[j0]
                    for j1=j0+1:trlen
                        @inbounds cii_id = tr[j1]
                        push_edge2(tape,cii_id, ci_id)
                    end
                end 
            elseif (op_sym == :/) # binary operator /
                assert(n==2)
                @inbounds ri_id = tr[trlen]
                @inbounds li_id = tr[trlen-1]
                push_edge2(tape,ri_id,li_id)
                push_edge2(tape,ri_id,ri_id)
            else # other binary
                assert(n==2)
                @inbounds ri_id = tr[trlen]
                @inbounds li_id = tr[trlen-1]
                push_edge2(tape,li_id, li_id)
                push_edge2(tape,ri_id, li_id)
                push_edge2(tape,ri_id, ri_id)
            end
            # @show "after creating", tape.live_vars
            trlen -= n
        end
    end #end while loop

    for i=1:tape.nvar
        @inbounds lvi = tape.bh[i] 
        # @show i, length(lvi)
        for j = 1:length(lvi)
            # @inbounds (v_idx,w) = lvi[j]
            @inbounds v_idx = lvi[j].i
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


    # @show tape.t_push_edge
    # writedlm(open("log4000_1_bh_length.txt","w"),tape.bh_length)
    return tape.nzh
end

function forward_pass2_2ord{I,V}(tape::Tape{I,V}, vvals::Array{V,1}, pvals::Array{V,1})
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
            idx += 1 #skip ID
            stklen += 1
            @inbounds stk[stklen] = pvals[tt[idx]]
            idx += 2 #skip TYPE_P
        elseif(ntype == TYPE_V)
            idx += 1 #skip ID
            stklen += 1
            @inbounds stk[stklen] = vvals[tt[idx]]
            idx += 2 #skip TYPE_V
        elseif(ntype == TYPE_O)
            idx += 1 #skip ID
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
        # else 
        #     @assert false
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

@inline function update2(tape,to,from,w)
    # @show "update2 - ", to, "<-- ", from, w
    # assert(tape.bh[to][tape.bh_idxes[to]].i == from)
    @inbounds tape.bh_idxes[to] += 1  
    @inbounds tape.bh[to][tape.bh_idxes[to]].w = w
end

function reverse_pass2_2ord{I,V}(tape::Tape{I,V}, factor::V)
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
            idx -= 3
        elseif(ntype == TYPE_V)
            @inbounds v_idx = tt[idx]
            idx -= 1
            # i_idx = idx + tape.nvar
            @inbounds i_id = tt[idx]
            idx -= 1
            idx -= 1 #skip TYPE_V
            @inbounds lvi = tape.bh[i_id]
            for j = 1:length(lvi)
                @inbounds p = lvi[j]
                p_id = p.i
                w = p.w
                if p_id == i_id
                    update2(tape,v_idx,v_idx,w)
                elseif p_id == v_idx  #p belong to i's child
                    update2(tape,v_idx,v_idx,2.0*w)
                else
                    update2(tape,p_id,v_idx,w)
                end
            end
        elseif(ntype == TYPE_O)
            @inbounds n = tt[idx]
            idx -= 1
            @inbounds oc = tt[idx]
            @inbounds op_sym = OP[oc]
            idx -= 1

            @inbounds i_id = tt[idx]
            idx -= 1
            idx -= 1 #skip TYPE_O
            @inbounds lvi = tape.bh[i_id]
            @inbounds tr0_id = tr[trlen]

            imm_counter = zero(I)

            if n==1
                #pushing
                @inbounds t0 = imm[immlen-1]
                for j=1:length(lvi)
                    @inbounds p = lvi[j]
                    p_id = p.i
                    w = p.w
                    if p_id == i_id                     
                        w_bar = t0*t0*w
                        update2(tape,tr0_id,tr0_id,w_bar)     
                    else
                        w_bar = t0*w
                        update2(tape,tr0_id,p_id,w_bar)
                    end
                end
                #creating
                if(op_sym != :-)
                    @inbounds w_bar = adj*imm[immlen]
                    update2(tape,tr0_id,tr0_id,w_bar)
                end

                #updating 
                adjlen += 1
                @inbounds adjs[adjlen] = t0*adj
                imm_counter = 2
            elseif n == 2
                @inbounds li_id = tr[trlen-1]

                if op_sym == :/
                    #pushing
                    @inbounds dl = imm[immlen-4]
                    @inbounds dr = imm[immlen-3]
                    
                    for j = 1:length(lvi)
                        @inbounds p = lvi[j]
                        p_id = p.i
                        w = p.w
                    
                        dlw = dl*w
                        drw = dr*w
                        if p_id == i_id
                            update2(tape,li_id,li_id,dl*dlw)
                            update2(tape,tr0_id,li_id, dl*drw)
                            update2(tape,tr0_id,tr0_id, dr*drw)
                        else 
                            update2(tape,li_id,p_id,dlw)
                            update2(tape,tr0_id,p_id,drw)
                        end
                    end
                    #creating
                    @inbounds dll = imm[immlen-2]
                    @inbounds dlr = imm[immlen-1]
                    @inbounds drr = imm[immlen]
                    update2(tape,tr0_id,li_id,adj*dlr)
                    update2(tape,tr0_id,tr0_id,adj*drr)
                    
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
                    
                    for j = 1:length(lvi)
                        @inbounds p = lvi[j]
                        # p_idx = p.i
                        p_id = p.i
                        w = p.w
                    
                        dlw = dl*w
                        drw = dr*w
                        if p_id == i_id
                            update2(tape,li_id,li_id,dl*dlw)
                            update2(tape,tr0_id,li_id,dl*drw)
                            update2(tape,tr0_id, tr0_id,dr*drw)
                        else 
                            update2(tape,li_id,p_id,dlw)
                            update2(tape,tr0_id,p_id,drw)
                        end
                    end
                    #creating
                    update2(tape,tr0_id,li_id,adj)
                    
                    #updating
                    imm_counter = 3
                    adjlen += 1
                    @inbounds adjs[adjlen] = dl*adj
                    adjlen += 1
                    @inbounds adjs[adjlen] = dr*adj
                elseif  op_sym == :-
                    #pushing
                    for j = 1:length(lvi)
                        @inbounds p = lvi[j]
                        p_id = p.i
                        w = p.w
                        if p_id == i_id
                            update2(tape,li_id,li_id,w)
                            update2(tape,tr0_id,li_id,-1.0*w)
                            update2(tape,tr0_id,tr0_id,w)
                        else
                            update2(tape,li_id,p_id,w)
                            update2(tape,tr0_id,p_id,-1.0*w)
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
                    for j = 1:length(lvi)
                        @inbounds p = lvi[j]
                        # p_idx = p.i
                        p_id = p.i
                        w = p.w
                        if p_id == i_id
                            update2(tape,li_id,li_id,w)
                            update2(tape,tr0_id,li_id,w)
                            update2(tape,tr0_id,tr0_id,w)
                        else
                            update2(tape,li_id,p_id,w)
                            update2(tape,tr0_id,p_id,w)
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
                    
                    for j = 1:length(lvi)
                        @inbounds p = lvi[j]
                        p_id = p.i
                        w = p.w
                    
                        dlw = dl*w
                        drw = dr*w
                        if p_id == i_id
                            update2(tape,li_id,li_id,dl*dlw)
                            update2(tape,tr0_id,li_id,dl*drw)
                            update2(tape,tr0_id,tr0_id,dr*drw)
                        else 
                            update2(tape,li_id,p_id,dlw)
                            update2(tape,tr0_id,p_id,drw)
                        end
                    end
                    #creating
                    @inbounds dll = imm[immlen-2]
                    @inbounds dlr = imm[immlen-1]
                    @inbounds drr = imm[immlen]
                    update2(tape,li_id,li_id,adj*dll)
                    update2(tape,tr0_id,li_id,adj*dlr)
                    update2(tape,tr0_id,tr0_id,adj*drr)

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
                    for j = 1:length(lvi)
                        @inbounds p = lvi[j]
                        p_id = p.i
                        w = p.w
                        if p_id == i_id  
                            for j0=trlen-n+1:trlen
                                @inbounds ci_id = tr[j0]
                                update2(tape,ci_id,ci_id,w)
                                for j1=j0+1:trlen
                                    @inbounds cii_id = tr[j1]
                                    update2(tape,cii_id,ci_id,w)
                                end #j1 +=1
                            end #j0+=1
                        else #p_idx != i_idx
                            for k=trlen-n+1:trlen
                                @inbounds ci_id = tr[k]
                                update2(tape,ci_id,p_id,w)
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
                    for j = 1:length(lvi)
                        @inbounds p = lvi[j]
                        p_id = p.i
                        w = p.w
                        
                        k = r
                        if p_id == i_id
                            for j0=trlen-n+1:trlen
                                @inbounds ci_id = tr[j0]
                                @inbounds t0 = imm[k]
                                t1 = t0*w
                                @inbounds w_bar0 = t0 * t1
                                update2(tape,ci_id,ci_id,w_bar0)
                                k0 = k + 1
                                for j1=j0+1:trlen
                                    @inbounds cii_id = tr[j1]
                                    @inbounds w_bar1 = t1*imm[k0]
                                    update2(tape,cii_id,ci_id,w_bar1)
                                    k0 += 1
                                end
                                k += 1
                            end
                        else #p_idx != i_idx 
                            for j=trlen -n+1:trlen
                                @inbounds ci_id = tr[j]
                                @inbounds w_bar = imm[k] * w
                                update2(tape,ci_id,p_id,w_bar)
                                k += 1
                            end
                        end
                    end
                    #creating
                    k = r + n
                    for j0=trlen-n+1:trlen
                        @inbounds ci_id = tr[j0]
                        for j1=j0+1:trlen
                            @inbounds cii_id = tr[j1]
                            @inbounds w = adj*imm[k]
                            update2(tape,cii_id,ci_id,w)
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
        # else
        #     @assert false
        end #end TYPE_O

        # println("++++++++++++++++++++++++++++++++++++")
    end  #end while
    
    # tic()
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
    # t_recover = toq()
    # @show t_recover
    return tape.nzh
end



#Interface function
function hess_structure2{I,V}(tape::Tape{I,V})
    # @time hess_struct2(tape,0)
    nz = hess_struct2(tape)
    return nz
end

function hess_reverse2{I,V}(tape::Tape{I,V},vvals::Vector{V},pvals::Vector{V})
    hess_reverse2(tape,vvals,pvals,1.0)
end

function hess_reverse2{I,V}(tape::Tape{I,V},vvals::Vector{V},pvals::Vector{V}, factor::V)
    # @time forward_pass2_2ord(tape,vvals,pvals)
    # @time reverse_pass2_2ord(tape,factor)

    forward_pass2_2ord(tape,vvals,pvals)
    reverse_pass2_2ord(tape,factor)
end
