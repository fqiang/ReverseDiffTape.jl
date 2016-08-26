#edge pusing algorithm for Hessian reverse AD

function reset_hess2{I,V}(tape::Tape{I,V})
    for i =1:length(tape.bhi)
        @inbounds tape.bhi[i] = Vector{Int}()
        @inbounds tape.bhv[i] = Vector{Float64}()
    end

    fill!(tape.bh_idxes,zero(I))

    tape.h_I = Vector{I}() #hess_I
    tape.h_J = Vector{I}() #hess_J
    tape.hess = Vector{V}() #hess value
    tape.nzh = -one(I)     #hess indicator
end

function prepare_reeval_hess2{I,V}(tape::Tape{I,V})
    fill!(tape.bh_idxes,zero(I))
end

@inline function push_edge2{I,V}(tape::Tape{I,V},to::I,from::I)
    # @show "push_edge2 - ",to," <--- ", from
    @inbounds push!(tape.bhi[to],from)
    @inbounds push!(tape.bhv[to],0.0)
    if to <= tape.nvar && from > tape.nvar
        push_edge2(tape,from,to)
    end
end

function hess_struct2{I,V}(tape::Tape{I,V})
    @timing tape.enable_timing_stats tic()

    if(tape.nzh != -1)
        return 
    end

    tt = tape.tt
    tr = tape.tr
    bhi = tape.bhi
    idx = length(tt)
    trlen = length(tr)
    assert(tape.nnode-1 == length(tr))
    
    while (idx > 0)
        @inbounds ntype = tt[idx]
        idx -= 1
        if(ntype == TYPE_P)
            idx -= 3
        elseif(ntype == TYPE_V)
            @inbounds v_id = tt[idx]
            # @inbounds lvi = bh[v_id]
            # @inbounds e = length(lvi)
            # @inbounds s = bh_v_idxes[v_id] + 1

            # for i = s:e
            #     @inbounds p_id = lvi[i].i
            #     if p_id <= tape.nvar
            #         if p_id < v_id 
            #             push!(tape.h_I, v_id)
            #             push!(tape.h_J, p_id)
            #         else
            #             push!(tape.h_J, v_id)
            #             push!(tape.h_I, p_id)
            #         end
            #     end
            # end
            # @inbounds bh_v_idxes[v_id] = e

            idx -= 2 #skip TYPE_V
        elseif(ntype == TYPE_O)
            @inbounds n = tt[idx]
            @inbounds oc = tt[idx-1]

            #pushing
            @inbounds i_id = tt[idx-2]
            idx -= 4
            # @inbounds lvi = bh[i_id]
            @inbounds lvi = bhi[i_id]
            for j = 1:length(lvi)
                # @inbounds p_id = lvi[j].i
                @inbounds p_id = lvi[j]
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
            trlen -= n
        end
    end #end while loop

    @timing tape.enable_timing_stats tape.ep_structure += toq()
    nothing
end

function hess_findnz2{I,V}(tape::Tape{I,V})
    @timing tape.enable_timing_stats tic()
    bhi = tape.bhi
    nvar = tape.nvar
    nz = zero(I)
    for i = 1:nvar
        @inbounds lvi = bhi[i]
        for j = 1:length(lvi)
            @inbounds v_id = lvi[j]
            if v_id <= nvar
                nz += 1
            end
        end
    end
    tape.nzh = nz
    @timing tape.enable_timing_stats tape.ep_findnz += toq()
    nothing
end

function hess_struct2_recovery{I,V}(tape::Tape{I,V}) 
    @timing tape.enable_timing_stats tic()
    resize!(tape.h_I, tape.nzh)
    resize!(tape.h_J, tape.nzh)
    resize!(tape.hess, tape.nzh)

    bhi = tape.bhi
    nvar = tape.nvar
    h_I = tape.h_I
    h_J = tape.h_J
    nz = zero(I)
    
    for i=1:nvar
        # @inbounds lvi = bh[i]
        @inbounds lvi = bhi[i]
        for j = 1:length(lvi)
            # @inbounds v_id = lvi[j].i
            @inbounds v_id = lvi[j]
            if v_id <= nvar
                if v_id <= i  #assert lower trangular
                    nz += 1
                    @inbounds h_I[nz] = i
                    @inbounds h_J[nz] = v_id
                else 
                    nz += 1
                    @inbounds h_I[nz] = v_id
                    @inbounds h_J[nz] = i
                end
            end
        end
    end

    @assert tape.nzh == nz
    @timing tape.enable_timing_stats tape.ep_structure_recovery += toq()
    nothing
end

function forward_pass2_2ord{I,V}(tape::Tape{I,V}, vvals::Array{V,1}, pvals::Array{V,1})
    @timing tape.enable_timing_stats tic()

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
    if length(tape.imm) < immlen
        resize!(tape.imm,immlen) #if not resize memory will be hold in julia , and not claimed by gc
        # @show "ep - imm length $immlen"
    end

    @timing tape.enable_timing_stats tape.ep_forward_time += toq()

    return @inbounds stk[1]
end

@inline function update2(tape,to,from,w)
    # @show "update2 - ", to, "<-- ", from, w
    @inbounds tape.bh_idxes[to] += 1  
    # assert(tape.bhi[to][tape.bh_idxes[to]] == from)
    @inbounds tape.bhv[to][tape.bh_idxes[to]] = w

    if to <= tape.nvar && from > tape.nvar
        update2(tape,from,to,w)
    end
end

function reverse_pass2_2ord{I,V}(tape::Tape{I,V})
    tr = tape.tr
    tt = tape.tt
    bhi = tape.bhi
    bhv = tape.bhv
    bh_idxes = tape.bh_idxes
    idx = length(tt)
    trlen = length(tr)
    imm = tape.imm
    immlen = tape.imm2ordlen
    nnz = zero(I)
    # assert(length(imm) == immlen)  

    adjs = tape.stk
    adjlen = 1
    @inbounds adjs[1] = one(V)  #initialize adjoint = 1

    while(idx > 0)
        # println("++++++++++++++++++++++++++++++++++++")
        # @show idx
        # @show trlen,immlen, adjlen
        @timing tape.enable_timing_stats tic()   #timing
        @timing tape.enable_timing_stats node_id = 0  #timing

        @inbounds ntype = tt[idx]
        idx -= 1
        
        #adjoints
        @inbounds adj = adjs[adjlen]
        adjlen -= 1
        # @show adj

        if ntype == TYPE_P
            @timing tape.enable_timing_stats node_id = tt[idx-1]  #timing
            idx -= 3
        elseif ntype == TYPE_V
            @inbounds v_id = tt[idx]

            # @inbounds lvi = bh[v_id]
            # @inbounds e = bh_idxes[v_id]
            # @inbounds s = bh_v_idxes[v_id] + 1

            # for i = s:e
            #     @inbounds p = lvi[i]
            #     if p.i <= tape.nvar
            #         nnz += 1
            #         @inbounds tape.hess[nnz] = p.w
            #     end
            # end
            # @inbounds bh_v_idxes[v_id] = e

            idx -= 2 #skip TYPE_V
            @timing tape.enable_timing_stats node_id = v_id  #timing
        elseif ntype == TYPE_O 
            @inbounds n = tt[idx]
            @inbounds oc = tt[idx-1]
            @inbounds op_sym = OP[oc]
            @inbounds i_id = tt[idx-2]
            @timing tape.enable_timing_stats node_id = i_id #timing
            idx -= 4 #skip TYPE_O

            @inbounds lvi = bhi[i_id]
            @inbounds lvv = bhv[i_id]
            @inbounds tr0_id = tr[trlen]
            # @show op_sym, lvi
            imm_counter = zero(I)

            if n==1
                #pushing
                @inbounds t0 = imm[immlen-1]
                for j=1:length(lvi)
                    # @inbounds p = lvi[j]
                    # p_id = p.i
                    # w = p.w
                    @inbounds p_id = lvi[j]
                    @inbounds w = lvv[j]
                    if p_id == i_id                     
                        w_bar = t0*t0*w
                        update2(tape,tr0_id,tr0_id,w_bar)     
                    else
                        w_bar = t0*w
                        if tr0_id == p_id
                            w_bar = 2.0*w_bar
                        end
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
                        # @inbounds p = lvi[j]
                        # p_id = p.i
                        # w = p.w
                        @inbounds p_id = lvi[j]
                        @inbounds w = lvv[j]
                    
                        dlw = dl*w
                        drw = dr*w
                        if p_id == i_id
                            update2(tape,li_id,li_id,dl*dlw)
                            update2(tape,tr0_id,li_id, dl*drw)
                            update2(tape,tr0_id,tr0_id, dr*drw)
                        else
                            if li_id == p_id
                                dlw = 2.0*dlw
                            end
                            if tr0_id == p_id
                                drw = 2.0*drw
                            end
                            update2(tape,li_id,p_id,dlw)
                            update2(tape,tr0_id,p_id,drw)
                        end
                    end
                    #creating
                    @inbounds dll = imm[immlen-2]
                    @inbounds dlr = imm[immlen-1]
                    @inbounds drr = imm[immlen]
                    if tr0_id == li_id 
                        update2(tape,tr0_id,li_id,2.0*adj*dlr)
                    else
                        update2(tape,tr0_id,li_id,adj*dlr)
                    end
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
                        # @inbounds p = lvi[j]
                        # p_id = p.i
                        # w = p.w
                        @inbounds p_id = lvi[j]
                        @inbounds w = lvv[j]
                    
                        dlw = dl*w
                        drw = dr*w
                        if p_id == i_id
                            update2(tape,li_id,li_id,dl*dlw)
                            update2(tape,tr0_id,li_id,dl*drw)
                            update2(tape,tr0_id, tr0_id,dr*drw)
                        else 
                            if li_id == p_id
                                dlw = 2.0*dlw
                            end
                            if tr0_id == p_id
                                drw = 2.0*drw
                            end
                            update2(tape,li_id,p_id,dlw)
                            update2(tape,tr0_id,p_id,drw)
                        end
                    end
                    #creating
                    if tr0_id == li_id
                        update2(tape,tr0_id,li_id,2.0*adj)
                    else
                        update2(tape,tr0_id,li_id,adj)
                    end
                    
                    #updating
                    imm_counter = 3
                    adjlen += 1
                    @inbounds adjs[adjlen] = dl*adj
                    adjlen += 1
                    @inbounds adjs[adjlen] = dr*adj
                elseif  op_sym == :-
                    #pushing
                    for j = 1:length(lvi)
                        # @inbounds p = lvi[j]
                        # p_id = p.i
                        # w = p.w
                        @inbounds p_id = lvi[j]
                        @inbounds w = lvv[j]

                        if p_id == i_id
                            update2(tape,li_id,li_id,w)
                            update2(tape,tr0_id,li_id,-1.0*w)
                            update2(tape,tr0_id,tr0_id,w)
                        else
                            dlw = w
                            drw = -1.0*w
                            if li_id == p_id
                                dlw = 2.0*dlw
                            end
                            if tr0_id == p_id
                                drw = 2.0*drw
                            end
                            update2(tape,li_id,p_id,dlw)
                            update2(tape,tr0_id,p_id,drw)
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
                        # @inbounds p = lvi[j]
                        # p_id = p.i
                        # w = p.w
                        @inbounds p_id = lvi[j]
                        @inbounds w = lvv[j]

                        if p_id == i_id
                            update2(tape,li_id,li_id,w)
                            update2(tape,tr0_id,li_id,w)
                            update2(tape,tr0_id,tr0_id,w)
                        else
                            dlw = w
                            drw = w
                            if li_id == p_id
                                dlw = 2.0*w
                            end
                            if tr0_id == p_id
                                drw = 2.0*w
                            end
                            update2(tape,li_id,p_id,dlw)
                            update2(tape,tr0_id,p_id,drw)
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
                        # @inbounds p = lvi[j]
                        # p_id = p.i
                        # w = p.w
                        @inbounds p_id = lvi[j]
                        @inbounds w = lvv[j]

                        dlw = dl*w
                        drw = dr*w
                        if p_id == i_id
                            update2(tape,li_id,li_id,dl*dlw)
                            update2(tape,tr0_id,li_id,dl*drw)
                            update2(tape,tr0_id,tr0_id,dr*drw)
                        else 
                            if li_id == p_id 
                                dlw = 2.0*dlw
                            end
                            if tr0_id == p_id
                                drw = 2.0*drw
                            end
                            update2(tape,li_id,p_id,dlw)
                            update2(tape,tr0_id,p_id,drw)
                        end
                    end
                    #creating
                    @inbounds dll = imm[immlen-2]
                    @inbounds dlr = imm[immlen-1]
                    @inbounds drr = imm[immlen]
                    update2(tape,li_id,li_id,adj*dll)
                    
                    if tr0_id == li_id
                        update2(tape,tr0_id,li_id,2.0*adj*dlr)
                    else
                        update2(tape,tr0_id,li_id,adj*dlr)
                    end
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
                        # @inbounds p = lvi[j]
                        # p_id = p.i
                        # w = p.w
                        @inbounds p_id = lvi[j]
                        @inbounds w = lvv[j]

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
                                w_bar = w
                                if ci_id == p_id
                                    w_bar = 2.0*w
                                end
                                update2(tape,ci_id,p_id,w_bar)
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
                        # @inbounds p = lvi[j]
                        # p_id = p.i
                        # w = p.w
                        p_id = lvi[j]
                        w = lvv[j]

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
                                if ci_id == p_id
                                    w_bar = 2.0*w_bar
                                end
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
                            if cii_id == ci_id 
                                update2(tape,cii_id, ci_id, 2.0*w)
                            else
                                update2(tape,cii_id,ci_id,w)
                            end
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

        @timing tape.enable_timing_stats tape.ep_reverse_times[node_id] += toq()
        # println("++++++++++++++++++++++++++++++++++++")
    end  #end while
    nothing
end

# function separate(bh::Vector{Vector{mPair{Int,Float64}}})
#     bhi = Vector{Vector{Int}}(length(bh))
#     bhv = Vector{Vector{Float64}}(length(bh))

#     for i = 1:length(bh)
#         @inbounds lvi = bh[i]
#         @inbounds bhi[i] = Vector{Int}(length(lvi))
#         @inbounds bhv[i] = Vector{Float64}(length(lvi))
#         @inbounds bhii = bhi[i]
#         @inbounds bhvi = bhv[i]
#         for j = 1:length(bhii)
#             @inbounds bhii[j] = lvi[j].i
#             @inbounds bhvi[j] = lvi[j].w
#         end
#     end
#     return bhi, bhv
# end

# function recovery2(tape::Tape{Int,Float64}, bhi::Vector{Vector{Int}}, bhv::Vector{Vector{Float64}}, h::Vector{Float64}, nvar::Int, factor::Float64)
#     @timing tape.enable_timing_stats tic()
#     nnz = zero(Int)
#     for i = 1:nvar
#         @inbounds bhii = bhi[i]
#         @inbounds bhvi = bhv[i]
#         for j = 1:length(bhii)
#             @inbounds id = bhii[j]
#             if id <= nvar
#                 nnz += 1
#                 @inbounds h[nnz] = factor*bhvi[j]
#             end
#         end
#     end 
#     @assert (length(h) == nnz) length(h),nnz
#     @timing tape.enable_timing_stats tape.ep_recovery_time += toq()
# end

@inline function hess_recovery2(tape::Tape{Int,Float64}, factor::Float64)
    @timing tape.enable_timing_stats tic()
    h = tape.hess
    nvar = tape.nvar
    bhi = tape.bhi
    bhv = tape.bhv

    nnz = zero(Int)
    for i = 1:nvar
        @inbounds bhii = bhi[i]
        @inbounds bhvi = bhv[i]
        for j = 1:length(bhii)
            @inbounds id = bhii[j]
            if id <= nvar
                nnz += 1
                @inbounds h[nnz] = factor*bhvi[j]
            end
        end
    end 
    @assert (length(h) == nnz) length(h),nnz
    @timing tape.enable_timing_stats tape.ep_recovery_time += toq()
    nothing
end

# function handle_sum_node(tape::Tape{Int,Float64}, i_id::Int, n::Int, trlen::Int)
#     @inbounds lvi = tape.bhi[i_id]
#     @inbounds lvv = tape.bhv[i_id]
#     tr = tape.tr

#     for j = 1:length(lvi)
#         @inbounds p_id = lvi[j]
#         @inbounds w = lvv[j]

#         if p_id == i_id  
#             for j0=trlen-n+1:trlen
#                 @inbounds ci_id = tr[j0]
#                 update2(tape,ci_id,ci_id,w)
#                 for j1=j0+1:trlen
#                     @inbounds cii_id = tr[j1]
#                     update2(tape,cii_id,ci_id,w)
#                 end #j1 +=1
#             end #j0+=1
#         else #p_idx != i_idx
#             for k=trlen-n+1:trlen
#                 @inbounds ci_id = tr[k]
#                 w_bar = w
#                 if ci_id == p_id
#                     w_bar = 2.0*w
#                 end
#                 update2(tape,ci_id,p_id,w_bar)
#             end
#         end
#     end
# end

#Interface function
function hess_structure2{I,V}(tape::Tape{I,V})
    @assert tape.bh_type == 1
    hess_struct2(tape)
    hess_findnz2(tape)
    hess_struct2_recovery(tape)
end

function hess_reverse2{I,V}(tape::Tape{I,V},vvals::Vector{V},pvals::Vector{V})
    hess_reverse2(tape,vvals,pvals,1.0)
end

function hess_reverse2{I,V}(tape::Tape{I,V},vvals::Vector{V},pvals::Vector{V}, factor::V)
    @timing tape.enable_timing_stats tape.ep_n += 1

    forward_pass2_2ord(tape,vvals,pvals)
    reverse_pass2_2ord(tape)
    # @time recovery2(tape, tape.bhi,tape.bhv,tape.hess,tape.nvar, factor)
    hess_recovery2(tape,factor)
end
