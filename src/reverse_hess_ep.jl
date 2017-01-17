#edge pusing algorithm for Hessian reverse AD

function reset_hess{I,V}(tape::Tape{I,V})
    for i = 1:length(tape.bhi)
        @inbounds tape.bhi[i] = Vector{Int}()
        @inbounds tape.bhv[i] = Vector{Float64}()
    end
    fill!(tape.bh_idxes,zero(I))
    tape.nzh = -one(I)     #hess indicator
    nothing
end

@inline function push_edge{I}(bhi::Vector{Vector{I}},to::I,from::I, nvar::I)
    # @show "push_edge - ",to," <--- ", from
    @inbounds push!(bhi[to],from)
    if to <= nvar && from > nvar
        push_edge(bhi,from,to,nvar)
    end
end

function hess_struct{I}(ttstart::I, ttend::I, tt::Vector{I}, trend::I, tr::Vector{I}, nvar::I, bhi::Vector{Vector{I}})
    #@timing tape.enable_timing_stats tic()
    idx = ttend
    trlen = trend - 1
    
    @inbounds while idx > ttstart
        @inbounds ntype = tt[idx]
        # @show ntype, trlen, idx
        if ntype == TYPE_P
            #@assert false
        elseif ntype == TYPE_V
            idx -= 3 
        elseif ntype == TYPE_O
            @inbounds n = tt[idx-1]
            @inbounds oc = tt[idx-2]

            #pushing
            @inbounds i_id = tt[idx-3]
            @inbounds lvi = bhi[i_id]
            for j = 1:length(lvi)
                # @inbounds p_id = lvi[j].i
                @inbounds p_id = lvi[j]
                if p_id == i_id
                    for j0=trlen-n+1:trlen 
                        @inbounds ci_id = tr[j0]
                        push_edge(bhi,ci_id, ci_id, nvar)
                        for j1=j0+1:trlen
                            @inbounds cii_id = tr[j1]
                            push_edge(bhi,cii_id, ci_id, nvar)
                        end
                    end
                else  #when i_id != p_id
                    for j0 = trlen -n + 1:trlen
                        @inbounds ci_id = tr[j0]
                        push_edge(bhi,ci_id, p_id, nvar)
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
                push_edge(bhi,ci_id, ci_id, nvar)
            elseif (op_sym == :*)
                # @show "times ", n
                for j0 = trlen -n + 1:trlen
                    @inbounds ci_id = tr[j0]
                    for j1=j0+1:trlen
                        @inbounds cii_id = tr[j1]
                        push_edge(bhi,cii_id, ci_id, nvar)
                    end
                end 
            elseif op_sym == :/ # binary operator /
                assert(n==2)
                @inbounds ri_id = tr[trlen]
                @inbounds li_id = tr[trlen-1]
                push_edge(bhi,ri_id,li_id, nvar)
                push_edge(bhi,ri_id,ri_id, nvar)
            else # other binary
                assert(n==2)
                @inbounds ri_id = tr[trlen]
                @inbounds li_id = tr[trlen-1]
                push_edge(bhi,li_id, li_id, nvar)
                push_edge(bhi,ri_id, li_id, nvar)
                push_edge(bhi,ri_id, ri_id, nvar)
            end
            trlen -= n
            idx -= 5
        elseif ntype == TYPE_O1
            @inbounds oc = tt[idx-2]
            #pushing 
            @inbounds i_id = tt[idx-3]
            @inbounds lvi = bhi[i_id]
            @inbounds c_id = tr[trlen]

            for j = 1:length(lvi)
                @inbounds p_id = lvi[j]
                if p_id == i_id
                    push_edge(bhi,c_id, c_id, nvar)
                else
                    push_edge(bhi,c_id,p_id, nvar)
                end
            end
            #creating
            @inbounds op_sym = OP[oc]
            if op_sym == :/ || op_sym == :^
                push_edge(bhi,c_id,c_id, nvar)
            end
            trlen -= 1
            idx -= 5
        elseif ntype == TYPE_O2
            @inbounds oc = tt[idx-2]
            #pushing 
            @inbounds i_id = tt[idx-3]
            @inbounds lvi = bhi[i_id]
            @inbounds c_id = tr[trlen]

            for j = 1:length(lvi)
                @inbounds p_id = lvi[j]
                if p_id == i_id
                    push_edge(bhi,c_id, c_id, nvar)
                else
                    push_edge(bhi,c_id,p_id, nvar)
                end
            end
            #creating
            @inbounds op_sym = OP[oc]
            if op_sym == :^
                push_edge(bhi,c_id,c_id, nvar)
            end
            trlen -= 1
            idx -= 5
        elseif ntype == TYPE_O3
            @inbounds oc = tt[idx-3]
            @inbounds i_id = tt[idx-4]
            @inbounds vidx = tt[idx-1]

            @inbounds lvi = bhi[i_id]
            #pushing 
            for j = 1:length(lvi)
                @inbounds p_id = lvi[j]
                if p_id == i_id
                    push_edge(bhi,vidx, vidx, nvar)
                else
                    push_edge(bhi,vidx,p_id, nvar)
                end
            end
            #creating
            @inbounds op_sym = OP[oc]
            if op_sym == :/ || op_sym == :^
                push_edge(bhi,vidx,vidx, nvar)
            end
            idx -= 6
        elseif ntype == TYPE_O4
            @inbounds oc = tt[idx-3]
            @inbounds i_id = tt[idx-4]
            @inbounds vidx = tt[idx-1]
            
            @inbounds lvi = bhi[i_id]
            #pushing 
            for j = 1:length(lvi)
                @inbounds p_id = lvi[j]
                if p_id == i_id
                    push_edge(bhi,vidx, vidx, nvar)
                else
                    push_edge(bhi,vidx,p_id, nvar)
                end
            end
            #creating
            @inbounds op_sym = OP[oc]
            if op_sym == :^
                push_edge(bhi,vidx,vidx, nvar)
            end
            idx -= 6
        elseif ntype == TYPE_O5
            @inbounds oc = tt[idx-2]
            @inbounds i_id = tt[idx-3]
            @inbounds vidx = tt[idx-1]

            @inbounds lvi = bhi[i_id]
            @inbounds c_id = tr[trlen]
            #pusing 
            for j = 1:length(lvi)
                @inbounds p_id = lvi[j]
                if p_id == i_id
                    push_edge(bhi, vidx, vidx, nvar)
                    push_edge(bhi, c_id, vidx, nvar)
                    push_edge(bhi, c_id, c_id, nvar)
                else
                    push_edge(bhi, vidx, p_id, nvar)
                    push_edge(bhi, c_id, p_id, nvar)
                end
            end
            #creating
            @inbounds op_sym = OP[oc]
            if op_sym ==:*
                push_edge(bhi,c_id,vidx, nvar)
            elseif op_sym ==:/
                push_edge(bhi,c_id,vidx, nvar)
                push_edge(bhi,c_id,c_id, nvar)
            elseif op_sym ==:^
                push_edge(bhi,vidx,vidx, nvar)
                push_edge(bhi,c_id,vidx, nvar)
                push_edge(bhi,c_id,c_id, nvar)
            else
                #@assert op_sym == :+ || op_sym ==:-
            end
            trlen -= 1
            idx -= 5
        elseif ntype == TYPE_O6
            @inbounds oc = tt[idx-2]
            @inbounds i_id = tt[idx-3]
            @inbounds vidx = tt[idx-1]

            @inbounds lvi = bhi[i_id]
            @inbounds c_id = tr[trlen]
            #pusing 
            for j = 1:length(lvi)
                @inbounds p_id = lvi[j]
                if p_id == i_id
                    push_edge(bhi, c_id, c_id, nvar)
                    push_edge(bhi, vidx, c_id, nvar)
                    push_edge(bhi, vidx, vidx, nvar)
                else
                    push_edge(bhi, c_id, p_id, nvar)
                    push_edge(bhi, vidx, p_id, nvar)
                end
            end
            #creating
            @inbounds op_sym = OP[oc]
            if op_sym ==:*
                push_edge(bhi,vidx,c_id, nvar)
            elseif op_sym ==:/
                push_edge(bhi,vidx,c_id, nvar)
                push_edge(bhi,vidx,vidx, nvar)
            elseif op_sym ==:^
                push_edge(bhi,c_id,c_id, nvar)
                push_edge(bhi,vidx,c_id, nvar)
                push_edge(bhi,vidx,vidx, nvar)
            else
                #@assert op_sym == :+ || op_sym ==:-
            end
            trlen -= 1
            idx -= 5

        elseif ntype == TYPE_O7
            @inbounds oc = tt[idx-2]
            @inbounds i_id = tt[idx-3]
            @inbounds vidx = tt[idx-1]

            @inbounds lvi = bhi[i_id]
            #pushing 
            for j = 1:length(lvi)
                @inbounds p_id = lvi[j]
                if p_id == i_id 
                    push_edge(bhi,vidx,vidx, nvar)
                else
                    push_edge(bhi,vidx,p_id, nvar)
                end
            end
            #creating
            @inbounds op_sym = OP[oc]
            if op_sym != :-
                push_edge(bhi,vidx,vidx, nvar)
            end
            idx -= 5
        else 
            #@assert false
        end
    end #end while loop
    
    @assert trlen == 0 && idx == (ttstart-1)
    
    #@timing tape.enable_timing_stats tape.ep_structure += toq()
    nothing
end

@inline function hess_findnz{I,V}(bhi::Vector{Vector{I}}, idxmap::Vector{I}, HH::Vector{V}, nvar::I)
    #@timing tape.enable_timing_stats tic()
    nnz = zero(I)
    for i = 1:nvar
        @inbounds lvi = bhi[i]
        for j = 1:length(lvi)
            @inbounds v_id = lvi[j]
            if v_id <= nvar
                nnz += 1
            end
        end
    end
    resize!(idxmap, nnz)
    fill!(idxmap,0)
    resize!(HH, nnz)
    #@timing tape.enable_timing_stats tape.ep_findnz += toq()
    return nnz
end

function hess_struct_recovery{I}(bhi::Vector{Vector{I}}, idxmap::Vector{I}, nvar::I, nnz_original::I, h_I::Vector{I}, h_J::Vector{I}) 
    #@timing tape.enable_timing_stats tic()
    @assert length(h_I) == length(h_J) 
    @assert length(idxmap) == nnz_original
    
    hh_I = zeros(I,nnz_original)
    hh_J = zeros(I,nnz_original)
    start = one(I)
    for i=1:nvar
        @inbounds lvi = bhi[i]
        for j = 1:length(lvi)
            @inbounds v_id = lvi[j]
            if v_id <= nvar
                if v_id <= i  #assert lower trangular
                    @inbounds hh_I[start] = i
                    @inbounds hh_J[start] = v_id
                else 
                    @inbounds hh_I[start] = v_id
                    @inbounds hh_J[start] = i
                end
                start += 1
            end
        end
    end
    @assert (start - 1) == nnz_original "$start $(nnz_original)"

    #postprocess non-repeated
    mergedindices = zeros(I,length(hh_I))
    @inbounds mergednnz = [0]
    mergedmap = zeros(I,length(hh_I))
    function combine(idx1,idx2)
        @inbounds if mergedmap[idx1] == 0 && mergedmap[idx2] != 0
            @inbounds mergednnz[1] += 1
            @inbounds mergedmap[idx1] = idx2
            @inbounds mergedindices[mergednnz[1]] = idx1
            return idx2
        else
            @inbounds @assert mergedmap[idx2] == 0
            @inbounds mergednnz[1] += 1
            @inbounds mergedmap[idx2] = idx1
            @inbounds mergedindices[mergednnz[1]] = idx2
            return idx1
        end
    end
    Hmat = sparse(hh_I, hh_J, [i for i in 1:length(hh_I)], nvar, nvar, combine)
    nnz_non_repeat = length(Hmat.nzval)
    pre = length(h_I)
    start = pre + 1
    append!(h_I, zeros(nnz_non_repeat))
    append!(h_J, zeros(nnz_non_repeat))

    for row in 1:nvar
        @inbounds for pos in Hmat.colptr[row]:(Hmat.colptr[row+1]-1)
            @inbounds col = Hmat.rowval[pos]
            @inbounds origidx = Hmat.nzval[pos] # this is the original index of this element
            @inbounds idxmap[origidx] = pos
            @inbounds h_I[start] = row
            @inbounds h_J[start] = col
            start += 1
        end
    end

    @inbounds for k in 1:mergednnz[1]
        @inbounds origidx = mergedindices[k]
        @inbounds mergedwith = mergedmap[origidx]
        @inbounds @assert idxmap[origidx] == 0
        @inbounds @assert idxmap[mergedwith] != 0
        @inbounds idxmap[origidx] = idxmap[mergedwith]
    end

    for i in 1:nnz_original
         @inbounds @assert idxmap[i] != 0
    end

    # @assert start - pre  - 1 == nz "$start $pre $(nz)"
    #@timing tape.enable_timing_stats tape.ep_structure_recovery += toq()

    nothing
end



function forward_pass_2ord{I,V}(ttstart::I, ttend::I, tt::Vector{I}, vvals::Array{V,1}, pvals::Array{V,1}, stk::Vector{V}, vals::Vector{V})
    #@timing tape.enable_timing_stats tic()

    idx = ttstart
    vallen = zero(I)
    stklen = zero(I)

    @inbounds while idx <= ttend
        ntype = tt[idx]
        # @show ntype
        # if ntype == TYPE_P
        #     #@assert false
        # else
        if ntype == TYPE_V
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
            vals[vallen] = stk[stklen]
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
        #     #@assert false
        end
    end
    @assert stklen == 1 && idx == (ttend + 1)
    vallen += 1
    @inbounds vals[vallen] = stk[stklen]
    #@timing tape.enable_timing_stats tape.ep_forward_time += toq()
    return vallen
end

@inline function update(bhv,to,from,w,bh_idxes,nvar)
    # @show "update - ", to, "<-- ", from, w
    @inbounds bh_idxes[to] += 1  
    # @assert tape.bhi[to][tape.bh_idxes[to]] == from
    @inbounds bhv[to][bh_idxes[to]] = w

    if to <= nvar && from > nvar
        update(bhv,from,to,w,bh_idxes,nvar)
    end
end

function reverse_pass_2ord{I,V}(ttstart::I, ttend::I, tt::Vector{I}, trend::I, tr::Vector{I}, 
    bhi::Vector{Vector{I}}, bhv::Vector{Vector{V}}, bh_idxes::Vector{I}, nvar::I,
    vvals::Vector{V}, pvals::Vector{V}, vallen::I, vals::Vector{V}, stk::Vector{V}, imm::Vector{V})
    
    idx = ttend
    trlen = trend - 1
    vallen -= 1

    adjs = stk
    adjlen = 1
    @inbounds adjs[adjlen] = one(V)  #initialize adjoint = 1

    while idx > ttstart
        #@timing tape.enable_timing_stats tic()   #timing
        #@timing tape.enable_timing_stats node_id = 0  #timing

        @inbounds ntype = tt[idx]
        @inbounds adj = adjs[adjlen]
        
        # @show ntype,trlen,idx
        
        # if ntype == TYPE_P
        #     #@assert false
        # else
        if ntype == TYPE_V
            #@timing tape.enable_timing_stats node_id = tt[idx-1]  #timing
            adjlen -= 1
            idx -= 3 #skip TYPE_V
            
        elseif ntype == TYPE_O 
            @inbounds n = tt[idx-1]
            @inbounds oc = tt[idx-2]

            @inbounds i_id = tt[idx-3]
            @inbounds lvi = bhi[i_id]
            @inbounds lvv = bhv[i_id]

            @inbounds op_sym = OP[oc]

            if n == 1
                @inbounds eval_2ord(OP[oc],vals[vallen],imm)
                @inbounds ci_id = tr[trlen]
                #pushing
                for j = 1:length(lvi)
                    @inbounds p_id = lvi[j]
                    @inbounds w = lvv[j]
                    if p_id == i_id
                        @inbounds w_bar = imm[1]*imm[1]*w
                        update(bhv,ci_id, ci_id, w_bar,bh_idxes,nvar)
                    else
                        @inbounds w_bar = imm[1]*w
                        update(bhv,ci_id, p_id, ci_id == p_id? 2.0*w_bar:w_bar,bh_idxes,nvar)
                    end
                end
                #creating
                if op_sym!=:-
                    @inbounds w_bar = adj*imm[2]
                    update(bhv,ci_id, ci_id, w_bar,bh_idxes,nvar)
                end
            elseif n == 2
                @inbounds eval_2ord(OP[oc],vals[vallen-1],vals[vallen],imm)
                @inbounds dl = imm[1]
                @inbounds dr = imm[2]
                @inbounds li_id = tr[trlen-1]
                @inbounds ri_id = tr[trlen]

                #pushing
                for j = 1:length(lvi)
                    @inbounds p_id = lvi[j]
                    @inbounds w = lvv[j]
                    dlw = dl*w
                    drw = dr*w
                    if p_id == i_id
                        update(bhv,li_id, li_id, dl*dlw,bh_idxes,nvar)
                        dlrw = dl*drw
                        update(bhv,ri_id, li_id, ri_id == li_id? 2.0*dlrw:dlrw,bh_idxes,nvar)
                        update(bhv,ri_id, ri_id, dr*drw,bh_idxes,nvar)
                    else
                        update(bhv,li_id, p_id, li_id == p_id? 2.0*dlw:dlw,bh_idxes,nvar)
                        update(bhv,ri_id, p_id, ri_id == p_id? 2.0*drw:drw,bh_idxes,nvar)
                    end
                end

                #creating
                if op_sym == :+ || op_sym==:-
                    #nothing
                elseif op_sym == :/
                    #@assert n==2
                    @inbounds dlr = imm[4]
                    @inbounds drr = imm[5]
                    w_bar = adj*dlr
                    update(bhv,ri_id,li_id, ri_id == li_id?2.0*w_bar:w_bar,bh_idxes,nvar)
                    update(bhv,ri_id,ri_id,adj*drr,bh_idxes,nvar)

                elseif op_sym == :*
                    #@assert n==2
                    @inbounds dlr = imm[4]
                    update(bhv,ri_id, li_id,ri_id == li_id? 2.0*adj:adj,bh_idxes,nvar)
                     
                else
                    @assert op_sym == :^
                    @inbounds dll = imm[3]
                    @inbounds dlr = imm[4]
                    @inbounds drr = imm[5]
                    update(bhv,li_id,li_id,adj*dll,bh_idxes,nvar)
                    w_bar = adj*dlr
                    update(bhv,ri_id,li_id,ri_id==li_id?2.0*w_bar:w_bar,bh_idxes,nvar)
                    update(bhv,ri_id,ri_id,adj*drr,bh_idxes,nvar)

                end
            else
                @assert n>2 
                @inbounds eval_2ord(OP[oc],vals,vallen-n+1,vallen,imm)
                #pushing
                if op_sym == :+
                    for j = 1:length(lvi)
                        @inbounds p_id = lvi[j]
                        @inbounds w = lvv[j]
                        if p_id == i_id
                            for j0=trlen-n+1:trlen 
                                @inbounds ci_id = tr[j0]
                                w_bar = w
                                update(bhv, ci_id, ci_id, w_bar,bh_idxes,nvar)
                                for j1=j0+1:trlen
                                    @inbounds cii_id = tr[j1]
                                    update(bhv,cii_id, ci_id,w_bar,bh_idxes,nvar)
                                end
                            end
                        else  
                            for j0 = trlen -n + 1:trlen
                                @inbounds ci_id = tr[j0]
                                w_bar = w
                                update(bhv,ci_id, p_id, ci_id == p_id? 2.0*w_bar:w_bar,bh_idxes,nvar)
                            end
                        end
                    end
                    #creating 
                        #nothing
                else
                    @assert op_sym == :*
                    for j = 1:length(lvi)
                        @inbounds p_id = lvi[j]
                        @inbounds w = lvv[j]
                        if p_id == i_id
                            k=1
                            for j0=trlen-n+1:trlen 
                                @inbounds ci_id = tr[j0]
                                @inbounds t0 = imm[k]
                                t1 = t0 * w
                                w_bar = t0 * t1
                                update(bhv, ci_id, ci_id, w_bar,bh_idxes,nvar)
                                k0 = k + 1
                                for j1=j0+1:trlen
                                    @inbounds cii_id = tr[j1]
                                    @inbounds w_bar1 = t1*imm[k0]
                                    update(bhv,cii_id, ci_id,cii_id==ci_id? 2.0*w_bar1:w_bar1,bh_idxes,nvar)
                                    k0 += 1
                                end
                                k += 1
                            end
                        else  
                            k=1
                            for j0 = trlen -n + 1:trlen
                                @inbounds ci_id = tr[j0]
                                @inbounds w_bar = imm[k]*w
                                update(bhv,ci_id, p_id, ci_id == p_id? 2.0*w_bar:w_bar,bh_idxes,nvar)
                                k += 1
                            end
                        end
                    end
                    #creating
                    # @show "creating", adj
                    k=n+1
                    @simd for j0 = trlen -n + 1:trlen
                        @inbounds ci_id = tr[j0]
                        for j1=j0+1:trlen
                            @inbounds cii_id = tr[j1]
                            @inbounds w_bar = adj*imm[k]
                            update(bhv,cii_id, ci_id, cii_id == ci_id? 2.0*w_bar:w_bar,bh_idxes,nvar)
                            k += 1
                        end
                    end
                end
            end

            #updating
            if op_sym == :+
                @simd for i=1:n
                    @inbounds adjs[adjlen+i-1] = adj
                end
            else
                @simd for i=1:n
                    @inbounds adjs[adjlen+i-1] = imm[i]*adj
                end
            end
            adjlen += n-1
            trlen -= n
            vallen -= n
            idx -= 5
        elseif ntype == TYPE_O1
            @inbounds pval = pvals[tt[idx-1]]
            @inbounds oc = tt[idx-2]
            @inbounds i_id = tt[idx-3]
            @inbounds lvi = bhi[i_id]
            @inbounds lvv = bhv[i_id]

            @inbounds eval_2ord(OP[oc],pval,vals[vallen],imm)
            @inbounds c_id = tr[trlen]
            #pushing
            for j = 1:length(lvi)
                @inbounds p_id = lvi[j]
                @inbounds w = lvv[j]
                @inbounds t0 = imm[2]*w
                if p_id == i_id
                    @inbounds w_bar = imm[2]*t0
                    update(bhv,c_id,c_id,w_bar,bh_idxes,nvar)
                else
                    w_bar = t0
                    update(bhv,c_id,p_id,c_id==p_id?2.0*w_bar:w_bar,bh_idxes,nvar)
                end
            end
            #creating
            @inbounds op_sym = OP[oc]
            if op_sym == :/ || op_sym == :^
                @inbounds update(bhv,c_id,c_id,adj*imm[5],bh_idxes,nvar)
            end

            #updating
            @inbounds adjs[adjlen] = imm[2]*adj
            trlen -= 1
            vallen -= 1
            idx -= 5


        elseif ntype == TYPE_O2
            @inbounds pval = pvals[tt[idx-1]]
            @inbounds oc = tt[idx-2]
            @inbounds i_id = tt[idx-3]
            @inbounds lvi = bhi[i_id]
            @inbounds lvv = bhv[i_id]

            @inbounds eval_2ord(OP[oc],vals[vallen],pval,imm)
            @inbounds c_id = tr[trlen]
            #pushing
            for j = 1:length(lvi)
                @inbounds p_id = lvi[j]
                @inbounds w = lvv[j]
                @inbounds t0 = imm[1]*w
                if p_id == i_id
                    @inbounds w_bar = imm[1]*t0
                    update(bhv,c_id,c_id,w_bar,bh_idxes,nvar)
                else
                    w_bar = t0
                    update(bhv,c_id,p_id,c_id==p_id?2.0*w_bar:w_bar,bh_idxes,nvar)
                end
            end
            #creating
            @inbounds op_sym = OP[oc]
            if op_sym == :^
                @inbounds update(bhv,c_id,c_id,adj*imm[3],bh_idxes,nvar)
            end

            #updating
            @inbounds adjs[adjlen] = imm[1]*adj
            trlen -= 1
            vallen -= 1
            idx -= 5

        elseif ntype == TYPE_O3
            @inbounds vidx = tt[idx-1]
            @inbounds vval = vvals[vidx]
            @inbounds pval = pvals[tt[idx-2]]
            @inbounds oc = tt[idx-3]
            @inbounds i_id = tt[idx-4]
            @inbounds lvi = bhi[i_id]
            @inbounds lvv = bhv[i_id]

            @inbounds eval_2ord(OP[oc],pval,vval,imm)
            #pushing
            for j = 1:length(lvi)
                @inbounds p_id = lvi[j]
                @inbounds w = lvv[j]
                @inbounds t0 = imm[2]*w
                if p_id == i_id
                    @inbounds w_bar = imm[2]*t0
                    update(bhv,vidx,vidx,w_bar,bh_idxes,nvar)
                else
                    w_bar = t0
                    update(bhv,vidx,p_id,vidx==p_id?2.0*w_bar:w_bar,bh_idxes,nvar)
                end
            end
            #creating
            @inbounds op_sym = OP[oc]
            if op_sym == :/ || op_sym == :^
                @inbounds update(bhv,vidx,vidx,adj*imm[5],bh_idxes,nvar)
            end

            #updating
            adjlen -= 1
            idx -= 6

        elseif ntype == TYPE_O4
            @inbounds vidx = tt[idx-1]
            @inbounds vval = vvals[vidx]
            @inbounds pval = pvals[tt[idx-2]]
            @inbounds oc = tt[idx-3]
            @inbounds i_id = tt[idx-4]
            @inbounds lvi = bhi[i_id]
            @inbounds lvv = bhv[i_id]

            @inbounds eval_2ord(OP[oc],vval,pval,imm)
            #pushing
            for j = 1:length(lvi)
                @inbounds p_id = lvi[j]
                @inbounds w = lvv[j]
                @inbounds t0 = imm[1]*w
                if p_id == i_id
                    @inbounds w_bar = imm[1]*t0
                    update(bhv,vidx,vidx,w_bar,bh_idxes,nvar)
                else
                    w_bar = t0
                    update(bhv,vidx,p_id,vidx==p_id?2.0*w_bar:w_bar,bh_idxes,nvar)
                end
            end
            #creating
            @inbounds op_sym = OP[oc]
            if op_sym == :^
                @inbounds update(bhv,vidx,vidx,adj*imm[3],bh_idxes,nvar)
            end

            #updating
            adjlen -= 1
            idx -= 6

        elseif ntype == TYPE_O5
            @inbounds vidx = tt[idx-1]
            @inbounds vval = vvals[vidx]
            @inbounds oc = tt[idx-2]
            @inbounds i_id = tt[idx-3]

            @inbounds lvi = bhi[i_id]
            @inbounds lvv = bhv[i_id]
            
            @inbounds eval_2ord(OP[oc],vval,vals[vallen],imm)
            @inbounds c_id = tr[trlen]
            #pushing
            for j = 1:length(lvi)
                @inbounds p_id = lvi[j]
                @inbounds w = lvv[j]
                @inbounds dlw = imm[1]*w
                @inbounds drw = imm[2]*w
                if p_id == i_id
                    @inbounds update(bhv,vidx,vidx,imm[1]*dlw,bh_idxes,nvar)
                    @inbounds update(bhv,c_id,vidx,imm[1]*drw,bh_idxes,nvar)
                    @inbounds update(bhv,c_id,c_id,imm[2]*drw,bh_idxes,nvar)
                else
                    update(bhv,vidx,p_id,vidx==p_id?2.0*dlw:dlw,bh_idxes,nvar)
                    update(bhv,c_id,p_id,c_id==p_id?2.0*drw:drw,bh_idxes,nvar)
                end
            end
            #creating
            @inbounds op_sym = OP[oc]
            @inbounds dll = imm[3]
            @inbounds dlr = imm[4]
            @inbounds drr = imm[5]
            adlr = adj*dlr
            if op_sym == :*
                update(bhv,c_id,vidx,c_id == vidx?2.0*adlr:adlr,bh_idxes,nvar)
            elseif op_sym == :/
                update(bhv,c_id,vidx,c_id == vidx?2.0*adlr:adlr,bh_idxes,nvar)
                update(bhv,c_id,c_id,adj*drr)
            elseif op_sym == :^
                update(bhv,vidx,vidx,adj*dll)
                update(bhv,c_id,vidx,c_id == vidx?2.0*adlr:adlr,bh_idxes,nvar)
                update(bhv,c_id,c_id,adj*drr,bh_idxes,nvar)
            else
                #@assert op_sym == :+ || op_sym ==:-
            end

            #updating
            # @inbounds adjs[adjlen] = imm[1]*adj
            @inbounds adjs[adjlen] = imm[2]*adj
            trlen -= 1
            vallen -= 1
            idx -= 5

        elseif ntype == TYPE_O6
            @inbounds vidx = tt[idx-1]
            @inbounds vval = vvals[vidx]
            @inbounds oc = tt[idx-2]
            @inbounds i_id = tt[idx-3]

            @inbounds lvi = bhi[i_id]
            @inbounds lvv = bhv[i_id]
            
            @inbounds eval_2ord(OP[oc],vals[vallen],vval,imm)
            @inbounds c_id = tr[trlen]
            #pushing
            for j = 1:length(lvi)
                @inbounds p_id = lvi[j]
                @inbounds w = lvv[j]
                @inbounds dlw = imm[1]*w
                @inbounds drw = imm[2]*w
                if p_id == i_id
                    @inbounds update(bhv,c_id,c_id,imm[1]*dlw,bh_idxes,nvar)
                    @inbounds update(bhv,vidx,c_id,imm[1]*drw,bh_idxes,nvar)
                    @inbounds update(bhv,vidx,vidx,imm[2]*drw,bh_idxes,nvar)
                else
                    update(bhv,c_id,p_id,c_id==p_id?2.0*dlw:dlw,bh_idxes,nvar)
                    update(bhv,vidx,p_id,vidx==p_id?2.0*drw:drw,bh_idxes,nvar)
                end
            end
            #creating
            @inbounds op_sym = OP[oc]
            @inbounds dll = imm[3]
            @inbounds dlr = imm[4]
            @inbounds drr = imm[5]
            adlr = adj*dlr
            if op_sym == :*
                update(bhv,vidx,c_id,vidx == c_id?2.0*adlr:adlr,bh_idxes,nvar)
            elseif op_sym == :/
                update(bhv,vidx,c_id,vidx == c_id?2.0*adlr:adlr,bh_idxes,nvar)
                update(bhv,vidx,vidx,adj*drr,bh_idxes,nvar)
            elseif op_sym == :^
                update(bhv,c_id,c_id,adj*dll,bh_idxes,nvar)
                update(bhv,vidx,c_id,vidx == c_id?2.0*adlr:adlr,bh_idxes,nvar)
                update(bhv,vidx,vidx,adj*drr,bh_idxes,nvar)
            else
                #@assert op_sym == :+ || op_sym ==:-
            end

            #updating
            @inbounds adjs[adjlen] = imm[1]*adj
            trlen -= 1
            vallen -= 1
            idx -= 5

        else
            @assert ntype == TYPE_O7
            @inbounds vidx = tt[idx-1]
            @inbounds oc = tt[idx-2]
            @inbounds i_id = tt[idx-3]
            @inbounds vval = vvals[vidx]

            @inbounds eval_2ord(OP[oc],vval,imm)

            @inbounds lvi = bhi[i_id]
            @inbounds lvv = bhv[i_id]
            #pushing 
            for j = 1:length(lvi)
                @inbounds p_id = lvi[j]
                @inbounds w = lvv[j]
                @inbounds t1 = imm[1]*w
                if p_id == i_id 
                    @inbounds w_bar = imm[1]*t1
                    update(bhv,vidx,vidx,w_bar,bh_idxes,nvar)
                else
                    update(bhv,vidx,p_id,vidx==p_id?2.0*t1:t1,bh_idxes,nvar)
                end
            end
            #creating
            @inbounds op_sym = OP[oc]
            if op_sym != :-
                @inbounds w_bar = adj*imm[2]
                update(bhv,vidx,vidx,w_bar,bh_idxes,nvar)
            end
            adjlen -= 1
            idx -= 5
        end

        #@timing tape.enable_timing_stats tape.ep_reverse_times[node_id] += toq()
        # println("++++++++++++++++++++++++++++++++++++")
    end  #end while

    # @show idx, vallen, adjlen, trlen
    @assert idx == (ttstart-1) && vallen == 0 && adjlen == 0 && trlen == 0
    nothing
end


@inline function hess_recovery{I,V}(bhi::Vector{Vector{I}},bhv::Vector{Vector{V}}, idxmap::Vector{I}, HH::Vector{V}, nvar::I, start::Int, H::Vector{V})
    #@timing tape.enable_timing_stats tic()
    pre = start 
    for i = 1:nvar
        @inbounds lvi = bhi[i]
        @inbounds lvv = bhv[i]
        for j = 1:length(lvi)
            @inbounds v_id = lvi[j]
            if v_id <= nvar
                @inbounds HH[start] = lvv[j]
                start += 1
            end
        end
    end 
    nnz =  start - pre 
    @assert nnz == length(HH)

    fill!(H,0)
    for i = 1:nnz
        @inbounds H[idxmap[i]] += HH[i]
    end

    #@timing tape.enable_timing_stats tape.ep_recovery_time += toq()
    nothing
end


#Interface function
function hess_structure{I,V}(ts::I, te::I,tt::Vector{I}, 
    re::I,tr::Vector{I},
    hs::HessStorage{I,V}, h_I::Vector{I}, h_J::Vector{I})

    hess_struct(ts,te,tt,re,tr,hs.gnvar,hs.bhi)
    @assert length(hs.bhi) == length(hs.bhv) == length(hs.bh_idxes)
    for i=1:length(hs.bhi)
        @inbounds resize!(hs.bhv[i],length(hs.bhi[i]))
    end
    nnz_repeat = hess_findnz(hs.bhi, hs.idxmap, hs.HH, hs.gnvar)
    hess_struct_recovery(hs.bhi, hs.idxmap, hs.gnvar, nnz_repeat, h_I, h_J)
    return length(h_I)
end

function hess_reverse{I,V}(ts::I, te::I, tt::Vector{I}, 
    re::I,tr::Vector{I},
    vvals::Vector{V},pvals::Vector{V}, stk::Vector{V}, vals::Vector{V}, imm::Vector{V}, hs::HessStorage{I,V}, start::I, H::Vector{V})
    #@timing tape.enable_timing_stats tape.ep_n += 1

    vallen = forward_pass_2ord(ts, te, tt, vvals, pvals, stk, vals)
    # @assert vallen == tape.nnode

    fill!(hs.bh_idxes,zero(I))
    reverse_pass_2ord(ts, te, tt, re, tr, hs.bhi, hs.bhv, hs.bh_idxes, hs.gnvar
        ,vvals,pvals, vallen, vals, stk, imm)
    # @time recovery2(tape, tape.bhi,tape.bhv,tape.hess,tape.nvar, factor)
    hess_recovery(hs.bhi, hs.bhv, hs.idxmap, hs.HH, hs.gnvar, start, H)
    nothing
end
