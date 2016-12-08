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

function prepare_reeval_hess{I,V}(tape::Tape{I,V})
    fill!(tape.bh_idxes,zero(I))
    nothing
end

@inline function push_edge{I,V}(tape::Tape{I,V},to::I,from::I)
    # @show "push_edge - ",to," <--- ", from
    @inbounds push!(tape.bhi[to],from)
    @inbounds push!(tape.bhv[to],0.0)
    if to <= tape.nvar && from > tape.nvar
        push_edge(tape,from,to)
    end
end

function hess_struct{I,V}(tape::Tape{I,V})
    #@timing tape.enable_timing_stats tic()

    if(tape.nzh != -1)
        return 
    end

    tt = tape.tt
    tr = tape.tr
    bhi = tape.bhi
    idx = length(tt)
    trlen = length(tr)
    #@assert tape.nnode-1 == length(tr)
    
    @inbounds while idx > 0
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
                        push_edge(tape,ci_id, ci_id)
                        for j1=j0+1:trlen
                            @inbounds cii_id = tr[j1]
                            push_edge(tape,cii_id, ci_id)
                        end
                    end
                else  #when i_id != p_id
                    for j0 = trlen -n + 1:trlen
                        @inbounds ci_id = tr[j0]
                        push_edge(tape,ci_id, p_id)
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
                push_edge(tape,ci_id, ci_id)
            elseif (op_sym == :*)
                # @show "times ", n
                for j0 = trlen -n + 1:trlen
                    @inbounds ci_id = tr[j0]
                    for j1=j0+1:trlen
                        @inbounds cii_id = tr[j1]
                        push_edge(tape,cii_id, ci_id)
                    end
                end 
            elseif op_sym == :/ # binary operator /
                assert(n==2)
                @inbounds ri_id = tr[trlen]
                @inbounds li_id = tr[trlen-1]
                push_edge(tape,ri_id,li_id)
                push_edge(tape,ri_id,ri_id)
            else # other binary
                assert(n==2)
                @inbounds ri_id = tr[trlen]
                @inbounds li_id = tr[trlen-1]
                push_edge(tape,li_id, li_id)
                push_edge(tape,ri_id, li_id)
                push_edge(tape,ri_id, ri_id)
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
                    push_edge(tape,c_id, c_id)
                else
                    push_edge(tape,c_id,p_id)
                end
            end
            #creating
            @inbounds op_sym = OP[oc]
            if op_sym == :/ || op_sym == :^
                push_edge(tape,c_id,c_id)
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
                    push_edge(tape,c_id, c_id)
                else
                    push_edge(tape,c_id,p_id)
                end
            end
            #creating
            @inbounds op_sym = OP[oc]
            if op_sym == :^
                push_edge(tape,c_id,c_id)
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
                    push_edge(tape,vidx, vidx)
                else
                    push_edge(tape,vidx,p_id)
                end
            end
            #creating
            @inbounds op_sym = OP[oc]
            if op_sym == :/ || op_sym == :^
                push_edge(tape,vidx,vidx)
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
                    push_edge(tape,vidx, vidx)
                else
                    push_edge(tape,vidx,p_id)
                end
            end
            #creating
            @inbounds op_sym = OP[oc]
            if op_sym == :^
                push_edge(tape,vidx,vidx)
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
                    push_edge(tape, vidx, vidx)
                    push_edge(tape, c_id, vidx)
                    push_edge(tape, c_id, c_id)
                else
                    push_edge(tape, vidx, p_id)
                    push_edge(tape, c_id, p_id)
                end
            end
            #creating
            @inbounds op_sym = OP[oc]
            if op_sym ==:*
                push_edge(tape,c_id,vidx)
            elseif op_sym ==:/
                push_edge(tape,c_id,vidx)
                push_edge(tape,c_id,c_id)
            elseif op_sym ==:^
                push_edge(tape,vidx,vidx)
                push_edge(tape,c_id,vidx)
                push_edge(tape,c_id,c_id)
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
                    push_edge(tape, c_id, c_id)
                    push_edge(tape, vidx, c_id)
                    push_edge(tape, vidx, vidx)
                else
                    push_edge(tape, c_id, p_id)
                    push_edge(tape, vidx, p_id)
                end
            end
            #creating
            @inbounds op_sym = OP[oc]
            if op_sym ==:*
                push_edge(tape,vidx,c_id)
            elseif op_sym ==:/
                push_edge(tape,vidx,c_id)
                push_edge(tape,vidx,vidx)
            elseif op_sym ==:^
                push_edge(tape,c_id,c_id)
                push_edge(tape,vidx,c_id)
                push_edge(tape,vidx,vidx)
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
                    push_edge(tape,vidx,vidx)
                else
                    push_edge(tape,vidx,p_id)
                end
            end
            #creating
            @inbounds op_sym = OP[oc]
            if op_sym != :-
                push_edge(tape,vidx,vidx)
            end
            idx -= 5
        else 
            #@assert false
        end
    end #end while loop

    #@assert trlen == 0
    #@assert idx == 0

    #@timing tape.enable_timing_stats tape.ep_structure += toq()
    nothing
end

function hess_findnz{I,V}(tape::Tape{I,V})
    #@timing tape.enable_timing_stats tic()
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
    #@timing tape.enable_timing_stats tape.ep_findnz += toq()
    nothing
end

function hess_struct_recovery{I,V}(tape::Tape{I,V}, h_I::Vector{I}, h_J::Vector{I}) 
    #@timing tape.enable_timing_stats tic()
    @assert length(h_I) == length(h_J)
    pre = length(h_I)
    start = pre + 1
    append!(h_I, zeros(tape.nzh))
    append!(h_J, zeros(tape.nzh))
    
    bhi = tape.bhi
    nvar = tape.nvar
    
    for i=1:nvar
        @inbounds lvi = bhi[i]
        for j = 1:length(lvi)
            @inbounds v_id = lvi[j]
            if v_id <= nvar
                if v_id <= i  #assert lower trangular
                    @inbounds h_I[start] = i
                    @inbounds h_J[start] = v_id
                else 
                    @inbounds h_I[start] = v_id
                    @inbounds h_J[start] = i
                end
                start += 1
            end
        end
    end

    # resize!(tape.hess, tape.nzh)

    @assert start - pre  - 1 == tape.nzh "$start $pre $(tape.nzh)"
    #@timing tape.enable_timing_stats tape.ep_structure_recovery += toq()
    nothing
end


function forward_pass_2ord{I,V}(tape::Tape{I,V}, vvals::Array{V,1}, pvals::Array{V,1})
    #@timing tape.enable_timing_stats tic()

    tt = tape.tt
    idx = one(I)
    stk = tape.stk
    vals = tape.vals
    vallen = zero(I)
    stklen = zero(I)

   @inbounds while(idx <= length(tt))
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
    @assert stklen == 1 && vallen + 1 == tape.nnode
    @inbounds ret = vals[vallen + 1] = stk[stklen]
    #@timing tape.enable_timing_stats tape.ep_forward_time += toq()
    return ret
end

@inline function update2(tape,to,from,w)
    # @show "update2 - ", to, "<-- ", from, w
    @inbounds tape.bh_idxes[to] += 1  
    # @assert tape.bhi[to][tape.bh_idxes[to]] == from
    @inbounds tape.bhv[to][tape.bh_idxes[to]] = w

    if to <= tape.nvar && from > tape.nvar
        update2(tape,from,to,w)
    end
end

function reverse_pass_2ord{I,V}(tape::Tape{I,V}, vvals::Vector{V}, pvals::Vector{V})
    tr = tape.tr
    tt = tape.tt
    imm = tape.imm
    vals = tape.vals
    bhi = tape.bhi
    bhv = tape.bhv
    bh_idxes = tape.bh_idxes
    idx = length(tt)
    trlen = length(tr)
    vallen = tape.nnode - 1

    adjs = tape.stk
    adjlen = 1
    @inbounds adjs[adjlen] = one(V)  #initialize adjoint = 1

    while idx > 0
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
                        update2(tape,ci_id, ci_id, w_bar)
                    else
                        @inbounds w_bar = imm[1]*w
                        update2(tape,ci_id, p_id, ci_id == p_id? 2.0*w_bar:w_bar)
                    end
                end
                #creating
                if op_sym!=:-
                    @inbounds w_bar = adj*imm[2]
                    update2(tape,ci_id, ci_id, w_bar)
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
                        update2(tape,li_id, li_id, dl*dlw)
                        dlrw = dl*drw
                        update2(tape,ri_id, li_id, ri_id == li_id? 2.0*dlrw:dlrw)
                        update2(tape,ri_id, ri_id, dr*drw)
                    else
                        update2(tape,li_id, p_id, li_id == p_id? 2.0*dlw:dlw)
                        update2(tape,ri_id, p_id, ri_id == p_id? 2.0*drw:drw)
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
                    update2(tape,ri_id,li_id, ri_id == li_id?2.0*w_bar:w_bar)
                    update2(tape,ri_id,ri_id,adj*drr)

                elseif op_sym == :*
                    #@assert n==2
                    @inbounds dlr = imm[4]
                    update2(tape,ri_id, li_id,ri_id == li_id? 2.0*adj:adj)
                     
                else
                    @assert op_sym == :^
                    @inbounds dll = imm[3]
                    @inbounds dlr = imm[4]
                    @inbounds drr = imm[5]
                    update2(tape,li_id,li_id,adj*dll)
                    w_bar = adj*dlr
                    update2(tape,ri_id,li_id,ri_id==li_id?2.0*w_bar:w_bar)
                    update2(tape,ri_id,ri_id,adj*drr)

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
                                update2(tape, ci_id, ci_id, w_bar)
                                for j1=j0+1:trlen
                                    @inbounds cii_id = tr[j1]
                                    update2(tape,cii_id, ci_id,w_bar)
                                end
                            end
                        else  
                            for j0 = trlen -n + 1:trlen
                                @inbounds ci_id = tr[j0]
                                w_bar = w
                                update2(tape,ci_id, p_id, ci_id == p_id? 2.0*w_bar:w_bar)
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
                                update2(tape, ci_id, ci_id, w_bar)
                                k0 = k + 1
                                for j1=j0+1:trlen
                                    @inbounds cii_id = tr[j1]
                                    @inbounds w_bar1 = t1*imm[k0]
                                    update2(tape,cii_id, ci_id,cii_id==ci_id? 2.0*w_bar1:w_bar1)
                                    k0 += 1
                                end
                                k += 1
                            end
                        else  
                            k=1
                            for j0 = trlen -n + 1:trlen
                                @inbounds ci_id = tr[j0]
                                @inbounds w_bar = imm[k]*w
                                update2(tape,ci_id, p_id, ci_id == p_id? 2.0*w_bar:w_bar)
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
                            update2(tape,cii_id, ci_id, cii_id == ci_id? 2.0*w_bar:w_bar)
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
                    update2(tape,c_id,c_id,w_bar)
                else
                    w_bar = t0
                    update2(tape,c_id,p_id,c_id==p_id?2.0*w_bar:w_bar)
                end
            end
            #creating
            @inbounds op_sym = OP[oc]
            if op_sym == :/ || op_sym == :^
                @inbounds update2(tape,c_id,c_id,adj*imm[5])
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
                    update2(tape,c_id,c_id,w_bar)
                else
                    w_bar = t0
                    update2(tape,c_id,p_id,c_id==p_id?2.0*w_bar:w_bar)
                end
            end
            #creating
            @inbounds op_sym = OP[oc]
            if op_sym == :^
                @inbounds update2(tape,c_id,c_id,adj*imm[3])
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
                    update2(tape,vidx,vidx,w_bar)
                else
                    w_bar = t0
                    update2(tape,vidx,p_id,vidx==p_id?2.0*w_bar:w_bar)
                end
            end
            #creating
            @inbounds op_sym = OP[oc]
            if op_sym == :/ || op_sym == :^
                @inbounds update2(tape,vidx,vidx,adj*imm[5])
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
                    update2(tape,vidx,vidx,w_bar)
                else
                    w_bar = t0
                    update2(tape,vidx,p_id,vidx==p_id?2.0*w_bar:w_bar)
                end
            end
            #creating
            @inbounds op_sym = OP[oc]
            if op_sym == :^
                @inbounds update2(tape,vidx,vidx,adj*imm[3])
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
                    @inbounds update2(tape,vidx,vidx,imm[1]*dlw)
                    @inbounds update2(tape,c_id,vidx,imm[1]*drw)
                    @inbounds update2(tape,c_id,c_id,imm[2]*drw)
                else
                    update2(tape,vidx,p_id,vidx==p_id?2.0*dlw:dlw)
                    update2(tape,c_id,p_id,c_id==p_id?2.0*drw:drw)
                end
            end
            #creating
            @inbounds op_sym = OP[oc]
            @inbounds dll = imm[3]
            @inbounds dlr = imm[4]
            @inbounds drr = imm[5]
            adlr = adj*dlr
            if op_sym == :*
                update2(tape,c_id,vidx,c_id == vidx?2.0*adlr:adlr)
            elseif op_sym == :/
                update2(tape,c_id,vidx,c_id == vidx?2.0*adlr:adlr)
                update2(tape,c_id,c_id,adj*drr)
            elseif op_sym == :^
                update2(tape,vidx,vidx,adj*dll)
                update2(tape,c_id,vidx,c_id == vidx?2.0*adlr:adlr)
                update2(tape,c_id,c_id,adj*drr)
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
                    @inbounds update2(tape,c_id,c_id,imm[1]*dlw)
                    @inbounds update2(tape,vidx,c_id,imm[1]*drw)
                    @inbounds update2(tape,vidx,vidx,imm[2]*drw)
                else
                    update2(tape,c_id,p_id,c_id==p_id?2.0*dlw:dlw)
                    update2(tape,vidx,p_id,vidx==p_id?2.0*drw:drw)
                end
            end
            #creating
            @inbounds op_sym = OP[oc]
            @inbounds dll = imm[3]
            @inbounds dlr = imm[4]
            @inbounds drr = imm[5]
            adlr = adj*dlr
            if op_sym == :*
                update2(tape,vidx,c_id,vidx == c_id?2.0*adlr:adlr)
            elseif op_sym == :/
                update2(tape,vidx,c_id,vidx == c_id?2.0*adlr:adlr)
                update2(tape,vidx,vidx,adj*drr)
            elseif op_sym == :^
                update2(tape,c_id,c_id,adj*dll)
                update2(tape,vidx,c_id,vidx == c_id?2.0*adlr:adlr)
                update2(tape,vidx,vidx,adj*drr)
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
                    update2(tape,vidx,vidx,w_bar)
                else
                    update2(tape,vidx,p_id,vidx==p_id?2.0*t1:t1)
                end
            end
            #creating
            @inbounds op_sym = OP[oc]
            if op_sym != :-
                @inbounds w_bar = adj*imm[2]
                update2(tape,vidx,vidx,w_bar)
            end
            adjlen -= 1
            idx -= 5
        end

        #@timing tape.enable_timing_stats tape.ep_reverse_times[node_id] += toq()
        # println("++++++++++++++++++++++++++++++++++++")
    end  #end while

    # @show idx, vallen, adjlen, trlen
    #@assert idx == 0
    #@assert vallen == 0 
    #@assert adjlen == 0
    #@assert trlen == 0
    nothing
end


@inline function hess_recovery(tape::Tape{Int,Float64}, start::Int, H::Vector{Float64})
    #@timing tape.enable_timing_stats tic()
    h = H
    nvar = tape.nvar
    bhi = tape.bhi
    bhv = tape.bhv

    pre = start
    for i = 1:nvar
        @inbounds lvi = bhi[i]
        @inbounds lvv = bhv[i]
        for j = 1:length(lvi)
            @inbounds v_id = lvi[j]
            if v_id <= nvar
                @inbounds h[start] = lvv[j]
                start += 1
            end
        end
    end 
    @assert start - pre == tape.nzh "$pre $start $(tape.nzh)"
    #@timing tape.enable_timing_stats tape.ep_recovery_time += toq()
    nothing
end


#Interface function
function hess_structure{I,V}(tape::Tape{I,V}, h_I::Vector{I}, h_J::Vector{I})
    #@assert tape.bh_type == 1
    hess_struct(tape)
    hess_findnz(tape)
    hess_struct_recovery(tape, h_I, h_J)
    nothing
end

function hess_reverse{I,V}(tape::Tape{I,V},vvals::Vector{V},pvals::Vector{V}, start::I, H::Vector{V})
    #@timing tape.enable_timing_stats tape.ep_n += 1

    forward_pass_2ord(tape,vvals,pvals)
    reverse_pass_2ord(tape,vvals,pvals)
    # @time recovery2(tape, tape.bhi,tape.bhv,tape.hess,tape.nvar, factor)
    hess_recovery(tape,start,H)
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
#     #@timing tape.enable_timing_stats tic()
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
#     #@assert (length(h) == nnz) length(h),nnz
#     #@timing tape.enable_timing_stats tape.ep_recovery_time += toq()
# end
