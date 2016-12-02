
import Lazy
using Base.Meta
using Calculus

#Operator overloading function for AD types
const B_OP_START = 1
const OP = Array{Symbol,1}([:+,:-,:*,:/,:^])
const B_OP_END = length(OP)
const U_OP_START = B_OP_END + 1

const S_TO_OC = Dict{Symbol,Int}()
const S_TO_DIFF = Dict{Symbol,Array{Any,1}}()
const S_TO_DIFF_FLAG = Dict{Symbol,Array{Bool,1}}()



##########################################################################################
#
# Differetitate rule building using Calculus 
#   https://github.com/johnmyleswhite/Calculus.jl/blob/master/src/differentiate.jl
#
##########################################################################################
import Calculus.differentiate

function differentiate(x::SymbolParameter{:abs},args,wrt)
    x = args[1]
    xp = differentiate(x,wrt)
    return :($x>0?$xp:-$xp)
end


for sym in OP
    (dx,dy) = differentiate("x $(sym) y",[:x,:y])
    # @show dx, dy
    dxx = differentiate(dx,:x)
    dxxi = isa(dxx,Number)
    # @show dxx
    dxy = differentiate(dx,:y)
    dxyi = isa(dxy,Number)
    # @show dxy
    dyy = differentiate(dy,:y)
    dyyi = isa(dyy,Number)
    # @show dyy
    # @show sym,dxx,dxy,dyy
    # @show sym,dxxi,dxyi,dyyi
    S_TO_DIFF[sym] = [dx,dy,dxx,dxy,dyy]
    S_TO_DIFF_FLAG[sym] = [dxxi,dxyi,dyyi]
end

for (sym, dx) = symbolic_derivatives_1arg()
    push!(OP,sym)
    try 
        dxx = differentiate(dx,:x)
        dxxi = isa(dxx,Number)
        S_TO_DIFF[sym] = [dx,dxx]
        S_TO_DIFF_FLAG[sym] = [dxxi]
    catch e
        pop!(OP)
        warn(e)
    end
end

const U_OP_END = length(OP)

##########################################################################################
#
# Operator overload for types to build tape
#   tape statistic is not initialized
#
##########################################################################################
for oc = B_OP_START:B_OP_END
    o = OP[oc]
    S_TO_OC[o] = oc
    @eval   ($o){I}(l::AD{I},r::AD{I}) = 
            begin
                return AD_O($(quot(o)),l,r)
            end
end

for oc = U_OP_START:U_OP_END
    o = OP[oc]
    S_TO_OC[o] = oc
    @eval   ($o){I}(l::AD{I}) = 
            begin
                return AD_O($(quot(o)),l)
            end
end

(+){I}(l::AD{I}, args...) = AD_O(:+,tuple(l,args...))
(*){I}(l::AD{I}, args...) = AD_O(:*,tuple(l,args...))

##########################################################################################

function replace_sym(s::Symbol,ex::AbstractString,r::AbstractString)
    return replace_sym(s,parse(ex),parse(r))
end
function replace_sym(s::Symbol,ex::Expr,r::AbstractString)
    return replace_sym(s,ex,parse(r))
end

function replace_sym(s::Symbol,ex::Number,r)
    return ex
end
function replace_sym(s::Symbol,ex::Symbol,r)
    if(ex==s) return r
    else return ex
    end
end
function replace_sym(s::Symbol,ex::Expr, r) 
    for i = 1:length(ex.args)
        ex.args[i] = replace_sym(s,ex.args[i],r)
    end
    return ex
end

##########################################################################################
#
# 0ord eval
#
##########################################################################################
## one argument functions
switchblock = Expr(:block)
for i = U_OP_START:U_OP_END
    o = OP[i]
    ex = :(return $(o)(v))
    push!(switchblock.args,quot(o),ex)
end
switchexpr = Expr(:macrocall, Expr(:.,:Lazy,quot(symbol("@switch"))), :s,switchblock)
@eval @inline function eval_0ord{V}(s::Symbol, v::V)
    if s==:-
        return -v
    else
        $switchexpr
    end
end

# 2 argument functions
switchblock = Expr(:block)

for i = B_OP_START:B_OP_END
    o = OP[i]
    ex = Expr(:block)
    push!(ex.args,parse("return $(o)(l,r)"))
    # ex = :(@inbounds return $(o)(l,r))
    push!(switchblock.args,quot(o),ex)
end

switchexpr = Expr(:macrocall, Expr(:.,:Lazy,quot(symbol("@switch"))), :s,switchblock)
@eval @inline function eval_0ord{V}(s::Symbol,l::V,r::V)
    $switchexpr
end

# 2+ argument functions
@eval @inline function eval_0ord{I,V}(s::Symbol, v::Vector{V}, i::I, e::I)
    # @show s,v
    ret = v[i]
    if s==:+
        @simd for j = i+1:e
            @inbounds ret +=  v[j]
        end
    else
        @simd for j = i+1:e
            @inbounds ret *=  v[j]
        end
    end
    return ret
end


##########################################################################################
#
# 1ord eval
#
##########################################################################################
# one argument functions
switchblock = Expr(:block)
for i = U_OP_START:U_OP_END
    o = OP[i]
    dx = S_TO_DIFF[o][1]
    dx = replace_sym(:x,dx,"v")
    # @show dx
    ex = Expr(:block)
    push!(ex.args,parse("@inbounds imm[1]=$(dx)"))
    push!(switchblock.args,quot(o),ex)
end
switchexpr = Expr(:macrocall,Expr(:.,:Lazy,quot(symbol("@switch"))),:s,switchblock)
@eval @inline function eval_1ord{V}(s::Symbol,v::V,imm::Vector{V})
    if s==:-
        @inbounds imm[1] = -one(V)
    else
        $switchexpr
    end
    nothing
end

# 2 argument functions
switchblock = Expr(:block)

for i = B_OP_START:B_OP_END
    o = OP[i]
    ex = Expr(:block)
    
    (dx,dy) = S_TO_DIFF[o]
    dx = replace_sym(:y,replace_sym(:x,dx,"l"),"r")
    dy = replace_sym(:y,replace_sym(:x,dy,"l"),"r")
    # @show dx, dy
    ex = Expr(:block)
    push!(ex.args,parse("@inbounds imm[1]=$(dx)"))
    push!(ex.args,parse("@inbounds imm[2]=$(dy)"))
    push!(switchblock.args,quot(o),ex)
end

switchexpr = Expr(:macrocall, Expr(:.,:Lazy,quot(symbol("@switch"))), :s,switchblock)
# @show switchexpr
@eval @inline function eval_1ord{V}(s::Symbol,l::V, r::V, imm::Vector{V})
    if s==:^
        exponent::V = r
        base::V = l
        t0 = l<0.0?0.0:log(l)
        t = l^r
        if exponent == 2.0
            @inbounds imm[1] = 2.0*l
            @inbounds imm[2] = t*t0 
        else
            @inbounds imm[1] = r*l^(r-1.0)
            @inbounds imm[2] = t*t0
        end
    else
        $switchexpr
    end
    nothing
end

@eval @inline function eval_1ord{I,V}(s::Symbol, v::Vector{V},i::I,e::I,imm::Vector{V})
    if s==:+
        # @inbounds imm[1:(e-i)+1] = one(V)
    else
        imm_i = zero(I)
        @simd for j0=i:e
            c=one(V)
            @simd for j1=i:e
                if j0!=j1
                    @inbounds c *= v[j1]
                end
            end
            imm_i += 1
            @inbounds imm[imm_i] = c
        end
    end
    nothing
end

##########################################################################################
#
# 2ord eval
#
##########################################################################################
#one argument funtion
switchblock = Expr(:block)
for i = U_OP_START:U_OP_END
    o = OP[i]
    (dx,dxx) = S_TO_DIFF[o]
    dx = replace_sym(:x,dx,"v")
    dxx = replace_sym(:x,dxx,"v")
    # @show dx, dxx
    ex = Expr(:block)
    push!(ex.args,parse("@inbounds imm[1] = $(dx)"))
    push!(ex.args,parse("@inbounds imm[2] = $(dxx)"))
    push!(switchblock.args,quot(o),ex)
end
switchexpr = Expr(:macrocall,Expr(:.,:Lazy,quot(symbol("@switch"))),:s,switchblock)
@eval @inline function eval_2ord{V}(s::Symbol,v::V,imm::Vector{V})
    if s==:-
        @inbounds imm[1] = -one(V)
        @inbounds imm[2] = zero(V)
    else
        $switchexpr
    end
    nothing
end

# 2 argument function
switchblock = Expr(:block)
for i = B_OP_START:B_OP_END
    o = OP[i]
    ex = Expr(:block)
    (dx,dy,dxx,dxy,dyy) = S_TO_DIFF[o]
    dx = replace_sym(:y,replace_sym(:x,dx,"l"),"r")
    dy = replace_sym(:y,replace_sym(:x,dy,"l"),"r")
    dxx = replace_sym(:y,replace_sym(:x,dxx,"l"),"r")
    dxy = replace_sym(:y,replace_sym(:x,dxy,"l"),"r")
    dyy = replace_sym(:y,replace_sym(:x,dyy,"l"),"r")
    # @show dx
    # @show dy
    # @show dxx
    # @show dyy
    # @show dxy
    ex = Expr(:block)
    push!(ex.args,parse("@inbounds imm[1] = $(dx)"))
    push!(ex.args,parse("@inbounds imm[2] = $(dy)"))
    push!(ex.args,parse("@inbounds imm[3] = $(dxx)"))
    push!(ex.args,parse("@inbounds imm[4] = $(dxy)"))
    push!(ex.args,parse("@inbounds imm[5] = $(dyy)"))
    push!(switchblock.args,quot(o),ex)
end
switchexpr = Expr(:macrocall, Expr(:.,:Lazy,quot(symbol("@switch"))), :s,switchblock)
@eval @inline function eval_2ord{V}(s::Symbol,l::V, r::V, imm::Vector{V})
    if s == :^
        exponent = r
        base = l
        if exponent == 2.0
            t = base*base
            t2 = base<0.0?0.0:log(base)
            @inbounds imm[1] = 2.0*base
            @inbounds imm[2] = t*t2
            @inbounds imm[3] = 2.0
            @inbounds imm[4] = base + 2.0*base*t2
            @inbounds imm[5] = t*t2*t2
        else
            t = base^exponent
            t2 = base<0.0?0.0:log(base)
            t3 = base^(exponent-1.0)
            @inbounds imm[1] = exponent*base^t3
            @inbounds imm[2] = t*t2
            @inbounds imm[3] = exponent*(exponent-1.0)*base^(exponent-2.0)
            @inbounds imm[4] = t3 + exponent*t3*t2
            @inbounds imm[5] = t*t2*t2
        end
    else
        $switchexpr
    end
    nothing
end

@eval @inline function eval_2ord{I,V}(s::Symbol,v::Vector{V},i::I,e::I,imm::Vector{V})
    if s==:+
        # @inbounds imm[1:(e-i)+1] = one(V)
    else
        #first order    
        imm_i = zero(I)
        for j0=i:e
            f_ord = one(V)
            for j1=i:j0-1
                @inbounds f_ord *= v[j1]
            end
            for j1=j0+1:e
                @inbounds f_ord *= v[j1]
            end
            imm_i += 1
            @inbounds imm[imm_i] = f_ord
        end
        #second order
        for j0=i:e
            for j1=j0+1:e
                s_ord = one(V)
                for j3=i:j0-1
                    @inbounds s_ord *= v[j3]
                end
                for j3=j0+1:j1-1
                    @inbounds s_ord *= v[j3]
                end
                for j3=j1+1:e
                    @inbounds s_ord *= v[j3]
                end
                imm_i+=1
                @inbounds imm[imm_i] = s_ord
            end
        end
    end
    nothing
end
##########################################################################################

