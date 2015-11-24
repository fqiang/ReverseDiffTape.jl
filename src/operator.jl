
import Lazy
using Base.Meta
using Calculus

#Operator overloading function for AD types
const B_OP_START = 1
const OP = Array{Symbol,1}([:+,:-,:*,:/,:^])
const S_TO_OC = Dict{Symbol,OP_TYPE}()
const S_TO_DIFF = Dict{Symbol,Array{Any,1}}()
const S_TO_DIFF_FLAG = Dict{Symbol,Array{Bool,1}}()
const B_OP_END = length(OP)
const U_OP_START = B_OP_END + 1


##########################################################################################
#
# Differetitate rule building using Calculus 
#	https://github.com/johnmyleswhite/Calculus.jl/blob/master/src/differentiate.jl
#
##########################################################################################
import Calculus.differentiate

function differentiate(x::SymbolParameter{:abs},args,wrt)
	x = args[1]
	xp = differentiate(x,wrt)
	return :($x>0?$xp:-$xp)
end

function differentiate(x::SymbolParameter{:abs2},args,wrt)
	x = args[1]
	xp = differentiate(x,wrt)
	return :($x>0?$xp:-$xp)
end

for sym in OP
	(dx,dy) = differentiate("x $(sym) y",[:x,:y])
	# @show dx, dy
	dxx = differentiate(dx,:x)
	dxxi = dxx == 0? true:false
	# @show dxx
	dxy = differentiate(dx,:y)
	dxyi = dxy == 0? true:false
	# @show dxy
	dyy = differentiate(dy,:y)
	dyyi = dyy == 0? true:false
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
		dxxi = dxx == 0? true:false
		S_TO_DIFF_FLAG[sym] = [dxxi]
		S_TO_DIFF[sym] = [dx,dxx]
	catch e
		pop!(OP)
		warn(e)
	end
	
end

const U_OP_END = length(OP)

##########################################################################################
#
# Operator overload for types to build tape
# 	tape statistic is not initialized
#
##########################################################################################
for oc = B_OP_START:B_OP_END
	o = OP[oc]
	S_TO_OC[o] = oc
	@eval 	($o){I}(l::AD{I},r::AD{I}) = 
			begin
				return AD_O($(quot(o)),l,r)
			end
end

for oc = U_OP_START:U_OP_END
	o = OP[oc]
	S_TO_OC[o] = oc
	@eval 	($o){I}(l::AD{I}) = 
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
# eval_0ord, eval_1ord, eval_2ord function gen
#
##########################################################################################
## one argument functions
switchblock = Expr(:block)
for i = U_OP_START:U_OP_END
	ex = :(return $(OP[i])(x))
    push!(switchblock.args,quot(OP[i]),ex)
end
switchexpr = Expr(:macrocall, Expr(:.,:Lazy,quot(symbol("@switch"))), :s,switchblock)
@eval function eval_0ord{V}(s::Symbol, x::V)
    $switchexpr
end

# 2+ argument functions
switchblock = Expr(:block)

for i = B_OP_START:B_OP_END
	o = OP[i]
	ex = Expr(:block)
	if(o==:+ || o==:* || o ==:^)
        continue
	else
		ex = :(@inbounds return $(o)(v[i],v[i+1]))
	end
	push!(switchblock.args,quot(o),ex)
end

switchexpr = Expr(:macrocall, Expr(:.,:Lazy,quot(symbol("@switch"))), :s,switchblock)
@eval @inline function eval_0ord{I,V}(s::Symbol, v::Array{V,1}, i::I, e::I)
    # @show s,v
    if s == :+
    	counter = zero(V)
        for j in i:e
            @inbounds counter += v[j]
        end
        return counter
    elseif s == :*
        counter = v[i]
        for j = i+1:e
            @inbounds counter *= v[j]
        end
        return counter
    elseif s == :^
        @inbounds exponent = v[i+1]
        @inbounds base = v[i]
        if exponent == 2
            return base*base
        else
            return base^exponent
        end
    end
    $switchexpr
end

##########################################################################################
# one argument functions
switchblock = Expr(:block)
for i = U_OP_START:U_OP_END
	o = OP[i]
	dx = S_TO_DIFF[o][1]
	dx = replace_sym(:x,dx,"v")
	# @show df
	ex = Expr(:block)
	push!(ex.args,parse("@inbounds imm[imm_i]=$(dx)"))
	push!(ex.args,:(@inbounds return $(o)(v)))
	push!(switchblock.args,quot(o),ex)
end
switchexpr = Expr(:macrocall,Expr(:.,:Lazy,quot(symbol("@switch"))),:s,switchblock)
@eval @inline function eval_1ord{I,V}(s::Symbol,v::V,imm::Array{V,1}, imm_i::I)
	$switchexpr
end

# 2+ argument functions
switchblock = Expr(:block)

for i = B_OP_START:B_OP_END
	o = OP[i]
	ex = Expr(:block)
	if(o==:+ || o==:* || o == :^)
		continue
	else
		(dx,dy) = S_TO_DIFF[o]
		dx = replace_sym(:y,replace_sym(:x,dx,"v[i]"),"v[i+1]")
		dy = replace_sym(:y,replace_sym(:x,dy,"v[i]"),"v[i+1]")
		# @show dx, dy
		ex = Expr(:block)
		push!(ex.args,parse("@inbounds imm[imm_i]=$(dx)"))
		push!(ex.args,parse("@inbounds imm[imm_i+1]=$(dy)"))
		push!(ex.args,:(@inbounds return $(o)(v[i],v[i+1])))
	end
	push!(switchblock.args,quot(o),ex)
end

switchexpr = Expr(:macrocall, Expr(:.,:Lazy,quot(symbol("@switch"))), :s,switchblock)

#should be auto generated code
@eval @inline function eval_1ord{I,V}(s::Symbol,v::Array{V,1}, i::I, e::I, imm::Array{V,1}, imm_i::I)
	if s==:+
		ret = zero(V)
		@simd for j=i:e
			@inbounds imm[imm_i] = one(V)
			imm_i += 1
			@inbounds ret += v[j]
		end
		return ret
	elseif s==:*
		@inbounds ret = v[i]
		@simd for j=i+1:e
			@inbounds ret *= v[j] 
		end
		# @simd for j=i:e  #unstable numerics for zero valued 
		# 	@inbounds imm[imm_i]=ret/v[j]
		# 	imm_i += 1
		# end
		@simd for j0=i:e
			c=one(V)
			@simd for j1=i:e
				if(j0!=j1)
					c*=v[j1]
				end
			end
			imm[imm_i] = c
			imm_i+=1
		end
		return ret
	elseif s==:^
		@inbounds exponent = v[i+1]
		@inbounds base = v[i]
		if exponent == 2
			imm[imm_i] = 2*base
			t = base*base
			imm[imm_i+1] = t==0?0:t*log(base) #not sure if this is way to handle log(-negative)
			return t
		else
			imm[imm_i] = exponent*base^(exponent-1)
			t = base^exponent
			imm[imm_i+1] = t*log(base)
			return t
		end
	end
	$switchexpr
end

##########################################################################################
#one argument funtion
switchblock = Expr(:block)
# @show switchblock
for i = U_OP_START:U_OP_END
	o = OP[i]
	(dx,dxx) = S_TO_DIFF[o]
	dx = replace_sym(:x,dx,"v")
	dxx = replace_sym(:x,dxx,"v")
	# @show dx, dxx
	ex = Expr(:block)
	push!(ex.args,parse("@inbounds imm[imm_i] = $(dx)"))
	push!(ex.args,parse("@inbounds imm[imm_i+1] = $(dxx)"))
	push!(ex.args,:(@inbounds return 2,$(o)(v)))
	push!(switchblock.args,quot(o),ex)
end
switchexpr = Expr(:macrocall,Expr(:.,:Lazy,quot(symbol("@switch"))),:s,switchblock)
@eval @inline function eval_2ord{I,V}(s::Symbol,v::V,imm::Array{V,1}, imm_i::I)
	$switchexpr
end

# 2+ argument function
switchblock = Expr(:block)
for i = B_OP_START:B_OP_END
	o = OP[i]
	ex = Expr(:block)
	if(o==:+ || o==:- || o==:* || o == :^)
		continue
	else
		(dx,dy,dxx,dxy,dyy) = S_TO_DIFF[o]
		dx = replace_sym(:y,replace_sym(:x,dx,"v[i]"),"v[i+1]")
		dy = replace_sym(:y,replace_sym(:x,dy,"v[i]"),"v[i+1]")
		dxx = replace_sym(:y,replace_sym(:x,dxx,"v[i]"),"v[i+1]")
		dxy = replace_sym(:y,replace_sym(:x,dxy,"v[i]"),"v[i+1]")
		dyy = replace_sym(:y,replace_sym(:x,dyy,"v[i]"),"v[i+1]")
		# @show dx
		# @show dy
		# @show dxx
		# @show dyy
		# @show dxy
		ex = Expr(:block)
		push!(ex.args,parse("@inbounds imm[imm_i] = $(dx)"))
		push!(ex.args,parse("@inbounds imm[imm_i+1] = $(dy)"))
		push!(ex.args,parse("@inbounds imm[imm_i+2] = $(dxx)"))
		push!(ex.args,parse("@inbounds imm[imm_i+3] = $(dxy)"))
		push!(ex.args,parse("@inbounds imm[imm_i+4] = $(dyy)"))
		push!(ex.args,:(@inbounds return 5,$(o)(v[i],v[i+1])))
	end
	push!(switchblock.args,quot(o),ex)
end
switchexpr = Expr(:macrocall, Expr(:.,:Lazy,quot(symbol("@switch"))), :s,switchblock)

#should be auto generated code
@eval @inline function eval_2ord{I,V}(s::Symbol,v::Array{V,1},i::I,e::I,imm::Array{V,1},imm_i::I)
	
	if(s==:+)
		@inbounds ret = v[i]
		@simd for j=i+1:e
			@inbounds ret += v[j]
		end
		return 0,ret
		#do not put second order , since all zeros
		#do not put first order , since all one
	elseif(s==:-)
		@inbounds ret = v[i]
		@inbounds ret -= v[i+1]
		return 0, ret
	elseif(s==:*)
		ret = one(V)
		@simd for j=i:e
			@inbounds ret *= v[j]
		end
		# non stable for zero values
		# n = e - i + 1
		# for j=0:n-1
		# 	dxj = ret/v[j+i]
		# 	@inbounds imm[j+imm_i] = dxj  #first order
		# 	imm_offset = round(I,(n+(n-j))*(j+1)/2)
		# 	@simd for k=j+1:n-1
		# 		@inbounds dxjk = dxj/v[k]
		# 		@inbounds imm[imm_offset+imm_i]=dxjk  #every cross second order, diagonal is zero not pushed. 
		# 		imm_offset += 1
		# 	end
		# end
		# return round(I,n+(n-1)*n/2), ret
		imm_off = zero(I)
		for j0=i:e
			f_ord = one(V)
			for j1=i:j0-1
				f_ord *= v[j1]
			end
			for j1=j0+1:e
				f_ord *= v[j1]
			end
			imm[imm_off+imm_i] = f_ord
			imm_off += 1
		end
		for j0=i:e
			for j1=j0+1:e
				s_ord = one(V)
				for j3=i:j0-1
					s_ord *= v[j3]
				end
				for j3=j0+1:j1-1
					s_ord *= v[j3]
				end
				for j3=j1+1:e
					s_ord *= v[j3]
				end
				imm[imm_off+imm_i] = s_ord
				imm_off+=1
			end
		end
		# n = e-i+1
		# @show n, imm_off, round(I,n+(n-1)*n/2)
		# assert(imm_off==round(I,n+(n-1)*n/2))
		return imm_off, ret
	elseif(s==:^)
		@inbounds exponent = v[i+1]
		@inbounds base = v[i]
		if exponent == 2
			t = base*base
			t2 = log(base)
			imm[imm_i] = 2*base
			imm[imm_i+1] = t*t2
			imm[imm_i+2] = 2.0
			imm[imm_i+3] = base + 2*base*t2
			imm[imm_i+4] = t*t2*t2
			return 5,t
		else
			t = base^exponent
			t2 = log(base)
			t3 = base^(exponent-1.0)
			imm[imm_i] = exponent*base^t3
			imm[imm_i+1] = t*t2
			imm[imm_i+2] = exponent*(exponent-1.0)*base^(exponent-2.0)
			imm[imm_i+3] = t3 + exponent*t3*t2
			imm[imm_i+4] = t*t2*t2
			return 5,t
		end
	end
	$switchexpr
end

##########################################################################################
#
# tape builder from a Julia Expr type
#	tape memory property is initialized 
#
##########################################################################################

#building tape with Julia expression
function tapeBuilder{I,V}(expr::Expr,tape::Tape{I,V}, pvals::Array{V,1})
	vset = Set{I}()
	istk = Vector{I}()
	tapeBuilder(expr,tape, pvals, vset, istk)
	tape.nvar = length(vset)
	assert(length(tape.tr) == tape.nnode-1)
	# @show tape
end

function tapeBuilder{I,V}(expr::Expr,tape::Tape{I,V}, pvals::Array{V,1},vset::Set{I},istk::Vector{I})
	tt = tape.tt
	head = expr.head
	if(head == :ref)  #a JuMP variable
		assert(length(expr.args) == 2)
		vidx = expr.args[2]
		push!(tt,TYPE_V)
		push!(tt,vidx)
		push!(tt,TYPE_V)
		
		push!(istk,length(tt)-2)
		push!(vset,vidx)
		tape.nvnode += 1
		tape.nnode += 1
		tape.eset[length(tt)-2] = Dict{I,V}()
		tape.liveVar[length(tt)-2] = Set{I}()
	elseif(head == :call)
		# @show expr.args[2]
		op = expr.args[1]
		assert(typeof(op)==Symbol)
		for i in 2:length(expr.args)
			tapeBuilder(expr.args[i],tape,pvals,vset,istk)
		end
		n = length(expr.args)-1
		push!(tt,TYPE_O)
		push!(tt,S_TO_OC[op])
		push!(tt,n)
		push!(tt,TYPE_O)

		t = Vector{I}()
		for i=1:n
			push!(t,pop!(istk))
		end
		append!(tape.tr,reverse!(t))
		push!(istk,length(tt)-3)
		tape.nnode += 1
		tape.maxoperands < length(expr.args)-1? tape.maxoperands = length(expr.args)-1:nothing
		tape.imm2ord += n + round(I,n*(n+1)/2)
		tape.eset[length(tt)-3] = Dict{I,V}()
		tape.liveVar[length(tt)-3] = Set{I}()
		# @show length(tt) - 3
    else
    	println("error !")
    	dump(expr)
    	assert(false)
    end
    nothing
end

function tapeBuilder{I,V}(expr::Real, tape::Tape{I,V}, pvals::Array{V,1},vset::Set{I},istk::Vector{I}) #a JuMP parameter
	tt = tape.tt
	push!(tt,TYPE_P)
	push!(tt,length(pvals)+1)
	push!(tt,TYPE_P)

	push!(istk,length(tt)-2)
	push!(pvals,expr)
	tape.nnode += 1
	tape.eset[length(tt)-2] = Dict{I,V}()
	tape.liveVar[length(tt)-2] = Set{I}()
end

##########################################################################################
#
# tape builder from types
#
##########################################################################################
function tapeBuilder{I}(data::Array{I,1})
	tape = Tape{I,Float64}(data)
	return tape
end


