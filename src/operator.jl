
import Lazy
using Base.Meta
using Calculus

#Operator overloading function for AD types
const B_OP_START = 1
const OP = Array{Symbol,1}([:+,:-,:*,:/,:^])
const S_TO_OC = Dict{Symbol,OP_TYPE}()
const S_TO_DIFF = Dict{Symbol,Array{Any,1}}()
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
	# @show dxx
	dxy = differentiate(dx,:y)
	# @show dxy
	dyy = differentiate(dy,:y)
	# @show dyy
	S_TO_DIFF[sym] = [dx,dy,dxx,dxy,dyy]
end

for (sym, dx) = symbolic_derivatives_1arg()
	push!(OP,sym)
	try 
		dxx = differentiate(dx,:x)
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
	if(o==:+ || o==:*)
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
		@simd for j=i:e
			@inbounds imm[imm_i]=ret/v[j]
			imm_i += 1
		end
		return ret
	end
	$switchexpr
end

##########################################################################################

switchblock = Expr(:block)
# @show switchblock
for i = U_OP_START:U_OP_END
	o = OP[i]
	(dx,dxx) = S_TO_DIFF[o]
	dx = replace_sym(:x,dx,"v[i]")
	dxx = replace_sym(:x,dxx,"v[i]")
	# @show dx, dxx
	ex = Expr(:block)
	push!(ex.args,Expr(:call,:push!,:imm,dx))
	push!(ex.args,Expr(:call,:push!,:imm,dxx))
	push!(ex.args,parse("r[1]=$(o)(v[i])"))
	push!(switchblock.args,quot(o),ex)
end

for i = B_OP_START:B_OP_END
	o = OP[i]
	ex = Expr(:block)
	if(o==:+)
		push!(ex.args,parse("ret = zero(V)"))
		#do not put second order , since all zeros
		#do not put first order , since all one
		push!(ex.args,parse("@simd for j=i:length(v) \n ret += v[j] \n end"))
		push!(ex.args,parse("r[1] = ret"))
	elseif(o==:*)
		push!(ex.args,parse("ret = one(V)"))
		push!(ex.args,parse("@simd for j=i:length(v) \n ret *= v[j] \n end"))
		push!(ex.args,parse("r[1] = ret"))
		push!(ex.args,parse("@simd for j=i:length(v) \n dxj = r[1]/v[j] \n push!(imm,dxj) \n for k=j+1:length(v) \n dxjk = dxj/v[k] \n push!(imm,dxjk) \n end \n end"))
				# for j=i:length(v)
				# 	dxj = r[1]/v[j]
				# 	push!(imm,dxj)  #first order
				# 	for k=j+1:length(v)
				# 		dxjk = dxj/v[k]
				# 		push!(imm,dxjk)  #every cross second order, diagonal is zero not pushed. 
				# 	end
				# end
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
		push!(ex.args,Expr(:call,:push!,:imm,dx))
		push!(ex.args,Expr(:call,:push!,:imm,dy))
		push!(ex.args,Expr(:call,:push!,:imm,dxx))
		push!(ex.args,Expr(:call,:push!,:imm,dxy))
		push!(ex.args,Expr(:call,:push!,:imm,dyy))
		push!(ex.args,parse("r[1]=$(o)(v[i],v[i+1])"))
	end
	push!(switchblock.args,quot(o),ex)
end

switchexpr = Expr(:macrocall, Expr(:.,:Lazy,quot(symbol("@switch"))), :s,switchblock)

#should be auto generated code
@eval function eval_2ord{I,V}(s::Symbol,v::Array{V,1},i::I,imm::Array{V,1},r::Array{V,1})
	@inbounds $switchexpr
end

##########################################################################################
#
# tape builder from a Julia Expr type
#	tape memory property is initialized 
#
##########################################################################################

#building tape with Julia expression
function tapeBuilder{I,V}(expr::Expr,tape::Tape{I}, pvals::Array{V,1})
	vset = Set{I}()
	tapeBuilder(expr,tape, pvals, vset)
	tape.nvar = length(vset)
	# @show tape
end

function tapeBuilder{I,V}(expr::Expr,tape::Tape{I}, pvals::Array{V,1}, vset::Set{I})
	tt = tape.tt
	head = expr.head
	if(head == :ref)  #a JuMP variable
		assert(length(expr.args) == 2)
		vidx = expr.args[2]
		push!(tt,TYPE_V)
		push!(tt,vidx)
		push!(tt,TYPE_V)

		push!(vset,vidx)
		tape.nvnode += 1
		tape.nnode += 1
	elseif(head == :call)
		# @show expr.args[2]
		op = expr.args[1]
		assert(typeof(op)==Symbol)
		for i in 2:length(expr.args)
			tapeBuilder(expr.args[i],tape,pvals,vset)
		end

		push!(tt,TYPE_O)
		push!(tt,S_TO_OC[op])
		push!(tt,length(expr.args)-1)
		push!(tt,TYPE_O)
		tape.nnode += 1
		tape.maxoperands < length(expr.args)-1? tape.maxoperands = length(expr.args)-1:nothing
    else
    	println("error !")
    	dump(expr)
    	assert(false)
    end
    nothing
end

function tapeBuilder{I,V}(expr::Real, tape::Tape{I}, pvals::Array{V,1},vset::Set{I}) #a JuMP parameter
	tt = tape.tt
	push!(tt,TYPE_P)
	push!(pvals,expr)
	push!(tt,length(pvals))
	push!(tt,TYPE_P)
	tape.nnode += 1
end

##########################################################################################
#
# tape builder from types
#
##########################################################################################
function tapeBuilder{I}(data::Array{I,1})
	tape = Tape{I}(data)
	return tape
end


