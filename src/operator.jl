
import Lazy
using Base.Meta

#Operator overloading function for AD types
const B_OP_START = 1
const OP = Array{Symbol,1}([:+,:-,:*,:/,:^])
const S_TO_OC = Dict{Symbol,OP_TYPE}()
const S_TO_DIFF = Dict{Symbol,Array{Any,1}}()
const B_OP_END = length(OP)
const U_OP_START = B_OP_END + 1


for sym in OP
	(dx,dy) = differentiate("x $(sym) y",[:x,:y])
	@show dx, dy
	S_TO_DIFF[sym] = [dx,dy]
end

for (sym, exp) = symbolic_derivatives_1arg()
	push!(OP,sym)
	S_TO_DIFF[sym] = [exp]
end

const U_OP_END = length(OP)


function replace_sym(s::Symbol,ex::AbstractString,r::AbstractString)
	return replace_sym(s,parse(ex),parse(r))
end
function replace_sym(s::Symbol,ex::Expr,r::AbstractString)
	return replace_sym(s,ex,parse(r))
end
function replace_sym(s::Symbol,ex::Number,r)
	if(ex == 1) return parse("one(V)")
	elseif(ex == 0) return parse("zero(V)") 
	else return ex
	end
end
function replace_sym(s::Symbol,ex::Symbol,r::Expr)
	if(ex==s) return r
	else return ex
	end
end
function replace_sym(s::Symbol,ex::Expr, r::Expr) 
	for i = 1:1:length(ex.args)
		ex.args[i] = replace_sym(s,ex.args[i],r)
	end
	return ex
end

##########################################################################################
#
# Operator overload for types to build tape
# 	tape statistic is not initialized
#
##########################################################################################
for oc = B_OP_START:1:B_OP_END
	o = OP[oc]
	S_TO_OC[o] = oc
	@eval 	($o){I}(l::AD{I},r::AD{I}) = 
			begin
				return AD_O($(quot(o)),l,r)
			end
end

for oc = U_OP_START:1:U_OP_END
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

switchblock = Expr(:block)
for i = U_OP_START:1:U_OP_END
	ex = parse("@inbounds v[1]=$(OP[i])(vals[i])")
    push!(switchblock.args,quot(OP[i]),ex)
end

for i = B_OP_START:1:B_OP_END
	o = OP[i]
	ex = Expr(:block)
	if(o==:+)
		push!(ex.args,parse("@inbounds v[1] = zero(V)"))
		push!(ex.args,parse("@simd for j=i:1:length(vals) \n @inbounds v[1]+=vals[j] \n end"))
	elseif(o==:*)
		push!(ex.args,parse("@inbounds v[1] = vals[i]"))
		push!(ex.args,parse("@simd for j=i+1:1:length(vals) \n @inbounds v[1] *= vals[j] \n end"))
	else
		ex = parse("@inbounds  v[1]=$(o)(vals[i],vals[i+1])")
	end
	push!(switchblock.args,quot(o),ex)
end

switchexpr = Expr(:macrocall, Expr(:.,:Lazy,quot(symbol("@switch"))), :s,switchblock)

@eval function evaluate{V,I}(s::Symbol, vals::Array{V,1}, i::I, v::Array{V,1})
    @inbounds $switchexpr
end

##########################################################################################

switchblock = Expr(:block)
# @show switchblock
for i = U_OP_START:1:U_OP_END
	o = OP[i]
	df = S_TO_DIFF[o][1]
	df = replace_sym(:x,df,"v[i]")
	@show df
	ex = Expr(:block)
	push!(ex.args,Expr(:call,:push!,:imm,df))
	push!(ex.args,parse("r[1]=$(o)(v[i])"))
	push!(switchblock.args,quot(o),ex)
end

for i = B_OP_START:1:B_OP_END
	o = OP[i]
	ex = Expr(:block)
	if(o==:+)
		push!(ex.args,parse("ret = zero(V)"))
		push!(ex.args,parse("@simd for j=i:1:length(v) \n push!(imm,one(V)) \n ret += v[j] \n end"))
		push!(ex.args,parse("r[1] = ret"))
	elseif(o==:*)
		push!(ex.args,parse("ret = one(V)"))
		push!(ex.args,parse("@simd for j=i:1:length(v) \n ret *= v[j] \n end"))
		push!(ex.args,parse("r[1] = ret"))
		push!(ex.args,parse("@simd for j=i:1:length(v) \n push!(imm,r[1]/v[j]) \n end"))
	else
		(dx,dy) = S_TO_DIFF[o]
		dx = replace_sym(:x,dx,"v[i]")
		dx = replace_sym(:y,dx,"v[i+1]")
		dy = replace_sym(:x,dy,"v[i]")
		dy = replace_sym(:y,dy,"v[i+1]")
		@show dx, dy
		ex = Expr(:block)
		push!(ex.args,Expr(:call,:push!,:imm,dx))
		push!(ex.args,Expr(:call,:push!,:imm,dy))
		push!(ex.args,parse("r[1]=$(o)(v[i],v[i+1])"))
	end
	push!(switchblock.args,quot(o),ex)
end

switchexpr = Expr(:macrocall, Expr(:.,:Lazy,quot(symbol("@switch"))), :s,switchblock)

#should be auto generated code
@eval function eval_idd{V,I}(s::Symbol,v::Array{V,1},i::I,imm::Array{V,1},r::Array{V,1})
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
		for i in 2:1:length(expr.args)
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


