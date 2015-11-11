
import Lazy
using Base.Meta


#Operator overloading function for AD types
const B_OP_START = 1
const OP = Array{Symbol,1}([:+,:-,:*,:/,:^])
const S_TO_OC = Dict{Symbol,OP_TYPE}()
const S_TO_DER_1ARG = Dict{Symbol,Expr}()
const B_OP_END = length(OP)
const U_OP_START = B_OP_END + 1
for (sym, exp) = symbolic_derivatives_1arg()
	push!(OP,sym)
	S_TO_DER_1ARG[sym] = exp
end
const U_OP_END = length(OP)


##########################################################################################
#
# Operator overload for types
#
##########################################################################################
for oc = B_OP_START:1:B_OP_END
	o = OP[oc]
	S_TO_OC[o] = oc
	@eval 	($o){I,V}(l::Placeholder{I,V},r::Placeholder{I,V}) = 
			begin
				assert(l.tape == r.tape)
				tape = l.tape
				return AD_O{I,V}(tape,$(quot(o)),l,r)
			end
end


for oc = U_OP_START:1:U_OP_END
	o = OP[oc]
	S_TO_OC[o] = oc

	@eval 	($o)(l::Placeholder) = 
			begin
				tape = l.tape
				return AD_O(tape,$(quot(o)),l)
			end
end
##########################################################################################

switchblock = Expr(:block)

for i = U_OP_START:1:U_OP_END
	ex = parse("@inbounds v[1]=$(OP[i])(vals[i])")
    push!(switchblock.args,quot(OP[i]),ex)
end

for i = B_OP_START:1:B_OP_END
	o = OP[i]
	if(o==:+)
		ex = Expr(:block)
		push!(ex.args,parse("@inbounds v[1] = zero(V)"))
		push!(ex.args,parse("@simd for j=i:1:length(vals) \n @inbounds v[1]+=vals[j] \n end"))
	elseif(o==:*)
		ex = Expr(:block)
		push!(ex.args,parse("@inbounds v[1] = vals[i]"))
		push!(ex.args,parse("@simd for j=i+1:1:length(vals) \n @inbounds v[1] *= vals[j] \n end"))
	else
		ex = parse("@inbounds  v[1]=$(o)(vals[i],vals[i+1])")
	end
	push!(switchblock.args,quot(o),ex)
end

switchexpr = Expr(:macrocall, Expr(:.,:Lazy,quot(symbol("@switch"))), :s,switchblock)

@eval function evaluate{V,I}(s::Symbol, vals::Array{V,1}, i::I, v::Array{V,1})
    $switchexpr
end

swithblock = Expr(:block)
for i = U_OP_START:1:U_OP_END
	ex = parse("")
end

##########################################################################################
#
# tape builder from a Julia Expr type
#
##########################################################################################

#building tape with Julia expression
function tapeBuilder(expr::Expr,tape::Tape{IDX_TYPE}, pvals::TV_TYPE)
	vset = Set{IDX_TYPE}()
	tapeBuilder(expr,tape, pvals, vset)
	tape.nvar = length(vset)
	# @show tape
end

function tapeBuilder(expr::Expr,tape::Tape{IDX_TYPE}, pvals::TV_TYPE, vset::Set{IDX_TYPE})
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

function tapeBuilder(expr::Real, tape::Tape{IDX_TYPE}, pvals::TV_TYPE,vset::Set{IDX_TYPE}) #a JuMP parameter
	tt = tape.tt
	push!(tt,TYPE_P)
	push!(pvals,expr)
	push!(tt,length(pvals))
	push!(tt,TYPE_P)
	tape.nnode += 1
end
##########################################################################################

