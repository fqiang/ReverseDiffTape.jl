#Operator overloading function for AD types
const B_OP_START = 1
const OP = Array{Symbol,1}([:+,:-,:*,:/,:^])
const OC_TO_OP = Dict{Symbol,OP_TYPE}()
const B_OP_END = length(OP)
const U_OP_START = B_OP_END + 1
for (syms, exp) = symbolic_derivatives_1arg()
	push!(OP,syms)
end
const U_OP_END = length(OP)

for oc = B_OP_START:1:B_OP_END
	o = OP[oc]
	OC_TO_OP[o] = oc

	# println("setup operator ",o)
	for LT in AD_TYPES
		for RT in AD_TYPES
			# println("setup ",LT," ",RT)
			if ((LT == Type{AD_V} || LT == Type{AD_P}) && (RT == Type{AD_V} || RT == Type{AD_P}))
				eval(quote
						($o)(::($LT),::($RT), l::Placeholder, r::Placeholder) =
						begin
							assert(l.tt == r.tt)
							tt = l.tt
							push!(tt,l.idx)
							push!(tt,l.t)
							lidx = length(tt)
							push!(tt,r.idx)
							push!(tt,r.t)
							ridx = length(tt)
							this = AD_O(tt,$(oc),lidx,ridx)
							return this
						end
					end)
			elseif ((LT == Type{AD_V} || LT == Type{AD_P}) && RT == Type{AD_O})
				eval(quote
						($o)(::($LT),::($RT), l::Placeholder, r::Placeholder) =
						begin
							assert(l.tt == r.tt)
							tt = l.tt
							push!(tt,l.idx)
							push!(tt,l.t)
							lidx = length(tt)
							ridx = r.idx
							this = AD_O(tt,$(oc),lidx,ridx)
							return this
						end
				end)
			elseif ((LT == Type{AD_O}) && (RT == Type{AD_V} || RT == Type{AD_P}))
				eval(quote
						($o)(::($LT),::($RT), l::Placeholder, r::Placeholder) =
						begin
							assert(l.tt == r.tt)
							tt = l.tt
							lidx = l.idx
							push!(tt,r.idx)
							push!(tt,r.t)
							ridx = length(tt)
							this = AD_O(tt,$(oc),lidx,ridx)
							return this
						end
				end)
			elseif ((LT == Type{AD_O}) && (RT == Type{AD_O}))
				eval(quote
						($o)(::($LT),::($RT), l::Placeholder, r::Placeholder) =
						begin
							assert(l.tt == r.tt)
							tt = l.tt
							lidx = l.idx
							ridx = r.idx
							this = AD_O(tt,$(oc),lidx,ridx)
							return this
						end
				end)			
			else
				println(LT,"  ",RT)
				assert(false)
			end
		end
	end

	eval(quote
			($o)(l::Placeholder,r::Placeholder) = ($o)(typeof(l),typeof(r),l,r)
		end)
end


for oc = U_OP_START:1:U_OP_END
	o = OP[oc]
	OC_TO_OP[o] = oc
	# println("setup operator ",o)
	for LT in AD_TYPES
		# println("setup ",LT)
		if (LT == Type{AD_P} || LT == Type{AD_V})
			eval(quote
					($o)(::($LT), l::Placeholder) =
					begin
						tt = l.tt
						push!(tt,l.idx)
						push!(tt,l.t)
						lidx = length(tt)
						this = AD_O(tt,$(oc),lidx)
						return this
					end
			end)			
		elseif (LT == Type{AD_O})
			eval(quote
					($o)(::($LT), l::Placeholder) =
					begin
						tt = l.tt
						lidx = l.idx
						this = AD_O(tt,$(oc),lidx)
						return this
					end
			end)	
		else
			println(LT,"  ",RT)
			assert(false)
		end
	end

	eval(quote
			($o)(l::Placeholder) = ($o)(typeof(l),l)
		end)


end

##########################################################################################
import Lazy
using Base.Meta

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
##########################################################################################

#building tape with Julia expression
function tapeBuilder(expr::Expr,tape::Tape{IDX_TYPE}, pvals::TV_TYPE)
	vset = Set{IDX_TYPE}()
	tapeBuilder(expr,tape, pvals, vset)
	tape.nvar = length(vset)
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
		# if(length(expr.args) >= 3) #binary should has 3 arguments, and we should also support multiple argument such as sum/prod list
			# @show expr.args[2]
			op = expr.args[1]
			assert(typeof(op)==Symbol)
			for i in 2:1:length(expr.args)
				tapeBuilder(expr.args[i],tape,pvals,vset)
			end

			push!(tt,TYPE_O)
			push!(tt,OC_TO_OP[op])
			push!(tt,length(expr.args)-1)
			push!(tt,TYPE_O)
			tape.nnode += 1
			tape.maxoperands < length(expr.args)-1? tape.maxoperands = length(expr.args)-1:nothing
			######################### -- old code for only support binary op
			# l = tapeBuilder(expr.args[2],tt,pvals)
			# # @show l
			# # @show expr.args[3]
			# r = tapeBuilder(expr.args[3],tt,pvals)
			# # @show r
			# ret = ReverseDiffTape.(expr.args[1])(l,r)
			# for i=4:1:length(expr.args)
			# 	l = tapeBuilder(expr.args[i],tt,pvals)
			# 	ret = ReverseDiffTape.(expr.args[1])(l,ret)
			# end
			# return ret
			####################### 
		# elseif(length(expr.args) == 2) #unary has only 2 arguments
		# 	op = expr.args[1]
		# 	assert(typeof(op)==Symbol)
		# 	push!(tt,TYPE_O)
		# 	push!(tt,OC_TO_OP[op])
		# 	push!(tt,1)
		# 	push!(tt,TYPE_O)
		# else
		# 	dump(expr)
		# 	assert(false)
		# end
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

