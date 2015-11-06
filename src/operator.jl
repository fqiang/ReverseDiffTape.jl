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

	eval(quote
			evaluate(::Type{Val{OP[$(oc)]}},lval::VV_TYPE, rval::VV_TYPE) =
			begin
				# ex = Expr(:call,OP[$(oc)],lval,rval)
				# debug("evaluate:",ex)
				return ($o)(lval,rval)
			end

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

	eval(quote
			evaluate(::Type{Val{OP[$(oc)]}},lval::VV_TYPE) =
			begin
				# ex = Expr(:call,OP[$(oc)],lval)
				# debug("evaluate:",ex)
				return ($o)(lval)
			end
		end)
end

evaluate(s::Symbol, lval::VV_TYPE) = evaluate(Val{s},lval)::Float64
evaluate(s::Symbol, lval::VV_TYPE, rval::VV_TYPE) = evaluate(Val{s},lval,rval)::Float64
evaluate(oc::OP_TYPE, lval::VV_TYPE) = evaluate(Val{OP[oc]},lval)::Float64
evaluate(oc::OP_TYPE, lval::VV_TYPE, rval::VV_TYPE) = evaluate(Val{OP[oc]},lval,rval)::Float64


#building tape with Julia expression
function tapeBuilder(expr::Expr,tt::TT_TYPE, pvals::TV_TYPE)
	head = expr.head
	if(head == :ref)
		assert(length(expr.args) == 2)
		vidx = expr.args[2]
		return AD_V(tt,vidx)
	elseif(head == :call)
		if(length(expr.args) >= 3) #as binary
			# @show expr.args[2]
			l = tapeBuilder(expr.args[2],tt,pvals)
			# @show l
			# @show expr.args[3]
			r = tapeBuilder(expr.args[3],tt,pvals)
			# @show r
			ret = ReverseDiffTape.(expr.args[1])(l,r)
			for i=4:1:length(expr.args)
				l = tapeBuilder(expr.args[i],tt,pvals)
				ret = ReverseDiffTape.(expr.args[1])(l,ret)
			end
			return ret
		elseif(length(expr.args) == 2) #unary 
			l = tapeBuilder(expr.args[2],tt,pvals)
            return ReverseDiffTape.(expr.args[1])(l)
		else
			dump(expr)
			assert(false)
		end
    else
    	println("error !")
    	dump(expr)
    	assert(false)
    end

end

function tapeBuilder(expr::Real, tt::TT_TYPE, pvals::TV_TYPE)
	return AD_P(tt,pvals,expr)
end

