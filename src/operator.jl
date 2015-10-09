#Operator overloading function for AD types
const OP = (:+,:-,:*,:/,:^,:sin,:cos)

const B_OP_START = 1
const B_OP_END = 5
for oc = B_OP_START:1:B_OP_END
	o = OP[oc]
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
							this = AD_O(tt,$(oc),convert(UInt,lidx),convert(UInt,ridx))
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
							this = AD_O(tt,$(oc),convert(UInt,lidx),convert(UInt,ridx))
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
							this = AD_O(tt,$(oc),convert(UInt,lidx),convert(UInt,ridx))
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
							this = AD_O(tt,$(oc),convert(UInt,lidx),convert(UInt,ridx))
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

const U_OP_START = 6
const U_OP_END = 7
for oc = U_OP_START:1:U_OP_END
	o = OP[oc]
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
						this = AD_O(tt,$(oc),convert(UInt,lidx))
						return this
					end
			end)			
		elseif (LT == Type{AD_O})
			eval(quote
					($o)(::($LT), l::Placeholder) =
					begin
						tt = l.tt
						lidx = l.idx
						this = AD_O(tt,$(oc),convert(UInt,lidx))
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

evaluate(s::Symbol, lval::VV_TYPE) = evaluate(Val{s},lval)
evaluate(s::Symbol, lval::VV_TYPE, rval::VV_TYPE) = evaluate(Val{s},lval,rval)
evaluate(oc::OP_TYPE, lval::VV_TYPE) = evaluate(Val{OP[oc]},lval)
evaluate(oc::OP_TYPE, lval::VV_TYPE, rval::VV_TYPE) = evaluate(Val{OP[oc]},lval,rval)

