facts("Hessian EP algorithm x1^2*x2^2") do
    p = Vector{Float64}()
    x = Vector{Float64}()
    x1 = AD_V(x,1.1)
    x2 = AD_V(x,2.2)
    p2 = AD_P(p,2.0)
    c = x1^p2*x2^p2
    tt = tapeBuilder(c.data)

    hess_structure2(tt)
    hess_reverse2(tt,x,p)
    
    h = sparse(tt.h_I,tt.h_J,tt.hess)

    @fact length(h.nzval) --> 3
    @fact h[1 , 1] --> 2.2*2.2*2
    @fact h[2 , 1] --> 2*1.1*2*2.2
    @fact h[2 , 2] --> 1.1*1.1*2 
end

facts("Hessian EP algorithm sin(x1)") do
    p = Vector{Float64}()
    x = Vector{Float64}()
    x1 = AD_V(x,1.1)
    c = sin(x1)
    tt = tapeBuilder(c.data)
    hess_structure2(tt)
    hess_reverse2(tt,x,p)

    h = sparse(tt.h_I,tt.h_J,tt.hess)

    @fact length(h.nzval) --> 1
    @fact h[1 , 1] --> -sin(1.1) 
end


facts("Hessian EP algorithm cos(x1)") do
    p = Vector{Float64}()
    x = Vector{Float64}()
    x1 = AD_V(x,1.1)
    c = cos(x1)
    tt = tapeBuilder(c.data)
    hess_structure2(tt)
    hess_reverse2(tt,x,p)
    h = sparse(tt.h_I,tt.h_J,tt.hess)

    @fact length(h.nzval) --> 1
    @fact h[1 , 1] --> -cos(1.1) 
end

facts("Hessian EP algorithm x1*x2") do
    p = Vector{Float64}()
    x = Vector{Float64}()
    x1 = AD_V(x,1.1)
    x2 = AD_V(x,2.2)
    c = x1*x2
    tt = tapeBuilder(c.data)
    hess_structure2(tt)
    hess_reverse2(tt,x,p)
    h = sparse(tt.h_I,tt.h_J,tt.hess)
    @fact length(h.nzval) --> 1
    @fact h[2 , 1] --> 1.0
end


facts("Hessian EP algorithm sin(x1*x2)") do
    p = Vector{Float64}()
    x = Vector{Float64}()
    x1 = AD_V(x,1.1)
    x2 = AD_V(x,2.2)
    c = sin(x1*x2)
    tt = tapeBuilder(c.data)
    tt = tapeBuilder(c.data)
    hess_structure2(tt)
    hess_reverse2(tt,x,p)
    h = sparse(tt.h_I,tt.h_J,tt.hess)
    @fact length(h.nzval) --> 3
    @fact h[2 , 1] --> cos(1.1*2.2) - 1.1*2.2*sin(1.1*2.2)
    @fact h[1 , 1] --> -2.2^2*sin(1.1*2.2)
    @fact h[2 , 2] --> -1.1^2*sin(1.1*2.2)
end

facts("Hessian EP algorithm cos(sin(x1))") do
    p = Vector{Float64}()
    x = Vector{Float64}()
    x1 = AD_V(x,1.1)
    c = cos(sin(x1))
    tt = tapeBuilder(c.data)
    hess_structure2(tt)
    hess_reverse2(tt,x,p)
    h = sparse(tt.h_I,tt.h_J,tt.hess)
    @fact length(h.nzval) --> 1
    @fact h[1 , 1] --> sin(sin(1.1))*sin(1.1) - cos(sin(1.1))*(cos(1.1)^2)
end


facts("Hessian EP algorithm cos(sin(x1*x2))") do
    p = Vector{Float64}()
    x = Vector{Float64}()
    x1 = AD_V(x,1.1)
    x2 = AD_V(x,2.2)
    c = cos(sin(x1*x2))
    tt = tapeBuilder(c.data)
    hess_structure2(tt)
    hess_reverse2(tt,x,p)
    h = sparse(tt.h_I,tt.h_J,tt.hess)
    @fact length(h.nzval) --> 3
    @fact h[1 , 1] --> roughly(2.2^2*sin(sin(1.1*2.2))*sin(1.1*2.2)-2.2^2*cos(sin(1.1*2.2))*cos(1.1*2.2)^2)
end

facts("Hessian EP algorithm x1*x1") do
    p = Vector{Float64}()
    x = Vector{Float64}()
    x1 = AD_V(x,1.1)
    c = x1*x1
    tt = tapeBuilder(c.data)
    hess_structure2(tt)
    hess_reverse2(tt,x,p)
    h = sparse(tt.h_I,tt.h_J,tt.hess)
    @fact length(h.nzval) --> 1
    @fact h[1 , 1] --> 2.0
end

facts("Hessian EP algorithm x1*x1*x2*x2") do
    p = Vector{Float64}()
    x = Vector{Float64}()
    x1 = AD_V(x,1.1)
    x2 = AD_V(x,2.2)
    p2 = AD_P(p,2)
    c = x1*x1*x2*x2
    tt = tapeBuilder(c.data)
    hess_structure2(tt)
    hess_reverse2(tt,x,p)
    h = sparse(tt.h_I,tt.h_J,tt.hess)
    @fact length(h.nzval) --> 3
    @fact h[1 , 1] --> 2.2*2.2*2
    @fact h[2 , 1] --> 2*1.1*2*2.2
    @fact h[2 , 2] --> 1.1*1.1*2 
end

facts("Hessian EP algorithm x1*x1*x1") do
    p = Vector{Float64}()
    x = Vector{Float64}()
    x1 = AD_V(x,1.1)
    p2 = AD_P(p,2)
    c = x1*x1*x1
    tt = tapeBuilder(c.data)
    hess_structure2(tt)
    hess_reverse2(tt,x,p)
    h = sparse(tt.h_I,tt.h_J,tt.hess)
    @fact length(h.nzval) --> 1
    @fact h[1 , 1] --> 6*1.1
end

facts("Hessian EP algorithm cos(x1*x2)") do
    p = Vector{Float64}()
    x = Vector{Float64}()
    x1 = AD_V(x,1.1)
    x2 = AD_V(x,2.2)
    c = cos(x1*x2)
    tt = tapeBuilder(c.data)
    hess_structure2(tt)
    hess_reverse2(tt,x,p)
    h = sparse(tt.h_I,tt.h_J,tt.hess)
    @fact length(h.nzval) --> 3
    @fact h[1 , 1] --> -cos(1.1*2.2)*2.2*2.2
    @fact h[2 , 1] --> -(sin(1.1*2.2)+cos(1.1*2.2)*1.1*2.2)
    @fact h[2 , 2] -->  -cos(1.1*2.2)*1.1*1.1 
end

facts("Hessian EP algorithm y=sin(x1) + cos(x2^p2) ") do
    p = Vector{Float64}()
    x = Vector{Float64}()
    x1 = AD_V(x,1.1)
    x2 = AD_V(x,2.2)
    p2 = AD_P(p,2.0)
    c = sin(x1) + cos(x2^p2)
    tt = tapeBuilder(c.data)
    hess_structure2(tt)
    hess_reverse2(tt,x,p)
    h = sparse(tt.h_I,tt.h_J,tt.hess)
    @fact length(h.nzval) --> 2
    @fact h[1,1] --> -sin(1.1)
    @fact h[2,2] --> -( 2*2*2.2*2.2*cos(2.2*2.2) + 2*sin(2.2*2.2) )
end

facts("Hessian EP algorithm y=sin(x1) + cos(x2^p2) ") do
    p = Vector{Float64}()
    x = Vector{Float64}()
    x1 = AD_V(x,1.1)
    x2 = AD_V(x,2.2)
    p2 = AD_P(p,2.0)
    c = x1*(x1+x2)
    tt = tapeBuilder(c.data)
    hess_structure2(tt)
    hess_reverse2(tt,x,p)
    h = sparse(tt.h_I,tt.h_J,tt.hess)
    @fact length(h.nzval) --> 2
    @fact h[1,1] --> 2.0
    @fact h[2,1] --> 1.0
end

facts("Hessian EP algorithm y=x1/x2 + x1/(x1+x2)") do
    p = Vector{Float64}()
    x = Vector{Float64}()
    x1 = AD_V(x,1.1)
    x2 = AD_V(x,2.2)
    p2 = AD_P(p,2.0)
    c = x1/x2+x1/(x1+x2)
    tt = tapeBuilder(c.data)
    hess_structure2(tt)
    hess_reverse2(tt,x,p)
    h = sparse(tt.h_I,tt.h_J,tt.hess)

    #symbolic
    y1 = 1.1
    y2 = 2.2
    (dy1,dy2) = differentiate("y1/y2+y1/(y1+y2)",[:y1,:y2])
    hess = Vector{Float64}()
    push!(hess, @eval $(differentiate(dy1,[:y1])[1]))
    push!(hess, @eval $(differentiate(dy1,[:y2])[1]))
    push!(hess, @eval $(differentiate(dy2,[:y2])[1]))

    @fact length(h.nzval) --> 3
    @fact h[1,1] --> roughly(hess[1])
    @fact h[2,1] --> roughly(hess[2])
    @fact h[2,2] --> roughly(hess[3])
end


facts("Hessian EP algorithm sin(x1)+cos(x2^2)*1-x3*2") do
    p = Vector{Float64}()
    x = Vector{Float64}()
    x1 = AD_V(x, 1.1)
    x2 = AD_V(x, 2.2)
    x3 = AD_V(x, 3.3)
    p1 = AD_P(p, 1.0)
    p2 = AD_P(p, 2.0)
    c = sin(x1)+cos(x2^p2) * p1 - x3*p2
    tt = tapeBuilder(c.data)
    hess_structure2(tt)
    hess_reverse2(tt,x,p)
    h = sparse(tt.h_I,tt.h_J,tt.hess)

    @fact length(h.nzval) --> 2
    @fact h[1 , 1] --> -sin(1.1)
    @fact h[2 , 2] --> -2*sin(2.2^2)-4*2.2^2*cos(2.2^2)
end