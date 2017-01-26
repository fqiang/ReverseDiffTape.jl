function write_log_dat_ampl(lambda,M,N)
	theta = Vector{Float64}(N)	
	rand!(theta)
	theta = theta.*7

	f = open(string("log",M,"_",N,".dat"),"w")
	
	write(f, string("# true theta \n"))
	for i=1:N
		write(f,@sprintf("# \t %0.2f \n", theta[i]))
	end

	write(f, "reset data; \n")
	write(f, "data; \n")
	write(f, string("param m := ",M,"; \n"))
	write(f, string("param n := ",N,"; \n"))
	write(f, string("param lambda := ", lambda ,"; \n"))
	write(f, "param \t x: ")
	for i=1:N
		write(f, string("\t ",i))
	end
	write(f, "\t := \n")

	x = randn(M,N)
	x[:,N] = ones(M)
	for i=1:M
		write(f, string("\t",i))
		xx = x[i,:]
		for j=1:N
			write(f, string("\t",xx[j]))
		end
		write(f, "\n");
	end
	write(f, "\t ; \n")
	write(f, "param \t y \t := \n")
	
	p = 1./(1+exp(-x*theta))
	y = Vector{Float64}(M)
	for i=1:M
		write(f, string("\t",i))
		y[i] = generate_bernoulli(p[i])
		write(f, string("\t", y[i], "\n"))
	end
	write(f, "\t ; \n")

	write(f, "option solver ipopt; \n");
	write(f, "solve; \n");

	close(f)
	return theta, x, y
end

function generate_bernoulli(p)
	a = rand()
	if(a<=p)
		return 1.0
	else
		return -1.0
	end
end

function write_log_dat_julia(lambda,t,x,y)
	M, N = size(x)
	f = open(string("log",M,"_",N,"_dat.jl"),"w")
	write(f, string("lambda = ",lambda,"\n"))
	write(f, string("M=",M,"\n"))
	write(f, string("N=",N,"\n"))
	write(f,string("x=Array{Float64,2}(M,N) \n"))
	for i = 1:M
		write(f, string("x[",i,",:]=",x[i,:], "\n"))
	end
	write(f,string("y=",y))
	write(f,"\n")

	close(f)
end

function genLogData(M,N)
	lambda = 1.0
	t, x, y = write_log_dat_ampl(lambda, M, N )
	write_log_dat_julia(lambda, t , x , y)
end

# ---------------------------
function genRandData(N,K)
	A = write_random_dat_ampl(N,K)
	writedlm(open(string("randset",N,"_",K,".txt"),"w"),A)
end

function rand_set(n,k)
    ret = []
    for i=1:k
        push!(ret,rand(1:n))
    end
    return ret
end


function write_random_dat_ampl(N, K)
	f = open(string("randset",N,"_",K,".dat"),"w")

	write(f, string("data; \n"))
	write(f, string("param n:=",N,"; \n"))
	write(f, string("param k:=",K,"; \n"))
	write(f, string("param rand: "));
	for i =1:K
		write(f, string("\t ",i))
	end
	write(f, string(":= \n"))
	
	A = Array{Int, 2}(N,K);
	for i =1:N
		A[i,:] = rand_set(N,K)
	end

	for i=1:N
		AA = A[i,:]
		write(f, string("\t",i,"\t"))
		for aa in AA
			write(f, string("\t",aa))
		end
		write(f,"\n")
	end
	write(f,string("; \n"))
	write(f,string("option solver ipopt; \n"))
	write(f,string("solve; \n"))

	write(f, "\n")
	close(f)
	return A
# data; 
# param n:=4000;
# param k:=2;

# #display N;
# #display K;

# for {i in 1..n}
# {
#     let {b in 1..k} rand[i, b] := round(Uniform(1,n));
# }

# #display rand;

# option solver ipopt;
# solve;


end

