using Random
using LinearAlgebra
using Convex
using ECOS
using Plots
using JLD
using CSV
using DataFrames

# parts b, c, d
cd("/Users/babhru/Documents/GitHub/20T2-406/docs/homework/hw2")
# load data
A = load("hw2_p2_sparse_A.jld")["data"]
b = load("hw2_p2_sparse_b.jld")["data"]
x0 = load("hw2_p2_sparse_signal.jld")["data"]

(m,n) = size(A)
x1 = Variable(n);

loss = sumsquares(A*x1-b)
reg = norm(x1,1)

λ = 10
problem = minimize(loss+λ*reg)
solve!(problem,ECOS.Optimizer(verbose = false))
free!(x1)
xa = x1.value;

λ = .001
problem = minimize(loss+λ*reg)
solve!(problem,ECOS.Optimizer(verbose = false))

xb = x1.value;

# generate logspace between -3,3

Lambda = exp10.(range(-3, stop=3, length=50))


plot(x0)
plot!(xa)
plot!(xb)


