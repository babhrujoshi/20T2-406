using Random
using LinearAlgebra
using Convex
using ECOS
using Plots
using JLD
using CSV
using DataFrames

cd("/Users/babhru/Documents/GitHub/20T2-406/docs/homework/hw2")

# load data
A = load("hw2_p2_smooth_A.jld")["data"]
b = load("hw2_p2_smooth_b.jld")["data"]
x0 = load("hw2_p2_smooth_signal.jld")["data"]

(m,n) = size(A)
x1 = Variable(n);


loss = sumsquares(A*x1-b)
reg = sum(abs(x1[2:end]-x1[1:end-1]))

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


# plot(x0)
plot!(xa)
plot!(xb)


plot_temp = plot(1:n, [x0,xb,xa], label = ["original" "lambda = 0.01" "lambda = 10"], xlabel="index, i", ylabel="x(i)",
title = "Total variation, solution comparision")


# generate logspace between -3,3
Lambda = exp10.(range(-3, stop=3, length=100))
f1_values = []
f2_values = []
for i in 1:100
    local λ = Lambda[i]
    local problem = minimize(loss+λ*reg)
    solve!(problem,ECOS.Optimizer(verbose = false))
    xd = x1.value;
    free!(x1)
    push!(f1_values, 0.5*norm(A*xd-b,2)^2)
    push!(f2_values, sum(abs.(xd[2:end]-xd[1:end-1]) ))

end
plot_temp = plot(f1_values, f2_values, label = "Pareto optimal curve", xlabel="f1", ylabel="f2",title = "Total Variation, pareto curve")