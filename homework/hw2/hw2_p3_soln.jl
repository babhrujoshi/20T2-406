using Random
using LinearAlgebra
using Convex
using ECOS
using Plots
using JLD
using CSV
using DataFrames

cd("/Users/babhru/Documents/GitHub/20T2-406/docs/homework/hw2")

C = load("hw2_p3_C.jld")["data"]
d = load("hw2_p3_d.jld")["data"]
# load signal
x0 = load("hw2_p3_signal.jld")["data"]

s = 1
α = 0.5
β=0.25
epsilon = 10^(-4)

#Gradient descent with constant step size
x = []
obj_gd = []
push!(x, [1000;-500])

for iter in 1:100
    x_temp = x[iter]
    r = sum((C.-x_temp).^2,dims=1)' - d.^2;
    push!(obj_gd, norm(r,2)^2/2)
    g = 4*sum((ones(2,1)*r').*(x_temp.-C),dims = 2)
    push!(x, x_temp - g*1e-10);
end
plot(obj_gd,yaxis=:log, label = "Grad. desc. constant")

#Gradient descent with line search
x = []
obj_gd_line = []
push!(x, [1000;-500])

for iter in 1:100
    x_temp = x[iter]
    r = sum((C.-x_temp).^2,dims=1)' - d.^2;
    push!(obj_gd_line, norm(r,2)^2/2)
    g = 4*sum((ones(2,1)*r').*(x_temp.-C),dims = 2)

    global t = s
    global obj_cand = obj_gd_line[iter]
    counter = 1
    while (obj_gd_line[iter] .- obj_cand .- α*t*(g'*g))[1] < 1e-10 && counter < 1000 
        global t = t*β
        xcand = x_temp .- t*g
        rcand = sum((C.-xcand).^2,dims=1)' .- d.^2;
        global obj_cand = norm(rcand,2)/2
        counter +=1
    end   
    push!(x, x_temp - t*g);   
end
plot!(obj_gd_line,yaxis=:log, label = "Grad. desc. line search")

#Gauss newton with constant step size
x = []
obj_gd_gauss = []
push!(x, [1000;-500])

for iter in 1:100
    x_temp = x[iter]
    r = sum((C.-x_temp).^2,dims=1)' - d.^2
    push!(obj_gd_gauss, norm(r,2)^2/2)
    J = (x_temp.-C)'
    z = J \ r
    push!(x, x_temp-z)
end
plot!(obj_gd_gauss,yaxis=:log, label = "Gauss-Newton constant")

#Gauss newton with line search
x = []
obj_gd_gauss_line = []
push!(x, [1000;-500])

for iter in 1:100
    x_temp = x[iter]
    r = sum((C.-x_temp).^2,dims=1)' - d.^2
    push!(obj_gd_gauss_line, norm(r,2)^2/2)
    J = (x_temp.-C)'
    z = J \ r
    global t = s
    global obj_cand = norm(r,2)^2/2
    g = 4*sum((ones(2,1)*r').*(x_temp.-C),dims = 2)
    counter = 1
    while (obj_gd_gauss_line[iter] .- obj_cand .- α*t*(g'*z))[1] < 1e-10 && counter < 1000
        global t = t*β
        xcand = x_temp .- t*z
        rcand = sum((C.-xcand).^2,dims=1)' .- d.^2;
        global obj_cand = norm(rcand,2)/2
        counter += 1
    end   

    push!(x, x_temp-t*z)
end
plot!(obj_gd_gauss_line,yaxis=:log,label= "Gauss-Newton line search")

