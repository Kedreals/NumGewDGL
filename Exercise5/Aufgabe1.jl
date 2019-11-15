using Plots
pyplot()
include("../integration_methods.jl")

function convergence_order(y_numeric, y_analytic, hs)
    err = log.(abs.(y_numeric .- y_analytic))
    hlog = log.(hs)

    co = (err[1:end-1] .- err[2:end]) ./ (hlog[1:end-1] .- hlog[2:end])
    return co
end

f1(x, y) = x .+ 2.0 .* y
f2(x, y) = 1.0 .+ 2.0 .* (x .+ 2.0 .* y)

y(x) = ((5.0/4.0) .* ℯ.^(2.0 .* x .- 2)) - (x ./ 2.0) .- (1.0/4.0)
h = 0.1

x0 = 1
x1 = 2
y0 = 1.0/2.0

x_n = collect(x0:h:x1)
x_a = collect(x0:((x1-x0)/1000.0):x1)

y_n = explicit_euler_all_values(f1, y0, x0, h, x1)
y_a = y(x_a)
y_ni = explicit_euler_modified(f1, f2, y0, x0, h, x1)

p = plot(x_a, y_a, label="analytic solution", yticks=0:0.5:10)
plot!(x_n, y_n, label="explicit euler h="*string(h))
plot!(x_n, y_ni, label="modified explicit euler h="*string(h))

savefig("Aufgabe1.pdf")

hs = [i^(-1) for i in 1:10000]
y_n = [explicit_euler_all_values(f1, y0, x0, hs[i], x1)[end] for i in 1:size(hs,1)]
y_ni = [explicit_euler_modified(f1, f2, y0, x0, hs[i], x1)[end] for i in 1:size(hs,1)]

con = convergence_order(y_n, y_a[end], hs)
coni= convergence_order(y_ni,y_a[end], hs)

plot(1 ./hs[1:end-1], con, label="explicit euler", title="convergence_order", xlabel="steps per unit", yticks=0:0.1:2)
plot!(1 ./hs[1:end-1],coni,label="modified explicit euler")
savefig("ConvergenceOrder.pdf")

print(con[end], "\n", coni[end], "\n")
