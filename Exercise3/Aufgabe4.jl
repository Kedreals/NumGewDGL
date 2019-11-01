using Plots
pyplot()
include("../integration_methods.jl")

f(x, y) = y^(3/2)
y0 = 1
x0 = 0
x1 = 1
h = 0.2

cp(ξ,η) = -ξ+2/sqrt(η)
cn(ξ,η) = -ξ-2/sqrt(η)

Yp(x,ξ,η)=4 ./ ((x .+ cp(ξ,η)) .^2)
Yn(x,ξ,η)=4 ./ ((x .+ cn(ξ,η)) .^2)

yh = explicit_euler_vec(f, y0, x0, h, x1)
x = collect(x0:h:x1)

print("yh = ", yh, "\n")
print("cp = ", [cp(x[i], yh[i]) for i in 1:size(yh,1)-1], "\n")
print("cn = ", [cn(x[i], yh[i]) for i in 1:size(yh,1)-1], "\n")

p1 = plot(x, yh, title = "negative c", label="explicit euler", xticks=x0:h:x1, linewidth=2.5)

for i in 1:size(yh,1)-1
    plot!(x, Yn(x, x[i], yh[i]), label="i="*string(i-1))
end

p2 = plot(x, yh, title = "positive c", label= "explicit euler", xticks=x0:h:x1, linewidth=2.5)

for i in 1:size(yh,1)-1
    plot!(x, Yp(x, x[i], yh[i]), label="i="*string(i-1))
end

plot(p2, p1, size=(1720, 890))

savefig("Aufgabe4.pdf")
