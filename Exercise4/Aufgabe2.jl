using LinearAlgebra
include("../integration_methods.jl")

f(t, x)=[-x[1]+x[2], x[1]-x[2]]
fv(t,x)=[-x[:,1].+x[:,2] x[:,1].-x[:,2]]
x(t) = 1.0/2.0 .* [1+ℯ^(-2*t), 1-ℯ^(-2*t)]

t0 = 0
x0 = [1,0]
h = 1.0/10.0
t = 1.0

x_h=explicit_euler_vec(f, x0, t0, h, t)

for i in 1:size(x_h,1)-1
    print("t = ", (i-1)*h, ":   x_h(t) = ", x_h[i,:], ",  x(t) = ", x((i-1)*h), "\n")
end

t_i = collect(0:h:1)

error = 0

for (i,t) in enumerate(t_i)
    n = norm(f(t, x_h[i,:]))
    int = (ℯ^(2*(1-t-t_i[i+1])) * (ℯ^(2*t_i[i+1])+ℯ^(2*t)*(-1+2*(t-t_i[i+1])))) / 4.0

    global error += n*int
    if i+1 == size(t_i,1)
        break
    end
end

error_besser = 0

for (i,t) in enumerate(t_i)
    n = norm(f(t,x_h[i,:]))

    global error_besser += n*h^2/2
    if i+1 == size(t_i,1)
        break
    end
end

print("\ne_h(1) =", error, "\n")
print("bessere Fehler Abschätzung = ", error_besser, "\nreal error = ", norm(x_h[end,:] - x(1)), "\n")
