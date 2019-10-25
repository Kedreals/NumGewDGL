using Plots
gr()

function explicit_euler(f, y0, x0, h, x)
    res = y0
    for t in x0:h:x
        res += f(res)*h
    end
    return res
end

f(y) = y^2
y(x) = -1/(x-1)
errorbound(h) = 4*(ℯ^2 - 1)*h

y0 = 1.0
x0 = 0.0
x = 1.0/2.0

k_min = 1
k_max = 12

results=zeros(k_max)
errorbounds = zeros(k_max)

for k in k_min:k_max
    h = 2.0^(-k)
    results[k] = explicit_euler(f, y0, x0, h, x)
    errorbounds[k] = errorbound(h)
end

errors = abs.(results .- y(x))

plot(k_min:k_max, [errors, errorbounds], title="Error of explicit euler with h=2^(-k)", xticks=k_min:k_max, xlabel="k", labels=["nummerical Error", "Error bound"])

savefig("errors_Aufgabe6.pdf")
