using Plots: gr, plot, savefig
gr()
# exact solution
y(x) = (4.0/15.0)*x.^(5/2) - (19/15)*x .+ 1
# second derivative
ddy(x)= sqrt.(x)

#sample point distance
h = 0.01

# get all sample points
xi = collect(0:h:1)
# get the coerrect results
yi = y(xi)

# set up the right side of the linear equation system
b = ddy(xi)
b[1] = y(0)
b[end] = y(1)

# set up the system matrix
Mat = zeros(size(b,1), size(b,1))
Mat[1,1] = h^2
Mat[end,end] = h^2
for i = 2:size(b,1)-1
    Mat[i, i-1] = 1
    Mat[i,i] = -2
    Mat[i, i+1] = 1
end
Mat = Mat ./ h^2

# calculate the approximated function values
nyi = Mat\b

# calculate the error
error = abs.(nyi-yi)

plot(xi, error, title="Error over x", xlabel="x")
savefig("error_Aufgabe4.pdf")
