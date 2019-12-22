include("../integration_methods.jl")

f(t,x)=[-0.04*x[1]+10.0^4.0*x[2]*x[3];
        0.04*x[1]-10.0^4.0*x[2]*x[3]-3.0*10.0^7.0*x[2]*x[2];
        3.0*10.0^7.0*x[2]*x[2]]

df(t, x)=[(-0.04)                    (10.0^4.0*x[3])  (10.0^x[2]);
           (0.04) (-10.0^4.0*x[3]-6.0*10.0^7.0*x[2]) (-10.0^x[2]);
              (0)                (6.0*10.0^7.0*x[2])  (0)]

x_s = [1;0;0]
x_e = [9.886739393819*10.0^-1.0;3.447715743689*10.0^-5.0;1.129158346063*10.0^-2.0]
t_s = 0

h=10.0^-3.0
t_e = 0.3

#first call of @time is a dump
print("first time @time call")
@time ERK_vec(f, x_s, t_s, h, t_e, c=[0], A=[0], b=[1])
@time implicit_euler_vec(f, df, x_s, t_s, h, t_e)
print("End of test time measurements!-----------------------\n\n")

print("Explicit euler with h=$h:")
x_explicit = @time ERK_vec(f, x_s, t_s, h, t_e, c=[0], A=[0], b=[1])
h *= 10.0
print("Implicit euler with h=$h:")
print("1 Newton step:")
x_implicit_1 = @time implicit_euler_vec(f, df, x_s, t_s, h, t_e, newton_terminate=1)
print("10 Newton steps:")
x_implicit_10 = @time implicit_euler_vec(f, df, x_s, t_s, h, t_e, newton_terminate=10)
print("20 Newton steps:")
x_implicit_20 = @time implicit_euler_vec(f, df, x_s, t_s, h, t_e, newton_terminate=20)

print("At t=$t_e\nExplicit Error=$(x_explicit[end,:]-x_e)\nImplicit Error 01 steps: $(x_implicit_1[end,:]-x_e)\nImplicit Error 10 steps: $(x_implicit_10[end,:]-x_e)\nImplicit Error 20 steps: $(x_implicit_20[end,:]-x_e)\n")
