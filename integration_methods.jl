using LinearAlgebra

function explicit_euler(f, y0, x0, h, x)
    res = y0
    for t in x0:h:(x-h)
        res += f(t, res)*h
    end
    return res
end

function ooc(hs, errors)
    return (log.(abs.(errors[1:end-1]))-log.(abs.(errors[2:end])))./(log.(hs[1:end-1])-log.(hs[2:end]))
end

function explicit_euler_all_values(f, y0, x0, h, x_end)
    x = collect(x0:h:x_end)
    y = zeros(size(x,1))
    y[1] = y0

    for i in 1:size(y,1)-1
        y[i+1] = y[i] + f(x[i], y[i])*h
    end

    return y
end

function special_linear_2_schritt_verfahren(f, df, x0s, t0, h, t_end, α; newton_iterationen=5)
    t = collect(t0:h:t_end)
    x = zeros(size(t, 1))
    x[1:2] = x0s
    αs = [1.0, -1.0-α, α]
    βs = [α/12.0+5.0/12.0, -2.0/3.0*α+2.0/3.0, -5.0*α/12.0-1.0/12.0]

    for i in 1:size(t,1)-2
        x[i+2] = x[i+1]
        g(T,X)=αs[1]*X + αs[2]*x[i+1] + αs[3]*x[i] -
            h*(βs[1]*f(T, X) + βs[2]*f(t[i+1], x[i+1]) + βs[3]*f(t[i], x[i]))

        dg(T,X)= 1.0-βs[1]*h*df(T,X)
        for j in 1:newton_iterationen
            s = g(t[i+2], x[i+2])/dg(t[i+2], x[i+2])
            x[i+2] = x[i+2] - s
        end
    end
    return x
end

function explicit_euler_modified(f1, f2, y0, x0, h, x_end)
    x = collect(x0:h:x_end)
    y = zeros(size(x,1))
    y[1] = y0

    for i in 1:size(y,1)-1
        y[i+1] = y[i] + h*f1(x[i], y[i]) + (h^2.0)/2.0 * f2(x[i], y[i])
    end
    return y
end

function explicit_euler_vec(f, x0, t0, h, t_end)
    t = collect(t0:h:t_end)
    x_h = zeros((size(t,1), size(x0,1)))
    x_h[1,:] = x0

    for i in 1:size(t,1)-1
        x_h[i+1,:] = x_h[i,:] .+ h.*f(t[i], x_h[i,:])
    end
    return x_h
end

function ERK_vec(f, x0, t0, h, t_end; c=[0; 1/2; 1/2; 1], A=[0 0 0 0; 1/2 0 0 0; 0 1/2 0 0; 0 0 1 0], b=[1/6 1/3 1/3 1/6])
    t = collect(t0:h:t_end)
    x_h = zeros((size(t,1), size(x0,1)))
    x_h[1,:] = x0

    for i in 1:size(t,1)-1
        v = zeros((size(c,1), size(x0,1)))
        v[1,:] = f(t[i], x_h[i,:])
        for j in 2:size(v,1)
            v[j,:] = f(t[i]+c[j]*h, x_h[i,:] + h.*sum(A[j,1:(j-1)].*v[1:(j-1),:], dims =1)')
        end
        x_h[i+1,:] = x_h[i,:] + h .* sum(b[1:size(v,1)] .* v[1:size(v,1),:], dims=1)'
    end
    return x_h
end

function _step(f, x, t, h, c, A, b)
    k = zeros((size(c,1), size(x0,1)))
    k[1,:] = f(t, x)
    for j in 2:size(k,1)
        k[j,:] = f(t+c[j]*h, x + h.*sum(A[j,1:(j-1)].*k[1:(j-1),:], dims=1)')
    end
    x_b = x + h.* sum(b[2,1:size(k,1)] .* k[1:size(k,1),:], dims=1)'
    x_t = x + h.* sum(b[1,1:size(k,1)] .* k[1:size(k,1),:], dims=1)'
    return x_b, x_t
end

function _errApprox(x_b, x_t, x_h)
    err = 0
    for j in 1:size(x_b, 1)
        err = max(err, abs(x_b[j]-x_t[j])/(1.0+abs(x_h[j])))
    end
    return err
end

function ARK_vec(f, x0, t0, h_start, t_end; c=[0; 1], A = [0 0; 1 0], b=[1 0;1/2 1/2], tol=10.0^-4.0, q=[1,2])
    h = h_start; t = t0; x_b = x0; x_t = x0;
    x_h = []; t_h = [];
    err = 0;
    adaptive_h(H, E) = min(2.0, max(1.0/2.0, 9.0/10.0*(tol/E)^(1.0/(min(q...)+1.0))))*H
    push!(x_h, x_b)
    push!(t_h, t)
    while t<= t_end && t_h[end] != t_end
        #h = h_start
        if t+h > t_end
            h = t_end-t
        end
        x_b, x_t = _step(f, x_h[end], t, h, c, A, b)

        err = _errApprox(x_b, x_t, x_h[end])

        while err > tol
            h = adaptive_h(h, err)
            if t+h > t_end
                h = t_end-t
                x_b, x_t = _step(f, x_h[end], t, h, c, A, b)
                break;
            end
            x_b, x_t = _step(f, x_h[end], t, h, c, A, b)
            err = _errApprox(x_b, x_t, x_h[end])
        end
        t+=h;
        h = adaptive_h(h, err)
        push!(x_h,x_b)
        push!(t_h, t)
    end
    return x_h, t_h
end

function RK3_vec(f, x0, t0, h, t_end)
    t = collect(t0:h:t_end)
    x_h = zeros((size(t,1), size(x0,1)))
    x_h[1,:] = x0

    for i in 1:size(t,1)-1
        v1 = f(t[i], x_h[i,:])
        v2 = f(t[i] + 1.0/2.0 * h, x_h[i,:] + h/2.0 * v1)
        v3 = f(t[i+1], x_h[i,:] + h*v2)
        v4 = f(t[i+1], x_h[i,:] + h*v3)

        x_h[i+1,:] = x_h[i,:] + h*(1.0/6.0 * v1 + 2.0/3.0 * v2 + 1.0/6.0 * v4)
    end
    return x_h
end

function implicit_euler(f, df, y0, x0, h, x)
    res = y0

    for t in x0+h:h:x
        y_new = res
        for j in 1:5
            y_new = y_new - (y_new - res - h*f(x, y_new))/(1-h*df(x,y_new))
        end
        res = y_new
    end

    return res
end

function implicit_euler_vec(f, J, x0, t0, h, te; newton_terminate=5)
    t_val = collect(t0:h:te)
    x_val = zeros(size(t_val, 1), size(x0, 1))
    x_val[1,:] = x0
    identity = I(size(x0, 1))

    for i in 2:size(t_val, 1)
        x_val[i,:] = x_val[i-1,:]
        #use Newton Ralphs to find root of x_n+1 - x_n-h*f(t_n+1, x_n+1)
        Dg(t,X)= identity-h*J(t, X)
        g(t,X) = X-x_val[i-1,:]-h*f(t,X)
        for j in 1:newton_terminate

            step = Dg(t_val[i],x_val[i,:])\-g(t_val[i], x_val[i,:]);
            x_val[i,:] = x_val[i,:] + step
        end
    end
    return x_val
end

function implicit_euler_all_values(f, df, y0, x0, h, x_end; newton_terminate=5)
    x_val = collect(x0:h:x_end)
    y_val = zeros(size(x_val, 1))
    y_val[1] = y0

    for i in 2:size(x_val, 1)
        y_new = y_val[i-1]
        # use Newton-Ralphs to find root of y_n+1 - y_n - h*f(x_n+1, y_n+1)
        # that root is the next y value
        for j in 1:newton_terminate
            y_new = y_new - (y_new - y_val[i-1] - h*f(x_val[i], y_new))/(1-h*df(x_val[i], y_new))
        end
        y_val[i] = y_new
    end
    return y_val
end

function milne_simpson_two(f, df, x_start, t0, h, t_end)
    t = collect(t0:h:t_end)
    res = zeros((size(t,1), size(x_start, 1)))
    res[1:2,:] = x_start[1:2]

    for n in 2:size(t,1)-1
        res[n+1,:] = res[n,:]
        s = (1.0/3.0)*h*(4.0*f(t[n], res[n,:]) + f(t[n-1], res[n-1,:]))
        #use Newton-Ralphs to find root of x_{n+1}-x_{n-1}-1/3*h*sum_{k=0}^2\beta_{3,k}f(t_{j+1-k},x_{j+1-k})
        for j in 1:5
            res[n+1,:]= res[n+1,:] - (res[n+1,:] - res[n-1,:]-(1.0/3.0)*h*f(t[n+1],res[n+1,:])-s)/(1.0 .-h*df(t[n+1],res[n+1,:]))
        end
    end
    return res
end

function Nyström_three(f, x_start, t0, h, t_end)
    t = collect(t0:h:t_end)
    res = zeros((size(t,1), size(x_start, 1)))
    res[1:3, :] = x_start[1:3]

    for n in 3:size(t,1)-1
        res[n+1,:] = res[n-1,:] + h/3.0 * (7*f(t[n], res[n,:]) - 2*f(t[n-1], res[n-1,:]) + f(t[n-2],res[n-2, :]))
    end
    return res
end

function adams_moulton_three(f, df, x_start, t0, h, t_end)
    t = collect(t0:h:t_end)
    res = zeros((size(t,1), size(x_start, 1)))
    res[1:3,:] = x_start[1:3]

    for n in 3:size(t,1)-1
        res[n+1,:] = res[n,:]
        s = h*((8.0/12.0)*f(t[n],res[n,:])-(1.0/12.0)*f(t[n-1],res[n-1,:]))
        #use Newton-Ralphs to find root of x_n+1-x_n-h*sum_{k=0}^2\beta_{3,k}f(t_{j+1-k},x_{j+1-k})
        for j in 1:5
            res[n+1,:]= res[n+1,:] .- (res[n+1,:] .- res[n,:] .-h*(5.0/12.0)*f(t[n+1],res[n+1,:]).-s)/(1.0 .-h*df(t[n+1],res[n+1,:]))
        end
    end
    return res
end

function adams_bashforth_three(f, x_start, t0, h, t_end)
    t = collect(t0:h:t_end);
    res = zeros((size(t,1), size(x_start,1)))
    res[1:3,:] = x_start[1:3];

    for n in 3:size(t,1)-1
        res[n+1,:] = res[n,:] + h*(23.0/12.0*f(t[n], res[n,:])-16.0/12.0*f(t[n-1],res[n-1,:])+5.0/12.0*f(t[n-2],res[n-2,:]))
    end


    return res
end

function adams_bashforth_five(f, x0, t0, h, t_end)
    t = collect(t0:h:t_end);
    res = zeros((size(t,1), size(x0,1)))
    res[1,:] = x0;

    if size(t,1)< 2
        return res;
    end
    res[2,:] = res[1,:] + h*f(t[1], res[1,:])

    if size(t,1) < 3
        return res;
    end
    res[3,:] = res[2,:] + h*(3.0/2.0*f(t[2],res[2,:])-1.0/2.0*f(t[1],res[1,:]));

    if size(t,1) < 4
        return res;
    end
    res[4,:] = res[3,:] + h*(23.0/12.0*f(t[3], res[3,:])-16.0/12.0*f(t[2],res[2,:])+5.0/12.0*f(t[1],res[1,:]));

    if size(t,1) < 5
        return res;
    end
    res[5,:] = res[4,:] + h*(55.0/24.0*f(t[4], res[4,:])-59.0/24.0*f(t[3],res[3,:])+37.0/24.0*f(t[2],res[2,:])-9.0/24.0*f(t[1],res[1,:]));

    for n in 1:size(t,1)-5
        res[n+5,:] = res[n+4,:] +
            h*(1901.0/720.0*f(t[n+4],res[n+4,:])-
                2774.0/720.0*f(t[n+3],res[n+3,:])+
                2616.0/720.0*f(t[n+2],res[n+2,:])-
                1274.0/720.0*f(t[n+1],res[n+1,:])+
                251.0/720.0*f(t[n],res[n,:]));
    end

    return res
end
