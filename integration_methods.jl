function explicit_euler(f, y0, x0, h, x)
    res = y0
    for t in x0:h:(x-h)
        res += f(t, res)*h
    end
    return res
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

function implicit_euler_all_values(f, df, y0, x0, h, x_end)
    x_val = collect(x0:h:x_end)
    y_val = zeros(size(x_val, 1))
    y_val[1] = y0

    for i in 2:size(x_val, 1)
        y_new = y_val[i-1]
        # use Newton-Ralphs to find root of y_n+1 - y_n - h*f(x_n+1, y_n+1)
        # that root is the next y value
        for j in 1:5
            y_new = y_new - (y_new - y_val[i-1] - h*f(x_val[i], y_new))/(1-h*df(x_val[i], y_new))
        end
        y_val[i] = y_new
    end
    return y_val
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
