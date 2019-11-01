function explicit_euler(f, y0, x0, h, x)
    res = y0
    for t in (x0+h):h:x
        res += f(t, res)*h
    end
    return res
end

function explicit_euler_vec(f, y0, x0, h, x_end)
    x = collect(x0:h:x_end)
    y = zeros(size(x,1))
    y[1] = y0

    for i in 1:size(y,1)-1
        y[i+1] = y[i] + h*f(x[i], y[i])
    end

    return y
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

function implicit_euler_vec(f, df, y0, x0, h, x_end)
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
