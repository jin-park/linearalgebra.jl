export sys_of_linear_eqs
function sys_of_linear_eqs(A, b)
    L, U = LU(A)
    # L(Ux) = b
    # L * temp = b
    # Ux = temp
    temp = zeros(size(A)[1])
    x = zeros(size(A)[1])
    temp[1] = b[1]
    for i in 2:size(A)[1]
        temp[i] = b[i] - sum(L[i, 1:i-1] .* temp[1:i-1])
    end
    x[end] = temp[end] / U[end, end]
    for i in size(A)[1]-1:-1:1
        if U[i, i] == 0
            continue
        end
        x[i] = (temp[i] - sum(U[i, i+1:end] .* x[i+1:end])) / U[i, i]
    end
    return x
end