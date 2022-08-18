export gaussian_elimination
# switches rows a and b
function switch!(A, a, b)
    for i in axes(A)[2]
        A[a, i], A[b, i] = A[b, i], A[a, i]
    end
end

# multiplies row a by multiplier and adds it to row b
function multiply_then_add!(A, a, b, multiplier)
    for i in axes(A)[2]
        A[b, i] += A[a, i] * multiplier
    end
end

# multiplies row by multiplier
function multiply!(A, a, multiplier)
    for i in axes(A)[2]
        A[a, i] *= multiplier
    end
end

function trailing_zeros(A, a)
    count = 0
    for i in size(A)[2]:-1:1
        if A[a, i] != 0
            break
        end
        count += 1
    end
    return count
end

function leading_zeros(A, a)
    count = 0
    for i in axes(A)[2]
        if A[a, i] != 0
            break
        end
        count += 1
    end
    return count
end

function sort_by_zeros!(A)
    num_zeros = zeros(size(A)[1])
    for i in axes(A)[1]
        num_zeros[i] = leading_zeros(A, i)
    end
    zeros_sorted = sort(num_zeros)
    for i in axes(A)[1]
        if num_zeros[i] == zeros_sorted[i]
            continue
        end

        for j in i+1:size(A)[2]
            if num_zeros[j] == zeros_sorted[i]
                switch!(A, i, j)
                num_zeros[i], num_zeros[j] = num_zeros[j], num_zeros[i]
            end
        end
    end
end

function gaussian_elimination(A, b) 
    s = size(A)
    if s[1] != size(b)[1]
        error("Dimensions of A and b don't match")
    elseif s[1] < size(b)[1]
        error("Not enough equations compared to the number of variables")
    end

    augmented = Matrix{Float64}(hcat(A, b))

    sort_by_zeros!(augmented)
    if leading_zeros(A, s[1]) == s[2] && count(==(s[2]), num_zeros)
        error("There is no unique solution")
    end

    for i in axes(A)[1]
        for j in i+1:s[1]
            if leading_zeros(augmented, j) == j-1
                continue
            end

            multiply_then_add!(augmented, i, j, -augmented[j, i]/augmented[i, i])
        end
    end

    for i in s[1]:-1:1
        for j in i-1:-1:1
            if trailing_zeros(augmented, j) == s[2]-j
                continue
            end

            multiply_then_add!(augmented, i, j, -augmented[j, i]/augmented[i, i])
        end
    end

    solution = Vector{Float64}()

    for i in axes(A)[1]
        push!(solution, augmented[i, end]/augmented[i, i])
    end

    return solution
end