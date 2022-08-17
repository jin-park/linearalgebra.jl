# solves the equation Ax = b

# 1. form the augmented matrix [A b]
# 2. perform row operations -> row echelon form
# 3. determine if the equation has a trivial or non trivial solution
# 4. 

# operations are too expensive
# function row_switching_matrix(n, a, b) # rows a and b, number of rows n
#     m = zeros(n, n)
#     for i in 1:n
#         m[i, i] = 1
#     end
#     m[a, a] = 0
#     m[a, b] = 1
#     m[b, b] = 0
#     m[b, a] = 1
#     return m
# end
export switch!, gaussian_elimination, multiply_then_add!, multiply!, sort_by_zeros!

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
    size = size(A)
    if size[1] != size(b)[1]
        error("Dimensions of A and b don't match")
    end
    augmented = hcat(A, b)

end