# LU Decomposition for square matrices


export LU
function LU(A::AbstractMatrix) # iterative method
    # only square matrices
    # L is going to be a unit lower triangle matrix i.e. the diagonals will all be 1
    s = size(A)[1]
    L = Matrix{Float64}(undef, s, s)
    U = Matrix{Float64}(undef, s, s)
    
    # initialize all the 0 and 1 values
    for i in 1:s
        L[i, i] = 1
        L[1:i-1, i] .= 0
        U[i+1:s, i] .= 0
    end
    
    # order of procedure: figure out row 1 of U
    #                     then repeat i=2 to s:
    #                                   figure out position i, n of L
    #                                   figure out row i, n+1 of U
    U[1, :] = A[1, :]
    for i in 2:s
        # row i of L
        for j in 1:s
            # figure out L
            # from L, [i, 1:j-1] .* from U, [1:i, j] + [i, j+1:s] .*  [j+1:s, i]
            # L[i, j] = (A[i, j] - sum from above) / U[j, i]
            if j < i
                left = (L[i, 1:j-1] .* U[1:i-1, j]) |> sum
                right = (L[i, j+1:end] .* U[i+1:end, j]) |> sum
                L[i, j] = (A[i, j] - (left + right)) / U[j, j]
            else
                # figure out U
                up = (L[i, 1:j-1] .* U[1:i-1, j]) |> sum
                down = (L[i, j+1:end] .* U[i+1:end, j]) |> sum
                U[i, j] = (A[i, j] - (up + down)) / L[j, j]
            end
        end 
    end
    print("L: ")
    println(L)
    print("U: ")
    println(U)
end























    # a = A[1, 1]
    # w_transpose = A[1, 2:end]
    # v = A[2:end, 1]
    # A′ = A[2:end, 2:end] # A \prime<TAB>
