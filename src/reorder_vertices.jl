function reorder(t,s, strat::CommonVertex)

    # Find the permutation P of t and s that make
    # Pt = [A1, A2, P]
    # Ps = [B1, B2, P]
    I = zeros(Int,1)
    J = zeros(Int,1)
    e = 1
    for i in 1:3
        v = t[i]
        for j in 1:3
            w = s[j]
            if norm(w-v) < eps(eltype(v)) * 1.0e3
                I[e] = i
                J[e] = j
                e += 1
                break
            end
        end
        e == 2 && break
    end

    prepend!(I, setdiff([1,2,3], I))
    prepend!(J, setdiff([1,2,3], J))

    # inverse permutations
    K = indexin([1,2,3], I)
    L = indexin([1,2,3], J)

    return I, J, K, L
end


function reorder(t, s, strat::CommonEdge)

    I = zeros(Int,3)
    J = zeros(Int,3)
    e = 1
    for i in 1:3
        v = t[i]
        for j in 1:3
            w = s[j]
            if norm(w-v) < eps(eltype(v)) * 1.0e3
                I[e] = i
                J[e] = j
                e += 1
                break
            end
        end
    end
    I[3] = setdiff([1,2,3], I[1:2])[1]
    J[3] = setdiff([1,2,3], J[1:2])[1]

    # inverse permutations
    K = indexin([1,2,3], I)
    L = indexin([1,2,3], J)

    return I, J, K, L
end
