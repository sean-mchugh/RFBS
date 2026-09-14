using ExponentialUtilities

struct RFCladoEvent
    left::Int
    right::Int
    prob::Float64
end

struct RFSparseCladoEvents
    by_parent::Vector{Vector{RFCladoEvent}}
    nnz::Int
    nstates::Int
end

struct RFBranchTransitionCache
    by_child::Vector{Union{Nothing, Matrix{Float64}}}
end

function make_sparse_clado_events(cladoPmat::Array{Float64, 3})
    nparents, nleft, nright = size(cladoPmat)
    nleft == nparents || error("cladoPmat left-child dimension does not match parent dimension")
    nright == nparents || error("cladoPmat right-child dimension does not match parent dimension")

    by_parent = [RFCladoEvent[] for _ in 1:nparents]
    nnz = 0

    for parent in 1:nparents
        parent_events = by_parent[parent]
        for left in 1:nleft
            for right in 1:nright
                prob = cladoPmat[parent, left, right]
                if prob != 0.0
                    push!(parent_events, RFCladoEvent(left, right, prob))
                    nnz += 1
                end
            end
        end
    end

    return RFSparseCladoEvents(by_parent, nnz, nparents)
end

count_sparse_clado_events(events::RFSparseCladoEvents) = events.nnz

function node_path_max_node(node_path::Vector{Vector{Int64}})
    max_node = 0
    for trav_node in node_path
        for node in trav_node
            max_node = max(max_node, node)
        end
    end
    return max_node
end

function make_child_branch_lengths(brs::Matrix{Float64}, node_path::Vector{Vector{Int64}})
    max_node = max(node_path_max_node(node_path), Int(maximum(brs[:, 2])))
    branch_lengths = Vector{Union{Nothing, Float64}}(nothing, max_node)

    for row in axes(brs, 1)
        child = Int(brs[row, 2])
        branch_lengths[child] = brs[row, 3]
    end

    return branch_lengths
end

function make_branch_transition_cache(Q::Matrix{Float64},
                                      brs::Matrix{Float64},
                                      node_path::Vector{Vector{Int64}})
    branch_lengths = make_child_branch_lengths(brs, node_path)
    return make_branch_transition_cache(Q, branch_lengths, node_path)
end

function make_branch_transition_cache(Q::Matrix{Float64},
                                      branch_lengths::Vector{Union{Nothing, Float64}},
                                      node_path::Vector{Vector{Int64}})
    transition_mats = Vector{Union{Nothing, Matrix{Float64}}}(nothing, length(branch_lengths))

    for trav_node in node_path
        left_child = trav_node[2]
        right_child = trav_node[3]

        if transition_mats[left_child] === nothing
            left_length = branch_lengths[left_child]
            left_length === nothing && error("No branch length found for child node $(left_child)")
            transition_mats[left_child] = exponential!(Q * left_length, ExpMethodHigham2005())
        end

        if transition_mats[right_child] === nothing
            right_length = branch_lengths[right_child]
            right_length === nothing && error("No branch length found for child node $(right_child)")
            transition_mats[right_child] = exponential!(Q * right_length, ExpMethodHigham2005())
        end
    end

    return RFBranchTransitionCache(transition_mats)
end
