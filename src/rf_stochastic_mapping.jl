import Random

const SimmapBranch = Vector{Tuple{Int, Float64}}

struct SimmapEdge
    parent::Int
    child::Int
    length::Float64
end

function simmap_exit_rate(Q::AbstractMatrix{<:Real}, state::Int)
    total = 0.0
    for next_state in axes(Q, 2)
        if next_state != state && Q[state, next_state] > 0.0
            total += Float64(Q[state, next_state])
        end
    end
    return total
end

function simmap_reachable(Q::AbstractMatrix{<:Real}, start_state::Int, end_state::Int)
    start_state == end_state && return true

    seen = falses(size(Q, 1))
    stack = [start_state]
    seen[start_state] = true

    while !isempty(stack)
        state = pop!(stack)
        for next_state in axes(Q, 2)
            if next_state != state && Q[state, next_state] > 0.0 && !seen[next_state]
                next_state == end_state && return true
                seen[next_state] = true
                push!(stack, next_state)
            end
        end
    end

    return false
end

function draw_simmap_next_state(Q::AbstractMatrix{<:Real},
                                state::Int,
                                total_rate::Float64,
                                rng::Random.AbstractRNG)
    draw = rand(rng) * total_rate
    cumulative = 0.0

    for next_state in axes(Q, 2)
        if next_state != state && Q[state, next_state] > 0.0
            cumulative += Float64(Q[state, next_state])
            if draw <= cumulative
                return next_state
            end
        end
    end

    error("could not draw next SIMMAP state from state $state")
end

function push_simmap_segment!(history::SimmapBranch,
                              state::Int,
                              duration::Float64;
                              atol::Float64 = 1e-12)
    if duration < -atol
        error("negative SIMMAP segment duration $duration for state $state")
    end

    clipped_duration = duration < 0.0 ? 0.0 : duration
    if isempty(history)
        push!(history, (state, clipped_duration))
    elseif history[end][1] == state
        previous_state, previous_duration = history[end]
        history[end] = (previous_state, previous_duration + clipped_duration)
    elseif clipped_duration > atol
        push!(history, (state, clipped_duration))
    end

    return history
end

function check_simmap_inputs(Q::AbstractMatrix{<:Real},
                             start_state::Int,
                             end_state::Int,
                             branch_length::Real)
    size(Q, 1) == size(Q, 2) || error("Q must be square")
    start_state in axes(Q, 1) || error("start_state $start_state is outside Q")
    end_state in axes(Q, 1) || error("end_state $end_state is outside Q")
    branch_length >= 0.0 || error("branch_length must be nonnegative")
end

function gillespie_rejection_map!(history::SimmapBranch,
                                  Q::AbstractMatrix{<:Real},
                                  start_state::Int,
                                  end_state::Int,
                                  branch_length::Real;
                                  max_attempts::Int = 10_000,
                                  rng::Random.AbstractRNG = Random.GLOBAL_RNG,
                                  atol::Float64 = 1e-12)
    check_simmap_inputs(Q, start_state, end_state, branch_length)
    t = Float64(branch_length)

    if t <= atol
        start_state == end_state || error("zero-length branch cannot map $start_state to $end_state")
        empty!(history)
        push!(history, (start_state, 0.0))
        return history
    end

    simmap_reachable(Q, start_state, end_state) ||
        error("impossible SIMMAP branch: no positive-rate path from $start_state to $end_state")

    for _ in 1:max_attempts
        empty!(history)
        current_state = start_state
        elapsed = 0.0

        while elapsed < t - atol
            total_rate = simmap_exit_rate(Q, current_state)

            if total_rate <= atol
                push_simmap_segment!(history, current_state, t - elapsed; atol = atol)
                elapsed = t
                break
            end

            wait = -log(rand(rng)) / total_rate
            if elapsed + wait >= t - atol
                push_simmap_segment!(history, current_state, t - elapsed; atol = atol)
                elapsed = t
                break
            end

            push_simmap_segment!(history, current_state, wait; atol = atol)
            elapsed += wait
            current_state = draw_simmap_next_state(Q, current_state, total_rate, rng)
        end

        if current_state == end_state
            return history
        end
    end

    empty!(history)
    error("failed to sample SIMMAP branch from $start_state to $end_state over length $t after $max_attempts attempts")
end

function gillespie_rejection_map(Q::AbstractMatrix{<:Real},
                                 start_state::Int,
                                 end_state::Int,
                                 branch_length::Real;
                                 max_attempts::Int = 10_000,
                                 rng::Random.AbstractRNG = Random.GLOBAL_RNG,
                                 atol::Float64 = 1e-12)
    history = Tuple{Int, Float64}[]
    return gillespie_rejection_map!(history,
                                    Q,
                                    start_state,
                                    end_state,
                                    branch_length;
                                    max_attempts = max_attempts,
                                    rng = rng,
                                    atol = atol)
end

function build_simmap_buffer(tree::rtree)
    return [Tuple{Int, Float64}[] for _ in axes(tree.ed, 1)]
end

function simmap_branch_length_lookup(brs::AbstractMatrix{<:Real})
    lookup = Dict{Tuple{Int, Int}, Float64}()
    for row in axes(brs, 1)
        lookup[(Int(brs[row, 1]), Int(brs[row, 2]))] = Float64(brs[row, 3])
    end
    return lookup
end

function build_simmap_edges(tree::rtree, brs::AbstractMatrix{<:Real})
    branch_lengths = simmap_branch_length_lookup(brs)
    edges = SimmapEdge[]

    for edge_index in axes(tree.ed, 1)
        parent = tree.ed[edge_index, 1]
        child = tree.ed[edge_index, 2]
        length = get(branch_lengths, (parent, child), tree.el[edge_index])
        push!(edges, SimmapEdge(parent, child, length))
    end

    return edges
end

function fill_simmap_buffer!(buffer::Vector{SimmapBranch},
                             Q::AbstractMatrix{<:Real},
                             anc_states::AbstractVector{<:Integer},
                             simmap_edges::Vector{SimmapEdge};
                             max_attempts::Int = 10_000,
                             rng::Random.AbstractRNG = Random.GLOBAL_RNG,
                             atol::Float64 = 1e-12)
    n_nodes = length(anc_states) ÷ 2
    length(anc_states) == 2 * n_nodes || error("anc_states must be vcat(Nst, Cst)")
    length(buffer) == length(simmap_edges) || error("SIMMAP buffer and edge metadata lengths differ")

    for edge_index in eachindex(simmap_edges)
        edge = simmap_edges[edge_index]
        child = edge.child
        child <= n_nodes || error("child node $child is outside ancestral state record")

        start_state = Int(anc_states[n_nodes + child])
        end_state = Int(anc_states[child])
        start_state > 0 || error("missing post-cladogenesis start state for child node $child")
        end_state > 0 || error("missing end state for child node $child")

        gillespie_rejection_map!(buffer[edge_index],
                                 Q,
                                 start_state,
                                 end_state,
                                 edge.length;
                                 max_attempts = max_attempts,
                                 rng = rng,
                                 atol = atol)
    end

    return buffer
end

function format_simmap_float(x::Real)
    return string(round(Float64(x), sigdigits = 12))
end

function format_simmap_branch(history::SimmapBranch)
    isempty(history) && error("cannot format empty SIMMAP branch history")
    return "{" * join([string(state) * "," * format_simmap_float(duration) for (state, duration) in history], ":") * "}"
end

function simmap_children_by_parent(tree::rtree)
    children = Dict{Int, Vector{Int}}()
    for edge_index in axes(tree.ed, 1)
        parent = tree.ed[edge_index, 1]
        child = tree.ed[edge_index, 2]
        if !haskey(children, parent)
            children[parent] = Int[]
        end
        push!(children[parent], child)
    end
    return children
end

function simmap_root_node(tree::rtree)
    parents = Set(tree.ed[:, 1])
    children = Set(tree.ed[:, 2])
    roots = collect(setdiff(parents, children))
    length(roots) == 1 || error("expected one root node, found $(length(roots))")
    return first(roots)
end

function simmap_edge_index_by_child(simmap_edges::Vector{SimmapEdge})
    by_child = Dict{Int, Int}()
    for edge_index in eachindex(simmap_edges)
        by_child[simmap_edges[edge_index].child] = edge_index
    end
    return by_child
end

function format_simmap_newick(tree::rtree,
                              buffer::Vector{SimmapBranch},
                              simmap_edges::Vector{SimmapEdge})
    children_by_parent = simmap_children_by_parent(tree)
    edge_by_child = simmap_edge_index_by_child(simmap_edges)
    ntips = length(tree.tlab)
    root = simmap_root_node(tree)

    function render_node(node::Int)
        if haskey(children_by_parent, node)
            rendered_children = [render_node(child) for child in children_by_parent[node]]
            node_text = "(" * join(rendered_children, ",") * ")"
        elseif node <= ntips
            node_text = tree.tlab[node]
        else
            node_text = ""
        end

        if node == root
            return node_text
        end

        edge_index = edge_by_child[node]
        return node_text * ":" * format_simmap_branch(buffer[edge_index])
    end

    return render_node(root) * ";"
end

function simmap_log_header(tree::rtree)
    edge_columns = join(string.(1:size(tree.ed, 1)), "\t")
    return "Iteration\t" * edge_columns * "\tsimmap"
end

function simmap_log_row(iter::Integer,
                        tree::rtree,
                        buffer::Vector{SimmapBranch},
                        simmap_edges::Vector{SimmapEdge})
    branch_columns = join([format_simmap_branch(buffer[i]) for i in eachindex(buffer)], "\t")
    return string(iter) * "\t" * branch_columns * "\t" * format_simmap_newick(tree, buffer, simmap_edges)
end
