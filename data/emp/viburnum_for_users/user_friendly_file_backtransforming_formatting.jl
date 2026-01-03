



function one_hot_encode(mat::Matrix{Int})
    n_rows = size(mat, 1)
    result = zeros(Int, n_rows, 3)
    for i in 1:n_rows
        for val in mat[i, :]
            if 1 ≤ val ≤ 3
                result[i, val] = 1
            end
        end
    end
    return result
end


function process_aff_to_user_files(tip_rf_ranges, forbidden_aff, tree)
    # Convert the tip_rf_ranges to a matrix
    rf_mat = reduce(vcat, permutedims(tip_rf_ranges))

    # Create a one-hot encoded matrix for the forbidden affinities
    zeros_mat = one_hot_encode(forbidden_aff)

    # Create matrices for twos and ones
    twos_mat = [x .== 2 for x in tip_rf_ranges]
    twos_mat = 1 .* (reduce(vcat, permutedims.(twos_mat)))

    ones_mat = [x .== 1 for x in tip_rf_ranges]
    ones_mat = 1 .* (reduce(vcat, permutedims.(ones_mat)))

    zeros_mat= one_hot_encode( forbidden_aff )


    # Example input (substitute your full data)
    data = tip_rf_ranges


    twos_mat = [x .== 2 for x in data]
    twos_mat = 1 .* (reduce(vcat, permutedims.(twos_mat)))

    # Similarly, for 1s
    ones_mat = [x .== 1 for x in data]
    ones_mat = 1 .*  reduce(vcat, permutedims.(ones_mat))


    function write_matrix_with_rowcol_names(filename, matrix, colnames::Vector{String}, rownames::Vector{String})
        open(filename, "w") do io
            # Write header: empty cell + column names
            println(io, join(["ID"; colnames], ","))

            # Write each row with row name prepended
            for (i, row) in enumerate(eachrow(matrix))
                println(io, join([rownames[i]; row], ","))
            end
        end
    end

    # Example column names
    colnames = ["Tropical", "Warm", "Cold"]


    rownames = tree.tlab  # Example row names, adjust as needed
    write_matrix_with_rowcol_names("data/emp/viburnum_for_vignette/"* "non_affinities.csv", Int.(zeros_mat), colnames, rownames)
    write_matrix_with_rowcol_names("data/emp/viburnum_for_vignette/" "enabled_affinities.csv", Int.(ones_mat), colnames, rownames)
    write_matrix_with_rowcol_names("data/emp/viburnum_for_vignette/" "established_affinities.csv", Int.(twos_mat), colnames, rownames)

    # For combined matrix:
    combined_array = Matrix{Union{Missing, String}}(missing, size(zeros_mat))

    for i in 1:size(zeros_mat, 1), j in 1:size(zeros_mat, 2)
        if zeros_mat[i, j] == 1
            combined_array[i, j] = "0"
        elseif ones_mat[i, j] == 1
            combined_array[i, j] = "1"
        elseif twos_mat[i, j] == 1
            combined_array[i, j] = "2"
        end
    end

    na_str_array =  [ismissing(x) ? "NA" : x for x in combined_array]

    write_matrix_with_rowcol_names("data/emp/viburnum_for_vignette/combined.csv", na_str_array, colnames, rownames)
