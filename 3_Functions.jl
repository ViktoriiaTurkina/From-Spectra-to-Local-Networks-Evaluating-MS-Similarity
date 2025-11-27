"function to get Vector from a string"
    function getVec(matStr)
        if matStr[1] .== '['
            if contains(matStr, ", ")
                str = split(matStr[2:end-1],", ")
            else
                str = split(matStr[2:end-1]," ")
            end
        elseif matStr[1] .== 'A'
            if contains(matStr, ", ")
                str = split(matStr[5:end-1],", ")
            else
                str = split(matStr[5:end-1]," ")
            end
        elseif matStr[1] .== 'F'
            if matStr .== "Float64[]"
                return []
            else
                str = split(matStr[9:end-1],", ")
            end
        else
            println("New vector start")
        end
        if length(str) .== 1 && cmp(str[1],"") .== 0
            return []
        else
            str = parse.(Float64, str)
            return str
        end
    end

""" spectral database_filtering"""
    function first_filtring(data_set,InChIKey,Resolution,MZ_VALUES, PRECURSOR_ION)
        #filtering N/A and missing values of InCIkeys and Resolution columns
        data_set = data_set[.!ismissing.(data_set[!,InChIKey]),:]
        data_set = data_set[(data_set[!,InChIKey].!= "N/A"),:]
        data_set = data_set[(data_set[!,InChIKey].!= "NA"),:]

        #filtering N/A and missing values of Resolution columns
        data_set = data_set[.!ismissing.(data_set[!,MZ_VALUES]),:]
        data_set = data_set[(data_set[!,MZ_VALUES].!= "[]"),:]
        
        #filtering N/A and missing values of PrecurcorIon 
        data_set = data_set[.!ismissing.(data_set[!,PRECURSOR_ION]),:]
        data_set = data_set[(data_set[!,PRECURSOR_ION].!= "N/A"),:]
        data_set = data_set[(data_set[!,PRECURSOR_ION].!= "NA"),:]
    

        #filtering resolution < 5000
        data_set= data_set[.!ismissing.(data_set[!,Resolution]),:]

        temp =[]
        for i = 1:length(data_set[!,Resolution])
            try
                if  typeof(data_set[i,Resolution]) == String
                    if parse(Int64,data_set[i,Resolution]) < 5000
                        push!(temp,i)
                    end 
                elseif typeof(data_set[i,Resolution]) == Int64
                    if data_set[i,Resolution] < 5000
                        push!(temp,i)
                    end  
                end


            catch
                push!(temp,i)
            end
        end

        deleteat!(data_set,temp)

        return data_set
    end

    """ function to 1.remove prec ions and 2. Remove spectra with less than X frags"""

    function remove_prec_ion(data_set, n_frags, tol, PRECURSOR_ION, MZ_VALUES, MZ_INT, MZ_INT_REL)
        
        locs2del = []
        for i in ProgressBar(1:size(data_set,1))
            prec_ion = data_set[i, PRECURSOR_ION]
            mz_vals = getVec(data_set[i, MZ_VALUES])
            mz_int = getVec(data_set[i, MZ_INT])

            if length(mz_vals) !== length(mz_int)
                push!(locs2del, i)
            else

                mz_vals_rm_prec = mz_vals[abs.(prec_ion.-mz_vals).>= tol]
                mz_int_rm_prec = mz_int[abs.(prec_ion.-mz_vals).>= tol]
                mz_rel_int_rm_prec = normalization(mz_int_rm_prec)

                if length(mz_vals_rm_prec) > n_frags
                    data_set[i, MZ_VALUES] = "Any[" * join(mz_vals_rm_prec, ", ") * "]"
                    data_set[i, MZ_INT] = "Any[" * join(mz_int_rm_prec, ", ") * "]"
                    data_set[i, MZ_INT_REL] = "Any[" * join(mz_rel_int_rm_prec, ", ") * "]"
                else
                    push!(locs2del, i)
                end
            end

        end

        return data_set[Not(locs2del),:]

    end

    function filtering_IonModeANDproton(data_set, Ion_mode,Pmass, Exact_mass)
        #filter positive (H+) and negative modes
        negative_ind = Int64[]
        positive_ind = Int64[]
        @views ion_mode::Vector{String} = data_set[!, Ion_mode]
        Pmass::Vector{Float64} = data_set[!,Pmass] 
        Exact_mass::Vector{Float64} = data_set[!,Exact_mass]
        Precursor_type = data_set[!,"PRECURSOR_TYPE"]

        for i::Int64 in ProgressBar(1:length(ion_mode))
            if (ion_mode[i] == "POSITIVE" || ion_mode[i] == "P") && (abs(Pmass[i] - Exact_mass[i]) < 1.01) && (Precursor_type[i] == "[M+H]+" || Precursor_type[i] == "[Cat-H]+")
                push!(positive_ind, i)
            elseif ion_mode[i] == "NEGATIVE" || ion_mode[i] == "N"
                push!(negative_ind, i)
            end
        end

        negative_set = data_set[negative_ind,:]
        positive_set = data_set[positive_ind,:]

        return positive_set, negative_set
    end

    function filtering_IonMode(data_set, Ion_mode)
        #filter positive (H+) and negative modes
        negative_ind = Int64[]
        positive_ind = Int64[]
        @views ion_mode::Vector{String} = data_set[!, Ion_mode]

        Precursor_type = data_set[!,"PRECURSOR_TYPE"]

        for i::Int64 in ProgressBar(1:length(ion_mode))
            if ion_mode[i] == "POSITIVE" || ion_mode[i] == "P"
                push!(positive_ind, i)
            elseif ion_mode[i] == "NEGATIVE" || ion_mode[i] == "N"
                push!(negative_ind, i)
            end
        end

        negative_set = data_set[negative_ind,:]
        positive_set = data_set[positive_ind,:]

        return positive_set, negative_set
    end
""" """


"""local networks"""

    function screen_dataset(INCHIKEY,PRECURSOR_ION, mz_tol)
            valid_entries = DataFrame(N = collect(1:length(INCHIKEY)),INCHIKEY = INCHIKEY, PRECURSOR_ION = PRECURSOR_ION)

        # Sort by mass
            sorted_entries = valid_entries[sortperm(valid_entries.PRECURSOR_ION), :]

        #generate empty output file
            output = DataFrame(PRECURSOR_ION = [], INCHIKEYS = [], n_INCHIKEY = [], n_SPECTRA = [], total_n_SPECTRA = [], locs = [])
            count_spectra = 0

            while count_spectra < length(sorted_entries.PRECURSOR_ION)

                mz = sorted_entries.PRECURSOR_ION[count_spectra+1]
                ind = (abs.(sorted_entries.PRECURSOR_ION .- mz) .<= mz_tol)
                network_i = sorted_entries[ind,:]
                
                n_SPECTRA  = [count(==(i),network_i.INCHIKEY) for i in unique(network_i.INCHIKEY)]

                output_i = DataFrame(
                    PRECURSOR_ION = mz, 
                    INCHIKEYS = join(unique(network_i.INCHIKEY),"; "), 
                    n_INCHIKEY = length(unique(network_i.INCHIKEY)), 
                    n_SPECTRA = join(n_SPECTRA,"; "), 
                    total_n_SPECTRA = length(network_i.INCHIKEY),
                    locs = join( network_i.N , "; ")
                    )

                append!(output,  output_i) 
                count_spectra += length(network_i.INCHIKEY)

            end
            
            return output

    end

    function normalization(int_spectra)
        rel_int_spectra = []

        for i in eachindex(int_spectra)
            rel_int_spectra_i = round.(int_spectra[i]./maximum(int_spectra), digits = 3)
            push!(rel_int_spectra, rel_int_spectra_i)
        end

        return rel_int_spectra

    end


    function star_network(similarity_matrix, inchikeys, central_nodes, thr)
        
        n = length(inchikeys)
        g_star = SimpleGraph(n)
            TP = 0
            FP = 0
            TN = 0
            FN = 0
            # Star network
            for central_node in central_nodes
                for j in 1:n
                    if j != central_node && similarity_matrix[central_node, j] ≥ thr
                        add_edge!(g_star, central_node, j)
                        if inchikeys[central_node] == inchikeys[j] 
                            TP +=1
                        elseif inchikeys[central_node] != inchikeys[j]
                            FP +=1
                        end
                    elseif j != central_node && similarity_matrix[central_node, j] < thr

                        if inchikeys[central_node] == inchikeys[j] 
                            FN +=1
                        elseif inchikeys[central_node] != inchikeys[j]
                            TN +=1
                        end
                    end


                end
            end

            #connect central nodes
            for i in 1:length(central_nodes)
                node_i = central_nodes[i]
                for j in i+1:length(central_nodes)  # only central-to-central without repeats
                    node_j = central_nodes[j]
                    if similarity_matrix[node_i, node_j] ≥ thr
                        add_edge!(g_star, node_i, node_j)
                        FP +=1
                    else 
                        TN +=1
                    end
                end
            end

            return g_star,  Dict("TP"=>TP, "FP"=>FP, "TN" =>TN, "FN"=>FN)
    end

    function star_network_simulation(similarity_matrix, inchikeys, central_nodes, thr)
        
        n = length(inchikeys)

            TP = 0
            FP = 0
            TN = 0
            FN = 0
            # Star network
            for central_node in central_nodes
                for j in 1:n
                    if j != central_node && similarity_matrix[central_node, j] ≥ thr

                        if inchikeys[central_node] == inchikeys[j] 
                            TP +=1
                        elseif inchikeys[central_node] != inchikeys[j]
                            FP +=1
                        end
                    elseif j != central_node && similarity_matrix[central_node, j] < thr

                        if inchikeys[central_node] == inchikeys[j] 
                            FN +=1
                        elseif inchikeys[central_node] != inchikeys[j]
                            TN +=1
                        end
                    end


                end
            end

            #connect central nodes
            for i in 1:length(central_nodes)
                node_i = central_nodes[i]
                for j in i+1:length(central_nodes)  # only central-to-central without repeats
                    node_j = central_nodes[j]
                    if similarity_matrix[node_i, node_j] ≥ thr

                        FP +=1
                    else 
                        TN +=1
                    end
                end
            end

            return  Dict("TP"=>TP, "FP"=>FP, "TN" =>TN, "FN"=>FN)
    end

    function select_central_nodes(similarity_matrix, inchikeys, thr)
        inchikey_indices = [findall(x -> x == i, inchikeys) for i in unique(inchikeys)]
        central_nodes = []

        for i in eachindex(inchikey_indices)

            inchikey_indeces = inchikey_indices[i]
            connected_nodes_i = [count(j.>=thr) for j in eachcol(similarity_matrix[:, inchikey_indeces])]
            max_val_i = maximum(connected_nodes_i )
            Random.seed!(1234) 
            central_node_i = inchikey_indeces[rand(findall(x -> x == max_val_i, connected_nodes_i))]

            push!(central_nodes, central_node_i)
        end
        
        return central_nodes
    end

    #3. build molecular network 
    function generic_network(similarity_matrix, n,thr)
        g_full = SimpleGraph(n)
        
            # Full network
            for i in 1:n-1
                for j in i+1:n
                    if similarity_matrix[i, j] ≥ thr
                        add_edge!(g_full, i, j)
                    end
                end
            end

            return g_full
    end


    function similarity_matrix_calc(mz_spectra_i, int_spectra_i, similarity; not_distance = true, need_normalize_result = true )
        #if not disctance (default)
        if not_distance
            similarity_functions = Dict(
                "cosine" => (similarity_function = cosine_similarity, is_index = true),
                "bi-cosine" => (similarity_function = binary_cosine_similarity, is_index = true),
                "jaccard" => (similarity_function = jaccard_similarity, is_index = true),
                "bi-jaccard" => (similarity_function = binary_jaccard_similarity, is_index = true),
                "entropy" => (similarity_function = unweighted_entropy_similarity_new, is_index = true),
                "w-entropy" => (similarity_function = entropy_similarity_new, is_index = true))   

            if !(similarity in keys(similarity_functions))
                error("Unsupported similarity type: $similarity_type. Supported types are: $(keys(similarity_functions))")
            end

            similarity_function = similarity_functions[similarity]

            n = length(mz_spectra_i)
            similarity_matrix = zeros(Float64, n, n)
            for i in 1:n-1

                for j in i:n
                
                    s1 = [mz_spectra_i[i] int_spectra_i[i]]
                    s2 = [mz_spectra_i[j] int_spectra_i[j]]
                    sim1 = similarity_function[1](s1, s2)
                    similarity_matrix[i, j] = sim1
                    similarity_matrix[j, i] = sim1
                end
            end
            return similarity_matrix

        #if distance
        elseif !not_distance
            n = length(mz_spectra_i)
            similarity_matrix = zeros(Float64, n, n)
            for i in 1:n-1
                for j in i:n
                    s1 = [mz_spectra_i[i] int_spectra_i[i]]
                    s2 = [mz_spectra_i[j] int_spectra_i[j]]
                    sim1 = distance(s1, s2, similarity,  ms2_da = 0.01; need_normalize_result = need_normalize_result)
                    similarity_matrix[i, j] = sim1
                    similarity_matrix[j, i] = sim1
                end
            end
            return similarity_matrix


        end

    end



    function match_peaks_in_spectra(spec_a::AbstractArray{T,2}, spec_b::AbstractArray{T,2}; ms2_ppm=nothing, ms2_da=nothing) where T <: Real
        # Convert to Float64 for calculations
    
        a_mz = Float64.(spec_a[:,1])
        a_i  = Float64.(spec_a[:,2])
        b_mz = Float64.(spec_b[:,1])
        b_i  = Float64.(spec_b[:,2])

        # determine tolerance function
        if ms2_da !== nothing
            tolfun = (mz1,mz2)->abs(mz1-mz2) <= ms2_da
        elseif ms2_ppm !== nothing
            tolfun = (mz1,mz2)->abs(mz1-mz2) <= mz2ppm(mz1, ms2_ppm)
        else
            error("MS2 tolerance need to be defined!")
        end

        # Build union of m/zs (sorted unique of all peaks)
        mz_union = sort(unique(vcat(a_mz,b_mz)))
        N = length(mz_union)
        res = Array{Float32,2}(undef, N, 3)

        for (idx, mzu) in enumerate(mz_union)
            # find nearest in a and b within tolerance
            ai = find_nearest_intensity(mzu, a_mz, a_i, ms2_ppm, ms2_da)
            bi = find_nearest_intensity(mzu, b_mz, b_i, ms2_ppm, ms2_da)
            res[idx,1] = Float32(mzu)
            res[idx,2] = Float32(ai)
            res[idx,3] = Float32(bi)
        end
        return res
    end

    # helper: find nearest intensity for a given mz in arrays, with tolerance
    function find_nearest_intensity(mz::Float64, mz_array::Vector{Float64}, inten_array::Vector{Float64}, ms2_ppm, ms2_da)
        if length(mz_array) == 0
            return 0.0
        end
        best = 1
        bestd = abs(mz - mz_array[1])
        for (i, mzv) in enumerate(mz_array)
            d = abs(mz - mzv)
            if d < bestd
                bestd = d; best = i
            end
        end
        if ms2_da !== nothing
            tol = ms2_da
        else
            tol = mz2ppm(mz, ms2_ppm)
        end
        if bestd <= tol
            return inten_array[best]
        else
            return 0.0
        end
    end



    function visualize_network(g::Graphs.SimpleGraph, node_names::Vector{String31}, inchi_keys::Vector{String31}, reference_indices::Vector{Any})
        @assert length(node_names) == length(inchi_keys) "node_names and inchi_keys must have the same length"

        nx, ny = 100, 100
        x = range(0, 1, length=nx)
        y = range(0, 1, length=ny)

        # Define a function that determines color based on x and y (creating a linear gradient)
        gradient_matrix = [xv + yv for xv in x, yv in y]  # Linear gradient

        # Normalize values to [0, 1] for proper color mapping
        gradient_matrix .= (gradient_matrix .- minimum(gradient_matrix)) ./ (maximum(gradient_matrix) - minimum(gradient_matrix))

        # Generate colormap for later color lookup
        colormap = cgrad([RGB(0.85,0.70,0.90),
                  RGB(0.70,0.85,0.95),
                  RGB(0.95,0.85,0.70),
                  RGB(0.75,0.90,0.75),
                  RGB(0.90,0.75,0.75)], categorical=true)

        # Sample unique InChI keys
        unique_inchikeys = unique(inchi_keys)

        # Create a mapping from InChI key to a unique number
        inchi_key_numbers = Dict(unique_inchikeys[i] => i for i in 1:length(unique_inchikeys))

        # Compute Levenshtein similarity matrix for unique InChI keys
        dist = Levenshtein()
        n_unique = length(unique_inchikeys)
        pairwise_similarity_unique = [compare(unique_inchikeys[i], unique_inchikeys[j], dist) for i in 1:n_unique, j in 1:n_unique]

        # Convert similarity to distance (1 - similarity)
        pairwise_distance = 1 .- pairwise_similarity_unique

        # Perform classical MDS to reduce the dimensionality to 2D
        coordinates = classical_mds(pairwise_distance, 2)

        # Normalize coordinates to [0, 1]
        min_val, max_val = minimum(coordinates), maximum(coordinates)
        normalized_coordinates = (coordinates .- min_val) ./ (max_val - min_val)

        # Create mapping from InChI keys to 2D coordinates
        inchi_to_coordinates = Dict(unique_inchikeys[i] => normalized_coordinates[:, i] for i in 1:n_unique)

        # Map InChI keys to colors
        # inchi_to_color = Dict()
        # for inchi in unique_inchikeys
        #     if haskey(inchi_to_coordinates, inchi)
        #         coord = inchi_to_coordinates[inchi]
        #         gradient_value = gradient_matrix[Int(round(coord[1] * (nx - 1))) + 1, Int(round(coord[2] * (ny - 1))) + 1]
        #         inchi_to_color[inchi] = get(colormap, gradient_value)
        #     else
        #         inchi_to_color[inchi] = :gray  # Fallback color if something goes wrong
        #     end
        # end
        
        inchi_to_color = Dict()
        for inchi in unique_inchikeys
            coord = inchi_to_coordinates[inchi]
            gradient_value = gradient_matrix[Int(round(coord[1] * (nx - 1))) + 1, Int(round(coord[2] * (ny - 1))) + 1]
            # Slightly desaturate colors for a "pastel" look
            inchi_to_color[inchi] = RGBA(get(colormap, gradient_value).r * 0.8 + 0.2,
                                        get(colormap, gradient_value).g * 0.8 + 0.2,
                                        get(colormap, gradient_value).b * 0.8 + 0.2,
                                        1.0)
        end


        # Assign colors to nodes based on their InChI key
        marker_colors = [get(inchi_to_color, inchi, :gray) for inchi in inchi_keys]

        # Get layout positions and apply scaling
        x_positions, y_positions = spring_layout_seed(g; seed=42)
        scaling_factor = 0.95
        x_positions .*= scaling_factor
        y_positions .*= scaling_factor

        # Plot edges
        p = Plots.plot(xlim=(-1.05, 1.05), ylim=(-1.05, 1.05), legend=false, grid=false, axis=false, dpi = 600)
        for e in edges(g)
            src_idx, dst_idx = src(e), dst(e)
            plot!([x_positions[src_idx], x_positions[dst_idx]], 
                [y_positions[src_idx], y_positions[dst_idx]], 
                line=:path, lw=0.5, color=:black, legend=false)
        end

        # Scatter plot for nodes
        # scatter!(x_positions, y_positions, 
        #         marker=:circle, markersize=10, 
        #         markercolor=marker_colors, 
        #         markerstrokewidth=0.5, hover=node_names)

        scatter!(x_positions, y_positions, 
        marker=:circle, markersize=10, 
        markercolor=marker_colors, 
        markerstrokewidth=0.5, legend = false)

        # Add text annotations for reference indices
        # for i in reference_indices
        #     inchi_key = inchi_keys[i]
        #     node_name = node_names[i]
        #     first_word = split(node_name)[1]

        #     centroid_x = x_positions[i]
        #     centroid_y = y_positions[i]

        #     annotation_x = centroid_x
        #     annotation_y = centroid_y + 0.10

        #     max_adjustments = 20
        #     adjustment_count = 0
        #     while is_overlapping(annotation_x, annotation_y, x_positions, y_positions) && adjustment_count < max_adjustments
        #         annotation_y += 0.01
        #         adjustment_count += 1
        #     end

        #     if annotation_y <= centroid_y
        #         annotation_y = centroid_y + 0.10
        #     end

        #     Plots.annotate!(p, annotation_x, annotation_y, Plots.text(first_word, :black, 8))
        # end

        # Add number annotations directly on the node
        text_size = 6
        for (i, key) in enumerate(inchi_keys)
            if haskey(inchi_key_numbers, key)
                number = inchi_key_numbers[key]
            else
                number = "?"  # Fallback in case of missing key
            end
            annotation_x = x_positions[i]
            annotation_y = y_positions[i]
            Plots.annotate!(p, annotation_x, annotation_y, Plots.text(string(number), :black, text_size, :center))
        end

        #  Add a legend mapping numbers → InChI keys
        legend_text = ["$(inchi_key_numbers[k]) → $(k)" for k in unique_inchikeys]
        Plots.plot!(p, legend=false)
        Plots.annotate!(p, -1.0, -1.1, Plots.text(join(legend_text, "\n"), :black, 6, :left))

        display(p)

        return p
    end


    #5 calculate jaccard index, Sensitivity and Specificity
    function calculate_metrics(conf_mx_res)
        #input is a dictinoary eith keys TP, FP, TN, FN
        jaccard_index = []
        speceficity = []
        sensitivity = []
        precision = []
        for i in eachindex(conf_mx_res["TP"])
            TP = conf_mx_res["TP"][i]
            TN = conf_mx_res["TN"][i]
            FP = conf_mx_res["FP"][i]
            FN = conf_mx_res["FN"][i]
            # Calculate jaccard index
            if TN + FP + FN > 0
                jaccard_index_i = TP / (TP + FP + FN)
            else
                jaccard_index_i = 0.0
            end
            push!(jaccard_index, round(jaccard_index_i, digits = 3))


            # Calculate Speceficity
            if TN + FP > 0
                speceficity_i = TN / (TN + FP)
            else
                speceficity_i = 0.0
            end
            push!(speceficity, round(speceficity_i, digits = 3))

            # Calculate Sensitivity
            if TP + FN > 0
                sensitivity_i = TP / (TP + FN)
            else
                sensitivity_i = 0.0
            end
            push!(sensitivity, round(sensitivity_i, digits = 3))

            # Calculate Precision
            if TP + FP > 0
                precision_i = TP / (TP + FP)
            else
                precision_i = 0.0
            end
            push!(precision, round(precision_i,digits = 3))
        end
        return jaccard_index, speceficity, sensitivity,  precision 
    end

    #4. iterate mol networks through thresholds
    function thr_iters_star_network(thr, similarity_matrix, inchikeys; vis = nothing)
        networks = []
        conf_mx_res = Array{Float64}(undef, length(thr), 4) 
        for i in eachindex(thr)

            central_nodes = select_central_nodes(similarity_matrix, inchikeys, thr[i])
            g_star, conf_mx = star_network(similarity_matrix, inchikeys, central_nodes, thr[i])
            conf_mx_res[i,:] = [conf_mx["TP"], conf_mx["FP"], conf_mx["TN"], conf_mx["FN"]]

            push!(networks, g_star)

            if !isnothing(vis)
                visualize_network(g_star, inchikeys, inchikeys,  central_nodes)
            end
        end
        conf_mx_res_dic = Dict(
            "TP"=>conf_mx_res[:,1],
            "FP"=>conf_mx_res[:,2],
            "TN"=>conf_mx_res[:,3],
            "FN"=>conf_mx_res[:,4])

        return  networks, conf_mx_res_dic

    end

    function thr_iters_star_network_simulation(thr, similarity_matrix, inchikeys; vis = nothing)
      
        conf_mx_res = Array{Float64}(undef, length(thr), 4) 
        for i in eachindex(thr)

            central_nodes = select_central_nodes(similarity_matrix, inchikeys, thr[i])
            conf_mx = star_network_simulation(similarity_matrix, inchikeys, central_nodes, thr[i])
            conf_mx_res[i,:] = [conf_mx["TP"], conf_mx["FP"], conf_mx["TN"], conf_mx["FN"]]

            if !isnothing(vis)
                visualize_network(g_star, inchikeys, inchikeys,  central_nodes)
            end
        end
        conf_mx_res_dic = Dict(
            "TP"=>conf_mx_res[:,1],
            "FP"=>conf_mx_res[:,2],
            "TN"=>conf_mx_res[:,3],
            "FN"=>conf_mx_res[:,4])

        return conf_mx_res_dic

    end

""" similarity functions """

    #__________________________________________________
    # Function to calculate Jaccard Index
    # Jaccard index or Tanimoto coefficient is utilized to identify the degree of overlap between two spectra
    # 1. Done through a binary fashion as it checks the total pool of available/present m/z fragments
    # 2. Then checks the presence of fragments in both spectra
    # 3. Divides the amount of fragments present in both by the amount of unique fragments.
    # Input: 2 spectra that require to be compared
    # Output: the similarity value between the compared spectra

    function binary_jaccard_similarity(spec1::AbstractArray{T,2}, spec2::AbstractArray{T,2},  tolerance = 0.01) where T


        # Extract m/z values from both spectra
        mz1 = spec1[:,1]
        mz2 = spec2[:,1]

        # Combine m/z values from both spectra and round to 2 decimal places
        all_mz = unique(round.(vcat(mz1, mz2), digits=2))

        # Return 0.0 if there are no unique m/z channels
        if isempty(all_mz)
            return 0.0
        end

        # Create binary vectors for presence/absence of each m/z in both spectra
        binary1 = [any(abs(round(mz - mz1_val, digits = 2)) <= tolerance for mz1_val in mz1) for mz in all_mz]
        binary2 = [any(abs(round(mz - mz2_val, digits = 2)) <= tolerance for mz2_val in mz2) for mz in all_mz]

        # Calculate intersection and union of the binary vectors
        intersection = sum(binary1[i] && binary2[i] for i in 1:length(all_mz))
        union = sum(binary1[i] || binary2[i] for i in 1:length(all_mz))

        # Return the Jaccard Index (intersection/union)
        return intersection / union
    end

    #__________________________________________________
    # Jaccard with intensities
    # Works similar as jaccard index by identifying the overlapping area/total area.
    # Overlapping area is found by the minimal intensities of all common m/z.
    # Total area is calculated as the maximum intensities for every m/z.
    # Note: dividing 0 by any value is mathematically impossible and taken into account to prevent erros.
    # Input: 2 spectra that require to be compared
    # Output: the similarity value between the compared spectra

    function jaccard_similarity(spec1::AbstractArray{T,2}, spec2::AbstractArray{T,2}; tolerance = 0.01) where T
        # Define tolerance for m/z matching

        # Extract m/z and intensity values
        mz1, intensity1 = spec1[:,1], spec1[:,2]
        mz2, intensity2 = spec2[:,1], spec2[:,2]

        # Sort by m/z for efficient searching
        sorted_indices1 = sortperm(mz1)
        sorted_indices2 = sortperm(mz2)
        
        mz1, intensity1 = mz1[sorted_indices1], intensity1[sorted_indices1]
        mz2, intensity2 = mz2[sorted_indices2], intensity2[sorted_indices2]

        # Initialize intersection and union sums
        intersection = 0.0
        union = 0.0

        # Use two-pointer technique for efficient matching
        i, j = 1, 1
        while i ≤ length(mz1) && j ≤ length(mz2)
            diff = round(mz1[i] - mz2[j],digits = 2)
            if abs(diff) <= tolerance
                # Peaks match within tolerance
                min_intensity = min(intensity1[i], intensity2[j])
                max_intensity = max(intensity1[i], intensity2[j])
                
                intersection += min_intensity
                union += max_intensity
                
                i += 1
                j += 1  # Move both pointers
            elseif diff < 0
                # mz1[i] is smaller, so process it separately
                union += intensity1[i]
                i += 1
            else
                # mz2[j] is smaller, so process it separately
                union += intensity2[j]
                j += 1
            end
        end

        # Add remaining intensities from unmatched peaks
        while i ≤ length(mz1)
            union += intensity1[i]
            i += 1
        end
        while j ≤ length(mz2)
            union += intensity2[j]
            j += 1
        end

        # Return the Jaccard Index (weighted intersection/union)
        return union == 0.0 ? 0.0 : intersection / union
    end

    function normalize_distance(dist::Real, dist_range)
        lo = float(dist_range[1])
        hi = float(dist_range[2])
        if !isfinite(dist)
            return 1.0
        end
        if isfinite(hi)
            if hi == lo
                return 0.0
            end
            x = (dist - lo) / (hi - lo)
            return clamp(x, 0.0, 1.0)
        else
            # hi is Inf: map using a normalized saturating function
            if dist <= lo
                return 0.0
            else
                # use 1 - lo/dist as a smooth mapping (bounded in (0,1])
                val = 1.0 - lo / (dist + lo + EPS)
                return clamp(val, 0.0, 1.0)
            end
        end
    end


    function entropy_similarity_new(spec_a, spec_b; tolerance = 0.01)
        #spec_a = Float32[69.071 7.917962; 86.066 1.021589; 86.0969 100.0]
        #spec_b   = Float32[41.04 37.16; 69.07 66.83; 86.1 999.0]

        # Calculate entropy similarity.
        return me.calculate_entropy_similarity(spec_a, spec_b, ms2_tolerance_in_da = tolerance)
    end
    
    function unweighted_entropy_similarity_new(spec_a, spec_b; tolerance = 0.01)
        #input : Array 1st column mz 2nd column Intentsities 
        #spec_a = Float32[69.071 7.917962; 86.066 1.021589; 86.0969 100.0]
        #spec_b   = Float32[41.04 37.16; 69.07 66.83; 86.1 999.0]
        # Calculate unweighted entropy similarity.
        return me.calculate_unweighted_entropy_similarity(spec_a, spec_b, ms2_tolerance_in_da = tolerance)
    end


    function ms2deepscore_similarity_matrix(mz_networks_i, similarity_measure)


        mz_spectra_i = getVec.(mz_networks_i.MZ_VALUES)
        int_spectra_i = getVec.(mz_networks_i.MZ_INT_REL)
        inchikeys = mz_networks_i.INCHIKEY
        prec_ion_spectra_i = mz_networks_i.PRECURSOR_ION
        ion_mode_i = mz_networks_i.ION_MODE
        failed = []
        spectra = []
        for i in eachindex(mz_spectra_i)
            try
                spectrum_i =  matchms_spectrum.Spectrum(mz=np.array(mz_spectra_i[i]),
                                    intensities=np.array(int_spectra_i[i]),
                                    metadata=Dict("id" => inchikeys[i],
                                    "precursor_mz" => prec_ion_spectra_i[i],
                                    "IONMODE"=>ion_mode_i[i]))

                push!(spectra, spectrum_i)
            catch 
                println("error spectrum_i i = $i")
                push!(failed,i)
                continue
            end

        end
        
        scores = matchms_spectrum.calculate_scores(spectra, spectra, similarity_measure,
                                    is_symmetric=true # Set to False if you match to a library
                                    )
        return scores.to_array(), failed


    end


    # Function to calculate the cosine score similarity
    # Used in a for loop to compare multiple spectra
    # 1. Takes m/z and intensities of the cleaned spectra and puts these in a dictionary
    # 2. Finds the common m/z which are present in both spectra by a 10 ppm tolerance
    # 3. Whenever no common m/z are found the score is returned with a 0.0 similarity
    # 4. Calculating the dot product of the found common fragments based on their intensities
    # 5. Scales the cosine similarity by the magnitude of the vectors such that it always returns between 0 and 1
    # Input: 2 spectra that require to be compared
    # Output: the similarity value between the compared spectra

    # The only similarity metric that prevents fragments to be double matched. This was required as in one network there were spectra with
    # roughly 100 fragments in a spectra. As a result of this some fragments were matched twice, hence similarities >1 were obtained.
    # function is modified as such that this is now fully prevented. This modification seemed unnecesarry for other similarity measures.
        # spec1 = [ 31.02  0.403
        #         39.02  0.052
        #         41.04  1.0
        #         43.01  0.047
        #         44.02  0.116]

        #                 spec2 = [ 31.02  0.403
        #         39.02  0.052
        #         41.04  1.0
        #         43.01  0.047
        #         44.02  0.116]
                

        # spec2 =  [31.02  0.053
        #         39.02  0.189
        #         41.04  0.057
        #         42.03  0.02
        #         43.02  1.0
        #         44.03  0.258]
        # tolerance = 0.01
    function spectra_alignment(spec1::AbstractArray{T,2}, spec2::AbstractArray{T,2}; tolerance = 0.01)  where T

        # Convert to Float64 for calculations
    
        a_mz = Float64.(spec1[:,1])
        a_i  = Float64.(spec1[:,2])
        b_mz = Float64.(spec2[:,1])
        b_i  = Float64.(spec2[:,2])

        # Build union of m/zs (sorted unique of all peaks)

        mz_all = sort(unique(vcat(a_mz,b_mz)))
        used_mz2 = []
        mz_union = []
        for mz_val in mz_all
            if mz_val ∉ used_mz2
                diff = abs.(round.(mz_val .- mz_all,digits = 2))
                close_matches =  diff .<= tolerance
            
                if !isempty(close_matches)
                    best_match = mz_all[argmin(diff)]
                    mz_all[close_matches]

                    used_mz2 = [used_mz2; mz_all[close_matches]]  # Prevent double matching
                    push!(mz_union, best_match)
                    
                end
            else
                continue
            end

        end

        N = length(mz_union)
        res = Array{Float64,2}(undef, N, 3)

        for (idx, mzu) in enumerate(mz_union)
            # find nearest in a and b within tolerance
            ai = find_nearest_intensity(mzu, a_mz, a_i, tolerance)
            bi = find_nearest_intensity(mzu, b_mz, b_i, tolerance)
            res[idx,1] = Float64(mzu)
            res[idx,2] = Float64(ai)
            res[idx,3] = Float64(bi)
        end


        return res
    end

    function find_nearest_intensity(mz::Float64, mz_array::Vector{Float64}, inten_array::Vector{Float64}, tol)
        if length(mz_array) == 0
            return 0.0
        end
        best = 1
        bestd = abs(mz - mz_array[1])
        for (i, mzv) in enumerate(mz_array)
            d = abs(mz - mzv)
            if d < bestd
                bestd = d; best = i
            end
        end

        if bestd <= tol
            return inten_array[best]
        else
            return 0.0
        end
    end


    function cosine_similarity(spec1::AbstractArray{T,2}, spec2::AbstractArray{T,2}; tolerance = 0.01) where T

        aligned_spectra = spectra_alignment(spec1, spec2; tolerance = tolerance)

        # Compute dot product with strict one-to-one matching
        dot_product = sum( aligned_spectra[:,2] .* aligned_spectra[:,3])

        norm1 = sqrt(sum(v^2 for v in aligned_spectra[:,2]))
        norm2 = sqrt(sum(v^2 for v in aligned_spectra[:,3]))

        return round(dot_product / max(norm1* norm2), digits = 3)  # Prevent division by zero
    end



    #__________________________________________________
    # Function to calculate cosine similarity in a binary fashion 
    # 1. Utilizes a master vector containing all unique m/z fragments
    # 2. Makes two binary vectors for each spectrum specific
    # 3. Calculates the DOT product between the two binary vectors
    # Input: 2 spectra that require to be compared
    # Output: the similarity value between the compared spectra

    function binary_cosine_similarity(spec1::AbstractArray{T,2}, spec2::AbstractArray{T,2}; tolerance = 0.01) where T
        # Define tolerance for m/z matching

        # Extract m/z values from both spectra
        mz1 = spec1[:,1]
        mz2 = spec2[:,1]

        # Combine m/z values from both spectra and round to 2 decimal places
        all_mz = unique(round.(vcat(mz1, mz2), digits=2))

        # Return 0.0 if there are no unique m/z channels
        if isempty(all_mz)
            return 0.0
        end

        # Create binary vectors for presence/absence of each m/z in both spectra
        binary1 = [any(abs(round(mz - mz1_val,digits = 2)) <= tolerance for mz1_val in mz1) for mz in all_mz]
        binary2 = [any(abs(round(mz - mz2_val,digits = 2)) <= tolerance for mz2_val in mz2) for mz in all_mz]

        # Calculate the dot product of the binary vectors
        dot_product = sum(binary1[i] * binary2[i] for i in 1:length(all_mz))

        # Calculate the norms of the binary vectors
        norm1 = sqrt(sum(binary1[i]^2 for i in 1:length(all_mz)))
        norm2 = sqrt(sum(binary2[i]^2 for i in 1:length(all_mz)))

        # Return cosine similarity
        return dot_product / (norm1 * norm2)
    end


""" plotting result"""

    function res_vec2mat(res)
        JI_mat = reduce(hcat, [parse.(Float64, vcat(split.(res.JI[i], "; ")...)) for i in eachindex(res.JI)])'

        specificity_mat = reduce(hcat, [parse.(Float64, vcat(split.(res.specificity[i], "; ")...)) for i in eachindex(res.specificity)])'

        sensitivity_mat = reduce(hcat, [parse.(Float64, vcat(split.(res.sensitivity[i], "; ")...)) for i in eachindex(res.sensitivity)])'

        precision_mat = reduce(hcat, [parse.(Float64, vcat(split.(res.precision[i], "; ")...)) for i in eachindex(res.precision)])'

        accuracy_mat = reduce(hcat, [parse.(Float64, vcat(split.(res.Accuracy[i], "; ")...)) for i in eachindex(res.Accuracy)])'

        f1score_mat = reduce(hcat, [parse.(Float64, vcat(split.(res.F1score[i], "; ")...)) for i in eachindex(res.F1score)])'

        FDR_mat = reduce(hcat, [parse.(Float64, vcat(split.(res.FDR[i], "; ")...)) for i in eachindex(res.FDR)])'

        return JI_mat, specificity_mat, sensitivity_mat, precision_mat, accuracy_mat, f1score_mat, FDR_mat
    end
    

    function plot_bar_heatmap_all(path2res, path2save)
        
        res = CSV.read( path2res ,DataFrame)
        metric = split(path2res,"_")[end][1:end-4]

        JI_mat, specificity_mat, sensitivity_mat, precision_mat = res_vec2mat(res)

        plot_bar_heatmap(JI_mat, "Jaccard_Index", metric,path2save)
        plot_bar_heatmap(specificity_mat, "Specificity", metric,path2save)
        plot_bar_heatmap(sensitivity_mat, "Sensitivity", metric,path2save)
        plot_bar_heatmap(precision_mat,"Precision", metric,path2save)

    end

    function plot_bar_heatmap(mat, mat_name, metric,path2save)
        n_bin = 10
        ybins = collect(0:1/15:1.1)
        #boxplot
        display(sp.boxplot(mat, label = false, title = "$(mat_name)_$metric",dpi = 600))
        sp.savefig(joinpath(path2save,"boxplot_$(mat_name)_$metric.png"))

        # Count values per bin for each column
        counts = zeros(Int, length(ybins)-1, size(mat, 2))
        for j in 1:size(mat, 2)
            h = fit(Histogram, mat[:, j], ybins)
            counts[:, j] .= h.weights
        end

        bin_centers = (ybins[1:end-1] .+ ybins[2:end]) ./ 2
        display(sp.heatmap(1:size(mat, 2), bin_centers, counts,
                xlabel = "Column", ylabel = "Value bins", colorbar_title = "Count", title = "$(mat_name)_$metric", dpi = 600))

        sp.savefig(joinpath(path2save,"heatmap_$(mat_name)_$metric.png"))

    end
        