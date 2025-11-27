
include("1_SpectralSimilarity.jl")
include("2_imports.jl")
include("3_Functions.jl")




""" filter initial database"""

    path2dtb = "C:/Users/vturkin/OneDrive - UvA/Desktop/Project/2_Projects/Databases"
    internal_dtb = CSV.read(joinpath(path2dtb,"Database_INTERNAL_2022-11-17.csv"), DataFrame)

    data_set = deepcopy(internal_dtb)
    InChIKey = :INCHIKEY
    Resolution = :RESOLUTION
    MZ_VALUES = :MZ_VALUES
    PRECURSOR_ION = :PRECURSOR_ION
    MZ_INT = :MZ_INT
    MZ_INT_REL = :MZ_INT_REL
    Ion_mode = :ION_MODE


    #filter on resolution and missing values
    data_set_first_filt = first_filtring(data_set,InChIKey,Resolution,MZ_VALUES, PRECURSOR_ION)
    positive_data_set, negative_data_set = filtering_IonMode(data_set_first_filt, Ion_mode)

    positive_data_set.MZ_INT_REL = Vector{String}(positive_data_set.MZ_INT_REL)
    data_set = deepcopy(positive_data_set)
    # filter on n of fragments and remove precursor_ion
    n_frags = 3
    tol = 0.01
    data_set_precion_filter = remove_prec_ion(data_set, n_frags, tol, PRECURSOR_ION, MZ_VALUES, MZ_INT, MZ_INT_REL)

    data_set_precion_filter.ION_MODE[startswith.("P",data_set_precion_filter.ION_MODE)] .= "positive"
    data_set_precion_filter.ION_MODE[startswith.("N",data_set_precion_filter.ION_MODE)] .= "negative"
    data_set_precion_filter.ION_MODE[startswith.("POSITIVE",data_set_precion_filter.ION_MODE)] .= "positive"
    data_set_precion_filter.ION_MODE[startswith.("NEGATIVE",data_set_precion_filter.ION_MODE)] .= "negative"

    #save filtered file
    CSV.write("Database_positive_filtered.csv",data_set_precion_filter)
"""  """
"""n_mol_network in dataset and mz (from here iterations on mx tol can start)"""

    data_set =  CSV.read("Database_positive_filtered.csv",DataFrame)
    # Define input values
    INCHIKEY = data_set[:,:INCHIKEY]     # Input for Inhcl keys based on filtered high res data
    PRECURSOR_ION = data_set[:, :PRECURSOR_ION]    # Input for Exact masses of molecules within the database

    # Search for masses in the database for networks
    mz_tol = 0.01             # Define the mass window -> 10mDa
    n_INCHIKEY = 2                       # Assings minimal amount (> or =) of Inchl keys related to the masses
    n_spec = 5                  # Assings minimal amount (> or =) of spectra related to the masses

    screened_dataset = screen_dataset(INCHIKEY,PRECURSOR_ION, mz_tol)

    CSV.write("Screened_dataset_0.01tol_positive.csv", screened_dataset)
    screened_dataset = CSV.read("Screened_dataset_0.01tol_positive.csv", DataFrame)
    #13638 - 8290 = 5348

    mz_networks_singleINCHIKEY = screened_dataset[[(screened_dataset.n_INCHIKEY[i]==1) for i in eachindex(screened_dataset.n_INCHIKEY)],:]
    minimum(parse.(Int64,mz_networks_singleINCHIKEY.n_SPECTRA))
    maximum(parse.(Int64,mz_networks_singleINCHIKEY.n_SPECTRA))

    mz_networks = screened_dataset[[(screened_dataset.n_INCHIKEY[i]>=n_INCHIKEY && screened_dataset.total_n_SPECTRA[i]>=n_spec) for i in eachindex(screened_dataset.n_INCHIKEY)],:]
    CSV.write("mz_networks_positive_4inchi_9spcs.csv", mz_networks)



""" investigating subdatasets """
    mz_networks_full = CSV.read("mz_networks_positive.csv", DataFrame)

    mz_networks_full

    coount_n_inshikey = [count(mz_networks_full.n_INCHIKEY.== i) for i in unique(mz_networks_full.n_INCHIKEY)]

    length(mz_networks_full.n_INCHIKEY)

    count_n = coount_n_inshikey[sortperm(coount_n_inshikey, rev = true)]
    un_inchikey = unique(mz_networks_full.n_INCHIKEY)[sortperm(coount_n_inshikey, rev = true)]

    n_spectra_n_inchikey = mz_networks_full.total_n_SPECTRA[mz_networks_full.n_INCHIKEY.==2]
    precion_n_inchikey = mz_networks_full.PRECURSOR_ION[mz_networks_full.n_INCHIKEY.==3]



    minimum(count_n)
    mean(count_n)
    median(count_n)
    maximum(un_inchikey)

    sp.bar(
        coount_n_inshikey[sortperm(coount_n_inshikey, rev = true)],
        xticks = (1:length(unique(mz_networks_full.n_INCHIKEY)),
                unique(mz_networks_full.n_INCHIKEY)[sortperm(coount_n_inshikey, rev = true)]),
        xrotation = 90,
        size = (1200, 800),
        xlabel = "Number of InChIKeys per network",
        label = false,
        ylabel = "Occurence",
        bottom_margin = 10mm, 
        left_margin = 10mm, 
        dpi = 600
    )

    sp.savefig("5_results/figures_tol_mz001_positive/inchikey_barplot_positive.png")


    count(count_n.==1)
    maximum(un_inchikey)



""" build single mol network (iteration on single mz tol start)""" 
    data_set =  CSV.read("Database_INTERNAL_2022-11-17_filtered.csv",DataFrame) 



    mz_networks = CSV.read("mz_networks.csv", DataFrame)

    locs_spectra = parse.(Int64, vcat(split.(mz_networks.locs[1], "; ")...))
    mz_networks_i = data_set[locs_spectra,:]
    mz_spectra_i = getVec.(mz_networks_i.MZ_VALUES)
    int_spectra_i = getVec.(mz_networks_i.MZ_INT_REL)

    #find centroid:
    #1. first calculate similarity between all spectra
        # similarity matrix calculations
        similarity = "jaccard"
        similarity_matrix = similarity_matrix_calc(mz_spectra_i, int_spectra_i, similarity; not_distance = false)
        #similarity_matrix = similarity_matrix_calc(mz_spectra_i, int_spectra_i, similarity; not_distance = false, need_normalize_result = true)
        #similarity_matrix = similarity_matrix_calc1(mz_spectra_i, int_spectra_i, similarity; not_distance = true, need_normalize_result = false)

    # spec_a = [31.02 39.02 41.04 43.01 44.02
    #  0.403 0.052 1.0 0.047 0.116]'

    #  spec_b = [ 31.02 39.02 41.04 42.03 43.02 44.03
    #  0.053 0.189 0.057 0.02 1.0 0.258]'
    #  ms2_da = 0.01


    #2. find the most connected ones based on thr ans the number of connected nodes
        thr = 0.1
        inchikeys = mz_networks_i.INCHIKEY 
        central_nodes = select_central_nodes(similarity_matrix, inchikeys, thr)

    #3. build start network
        g_star, conf_mx = star_network(similarity_matrix, inchikeys, central_nodes, thr)
        visualize_network(g_star, inchikeys, inchikeys,  central_nodes)
        n = length(inchikeys)
        gplot(g_star, nodelabel=1:n)


    #4. iterate mol networks through thresholds
        if not_distance
            thr = collect(0:0.1:1)
            networks, conf_mx_res = thr_iters_star_network(thr, similarity_matrix, inchikeys)
        else
            thr = collect(0:0.1:1) .-1
            networks, conf_mx_res = thr_iters_star_network(thr, (similarity_matrix).*-1, inchikeys)
        end

    #5 calculate jaccard index, Sensitivity and Specificity
       
        jaccard_index, specificity, sensitivity,  precision  = calculate_metrics(conf_mx_res)

        Plots.plot(thr, jaccard_index, label = "jaccard_index", marker = :o)
        Plots.plot!(thr, specificity, label = "specificity", marker = :o)
        Plots.plot!(thr, sensitivity, label = "sensitivity", marker = :o)
        Plots.plot!(thr, precision, label = "precision", marker = :o)
        Plots.plot!(legend = :outerright)
    
""" batch calculations jaccard, cosine, entropy"""
    data_set =  CSV.read("Database_positive_filtered.csv",DataFrame)
    mz_networks_full = CSV.read("mz_networks_positive.csv", DataFrame)
    mz_networks = mz_networks_full[:,:]

    similarity_list = [ "entropy", "w-entropy"]
    similarity_list = ["jaccard", "bi-jaccard", "cosine", "bi-cosine"]
    similarity_list = ["bi-jaccard"]
    thr = collect(0:0.1:1)

    for similarity in similarity_list
        println(similarity)
        res = DataFrame()
        for i in ProgressBar(1:size(mz_networks,1))
         
                locs_spectra_i = parse.(Int64, vcat(split.(mz_networks.locs[i], "; ")...))
                mz_networks_i = data_set[locs_spectra_i,:]
                mz_spectra_i = getVec.(mz_networks_i.MZ_VALUES)
                int_spectra_i = getVec.(mz_networks_i.MZ_INT_REL)
                inchikeys = mz_networks_i.INCHIKEY


                #1. first calculate similarity between all spectra
                    similarity_matrix = similarity_matrix_calc(mz_spectra_i, int_spectra_i, similarity, not_distance = true)

                #2. iterate mol networks through thresholds
                    conf_mx_res = thr_iters_star_network_simulation(thr, similarity_matrix, inchikeys)

                #5 calculate jaccard index, Sensitivity and Specificity
                
                    jaccard_index, specificity, sensitivity,  precision  = calculate_metrics(conf_mx_res)

                #6 write results

                    res_i = DataFrame( 
                        TP = join(conf_mx_res["TP"], "; "),
                        FP = join(conf_mx_res["FP"], "; "),
                        TN = join(conf_mx_res["TN"], "; "),
                        FN = join(conf_mx_res["FN"], "; "),
                        JI = join(jaccard_index,"; "), 
                        specificity = join(specificity,"; "), 
                        sensitivity = join(sensitivity,"; "),  
                        precision = join(precision,"; "))
                    
                    append!(res, res_i)

        end

        CSV.write("res_$similarity.csv", res; bufsize = 16777216)

    end
    
        Plots.plot(thr, jaccard_index, label = "jaccard_index", marker = :o)
        Plots.plot!(thr, speceficity, label = "specificity", marker = :o)
        Plots.plot!(thr, sensitivity, label = "sensitivity", marker = :o)
        Plots.plot!(thr, precision, label = "precision", marker = :o)
        Plots.plot!(legend = :outerright)
    


""" plotting results """

    res_files = readdir("5_results/tol_mz001_allsims")[contains.(readdir("5_results/tol_mz001_allsims"),".csv")]
    path2save = "5_results/tol_mz001_allsims/figures"

    for i in ProgressBar(1:length(res_files)) 
        path2res = joinpath("5_results/tol_mz001_allsims",res_files[i])
        plot_bar_heatmap_all(path2res, path2save)

    end


""" """
""" batch calculations ms2deepscore"""
    data_set =  CSV.read("Database_positive_filtered.csv",DataFrame)
    mz_networks_full = CSV.read("mz_networks_positive.csv", DataFrame)
    mz_networks = mz_networks_full[:,:]
        # Load pretrained model
    model = ms2deepscore_model.load_model("ms2deepscore_model.pt")
    similarity_measure = sim_measure.MS2DeepScore(model)
   
   
    similarity = "ms2deepscore"
    thr = collect(0:0.1:1)
    

    res = DataFrame()


    for i in ProgressBar(1:size(mz_networks,1))
            locs_spectra_i = parse.(Int64, vcat(split.(mz_networks.locs[i], "; ")...))
            mz_networks_i = data_set[locs_spectra_i,:]
            mz_spectra_i = getVec.(mz_networks_i.MZ_VALUES)
            int_spectra_i = getVec.(mz_networks_i.MZ_INT_REL)
            inchikeys = mz_networks_i.INCHIKEY

        #1. first calculate similarity between all spectra
            similarity_matrix,failed = ms2deepscore_similarity_matrix(mz_networks_i, similarity_measure)

        #2. iterate mol networks through thresholds
            networks, conf_mx_res = thr_iters_star_network(thr, similarity_matrix, inchikeys[Not(failed)])

        #5 calculate jaccard index, Sensitivity and Specificity
        
            jaccard_index, specificity, sensitivity,  precision  = calculate_metrics(conf_mx_res)

        #6 write results

            res_i = DataFrame(
                MolNet = join(networks, "; "), 
                TP = join(conf_mx_res["TP"], "; "),
                FP = join(conf_mx_res["FP"], "; "),
                TN = join(conf_mx_res["TN"], "; "),
                FN = join(conf_mx_res["FN"], "; "),
                JI = join(jaccard_index,"; "), 
                specificity = join(specificity,"; "), 
                sensitivity = join(sensitivity,"; "),  
                precision = join(precision,"; "))
            
            append!(res, res_i)

    end

    CSV.write("res_jaccard_$similarity.csv", res; bufsize = 16777216)


    n_bin = 10
    ybins = collect(0:1/15:1.1)
    res
    # Count values per bin for each column
    counts = zeros(Int, length(ybins)-1, size(precision_mat, 2))
    for j in 1:size(precision_mat, 2)
        h = fit(Histogram, precision_mat[:, j], ybins)
        counts[:, j] .= h.weights
    end

    bin_centers = (ybins[1:end-1] .+ ybins[2:end]) ./ 2
    sp.heatmap(1:size(mat, 2), bin_centers, counts,
            xlabel = "Column", ylabel = "Value bins", colorbar_title = "Count", title = "precision")




""" batch calculations distances"""
    data_set =  CSV.read("Database_positive_filtered.csv",DataFrame);
    mz_networks_full = CSV.read("mz_networks_positive.csv", DataFrame);
    mz_networks = mz_networks_full[:,:];

    similarity_list = collect(keys(math_funcs));

    thr = round.(collect(0:0.1:1) .-1,digits = 2);
    similarity = similarity_list[12];

    for similarity in similarity_list
        @time begin
            res = DataFrame()
            for i in ProgressBar(1:10)
            #for i in ProgressBar(1:size(mz_networks,1))
                
                    locs_spectra_i = parse.(Int64, vcat(split.(mz_networks.locs[i], "; ")...))
                    mz_networks_i = data_set[locs_spectra_i,:]
                    mz_spectra_i = getVec.(mz_networks_i.MZ_VALUES)
                    int_spectra_i = getVec.(mz_networks_i.MZ_INT_REL)
                    inchikeys = mz_networks_i.INCHIKEY

                #1. first calculate similarity between all spectra
                    similarity_matrix = similarity_matrix_calc(mz_spectra_i, int_spectra_i, similarity, not_distance = false)

                #2. iterate mol networks through thresholds
                    networks, conf_mx_res = thr_iters_star_network(thr, (similarity_matrix)*-1, inchikeys)

                #5 calculate jaccard index, Sensitivity and Specificity
                
                    jaccard_index, specificity, sensitivity,  precision  = calculate_metrics(conf_mx_res)

                #6 write results

                    res_i = DataFrame(
                        MolNet = join(networks, "; "), 
                        TP = join(conf_mx_res["TP"], "; "),
                        FP = join(conf_mx_res["FP"], "; "),
                        TN = join(conf_mx_res["TN"], "; "),
                        FN = join(conf_mx_res["FN"], "; "),
                        JI = join(jaccard_index,"; "), 
                        specificity = join(specificity,"; "), 
                        sensitivity = join(sensitivity,"; "),  
                        precision = join(precision,"; "))
                    
                    append!(res, res_i)
                
            end
        end
        CSV.write("res_$similarity.csv", res; bufsize = 54359802)

    end
    
        Plots.plot(thr, jaccard_index, label = "jaccard_index", marker = :o)
        Plots.plot!(thr, speceficity, label = "specificity", marker = :o)
        Plots.plot!(thr, sensitivity, label = "sensitivity", marker = :o)
        Plots.plot!(thr, precision, label = "precision", marker = :o)
        Plots.plot!(legend = :outerright)

""" """

""" batch calculations distances_simulation (with no network saved)"""


    data_set =  CSV.read("Database_INTERNAL_2022-11-17_filtered.csv",DataFrame);
    mz_networks_full = CSV.read("mz_networks.csv", DataFrame);
    mz_networks = mz_networks_full[:,:];

    similarity_list = collect(keys(math_funcs))

    thr = round.(collect(0:0.1:1) .-1,digits = 2)
    similarity = similarity_list[22]

    for similarity in similarity_list[26:28]
        @time begin
            println(similarity)
            res = DataFrame()
            #for i in ProgressBar(1:10)
            for i in ProgressBar(1:size(mz_networks,1))

                    locs_spectra_i = parse.(Int64, vcat(split.(mz_networks.locs[i], "; ")...))
                    mz_networks_i = data_set[locs_spectra_i,:]
                    mz_spectra_i = getVec.(mz_networks_i.MZ_VALUES)
                    int_spectra_i = getVec.(mz_networks_i.MZ_INT_REL)
                    inchikeys = mz_networks_i.INCHIKEY

                #1. first calculate similarity between all spectra
                    similarity_matrix = similarity_matrix_calc(mz_spectra_i, int_spectra_i, similarity, not_distance = false)

                #2. iterate mol networks through thresholds
                    conf_mx_res = thr_iters_star_network_simulation(thr, (similarity_matrix)*-1, inchikeys)

                #5 calculate jaccard index, Sensitivity and Specificity
                
                    jaccard_index, specificity, sensitivity,  precision  = calculate_metrics(conf_mx_res)

                #6 write results

                    res_i = DataFrame( 
                        TP = join(conf_mx_res["TP"], "; "),
                        FP = join(conf_mx_res["FP"], "; "),
                        TN = join(conf_mx_res["TN"], "; "),
                        FN = join(conf_mx_res["FN"], "; "),
                        JI = join(jaccard_index,"; "), 
                        specificity = join(specificity,"; "), 
                        sensitivity = join(sensitivity,"; "),  
                        precision = join(precision,"; "))
                    
                    append!(res, res_i)
            end
        end

        CSV.write("res_$similarity.csv", res; bufsize = 54359802)

    end
    
        Plots.plot(thr, jaccard_index, label = "jaccard_index", marker = :o)
        Plots.plot!(thr, speceficity, label = "specificity", marker = :o)
        Plots.plot!(thr, sensitivity, label = "sensitivity", marker = :o)
        Plots.plot!(thr, precision, label = "precision", marker = :o)
        Plots.plot!(legend = :outerright)


""" """


""" comparison between different metrics"""

    #calculate median value
    x = collect(0:0.1:1.0)
    #path2files = "5_results/tol_mz001_allsims/res_mol_network"
    path2files = "5_results/tol_mz001_positive"
    res_files = readdir(path2files )[contains.(readdir(path2files ),".csv")]
 
       
    """ calculate additional metrics"""
        path2files_old = "5_results/res_restore"
        res_filesold = readdir(path2files_old  )[contains.(readdir(path2files_old  ),".csv")]
        readdir("5_results/res_2810")
        
        res_i = CSV.read(joinpath( path2files_old ,res_filesold[i]),DataFrame)

    
            tp_i = [parse.(Float64, split(x, "; ")) for x in res_i.TP]
            tp_matrix = hcat(tp_i...)'

            fp_i =  [parse.(Float64, split(x, "; ")) for x in res_i.FP]
            fp_matrix = hcat(fp_i...)'

            tn_i = [parse.(Float64, split(x, "; ")) for x in res_i.TN]
            tp_matrix = hcat(tn_i...)'

            fn_i =  [parse.(Float64, split(x, "; ")) for x in res_i.FN]
            fp_matrix = hcat(fn_i...)'

            specificity= [parse.(Float64, split(x, "; ")) for x in res_i.specificity]
            specificity_matrix = hcat(specificity...)'

            sensitivity_i = [parse.(Float64, split(x, "; ")) for x in res_i.sensitivity]
            sensitivity_matrix = hcat(sensitivity_i...)'

            JaccardScore_vector = tp_i[1]./(tp_i[1].+fp_i[1] .+ fn_i[1])
            JaccardScore_matrix = tp_matrix./(tp_matrix.+fp_matrix .+ fn_matrix)

            path2filesold = "5_results/tol_mz001_positive - Copy"
            res_i_old = CSV.read(joinpath(path2filesold ,res_files[i]),DataFrame)

       

        #

        for i in ProgressBar(1: length(res_filesold))
            
            res_i = CSV.read(joinpath(path2files_old,res_filesold[i]),DataFrame)

            tp_i = [parse.(Float64, split(x, "; ")) for x in res_i.TP]
            tp_matrix = hcat(tp_i...)'

            fp_i =  [parse.(Float64, split(x, "; ")) for x in res_i.FP]
            fp_matrix = hcat(fp_i...)'

            tn_i = [parse.(Float64, split(x, "; ")) for x in res_i.TN]
            tn_matrix = hcat(tn_i...)'

            fn_i = [parse.(Float64, split(x, "; ")) for x in res_i.FN]
            fn_matrix = hcat(fn_i...)'

            specificity= [parse.(Float64, split(x, "; ")) for x in res_i.specificity]
            specificity_matrix = hcat(specificity...)'

            sensitivity_i = [parse.(Float64, split(x, "; ")) for x in res_i.sensitivity]
            sensitivity_matrix = hcat(sensitivity_i...)'

            fn_matrix = (tp_matrix.* (sensitivity_matrix.*(-1) .+1)) ./ sensitivity_matrix
            fn_matrix .= ifelse.(isnan.(fn_matrix ), 0.0, fn_matrix )
            fn_matrix  .= round.(fn_matrix , digits = 0)
            fn_vecs = [fn_matrix[i, :] for i in 1:size(fn_matrix, 1)]

            tn_matrix = (fp_matrix.* specificity_matrix)./ (specificity_matrix .*(-1) .+1)
            tn_matrix .= ifelse.(isnan.(tn_matrix), 0.0, tn_matrix)
            tn_matrix .= round.(tn_matrix, digits = 0)
            tn_vecs = [tn_matrix[i, :] for i in 1:size(tn_matrix, 1)]
            

            #sensitivity
            sensitivity = tp_matrix./(tp_matrix.+ fn_matrix)
            sensitivity .= ifelse.(isnan.(sensitivity), 0.0, sensitivity)
            sensitivity .= round.(sensitivity, digits = 2)
            sensitivity_vecs = [sensitivity[i, :] for i in 1:size(sensitivity, 1)]

            #specificity
            specificity = tn_matrix./(tn_matrix.+ fp_matrix)
            specificity .= ifelse.(isnan.(specificity), 0.0, specificity)
            specificity .= round.(specificity, digits = 2)
            specificity_vecs = [specificity[i, :] for i in 1:size(specificity, 1)]

            # #jaccard score
            JaccardScore = tp_matrix./(tp_matrix.+fp_matrix .+ fn_matrix)
            JaccardScore .= ifelse.(isnan.(JaccardScore), 0.0, JaccardScore)
            JaccardScore .= round.(JaccardScore, digits = 2)
            JaccardScore_vecs = [JaccardScore[i, :] for i in 1:size(JaccardScore, 1)]

            #precision
            precision = tp_matrix./(tp_matrix.+fp_matrix)
            precision .= ifelse.(isnan.(precision), 0.0, precision)
            precision .= round.(precision, digits = 2)
            precision_vecs = [precision[i, :] for i in 1:size(precision, 1)]
             
            #fdr = fp/(tp+Fp) NaN  = all positives -> disconnected
            FDR = fp_matrix./(tp_matrix.+fp_matrix)
            FDR .= ifelse.(isnan.(FDR), 0.0, FDR)
            FDR .= round.(FDR, digits = 3)
            FDR_vecs = [FDR[i, :] for i in 1:size(FDR, 1)]

            #accuracy = (TP+Tn)/(TP+Tn+FP+FN) NaN all disconnected therefore accuracy 0
            accuracy = (tp_matrix .+ tn_matrix)./(tp_matrix.+ tn_matrix .+ fp_matrix .+ fn_matrix)
            accuracy .= ifelse.(isnan.(accuracy), 0.0, accuracy)
            accuracy .= round.(accuracy, digits = 3)
            accuracy_vecs = [accuracy[i, :] for i in 1:size(accuracy, 1)]

            #f1-score  = 2*TP/(2*TP+FP+FN)
            f1score = (tp_matrix .*2)./(tp_matrix.*2 .+ fp_matrix .+ fn_matrix)
            f1score .= ifelse.(isnan.(f1score), 0.0, f1score)
            f1score .= round.(f1score, digits = 3)
            f1score_vecs = [f1score[i, :] for i in 1:size(f1score, 1)]

            try
                res_new = [res_i[:,1:3] DataFrame(TN = join.(tn_vecs,"; "), FN = join.(fn_vecs,"; ")) res_i[:,6:9] DataFrame(FDR = join.(FDR_vecs,"; "), Accuracy = join.(accuracy_vecs,"; "), F1score = join.(f1score_vecs,"; "))]

                CSV.write(joinpath("5_results/res_2810", res_filesold[i]),res_new)
            catch
                res_new = [res_i[:,1:2] DataFrame(TN = join.(tn_vecs,"; "), FN = join.(fn_vecs,"; ")) res_i[:,5:8] DataFrame(FDR = join.(FDR_vecs,"; "), Accuracy = join.(accuracy_vecs,"; "), F1score = join.(f1score_vecs,"; "))]
                CSV.write(joinpath("5_results/res_2810", res_filesold[i]),res_new)
            end

        end
    """ """
    path2files = "5_results/res_2810"
    res_files = readdir("5_results/res_2810")

    JI_meadian_all = Matrix{Float64}(undef, length(res_files), length(x))
    specificity_all = Matrix{Float64}(undef, length(res_files), length(x))
    sensitivity_all = Matrix{Float64}(undef, length(res_files), length(x))
    precision_all = Matrix{Float64}(undef, length(res_files), length(x))
    FDR_all = Matrix{Float64}(undef, length(res_files), length(x))
    accuracy_all = Matrix{Float64}(undef, length(res_files), length(x))
    f1score_all = Matrix{Float64}(undef, length(res_files), length(x))
      

    for i in ProgressBar(1:length(res_files))

        res_i = CSV.read(joinpath(path2files ,res_files[i]),DataFrame)

        JI_mat_i, specificity_mat_i, sensitivity_mat_i, precision_mat_i, accuracy_mat_i, f1score_mat_i, FDR_mat_i = res_vec2mat(res_i)

        JI_meadian_all[i,:]= median(JI_mat_i, dims=1)
        specificity_all[i,:]= median(specificity_mat_i, dims=1)
        sensitivity_all[i,:]= median(sensitivity_mat_i, dims=1)
        precision_all[i,:]= median(precision_mat_i, dims=1)
        FDR_all[i,:]= median(FDR_mat_i, dims=1)
        accuracy_all[i,:]= median(accuracy_mat_i, dims=1)
        f1score_all[i,:]= median(accuracy_mat_i, dims=1)


    end

    TP_median_all = Matrix{Float64}(undef, length(res_files), length(x))
    FP_median_all = Matrix{Float64}(undef, length(res_files), length(x))
    TN_median_all = Matrix{Float64}(undef, length(res_files), length(x))
    FN_median_all = Matrix{Float64}(undef, length(res_files), length(x))

    #for i in ProgressBar(1:length(res_files))
    for i in ProgressBar(1:length(res_files))

        res_i = CSV.read(joinpath(path2files ,res_files[i]),DataFrame)
        TP_mat_i = reduce(hcat, [parse.(Float64, vcat(split.(res_i.TP[i], "; ")...)) for i in eachindex(res_i.TP)])'
        FP_mat_i = reduce(hcat, [parse.(Float64, vcat(split.(res_i.FP[i], "; ")...)) for i in eachindex(res_i.FP)])'
        TN_mat_i = reduce(hcat, [parse.(Float64, vcat(split.(res_i.TN[i], "; ")...)) for i in eachindex(res_i.TN)])'
        FN_mat_i = reduce(hcat, [parse.(Float64, vcat(split.(res_i.FN[i], "; ")...)) for i in eachindex(res_i.FN)])'

        #JI_mat_i, specificity_mat_i, sensitivity_mat_i, precision_mat_i, accuracy_mat_i, f1score_mat_i, FDR_mat_i = res_vec2mat(res_i)

        TP_median_all[i,:]= median(TP_mat_i, dims=1)
        FP_median_all[i,:]= median(FP_mat_i, dims=1)
        TN_median_all[i,:]= median(TN_mat_i, dims=1)
        FN_median_all[i,:]= median(FN_mat_i, dims=1)

    end

    sim_label_initial = split.(res_files,"res_")
    sim_label =[sim_label_initial[i][2][1:end-4] for i in eachindex(sim_label_initial)]
    # Dendrogram + heatmap
    xlabelplot  = collect(0.0:0.1:1)
    # Rotate: Transpose matrix and flip hcorder

    sp.heatmap(TP_median_all,
                yflip=true,
                yticks=(1:length(res_files), sim_label),
                xticks = (1:length(xlabelplot), string.(xlabelplot)),
                color=:viridis,
                xlabel = "Threshold",
                title = "TP", dpi = 600)
                sp.savefig("5_results/figures_tol_mz001_positive/tp.png")

    sp.heatmap(FP_median_all,
                yflip=true,
                yticks=(1:length(res_files), sim_label),
                xticks = (1:length(xlabelplot), string.(xlabelplot)),
                color=:viridis,
                xlabel = "Threshold",
                title = "FP", dpi = 600)
                sp.savefig("5_results/figures_tol_mz001_positive/fp.png")

    sp.heatmap(TN_median_all,
                yflip=true,
                yticks=(1:length(res_files), sim_label),
                xticks = (1:length(xlabelplot), string.(xlabelplot)),
                color=:viridis,
                xlabel = "Threshold",
                title = "TN", dpi = 600)
                sp.savefig("5_results/figures_tol_mz001_positive/tn.png")



    """ plotting """

        # Make a DataFrame with one column per threshold
        df_JI = DataFrame(JI_meadian_all, :auto)
        # Rename columns to match threshold values in x
        rename!(df_JI, string.(x))
        res_all_JI = [DataFrame(Similarity = res_files) df_JI]
        
        CSV.write(joinpath("5_results","res_all_JI_2810.csv"),res_all_JI)

        my_palette = palette(:tab20,24)              # distinct colors
        linestyles = [:solid, :dash, :dot, :dashdot]

        p = sp.plot(title="Jaccard Score median", size=(1000,670), dpi=600,     
        xtickfont = font(14),            # X-axis tick font size
        ytickfont = font(14),            # Y-axis tick font size
        guidefont = font(16),            # Axis label font size
        titlefont = font(18),
        xlabel = "Threshold value" )           # Title font size)

            for i in 1:size(JI_meadian_all, 1)
                display(
                    sp.plot!(x, JI_meadian_all[i, :], 
                    label = "$(res_all_JI.Similarity[i][5:end-4])",
                    legendfont = font(12),
                    legend = :outertopright,
                    color = my_palette[mod1(i, length(my_palette))],
                    linestyle = linestyles[mod1(i, length(linestyles))],
                    lw = 2.5))
            end
        sp.savefig(joinpath("5_results/figures_tol_mz001_positive", "JI_median2810.png"))

        
        # Make a DataFrame with one column per threshold
        df_specificity_ = DataFrame(specificity_all, :auto)
        # Rename columns to match threshold values in x
        rename!(df_specificity_, string.(x))
        res_all_specificity_ = [DataFrame(Similarity = res_files) df_specificity_]
        CSV.write(joinpath("5_results","res_all_specificity_2810.csv"),res_all_specificity_)

        my_palette = palette(:tab20,23)               # distinct colors
        linestyles = [:solid, :dash, :dot, :dashdot]

        p = sp.plot(title="Specificity median", size=(1000,670), dpi=600,
        xtickfont = font(14),            # X-axis tick font size
        ytickfont = font(14),            # Y-axis tick font size
        guidefont = font(16),            # Axis label font size
        titlefont = font(18),
        xlabel = "Threshold value" )

            for i in 1:size(specificity_all, 1)
                display(
                    sp.plot!(x, specificity_all[i, :], 
                    label = "$(res_all_specificity_.Similarity[i][5:end-4])",
                    legendfont = font(12),
                    legend = :outertopright,
                    color = my_palette[mod1(i, length(my_palette))],
                    linestyle = linestyles[mod1(i, length(linestyles))],
                    lw = 2.5))
            end
        sp.savefig(joinpath("5_results/figures_tol_mz001_positive", "Specificity_median2810.png"))


        # Make a DataFrame with one column per threshold
        df_sensitivity = DataFrame(sensitivity_all, :auto)
        # Rename columns to match threshold values in x
        rename!(df_sensitivity, string.(x))
        res_all_sensitivity = [DataFrame(Similarity = res_files) df_sensitivity]
        CSV.write(joinpath("5_results","res_all_sensitivity2810.csv"), res_all_sensitivity)


        my_palette = palette(:tab20,23)               # distinct colors
        linestyles = [:solid, :dash, :dot, :dashdot]

        p = sp.plot(title="Sensitivity median", size=(1000,670), dpi=600,
        xtickfont = font(14),            # X-axis tick font size
        ytickfont = font(14),            # Y-axis tick font size
        guidefont = font(16),            # Axis label font size
        titlefont = font(18),
        xlabel = "Threshold value" )

            for i in 1:size(sensitivity_all, 1)
                display(
                    sp.plot!(x, sensitivity_all[i, :], 
                    label = "$(res_all_sensitivity.Similarity[i][5:end-4])",
                    legendfont = font(12),
                    legend = :outertopright,
                    color = my_palette[mod1(i, length(my_palette))],
                    linestyle = linestyles[mod1(i, length(linestyles))],
                    lw = 2.5))
            end
        sp.savefig(joinpath("5_results/figures_tol_mz001_positive", "Sensitivity_median2810.png"))

        # Make a DataFrame with one column per threshold
        df_precision = DataFrame(precision_all, :auto)
        # Rename columns to match threshold values in x
        rename!(df_precision, string.(x))
        res_all_precision = [DataFrame(Similarity = res_files) df_precision]
        CSV.write(joinpath("5_results","res_all_precision2810.csv"), res_all_precision)


        my_palette = palette(:tab20,23)               # distinct colors
        linestyles = [:solid, :dash, :dot, :dashdot]

        p = sp.plot(title="Precision median", size=(1000,670), dpi=600,
        xtickfont = font(14),            # X-axis tick font size
        ytickfont = font(14),            # Y-axis tick font size
        guidefont = font(16),            # Axis label font size
        titlefont = font(18),
        xlabel = "Threshold value" )

            for i in 1:size(precision_all, 1)
                display(
                    sp.plot!(x, precision_all[i, :], 
                    label = "$(res_all_precision.Similarity[i][5:end-4])",
                    legendfont = font(12),
                    legend = :outertopright,
                    color = my_palette[mod1(i, length(my_palette))],
                    linestyle = linestyles[mod1(i, length(linestyles))],
                    lw = 2.5))
            end
        sp.savefig(joinpath("5_results/figures_tol_mz001_positive", "Precision_median2810.png"))

        
        # Make a DataFrame with one column per threshold
        df_accuracy = DataFrame(accuracy_all, :auto)
        # Rename columns to match threshold values in x
        rename!(df_accuracy, string.(x))
        res_all_accuracy = [DataFrame(Similarity = res_files) df_accuracy]
        CSV.write(joinpath("5_results","res_all_accuracy2810.csv"), res_all_accuracy)


        my_palette = palette(:tab20,23)               # distinct colors
        linestyles = [:solid, :dash, :dot, :dashdot]

        p = sp.plot(title="Accuaracy median", size=(1000,670), dpi=600,
        xtickfont = font(14),            # X-axis tick font size
        ytickfont = font(14),            # Y-axis tick font size
        guidefont = font(16),            # Axis label font size
        titlefont = font(18),
        xlabel = "Threshold value" )

            for i in 1:size(accuracy_all, 1)
                display(
                    sp.plot!(x, accuracy_all[i, :], 
                    label = "$(res_all_sensitivity.Similarity[i][5:end-4])",
                    legendfont = font(12),
                    legend = :outertopright,
                    color = my_palette[mod1(i, length(my_palette))],
                    linestyle = linestyles[mod1(i, length(linestyles))],
                    lw = 2.5))
            end
        sp.savefig(joinpath("5_results/figures_tol_mz001_positive", "Accuracy_median2810.png"))


                
        # Make a DataFrame with one column per threshold
        df_FDR = DataFrame(FDR_all, :auto)
        # Rename columns to match threshold values in x
        rename!(df_FDR, string.(x))
        res_all_FDR = [DataFrame(Similarity = res_files) df_FDR]
        CSV.write(joinpath("5_results","res_all_FDR2810.csv"), res_all_FDR)


        my_palette = palette(:tab20,23)               # distinct colors
        linestyles = [:solid, :dash, :dot, :dashdot]

        p = sp.plot(title="FDR median", size=(1000,670), dpi=600,
        xtickfont = font(14),            # X-axis tick font size
        ytickfont = font(14),            # Y-axis tick font size
        guidefont = font(16),            # Axis label font size
        titlefont = font(18),
        xlabel = "Threshold value" )

            for i in 1:size(FDR_all, 1)
                display(
                    sp.plot!(x, FDR_all[i, :], 
                    label = "$(res_all_FDR.Similarity[i][5:end-4])",
                    legendfont = font(12),
                    legend = :outertopright,
                    color = my_palette[mod1(i, length(my_palette))],
                    linestyle = linestyles[mod1(i, length(linestyles))],
                    lw = 2.5))
            end
        sp.savefig(joinpath("5_results/figures_tol_mz001_positive", "FDR_median2810.png"))

        # Make a DataFrame with one column per threshold
        df_f1score = DataFrame(f1score_all, :auto)
        # Rename columns to match threshold values in x
        rename!(df_f1score, string.(x))
        res_all_f1score = [DataFrame(Similarity = res_files) df_f1score]


        CSV.write(joinpath("5_results","res_all_f1score2810.csv"), res_all_f1score)


        my_palette = palette(:tab20,23)               # distinct colors
        linestyles = [:solid, :dash, :dot, :dashdot]

        p = sp.plot(title="F1-score median", size=(1000,670), dpi=600,
        xtickfont = font(14),            # X-axis tick font size
        ytickfont = font(14),            # Y-axis tick font size
        guidefont = font(16),            # Axis label font size
        titlefont = font(18),
        xlabel = "Threshold value" )

            for i in 1:size(f1score_all, 1)
                display(
                    sp.plot!(x, f1score_all[i, :], 
                    label = "$(res_all_f1score.Similarity[i][5:end-4])",
                    legendfont = font(12),
                    legend = :outertopright,
                    color = my_palette[mod1(i, length(my_palette))],
                    linestyle = linestyles[mod1(i, length(linestyles))],
                    lw = 2.5))
            end
        sp.savefig(joinpath("5_results/figures_tol_mz001_positive", "f1score_median2810.png"))

""" set thr 0.7 """

    #path2files = "5_results/tol_mz001_allsims/res_mol_network"
    path2files = "5_results/res_2810"
    res_files = readdir(path2files )[contains.(readdir(path2files ),".csv")]

    JI_07 = DataFrame()
    specificity_07 = DataFrame()
    sensitivity_07 = DataFrame()
    precision_07 = DataFrame()
    accuracy_07 = DataFrame()
    f1score_07 = DataFrame()
    FDR_07 = DataFrame()



    for i in ProgressBar(1:length(res_files))
        res_i = CSV.read(joinpath(path2files ,res_files[i]),DataFrame)

        JI_mat, specificity_mat, sensitivity_mat, precision_mat, accuracy_mat, f1score_mat, FDR_mat = res_vec2mat(res_i)

        col_name = res_files[i][5:end-4]
        JI_07  = [JI_07 DataFrame("$col_name" => JI_mat[:,8])]
        specificity_07 = [specificity_07 DataFrame("$col_name" => specificity_mat[:,8])]
        sensitivity_07 = [sensitivity_07 DataFrame("$col_name" => sensitivity_mat[:,8])]
        precision_07 = [precision_07 DataFrame("$col_name" => precision_mat[:,8])]
        accuracy_07 = [accuracy_07 DataFrame("$col_name" => accuracy_mat[:,8])]
        f1score_07 = [f1score_07 DataFrame("$col_name" => f1score_mat[:,8])]
        FDR_07 = [FDR_07 DataFrame("$col_name" =>FDR_mat[:,8])]


    end


    sp.boxplot(
        Matrix(specificity_07  ),
        orientation = :horizontal,
        yticks = (1:ncol(specificity_07  ), names(specificity_07  )),
        label = "", title = "Specificity, threshold = 0.7", dpi = 600
    )

    sp.savefig(joinpath("5_results/figures_tol_mz001_positive", "Specificityy_ver2.png"))

    median(Matrix(JI_07))
    maximum(Matrix(JI_07))
    minimum(Matrix(JI_07))
    count(Matrix(JI_07).==1)

    no_helliner = names(JI_07)[names(JI_07).!="hellinger"]
    no_clark = no_helliner[no_helliner.!="clark"]

    median(Matrix(JI_07[:,no_clark]))

    locs_only1 = []
    for i  = 1: size(JI_07,1)
        n = count(Vector(JI_07[i,:]).==1.0)
        if n==1
            push!(locs_only1,i)
        end
    end
    
    length(locs)/size(JI_07,1)
    JI_07[locs,:]
    JI_07[locs_only1,:]

    x = []
    for i in 1:size(JI_07,2)
        n = count(Vector(JI_07[locs_only1,i]).==1.0)
        push!(x, n)

    end
    

    # assuming you have:
    # x :: Vector{<:Number}  (length 23)
    # names(JI_07) :: Vector{String} (also length 23)
    x = vec(median(Matrix(JI_07),dims = 1))

    sorted_idx = sortperm(x; rev=true)
    x_sorted = x[sorted_idx]
    labels_sorted = names(JI_07)[sorted_idx]



    sp.bar(
        labels_sorted,
        x_sorted;
        legend = false,
        xlabel = "Similarity measure",
        ylabel = "Count (resolved networks)",
        permute=(:x,:y),
        xticks = (1:length(labels_sorted), labels_sorted),      # smaller font for long names
        left_margin = 15mm,       # extra space for labels
        fillcolor = :steelblue,
        linecolor = :black,
        ylims = (0, maximum(x_sorted) + 0.2),
        dpi = 600
    )

    sp.savefig(joinpath("5_results/figures_tol_mz001_positive", "median_07.png"))

""" investigating behavior"""

    path2res = "C:/Users/vturkin/OneDrive - UvA/Desktop/Project/2_Projects/13_mol_network/5_results/tol_mz001_allsims"
    readdir(path2res)

    JI_meadian_all
    res_all_JI.Similarity[maximum(JI_meadian_all, dims = 2)[:,1].<0.2]


""" hierarchical clustering of different metrics"""

    res_all = [JI_meadian_all specificity_all sensitivity_all precision_all]

    D = pairwise(Euclidean(), res_all')    # distance matrix

    hc = hclust(D, linkage=:average)
    hcorder = hc.order

    sim_label_initial = split.(res_all_JI.Similarity[hcorder],"res_")
    sim_label =[sim_label_initial[i][2][1:end-4] for i in eachindex(sim_label_initial)][hcorder]
    # Dendrogram + heatmap
    xlabelplot  = [collect(0.0:0.1:1); collect(0.0:0.1:1); collect(0.0:0.1:1); collect(0.0:0.1:1)]
    # Rotate: Transpose matrix and flip hcorder
    mat_rotated = res_all[hcorder,:]
    p = sp.heatmap(mat_rotated, 
                yflip=true, 
                colorbar=false,
                yticks=(1:length(res_all_JI.Similarity), sim_label),
                xticks = (1:length(xlabelplot), string.(xlabelplot)),
                xrotation=90,
                color=:viridis,
                xlabel = "Jaccard Index      Specificity    Sensitivity    Precision")
        display(p)

    # Standalone colorbar
    p_cbar = sp.heatmap([0 1],  color=:viridis, colorbar=true,
                 axis=nothing, framestyle=:none, lims=(-1,0), clim=(minimum(mat_rotated), maximum(mat_rotated)))
    # Plot: left dendrogram + right heatmap
    sp.plot(p,
        sp.plot(hc, yticks=false, xticks=false,xlims = (0,5), orientation=:horizontal, framestyle = :none), 
        p_cbar,
        layout=sp.grid(1, 3, widths=[0.65, 0.3, 0.05]),
        size=(1100, 800), dpi = 600
    )  


    sp.savefig("5_results/figures_tol_mz001_positive/all_res_clustering_dendogram_3010.png")



    res_all = [JI_meadian_all specificity_all sensitivity_all precision_all]

    D = pairwise(Euclidean(), res_all')    # distance matrix
    hc = hclust(D, linkage=:ward)
    hcorder = hc.order

    # cluster assignments
    clusters = cutree(hc, k=5)

    # Reorder labels
    sim_label_initial = split.(res_all_JI.Similarity[hcorder],"res_")
    sim_label = [sim_label_initial[i][2][1:end-4] for i in eachindex(sim_label_initial)][hcorder]

    xlabelplot  = [collect(0.0:0.1:1); collect(0.0:0.1:1); collect(0.0:0.1:1); collect(0.0:0.1:1)]
    mat_rotated = res_all[hcorder,:]

    # assign a color to each cluster
    my_palette = distinguishable_colors(20)
    row_colors = [my_palette[c] for c in clusters[hcorder]]

    # Heatmap with colored tick labels
    p = sp.heatmap(mat_rotated, 
                yflip=true,
                colorbar=false,
                yticks=(1:length(sim_label), sim_label),
                xticks=(1:length(xlabelplot), string.(xlabelplot)),
                xrotation=90,
                color=:viridis,
                xlabel="Jaccard Index   Specificity   Sensitivity   Precision",
                foreground_color_axis =row_colors)  # <-- cluster coloring of labels

    # Dendrogram
    p_dendro = sp.plot(hc, 
                    yticks=false, xticks=false, xlims=(0,5),
                    orientation=:horizontal, framestyle=:none)

    # Colorbar
    p_cbar = sp.heatmap([0 1], color=:viridis, colorbar=true,
                    axis=nothing, framestyle=:none, lims=(-1,0),
                    clim=(minimum(mat_rotated), maximum(mat_rotated)))

    # Combine
    sp.plot(p, p_dendro, p_cbar,
        layout=sp.grid(1, 3, widths=[0.65, 0.3, 0.05]),
        size=(1200, 800), dpi=600)

        

 # for `findfirst`





""" pca """
    
    pca = pyimport("sklearn.decomposition")
    pca_test = pca.PCA(n_components=3)

    res_all = [JI_meadian_all specificity_all sensitivity_all precision_all]
    #scores
    X_reduced = pca_test.fit_transform(res_all)
    #loadings
    loadings = pca_test.components_
    var = round.(pca_test.explained_variance_ratio_, digits = 2)


    sp.scatter(X_reduced[:,1],X_reduced[:,2],
        marker=:circle, 
        linewidth=0,
        color = my_palette[1:length(labels)],  # assign colors
        label = labels,
        markerstrokewidth=0)

    # color palette and marker shapes
    my_palette = palette(:tab20, 24)
    my_markers = [:circle, :utriangle, :star] 

    # labels from your data
    labels = [s[5:end-4] for s in res_all_JI.Similarity]

    p = sp.scatter(xlabel = "PC1 ($(var[1])%)", ylabel = "PC2 ($(var[2])%)" )
    for i in 1:length(labels)
        p = sp.scatter!(
            p,
            [X_reduced[i, 1]], [X_reduced[i, 2]],   # single point per series
            color = my_palette[mod1(i, length(my_palette))],
            marker = my_markers[mod1(i, length(my_markers))],
            label = labels[i],
            legend = :outertopright,
            markerstrokewidth=0,
            markersize = 8,
            dpi = 600
        )
    end
    
    display(p)
    sp.savefig("5_results/figures_tol_mz001_positive/all_res_pca_3010.png")


     clusters = cutree(hc, k=4)
 
    sp.scatter(X_reduced[ clusters.==1 ,1],X_reduced[ clusters.==1,2],
        marker=:circle, 
        linewidth=0,
        color = my_palette[1],  # assign colors
        label = false,
        markerstrokewidth=0)

     
    sp.scatter!(X_reduced[ clusters.==2 ,1],X_reduced[ clusters.==2,2],
        marker=:circle, 
        linewidth=0,
        color = my_palette[6],  # assign colors
        label = false,
        markerstrokewidth=0)

    sp.scatter!(X_reduced[ clusters.==3 ,1],X_reduced[ clusters.==3,2],
        marker=:circle, 
        linewidth=0,
        color = my_palette[3],  # assign colors
        label = false,
        markerstrokewidth=0)

    sp.scatter!(X_reduced[ clusters.==4 ,1],X_reduced[ clusters.==4,2],
        marker=:circle, 
        linewidth=0,
        color = my_palette[4],  # assign colors
        label = false,
        markerstrokewidth=0)

    sp.scatter!(X_reduced[ clusters.==5 ,1],X_reduced[ clusters.==5,2],
        marker=:circle, 
        linewidth=0,
        color = my_palette[5],  # assign colors
        label = false,
        markerstrokewidth=0)


""" end  """
""" optimal threshold statistics"""

    path2res = "5_results/res_2810"
    res_files = readdir(path2res)

    opt_thr_df = DataFrame()
    opt_JI_df = DataFrame()
    for i in eachindex(1: length(res_files))
        res_i = CSV.read(joinpath(path2res, res_files[i]), DataFrame)
        JI_mat_i = reduce(hcat, [parse.(Float64, vcat(split.(res_i.JI[i], "; ")...)) for i in eachindex(res_i.JI)])'
        locs = argmax(JI_mat_i, dims = 2)

        dummy_array = repeat(collect(0:0.1:1)', length(locs), 1)
        
        opt_thr_i =  vec(dummy_array[locs])
        name_i = res_files[i][5:end-4]

        opt_thr_df =[opt_thr_df DataFrame(name_i => opt_thr_i)]
        opt_JI_df =[opt_JI_df DataFrame(name_i => vec(JI_mat_i[locs]))]

    end
    


    my_palette = palette(:tab20,24)              # distinct colors
    linestyles = [:solid, :dash, :dot, :dashdot]


    p = sp.density(opt_thr_df[:, 1], 
        bandwidth=0.1, 
        label=names(opt_thr_df)[1], 
        legend=:outerright,
        xlabel = "Threshold value",
        ylabel = "Density of occurence",
        dpi = 600)
    for i in 2:20
        sp.density!(opt_thr_df[:, i],
        bandwidth=0.1,
        label=names(opt_thr_df)[i],
        color = my_palette[mod1(i, length(my_palette))],
        linestyle = linestyles[mod1(i, length(linestyles))],
        lw = 1.5)
    end
    display(p)

    sp.savefig(joinpath("5_results/figures_tol_mz001_positive", "opt_threshold_density.png"))

    median(Matrix(opt_thr_df), dims = 1)
    sort(vec(median(Matrix(opt_thr_df), dims = 1)))

    median(Matrix(opt_JI_df), dims = 1)
    sort(vec(median(Matrix(opt_JI_df), dims = 1)))

    p = sp.density(opt_JI_df[:, 1], 
        bandwidth=0.1, 
        label=names(opt_JI_df)[1], 
        legend=:outerright,
        xlabel = "Threshold value",
        ylabel = "Density of occurence",
        dpi = 600)
    for i in 2:20
        sp.density!(opt_JI_df[:, i],
        bandwidth=0.1,
        label=names(opt_JI_df)[i],
        color = my_palette[mod1(i, length(my_palette))],
        linestyle = linestyles[mod1(i, length(linestyles))],
        lw = 1.5)
    end
    display(p)

    sp.savefig(joinpath("5_results/figures_tol_mz001_positive", "maximumJI_density.png"))

    n_pure_clusters = []
    for i = 1:size(opt_JI_df,1)
        if count(Vector(opt_JI_df[i,:]).==1.0)==1
            push!(n_pure_clusters,i)
        end
    end

    opt_thr_df[n_pure_clusters,:]


    
    x = []
    for i in 1:size(opt_JI_df,2)
        n = count(Vector(opt_JI_df[n_pure_clusters,i]).==1.0)
        push!(x, n)

    end
    

    # assuming you have:
    # x :: Vector{<:Number}  (length 23)
    # names(opt_JI_df) :: Vector{String} (also length 23)
    #x = vec(median(Matrix(opt_JI_df),dims = 1))

    sorted_idx = sortperm(x; rev=true)
    x_sorted = x[sorted_idx]
    labels_sorted = names(opt_JI_df)[sorted_idx]



    sp.bar(
        labels_sorted,
        x_sorted;
        legend = false,
        xlabel = "Similarity measure",
        ylabel = "Count (resolved networks)",
        permute=(:x,:y),
        xticks = (1:length(labels_sorted), labels_sorted),      # smaller font for long names
        left_margin = 15mm,       # extra space for labels
        fillcolor = :steelblue,
        linecolor = :black,
        ylims = (0, maximum(x_sorted) + 0.2),
        dpi = 600
    )

    sp.savefig(joinpath("5_results/figures_tol_mz001_positive", "median_optimal_thr.png"))
    CSV.write("5_results/res_optimal_thr.csv",opt_thr_df )
    CSV.write("5_results/res_optimal_thr_JI.csv",opt_JI_df )


""" resolved and not resolved vs the number of frags optimal thr"""

    data_set =  CSV.read("Database_positive_filtered.csv",DataFrame)
    mz_networks_full = CSV.read("mz_networks_positive.csv", DataFrame)

    opt_thr_df  = CSV.read("5_results/res_optimal_thr.csv", DataFrame)
    opt_JI_df = CSV.read("5_results/res_optimal_thr_JI.csv", DataFrame)

    all_resolved = []
    for i = 1:size(opt_JI_df,1)
        if count(Vector(opt_JI_df[i,:]).==1.0)>0
            push!(all_resolved,i)
        end
    end

    not_resolved = []
    for i = 1:size(opt_JI_df,1)
        if count(Vector(opt_JI_df[i,:]).==1.0)==0
            push!(not_resolved,i)
        end
    end


    mz_networks_full[all_resolved,:]

    sp.density(mz_networks_full[all_resolved,"n_INCHIKEY"], bandwidth=0.1)


    maximum(mz_networks_full[all_resolved,"n_INCHIKEY"])
    maximum(mz_networks_full[all_resolved,"PRECURSOR_ION"])
    maximum(mz_networks_full[all_resolved,"total_n_SPECTRA"])



    sp.density(mz_networks_full[not_resolved,"n_INCHIKEY"], bandwidth=0.1)

    maximum(mz_networks_full[not_resolved,"n_INCHIKEY"])
    

    limit_to8 = not_resolved[mz_networks_full[not_resolved,"n_INCHIKEY"].<=8]
    maximum(mz_networks_full[limit_to8,"PRECURSOR_ION"])

    limit_to8_and_1125mz = limit_to8[mz_networks_full[limit_to8,"PRECURSOR_ION"].<=1125]
    maximum(mz_networks_full[limit_to8_and_1125mz,"total_n_SPECTRA"])


    limit_to8_and_1125mz_and_201nspect = limit_to8_and_1125mz[mz_networks_full[limit_to8_and_1125mz,"total_n_SPECTRA"].<=201]

    mz_networks_full[limit_to8_and_1125mz_and_201nspect[6],"INCHIKEYS"]

    #record number of frags
    n_frags_vec = []
    n_frags_min = []
    n_frags_max = []
    n_frags_median = []

    
    for i in ProgressBar(1:length(limit_to8_and_1125mz_and_201nspect))
        locs_i = mz_networks_full[limit_to8_and_1125mz_and_201nspect[i],"locs"]
        vec_loc_i = parse.(Int64, split(locs_i,"; "))


        mzvalues_i = getVec.(data_set[vec_loc_i,"MZ_VALUES"])
        n_frags_vec_i = length.(mzvalues_i)

        push!(n_frags_vec, n_frags_vec_i)
        push!(n_frags_min, minimum(n_frags_vec_i))
        push!(n_frags_max, maximum(n_frags_vec_i))
        push!(n_frags_median, median(n_frags_vec_i))
    end


        #record number of frags
    n_frags_res_vec = []
    n_frags_res_min = []
    n_frags_res_max = []
    n_frags_res_median = []

    
    for i in ProgressBar(1:length(all_resolved))
        locs_i = mz_networks_full[all_resolved[i],"locs"]
        vec_loc_i = parse.(Int64, split(locs_i,"; "))


        mzvalues_i = getVec.(data_set[vec_loc_i,"MZ_VALUES"])
        n_frags_vec_i = length.(mzvalues_i)

        push!(n_frags_res_vec, n_frags_vec_i)
        push!(n_frags_res_min, minimum(n_frags_vec_i))
        push!(n_frags_res_max, maximum(n_frags_vec_i))
        push!(n_frags_res_median, median(n_frags_vec_i))
    end



    sp.density(n_frags_res_median, label = "resolved", title = "Density of median value of number of fragments",   linewidth = 3, dpi = 600)
    sp.density!(n_frags_median, label = "not resolved", linewidth = 3, legendfontsize = 12, tickfont  = font(10) )
    xlabel!("Number of fragments in MS2")
    ylabel!("Density")
    sp.savefig("5_results/figures_tol_mz001_positive/nfragsmedian.png")
    
    
    sp.density(n_frags_res_min, label = "resolved", title = "Density of minimum value of number of fragments",   linewidth = 3, dpi = 600)
    sp.density!(n_frags_min, label = "not resolved", linewidth = 3, legendfontsize = 12, tickfont  = font(10) )
    xlabel!("Number of fragments in MS2")
    ylabel!("Density")
    sp.savefig("5_results/figures_tol_mz001_positive/nfragsminimum.png")

    sp.density(n_frags_res_max, label = "resolved", title = "Density of maximum value of number of fragments",   linewidth = 3, dpi = 600)
    sp.density!(n_frags_max, label = "not resolved", linewidth = 3, legendfontsize = 12, tickfont  = font(10) )
    xlabel!("Number of fragments in MS2")
    ylabel!("Density")
     sp.savefig("5_results/figures_tol_mz001_positive/nfragsmaximum.png")



"""optimal threshols database"""
    data_set =  CSV.read("Database_positive_filtered.csv",DataFrame)
    mz_networks_full = CSV.read("mz_networks_positive.csv", DataFrame)


    path2res = "5_results/res_2810"
    res_files = readdir(path2res)

    opt_thr_res = DataFrame()

    for i in eachindex(1: length(res_files))
        res_i = CSV.read(joinpath(path2res, res_files[i]), DataFrame)
        JI_mat_i = reduce(hcat, [parse.(Float64, vcat(split.(res_i.JI[i], "; ")...)) for i in eachindex(res_i.JI)])'
        locs = argmax(JI_mat_i, dims = 2)

        dummy_array = repeat(collect(0:0.1:1)', length(locs), 1)
        
        opt_thr_i =  vec(dummy_array[locs])
        opt_JI_i = vec(JI_mat_i[locs])
        name_i = res_files[i][5:end-4]*" [THR,JI]"

        opt_i = [join([opt_thr_i[i], opt_JI_i[i]],"; ") for i = 1:8290]

        opt_thr_res = [opt_thr_res DataFrame(name_i => opt_i)]


    end


    #find ACCESSION

    accession = []
    for i  = 1:8290
        locs_i = parse.(Int64, split(mz_networks_full.locs[i],"; "))
        push!(accession, join(data_set.ACCESSION[locs_i], "; ")) 
    end

    final_df = [mz_networks_full DataFrame(ACCESSION= accession) opt_thr_res ]

    CSV.write("6_publication/database.csv", final_df)

""" end """