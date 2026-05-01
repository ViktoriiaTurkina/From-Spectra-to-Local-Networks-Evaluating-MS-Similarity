 module SpectralSimilarity

"""
This is a Julia version of previously published Python code: Li Y. et al. Nat Methods. 2021 Dec;18(12):1524-1531, doi: 10.1038/s41592-021-01331-z,  for spectral similarity and distance calculations:

- math_distance functions (many distance formulas)
- ms_distance functions (MS-specific aggregated distances)
- tools utilities: clean_spectrum, match_peaks_in_spectra, match_peaks_with_mz_info_in_spectra, normalize_distance
- spectral_similarity API: distance, all_distance, similarity, all_similarity, multiple_distance, multiple_similarity

Notes / limitations:
- The distance calculations in this code were implemented very closely from original Python implementations, 
  with only minimal changes needed to make them work in Julia.

"""

export distance, all_distance, similarity, all_similarity, multiple_distance, multiple_similarity, math_funcs

using Statistics

const EPS = 1e-12

# -----------------------------
# tools (helpers)
# -----------------------------

"""
clean_spectrum(spectrum; ms2_ppm=nothing, ms2_da=nothing, min_rel_intensity=0.01)

spectrum expected as Nx2 array [mz intensity].
Removes peaks with intensity < min_rel_intensity * max_intensity and normalizes intensities to sum=1.
"""
function clean_spectrum(spectrum::AbstractArray{T,2}; ms2_ppm=nothing, ms2_da=nothing, min_rel_intensity::Float64=0.01) where T <: Real
    if size(spectrum,1) == 0
        return Array{Float32,2}(undef,0,2)
    end
    mz = Float64.(spectrum[:,1])
    inten = Float64.(spectrum[:,2])
    max_i = maximum(inten)
    thr = min_rel_intensity * max_i
    mask = inten .>= thr
    mz_f = mz[mask]
    inten_f = inten[mask]
    if sum(inten_f) > 0
        inten_f .= inten_f ./ sum(inten_f)
    end
    out = hcat(Float32.(mz_f), Float32.(inten_f))
    return out
end

"""
match_peaks_in_spectra(spec_a, spec_b; ms2_ppm=nothing, ms2_da=nothing)

Return an Mx3 Float32 array where columns are [mz_reference, intensity_a, intensity_b].
It constructs a union of peak m/zs from both spectra and assigns intensities if a matching
peak is within tolerance, otherwise 0.
"""

function match_peaks_in_spectra_old(spec_a::AbstractArray{T,2}, spec_b::AbstractArray{T,2}; ms2_ppm=nothing, ms2_da=nothing) where T <: Real
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
    tolerance = ms2_da
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




"""
match_peaks_with_mz_info_in_spectra(spec_query, spec_reference; ms2_ppm=nothing, ms2_da=nothing)

Return Mx4 array: [mz_q, i_q, mz_r, i_r] with zeros where unmatched.
"""
function match_peaks_with_mz_info_in_spectra(spec_q::AbstractArray{T,2}, spec_r::AbstractArray{T,2}; ms2_ppm=nothing, ms2_da=nothing) where T <: Real
    # We'll attempt to match each query peak to the nearest reference peak within tol
    q_mz = Float64.(spec_q[:,1]); q_i = Float64.(spec_q[:,2])
    r_mz = Float64.(spec_r[:,1]); r_i = Float64.(spec_r[:,2])

    if ms2_da === nothing && ms2_ppm === nothing
        error("MS2 tolerance need to be defined!")
    end

    matched = Vector{Tuple{Float64,Float64,Float64,Float64}}()

    used_r = falses(length(r_mz))
    for j in eachindex(q_mz)
        mzq = q_mz[j]
        # find nearest r index within tolerance
        best_k = 0
        best_dist = Inf
        for k in eachindex(r_mz)
            d = abs(mzq - r_mz[k])
            if ms2_da !== nothing
                ok = d <= ms2_da
            else
                ok = d <= mz2ppm(mzq, ms2_ppm)
            end
            if ok && d < best_dist
                best_dist = d; best_k = k
            end
        end
        if best_k != 0
            push!(matched, (mzq, q_i[j], r_mz[best_k], r_i[best_k]))
            used_r[best_k] = true
        else
            push!(matched, (mzq, q_i[j], 0.0, 0.0))
        end
    end
    # also include unmatched r peaks
    for k in eachindex(r_mz)
        if !used_r[k]
            push!(matched, (0.0, 0.0, r_mz[k], r_i[k]))
        end
    end

    M = length(matched)
    out = Array{Float32,2}(undef, M, 4)
    for i in 1:M
        out[i,1] = Float32(matched[i][1])
        out[i,2] = Float32(matched[i][2])
        out[i,3] = Float32(matched[i][3])
        out[i,4] = Float32(matched[i][4])
    end
    return out
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

# ppm helper: absolute tolerance for an m/z at given ppm
mz2ppm(mz::Float64, ppm::Real) = abs(mz) * (ppm/1e6)

"""
normalize_distance(dist, dist_range)
Normalize a raw distance scalar into [0,1] according to the provided dist_range = (min,max).
If max is Inf, a monotonic transform is used to map larger distances closer to 1.
"""
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

# -----------------------------
# math_distance (port)
# -----------------------------

# Tsallis / entropy utilities
function calculate_tsallis_entropy(ints::AbstractVector{T}, entropy_dimension::Real) where T
    s = sum(ints)
    if s > 0
        ints2 = collect(float.(ints))
        ints2 ./= sum(ints2)
        return (sum(x->x^entropy_dimension, ints2) - 1.0) / (1.0 - entropy_dimension)
    else
        return 0.0
    end
end

function tsallis_entropy_distance(p::AbstractVector, q::AbstractVector, entropy_dimension::Real)
    if entropy_dimension <= 0
        error("entropy_dimension must be positive")
    elseif entropy_dimension == 1
        return entropy_distance(p, q) / log(4)
    else
        ent_p = calculate_tsallis_entropy(p, entropy_dimension)
        ent_q = calculate_tsallis_entropy(q, entropy_dimension)
        ent_pq = calculate_tsallis_entropy((p .+ q) ./ 2, entropy_dimension)
        N = sum(2 .*(p./2).^entropy_dimension .+ 2 .*(q./2).^entropy_dimension .- p.^entropy_dimension .- q.^entropy_dimension) / (1.0 - entropy_dimension)
        return (2.0 * ent_pq - ent_p - ent_q) / N
    end
end

using StatsBase: mean

function unweighted_entropy_distance(p::AbstractVector, q::AbstractVector)
    merged = p .+ q
    # approximate entropy using -sum(x*log(x)) with small eps handling
    function ent(x)
        s = sum(x)
        if s <= 0
            return 0.0
        end
        xx = x ./ s
        # avoid log(0)
        return -sum(ifelse.(xx .> 0.0, xx .* log.(xx), 0.0))
    end
    entropy_increase = 2 * ent(merged) - ent(p) - ent(q)
    return entropy_increase
end

function _weight_intensity_by_entropy(x::AbstractVector)
    WEIGHT_START = 0.25
    ENTROPY_CUTOFF = 3.0
    weight_slope = (1 - WEIGHT_START) / ENTROPY_CUTOFF
    s = sum(x)
    out = float.(x)
    if s > 0
        xx = out ./ s
        entropy_x = -sum(ifelse.(xx .> 0.0, xx .* log.(xx), 0.0))
        if entropy_x < ENTROPY_CUTOFF
            weight = WEIGHT_START + weight_slope * entropy_x
            out .= out .^ weight
            tot = sum(out)
            if tot > 0
                out ./= tot
            end
        end
    end
    return out
end

function entropy_distance(p::AbstractVector, q::AbstractVector)
    p_w = _weight_intensity_by_entropy(p)
    q_w = _weight_intensity_by_entropy(q)
    return unweighted_entropy_distance(p_w, q_w)
end

function _select_common_peaks(p::AbstractVector, q::AbstractVector)
    sel = q .> 0
    psel = p[sel]
    p_sum = sum(psel)
    if p_sum > 0
        psel ./= p_sum
    end
    qsel = q[sel]
    qsel ./= sum(qsel)
    return psel, qsel
end

# geometric distances
julia_euclidean_distance(p,q) = sqrt(sum((p .- q).^2))
function manhattan_distance(p,q) return sum(abs.(p .- q)) end
function chebyshev_distance(p,q) return maximum(abs.(p .- q)) end
function squared_euclidean_distance(p,q) return sum((p .- q).^2) end

# fidelity / matusita / etc
function fidelity_distance(p,q) return 1.0 - sum(sqrt.(p .* q)) end
function matusita_distance(p,q) return sqrt(sum((sqrt.(p) .- sqrt.(q)).^2)) end
function squared_chord_distance(p,q) return sum((sqrt.(p) .- sqrt.(q)).^2) end

function bhattacharya_1_distance(p,q)
    s = sum(sqrt.(p .* q))
    s = min(s, 1.0)
    return (acos(s))^2
end

function bhattacharya_2_distance(p,q)
    s = sum(sqrt.(p .* q))
    if s == 0
        return Inf
    else
        return -log(s)
    end
end

function harmonic_mean_distance(p,q) return 1.0 - 2.0 * sum((p .* q) ./ (p .+ q)) end

function probabilistic_symmetric_chi_squared_distance(p,q)
    return 0.5 * sum(((p .- q).^2) ./ (p .+ q .+ EPS))
end

function ruzicka_distance(p,q)
    return sum(abs.(p .- q)) / sum(max.(p, q))
end

function roberts_distance(p,q)
    return 1.0 - sum((p .+ q) ./ sum(p .+ q) .* min.(p, q) ./ max.(p, q))
end

function intersection_distance(p,q)
    return 1.0 - sum(min.(p,q)) / min(sum(p), sum(q))
end

function motyka_distance(p,q)
    return - (sum(min.(p,q)) / sum(p .+ q))
end

function canberra_distance(p,q)
    return sum(abs.(p .- q) ./ (abs.(p) .+ abs.(q) .+ EPS))
end

function baroni_urbani_buser_distance(p,q)
    # approximate implementation
    if maximum(p) < maximum(q)
        p, q = q, p
    end
    d1 = sqrt(sum(min.(p,q)) * abs.(sum(maximum.(p) .- max.(p,q))))
    return 1.0 - (sum(min.(p,q)) + d1) / (sum(max.(p,q)) + d1)
end

function penrose_size_distance(p,q)
    n = sum(p .> 0)
    return sqrt(n) * sum(abs.(p .- q))
end

function mean_character_distance(p,q)
    n = sum(p .> 0)
    return (1.0 / n) * sum(abs.(p .- q))
end

function lorentzian_distance(p,q) return sum(log.(1 .+ abs.(p .- q))) end

function penrose_shape_distance(p,q)
    p_avg = mean(p); q_avg = mean(q)
    return sqrt(sum(((p .- p_avg) .- (q .- q_avg)).^2))
end

function clark_distance(p,q)
    n = sum(p .> 0)
    return sqrt((1.0 / n) * sum(((p .- q) ./ (abs.(p) .+ abs.(q) .+ EPS)).^2))
end

function hellinger_distance(p,q)
    p_avg = mean(p); q_avg = mean(q)
    return sqrt(2.0 * sum((sqrt.(p ./ p_avg) .- sqrt.(q ./ q_avg)).^2))
end

function whittaker_index_of_association_distance(p,q)
    p_avg = mean(p); q_avg = mean(q)
    return 0.5 * sum(abs.(p ./ p_avg .- q ./ q_avg))
end

function symmetric_chi_squared_distance(p,q)
    p_avg = mean(p); q_avg = mean(q)
    n = sum(p .> 0)
    d1 = (p_avg + q_avg) / (n * (p_avg + q_avg)^2 + EPS)
    return sqrt(d1 * sum(((p * q_avg .- q * p_avg).^2) ./ (p .+ q .+ EPS)))
end

function pearson_correlation_distance(p,q)
    p_avg = mean(p); q_avg = mean(q)
    x = sum((q .- q_avg) .* (p .- p_avg))
    y = sqrt(sum((q .- q_avg).^2) * sum((p .- p_avg).^2))
    if x == 0 && y == 0
        return 0.0
    else
        return - x / (y + EPS)
    end
end

function improved_similarity_distance(p,q)
    n = sum(p .> 0)
    return sqrt((1.0 / n) * sum(((p .- q) ./ (p .+ q .+ EPS)).^2))
end

function absolute_value_distance(p,q)
    return sum(abs.(q .- p)) / (sum(p) + EPS)
end

function dot_product_distance(p,q)
    score = (sum(q .* p))^2 / ((sum(q.^2) + EPS) * (sum(p.^2) + EPS))
    return 1.0 - sqrt(score)
end

function cosine_distance(p,q) return dot_product_distance(p,q) end

function dot_product_reverse_distance(p,q)
    psel, qsel = _select_common_peaks(p,q)
    if sum(psel) == 0
        score = 0.0
    else
        score = (sum(qsel .* psel))^2 / ((sum(qsel.^2) + EPS) * (sum(psel.^2) + EPS))
    end
    return 1.0 - sqrt(score)
end

function spectral_contrast_angle_distance(p,q)
    return 1.0 - sum(q .* p) / (sqrt(sum(q.^2) * sum(p.^2)) + EPS)
end

function wave_hedges_distance(p,q)
    return sum(abs.(p .- q) ./ (max.(p,q) .+ EPS))
end

function jaccard_distance(p,q)
    return sum((p .- q).^2) / (sum(p.^2) + sum(q.^2) - sum(p .* q) + EPS)
end

function dice_distance(p,q)
    return sum((p .- q).^2) / (sum(p.^2) + sum(q.^2) + EPS)
end

function inner_product_distance(p,q) return 1.0 - sum(p .* q) end

function divergence_distance(p,q) return 2.0 * sum(((p .- q).^2) ./ ((p .+ q).^2 .+ EPS)) end

function avg_l_distance(p,q) return (sum(abs.(p .- q)) + maximum(abs.(p .- q))) end

function vicis_symmetric_chi_squared_3_distance(p,q) return sum(((p .- q).^2) ./ (max.(p,q) .+ EPS)) end

# -----------------------------
# ms_distance (port)
# -----------------------------

function weighted_dot_product_distance(spec_query::AbstractArray{T,2}, spec_reference::AbstractArray{T,2}; ms2_ppm=nothing, ms2_da=nothing) where T
    spec_matched = match_peaks_with_mz_info_in_spectra(spec_query, spec_reference; ms2_ppm=ms2_ppm, ms2_da=ms2_da)
    m_q = spec_matched[:,1]; i_q = spec_matched[:,2]
    m_r = spec_matched[:,3]; i_r = spec_matched[:,4]
    k = 0.6; l = 3.0
    w_q = (i_q .^ k) .* (m_q .^ l)
    w_r = (i_r .^ k) .* (m_r .^ l)
    return dot_product_distance(w_q, w_r)
end

function ms_for_id_v1_distance(spec_query::AbstractArray{T,2}, spec_reference::AbstractArray{T,2}; ms2_ppm=nothing, ms2_da=nothing) where T
    spec_matched = match_peaks_in_spectra(spec_query, spec_reference; ms2_ppm=ms2_ppm, ms2_da=ms2_da)
    i_q = spec_matched[:,2]; i_r = spec_matched[:,3]
    n_m = sum((i_q .> 0) .& (i_r .> 0))
    n_q = sum(i_q .> 0); n_r = sum(i_r .> 0)
    a = 0.25
    x = n_m^4
    y = n_q * n_r * (sum(abs.(i_q .- i_r))^a)
    if x == 0
        return Inf
    else
        return y / x
    end
end

function ms_for_id_distance(spec_query::AbstractArray{T,2}, spec_reference::AbstractArray{T,2}; ms2_ppm=nothing, ms2_da=nothing) where T
    if size(spec_query,1) == 0 || size(spec_reference,1) == 0
        return Inf
    end
    # filter intensity > 0.05
    spec_q_f = spec_query[spec_query[:,2] .> 0.05, :]
    spec_r_f = spec_reference[spec_reference[:,2] .> 0.05, :]
    spec_matched = match_peaks_with_mz_info_in_spectra(spec_q_f, spec_r_f; ms2_ppm=ms2_ppm, ms2_da=ms2_da)
    b = 4; c = 1.25; d = 2
    i_q = spec_matched[:,2]; i_r = spec_matched[:,4]
    matched_peak = (i_q .> 0) .& (i_r .> 0)
    n_m = sum(matched_peak)
    n_q = sum(i_q .> 0); n_r = sum(i_r .> 0)
    i_delta = (i_q .- i_r)[matched_peak]
    m_delta = (spec_matched[:,1] .- spec_matched[:,3])[matched_peak]
    s1 = (n_m^b) * (sum(i_q) + 2.0 * sum(i_r))^c
    s2 = (n_q + 2.0 * n_r)^d + sum(abs.(i_delta)) + sum(abs.(m_delta))
    if s2 == 0
        similarity = 0.0
    else
        similarity = s1 / s2
    end
    return -similarity
end

# -----------------------------
# spectral_similarity API
# -----------------------------

# methods name and ranges (partial translation)
const methods_name = Dict(
    "entropy" => "Entropy distance",
    "tsallis_entropy" => "Tsallis Entropy distance",
    "unweighted_entropy" => "Unweighted entropy distance",
    "euclidean" => "Euclidean distance",
    "manhattan" => "Manhattan distance",
    "chebyshev" => "Chebyshev distance",
    "squared_euclidean" => "Squared Euclidean distance",
    "fidelity" => "Fidelity distance",
    "matusita" => "Matusita distance",
    "squared_chord" => "Squared-chord distance",
    "bhattacharya_1" => "Bhattacharya 1 distance",
    "bhattacharya_2" => "Bhattacharya 2 distance",
    "harmonic_mean" => "Harmonic mean distance",
    "probabilistic_symmetric_chi_squared" => "Probabilistic symmetric χ2 distance",
    "ruzicka" => "Ruzicka distance",
    "roberts" => "Roberts distance",
    "intersection" => "Intersection distance",
    "motyka" => "Motyka distance",
    "canberra" => "Canberra distance",
    "baroni_urbani_buser" => "Baroni-Urbani-Buser distance",
    "penrose_size" => "Penrose size distance",
    "mean_character" => "Mean character distance",
    "lorentzian" => "Lorentzian distance",
    "penrose_shape" => "Penrose shape distance",
    "clark" => "Clark distance",
    "hellinger" => "Hellinger distance",
    "whittaker_index_of_association" => "Whittaker index of association distance",
    "symmetric_chi_squared" => "Symmetric χ2 distance",
    "pearson_correlation" => "Pearson/Spearman Correlation Coefficient",
    "improved_similarity" => "Improved Similarity",
    "absolute_value" => "Absolute Value Distance",
    "dot_product" => "Dot product distance",
    "cosine" => "Cosine distance",
    "dot_product_reverse" => "Reverse dot product distance",
    "spectral_contrast_angle" => "Spectral Contrast Angle",
    "wave_hedges" => "Wave Hedges distance",
    "jaccard" => "Jaccard distance",
    "dice" => "Dice distance",
    "inner_product" => "Inner product distance",
    "divergence" => "Divergence distance",
    "avg_l" => "Avg (L1, L∞) distance",
    "vicis_symmetric_chi_squared_3" => "Vicis-Symmetric χ2 3 distance",
    "ms_for_id_v1" => "MSforID distance version 1",
    "ms_for_id" => "MSforID distance",
    "weighted_dot_product" => "Weighted dot product distance",
)

const methods_range = Dict(
    "entropy" => (0.0, log(4.0)),
    "tsallis_entropy" => (0.0, 1.0),
    "unweighted_entropy" => (0.0, log(4.0)),
    "absolute_value" => (0.0, 2.0),
    "avg_l" => (0.0, 1.5),
    "bhattacharya_1" => (0.0, acos(0.0)^2),
    "bhattacharya_2" => (0.0, Inf),
    "canberra" => (0.0, Inf),
    "clark" => (0.0, Inf),
    "divergence" => (0.0, Inf),
    "euclidean" => (0.0, sqrt(2.0)),
    "hellinger" => (0.0, Inf),
    "improved_similarity" => (0.0, Inf),
    "lorentzian" => (0.0, Inf),
    "manhattan" => (0.0, 2.0),
    "matusita" => (0.0, sqrt(2.0)),
    "mean_character" => (0.0, 2.0),
    "motyka" => (-0.5, 0.0),
    "ms_for_id" => (-Inf, 0.0),
    "ms_for_id_v1" => (0.0, Inf),
    "pearson_correlation" => (-1.0, 1.0),
    "penrose_shape" => (0.0, sqrt(2.0)),
    "penrose_size" => (0.0, Inf),
    "probabilistic_symmetric_chi_squared" => (0.0, 1.0),
    "squared_chord" => (0.0, 2.0),
    "squared_euclidean" => (0.0, 2.0),
    "symmetric_chi_squared" => (0.0, 0.5 * sqrt(2.0)),
    "vicis_symmetric_chi_squared_3" => (0.0, 2.0),
    "wave_hedges" => (0.0, Inf),
    "whittaker_index_of_association" => (0.0, Inf)
)

# map method names to functions (for math_distance and ms_distance)
const math_funcs = Dict(
    "tsallis_entropy" => tsallis_entropy_distance,
    "unweighted_entropy" => unweighted_entropy_distance,
    "entropy" => entropy_distance,
    "euclidean" => squared_euclidean_distance, # note: original euclidean returns sqrt, but some code expects pairwise intensity arrays
    "manhattan" => manhattan_distance,
    "chebyshev" => chebyshev_distance,
    "squared_euclidean" => squared_euclidean_distance,
    "fidelity" => fidelity_distance,
    "matusita" => matusita_distance,
    "squared_chord" => squared_chord_distance,
    "bhattacharya_1" => bhattacharya_1_distance,
    "bhattacharya_2" => bhattacharya_2_distance,
    "harmonic_mean" => harmonic_mean_distance,
    "probabilistic_symmetric_chi_squared" => probabilistic_symmetric_chi_squared_distance,
    "ruzicka" => ruzicka_distance,
    "roberts" => roberts_distance,
    "intersection" => intersection_distance,
    "motyka" => motyka_distance,
    "canberra" => canberra_distance,
    "baroni_urbani_buser" => baroni_urbani_buser_distance,
    "penrose_size" => penrose_size_distance,
    "mean_character" => mean_character_distance,
    "lorentzian" => lorentzian_distance,
    "penrose_shape" => penrose_shape_distance,
    "clark" => clark_distance,
    "hellinger" => hellinger_distance,
    "whittaker_index_of_association" => whittaker_index_of_association_distance,
    "symmetric_chi_squared" => symmetric_chi_squared_distance,
    "pearson_correlation" => pearson_correlation_distance,
    "improved_similarity" => improved_similarity_distance,
    "absolute_value" => absolute_value_distance,
    "dot_product" => dot_product_distance,
    "cosine" => cosine_distance,
    "dot_product_reverse" => dot_product_reverse_distance,
    "spectral_contrast_angle" => spectral_contrast_angle_distance,
    "wave_hedges" => wave_hedges_distance,
    "jaccard" => jaccard_distance,
    "dice" => dice_distance,
    "inner_product" => inner_product_distance,
    "divergence" => divergence_distance,
    "avg_l" => avg_l_distance,
    "vicis_symmetric_chi_squared_3" => vicis_symmetric_chi_squared_3_distance
)

const ms_funcs = Dict(
    "ms_for_id_v1" => ms_for_id_v1_distance,
    "ms_for_id" => ms_for_id_distance,
    "weighted_dot_product" => weighted_dot_product_distance
)

# API functions

function distance(spectrum_query::AbstractArray{T,2}, spectrum_library::AbstractArray{T,2}, method::String; ms2_ppm=nothing, ms2_da=nothing, need_clean_spectra::Bool=false, need_normalize_result::Bool=true, entropy_dimension::Real=2.0) where T
    if ms2_ppm === nothing && ms2_da === nothing
        error("MS2 tolerance need to be defined!")
    end
    sq = Array{Float32,2}(spectrum_query)
    sl = Array{Float32,2}(spectrum_library)
    if need_clean_spectra
        sq = clean_spectrum(sq; ms2_ppm=ms2_ppm, ms2_da=ms2_da)
        sl = clean_spectrum(sl; ms2_ppm=ms2_ppm, ms2_da=ms2_da)
    end
    if size(sq,1) > 0 && size(sl,1) > 0
        fname = method
        if haskey(math_funcs, fname)
            f = math_funcs[fname]
            spec_matched = match_peaks_in_spectra(sq, sl; ms2_ppm=ms2_ppm, ms2_da=ms2_da)
            if fname == "tsallis_entropy"
                return f(spec_matched[:,2], spec_matched[:,3], entropy_dimension)
            else
                if need_normalize_result
                    dist_range = haskey(methods_range, method) ? methods_range[method] : (0.0, 1.0)
                    dist = f(spec_matched[:,2], spec_matched[:,3])
                    res_dist = normalize_distance(dist, dist_range)
                    return res_dist
                else 
                    res_dist = f(spec_matched[:,2], spec_matched[:,3])
                    return res_dist
                end
                
            end
        elseif haskey(ms_funcs, fname)
            f = ms_funcs[fname]
            return f(sq, sl; ms2_ppm=ms2_ppm, ms2_da=ms2_da)
        else
            error("Method name: $method error!")
        end
    else
        return need_normalize_result ? 1.0 : Inf
    end
end

function all_distance(spectrum_query::AbstractArray{T,2}, spectrum_library::AbstractArray{T,2}; ms2_ppm=nothing, ms2_da=nothing, need_clean_spectra::Bool=true, need_normalize_result::Bool=true, entropy_dimension::Real=2.0) where T
    if ms2_ppm === nothing && ms2_da === nothing
        error("MS2 tolerance need to be defined!")
    end
    sq = Array{Float32,2}(spectrum_query)
    sl = Array{Float32,2}(spectrum_library)
    if need_clean_spectra
        sq = clean_spectrum(sq; ms2_ppm=ms2_ppm, ms2_da=ms2_da)
        sl = clean_spectrum(sl; ms2_ppm=ms2_ppm, ms2_da=ms2_da)
    end
    result = Dict{String,Float64}()
    if size(sq,1) > 0 && size(sl,1) > 0
        spec_matched = match_peaks_in_spectra(sq, sl; ms2_ppm=ms2_ppm, ms2_da=ms2_da)
        for (method, _) in methods_name
            try
                if haskey(math_funcs, method)
                    f = math_funcs[method]
                    if method == "tsallis_entropy"
                        dist = f(spec_matched[:,2], spec_matched[:,3], entropy_dimension)
                    else
                        dist = f(spec_matched[:,2], spec_matched[:,3])
                    end
                elseif haskey(ms_funcs, method)
                    f = ms_funcs[method]
                    dist = f(sq, sl; ms2_ppm=ms2_ppm, ms2_da=ms2_da)
                else
                    error("Method name: $method error!")
                end
                if need_normalize_result
                    dist_range = haskey(methods_range, method) ? methods_range[method] : (0.0, 1.0)
                    dist = normalize_distance(dist, dist_range)
                end
                result[method] = float(dist)
            catch e
                # on failure, set Inf or 1 depending on normalization
                result[method] = need_normalize_result ? 1.0 : Inf
            end
        end
    else
        for (method, _) in methods_name
            result[method] = need_normalize_result ? 1.0 : Inf
        end
    end
    return result
end

function similarity(spectrum_query::AbstractArray{T,2}, spectrum_library::AbstractArray{T,2}, method::String; ms2_ppm=nothing, ms2_da=nothing, need_clean_spectra::Bool=true, need_normalize_result::Bool=true, entropy_dimension::Real=2.0) where T
    if need_normalize_result
        return 1.0 - distance(spectrum_query, spectrum_library, method; ms2_ppm=ms2_ppm, ms2_da=ms2_da, need_clean_spectra=need_clean_spectra, need_normalize_result=need_normalize_result, entropy_dimension=entropy_dimension)
    else
        return 0.0 - distance(spectrum_query, spectrum_library, method; ms2_ppm=ms2_ppm, ms2_da=ms2_da, need_clean_spectra=need_clean_spectra, need_normalize_result=need_normalize_result, entropy_dimension=entropy_dimension)
    end
end

function all_similarity(spectrum_query::AbstractArray{T,2}, spectrum_library::AbstractArray{T,2}; ms2_ppm=nothing, ms2_da=nothing, need_clean_spectra::Bool=true, need_normalize_result::Bool=true) where T
    all_similarity_score = all_distance(spectrum_query, spectrum_library; need_clean_spectra=need_clean_spectra, need_normalize_result=need_normalize_result, ms2_ppm=ms2_ppm, ms2_da=ms2_da)
    for (m, v) in collect(all_similarity_score)
        all_similarity_score[m] = need_normalize_result ? 1.0 - v : 0.0 - v
    end
    return all_similarity_score
end

function multiple_distance(spectrum_query::AbstractArray{T,2}, spectrum_library::AbstractArray{T,2}, methods::Union{Nothing,Vector{String}}=nothing; ms2_ppm=nothing, ms2_da=nothing, need_clean_spectra::Bool=true, need_normalize_result::Bool=true) where T
    if methods !== nothing
        result = Dict{String,Float64}()
        for m in methods
            d = distance(spectrum_query, spectrum_library, m; need_clean_spectra=need_clean_spectra, need_normalize_result=need_normalize_result, ms2_ppm=ms2_ppm, ms2_da=ms2_da)
            result[m] = float(d)
        end
        return result
    else
        return all_distance(spectrum_query, spectrum_library; ms2_ppm=ms2_ppm, ms2_da=ms2_da, need_clean_spectra=need_clean_spectra, need_normalize_result=need_normalize_result)
    end
end

function multiple_similarity(spectrum_query::AbstractArray{T,2}, spectrum_library::AbstractArray{T,2}, methods::Union{Nothing,Vector{String}}=nothing; ms2_ppm=nothing, ms2_da=nothing, need_clean_spectra::Bool=true, need_normalize_result::Bool=true) where T
    if methods !== nothing
        result = Dict{String,Float64}()
        for m in methods
            s = similarity(spectrum_query, spectrum_library, m; need_clean_spectra=need_clean_spectra, need_normalize_result=need_normalize_result, ms2_ppm=ms2_ppm, ms2_da=ms2_da)
            result[m] = float(s)
        end
        return result
    else
        return all_similarity(spectrum_query, spectrum_library; need_clean_spectra=need_clean_spectra, need_normalize_result=need_normalize_result, ms2_ppm=ms2_ppm, ms2_da=ms2_da)
    end
end

end # module
