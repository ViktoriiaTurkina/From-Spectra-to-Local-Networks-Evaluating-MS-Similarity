using CSV
using DataFrames
using Plots
using StatsBase
using PlotlyJS
using NetworkLayout
using ScikitLearn
using Clustering
using GraphPlot
using Graphs
using Distances
using MultivariateStats
using StatsPlots
using StringDistances
using PooledArrays
using Random
using BSON
using Conda
using PyCall
using ProgressBars
using Distances
using ColorSchemes
using .SpectralSimilarity


import StatsPlots as sp

me = pyimport("ms_entropy")
np = pyimport("numpy")


ms2deepscore = pyimport("matchms.calculate_scores")
matchms_spectrum = pyimport("matchms")
ms2deepscore_model = pyimport("ms2deepscore.models")
sim_measure = pyimport("ms2deepscore")




    # using Pkg, Conda, PyCall

    #     Conda.pip_interop(true,  :my_env)
    #     Conda.pip("install", "matchms")
    #     Conda.pip("install", "ms2deepscore")
        

    #     Pkg.build("PyCall")


