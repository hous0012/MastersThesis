using LinearAlgebra, Statistics, Dates, DataFrames

# Defining structure to save results in bFPCA procedure
struct bFPCA_struct
    mean::Vector{Vector}
    cov::Matrix
    eigenvalues::Vector
    eigenfunctions::Vector{Vector{Vector}}
    scores::Vector{Vector}
    FVE::Vector
    cumFVE::Vector
end

# Defining function to carry out bFPCA procedure
function bFPCA(fundata::DataFrame, n_gridpoints::Int; weights = false)
    # Ensuring that "fundata" is constructed correctly specified
    @assert eltype(fundata[:,1]) == Date
    @assert eltype(fundata[:,2]) == Z  

    # Misc
    n_samples = nrow(fundata);

    # Evaluating the supply/demand curves on grid of [0, 1], with grid specified by "n_gridpoints"
    t0 = LinRange(0, 1, n_gridpoints+1)
    eval_curve = [map(f -> f(t0[i]), fundata[:,2]) for i in eachindex(t0[2:end])]
    
    # Estimating mean
    μ = mean.(eval_curve);

    # Estimating covariance
    tildeX = [vcat([X[i] for X in eval_curve]...) - vcat(μ...) for i in 1:n_samples];
    boldX = Matrix(undef, n_samples, 2 * n_gridpoints);
    
    for i in 1:n_samples
        boldX[i,:] = tildeX[i]';
    end

    V = 1 / (n_samples - 1) .* boldX'boldX;

    # Computing partition length
    Δ = 1 / n_gridpoints;

    # Constructing weight matrix
    if(weights == true)
        # Recovering variance componentwise
        K11 = diag(V)[1:2:size(V, 1)];
        K22 = diag(V)[2:2:size(V, 1)];

        # Integrating (by discretization) componentwise
        intK11 = Δ * sum(K11);
        intK22 = Δ * sum(K22);

        # Inverting componentwise
        w1 = 1 / intK11;
        w2 = 1 / intK22;

        # Constructing relevant matrices
        W = Diagonal(repeat([w1, w2], outer = n_gridpoints));
        sqrtW = Diagonal(repeat([sqrt(w1), sqrt(w2)], outer = n_gridpoints));
        invsqrtW = Diagonal(repeat([1 / sqrt(w1), 1 / sqrt(w2)], outer = n_gridpoints));
    else
        W = sqrtW = invsqrtW = I    # Set weight matrices equal to identity
    end

    # Solving eigenproblem
    eig = eigen(sqrtW * V * sqrtW);

    # Recovering eigenfunctions
    all_eigf = [[(invsqrtW * real.(eig.vectors))[i:i+1, j]  for i in 1:2:size(V, 1)] for j in axes(V, 2)];

    # Recovering eigenvalues; recall must be non-negative.
    # Also, sorting eigenvalues from largest to smallest (descending).
    # Note that Julia follows LAPACK so eigenvalues are by 
    # default returned in an ascending order.
    all_eigv = (Δ * real(eig.values));
    ind_sorted_nn_eigv = reverse(findall(f -> f > 0, all_eigv));
    eigv = all_eigv[ind_sorted_nn_eigv];
    eigf = all_eigf[ind_sorted_nn_eigv];

    # Computing scores
    scores = [boldX * W * vcat(eigf[i]...) for i eachindex(ind_sorted_nn_eigv)];
    
    # FVE
    FVE = [eigv[i] / sum(eigv) for i in eachindex(eigv)];

    # Cumulative FVE
    cumFVE = cumsum(eigv) ./ sum(eigv);

    return bFPCA_struct(μ, V, eigv, eigf, scores, FVE, cumFVE)
end

############################################################# EXAMPLE OF USAGE #############################################################

# Generating artificial prices and cumulative quantities used to create (reparameterized) supply and demand curves
SupplyPrice1 = sort([0., 2., 1., 0.5, 0.4, 0.3, 0.8, 0.9, 2.1, 1.5, 1.3, 0.9, 1.5, 1.8])
SupplyPrice2 = sort([0.2, 1.4, 0.9, 0.7, 0.5, 1., 1.3, 1.2, 1.5, 1.8, 1.7, 1.9, 2, 1.8])
SupplyPrice3 = sort([0.5, 0.5, 0.3, 1., 1.2, 1., 0.6, 0.9, 0.5, 1.8, 1.4, 1.5, 1.2, 1.5])

SupplyCumulativeQuantity1 = [0., 0.14, 0.27, 0.35, 0.44, 0.53, 0.58, 0.63, 0.71, 0.73, 0.78, 0.87, 0.94, 0.99]
SupplyCumulativeQuantity2 = [0., 0.03, 0.1, 0.17, 0.25, 0.44, 0.58, 0.67, 0.84, 0.88, 0.95, 0.99, 1.22, 1.32]
SupplyCumulativeQuantity3 = [0., 0.28, 0.33, 0.38, 0.44, 0.59, 0.74, 0.91, 0.97, 1.12, 1.33, 1.49, 1.57, 1.73]

SupplyCurveExample1 = SupplyDemandCurve(SupplyPrice1, SupplyCumulativeQuantity1)
SupplyCurveExample2 = SupplyDemandCurve(SupplyPrice2, SupplyCumulativeQuantity2)
SupplyCurveExample3 = SupplyDemandCurve(SupplyPrice3, SupplyCumulativeQuantity3)

ReparameterizedSupplyCurveExample1 = ζ(SupplyCurveExample1)
ReparameterizedSupplyCurveExample2 = ζ(SupplyCurveExample2)
ReparameterizedSupplyCurveExample3 = ζ(SupplyCurveExample3)

ReparameterizedSupplyCurveExamples = [ReparameterizedSupplyCurveExample1, ReparameterizedSupplyCurveExample2, ReparameterizedSupplyCurveExample3]

# Creating an example of "fundata"
n_samples = 3
fundataExample = DataFrame(Date = Vector{Date}(undef, n_samples), SDCurve = Vector{Z}(undef, n_samples))

for i in 1:n_samples
    fundataExample[i,:]["Date"] = Date(2023, 6, i)
    fundataExample[i,:]["SDCurve"] = ReparameterizedSupplyCurveExamples[i]
end

# Performing bFPCA
fundata = fundataExample
n_gridpoints = 100
n_fpcs = 3

res = bFPCA(fundata, n_gridpoints, n_fpcs);

## Visualizing results
t0 = LinRange(0, 1, n_gridpoints+1)  
t0 = t0[2:end]

# Mean
let
    using Plots

    default(size = (625, 385))
    default(fontfamily = "Computer Modern", titlefontsize = 11, guidefontsize = 11, tickfontsize = 9, legendfontsize = 9, linewidth = 2)

    plot([], [], label = "")
    plot!(xlim = [0,2], ylim = [0,2])

    labels = ["SupplyExample$i" for i in 1:3]
    for i in 1:3
        plot!(map(fundataExample[!, "SDCurve"][i].η, t0), map(fundataExample[!, "SDCurve"][i].SDCurve, t0), linealpha = 1, label = labels[i])
    end

    plot!([v[1] for v in res.mean], [v[2] for v in res.mean], line = "black", label = "Mean")
    plot!(xlabel = "Quantity", ylabel = "Price", grid = false)   
end

# Covariance
let
    using Plots

    default(size = (600, 400))
    default(fontfamily = "Computer Modern", titlefontsize = 9, guidefontsize = 9, tickfontsize = 5, legendfontsize = 9, linewidth = 2)

    foo = size(res.cov, 1)

    Cov11 = res.cov[1:2:foo, 1:2:foo]
    Cov22 = res.cov[2:2:foo, 2:2:foo]
    Cov12 = res.cov[2:2:foo, 1:2:foo]
    Cov21 = res.cov[1:2:foo, 2:2:foo]

    Y = reshape(repeat(t0, outer = n_gridpoints), n_gridpoints, :);
    X = Y';

    p1 = surface(X, Y, Cov11, legend = 0, title = "Cov: Quantity")
    p2 = surface(X, Y, Cov12, legend = 0, title = "Cross-cov: Quantity and Price")
    p3 = surface(X, Y, Cov21, legend = 0, title = "Cross-cov: Price and Quantity")
    p4 = surface(Y, X, Cov22, legend = 0, title = "Cov: Price")

    plot(p1, p2, p3, p4, layout = (2, 2))
    plot!(colorbar = false)
end

# Eigenvalue and cumulative FVE
let
    using Plots

    default(size = (625, 385))
    default(fontfamily = "Computer Modern", titlefontsize = 11, guidefontsize = 11, tickfontsize = 9, legendfontsize = 9, linewidth = 2)
    
    p1 = plot(res.eigenvalues[1:10], linecolor = :black, markershape = :circle, markercolor = :white, markerstrokecolor = :black, markerstrokewidth = 1.5, xlabel = "Component Number", ylabel = "Eigenvalue")
    p2 = plot(res.cumFVE[1:10], linecolor = :black, markershape = :circle, markercolor = :white, markerstrokecolor = :black, markerstrokewidth = 1.5, xlabel = "Component Number", ylabel = "Cumulative FVE")

    plot(p1, p2, layout = (1, 2), grid = false, legend = false)
end

# Eigenfunctions
let
    using Plots

    default(size = (625, 385))
    default(fontfamily = "Computer Modern", titlefontsize = 11, guidefontsize = 11, tickfontsize = 9, legendfontsize = 9, linewidth = 2)

    p1 = plot([v[1] for v in res.eigenfunctions[1]], [v[2] for v in res.eigenfunctions[1]], linecolor = :black, xlabel = "Quantity", ylabel = "Price", label = "Eigenfunction1")
    p2 = plot([v[1] for v in res.eigenfunctions[2]], [v[2] for v in res.eigenfunctions[2]], linecolor = :black, xlabel = "Quantity", ylabel = "Price", label = "Eigenfunction2")

    plot(p1, p2, layout = (1, 2), grid = false)
end

# Scores
let
    using Plots

    default(size = (625, 385))
    default(fontfamily = "Computer Modern", titlefontsize = 11, guidefontsize = 11, tickfontsize = 9, legendfontsize = 9, linewidth = 2)

    labels = ["Score$i" for i in 1:n_fpcs]

    plot([], [], label = "")
    for i in 1:n_fpcs
        plot!(res.scores[i], label = labels[i])
    end
    plot!(xlabel = "Sample No.", grid = false)
end

# Reconstruction of data
let
    using Plots

    default(size = (600, 400))
    default(fontfamily = "Computer Modern", titlefontsize = 11, guidefontsize = 11, tickfontsize = 9, legendfontsize = 9, linewidth = 2)

    SampleNo = 3
    TruncationOrder = 10

    foo_res = bFPCA(fundata, n_gridpoints, TruncationOrder) # Run bFPCA procedure again so that "n_fpcs = TruncationOrder"

    KarhunenLoeve1 = [v[1] for v in foo_res.mean]
    KarhunenLoeve2 = [v[2] for v in foo_res.mean]

    for i in 1:TruncationOrder
        KarhunenLoeve1 += foo_res.scores[i][SampleNo]*[v[1] for v in foo_res.eigenfunctions[i]]
        KarhunenLoeve2 += foo_res.scores[i][SampleNo]*[v[2] for v in foo_res.eigenfunctions[i]]
    end

    plot(map(fundataExample[!, "SDCurve"][SampleNo].η, t0), map(fundataExample[!, "SDCurve"][SampleNo].SDCurve, t0), linecolor = :black, label = "SupplyExample$SampleNo", xlabel = "Quantity", ylabel = "Price")
    plot!(KarhunenLoeve1, KarhunenLoeve2, line = (:dash, "red"), label = "Reconstructed")
    plot!(grid = false)
end

