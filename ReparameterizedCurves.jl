begin
    # Define a structure encompasing all "η" functions
	struct H{Q <: Real} <: Function
		a::Q
		b::Q 
	end

    # This function "η" constructs an instance of "H" given a "SupplyDemandCurve". 
	# The interval of the function starts from the "startpoint" of the first interval
    # and ends in the "endpoint" of the last interval of the "SupplyDemandCurve"
	η(SDCurve::SupplyDemandCurve{R,Q}) where {R, Q <: Real} = begin
		a, b = startpoint(SDCurve.Intervals[1]), endpoint(SDCurve.Intervals[end])
		H(a,b)
    end

    # Defining how the "η" function should be evaluated at a given "t".
	function (η::H{Q})(t::Real) where {Q}
		a, b = η.a, η.b
		(b - a) * t + a
	end

    # Defining an "∘" operator to compose a "SupplyDemandCurve" and a function "η". 
	# The result is a new "SupplyDemandCurve" where the quantities have been scaled to the interval [0, 1].
    import Base:
			∘
    
	function ∘(SDCurve::SupplyDemandCurve{R,Q}, η::H{Q}) where {R, Q <: Real}
		Price, Intervals = SDCurve.Price, SDCurve.Intervals
		a, b = η.a, η.b

		if a == b
			return SupplyDemandCurve([Price[1], Price[1]], [0., 1.])
		else
		    CumulativeQuantity = [startpoint(Intervals[1]), map(endpoint, Intervals)...] |> unique
			𝑡 = map(x -> (x - a) / (b - a), CumulativeQuantity)
			return SupplyDemandCurve(Price, 𝑡)
		end
	end	
end

begin
    # Define a structure encompasing all "ζ" functions
	struct Z{R,Q} <: Function
		η::H{Q}
		SDCurve::SupplyDemandCurve{R,Q}
	end

    # Defining how the "ζ" should be evaluated at a given "t".
    # The result is a pair of real numbers, representing the quantity and price at 't' on the supply/demand curve.
	function (ζ::Z)(t::Real)
		@assert 0 ≤ t ≤ 1
		η, SDCurve = ζ.η, ζ.SDCurve

		return [η(t), SDCurve(t)]
	end

    # This function "ζ" constructs an instance of "Z" given a "SupplyDemandCurve".
    # The components are two-fold; a function "η" and a new "SupplyDemandCurve" function 
    # given by the composition of the function "η" and the original "SupplyDemandCurve".
	function ζ(SDCurve::SupplyDemandCurve{R,Q}) where {R,Q}
		Z{R,Q}(η(SDCurve), SDCurve ∘ η(SDCurve))
	end
end

############################################################# EXAMPLE OF USAGE #############################################################
let
    # Continuation of example in "SupplyDemandCurves.jl"
    SupplyPrice = sort([0., 2., 1., 0.5, 0.4, 0.3, 0.8, 0.9, 2.1, 1.5, 1.3, 0.9, 1.5, 1.8])
    DemandPrice = sort([0., 2., 1., 0.5, 0.4, 0.3, 0.8, 0.9, 2.1, 1.5, 1.3, 0.9, 1.5, 1.8], rev = true)

    SupplyCumulativeQuantity = [0., 0.14, 0.27, 0.35, 0.44, 0.53, 0.58, 0.63, 0.71, 0.73, 0.78, 0.87, 0.94, 0.99]
    DemandCumulativeQuantity = [0., 0.23, 0.35, 0.37, 0.41, 0.5, 0.52, 0.74, 0.88, 0.94, 1., 1.28, 1.43, 1.57]

    SupplyCurveExample = SupplyDemandCurve(SupplyPrice, SupplyCumulativeQuantity)
    DemandCurveExample = SupplyDemandCurve(DemandPrice, DemandCumulativeQuantity)

    # Reparameterized supply and demand curves
    ReparameterizedSupplyCurveExample = ζ(SupplyCurveExample)
    ReparameterizedDemandCurveExample = ζ(DemandCurveExample)

    # Visualizing example with "×" indicating sampling points
    using Plots

    default(size = (625, 385))
    default(fontfamily = "Computer Modern", titlefontsize = 11, guidefontsize = 11, tickfontsize = 9, legendfontsize = 9, linewidth = 2)

    n_samplingpoints = 81

    t = LinRange(0, 1, n_samplingpoints);
    eval_supply = map(ReparameterizedSupplyCurveExample, t)
    eval_demand = map(ReparameterizedDemandCurveExample, t)

    plot(t, [elem[1] for elem in eval_supply], [elem[2] for elem in eval_supply], linecolor = :green,  markershape = :xcross, markersize = 1, markerstrokewidth = 2, markercolor = :black, label = "Supply")
    plot!(t, [elem[1] for elem in eval_demand], [elem[2] for elem in eval_demand], linecolor = :red,  markershape = :xcross, markersize = 1, markerstrokewidth = 2, markercolor = :black, label = "Demand")

    plot!(xlabel = "'Time'", ylabel = "Quantity", zlabel = "Price")
end

