# Defining a structure for left-open right-closed intervals, that is, an interval 
# including all real numbers greater than "a" and less than or equal to "b"
struct LeftOpenRightClosedInterval{R}
    a::R
    b::R
	LeftOpenRightClosedInterval(a::R,b::R) where {R <: Real} = begin @assert a ≤ b; new{R}(a,b) end
end

# Defining an "∈" operator for a "LeftOpenRightClosedInterval" to check if a given 
# real number "x" lies within the considered interval
begin
	import Base: 
		∈

	∈(x::Q, I::LeftOpenRightClosedInterval{R}) where {Q <: Real, R <: Real} = begin
		return I.a < x ≤ I.b
	end
end

# Defining helper functions to retrieve the starting point, ending point, as well 
# as both the starting and ending point simultaneously, of a "LeftOpenRightClosedInterval".
begin 
	endpoint(I::LeftOpenRightClosedInterval) = I.b
	startpoint(I::LeftOpenRightClosedInterval) = I.a
	points(I::LeftOpenRightClosedInterval) = [I.a, I.b]
end

# Defining a structure that represents a supply/demand curve as a function of price and cumulative quantity
begin 
	struct SupplyDemandCurve{R,Q} <: Function
    	Price::AbstractVector{R}
    	Intervals::AbstractVector{LeftOpenRightClosedInterval{Q}}
	end

	function SupplyDemandCurve(Price::AbstractVector{R}, CumulativeQuantity::AbstractVector{Q}) where {R, Q <: Real} 
	    @assert length(Price) == length(CumulativeQuantity)
	    
        Intervals = LeftOpenRightClosedInterval{Q}[]
	    for i in eachindex(CumulativeQuantity[1:end-1])
			push!(Intervals, LeftOpenRightClosedInterval(CumulativeQuantity[i], CumulativeQuantity[i+1]))
	    end
	    SupplyDemandCurve{R,Q}(Price, Intervals)
	end
end

# Defining how the "SupplyDemandCurve" function should be evaluated at a given quantity "x".
# Specifically, it uses the intervals and prices of the SupplyDemandCurve to determine the appropriate price for the given quantity.
# If the quantity "x" lies within the "k"-th interval, the function returns the "k+1"-th price.
# If "x" is equal to the starting point of the first interval, the function returns the first price. 
# If "x" does not lie within any interval, the function returns zero to avoid errors.
# NOTE: The implementation below differs a little from how it is described in the thesis, particularly, in terms of the indices.
function (C::SupplyDemandCurve{R,Q})(x::S) where {R, Q <: Real, S <: Real}
    Price = C.Price
    Intervals = C.Intervals
	x0 = startpoint(Intervals[1])

	if x ≈ x0 
		return Price[1]
	end

    for k in eachindex(Intervals)
        if x ∈ Intervals[k]
            return Price[k+1] 
        end
    end

    return zero(R)
end

############################################################# EXAMPLE OF USAGE #############################################################
let
    # Generating artificial prices and cumulative quantities.
    # NOTE: the artificial prices only differ in terms of sorting; supply = ascending and demand = descending.
    # In practice, supply and demand prices of course also differ in terms of their actual value, however, this is 
    # merely meant to serve as an illustrative example for how supply/demand curves may be created using the "SupplyDemandCurve" function.
    SupplyPrice = sort([0., 2., 1., 0.5, 0.4, 0.3, 0.8, 0.9, 2.1, 1.5, 1.3, 0.9, 1.5, 1.8])
    DemandPrice = sort([0., 2., 1., 0.5, 0.4, 0.3, 0.8, 0.9, 2.1, 1.5, 1.3, 0.9, 1.5, 1.8], rev = true)
    
    SupplyCumulativeQuantity = [0., 0.14, 0.27, 0.35, 0.44, 0.53, 0.58, 0.63, 0.71, 0.73, 0.78, 0.87, 0.94, 0.99]
    DemandCumulativeQuantity = [0., 0.23, 0.35, 0.37, 0.41, 0.5, 0.52, 0.74, 0.88, 0.94, 1., 1.28, 1.43, 1.57]

    SupplyCurveExample = SupplyDemandCurve(SupplyPrice, SupplyCumulativeQuantity)
    DemandCurveExample = SupplyDemandCurve(DemandPrice, DemandCumulativeQuantity)
    
    # Visualizing example with "×" indicating sampling points
    using Plots
    
    default(size = (625, 385))
    default(fontfamily = "Computer Modern", titlefontsize = 11, guidefontsize = 11, tickfontsize = 9, legendfontsize = 9, linewidth = 2)

    n_samplingpoints = 81

    plot(LinRange(0, endpoint(SupplyCurveExample.Intervals[end]), n_samplingpoints), SupplyCurveExample, linecolor = :green,  markershape = :xcross, markersize = 2.5, markerstrokewidth = 2, markercolor = :black, label = "Supply")
    plot!(LinRange(0, endpoint(DemandCurveExample.Intervals[end]), n_samplingpoints), DemandCurveExample, linecolor = :red, markershape = :xcross, markersize = 2.5, markerstrokewidth = 2, markercolor = :black, label = "Demand")
    plot!(xlabel = "Quantity", ylabel = "Price", grid = false)
end

