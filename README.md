# Master's Thesis - Computational Resources

This repository hosts the scripts, developed in the Julia programming language, intended to supplement my master's thesis.

## Repository Structure

The repository contains three scripts, each mirroring a particular theoretical discourse presented in my master's thesis:
  - `SupplyDemandCurves.jl` allows for the construction of supply and demand curves based on user-specified market bids;
  - `ReparameterizedCurves.jl` facilitates the reparameterization of user-specified supply and demand curves. These curves need to be constructed using the `SupplyDemandCurves.jl` script;
  - `bFPCA.jl` may be used to perform bFPCA (Bivariate Functional Principal Components Analysis) based on user-specified bivariate functional data. The bivariate functional data needs to be constructed using the `ReparameterizedCurves.jl` script.

Each script includes an illustrative example demonstrating its usage.
  
