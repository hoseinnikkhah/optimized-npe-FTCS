# Electrokinetic Transport Model

## Overview
This MATLAB code implements an electrokinetic transport model for simulating the migration of contaminants through porous media under the influence of electroosmotic flow. The model focuses on the transport of Be-7 radionuclides in soil and incorporates various physical and chemical parameters to represent real-world conditions.

## Key Features
- Solves a one-dimensional advection-diffusion equation for contaminant transport.
- Considers electroosmotic flow, voltage gradients, zeta potential, and other electrokinetic effects.
- Calculates electroosmotic mobility using different models (Alshawabkeh, Vane, Shapiro).
- Incorporates temperature-dependent viscosity and permittivity of different media (oil, water, clay).
- Accounts for porosity, tortuosity, valency, and density of clay particles.

## Usage
The code initializes and solves the transport equations numerically using finite difference methods. It outputs concentration profiles for different contaminants over time. The results are then visualized using MATLAB plots to compare the original finite difference method (FDM) with an optimized FDM.

## Parameters
- Length of domain (L): 0.4 m
- Initial inventory of Be-7 (I0): 0.0000829 Bq m^-2 day^-1
- Decay constant (R1): 0.013 day^-1
- Simulation end time (tmax): 35 days
- Number of nodes (nx) and time steps (nt)
- Various physical and chemical properties of media (viscosity, permittivity, etc.)

## Applications
This model can be used to study the transport behavior of radionuclides and other contaminants in soil under the influence of electrokinetic processes. It finds applications in environmental science, radioactive waste management, and understanding subsurface contaminant migration.

## Results
The code generates plots comparing the concentration profiles of contaminants under original and optimized finite difference methods for different initial concentrations.

![FTCS_fig](https://github.com/hoseinnikkhah/optimized-npe-FTCS/assets/116885462/319034ea-d942-4b15-bf24-87b9f7db9be7)


## Contributors
- [Hosein Nikkkhah]

Feel free to contribute to the code and share your insights!

