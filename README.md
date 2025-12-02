# simple-energy-system-optimization
Simple MATLAB techno-economic optimisation model for a small power system (demo code).

This repository contains a small MATLAB script that builds a simple techno-economic
capacity expansion model for a power system with solar, wind and gas.

The model:
 Creates synthetic hourly profiles for demand, solar and wind.
 Defines investment and variable costs for each technology.
 Uses linear programming (`linprog`, Optimization Toolbox) to find the least-cost
  combination of capacities and hourly dispatch over one week.
 Outputs optimal capacities, system cost and an example dispatch plot.

This code is part of my portfolio and illustrates how I work with
modelling, optimisation and energy system analysis in MATLAB.
