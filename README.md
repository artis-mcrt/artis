# Artistools

ARTIS (Sim et al. 2007; Kromer & Sim 2009) is a 3D radiative transfer code for Type Ia supernovae using the Monte Carlo method with indivisible energy packets (Lucy 2002). This version includes non-LTE physics appropriate for the nebular phase described in Shingles et al. 2020.

## model.txt (1D case)
line 1: number of cells
line 2: time after explosion in days
lines >3: [cell number]   [velocity in km/s]  [log10(density in g/cm^3)]  [X_Fegroup]  [X_Ni56]  [X_Co56]  [X_Fe52]  [X_Cr48]