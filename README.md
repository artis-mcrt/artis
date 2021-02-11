# Artistools

ARTIS (Sim et al. 2007; Kromer & Sim 2009) is a 3D radiative transfer code for Type Ia supernovae using the Monte Carlo method with indivisible energy packets (Lucy 2002). This version includes non-LTE physics appropriate for the nebular phase described in Shingles et al. 2020.

## model.txt (1D case)
For a 1D model, we define the density and radioactive mass fractions for spherical shells as as a function of velocity.
```
line 1: [npts_model: number of cells]
line 2: [t_model_days: time after explosion in days]
lines >3: [mgi: cell number]   [shell outer velocity in km/s]  [log10(density in g/cm^3)]  [X_Fegroup]  [X_Ni56]  [X_Co56]  [X_Fe52]  [X_Cr48]
```

## abundances.txt
The elemental abundances are specified as mass fractions (or equivalently element masses or densities, as they will be normalised to 1. after input) in each of the model cells defined in the model.txt file.
```
lines: [mgi (cell number)]   [X_H]  [X_He]  ...  [X_Zn]
```
