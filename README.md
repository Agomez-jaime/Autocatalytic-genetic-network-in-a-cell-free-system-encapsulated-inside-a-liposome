# Autocatalytic Genetic Network in a Liposome

Monte Carlo simulations of an autocatalytic DNA replication network encapsulated inside a cell-free liposome, modeling the competition between functional DNA and parasitic DNA across successive bulk-replication cycles.

## What this simulates

A minimal protocell model in which:

- **Functional DNA (`dna`)** is replicated by a DNA polymerase (`dnap`) using an RNA polymerase template (`tp`) and a finite resource pool (`res`).
- **Parasitic DNA (`para`)** can emerge stochastically (probability `pb` per cycle) and competes for the same polymerase, hijacking replication resources without contributing useful function.
- After each cycle the contents of the liposome are **diluted** by a factor `dil` and propagated to the next cycle, mimicking serial-passage protocell experiments.

The simulation sweeps a logarithmic range of parasite-emergence probabilities and dilution factors, running 1000 Monte Carlo realizations per condition. For each run it records whether the functional DNA survives, goes extinct, or is overtaken by parasites — yielding phase diagrams of replicator persistence as a function of bulk parameters (volume, dilution, parasite-emergence rate, integration time).

## Files

| File | Description |
|---|---|
| `Monte_Carlo_Bulk_Simulation.m` | MATLAB driver: defines bulk parameters, runs the stochastic loop, and saves per-cycle concentrations of `dna`, `dnap`, `dnapdna`, `tp`, `para`, `dnappara`, `res`. |

## Running

Requires MATLAB. Open `Monte_Carlo_Bulk_Simulation.m` and run it (or hit `Ctrl+Enter` on each section to step through). Adjust the parameter sweep arrays at the top of the file (`probabilities`, the random `dil`/`tim` ranges, the number of Monte Carlo realizations) to explore different regimes.

## Credit

Underlying replicator model: **Andreea Stan, TU Delft**. This repository contains the bulk Monte Carlo extension and parameter-sweep driver used in the *Autocatalytic genetic networks in cell-free systems encapsulated inside liposomes* project.
