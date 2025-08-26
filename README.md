# Electric Machine Optimization with MATLAB-FEMM

This repository contains MATLAB codes and scripts developed for the **design and optimization of a Synchronous Permanent Magnet Motor (PMSM)** intended for **naval electric propulsion applications**.  
It was carried out as part of the course **Electrical Machines Design for Transportation Electrification (TE-1)** in the M.Sc. program *Technological & Management Advances in Intelligent Transportation Electrification Systems (TMAITES)* at Democritus University of Thrace (D.U.Th.).

We started from a base MATLAB–FEMM code, made significant modifications, and developed additional scripts to meet the **requirements and objectives** of the project.  
The optimization process explores both **single-objective (PSO)** and **multi-objective (MOPSO)** Particle Swarm Optimization approaches to identify optimal design parameters of the PMSM.

---
# Modified codes
The archives with the modified codes & the extra scripts which I create and use are:

PSO: 1) main_peirama.m, 2) ObjFunc_peirama.m, 3) PSO_peirama.m, 4) MASS_EFF_CODE.m, 5) Best_Values.m, 6) anaktisi_metavlitwn.m,
     7) CHECK_BAD_VALUES_OF_VARIABLES_CODE.m, 8) Convergence_CURVE_RECOVERY_CODE.m

MOPSO: 1) main_mopso_peirama.m, 2) MObjFunc_peirama.m, 3) MOPSO_peirama.m, 4) MASS_EFF_CODE.m, 5) Best_Values.m

# Project Objectives

- Optimize key **geometrical parameters** of the PMSM (e.g., pole width ratio, rotor/stator yoke thickness, magnet thickness).  
- Compare **single-objective optimization (PSO)** and **multi-objective optimization (MOPSO)** approaches.  
- Evaluate the trade-off between:
  - **Efficiency (η)**
  - **Motor Mass (Mmotor)**
- Provide convergence plots and Pareto fronts for detailed analysis of optimization results.  

---

# Repository Structure

- `main.m` → Main simulation script (entry point).  
- `PSO.m` → Implementation of Particle Swarm Optimization.  
- `MOPSO.m` → Implementation of Multi-Objective PSO.  
- `ObjFunc.m` → Objective function for single-objective optimization.  
- `MObjFunc.m` → Objective function for multi-objective optimization.  
- `geometry.fem` → FEMM model for electromagnetic simulations.  
- `Convergence_Curves.png` → Example convergence plot.  
- `Pareto_Fronts.png` → Example Pareto front (multi-objective results).  
- `/results/` → MATLAB `.mat` files with optimization runs.  

---

# Requirements

- **MATLAB** (R2021a or later recommended)  
- **FEMM** (Finite Element Method Magnetics) properly installed  
- **Windows OS** recommended (avoid Greek characters or spaces in file paths)  

---

# How to Run

1. Clone the repository:
   ```bash
   git clone https://github.com/NickosKokosal/Electric_Machine_Optimization_Matlab.git
   cd Electric_Machine_Optimization_Matlab
