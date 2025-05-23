# LAMMPS MD Analysis Tools

## Table of Contents
1. [OCLIMAX Script (oclimax.py)](#oclimax-script-oclimaxpy)
   - [Purpose](#purpose)
   - [Main Features](#main-features)
   - [How to Use](#how-to-use)
   - [Output Files](#output-files)
   - [Requirements](#requirements)

2. [MD Post-processing Script (MD_postprocessing.py)](#md-post-processing-script-md_postprocessingpy)
   - [Purpose](#purpose-1)
   - [Analysis Features](#analysis-features)
   - [Output Plots](#output-plots)
   - [Requirements](#requirements-1)

## OCLIMAX Script (oclimax.py)

### Purpose
The `oclimax.py` script processes molecular dynamics (MD) simulation data to calculate neutron scattering spectra. It takes the `velocity.dump` data from your MD simulation and converts it into a format suitable for neutron scattering analysis.

### Main Features
1. **Easy Setup**: Creates all necessary files and directories automatically
2. **Flexible Options**: 
   - Supports multiple calculation types:
     - Incoherent approximation
     - Coherent + incoherent analysis
     - Single crystal analysis
   - Multiple energy unit options (cm⁻¹, meV, THz)
   
### How to Use
1. Ensure the `elements` column in your `velocity.dump` file contains the correct elements
2. Place the `velocity.dump` file in your working directory
3. Run the script:
   ```bash
   python3 oclimax.py -t YOUR_TIMESTEP
   ```
   Replace YOUR_TIMESTEP = MD simulation timestep in femtoseconds * Dumping Frequency

### Output Files
- Creates directory: `oclimax_[timestep]fs`
- Generates:
  - `out_vis_inc_5K.csv`: Processed neutron scattering data
  - `out.tclimax`: Technical calculation details

### Requirements
- OCLIMAX program installed on your system
- Python 3.x

## MD Post-processing Script (MD_postprocessing.py)

### Purpose
A comprehensive analysis tool for molecular dynamics trajectories, providing various statistical and structural analyses of your simulation data.

### Analysis Features
1. **Root Mean Square Deviation (RMSD)**
   - Tracks structural changes over time
   - Results plotted in nanoseconds

2. **Root Mean Square Fluctuation (RMSF)**
   - Analyzes atomic mobility
   - Shows fluctuation by atom index

3. **Radius of Gyration**
   - Measures system size changes
   - Accounts for mass distribution

4. **Mean Square Displacement (MSD)**
   - Calculates atomic diffusion
   - Uses FFT for efficient computation

5. **Radial Distribution Function (RDF)**
   - Analyzes atomic structure
   - Optimized for large trajectories

### Output Plots
All plots are automatically saved in the `analysis_output` directory:
- `rmsd.png`: RMSD vs. Time
- `rmsf.png`: RMSF by Atom
- `rg.png`: Radius of Gyration vs. Time
- `msd.png`: Mean Square Displacement
- `rdf.png`: Radial Distribution Function

### Requirements
- Python 3.x
- Required packages:
  - MDAnalysis
  - NumPy
  - Matplotlib
  - tqdm

### Input Files Needed
- `system.data`: LAMMPS data file
- `<any traj file>.lammpstrj`: LAMMPS trajectory file