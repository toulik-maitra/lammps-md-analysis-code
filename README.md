# lammps-md-analysis-code

## OCLIMAX Script (oclimax.py)

The `oclimax.py` script is a helper tool that processes molecular dynamics (MD) simulation data to calculate neutron scattering spectra. Here's what it does in simple terms:

### What is it for?
- It takes the `velocity.dump` data from your molecular dynamics simulation 
- It converts this data into a format that can be used to calculate how neutrons would scatter off your simulated material

### Main Features:
1. **Easy Setup**: Creates all necessary files and directories automatically
2. **Flexible Options**: 
   - Can handle different types of calculations (incoherent approximation, coherent+incoherent, single crystal analysis)
   - Supports different energy units (cm⁻¹, meV, THz)
   
### How to Use:
0. Make sure that the 'elements' column in the `velocity.dump` file has the desired elements you are looking for. 
1. Make sure you have a `velocity.dump` file from your LAMMPS simulation in your working directory
2. Run the script using Python:
   ```bash
   python oclimax.py --timestep YOUR_TIMESTEP
   ```
   (Replace YOUR_TIMESTEP with the actual timestep used in your MD simulation in femtoseconds)

### What You Get:
- Creates a new directory named `oclimax_[timestep]fs`
- Generates two main output files:
  - `out_vis_inc_5K.csv`: Contains the processed neutron scattering data
  - `out.tclimax`: Contains technical details of the calculation

### Note:
This script requires the OCLIMAX program to be installed on your system to work properly.