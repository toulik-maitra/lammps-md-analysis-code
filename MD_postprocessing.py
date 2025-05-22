"""
Molecular Dynamics Trajectory Analysis Script
------------------------------------------
This script performs various analyses on MD trajectory data including RMSD, 
RMSF, radius of gyration, MSD, and RDF calculations.

Author: Toulik Maitra
Institution: University of California, Davis
Date: 2025
"""

import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis.rms import RMSD, RMSF
from MDAnalysis.analysis.rdf import InterRDF
from MDAnalysis.analysis.msd import EinsteinMSD
import matplotlib.pyplot as plt
from tqdm import tqdm
from pathlib import Path

# Global plotting parameters
plt.rcParams.update({
    'figure.figsize': (10, 6),
    'font.size': 12,
    'axes.labelsize': 12,
    'axes.titlesize': 14,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'lines.linewidth': 2,
    'axes.grid': True,
    'grid.alpha': 0.3
})

# Analysis parameters
DATA_FILE = 'system.data'
TRAJ_FILE = 'equil20.lammpstrj'
TRAJ_FORMAT = 'LAMMPSDUMP'
OUTPUT_DIR = 'analysis_output'
PLOT_DPI = 300

def save_figure(name):
    """Save figure with consistent formatting"""
    plt.tight_layout()
    plt.savefig(f'{OUTPUT_DIR}/{name}.png', dpi=PLOT_DPI, bbox_inches='tight')
    plt.close()

def calculate_rmsd(universe):
    """Calculate Root Mean Square Deviation"""
    print("\nCalculating RMSD...")
    rmsd = RMSD(universe, universe, select='all').run()
    plt.figure()
    plt.plot(rmsd.rmsd[:,1], rmsd.rmsd[:,2])
    plt.xlabel('Time (ps)')
    plt.ylabel('RMSD (Å)')
    plt.title('Root Mean Square Deviation vs Time')
    save_figure('rmsd')

def calculate_rmsf(universe):
    """Calculate Root Mean Square Fluctuation"""
    print("Calculating RMSF...")
    rmsf = RMSF(universe.select_atoms('all')).run()
    plt.figure()
    plt.plot(range(len(rmsf.results.rmsf)), rmsf.results.rmsf)
    plt.xlabel('Atom Index')
    plt.ylabel('RMSF (Å)')
    plt.title('Root Mean Square Fluctuation by Atom')
    save_figure('rmsf')

def calculate_radius_of_gyration(universe):
    """Calculate Radius of Gyration"""
    print("Calculating Radius of Gyration...")
    masses = universe.atoms.masses.copy()
    mass_sum = masses.sum()
    if mass_sum == 0:
        masses = np.ones_like(masses)
        mass_sum = masses.sum()

    rg_times, rg_vals = [], []
    for ts in tqdm(universe.trajectory, desc="Processing frames"):
        pos = universe.atoms.positions
        com = (masses[:, None] * pos).sum(axis=0) / mass_sum
        dpos = pos - com
        rg2 = (masses[:, None] * (dpos**2)).sum() / mass_sum
        rg_times.append(ts.time)
        rg_vals.append(np.sqrt(rg2))

    plt.figure()
    plt.plot(rg_times, rg_vals)
    plt.xlabel('Time (ps)')
    plt.ylabel('Rg (Å)')
    plt.title('Radius of Gyration vs Time')
    save_figure('rg')

def calculate_msd(universe):
    """Calculate Mean Square Displacement"""
    print("Calculating Mean Square Displacement...")
    msd = EinsteinMSD(universe, select='all', msd_type='xyz', fft=True).run()
    lag_times = np.arange(len(msd.results.timeseries)) * universe.trajectory.dt
    plt.figure()
    plt.plot(lag_times, msd.results.timeseries)
    plt.xlabel('Lag time (ps)')
    plt.ylabel('MSD (Å²)')
    plt.title('Mean Square Displacement vs Lag Time')
    save_figure('msd')

def calculate_rdf(universe):
    """Calculate Radial Distribution Function"""
    print("Calculating Radial Distribution Function...")
    try:
        rdf = InterRDF(universe.select_atoms('all'), universe.select_atoms('all')).run()
        plt.figure()
        plt.plot(rdf.bins, rdf.rdf)
        plt.xlabel('Distance (Å)')
        plt.ylabel('g(r)')
        plt.title('Radial Distribution Function')
        save_figure('rdf')
    except Exception as e:
        print(f"Note: RDF analysis skipped - {str(e)}")

def main():
    """Main execution function"""
    # Create output directory
    Path(OUTPUT_DIR).mkdir(exist_ok=True)
    
    # Load trajectory
    print("Loading trajectory...")
    universe = mda.Universe(DATA_FILE, TRAJ_FILE, format=TRAJ_FORMAT)
    
    # Perform analyses
    calculate_rmsd(universe)
    calculate_rmsf(universe)
    calculate_radius_of_gyration(universe)
    calculate_msd(universe)
    calculate_rdf(universe)
    
    print(f"\nAnalysis complete. Results saved in '{OUTPUT_DIR}' directory.")

if __name__ == "__main__":
    main()
