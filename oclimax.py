"""
OCLIMAX Script
------------------------------------------
This script is a helper tool that processes molecular dynamics (MD) 
simulation data to calculate neutron scattering spectra.

Edited by: Toulik Maitra
Institution: University of California, Davis
Date: 2025
"""


import os
import argparse
from shutil import copyfile


def write_params(task, e_unit):
    """
    Creates parameters file for OCLIMAX
    Args:
    task (int): Defines approximation method where 
            0:inc approx. 1:coh+inc. 2:single-xtal Q-E. 3:single-xtal Q-Q.
    e_unit (int): Defines energy units such that 0:cm-1 1:meV 2:THz.
    #################################################################
    Return:
    out.params in the folder
    """
    with open('out.params', 'w') as f:
        f.write('## General parameters \n'
                'TASK    =         {} # 0:inc approx. 1:coh+inc. 2:single-xtal Q-E. 3:single-xtal Q-Q\n'
                'INSTR   =         0  # 0:VISION 1:indirect traj 2:direct traj 3:Q-E or Q-Q mesh\n'
                'TEMP    =      5.00  # Temperature [K]\n'
                'E_UNIT  =         0 # Energy unit [eu] (0:cm-1,1:meV,2:THz)\n'
                'OUTPUT  =         0  # 0:standard, 1:restart, 2:SPE, 3:full, 4:DOS, 5:modes\n'

                '## Additional general parameters\n'
                'MAXO    =        10  # Maximum order of excitation\n'
                'CONV    =         2  # Start convolution from order=CONV (2 or 3)\n'
                'PHASE   =         0  # Phase factor of polarization vector (0 or 1)\n'
                'MASK    =         0  # Set 1 to apply mask on Q-E/Q-Q map (INSTR=3)\n'
                'ELASTIC =  -0.10E+01 -0.10E+01  # E Q, <0:no EL,0:cal res,>0:given res\n'

                '## E parameters\n'
                'MINE    =      8.00  # Energy range (minimum) to calculate [eu]\n'
                'MAXE    =   5000.00  # Energy range (maximum) to calculate [eu]\n'
                'dE      =      1.00  # Energy bin size [eu]\n'
                'ECUT    =     8.000  # Exclude modes below this cutoff energy [eu]\n'
                'ERES    =   0.25E+01  0.50E-02  0.10E-06  # E resolution coeff\n'

                '## Q parameters\n'
                'MINQ    =      0.50  # Q range (minimum) to calculate\n'
                'MAXQ    =     20.00  # Q range (maximum) to calculate\n'
                'dQ      =      0.50  # Q bin size\n'
                'QRES    =   0.10E+00 # Q resolution coeff (INSTR=3)\n'

                '## Instrument parameters\n'
                'THETA   =  135.0  45.0  # List of scattering angles [degree]\n'
                'Ef      =     32.00  # Final energy [eu] (INSTR=1)\n'
                'Ei      =   5000.00  # Incident energy [eu] (INSTR=2)\n'
                'L1      =     11.60  # L1 [m] for DGS (INSTR=2 or 3, ERES=0)\n'
                'L2      =      2.00  # L2 [m] for DGS (INSTR=2 or 3, ERES=0)\n'
                'L3      =      3.00  # L3 [m] for DGS (INSTR=2 or 3, ERES=0)\n'
                'dt_m    =      3.91  # dt_m [us] for DGS (INSTR=2 or 3, ERES=0)\n'
                'dt_ch   =      5.95  # dt_ch [us] for DGS (INSTR=2 or 3, ERES=0)\n'
                'dL3     =      3.50  # dL3 [cm] for DGS (INSTR=2 or 3, ERES=0)\n'
                '## Single crystal parameters\n'
                'HKL     =    0.0   0.0   0.0  # HKL (TASK=2 or 3)\n'
                'Q_vec   =    0.0   0.0   1.0  # Q vector dir (TASK=2 or 3)\n'
                'Q_vec_y =    1.0   0.0   0.0  # Q vector dir y-axis (TASK=3)\n'
                'MINQ_y  =      1.00  # Q range (minimum) y-axis (TASK=3)\n'
                'MAXQ_y  =      2.00  # Q range (maximum) y-axis (TASK=3)\n'
                'dQ_y    =      0.02  # Q bin size y-axis (TASK=3)\n'

                '## Wing parameters\n'
                'WING    =         0  # Wing calculation (0:no wing,1:isotropic,2:ST tensor)\n'
                'A_ISO   =    0.0350  # Isotropic A_external for wing calculation\n'
                'W_WIDTH =     150.0  # Energy width [eu] of initial wing)\n' .format(task, e_unit))
        
def check_oclimax_success(output_dir):
    """
    Simple check if OCLIMAX completed successfully
    """
    expected_file = os.path.join(output_dir, 'out_vis_inc_5K.csv')
    return os.path.exists(expected_file) and os.path.getsize(expected_file) > 0

def oclimax(dt, params=None, task=0, e_unit=0):
    """
    Creates oclimax_[timestep] directory, writes oclimax parameters file and runs oclimax simulation for Lammps MD.
    """
    print(" ----------------------------------------------------------------- ")
    print(" Starting OCLIMAX calculation... ")
    print(" ----------------------------------------------------------------- ")
    
    # Create directory with timestep in name
    output_dir = f'oclimax_{dt}fs'
    os.makedirs(output_dir, exist_ok=True)
    copyfile('velocity.dump', f'{output_dir}/velocity.dump')
    
    # Change to oclimax directory
    os.chdir(output_dir)

    if not params:
        write_params(task, e_unit)
        params = "out.params"
    
    # Run OCLIMAX
    os.system('oclimax convert -lt velocity.dump {} 1>> ocl.out 2>> ocl.err'.format(float(dt)))
    os.system('oclimax run out.tclimax {} 1>> ocl.out 2>> ocl.err'.format(params))

    # Check if run was successful
    success = check_oclimax_success('.')
    os.chdir('..')
    
    if success:
        print(" ----------------------------------------------------------------- ")
        print(" OCLIMAX completed successfully! ")
        print(f" Results are in: {output_dir}/out_vis_inc_5K.csv")
        print(" ----------------------------------------------------------------- ")
    else:
        print(" ----------------------------------------------------------------- ")
        print(" ERROR: OCLIMAX run failed! ")
        print(f" Check {output_dir}/ocl.err for details")
        print(" ----------------------------------------------------------------- ")
    
    return success

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--timestep", type=float, default=1., help="timestep*freq in MD simulation [fs]")
    args = parser.parse_args()
    
    if oclimax(dt=args.timestep):
        exit(0)  # Success
    else:
        exit(1)  # Failure