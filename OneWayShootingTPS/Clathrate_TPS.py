import numpy as np

from openmm.app import *
from openmm import *
from openmmtools.integrators import VVVRIntegrator
from openmm import unit
from openmm.unit import *

import mdtraj as md
from mdtraj.reporters import DCDReporter

import h5py

from OneWayShootingRangeTPS import TPS_sampler
import scipy.constants as constants
import os
import re
import fcntl
import time
import shutil
import sys
import glob
from sys import stdout
sys.path.append('../MCG')
import mcg


def cv(x_p, L):

    co2_coms, water_coms, box_length = mcg.CO2ClathrateCOM(x_p, L)
    MCG_OP, _, _ = mcg.MCG_optimized(co2_coms, water_coms, box_length)

    return MCG_OP


def stateFunction(x, box_length):
    c = cv(x, box_length)
    # print(s,c)

    if c <= 10:
        return "A", c 
    elif c >= 300:
        return "B", c
    else:
        return "0", c


## this is needed to safely name data folders squentially while starting multiple scripts at the same time
def get_next_run_folder(base_dir="data", lock_file="folder_lock.lock"):

    print(f"[PID {os.getpid()}] Calling get_next_run_folder")

    lock_path = os.path.join(base_dir, lock_file)

    # Open a lock file to coordinate access
    with open(lock_path, "w") as lock:
        # Wait for and acquire exclusive lock
        fcntl.flock(lock, fcntl.LOCK_EX)

        # Double-check once locked: get next available number
        existing = [d for d in os.listdir(base_dir) if re.match(r"run_\d+$", d)]
        numbers = [int(re.search(r"run_(\d+)", d).group(1)) for d in existing]
        next_number = max(numbers) + 1 if numbers else 0
        new_folder = os.path.join(base_dir, f"run_{next_number}")
        os.makedirs(new_folder)
        print(f"[PID {os.getpid()}] Created folder: {new_folder}")


        # Release lock automatically when closing
        return next_number, new_folder




def test_openmm_SPEx_sample():

    print(f"[PID {os.getpid()}] Starting test_openmm_SPEx_sample()")


    mu_gnd = 175
    alpha_gnd = 15
    beta_gnd = 2.5


    base_path = f"data_mu175_chain1_batch_4new"
    os.makedirs(base_path, exist_ok=True)

    #base_name = 'run'
    next_number, data_folder = get_next_run_folder(base_path)
    
    
    print(f"Created new folder: {data_folder}")

    total_trials = 50


    T = 260
    P = 500

    initial_traj_folder = 'InitialTrajectories/data_mu175_chain1_batch_3_last'
    initial_traj_idx = 20

    gro_folder = "gro_files"

    # Get the trajectory input from the matching folder
    traj_folder = os.path.join(initial_traj_folder, f"run_{next_number}")
    traj_files = glob.glob(os.path.join(traj_folder, "traj_*.dcd"))

    if len(traj_files) != 1:
        raise FileNotFoundError(f"Expected exactly one traj_*.dcd file in {traj_folder}, found {len(traj_files)}")

    input_file = traj_files[0]
    print(f"Using trajectory file: {input_file}")

    #shutil.copy(input_file, os.path.join(data_folder, os.path.basename(input_file)))

    input_coordinates=os.path.join(gro_folder, 'conf.gro')
    input_topology=os.path.join(gro_folder, 'topol.top')

    positions_mdtraj = md.load(input_file, top=input_coordinates)


    print(f'Loading input coordinates from: {input_coordinates}\nLoading input topology from: {input_topology}')

    #### LOADING OF FILES AND SIMULATION #####################################################

    try:
        gro = GromacsGroFile(input_coordinates)
        top = GromacsTopFile(input_topology, periodicBoxVectors=gro.getPeriodicBoxVectors())
    except:
        print("Problem loading files...")
        print("ending script...")

    
    # Define the target pressure and temperature
    pressure = P * bar
    temperature = T * kelvin

    # Define the integration timestep and collision rate
    timestep = 2.0 * femtoseconds
    collision_rate = 1.0 / picoseconds

    # Calculate the frequency in terms of integration steps
    barostat_frequency = int((4 * picoseconds) / timestep)  # This will be 2000 steps


    system = top.createSystem(nonbondedMethod = PME, nonbondedCutoff = 1 * nanometer, constraints = HBonds)  # Using particle Mesh Ewald, setting nonbonded radius to  1nm and constraining only Hydrogen bonds
    system.addForce(MonteCarloBarostat(pressure, temperature, barostat_frequency)) # Pressure, Temperature, frequency
    system.addForce(CMMotionRemover())
    integrator = VVVRIntegrator(temperature, collision_rate, timestep) 
    simulation_eq = Simulation(top.topology, system, integrator)


   # initial_traj = md.load_hdf5(input_file, stride = 1)
    initial_traj = md.load(input_file, top=input_coordinates)


    # Allocate velocity array
    velocities = np.zeros((initial_traj.n_frames, initial_traj.n_atoms, 3))
    
    for i in range(initial_traj.n_frames):
        # Set positions for this frame
        simulation_eq.context.setPositions(initial_traj.xyz[i])
        
        # Resample velocities at the desired temperature
        simulation_eq.context.setVelocitiesToTemperature(temperature)
    
        # Get the state and extract velocities
        state = simulation_eq.context.getState(getVelocities=True)
        vel = state.getVelocities(asNumpy=True)  # Quantity array (nm/ps)
    
        # Convert to ndarray in units of m/s or whatever you want
        velocities[i] = vel.value_in_unit(nanometers / picoseconds)
    
    # Assign to traj (note: mdtraj doesn't track units)
    initial_traj.velocities = velocities

    
    simulation_eq.context.setPositions(initial_traj.xyz[initial_traj_idx])

    #with h5py.File(input_file, "r") as f:
    #    velocities = f["velocities"][:]

    #initial_traj.velocities = velocities


    #simulation_eq.context.setVelocities(velocities[initial_traj_idx] * nanometers/picosecond)  # Ensure units
    simulation_eq.context.setVelocitiesToTemperature(temperature)
    simulation_eq.context.setPeriodicBoxVectors(*initial_traj.unitcell_vectors[initial_traj_idx])

    state = simulation_eq.context.getState()
    current_vectors = state.getPeriodicBoxVectors()

    simulation_eq.reporters.append(StateDataReporter(stdout, 1e5, step=True, potentialEnergy=True, temperature=True, volume = True))


    sampler = TPS_sampler(simulation_eq, data_folder, mu_gnd, alpha_gnd, beta_gnd)

    sampler.sample(initial_traj, stateFunction, cv, total_trials, stride=1e5, maxPathLength=3e8) #mpl in unit steps


test_openmm_SPEx_sample()
