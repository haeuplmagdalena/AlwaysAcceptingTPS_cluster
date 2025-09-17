from openmm.app import *
from openmm import *
from openmm.unit import *
from openmmtools import integrators
from sys import stdout
import numpy as np
import matplotlib.pyplot as plt
import mdtraj as md
import os
import argparse
from mdtraj.reporters import HDF5Reporter
import shutil


parser = argparse.ArgumentParser(description='folder inputs')
parser.add_argument("--input_file", type=str, required=True)
parser.add_argument('--N', dest='N', required = True)
parser.add_argument('--P', dest='P', default = 500)
parser.add_argument('--T', dest='T', default = 260)
parser.add_argument('--N_out', dest = 'N_out', default = 1)
args = parser.parse_args()


#### VARIABLE AND PATH SETUP ###############################################################

N = int(args.N)
N_out = int(args.N_out)
T = int(args.T)
P = int(args.P)
input_file = str(args.input_file)

data_folder = os.path.abspath(os.path.join(input_file, os.pardir))

print(f'Saving to folder: {data_folder}')

print(f'PRESSURE    = {P}\nTEMPERATURE = {T}')

input_coordinates=f'/leonardo_work/L-AUT_Coretti/CODE/VACF/gro_files/conf.gro'
input_topology=f'/leonardo_work/L-AUT_Coretti/CODE/VACF/gro_files/topol.top'

topology_mdtraj = md.load(input_coordinates)

positions_mdtraj = md.load(input_file, top=topology_mdtraj)


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


system = top.createSystem(nonbondedMethod = PME, nonbondedCutoff = 1 * nanometer, constraints = HBonds)  # Using particle Mesh Ewald, setting nonbonded radius to  1nm and constrianing only Hydrogen bonds
system.addForce(MonteCarloBarostat(pressure, temperature, barostat_frequency)) # Pressure, Temperature, frequency
system.addForce(CMMotionRemover())
integrator = integrators.VVVRIntegrator(temperature, collision_rate, timestep) # Temperature, collision rate, timestep
platform = Platform.getPlatformByName('CUDA')
simulation = Simulation(top.topology, system, integrator)

simulation.context.setPositions(positions_mdtraj.xyz[0])
simulation.context.setVelocitiesToTemperature(temperature)

hdf5_reporter = HDF5Reporter(f'{data_folder}/trajectory.h5', N_out, velocities=True)

# Append the reporter to the simulation
simulation.reporters.append(hdf5_reporter)
simulation.reporters.append(StateDataReporter(f'{data_folder}/energy.csv', N_out, step = True, totalEnergy=True,
        potentialEnergy=True, kineticEnergy=True, temperature=True, volume=True))
simulation.reporters.append(StateDataReporter(stdout, 100, step=True,
        potentialEnergy=True, temperature=True, volume = True))


simulation.step(N)

final_state = simulation.context.getState(getPositions=True, getVelocities=True)
with open(f'{data_folder}/final_state_file', 'w') as f:
    f.write(XmlSerializer.serialize(final_state))

positions = simulation.context.getState(getPositions=True).getPositions()
PDBFile.writeFile(simulation.topology, positions, open(f'{data_folder}/pdb_final.pdb', 'w'))

