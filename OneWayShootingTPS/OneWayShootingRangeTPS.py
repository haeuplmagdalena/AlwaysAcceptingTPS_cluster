import os
import numpy as np
import scipy.constants as constants
from scipy.special import logsumexp


from openmm.app import *
from openmm import *
from openmm import unit

from mdtraj.reporters import DCDReporter
import mdtraj

class TPS_sampler:
    """
    A class to perform path sampling using openmm molecular dynamics simulations.
    This class implements shooting moves and path acceptance/rejection.
    """
    
    def __init__(self, simulation_eq: Simulation, working_folder: str, mu: int, alpha: int, beta: int) -> None:
        """
        Initialize the sampler with equilibrium and metadynamics simulations,
        working folder, and bias force group.

        Args:
            simulation_eq (Simulation): The equilibrium simulation object.
            simulation_metad (Simulation): The metadynamics simulation object.
            working_folder (str): The directory to save output files.
        """
        self.mu_gnd = mu
        self.alpha_gnd = alpha
        self.beta_gnd = beta


        self.working_folder = working_folder
        self.simulation_eq = simulation_eq
        self.box_length = simulation_eq.topology.getPeriodicBoxVectors()[0][0]._value

        self.temperature = self.simulation_eq.integrator.getTemperature()
        self.beta = 1/(constants.R / 1000 * self.temperature._value)
        
        self.oldTrajectory = None  # Stores the previous trajectory
        self.trajectory_weights = []  # Weights for each trajectory
        self.frag_pn  = None  # New fragment of the trajectory
        self.frag_po = None  # Old fragment of the trajectory
        self.vel_po = None
        self.vel_pn = None

        self.velocities_pn = None

        self.sp_idx_old = None
        self.sp_idx_new = None
        self.cv_array_old = []
        self.cv_array_new = []

        self.target_state = None
    
        self.currentIndex = 1  # Current trajectory index

        self.acceptedPaths = 0  # Counter for accepted paths
        self.successful_SP_indices = []


    def simulate_to_state(self, x: np.ndarray, v: np.ndarray, sp_box_vectors: np.ndarray, 
                          stateFunction: callable, cv_list: list, vel_list: list, filename: str, 
                          maxSteps: int) -> str:
        """
        Simulate the system until it reaches a defined state (A or B) or maxSteps is reached.

        Args:
            x (np.ndarray): Initial positions.
            v (np.ndarray): Initial velocities.
            stateFunction (callable): Function to determine the state of the system.
            filename (str): File to save the trajectory.
            maxSteps (int): Maximum number of steps to simulate.
            self.stride (int): Frequency to save frames.

        Returns:
            str: The state ("A", "B", or "0" if no state is reached).
        """

        self.simulation_eq.context.setPositions(x)
        self.simulation_eq.context.setStepCount(0)
        self.simulation_eq.context.setTime(0)
        self.simulation_eq.context.setVelocities(v)

        # When using non-cubic box, change the following
        self.simulation_eq.context.setPeriodicBoxVectors(*sp_box_vectors)  #mdtraj gives an (frames,3,3) array and opnemm takes (a,b,c) 

        self.simulation_eq.reporters.clear()
        self.simulation_eq.reporters.append(DCDReporter(filename, self.stride))      #since the system is using PBC, this defaults to enforcePeriodicBox = True
        state = self.simulation_eq.context.getState(getPositions=True)

        self.simulation_eq.reporters[0].report(self.simulation_eq, state)

        steps = 0
        while(steps < maxSteps):        #this 'steps' variable lines up exactly with the openmm steps

            self.simulation_eq.step(self.stride)
            steps += self.stride
            
            state = self.simulation_eq.context.getState(
                getPositions=True,
                getVelocities=True
            )
            
            positions = state.getPositions(asNumpy=True)
            boxlength = state.getPeriodicBoxVectors()
            velocities = state.getVelocities(asNumpy=True) 
            vel_list.append(velocities)

            # making sure the positions and the box length are in nanometers
            positions_value = positions.value_in_unit(unit.nanometer).astype(np.float64)    # this used to be float32, changed for testing
            boxlength_value = boxlength[0][0].value_in_unit(unit.nanometer)
            

            inState, cv = stateFunction(positions_value, boxlength_value)
            cv_list.append((int(steps), cv))  # Save step count and CV as a tuple, appending onto the value at the SP


            # Testing CV discrepancy
            with open(os.path.join(self.working_folder, f"log_test.txt"), "a") as file:
                file.writelines(f"\nSTEPS: {steps}")
                file.writelines(f"\ncurrent_boxlength: {boxlength}")
                file.writelines(f"\CV: {cv}")

            if inState ==  "A" or inState == "B":
                self.simulation_eq.reporters.clear()
                return inState, cv_list, vel_list

        self.simulation_eq.reporters.clear()
        return "0", cv_list, vel_list


    def perform_shooting_move(self, x: np.ndarray, v: np.ndarray, stateFunction: callable, direction: int, sp_box_vectors: np.ndarray,  cv_po: np.ndarray,
                              maxSteps: int) -> tuple:
        """
        Perform a shooting move from given coordinates/velocities and return the resulting state.
        
        Note: This function updates the following global class variables:
            - self.sp_idx_new (int): Index of the shooting point in the trajectory
            - self.cv_array_new (np.ndarray): Updated collective variable array for the new trajectory
            - self.velocities_pn (list): Velocities of the newly generated trajectory segment
        
        Args:
            x (np.ndarray): Position of shooting point (shape: [n_atoms, 3])
            v (np.ndarray): Velocity of shooting point (shape: [n_atoms, 3])
            stateFunction (callable): Function to determine state ('A', 'B', or '0') from coordinates
            direction (int): Shooting direction (1=forward, -1=backward)
            sp_box_vectors (np.ndarray) : The box vectors at the shooting point
            cv_po (np.ndarray): Collective variable values for the "old" trajectory segment
            maxSteps (int): Maximum simulation steps for path generation
        
        Returns:
            state_pn (str): Resulting state after shooting ('A', 'B', or '0')
        """

        # Check initial state of shooting point
        init_state, cv_init = stateFunction(x.astype(np.float32), sp_box_vectors[0][0])

        # Initialising lists with shooting point data
        cv_pn = [[0, cv_init]]
        velocities_pn = [v]

        with open(os.path.join(self.working_folder, f"log_test.txt"), "a") as file:
            file.writelines(f"\nStarting from SP at {cv_init}")

        # If shooting point is already in target state, record nothing, since shooting point is saved in old traj fragment
        if(init_state == self.target_state):

            rep = DCDReporter(os.path.join(self.working_folder, "current_traj.dcd"), self.stride)
            rep.report(self.simulation_eq, self.simulation_eq.context.getState(getPositions=True))
            state_pn = init_state
        
        else:
            # Simulate until reaching A/B or maxSteps
            state_pn, cv_pn, velocities_pn = self.simulate_to_state(x, v, sp_box_vectors, stateFunction, cv_pn, velocities_pn, os.path.join(self.working_folder, "current_traj.dcd"), maxSteps)  

        cv_pn = np.array(cv_pn)

        self.velocities_pn = np.array(velocities_pn)

        if direction == 1:
            # Forward shooting: append new segment after old segment
            self.sp_idx_new = self.sp_idx_old   # sp index does not change in this case
            cv_pn[:, 0] += self.sp_idx_new # shift step by last step of old part
            cv_array_new = np.vstack(( cv_po, cv_pn[1:, :] )) 

        else:
            # Backward shooting: prepend reversed new segment before old segment
            shift = int(cv_pn[-1, 0] - cv_po[0, 0])
            cv_po[:, 0] += shift
            cv_pn[:, 1] = cv_pn[::-1, 1]
            self.sp_idx_new = int(cv_pn[-1, 0])
            cv_array_new = np.vstack(( cv_pn, cv_po[1:, :] )) 

        return state_pn, cv_array_new

    def pathAcceptance(self, state_pn: str, p_sel_old: float) -> bool:
        """
        Metropolis criterion for path acceptance. Returns True ONLY when:
        1. New path reaches target state (self.target_state)
        2. EITHER:
            a) p_sel_new > p_sel_old (always accept), OR
            b) Random number < p_sel_new/p_sel_old (probabilistic accept)
        
        Args:
            state_pn: State of new path segment ('A', 'B', or '0')
            p_sel_old: Shooting point probability from original path
        
        Returns:
            bool: Acceptance decision (True/False)

        """

        accepted = False

        # Only consider acceptance if new path reaches target state
        if state_pn == self.target_state: 

            # Calculate selection probabilities for new path
            self.p_sel_array_new = self.calculate_selection_probabilities(self.cv_array_new[:, 1], self.mu_gnd, self.alpha_gnd, self.beta_gnd)
            
            # Get selection probability at new shooting point
            p_sel_new = self.p_sel_array_new[int(self.sp_idx_new // self.stride)]

            # Metropolis acceptance criteria
            if p_sel_new > p_sel_old:
                accepted = True     # Always accept if probability increased

            elif np.random.random() < p_sel_new / p_sel_old:
                # Note: Correctly uses half-open interval [0,1) via np.random.random
                accepted = True

        return accepted


    def acceptPath(self, direction: int, trial: int) -> None:
        """
        Accepts the newly generated path and updates all trajectory references.
        
        Performs these key operations:
        1. Reconstructs the full trajectory by joining old and new fragments
        2. Updates velocities to match the new trajectory
        3. Saves the accepted path and associated metadata
        4. Increments counters and weights
        
        Updates these class variables:
            * self.oldTrajectory: The newly accepted full trajectory
            * self.frag_po: Old trajectory segment (pre/post shooting point)
            * self.frag_pn: New trajectory segment
            * self.vel_po: Velocities of old segment
            * self.vel_pn: Velocities of new segment

            * self.trajectory_weights: Appends weight=1 for new path
            * self.currentIndex: Increments trajectory version
            * self.acceptedPaths: Increments acceptance counter
        
        File Operations:
        - Saves these files in working_folder:
            1. New trajectory (traj_{n+1}.dcd)
            2. Shooting point info (sp_{n+1}.txt)
            3. Collective variables (cv_{n+1}.txt)
        """

        full_old_traj = mdtraj.Trajectory.load(os.path.join(self.working_folder, f"traj_{self.currentIndex}.dcd"), top=self.oldTrajectory.topology)
            
        if direction == 1:
            self.frag_po = full_old_traj[:self.sp_idx_old // self.stride + 1]
            self.frag_pn = mdtraj.Trajectory.load(os.path.join(self.working_folder, "current_traj.dcd"), top=self.oldTrajectory.topology)
            self.vel_po = self.oldTrajectory.velocities[:self.sp_idx_old // self.stride + 1]
            self.vel_pn = self.velocities_pn

            self.oldTrajectory = mdtraj.join([self.frag_po, self.frag_pn[1:]])   
            self.oldTrajectory.velocities = np.vstack((self.vel_po, self.vel_pn[1:]))
            
        else:
            self.frag_po = full_old_traj[self.sp_idx_old // self.stride:]
            self.frag_pn = mdtraj.Trajectory.load(os.path.join(self.working_folder, "current_traj.dcd"), top=self.oldTrajectory.topology)[::-1]
            self.vel_po = self.oldTrajectory.velocities[self.sp_idx_old // self.stride:]
            self.vel_pn = self.velocities_pn[::-1] * -1

            self.oldTrajectory = mdtraj.join([self.frag_pn, self.frag_po[1:]])     
            self.oldTrajectory.velocities = np.vstack((self.vel_pn, self.vel_po[1:]))


        np.save( os.path.join(self.working_folder, f'velocities_{self.currentIndex + 1}.npy'), self.oldTrajectory.velocities )
        
        self.oldTrajectory.save_dcd(os.path.join(self.working_folder, f"traj_{self.currentIndex + 1}.dcd"))
        

        with open(os.path.join(self.working_folder, f'log.txt'), "a") as spfile:  
            row = "{:<20} {:<20} {:<20} {:<20} {:<20} {:<20} {:<20}".format(trial, self.currentIndex + 1, "True", int(self.sp_idx_old), int(self.sp_idx_new), len(self.oldTrajectory), direction)
            spfile.write(row + "\n")


        np.savetxt(os.path.join(self.working_folder, f"cv_{self.currentIndex + 1}.txt"), self.cv_array_new, fmt='%d %.4f', header="frame CV")

        self.trajectory_weights.append(1)
        self.currentIndex += 1
        self.acceptedPaths += 1
        

    def calculate_selection_probabilities(self, X_cv: np.ndarray, mu_gnd: float, alpha_gnd: float, beta_gnd: float):
        """
        Calculates normalized selection probabilities for shooting points based on 
        collective variable values using a generalized Gaussian distribution.
        
        Args:
            X_cv (np.ndarray): Array of collective variable values
            mu_gnd (float): Center parameter of the distribution
            alpha_gnd (float): Scale parameter (controls width)
            beta_gnd (float): Shape parameter (controls tails)
            
        Returns:
            np.ndarray: Normalized probabilities (sum to 1.0)
            
        Note:
            Implements the generalized normal distribution formula:
            p(x) ~ exp( -(|x-μ|/α)^β )
            where μ=mu_gnd, α=alpha_gnd, β=beta_gnd
        """
        num = np.exp( -( np.abs(X_cv - mu_gnd) / alpha_gnd )**beta_gnd )
        p_sel_array = num / np.sum(num)

        return p_sel_array
        

    def sample(self, initialTrajectory: mdtraj.Trajectory, stateFunction: callable, cv: callable, 
            total_trials: int, stride: int, maxPathLength: int) -> None:
        """
        Perform transition path sampling using shooting moves to explore pathways between states.

        This method implements a shooting move algorithm to sample transition paths between 
        two states (A and B). It starts with an initial trajectory and iteratively modifies it
        by selecting shooting points, perturbing velocities, and accepting/rejecting new paths
        based on a Metropolis criterion.

        Parameters:
        -----------
        initialTrajectory : mdtraj.Trajectory
            The initial trajectory to start sampling from. Should connect states A and B.
        stateFunction : callable
            Function that takes coordinates and returns the state ('A', 'B', or None).
        cv : callable
            Collective variable function that takes coordinates and returns a CV value.
        total_trials : int
            Total number of shooting move trials to perform.
        stride : int
            Stride for saving frames and analyzing the trajectory.
        maxPathLength : int
            Maximum allowed length for generated paths.

        Returns:
        --------
        None
            Results are saved to files in the working folder, including:
            - Trajectory files (.dcd)
            - Collective variable files (.txt)
            - Log file with acceptance statistics
            - Trajectory weights and successful shooting point indices
        """

        # Store initial trajectory and its velocities
        self.oldTrajectory = initialTrajectory
        self.stride = int(stride)
        self.currentIndex = 0
        self.trajectory_weights = [1]  # Weight of initial trajectory is 1

        # Save initial trajectory
        self.oldTrajectory.save_dcd(os.path.join(self.working_folder, f"traj_{self.currentIndex}.dcd"))
        
        # Calculate collective variables for initial trajectory
        cv_list = []

        for steps in range(0, len(self.oldTrajectory.xyz)):
            L  = self.oldTrajectory.unitcell_vectors[steps][0][0]
            cv_list.append((self.stride*steps, cv(self.oldTrajectory.xyz[steps], L=L)))


        # Store and save CV values
        self.cv_array_old = np.array(cv_list)

        np.savetxt(os.path.join(self.working_folder, f"cv_{self.currentIndex}.txt"), 
                self.cv_array_old, fmt='%d %.4f', header="frame CV")
        
        # Calculate selection probabilities for shooting points
        p_sel_array_old = self.calculate_selection_probabilities(
            self.cv_array_old[:, 1], self.mu_gnd, self.alpha_gnd, self.alpha_gnd)
        
        # Initialize counters
        self.acceptedPaths = 0

        # Initialize log file
        with open(os.path.join(self.working_folder, f'log.txt'), "a") as spfile:  # Open file inside the function
            header = "{:<20} {:<20} {:<20} {:<20} {:<20} {:<20} {:<20}".format("Trial", "Index", "Successful", "SP frame on old traj", "SP frame on new traj", "Full Path Length", "Shooting Direction")
            spfile.write(header + "\n")
    

        # Main sampling loop
        for trial in range(1, total_trials):


            # Test logging
            with open(os.path.join(self.working_folder, f"log_test.txt"), "a") as file:
                file.writelines(f"\n\tTrial:\t{trial}/{total_trials}")

            print(f"\r{trial}/{total_trials}", end="", flush=True)
 
            # Store current state securely
            current_cv_snapshot = self.cv_array_old.copy()  # Explicit copy
            current_p_sel_snapshot = p_sel_array_old.copy()

            # 1. Select shooting point from current trajectory
            self.sp_idx_old = int(np.random.choice(self.cv_array_old[:, 0], p=p_sel_array_old))
            shooting_point = self.oldTrajectory.xyz[self.sp_idx_old // self.stride] 
            sp_box_vectors = self.oldTrajectory.unitcell_vectors[self.sp_idx_old // self.stride]
            shooting_point_vel_old = self.oldTrajectory.velocities[self.sp_idx_old // self.stride]
            p_sel_old = p_sel_array_old[self.sp_idx_old // self.stride]


            # 2. Choose shooting direction (forward or backward)
            rdm_idx = np.random.choice([0,1])
            shooting_direction = np.array([-1,1])[rdm_idx]
            shooting_point_vel = shooting_point_vel_old * shooting_direction

            # 3. Prepare old part of trajectory for acceptance calculation
            if shooting_direction == -1:
                # For backward shooting: old part is from SP to end
                self.target_state = "A"
                cv_po = self.cv_array_old[self.sp_idx_old // self.stride:]      
            else:
                # For forward shooting: old part is from start to SP
                self.target_state = "B"
                cv_po = self.cv_array_old[:self.sp_idx_old // self.stride + 1]

            newMaxPathLength = maxPathLength - len(cv_po)

            # 4. Perform shooting move and generate new trajectory segment
            state_pn, cv_array_new = self.perform_shooting_move(
                shooting_point, shooting_point_vel, stateFunction, 
                shooting_direction, sp_box_vectors, cv_po, newMaxPathLength)
        
            self.cv_array_new = cv_array_new.copy()

            # 5. Decide whether to accept the new path
            path_accepted = self.pathAcceptance(state_pn, p_sel_old)

            
            if path_accepted:               
                # Accept the new path and update all references
                self.acceptPath(shooting_direction, trial)
                self.cv_array_old = self.cv_array_new.copy()
                p_sel_array_old = self.p_sel_array_new

            else:

                # Write logdata
                with open(os.path.join(self.working_folder, f'log.txt'), "a") as spfile:  # Open file inside the function
                    row = "{:<20} {:<20} {:<20} {:<20} {:<20} {:<20} {:<20}".format(trial, self.currentIndex + 1, "False", int(self.sp_idx_old), int(self.sp_idx_new), len(self.cv_array_new), shooting_direction)
                    spfile.write(row + "\n")

                # Saving the generated but rejected snippets
                old_path = os.path.join(self.working_folder, "current_traj.dcd")
                new_path = os.path.join(self.working_folder, f"rejected_traj_trial_{trial}.dcd")
                os.rename(old_path, new_path)

                # Resetting cv and psel arrays to snapshots
                self.cv_array_old = current_cv_snapshot
                p_sel_array_old = current_p_sel_snapshot

                # Increment weight of current trajectory
                self.trajectory_weights[-1] += 1

                # Save CV file
                np.savetxt(os.path.join(self.working_folder, f"rejected_cv_trial_{trial}.txt"), self.cv_array_new, fmt='%d %.4f', header="frame CV")



            # Save final results
            np.save(os.path.join(self.working_folder, f"trajectory_weights.npy"), self.trajectory_weights)
            np.save(os.path.join(self.working_folder, f"successful_SP_indices.npy"), self.successful_SP_indices)
