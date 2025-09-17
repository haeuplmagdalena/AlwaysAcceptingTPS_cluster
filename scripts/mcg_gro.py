import numpy as np
import commons as params
from numba import jit
import mdtraj as md
import networkx as nx
import re
from scipy import linalg


def COM_calculation_gro(file):
    """Parse a .gro file and calculate the center of mass for CO2 (UNL) and water (SOL)."""
    C_M = params.C_M  # Carbon atomic mass
    O_M = params.O_M  # Oxygen atomic mass
    H_M = params.H_M  # Hydrogen atomic mass

    CAR_positions = []
    WATER_positions = []

    with open(file, 'r') as f:
        lines = f.readlines()
    
    atom_lines = lines[2:-1]
    box_line = lines[-1].strip()
    box_length = float(box_line.split()[0])  # Assuming a cubic box
    
    for line in atom_lines:
        residue_name = line[5:10].strip()
        atom_name = line[10:15].strip()
        x = float(line[20:28].strip())
        y = float(line[28:36].strip())
        z = float(line[36:44].strip())

        # Collect CO2 atoms (UNL residue)
        if residue_name == "UNL":
            CAR_positions.append((atom_name, np.array([x, y, z])))
        
        # Collect water atoms (SOL residue)
        elif residue_name == "SOL":
            WATER_positions.append((atom_name, np.array([x, y, z])))
    
    # Calculate COM for CO2 (UNL)
    CAR_COM = []
    for i in range(0, len(CAR_positions), 3):  # Groups of 3 atoms
        co_pos = CAR_positions[i][1]
        o1_pos = CAR_positions[i + 1][1]
        o2_pos = CAR_positions[i + 2][1]
        com = (co_pos * C_M + (o1_pos + o2_pos) * O_M) / (C_M + 2 * O_M)
        CAR_COM.append(com)
    
    # Calculate COM for water (SOL)
    WATER_COM = []
    for i in range(0, len(WATER_positions), 4):  # Groups of 4 atoms
        o_pos = WATER_positions[i][1]
        h1_pos = WATER_positions[i + 1][1]
        h2_pos = WATER_positions[i + 2][1]
        com = (o_pos * O_M + (h1_pos + h2_pos) * H_M) / (O_M + 2 * H_M)
        WATER_COM.append(com)
    
    return np.array(CAR_COM), np.array(WATER_COM), box_length

def COM_calculation_gro_multiframe(file):
    """Parse a multi-frame .gro file and calculate the COM and MCG_OP for each frame."""
    C_M = params.C_M  # Carbon atomic mass
    O_M = params.O_M  # Oxygen atomic mass
    H_M = params.H_M  # Hydrogen atomic mass

    # Data structures to store results
    all_car_com = []
    all_water_com = []
    all_box_lengths = []
    all_steps = []  # Store the 'step' value from the top line of each frame

    with open(file, 'r') as f:
        lines = f.readlines()

    # Identify the structure of the file
    num_atoms = int(lines[1].strip())  # Number of atoms (constant across frames)
    frame_size = num_atoms + 3        # 2 header lines + atom lines + 1 box line
    num_frames = len(lines) // frame_size

    for frame_idx in range(num_frames):
        start_line = frame_idx * frame_size
        end_line = start_line + frame_size
        frame_lines = lines[start_line:end_line]

        # Extract 'step' number from the first line of the frame
        first_line = frame_lines[0].strip()
        step_match = re.search(r'step=\s*(\d+)', first_line)  # Regex to extract number after 'step='
        if step_match:
            step_value = int(step_match.group(1))  # Extracted step value
            all_steps.append(step_value)
        else:
            print(f"Warning: Could not find step value in frame {frame_idx}")

        CAR_positions = []
        WATER_positions = []

        # Parse each frame
        for line in frame_lines[2:-1]:  # Skip header and box line
            residue_name = line[5:10].strip()
            atom_name = line[10:15].strip()
            x = float(line[20:28].strip())
            y = float(line[28:36].strip())
            z = float(line[36:44].strip())

            if residue_name == "UNL":
                CAR_positions.append((atom_name, np.array([x, y, z])))
            elif residue_name == "SOL":
                WATER_positions.append((atom_name, np.array([x, y, z])))

        # Calculate COM for CAR (CO2)
        CAR_COM = []
        for i in range(0, len(CAR_positions), 3):
            co_pos = CAR_positions[i][1]
            o1_pos = CAR_positions[i + 1][1]
            o2_pos = CAR_positions[i + 2][1]
            com = (co_pos * C_M + (o1_pos + o2_pos) * O_M) / (C_M + 2 * O_M)
            CAR_COM.append(com)

        # Calculate COM for WATER (SOL)
        WATER_COM = []
        for i in range(0, len(WATER_positions), 4):
            o_pos = WATER_positions[i][1]
            h1_pos = WATER_positions[i + 1][1]
            h2_pos = WATER_positions[i + 2][1]
            com = (o_pos * O_M + (h1_pos + h2_pos) * H_M) / (O_M + 2 * H_M)
            WATER_COM.append(com)

        # Get box length from the last line of the frame
        box_line = frame_lines[-1].strip()
        box_length = float(box_line.split()[0])  # Assuming a cubic box

        # Store the data for the current frame
        all_car_com.append(np.array(CAR_COM))
        all_water_com.append(np.array(WATER_COM))
        all_box_lengths.append(box_length)

    return np.array(all_steps), all_car_com, all_water_com, all_box_lengths

def COM_calculation_mdtraj(traj):
    """Calculate the center of mass (COM) for CO2 and water molecules from an mdtraj trajectory."""
    C_M = params.C_M  # Carbon atomic mass
    O_M = params.O_M  # Oxygen atomic mass
    H_M = params.H_M  # Hydrogen atomic mass
    
    # Data structures to store results
    all_car_com = []
    all_water_com = []
    all_box_lengths = []
    # all_steps = np.arange(traj.n_frames)  # Use frame index as step
    timesps = traj.time
    
    # Identify atom indices
    topology = traj.topology
    CAR_indices = []  # CO2 indices
    WATER_indices = []  # Water indices
    
    for residue in topology.residues:
        if residue.name == "UNL":  # CO2 molecule
            atoms = [atom.index for atom in residue.atoms]
            if len(atoms) == 3:  # Ensure CO2 has exactly 3 atoms
                CAR_indices.append(atoms)
        elif residue.name == "HOH":  # TIP4P/Ice water molecule
            oxygen = None
            hydrogens = []
            for atom in residue.atoms:
                if atom.element.symbol == "O":
                    oxygen = atom.index
                elif atom.element.symbol == "H":
                    hydrogens.append(atom.index)
            if oxygen is not None and len(hydrogens) == 2:  # Ignore MW/EP site
                WATER_indices.append([oxygen] + hydrogens)
    
    for frame in range(traj.n_frames):
        CAR_COM = []
        WATER_COM = []
        
        xyz = traj.xyz[frame]  # Positions of all atoms in this frame
        
        # Compute CO2 COM
        for indices in CAR_indices:
            co_pos = xyz[indices[0]]  # Carbon
            o1_pos = xyz[indices[1]]  # Oxygen 1
            o2_pos = xyz[indices[2]]  # Oxygen 2
            com = (co_pos * C_M + (o1_pos + o2_pos) * O_M) / (C_M + 2 * O_M)
            CAR_COM.append(com)
        
        # Compute Water COM
        for indices in WATER_indices:
            o_pos = xyz[indices[0]]  # Oxygen
            h1_pos = xyz[indices[1]]  # Hydrogen 1
            h2_pos = xyz[indices[2]]  # Hydrogen 2
            com = (o_pos * O_M + (h1_pos + h2_pos) * H_M) / (O_M + 2 * H_M)
            WATER_COM.append(com)
        
        # Extract box length (assuming cubic box)
        box_length = traj.unitcell_lengths[frame, 0]
        
        all_car_com.append(np.array(CAR_COM))
        all_water_com.append(np.array(WATER_COM))
        all_box_lengths.append(box_length)
    
    return timesps, all_car_com, all_water_com, all_box_lengths



def is_within_cone(guest1, guest2, water, box_length, angle = np.pi / 2):
    """
    Computes, where a given water lies in a cone spanned by two guest molecules guest 1 and guest 2.

    Parameters:
        guest 1 (int): Index of guest 1.
        guest 2 (int): Index of guest 2.
        water (int): Index of water.
        box_length (float) : Length of the box, used for the pbc.
        box_length (float): Length of the periodic box.
        angle (float) : Angle that cone spans from connecting line (g1 - g2) outwards.

    Returns:
        bool: Stating whether water is in (True) or outside (False) the cone.
    """

    # calculate distance vectors
    vec_g1_g2 = guest2 - guest1                                         # points from guest 1 to guest 2
    vec_g1_w = water - guest1                                           # points from guest 1 to water 
    vec_g2_w = guest2 - water                                           # points from water to guest 2 (so that angle calculation works as exptected)

    # periodic boundary conditions
    vec_g1_g2 -= np.rint(vec_g1_g2 / box_length) * box_length           # if distance is longer than half box length, distance gets decreased by box size (vector is preserved)
    vec_g1_w -= np.rint(vec_g1_w / box_length) * box_length
    vec_g2_w -= np.rint(vec_g2_w / box_length) * box_length
    
    # normalize vectors
    vec_g1_g2 /= np.linalg.norm(vec_g1_g2)                              # only unit-vectors are needed for angles
    vec_g1_w /= np.linalg.norm(vec_g1_w)
    vec_g2_w /= np.linalg.norm(vec_g2_w)
    
    # compute angle between vectors
    theta_g1 = np.arccos( np.dot(vec_g1_g2, vec_g1_w) )                 # |a| and |b| are 1 
    theta_g2 = np.arccos( np.dot(vec_g1_g2, vec_g2_w) )
    
    # return whether both angles are within cone
    return theta_g1 <= ( angle / 2 ) and theta_g2 <= ( angle / 2 )

def distance_matrix_pbc(pos1, pos2, box_length):
    """
    Computes the pairwise distance matrix between two position arrays
    using periodic boundary conditions.

    Parameters:
        pos1 (ndarray): A 2D array of shape (N, 3) representing N points.
        pos2 (ndarray): A 2D array of shape (M, 3) representing M points.
        box_length (float): Length of the periodic box.

    Returns:
        ndarray: A matrix of distances of shape (N, M).
    """

    # ensure 2d arrays
    pos1 = np.atleast_2d(pos1)                                          # (3,) input is converted into (1,3) when distances to only one guest are calculated
    pos2 = np.atleast_2d(pos2)
    
    #calculate distances
    disp = pos1[:, np.newaxis] - pos2                                   # newaxis ensures, that we can subtract the whole pos2 array from every single pos1 entry, every distance will be calculated
    
    #periodic boundary conditions
    disp -= box_length * np.rint(disp / box_length)                    
    
    # turn distance vectors into absolute distances
    distance_matrix = np.linalg.norm(disp, axis=2)

    if pos1.size == 0 or pos2.size == 0:
        print('error')
        return np.inf

    else:
        return distance_matrix

def distance_pbc(pos1, pos2, box_length):
    """
    Computes absolute distances between two points using periodic boundary conditions.

    Parameters:
        pos1 (ndarray): A 1D array of shape (, 3) representing 1 point.
        pos2 (ndarray): A 1D array of shape (, 3) representing 1 point.
        box_length (float): Length of the periodic box.

    Returns:
        disp (float): Scalar distance between pos1 and pos2.
    """

    disp = pos2 - pos1

    disp -= box_length * np.rint(disp / box_length)                    

    disp = linalg.norm(disp)

    return disp

def largest_cluster_size(adjacency_matrix):
    # If the adjacency matrix is all False (no edges), return 0
    if np.sum(adjacency_matrix) == 0:
        return 0,0
    
    # convert the adjacency matrix to a graph
    G = nx.from_numpy_array(adjacency_matrix)
    
    # find all connected components (each component is a cluster)
    connected_components = list(nx.connected_components(G))
    
    # calculate the size of each cluster
    cluster_sizes = np.array([len(component) for component in connected_components])
    # get the size of the largest cluster

    largest_cluster_idx = np.argmax(cluster_sizes)
    largest_cluster = list(connected_components[largest_cluster_idx])

    largest_cluster_count = len(connected_components[largest_cluster_idx])
    
    return largest_cluster_count, largest_cluster 

#@jit(nopython=True)

def MCG(CAR_COM, WATER_COM, box_length, order = 1, guest_cutoff = 0.9, water_cutoff = 0.6):

    """
    Computes the largest cluster of monomers Mutually Coordinated Guest (MCG) order parameter in a system based on pairwise distances,
    adjacency conditions, and geometric constraints involving water molecules. 

    Parameters:
        CAR_COM (ndarray): A 2D array of shape (N, 3) containing the center of mass (COM) positions
                           of N guest molecules (e.g., CO2) in 3D space.
        WATER_COM (ndarray): A 2D array of shape (M, 3) containing the center of mass (COM) positions
                             of M water molecules in 3D space.
        box_length (float): The length of the cubic simulation box for applying periodic boundary
                            conditions (PBC).
        order (int, optional): The minimum number of adjacencies required for a guest molecule to
                               be part of the MCG. Default is 1.
        guest_cutoff (float, optional): Maximum distance between two guest molecules
                                                  to consider them as neighbors. Default is 0.9.
        water_cutoff (float, optional): Maximum distance between a guest molecule and a water
                                              molecule for a valid bond. Default is 0.6.

    Returns:
        tuple:
            - MCG_OP (int): The size of the largest cluster of monomers satisfying the conditions.
            - MCG_monomers_position (ndarray): A 2D array of shape (K, 3) containing the COM
                                               positions of the monomers in the largest cluster.

    Steps:
        1. Compute the pairwise distance matrix between all guest molecules (CAR_COM) using
           `distance_matrix_pbc`, applying periodic boundary conditions (PBC).

        2. Construct an adjacency matrix (`adjacency_matrix_distance`) based on the condition
           that distances between guest molecules are less than or equal to `guest_cutoff`.

        3. Extract unique guest molecule pairs satisfying the adjacency condition using the
           upper triangular part of the adjacency matrix.

        4. For each pair of adjacent guest molecules, identify nearby water molecules that satisfy
           a relaxed distance condition (`guest_cutoff * 1.1`).

        5. For each nearby water molecule, check if it lies within the geometric cone defined by
           the two guest molecules and verify that its distance to both guest molecules is less than
           `water_cutoff`. If these conditions are satisfied, count the water molecule as a
           valid bond.

        6. If at least 5 valid bonds are formed for the pair of guest molecules, update the
           adjacency matrix to indicate a monomer connection between them.

        7. Remove guest molecules that do not satisfy the required number of adjacencies (1 for MCG-1).

        8. Perform cluster analysis on the final adjacency matrix to determine the largest cluster
           of connected guest molecules. This uses the `largest_cluster_size` helper function.

    """
    
    # build symmetric distance matrix between guests
    distance_matrix = distance_matrix_pbc(CAR_COM, CAR_COM, box_length)
    
    # bool matrix stating whether distance is smaller than or equal to minimum neighbour distance (default = 0.9)
    adjacency_matrix_distance = (distance_matrix <= guest_cutoff)

    # only take upper triangular because of symmetry and remove self-ajacency by k = 1
    adjacency_matrix_distance = np.triu(adjacency_matrix_distance, k = 1) 

    # set up final adjacency matrix, for bonds satifying MCG requirements
    adjacency_matrix = np.zeros_like(adjacency_matrix_distance)
    
    # build inhomogenous list of arrays, first index selects guest and array at that position gives indices of adjacent guests
    # Ex. adjacent_guests_idx[1] = (array([0, 2]),) means guest 1 is adjacent to guest 0 and 2
    adjacent_guests_idx = [np.where(arr == True) for arr in adjacency_matrix_distance[:]]    

    successful_pairs = []

    for guest1_idx, guest2_idx_vector in enumerate(adjacent_guests_idx):

        guest2_idx_vector = np.array(guest2_idx_vector)[0]                  # list is inhomogeneous so we extract a numpy array

        if guest2_idx_vector.size != 0:

            distance_matrix_g1_waters = distance_matrix_pbc(CAR_COM[guest1_idx], WATER_COM, box_length)[0]      # this is by function definition a matrix, but we only need a vector here
                                                                                                                # dimensions are len(waters) since its always the distance to g1
            distance_matrix_g2_waters = distance_matrix_pbc(CAR_COM[guest2_idx_vector], WATER_COM, box_length)  # dimensions are (len(guest2_idx_vector, len(waters)) so first index is guest index, second water index

            if np.sum(distance_matrix_g1_waters) < 5:                       #if g1 already doesn't have enough water neighbours it's invalid
                break
                
            g1_water_mask = (distance_matrix_g1_waters < water_cutoff) 
            g2_water_mask = (distance_matrix_g2_waters < water_cutoff) 

            mutual_waters_g2_mask = g1_water_mask & g2_water_mask           # waters that g2 has in common with g1 (index 1 indicates g2, index two indicated water index)

                                                                            # // TODO here i can directly take out g2s that don't have enough water

            for count, guest2_idx in enumerate(guest2_idx_vector):          # count is needed to not lose access to the position of the current index in the mask

                N_w = 0                                                     # water count
                mutual_waters_idx = np.where(mutual_waters_g2_mask[count] == True)[0]

                successful_waters = []                                      # this is only for the animation

                for water_idx in mutual_waters_idx:
                    
                    if is_within_cone(CAR_COM[guest1_idx], CAR_COM[guest2_idx], WATER_COM[water_idx], box_length) == True:

                        successful_waters.append(water_idx)
                        N_w += 1

                    if N_w == 5:

                        adjacency_matrix[guest1_idx, guest2_idx] = True     # connected components works for undirected, single edge connections

                        set = [guest1_idx, guest2_idx, *successful_waters]  # for animation
                        successful_pairs.append(set)

                        break


            
    # adjacent_indices = np.where(adjacency_matrix == True)
    # unique_monomers = unique_numbers = list(set(np.concatenate(adjacent_indices)))
    # MCG_monomers_position = CAR_COM[unique_monomers, :]
    
    # take out non MCG monomers (not more than given number of adjacencies, in this case one so its no issue really)

    # for guest_idx in range(adjacency_matrix.shape[1]): # should be symmetric
    #     #print(f'Guest:{guest_idx}')
    #     #print(f'SUM: {np.sum(adjacency_matrix[:, guest_idx]) + np.sum(adjacency_matrix[guest_idx, :])}')
    #     if np.sum(adjacency_matrix[:, guest_idx]) + np.sum(adjacency_matrix[guest_idx, :]) < order:
    #         adjacency_matrix[:, guest_idx], adjacency_matrix[guest_idx, :] = 0, 0
            
    #print(adjacency_matrix)
    # do cluster analysis on the adjacency matrix

    MCG_OP, monomer_idx = largest_cluster_size(adjacency_matrix)
    MCG_monomers_position = CAR_COM[monomer_idx]

    return MCG_OP, successful_pairs, MCG_monomers_position