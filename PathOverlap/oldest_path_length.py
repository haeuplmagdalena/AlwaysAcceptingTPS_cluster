def find_age_of_oldest_path(indices, directions, lengths, force_eval):
    '''
        Parameters:
            indices : list of tuples
                (old_index, new_index) are shooting point indices on the old and new path
            directions : list of int
                (+1 = forward shooting, -1 = backward shooting)
            first_path_length : int
                Length of the initial transition path (number of frames)

        Returns:
            AOPS_list : list of int
                Ages (in number of accepted paths) of the oldest path segment for each starting trajectory
            LOPS_list : list of lists
                Evolution of the oldest path segment lengths over successive accepted paths.
                Each sublist corresponds to one starting position in AOPS_list.
    '''

    LOPS_list = []  #lengths of oldest path segments
    force_eval_list = [] #force evaluations at specific path lengths 
    AOPS_list = []  #ages of oldest path segments
    
    # loop over all possible paths as starting paths
    for start in range(1, len(indices)):    #oldest path is only defined once another path has been accepted (so starts at idx 1)
        
        first_path_length = lengths[start - 1]
        
        # Setting up the first OPS after one new path has been generated
        _, first_path_n = indices[start]    #shooting point index on new path
        first_path_dir = directions[start]

        if first_path_dir == 1:
            SOPS = 0                    #SOPS (Start of oldest path segment) and EOPS (end of oldest path segment) are always defined in terms of the indices on the newest path
            EOPS = first_path_n
    
        if first_path_dir == -1:
            SOPS = first_path_n
            EOPS = first_path_length
    
        LOPS = EOPS - SOPS
        current_lengths = [LOPS]
        current_force_eval = [force_eval[start]]
        age = 1

        total_force_eval = [force_eval[start]]
    
    
        for idx, (o, n) in enumerate(indices[start + 1:], start = start + 1):       #(o, n) are sp indices on old/new path

    
            if directions[idx] == 1:        #easier version, no shift of indices necessary
                if n <= EOPS:
                    EOPS = n
                #otherwise the end stays where it was before
        
            elif directions[idx] == -1:     #indices need to be shifted
                shift = n - o
                if o >= SOPS:
                    SOPS = n
                    EOPS += shift
                else:                        # length stays the same but indices need to be shifted to fit current path 
                    SOPS += shift
                    EOPS += shift

    
            LOPS = EOPS - SOPS

            #when path dies
            if LOPS <= 0:
                AOPS_list.append(age)
                LOPS_list.append(current_lengths)
                force_eval_list.append(current_force_eval)
                break
                
            age += 1
            current_lengths.append(LOPS)
            #print(force_eval[idx])
            total_force_eval += force_eval[idx]
            current_force_eval.extend(total_force_eval)
            
        else:  # OPS survived all paths
            break  # Stop early since following starts will also survive

    return AOPS_list, LOPS_list, force_eval_list