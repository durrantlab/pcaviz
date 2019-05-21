import json
import csv

import MDAnalysis
from MDAnalysis.analysis import rms
from MDAnalysis.analysis.base import AnalysisFromFunction
from MDAnalysis.coordinates.memory import MemoryReader
import numpy as np
from sklearn.decomposition import PCA

import _Utils
from _UserParams import get_params

# Calculate PCA of the first trajectory using MDAnalysis' class.
def compress_PCA(traj, selection, cum_var):
    """
    This function performs PCA analysis using the MDAnalysis PCA analysis
    class.

    It obtains the information necessary for expansion back into a trajectory.

    :param MDAnalysis.Universe traj: The universe containing the trajectory to
                                     be compressed.

    :param str selection: The atom selection to be included in the compression
                          e.g. protein, name CA...

    :param float cum_var: The amount of variance to be explained cumulatively
                          by calculated components.

    :returns: The principal component vectors.
    :rtype: :class:'ndarray'

    :returns: The coordinates of selection in traj projected onto the
              calculated components.
    :rtype: :class:'ndarray'

    :returns: The coordinates of the average atomic positions of selection in
              traj in XYZ space.
    :rtype: :class:'ndarray'

    :returns: The first frame coordinates of selection in traj in XYZ space.
    :rtype: :class:'ndarray

    :returns: The atom type and number of the selection in traj.
    :rtype :class:'list'

    :returns: The residue number and name of the selection in traj.
    :rtype :class:'list'
    """

    # Get the indices of the atoms that match the selection.
    atom_group = traj.select_atoms(selection)
    idx_of_selected_atoms = atom_group.indices

    # Get the flattened coordinates of each frame.
    coor_data = []
    for frame in traj.trajectory:
        coor_data.append(frame.positions[idx_of_selected_atoms].ravel())
    coor_data = np.array(coor_data)

    # Append coordinates to python list then use ravel to flatten coordinates without
    # having to allocate new memory
    # Turning coord list into numpy array because it reduces memory consumption
    # You get a flatten numpy array containing the coordinates
    pca = PCA(n_components=min(coor_data.shape))
    pca.fit(coor_data)
    pca_vectors = pca.components_

    # Number of components should be equal to number of elements in coor_data array
    # Number of principal components = 3xs the number of atoms
    # pca.fit fits model to coordinate data
    # PCA vectors contain principal axes in PCA space, show max variance in data

    # Also get the coefficients for each frame.
    coords_project_onto_pca_space = pca.transform(coor_data)

    # coords_project_onto_pca_space applies PCA dimensionality reduction to coor_data

    # Calculate how many components are required to cumulatively explain the
    # desired variance.
    var = pca.explained_variance_ratio_
    cumulated_variance = np.array([np.sum(var[:i+1]) for i in range(var.size)])
    idxs_of_above_cum_var = np.where(cumulated_variance >= cum_var)[0]

    # var = % of variance explained by each of the selected components cum_var
    # sums up the % of variance explained by the components

    # If there are no components above threshold (because cum_var == 1.0),
    # keep all components. If not, keep only the number required to account
    # for the user-specified variance.
    if idxs_of_above_cum_var.size == 0:
        n_pcs = cumulated_variance.size
    else:
        n_pcs = idxs_of_above_cum_var[0] + 1

    # Keep only the appropriate number of components and coefficients.
    pca_vectors = pca_vectors[:n_pcs]
    coords_project_onto_pca_space = coords_project_onto_pca_space[:,:n_pcs]

    # Also calculate the mean structure.
    coords_avg_atoms = pca.mean_

    # Save the coordinates of the first structure too (so we can get the
    # correct bonds in the browser).
    coords_first_frame = coor_data[0]

    # Only select atoms of interest.
    atomgroup = traj.select_atoms(selection)
    # data_json will contain residue/atom info for user selection
    data_json = []

    # Save atom information.
    last_key = ""
    for a in atomgroup:
        at_name = a.name
        res_name = a.resname
        res_id = a.resid

        # If statement will check to see if residue id is present if it's not
        # present it will add the residue id, residue name, and the residue's
        # atoms.
        key = str(res_id) + res_name # Because list is not hashable.
        if key != last_key:
            data_json.append([int(res_id), res_name, at_name])
            last_key = key
        else:
            data_json.append(at_name)

    # Return only the information necessary for compression and expansion.
    return pca_vectors, coords_project_onto_pca_space, coords_avg_atoms, \
           coords_first_frame, data_json

def expand_PCA(pca_vectors, coords_project_onto_pca_space, coords_avg_atoms, precision):
    """
    This function decompresses a trajectory compressed with the compress_PCA function.

    It uses only the information necessary for expansion back into XYZ-space and 
    returns a list of the decompressed atomic coordinates.

    :param list pca_vectors: The principal components on which the trajectory
                             has been projected.

    :param list coords_project_onto_pca_space: The projected trajectory coordinates
                                               in PCA-space.

    :param list coords_avg_atoms: The average atomic coordinated of the
                                  original trajectory.

    :param int precision: How many decimal places to which the PCA vectors,
                          projection, and atomic coordinates were rounded.

    :returns: The decompressed trajectory. A list containing ndarrays of the atomic
              coordinates at each frame in the same format as an ordinary MDAnalysis
              Universe atom_selection.positions field.
    :rtype: :class:'list'
    """
    # print(coords_project_onto_pca_space[0])
    pca_vectors = np.true_divide(pca_vectors, (10**precision))
    coords_project_onto_pca_space = np.true_divide(coords_project_onto_pca_space,
                                                   (10**precision))
    coords_avg_atoms = np.true_divide(coords_avg_atoms, (10**precision))


    # Prepare list to hold expanded trajectory coordinates.
    expanded_coordinates = []
    # print(coords_project_onto_pca_space[0])
    # print(np.dot(coords_project_onto_pca_space[0], pca_vectors))
    # print(np.dot(coords_project_onto_pca_space[0], pca_vectors) + coords_avg_atoms)

    # Construct each frame's coordinates one at a time.
    for frame_num, frame_coefficients in enumerate(coords_project_onto_pca_space):

        # Calculate flattened XYZ coordinate delta values.
        XYZ_coor_delta = np.dot(frame_coefficients, pca_vectors)
        
        # Add the average atomic positions.
        XYZ_coor = XYZ_coor_delta + coords_avg_atoms

        expanded_coordinates.append(np.reshape(XYZ_coor, (len(XYZ_coor) // 3, 3)))
    
    # Return the expanded coordinates.
    return expanded_coordinates



# Main function.
if __name__ == '__main__':
    # Get user provided parameters and unpack frequently used parameters
    params = get_params()
    coor_file = params['coor_file']
    top_file = params['top_file']
    num_decimals = params['precision']

    # Load a trajectory for input file, aligning it if it is the first.
    traj = _Utils.load_traj(
        top_file, coor_file, align=False  # align_sel=params['align_sel']
    )
    # If a dual function file is being used, set coor_file to it for filename purposes.
    if coor_file is None:
        coor_file = top_file

    # Stride or trim the trajectory as desired, this transfers the trajectory
    # to memory.
    if params['stride'] > 1 or params['starting_frame'] > 0:
        sel = traj.select_atoms('all')
        analyze = AnalysisFromFunction(lambda ag: ag.positions.copy(), traj.atoms)
        coor = analyze.run().results
        traj = MDAnalysis.Merge(sel)
        traj = traj.load_new(coor[params['starting_frame']::params['stride']],
                             format=MemoryReader)

    # Calculate trajectory on principal components and retrieve PCA
    # information
    vect_components, PCA_coeff, coords_avg_atoms, coords_first_frame, data_json_info = \
        compress_PCA(traj, params['selection'], params['cum_var'])

    # Keep only the same number of top vect_components as the number of PCA_coeff
    vect_components = vect_components[:len(PCA_coeff[0])]

    # Compress information further by rounding and converting to ints.
    vect_components = _Utils.compress_list(vect_components, num_decimals)
    PCA_coeff = _Utils.compress_list(PCA_coeff, num_decimals)
    coords_avg_atoms = _Utils.compress_list(coords_avg_atoms,
                                            num_decimals,
                                            flatten=True)
    coords_first_frame = _Utils.compress_list(coords_first_frame,
                                              num_decimals,
                                              flatten=True)

    # If the user request to check the accuracy of the expanded trajectory.
    if params['check_accuracy']:
        # Prepare list to hold RMSD values.
        RMSD_lst = []

        # Calculate expanded coordinates.
        expanded_trajectory = expand_PCA(vect_components, PCA_coeff,
                                         coords_avg_atoms, num_decimals)


        # Calculate RMSD value of each frame between original and expanded trajectory.
        for frame_id, traj_fram in enumerate(traj.trajectory):
            original_positions = traj.select_atoms(params['selection']).positions
            RMSD_lst.append(rms.rmsd(original_positions, expanded_trajectory[frame_id]))
        
        # Create a list for outputting frame number.
        frame_numbers = list(range(len(RMSD_lst)))
        frame_numbers.append('mean')
        frame_numbers.append('standard deviation')

        # Also output the mean and standard deviation.
        mean_RMSD = np.mean(RMSD_lst)
        std_RMSD = np.std(RMSD_lst)
        RMSD_lst.extend([mean_RMSD, std_RMSD])
        print("The mean RMSD between the original and expanded trajectory is: " + \
              str(mean_RMSD))
        print("The standard deviation of the RMSD between the original and " + \
              "expanded trajectory is: " + str(std_RMSD))
        
        # Save a table containing results.
        tbl_flnm = _Utils.output_filename('starting_vs_expanded_RMSD',
							       'csv',
							       coor_file=coor_file,
							       output_dir=params['output_dir'])
        with open(tbl_flnm, 'w') as tf:
            table_writer = csv.writer(tf)
            table_writer.writerows(zip(frame_numbers, RMSD_lst))

    # dictionary to be outputted as a json.
    my_dict = {}
    my_dict['vecs'] = vect_components
    my_dict['coeffs'] = PCA_coeff
    my_dict['coors'] = coords_avg_atoms
    my_dict['first_coors'] = coords_first_frame
    my_dict['params'] = params
    my_dict['res_info'] = data_json_info

    # Write final json file.
    json_txt = json.dumps(my_dict) # , sort_keys=True, indent=4)

    # Remove spaces to compress further.
    json_txt = json_txt.replace(' ', '')

    # Obtain output file name.
    output_name = _Utils.output_filename('compressed',
                                         'json',
                                         coor_file=coor_file,
                                         output_dir=params['output_dir'])

    # Write the file.
    with open(output_name, 'w') as write_file:
        write_file.write(json_txt)
