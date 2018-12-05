import json

import MDAnalysis
from MDAnalysis.analysis import pca
import numpy as np

import _Utils
from _UserParams import get_params


# Calculate PCA of the first trajectory using MDAnalysis' class.
def get_PCA_trajectory(traj, selection, cum_var):

    # Calculate principal components and variation of specific atom selection.
    get_trajectory_pca = pca.PCA(traj, selection)
    pca_space = get_trajectory_pca.run()

    # Access the components themselves.
    pca_vectors = pca_space.p_components

    # Calculate how many components are required to cumulatively explain desired variance.
    n_pcs = np.where(get_trajectory_pca.cumulated_variance > cum_var)[0][0]

    # Only select atoms of interest.
    atomgroup = traj.select_atoms(selection)

    # Get average atomic coordinates.
    coords_atoms = atomgroup.positions

    # And project onto the principal components.
    coords_project_onto_pca_space = get_trajectory_pca.transform(atomgroup,n_components=n_pcs)
    
    # Return only the information necessary for compression and expansion
    return pca_vectors, coords_project_onto_pca_space, coords_atoms

# Main function.
if __name__ == '__main__':
    # Get user provided parameters and unpack frequently used parameters
    params = get_params()

    coor_file = params['coor_file']
    top_file = params['top_file']
    num_decimals = params['num_decimals']

    # Load a trajectory for coor_file, aligning it if it is the first.
    if coor_file is not None:
        traj = _Utils.load_traj(top_file,
                                coor_file,
                                align_sel=params['align_sel'],
                                ps_timestep=params['ps_timestep'])
    else:
        traj = _Utils.load_traj(top_file,
                                align_sel=params['align_sel'],
                                ps_timestep=params['ps_timestep'])

    # Stride or trim the trajectory as desired.
    if params['num_stride'] > 1 or params['starting_frame'] > 0:
        traj = _Utils.stride_and_trim_traj(traj, params['num_stride'],
                                           params['starting_frame'],
                                           params['ps_timestep'])

    # Calculate trajectory on principal components and retrieve PCA information
    vect_components, PCA_coeff, coords_atoms = get_PCA_trajectory(traj,
                                                                params['selection'],
                                                                params['cum_var'])

    # Compress information further by rounding and converting to ints.
    # Save to dictionary to be outputted as a json.
    my_dict = {}
    my_dict['vecs'] = _Utils.compress_list(vect_components, num_decimals)
    my_dict['coeffs'] = _Utils.compress_list(PCA_coeff, num_decimals)
    my_dict['coors'] = _Utils.compress_list(coords_atoms, num_decimals, is_coor=True)
    my_dict['params'] = params

    # Write final json file.
    json_txt = json.dumps(my_dict) # , sort_keys=True, indent=4)

    # Remove spaces to compress further.
    json_txt = json_txt.replace(' ', '')

    with open('PCA_file.json', 'w') as write_file:
        write_file.write(json_txt)
