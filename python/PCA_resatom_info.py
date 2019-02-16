import json

import MDAnalysis
from MDAnalysis.analysis.base import AnalysisFromFunction
from MDAnalysis.analysis import pca
from MDAnalysis.coordinates.memory import MemoryReader
import numpy as np

import _Utils
from _UserParams import get_params


# Calculate PCA of the first trajectory using MDAnalysis' class.
def get_PCA_trajectory(traj, selection, cum_var):
    """
    This function performs PCA analysis using the MDAnalysis PCA analysis class.

    It obtains the information necessary for expansion back into a trajectory.

    :param MDAnalysis.Universe traj: The universe containing the trajectory to be compressed.

    :param str selection: The atom selection to be included in the compression e.g. protein, name CA...

    :param float cum_var: The amount of variance to be explained cumulatively by calculated components.

    :returns: The principal component vectors.
    :rtype: :class:'ndarray'

    :returns: The coordinates of selection in traj projected onto the calculated components.
    :rtype: :class:'ndarray'

    :returns: The coordinates of the average atomic positions of selection in traj in XYZ space.
    :rtype: :class:'ndarray'

    :returns: The atom type and number of the selection in traj.
    :rtype :class:'list'

    :returns: The residue number and name of the selection in traj.
    :rtype :class:'list'
    """

    # Calculate principal components and variation of specific atom selection.
    get_trajectory_pca = pca.PCA(traj, selection)
    pca_space = get_trajectory_pca.run()

    # Access the components themselves.
    pca_vectors = pca_space.p_components

    # Calculate how many components are required to cumulatively explain
    # desired variance.
    n_pcs = np.where(get_trajectory_pca.cumulated_variance > cum_var)[0][0]

    # Only select atoms of interest.
    atomgroup = traj.select_atoms(selection)
    #data_json will contain residue/atom info for user selection
    data_json = []
    # loops over atom information turns the atom name, residue name, and resid
    # into a string

    for a in atomgroup:
        at_name = str(a.name)
        res_name = str(a.resname)
        res_id = str(a.resid)
        # if statement will check to see if residue id is present if it's not
        # present it will add the residue id, residue name, and the residue's
        # atoms
        if [res_id, res_name] not in data_json:
            data_json.append([res_id, res_name])
            data_json.append(at_name)
        elif [res_id, res_name] in data_json:
            data_json.append(at_name)


    # Get average atomic coordinates.
    coords_avg_atoms = traj.trajectory.timeseries(asel=atomgroup).mean(axis=1)

    # And project onto the principal components.
    coords_project_onto_pca_space = get_trajectory_pca.transform(atomgroup,n_components=n_pcs)


    # Return only the information necessary for compression and expansion.
    return pca_vectors, coords_project_onto_pca_space, coords_avg_atoms, data_json

# Main function.
if __name__ == '__main__':
    # Get user provided parameters and unpack frequently used parameters
    params = get_params()
    coor_file = params['coor_file']
    top_file = params['top_file']
    num_decimals = params['precision']

    # Load a trajectory for coor_file, aligning it if it is the first.
    if coor_file is not None:
        traj = _Utils.load_traj(top_file, coor_file, align_sel=params['align_sel'])
    else:
        traj = _Utils.load_traj(top_file, align_sel=params['align_sel'])

    # Stride or trim the trajectory as desired, this transfers the trajectory to memory.
    if params['stride'] > 1 or params['starting_frame'] > 0:
        sel = traj.select_atoms('all')
        analyze = AnalysisFromFunction(lambda ag: ag.positions.copy(), traj.atoms)
        coor = analyze.run().results
        traj = MDAnalysis.Merge(sel)
        traj = traj.load_new(coor[params['starting_frame']::params['stride']],
                                         format=MemoryReader)

    # Calculate trajectory on principal components and retrieve PCA information
    vect_components, PCA_coeff, coords_avg_atoms, data_json_info = get_PCA_trajectory(
        traj, params['selection'], params['cum_var']
    )

    #### JDD ADDITION:
    # Keep only the save number of top vect_components as the number of PCA_coeff
    vect_components = vect_components[:len(PCA_coeff[0])]

    # Compress information further by rounding and converting to ints.
    # Save to dictionary to be outputted as a json.

    my_dict = {}
    my_dict['vecs'] = _Utils.compress_list(vect_components, num_decimals)
    my_dict['coeffs'] = _Utils.compress_list(PCA_coeff, num_decimals)
    my_dict['coors'] = _Utils.compress_list(coords_avg_atoms, num_decimals, is_coor=True)
    my_dict['params'] = params
    my_dict['res_info'] = data_json_info

    # Write final json file.
    json_txt = json.dumps(my_dict) # , sort_keys=True, indent=4)

    # Remove spaces to compress further.
    json_txt = json_txt.replace(' ', '')

    # Obtain output file name.
    output_name = _Utils.output_filename('compressed', 'json', coor_file=coor_file,
                                         output_dir=params['output_dir'])

    # Write the file.
    with open(output_name, 'w') as write_file:
        write_file.write(json_txt)