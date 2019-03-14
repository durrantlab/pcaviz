import json
import MDAnalysis
from MDAnalysis.analysis.base import AnalysisFromFunction
from MDAnalysis.coordinates.memory import MemoryReader
import numpy as np
from sklearn.decomposition import PCA

import _Utils
from _UserParams import get_params

# Calculate PCA of the first trajectory using MDAnalysis' class.
def get_PCA_trajectory(traj, selection, cum_var):
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

    # Calculate all the principal components.
    pca = PCA(n_components=min(coor_data.shape))
    pca.fit(coor_data)
    pca_vectors = pca.components_

    # Also get the coefficients for each frame.
    coords_project_onto_pca_space = pca.transform(coor_data)

    # Calculate how many components are required to cumulatively explain the
    # desired variance.
    var = pca.explained_variance_ratio_
    cumulated_variance = np.array([np.sum(var[:i+1]) for i in range(var.size)])
    idxs_of_above_cum_var = np.where(cumulated_variance > cum_var)[0]

    # If there are no components above threshold (because cum_var == 1.0),
    # keep all components. If not, keep only the number required to account
    # for the user-specified variance.
    n_pcs = cumulated_variance.size if idxs_of_above_cum_var.size == 0 else idxs_of_above_cum_var[0]

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
    return pca_vectors, coords_project_onto_pca_space, coords_avg_atoms, coords_first_frame, data_json

def un_pca(vect_components, PCA_coeff, coords_first_frame, compressed_json_filename):
    """A temporary function that should "undo" to PCA (back to Cartesian space).

    :param vect_components: The PCA vectors.
    :type vect_components: numpy.array
    :param PCA_coeff: The PCA coefficients for each frame.
    :type PCA_coeff: numpy.array
    :param coords_first_frame: The first-frame coordinates.
    :type coords_first_frame: numpy.array
    :param compressed_json_filename: The name of the compressed json file.
    :type compressed_json_filename: string.
    """

    # Dummy variables.
    atomName = "X"
    resName = "XXX"
    chain = "X"
    resID = "0"
    element = "X"

    with open(compressed_json_filename + ".uncompressed.pdb", "w") as f:

        # Go through each of the frames.
        for frame_idx, frame_coeffs in enumerate(PCA_coeff):
            f.write("MODEL " + str(frame_idx) + "\n")

            # Get the coordinates for the current frame.
            frame_coors_delta = np.dot(frame_coeffs, vect_components)
            frame_coors = coords_first_frame + frame_coors_delta
            frame_coors_3d = np.reshape(frame_coors, (len(frame_coors) / 3,3))

            # Go through each of the coordinates in this frame.
            for x, y, z in frame_coors_3d:
                # Make a PDB line
                pdb_line = "ATOM  " + str(frame_idx).rjust(5) + \
                           atomName.rjust(5) + resName.rjust(4) + \
                           " " + chain + resID.rjust(4) + "    " + \
                           "{:8.3f}{:8.3f}{:8.3f}" + \
                           "  1.00  0.00          " + \
                           element + "  \n"

                pdb_line = pdb_line.format(x, y, z)

                f.write(pdb_line)

            f.write("ENDMDL\n")

# Main function.
if __name__ == '__main__':
    # Get user provided parameters and unpack frequently used parameters
    params = get_params()
    coor_file = params['coor_file']
    top_file = params['top_file']
    num_decimals = params['precision']

    # Load a trajectory for coor_file, aligning it if it is the first.
    if coor_file is not None:
        traj = _Utils.load_traj(
            top_file, coor_file, align=False  # align_sel=params['align_sel']
        )
    else:
        traj = _Utils.load_traj(
            top_file, align=False  # align_sel=params['align_sel']
        )

    # Stride or trim the trajectory as desired, this transfers the trajectory
    # to memory.
    if params['stride'] > 1 or params['starting_frame'] > 0:
        sel = traj.select_atoms('all')
        analyze = AnalysisFromFunction(
            lambda ag: ag.positions.copy(), traj.atoms
        )
        coor = analyze.run().results
        traj = MDAnalysis.Merge(sel)
        traj = traj.load_new(
            coor[params['starting_frame']::params['stride']],
            format=MemoryReader
        )

    # Calculate trajectory on principal components and retrieve PCA
    # information
    vect_components, PCA_coeff, coords_avg_atoms, coords_first_frame, data_json_info = get_PCA_trajectory(
        traj, params['selection'], params['cum_var']
    )

    # Keep only the same number of top vect_components as the number of
    # PCA_coeff
    vect_components = vect_components[:len(PCA_coeff[0])]

    # Compress information further by rounding and converting to ints. Save to
    # dictionary to be outputted as a json.
    my_dict = {}
    my_dict['vecs'] = _Utils.compress_list(vect_components, num_decimals)
    my_dict['coeffs'] = _Utils.compress_list(PCA_coeff, num_decimals)
    my_dict['coors'] = _Utils.compress_list(
        coords_avg_atoms, num_decimals, is_coor=True
    )
    my_dict['first_coors'] = _Utils.compress_list(
        coords_first_frame, num_decimals, is_coor=True
    )
    my_dict['params'] = params
    my_dict['res_info'] = data_json_info

    # Write final json file.
    json_txt = json.dumps(my_dict) # , sort_keys=True, indent=4)

    # Remove spaces to compress further.
    json_txt = json_txt.replace(' ', '')

    # Obtain output file name.
    output_name = _Utils.output_filename(
        'compressed', 'json', coor_file=coor_file,
        output_dir=params['output_dir']
    )

    # Write the file.
    with open(output_name, 'w') as write_file:
        write_file.write(json_txt)

    # Test it
    # un_pca(vect_components, PCA_coeff, coords_first_frame, "test.pdb") # output_name)
