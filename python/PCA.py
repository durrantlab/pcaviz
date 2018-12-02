import os.path
import copy
import json
import pandas as pd

import MDAnalysis
from MDAnalysis.analysis import pca
import numpy as np

import _Utils
from _UserParams import get_params, log_params

# Get params so we know if we are being run on a server.
params = get_params('PCA')

# We will always iterate, thus coor_file must be a list.
if params['coor_file'] is None:
    params['coor_file'] = [None]

# Get user parameters.
disp_frame = params['disp_frame']

# Create variables to allow breaking the trajectory loop later on.
trajs_projected = []
trajs = []

# Create variables for each first trajectory specific operation.
first_traj = None
first_frame = None

# Perform all PCA analyses for each provided coordinate file.
# Begin trajectory loop.
for coor_id, coor_file in enumerate(params['coor_file']):

    # Load a trajectory for coor_file, aligning it if it is the first.
    if first_traj is None:
        if coor_file is not None:
            first_traj = _Utils.load_traj(params['top_file'],
                                          coor_file,
                                          align_sel=params['align_sel'],
                                          ps_timestep=params['ps_timestep'])
        else:
            first_traj = _Utils.load_traj(params['top_file'],
                                          align_sel=params['align_sel'],
                                          ps_timestep=params['ps_timestep'])
        traj = first_traj
    else:
        traj = _Utils.load_traj(params['top_file'],
                                coor_file,
                                align=False,
                                ps_timestep=params['ps_timestep'])

    # Make sure provided frames exist in the trajectory.
    if np.amax(disp_frame) > len(traj.trajectory):
        raise IndexError("At least one of the frames you provided " + \
                         "to be displayed exceeds the length of trajectory " + \
                         str(coor_id))

    # Save the first frame of the primary, aligned trajectory. You'll use this
    # to align the additional trajectories (if any) in a bit.
    if first_frame is None:
        first_frame = _Utils.output_filename("_tmp.primary_traj.first_frame",
                                             "pdb",
                                             dont_randomize_output=True,
                                             output_dir=("." + os.sep),
                                             coor_file="",
                                             id=params['id'])
        _Utils.save_pdb_first_frame(first_traj,
                                             first_frame,
                                             params["align_sel"])
    else:
        # Align any additional trajectory to the first frame of the first.
        traj = _Utils.align_trajectory(traj,
                                                sel_str=params["align_sel"],
                                                pdb_filename=first_frame)
    trajs.append(traj)

# JDD: For all numerical lists (coors, vectors, etc.), several things should
# happen.
# 1. They should be rounded to a user-specified precision
# 2. They should be multiplied to eliminate the decimal (to compress). So if
#    the precision is 2, they would be multiplied by 100.
# Since these things need to happen to every list, it would be good to make a
# definition (function) dedicated to these tasks.

###
# Calculate pca of the first trajectory using MDAnalysis' class.
def PCA_variance(first_traj, params):
    if params["cum_var"]:
        # JDD: GOOD TO PUT COMMENTS HERE SAYING WHAT EACH STEP DOES.
        get_trajectory = first_traj
        get_trajectory_pca = pca.PCA(get_trajectory, params["selection"] )
        pca_space = get_trajectory_pca.run()
        pca_vectors = pca_space.p_components #vectors components
        n_pcs = np.where(get_trajectory_pca.cumulated_variance > params["cum_var"])[0][0]
        atomgroup = get_trajectory.select_atoms(params["selection"])
        coords_atoms = atomgroup.positions
        coords_project_onto_pca_space = get_trajectory_pca.transform(atomgroup,n_components=n_pcs) # coefficients
    else:
        print("JDD: WHAT HAPPENS IF cum_var NOT SPECIFIED? SHOULD BE REQUIRED.")
    return pca_vectors, coords_project_onto_pca_space, coords_atoms

# Calls PCA variance function
if first_traj is not None:
    # gets pca vector and coefficients from PCA_variance function
    vect_components, PCA_coeff, coords_atoms = PCA_variance(first_traj, params)
    #rounds and turns both PCA vectors and coefficients to 2 sig figs after decimal
    rounded_vect = np.around(vect_components, decimals=2)  # JDD: decimals is fixed at 2. Should be user parameter.
    # JDD: GOOD TO PUT MORE COMMENTS HERE TOO.
    rounded_coefficient = np.around(PCA_coeff, decimals=2)
    vect_list = np.array(rounded_vect).tolist() #PCA vector component as a list
    coefficient_list = np.array(rounded_coefficient).tolist() #PCA coefficients as a rounded list
    # JDD: WOULD BE GOOD IF ATOM_COORDINATES WERE ROUNDED
    atom_coordinates = np.array(coords_atoms).tolist()

# If we were using a dual function file, set coor_file to be the top_file.
if params['coor_file'][0] is None:
    params['coor_file'] = [params['top_file']]

# Remove first frame used for alignment.
_Utils.del_filename(first_frame)

#Creating a dictionary that contains PCA vectors, coefficients, and list of parameters
my_dict = {}
my_dict["vecs"] = vect_list  # JDD: I simplified the key names.
my_dict["coeffs"] = coefficient_list
my_dict["params"] = params
my_dict["coors"] = atom_coordinates

# writes final json file
json_txt = json.dumps(my_dict) # , sort_keys=True, indent=4)

# JDD: I removed spaces from the json_txt file to make it smaller
json_txt = json_txt.replace(" ", "")

with open("PCA_file.json", "w") as write_file:
    write_file.write(json_txt)
