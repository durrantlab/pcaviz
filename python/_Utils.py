import os
import glob

import matplotlib.pyplot as plt
from matplotlib.pyplot import figure, draw
import numpy as np

import MDAnalysis
from MDAnalysis.analysis import align


def output_filename(description, ext, coor_file, id,
                    dont_randomize_output=False, output_dir=None):
    # This function generates an output filename, respecting user parameters in
    # params. If temp_file is true, load a mini param. If it's None, use the
    # global params, unless params_to_use specified. If that's specified, use
    # it.
    
    # If the description is present, add "." to the beginning.
    if description != "":
        description = "." + description

    # Pick a random id, if any.
    if not dont_randomize_output:
        run_id_to_use = "." + id
    else:
        run_id_to_use = ""

    # Assemble the filename.
    filename = coor_file + description + run_id_to_use + "." + ext

    # Put it in the right path.
    if output_dir is not None: 
        filename = output_dir + os.path.basename(filename)
    
    print "Writing to file: " + filename

    return filename


def del_filename(filename):
    # Deletes the filename with path filename.
    print "Deleting file: " + filename
    os.unlink(filename)


def find_closest_points(data, codebook, labels):
    # This function takes in a codebook and labels array and identifies the
    # points closest to the codebook values within data. It returns an array
    # of the indices of the points in data and distances from closest cluster
    # center, as well as an array of the points themselves. Labels must be the
    # same length as data, and codebook must have the same width as data.

    # Prepare arrays for storing results.
    closest_frames = np.zeros((len(codebook), 2))
    closest_points = np.ndarray(closest_frames.shape)

    # Set the first frame as the closest point.
    for clust_index in range(0, len(codebook)):
        closest_frames[clust_index,1] = np.sqrt(
                            np.square(data[0,0] - codebook[clust_index,0]) + \
                            np.square(data[0,1] - codebook[clust_index,1]))
    
    # For every frame, find its distance from the closest codebook value and
    # keep track of the closest one
    for index, clust_index in enumerate(labels):
        dist = np.sqrt(np.square(data[index,0] - codebook[clust_index,0]) + \
                       np.square(data[index,1] - codebook[clust_index,1]))
        if dist < closest_frames[clust_index,1]:
            closest_frames[clust_index,0] = index
            closest_frames[clust_index,1] = dist
            closest_points[clust_index] = data[index]
    return closest_frames, closest_points


def begin_plot(title, ylabel, xlabel, figure_number=1):
    # This functions prepares plots with a title and axes labels.
    plt.figure(figure_number,figsize=[10,7])
    plt.title(title)
    plt.ylabel(ylabel)
    plt.xlabel(xlabel)


def find_top_and_coor_files():
    # This function finds and returns topology and coordinate files in the
    # working directory. Caution is advised as some MDAnalysis compatible
    # extensions are relatively common ones. Thus, please ensure your
    # working directory is clean and contains only the files you wish to use.

    # List of all MDAnalysis compatible dual function file formats.
    extensions_both = ('*.pdb', '*.data', '*.xyz', '*.txyz', '*.arc',
                        '*.gsd', '*.pdb/ent', '*.pdbqt', '*.pqr', '*.gro',
                        '*.crd', '*.dms', '*.config', '*.history', '*.mmtf',
                        '*.mol2', '*.tpr', '*.log')

    # Find first dual function file in working directory.
    list_of_both_files = []
    for ext in extensions_both:
        list_of_both_files.extend(glob.glob(ext))

    if len(list_of_both_files) > 0:
        top_file = list_of_both_files[0]
        coor_file = None
        if len(list_of_both_files) > 1:
            print "Warning: Multiple topology/coordinate files found. " +\
                  "Using first one found in directory: " + top_file + "."
    else:
        # List of all MDAnalysis compatible topology file formats,
        # excluding common extension names and formats that act as
        # both topology andcoordinate files.
        extensions_top = ('*.psf', '*.top', '*.prmtop', '*.parm7', '*.xml')

        # Find first topology file in working directory.
        list_of_top_files = []
        for ext in extensions_top:
            list_of_top_files.extend(glob.glob(ext))
        if len(list_of_top_files) > 0:
            top_file = list_of_top_files[0]
            if len(list_of_top_files) > 1:
                print "Warning: Multiple topology files found. " +\
                      "Using first one found in directory: " + top_file + "."

        # List of all MDAnalysis compatible coordinate file formats,
        # excluding common extension names and formats that act as
        # both topology and coordinate files.
        extensions = ('*.dcd', '*.xtc', '*.trr', '*.trj', '*.out', '*.trz',
                        '*.mdcrd', '*.inpcrd', '*.restrt', '*.ncdf', '*.nc')

        # Find first coordinate file in working directory.
        list_of_coor_files = []
        for ext in extensions:
            list_of_coor_files.extend(glob.glob(ext))
        if len(list_of_coor_files) > 0:
            coor_file = list_of_coor_files[0]
            if len(list_of_coor_files) > 1:
                print "Warning: Multiple coordinate files found. " +\
                      "Using first one found in directory: " + coor_file + "."

    return top_file, coor_file


def load_traj(top_file, coor_file=None, align=True, 
              align_sel='name CA', ps_timestep=.002):
    # Read in coordinate file and topology file as a universe.
    if coor_file is not None:
        traj = MDAnalysis.Universe(top_file, coor_file, dt=ps_timestep)
        traj.trajectory[0::100]
    else:
        traj = MDAnalysis.Universe(top_file, dt=ps_timestep)
        traj.trajectory[0::100]

    # Align the trajectory if requested.
    if align:
        traj = align_trajectory(traj, sel_str=align_sel)

    return traj


def check_axes(x_min, x_max, y_min, y_max):
    # Makes sure that if either x_min or x_max is specified,
    # the other is as well. Also checks y_min and y_max.
    # Inform user if one axis is specified while the other is not.
    # Make sure min values are less than max values.
    if (x_min is not None and x_max is None) or \
       (x_min is None and x_max is not None):
        raise ValueError('ERROR: If you specify --axis_x_min, you must also specify ',
                    '--axis_x_max, and visa versa.')
    if (y_min is not None and y_max is None) or \
       (y_min is None and y_max is not None):
        raise ValueError('ERROR: If you specify --axis_y_min, you must also specify ',
                    '--axis_y_max, and visa versa.')
    if y_min is not None:
        if y_min >= y_max:
            raise ValueError('ERROR: Please make sure your --axis_y_max value is larger ',
                            'than your --axis_y_min value.')
    if x_min is not None:
        if x_min >= x_max:
            raise ValueError('ERROR: Please make sure your --axis_x_max value is larger ',
                            'than your --axis_x_min value.')


def _save_dcd(traj, filename, sel_str):
    # This function saves a MDAnalysis universe, or a subset of atoms
    # therein, as a new file.

    sel = traj.select_atoms(sel_str)
    with MDAnalysis.Writer(filename,
                           multiframe=True,
                           n_atoms=traj.atoms.n_atoms) as W:
        for ts in traj.trajectory:
            W.write(sel)


def save_pdb_first_frame(traj, filename, sel_str):
    # This function saves the first frame of a MDAnalysis universe,
    # or a subset of atoms therein, as a new file.

    sel = traj.select_atoms(sel_str)
    with MDAnalysis.Writer(filename,
                           multiframe=True,
                           n_atoms=traj.atoms.n_atoms) as W:
        traj.trajectory[0]
        W.write(sel)


def align_trajectory(traj, sel_str="name CA", pdb_filename=None):
    # This function aligns a trajectory to its first frame.

    # If we are given a pdb to align to, load it as a universe
    # and align to it.
    if pdb_filename is not None:
        ref = MDAnalysis.Universe(pdb_filename)
        sel_ref = ref.select_atoms(sel_str)
        alignment = align.AlignTraj(traj, ref, in_memory=True, select=sel_str)
    else:
        # Get the reference (first frame).
        traj.trajectory[0]
        
        # Align the trajectory to current frame of trajectory
        # (in this case the first). Note that this is in memory.
        # Could run into problems with big trajectories...
        alignment = align.AlignTraj(traj, traj, in_memory=True, select=sel_str)
        alignment.run()

    # Aligned trajectory is in traj.
    return traj

