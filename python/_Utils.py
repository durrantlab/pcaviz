import os
import glob

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


def compress_list(to_round, num_decimals, is_coor=False):
    # This function compresses a list by rounding its elements and
    # losslessly converting them to integers
    rounded_list = np.around(to_round, decimals=num_decimals)
    multiplication_factor = 10 ** num_decimals
    rounded_list *= multiplication_factor
    compressed_list = np.array(rounded_list, dtype=np.int32)
    if is_coor:
        compressed_list = compressed_list.flatten().tolist()
    else:
        compressed_list = compressed_list.tolist()
    return compressed_list

# Prunes trajectory to a given number of frames.
# Trims the beginning of a trajectory keeping a given percentage.
def stride_and_trim_traj(traj, stride, start_frame, ps_timestep):
    
    # Make sure we are affecting every atom in the trajectory.
    sel = traj.select_atoms("all")

    # Find number of frames in trajectory
    num_frames = len(traj.trajectory)

    # Prune and trim as necessary. This also transfers the trajectory into memory.
    traj.transfer_to_memory(start, num_frames, stride)

    # Important for further RMSD analysis.
    print("If old dt was " + str(ps_timestep) + " ps, new dt is " + \
          str((ps_timestep * stride)) + " ps.")

    return traj
