import os
import glob

import numpy as np
import MDAnalysis
from MDAnalysis.analysis import align


def output_filename(description, ext, coor_file, output_dir=None):
    """
    This function generates an output filename of type ext.

    The file is saved to output_dir if provided. Defaults to the coor_file's directory.

    :param str description: A decription of what the file is.

    :param str ext: The file extension to attach to the output file.

    :param str coor_file: The coordinate file being handled.

    :param str output_dir: The directory to which the file should be saved.

    :return: A path and filename ordered as output_dir/coor_name.description.ext
    :rtype: :class:'str'
    """
    
    # If the description is present, add "." to the beginning.
    if description != "":
        description = "." + description

    # Assemble the filename.
    filename = os.path.splitext(coor_file)[0] + description + "." + ext

    # Put it in the right path.
    if output_dir is not None: 
        filename = output_dir + os.path.basename(filename)
   
    # Let the user know where it is going.
    print "Writing to file: " + filename
    
    # Return the filename.
    return filename


def find_top_and_coor_files():
    """
    This function finds and returns topology and coordinate files in the
    working directory. Caution is advised as some MDAnalysis compatible
    extensions are relatively common ones. Thus, please ensure your
    working directory is clean and contains only the files you wish to use.

    Priority is given to files that can serve as both topology and coordinate
    files. For all files, the first found within the directory will always be used.

    For more information on which files will be used, please see MDAnalysis' list of
    supported file readers.

    :return: The filenames of any topology/coordinate files found within the current working directory.
             If a dual function file is found, it is returned within top_file and coor_file is set to None.
    :rtype: :class:'str'
    """

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
              align_sel='name CA'):
    """
    Loads an MDAnalysis Universe from topology and/or coordinate files. Requires MDAnalysis.

    Also calls align_trajectory on the universe according to an atom selection if specified.

    :param str top_file: Filename of topology file to be used. If a dual function file is being used
                         it must be passed through top_file.

    :param str coor_file: Filename of coordinate file to be used. Defaults to None.

    :param bool align: If True, Universe is aligned. If False, Universe is not aligned. Defaults to True.

    :param str align_sel: Atom selection to use when aligning. Defaults to 'name CA'.

    :return: The Mdanalysis Universe constructed from provided files, potentially aligned to its first frame.
    :rtype: :class:'MDAnalysis.Universe'
    """

    # Read in coordinate file and topology file as a universe.
    if coor_file is not None:
        traj = MDAnalysis.Universe(top_file, coor_file)
        traj.trajectory[0::100]
    else:
        traj = MDAnalysis.Universe(top_file)
        traj.trajectory[0::100]

    # Align the trajectory if requested.
    if align:
        traj = align_trajectory(traj, align_sel=align_sel)

    return traj


def align_trajectory(traj, align_sel="name CA"):
    """
    This function aligns a trajectory to its first frame according to an atom selection.

    Uses MDAnalysis' align module.

    :param MDAnalysis.Universe traj: The universe containing the trajectory to be aligned.

    :param str sel_str: Atom selection to use when aligning. Defaults to 'name_CA'.

    :return: A new MDAnalysis Universe aligned to its first frame.
    :rtype: :class:'MDAnalysis.UNiverse'
    """

    # Get the reference (first frame).
    traj.trajectory[0]

    # Align the trajectory to current frame of trajectory
    # (in this case the first).
    alignment = align.AlignTraj(traj, traj, in_memory=True, select=align_sel)
    alignment.run()

    # Aligned trajectory is in traj.
    return traj


def compress_list(to_round, num_decimals, is_coor=False):
    """
    This function compresses a list by rounding its elements to num_decimals places.
    
    It then converts the elements to integers and flattens the output according to is_coor.

    :param array_like to_round: An array_like containing floats to be compressed.

    :param int num_decimals: How many decimal places to round to.

    :param bool is_coor: If True, output list will be flattened. If False list will retain its structure. Defaults to False.

    :return: A basic Python list of ints with shape identical to to_round, or flattened if is_coor set to True.
    :rtype: :class:'list'
    """

    # First round to specified number of decimal places.
    rounded_list = np.around(to_round, decimals=num_decimals)

    # Prepare for int conversion by moving decimals.
    multiplication_factor = 10 ** num_decimals
    rounded_list *= multiplication_factor

    # Convert to integers.
    compressed_list = np.array(rounded_list, dtype=np.int32)

    # We can flatten the coordinates to save a bit more space.
    # Only atomic coordinates can be flattened since other lists'
    # shapes cannot be inferred.
    if is_coor:
        compressed_list = compressed_list.flatten().tolist()
    else:
        compressed_list = compressed_list.tolist()

    # Return compressed list.
    return compressed_list

