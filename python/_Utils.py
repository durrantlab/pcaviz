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
    print("Writing to file: " + filename)

    # Return the filename.
    return filename


def load_traj(top_file, coor_file=None):
    """
    Loads an MDAnalysis Universe from topology and/or coordinate files. Requires MDAnalysis.

    :param str top_file: Filename of topology file to be used. If a dual function file is being used
                         it must be passed through top_file.

    :param str coor_file: Filename of coordinate file to be used. Defaults to None.

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

    return traj


def compress_list(to_round, num_decimals, flatten=False):
    """
    This function compresses a list by rounding its elements to num_decimals places.

    It then converts the elements to integers and flattens the output according to flatten.

    :param array_like to_round: An array_like containing floats to be compressed.

    :param int num_decimals: How many decimal places to round to.

    :param bool flatten: If True, output list will be flattened. If False list will retain its structure. Defaults to False.

    :return: A basic Python list of ints with shape identical to to_round, or flattened if flatten set to True.
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
    if flatten:
        compressed_list = compressed_list.flatten().tolist()
    else:
        compressed_list = compressed_list.tolist()

    # Return compressed list.
    return compressed_list
