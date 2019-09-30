#!/usr/bin/env python

import json
import csv
import os
import sys
import shutil
import argparse

import MDAnalysis
from MDAnalysis.analysis import rms
from MDAnalysis.analysis.base import AnalysisFromFunction
from MDAnalysis.coordinates.memory import MemoryReader
import numpy as np
from sklearn.decomposition import PCA

__version__ = "1.2"
FROM_PIP = False

class Utils:
    """A class (namespace) that contains utiliity definitions."""

    @staticmethod
    def get_output_filename(description, ext, coor_file, output_dir=None):
        """This definition generates an output filename of type ext.

        :param description: A decription of the file.
        :type description:  str
        :param ext:         The file extension to attach to the output file.
        :type ext:          str
        :param coor_file:   The coordinate file name being handled.
        :type coor_file:    str
        :param output_dir:  The directory to which the file should be saved,
                            defaults to the coor_file's directory.
        :type output_dir:   string, optional
        :return:            A path and filename ordered as
                            output_dir/coor_name.description.ext
        :rtype: str
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


    @staticmethod
    def load_traj(top_file, coor_file=None):
        """Loads an MDAnalysis Universe from topology and/or coordinate files.
        Requires MDAnalysis.

        :param top_file:  Filename of topology file to be used. If a dual
                          function file is being used it must be passed through
                          top_file.
        :type top_file:   str
        :param coor_file: Filename of coordinate file to be used. Defaults to
                          top_file.
        :type coor_file:  str, optional
        :return:          The MDAnalysis Universe built from the provided
                          files.
        :rtype:           MDAnalysis.Universe
        """

        # Read in coordinate file and topology file as a universe.
        if coor_file is not None:
            traj = MDAnalysis.Universe(top_file, coor_file)
            # Not sure why 100, but it doesn't change number of frames.
            # traj.trajectory[0::100]
            traj.trajectory[:]
        else:
            traj = MDAnalysis.Universe(top_file)
            # traj.trajectory[0::100]
            traj.trajectory[:]

        return traj


    @staticmethod
    def compress_list(to_round, precision, flatten=False):
        """This definition compresses a list by rounding its elements. It then
        converts the elements to integers and optionally flattens the output.

        :param to_round:  An array-like variable containing floats to be
                          compressed.
        :type to_round:   list
        :param precision: How many decimal places to round to.
        :type precision:  int
        :param flatten:   If True, the output list will be flattened. If False
                          list will retain its 2D structure. Defaults to False.
        :type flatten:    bool, optional
        :return:          A basic Python list of ints with shape identical to
                          to_round, or flattened if flatten is set to True.
        :rtype: list
        """

        # First round to specified number of decimal places.
        rounded_list = np.around(to_round, decimals=precision)

        # Prepare for int conversion by moving decimals.
        multiplication_factor = 10 ** precision
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


class PCAFuncs:
    """A class (namespace) that includes definitions related to PCA."""

    @staticmethod
    def compress_PCA(traj, selection, cum_var):
        """This definition performs PCA analysis using the MDAnalysis PCA
        analysis class. It saves the information necessary for expansion back
        into a trajectory.

        :param traj:      The universe containing the trajectory to be compressed.
        :type traj:       MDAnalysis.Universe
        :param selection: The atom selection to be included in the compression
                          e.g. protein, name CA...
        :type selection:  str
        :param cum_var:   The amount of variance to be explained cumulatively
                          by calculated components.
        :type cum_var:    float
        :return:          The principal component vectors.
        :rtype:           ndarray
        :returns:         The coordinates of the selected trajectory atoms
                          projected onto the calculated components.
        :rtype:           ndarray
        :returns:         The coordinates of the average positions of the
                          selected atoms in Cartesian space.
        :rtype:           ndarray
        :returns:         The first-frame coordinates of the selected
                          trajectory atoms in Cartesian space.
        :rtype:           ndarray
        :returns:         The atom type and number of the selected atoms.
        :rtype            list
        :returns:         The residue number and name of the selected atoms.
        :rtype            list
        """

        # Get the indices of the atoms that match the selection.
        atom_group = traj.select_atoms(selection)
        idx_of_selected_atoms = atom_group.indices

        # Get the flattened coordinates of each frame.
        coor_data = []
        for frame in traj.trajectory:
            coor_data.append(frame.positions[idx_of_selected_atoms].ravel())
        coor_data = np.array(coor_data)

        # Append coordinates to python list then use ravel to flatten
        # coordinates without having to allocate new memory

        # Turning coord list into numpy array because it reduces memory
        # consumption You get a flatten numpy array containing the coordinates
        pca = PCA(n_components=min(coor_data.shape))
        pca.fit(coor_data)
        pca_vectors = pca.components_

        # Number of components should be equal to number of elements in
        # coor_data array. Number of principal components = 3xs the number of
        # atoms. pca.fit fits model to coordinate data. PCA vectors contain
        # principal axes in PCA space, show max variance in data.

        # Also get the coefficients for each frame.
        coords_project_onto_pca_space = pca.transform(coor_data)

        # coords_project_onto_pca_space applies PCA dimensionality reduction
        # to coor_data

        # Calculate how many components are required to cumulatively explain
        # the desired variance.
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

            # If statement will check to see if residue id is present if it's
            # not present it will add the residue id, residue name, and the
            # residue's atoms.
            key = str(res_id) + res_name # Because list is not hashable.
            if key != last_key:
                data_json.append([int(res_id), res_name, at_name])
                last_key = key
            else:
                data_json.append(at_name)

        # Return only the information necessary for compression and expansion.
        return pca_vectors, coords_project_onto_pca_space, coords_avg_atoms, \
            coords_first_frame, data_json


    @staticmethod
    def expand_PCA(pca_vectors, coords_project_onto_pca_space, coords_avg_atoms, precision):
        """This definition decompresses a trajectory that was compressed with
        the compress_PCA definition. It uses only the information necessary for
        expansion back into Cartesian-space and returns a list of the
        decompressed atomic coordinates.

        :param pca_vectors:                   The principal components onto
                                              which the trajectory has been
                                              projected.
        :type pca_vectors:                    list
        :param coords_project_onto_pca_space: The projected trajectory
                                              coordinates in PCA-space.
        :type coords_project_onto_pca_space:  list
        :param coords_avg_atoms:              The average atomic coordinates of
                                              the original trajectory.
        :type coords_avg_atoms:               list
        :param precision:                     The number of decimal places
                                              after the PCA vectors,
                                              coefficients, and atomic
                                              coordinates were rounded.
        :type precision:                      int
        :return:                              The decompressed trajectory. A
                                              list containing ndarrays of the
                                              atomic coordinates at each frame,
                                              in the same format as an ordinary
                                              MDAnalysis Universe
                                              atom_selection.positions field.
        :rtype: list
        """

        pca_vectors = np.true_divide(pca_vectors, (10**precision))
        coords_project_onto_pca_space = np.true_divide(coords_project_onto_pca_space,
                                                    (10**precision))
        coords_avg_atoms = np.true_divide(coords_avg_atoms, (10**precision))


        # Prepare list to hold expanded trajectory coordinates.
        expanded_coordinates = []

        # Construct each frame's coordinates one at a time.
        for frame_num, frame_coefficients in enumerate(coords_project_onto_pca_space):

            # Calculate flattened XYZ coordinate delta values.
            XYZ_coor_delta = np.dot(frame_coefficients, pca_vectors)

            # Add the average atomic positions.
            XYZ_coor = XYZ_coor_delta + coords_avg_atoms

            expanded_coordinates.append(np.reshape(XYZ_coor, (len(XYZ_coor) // 3, 3)))

        # Return the expanded coordinates.
        return expanded_coordinates


    @staticmethod
    def save_uncompressed_pdb(coors, atom_inf, output_filename):
        """Saves an XYZ file of the uncompressed trajectory.

        :param coors:           A list of numpy.array coordinates.
        :type coors:            list
        :param atom_inf:        The residue index, name, and atom name of each
                                atom.
        :type atom_inf:         list
        :param output_filename: The filename of the output xyz file.
        :type output_filename:  str
        """

        atom_inf2 = []
        last_resindex = None
        last_resname = None
        for ai in atom_inf:
            if type(ai) == str:
                atom_inf2.append([last_resindex, last_resname, ai])
            else:
                atom_inf2.append(ai)
                last_resindex = ai[0]
                last_resname = ai[1]

        with open(output_filename, "w") as f:
            coor_len = len(coors[0])
            for i, coor in enumerate(coors):
                f.write(str(coor_len) + "\n")
                f.write(" pcaviz\n")
                for c, ai in zip(coor, atom_inf2):
                    x, y, z = c
                    residx, resname, atomname = ai
                    f.write(atomname + " " + str(x) + " " + str(y) + " " + str(z) + "\n")

class Tests:
    """A class/namespace that contains definitions related to running PCACViz
    in test mode."""

    @staticmethod
    def start():
        """Initial preparation step required to run PCAViz in testing mode.

        :return: The parameters to use for testing.
        :rtype: dict
        """

        output_dir = "tmp_testing" + os.sep
        print("Running in test mode. All other parameters will be ignored.")
        if os.path.exists(output_dir):
            print("Directory already exists: " + output_dir + "! Delete this directory before running the test.")
            sys.exit(0)
        os.mkdir(output_dir)
        examples_dir = os.path.dirname(os.path.realpath(__file__)) + os.sep + "examples" + os.sep
        params = {
            "top_file": examples_dir + "1J8K_example.pdb",
            "coor_file": examples_dir + "1J8K_example.pdb",
            "selection": "name C N CA",
            "output_dir": output_dir,
            "stride": 2,
            "cum_var": 0.8,
            "precision": 1,
            "check_accuracy": True,
            "test": True
        }
        return params


    @staticmethod
    def finish():
        """Runs after tests are complete. Checks if output values are as
        expected, cleans up temporary files, etc."""

        print("")
        print("Checking if calculated mean RMSD is 0.53...")
        mean = [l for l in open("tmp_testing" + os.sep + "1J8K_example.starting_vs_expanded_RMSD.csv")
                if "mean" in l][0].split(",")[1][:4]
        if (mean == "0.53"):
            print("\tPassed!")
        else:
            print("\tFailed!")
            shutil.rmtree('tmp_testing')
            sys.exit(0)

        PCA_inf = json.load(open("tmp_testing" + os.sep + "1J8K_example.compressed.json"))
        print("Checking 846 coordinates (282 atoms * 3)...")
        if len(PCA_inf["coors"]) == 846 and len(PCA_inf["first_coors"]) == 846:
            print("\tPassed!")
        else:
            print("\tFailed!")
            shutil.rmtree('tmp_testing')
            sys.exit(0)

        print("Checking 3 frames (5 original frames, stride 2)...")
        if len(PCA_inf["coeffs"]) == 3:
            print("\tPassed!")
        else:
            print("\tFailed!")
            shutil.rmtree('tmp_testing')
            sys.exit(0)

        print("Checking 2 principal components...")
        if len(PCA_inf["vecs"]) == 2:
            print("\tPassed!")
        else:
            print("\tFailed!")
            shutil.rmtree('tmp_testing')
            sys.exit(0)

        shutil.rmtree('tmp_testing')

    @staticmethod
    def check_accuracy(params, vect_components, PCA_coeff, coords_avg_atoms,
                       precision, coor_file, data_json_info, traj):
        """Checks the accuracy of a PCAViz-processed trajectory.

        :param params:           The user parameters.
        :type params:            dict
        :param vect_components:  The PCAViz-generated principal components.
        :type vect_components:   list
        :param PCA_coeff:        The PCA coefficients for each frame.
        :type PCA_coeff:         list
        :param coords_avg_atoms: The average atomic coordinates of the
                                 trajectory.
        :type coords_avg_atoms:  list
        :param precision:        The number of decimal places after the PCA
                                 vectors, coefficients, and atomic coordinates
                                 were rounded.
        :type precision:         int
        :param coor_file:        The coordinate filename.
        :type coor_file:         str
        :param data_json_info:   The residue index, name, and atom name of each
                                 atom.
        :type data_json_info:    dict
        :param traj:             The original trajectory.
        :type traj:              MDAnalysis.Universe
        """

        # Prepare list to hold RMSD values.
        RMSD_lst = []

        # Calculate expanded coordinates.
        expanded_trajectory = PCAFuncs.expand_PCA(vect_components, PCA_coeff,
                                                  coords_avg_atoms, precision)

        xyz_filename = Utils.get_output_filename('uncompressed',
                                                 'xyz',
                                                 coor_file=coor_file,
                                                 output_dir=params['output_dir'])
        PCAFuncs.save_uncompressed_pdb(expanded_trajectory, data_json_info, xyz_filename)

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
        print("Mean RMSD between original and expanded trajectories: " + \
                str(mean_RMSD))
        print("Standard deviation of RMSD between original and " + \
                "expanded trajectories: " + str(std_RMSD))


        print("Number top components required to account for >= " + \
              str(100 * params["cum_var"]) + "% of " + \
              "variance: " + str(len(vect_components)))

        # Save a table containing results.
        tbl_flnm = Utils.get_output_filename('starting_vs_expanded_RMSD',
                                    'csv',
                                    coor_file=coor_file,
                                    output_dir=params['output_dir'])
        with open(tbl_flnm, 'w') as tf:
            table_writer = csv.writer(tf)
            table_writer.writerows(zip(frame_numbers, RMSD_lst))

        # Report file size too
        json_file = tbl_flnm.replace("starting_vs_expanded_RMSD.csv", "compressed.json")
        import os
        print("Size of output JSON file: " + \
            str(
                np.round(
                    os.stat(json_file).st_size / 1024.0,  # Using this instead of 1000... definition is ambiguous
                    3
                )
            ) + "K"
        )

        # Temporarily compress the file to report that size as well.
        import gzip
        import shutil
        with open(json_file, 'rb') as f_in:
            with gzip.open(json_file + ".tmp.gz", 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        print("Approximate size if gz compressed: " + \
            str(
                np.round(
                    os.stat(json_file + ".tmp.gz").st_size / 1024.0,  # Using this instead of 1000... definition is ambiguous
                    3
                )
            ) + "K"
        )
        os.remove(json_file + ".tmp.gz")


def get_params():
    """This definition handles user command-line input. If the first argument
    is a json file, it will read the parameters from that file. Otherwise the
    parameters must be specified from the command line. Please use the -h flag
    for more information on available parameters.

    :return: A dictionary containing all the user-provided parameters as values
             associated with their respective keys.
    :rtype:  dict

    ::
        >>> python PCAViz.py -h
    """

    if FROM_PIP:
        exec_str = "pcaviz"
    else:
        exec_str = "python PCAViz.py"

    # If no parameters given, show the help file.
    if len(sys.argv) == 1:
        sys.argv.append("-h")

    # If a json file is provided,
    # read in parameters from json file.
    possible_json_file = sys.argv[1].lower()
    if os.path.splitext(possible_json_file)[1][1:] == "json":
        with open(sys.argv[1], 'r') as param_file:
            params = json.load(param_file)
    else:
        # Otherwise parse the parameters from the system arguments.
        # parser = argparse.ArgumentParser(description='Compress an MD trajectory, saving the output to a JSON file.')

        PARSER = argparse.ArgumentParser(
            formatter_class=argparse.RawDescriptionHelpFormatter,
            description="""
PCAViz """ + __version__ + """ is an open-source Python/JavaScript toolkit for sharing and
visualizing MD trajectories via a web browser. This Python script, called the
PCAViz Compressor, compresses an MD trajectory and saves the output to a JSON
file.""",
            epilog="""
EXAMPLES OF USE:

1. Create a compressed JSON file from a topology (1J8K_example.psf) and a
   trajectory/coordinate (1J8K_example.dcd) file.

""" + exec_str + """ --top_file examples/1J8K_example.psf --coor_file examples/1J8K_example.dcd

2. PDB files can also contain multiple frames. In this case, the same file
   serves as the topology and trajectory file.

""" + exec_str + """ --top_file examples/1J8K_example.pdb --coor_file examples/1J8K_example.pdb

3. By default, PCAViz includes only the backbone atoms in the output. These
   atoms should be enough for cartoon-style visualization. But you can select
   your own atoms to include in the output. See https://goo.gl/kVeQuN to learn
   how to construct an atom-selection string.

""" + exec_str + """ --top_file examples/1J8K_example.pdb --coor_file examples/1J8K_example.pdb --selection "name *"

4. Striding the trajectory frames can reduce file sizes. PCAViz will
   interpolate between the remaining frames to fill in the frames that are
   missing. Here we keep only every other frame:

""" + exec_str + """ --top_file examples/1J8K_example.pdb --coor_file examples/1J8K_example.pdb --stride 2

5. PCAViz allows users to control the compression settings. Two settings are
   available. First, the user can specify how much of the cumulative variance
   should be explained by the principal components. PCAViz will pick the
   number of components required to meet that goal. Second, the user can
   specify how precisely PCAViz should represent numeric values (e.g.,
   principal-component coefficients). For example, here we tell PCAViz to
   produce a trajectory that accounts for 80% of the variance, and to round all
   numbers in the output JSON file to the nearest hundredth (two decimal
   places):

""" + exec_str + """ --top_file examples/1J8K_example.pdb --coor_file examples/1J8K_example.pdb --cum_var 0.8 --precision 2

6. To find the ideal --cum_var and --precision parameters, you may wish to
   check how closely the PCAViz-compressed trajectory matches the original
   trajectory. You can instruct PCAViz to output a CSV file that provides a
   frame-by-frame RMSD comparison between the two trajectories. This option
   also outputs an XYZ trajectory file that you can visually compare to the
   original.

""" + exec_str + """ --top_file examples/1J8K_example.pdb --coor_file examples/1J8K_example.pdb --check_accuracy

7. By default, PCAViz saves the compressed JSON file to the same directory
   where the coordinate file is located. You can specify a different output
   directory if needed. The directory will be created if it doesn't exist.

""" + exec_str + """ --top_file examples/1J8K_example.pdb --coor_file examples/1J8K_example.pdb --output_dir "my_dir"

8. For debugging purposes, PCAViz also includes an option to test whether the
   code is fully functional.

""" + exec_str + """ --test

""")

        # First get all parameters.
        PARSER.add_argument('--top_file', metavar='t', type=str, nargs="?",
                            help='The topology filename (e.g., psf).')
        PARSER.add_argument('--coor_file', metavar='c', type=str, nargs="?",
                            help='The coordinate (trajectory) ' +\
                                 'filename (e.g., dcd).')
        PARSER.add_argument('--selection', metavar='s', type=str, nargs="?",
                            default="name CA C N O",
                            help='Which atoms to save to the JSON output ' +\
                                 'file (default: "name CA C N O"). See ' +\
                                 'https://goo.gl/kVeQuN to learn ' +\
                                 'how to construct an atom selection string.')
        PARSER.add_argument('--output_dir', metavar='od', type=str, nargs="?",
                            default=None,
                            help='The directory where output files should be saved. ' +\
                                 'The directory will be created if it does not exist. ' +\
                                 '(default: the directory where the coordinate file is located).')
        PARSER.add_argument('--stride', metavar='ns', type=int,
                            nargs="?", default=1,
                            help='How many frames to stride. For example, stride = 2 ' +\
                                 'means every other frame will be saved, and stride = 3 ' +\
                                 'means every third frame will be saved. (default: 1, no ' +\
                                 'stride).')
        PARSER.add_argument('--cum_var', metavar='var', type=float,
                            nargs='?', default=0.90,
                            help='The target cumulative variance, as a float. PCAViz will ' +\
                                 'use the minimum number of principal components required to ' +\
                                 'capture at least this variance (default: 0.90).')
        PARSER.add_argument('--precision', metavar='p', type=int,
                            nargs='?', default=2,
                            help='The number of decimal places to retain when rounding PCA ' +\
                                 'vectors, coefficients, and atomic coordinates ' +\
                                 '(default: 2, meaning round to the hundredths).')
        PARSER.add_argument('--check_accuracy', action='store_true',
                            help='Create a csv file containing the frame-to-frame RMSDs between ' +\
                                 'original- and decompressed-trajectory frames. Useful for testing ' +\
                                 'the impact of different settings on atom-position accuracy.')
        PARSER.add_argument('--test', action='store_true',
                            help='Tests PCAViz to make sure all components are functioning.')

        # Parse the arguments.
        args = PARSER.parse_args()

        # Ensure output_dir is formatted correctly.
        output_dir = args.output_dir
        if output_dir is not None and not output_dir.endswith(os.sep):
            output_dir = output_dir + os.sep
        if output_dir is not None and not os.path.exists(output_dir):
            os.mkdir(output_dir)

        # I want to convert it to a dictionary. Easier to use.
        params = vars(args)
        params['output_dir'] = output_dir

    # Return the dictionary of parameters.
    return params


def main():
    """The main definition."""
    # Get user provided parameters and unpack frequently used parameters
    params = get_params()

    # Redefine params if running in test mode.
    if params["test"]:
        params = Tests.start()

    coor_file = params['coor_file']
    top_file = params['top_file']
    precision = params['precision']

    # Load a trajectory for input file, aligning it if it is the first.
    traj = Utils.load_traj(top_file, coor_file)

    # If a dual function file is being used, set coor_file to it for filename
    # purposes.
    if coor_file is None:
        coor_file = top_file

    # Stride or trim the trajectory as desired, this transfers the trajectory
    # to memory.
    if params['stride'] > 1:  # or params['starting_frame'] > 0:
        sel = traj.select_atoms('all')
        analyze = AnalysisFromFunction(lambda ag: ag.positions.copy(), traj.atoms)
        coor = analyze.run().results
        traj = MDAnalysis.Merge(sel)
        traj.load_new(coor[0::params['stride']], format=MemoryReader)

    # Calculate trajectory on principal components and retrieve PCA
    # information
    vect_components, PCA_coeff, coords_avg_atoms, coords_first_frame, data_json_info = \
        PCAFuncs.compress_PCA(traj, params['selection'], params['cum_var'])

    # Keep only the same number of top vect_components as the number of PCA_coeff
    vect_components = vect_components[:len(PCA_coeff[0])]

    # Compress information further by rounding and converting to ints.
    vect_components = Utils.compress_list(vect_components, precision)
    PCA_coeff = Utils.compress_list(PCA_coeff, precision)
    coords_avg_atoms = Utils.compress_list(coords_avg_atoms,
                                            precision,
                                            flatten=True)
    coords_first_frame = Utils.compress_list(coords_first_frame,
                                              precision,
                                              flatten=True)

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
    output_name = Utils.get_output_filename('compressed',
                                         'json',
                                         coor_file=coor_file,
                                         output_dir=params['output_dir'])

    # Write the file.
    with open(output_name, 'w') as write_file:
        write_file.write(json_txt)

    # If the user request to check the accuracy of the expanded trajectory.
    if params['check_accuracy']:
        Tests.check_accuracy(params, vect_components, PCA_coeff, coords_avg_atoms,
                             precision, coor_file, data_json_info, traj)

    # If it's running in test mode, finish the test
    if params["test"]:
        Tests.finish()


if __name__ == '__main__':
    main()
