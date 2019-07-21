import sys
import argparse
import random
import os
import json

# from _Utils import find_top_and_coor_files


def get_params():
    """
    This function handles all user input directly using Python's built in argparse functionality.

    It directly accesses the system arguments when PCA.py is called. If the first argument is a json file
    it will read the parameters from that file, otherwise they must be specified directly.
    Please use the -h flag for more information on available parameters.

    :return: A dictionary containing all the user provided parameters as values associated with their respective keys.
    :rtype: :class:'dict'

    ::

        >>> python PCA.py -h
    """

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
        parser = argparse.ArgumentParser(description='Compress an MD trajectory, saving the output to a JSON file.')

        # First get all parameters.
        parser.add_argument('--top_file', metavar='t', type=str, nargs="?",
                            help='The topology filename (e.g., psf).')
        parser.add_argument('--coor_file', metavar='c', type=str, nargs="?",
                            help='The coordinate (trajectory)' +\
                                 ' filename (e.g., dcd).')
        parser.add_argument('--selection', metavar='s', type=str, nargs="?",
                            default="name CA C N O",
                            help='Which atoms to save to the JSON output ' +\
                                 'file (default: "name CA C N O"). See ' +\
                                 'https://goo.gl/kVeQuN to learn ' +\
                                 'how to construct an atom selection string.')
        # parser.add_argument('--align_sel', metavar='as', type=str, nargs="?",
        #                     default="NOT_SPECIFIED",
        #                     help='Which atoms to use for aligning the trajectorty (default, ' +\
        #                          'whatever --selection is). In some scenarios, ' +\
        #                          'you may wish to align by one set of atoms, but ' +\
        #                          'save a different set of atoms. For example, ' +\
        #                          'aligning by all CA, but saving all protein atoms.')
        parser.add_argument('--output_dir', metavar='od', type=str, nargs="?",
                            default=None,
                            help='The directory where output files should be saved. ' +\
                                 'The directory will be created if it does not exist. ' +\
                                 '(default: the directory where the coordinate file is located).')
        # parser.add_argument('--dir', action="store_true",
        #                     help='If present, no need to provide topology ' +\
        #                          'and coordinate files, scripts will select ' +\
        #                          'first trajectory in current directory.' +\
        #                          'Please only have the files you wish to load' +\
        #                          'in your working directory.')
        parser.add_argument('--stride', metavar='ns', type=int,
                            nargs="?", default=1,
                            help='How many frames to stride. For example, stride = 2 ' +\
                                 'means every other frame will be saved, and stride = 3 ' +\
                                 'means every third frame will be saved. (default: 1, no ' +\
                                 'stride).')
        # parser.add_argument('--starting_frame', metavar='sf',
        #                     type=int, nargs="?", default=0,
        #                     help='Which frame to start trimming from ' +\
        #                          '(default is 0).')
        parser.add_argument('--cum_var', metavar='var', type=float,
                            nargs='?', default=0.90,
                            help='The target cumulative variance, as a float. PCAViz will ' +\
                                 'use the minimum number of principal components required to ' +\
                                 'capture at least this variance (default: 0.90).')
        parser.add_argument('--precision', metavar='p', type=int,
                            nargs='?', default=2,
                            help='The number of decimal places to retain when rounding PCA ' +\
                                 'vectors, coefficients, and atomic coordinates ' +\
                                 '(default: 2, meaning round to the hundredths).')
        parser.add_argument('--check_accuracy', action='store_true',
                            help='Create a csv file containing the frame-to-frame RMSDs between ' +\
                                 'original- and decompressed-trajectory frames. Useful for testing ' +\
                                 'the impact of different settings on atom-position accuracy.')

        # Parse the arguments.
        args = parser.parse_args()

        # align_sel is selection if not specified.
        # if args.align_sel == "NOT_SPECIFIED":
        #     align_sel = args.selection
        # else:
        #     align_sel = args.align_sel

        # Ensure output_dir is formatted correctly.
        output_dir = args.output_dir
        if output_dir is not None and not output_dir.endswith(os.sep):
            output_dir = output_dir + os.sep
        if output_dir is not None and not os.path.exists(output_dir):
            os.mkdir(output_dir)

        # I want to convert it to a dictionary. Easier to use.
        params = vars(args)
        # params['align_sel'] = align_sel
        params['output_dir'] = output_dir

    # If we are told to find the files,
    # get first topology/coordinate/both file.
    # if params['dir']:
    #     top_file, coor_file = find_top_and_coor_files()
    #     params['top_file'] = top_file
    #     params['coor_file'] = coor_file

    print(params)

    # Return the dictionary of parameters.
    return params
