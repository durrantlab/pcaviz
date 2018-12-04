import sys
import argparse
import random
import os
import json

from _Utils import find_top_and_coor_files

#Obtains and formats all user provided parameters as a dictionary.
def get_params():

    run_id = "id" + str(random.randint(0,100000000))
    possible_json_file = sys.argv[1].lower()
    # If a json file is provided,
    # read in parameters from json file (web implementation).
    if os.path.splitext(possible_json_file)[1][1:] == "json":
        with open(sys.argv[1], 'r') as param_file:
            params = json.load(param_file)
    else:
        # Otherwise parse the parameters from the system arguments.
        parser = argparse.ArgumentParser(description='Process MD trajectory.')

        # First get all universal parameters.
        parser.add_argument('--top_file', metavar='t', type=str, nargs="?",
                            help='The topology filename (e.g., psf)')
        parser.add_argument('--coor_file', metavar='c', type=str, nargs="?",
                            help=('The coordinate (trajectory)',
                                  ' filename (e.g., dcd)'))
        parser.add_argument('--selection', metavar='s', type=str, nargs="?",
                            default="name CA",
                            help=('The atom selection for the region of ',
                                  'interest (default "name CA"). ',
                                  'See https://goo.gl/kVeQuN for details on ',
                                  'how to construct a selection.'))
        parser.add_argument('--align_sel', metavar='as', type=str, nargs="?",
                            default="NOT_SPECIFIED",
                            help=('The atom selection for aligning (default, ',
                                  'whatever --selection is). You may wish to ',
                                  'align by one selection, but analyze a ',
                                  'different selection. For example, aligning ',
                                  'by all CA, but clustering on the CA of a ',
                                  'specific flexible loop.'))
        parser.add_argument('--output_dir', metavar='od', type=str, nargs="?",
                            default=None,
                            help=('The directory to store output files ',
                                  '(coordinate file's directory by default).'))
        parser.add_argument('--dont_randomize_output', action="store_true",
                            help=('Enable file-output overwrite. If present, ',
                                  'do not include a random string ',
                                  'in the output files.'))
        parser.add_argument('--dir', action="store_true",
                            help=('If present, no need to provide topology ',
                                  'and coordinate files, scripts will select ',
                                  'first trajectory in current directory.'))
        parser.add_argument('--ps_timestep', metavar='p', type=float,
                            nargs="?", default=0.002,
                            help=('The number of picoseconds per frame ',
                                  '(default 0.002 = 2 fs).'))
        parser.add_argument('--num_stride', metavar='ns', type=int,
                            nargs="?", default=1,
                            help=('How many frames to stride by ',
                                  '(default 1, do nothing).'))
        parser.add_argument('--starting_frame', metavar='sf',
                            type=int, nargs="?", default=0,
                            help=('Which frame to start trimming from ',
                                '(default is 0).'))
        parser.add_argument('--cum_var', metavar='var', type=float,
                            nargs='?', default=0.90, 
                            help=('Cumulative variance type as decimal',
                                  '(default is 0.95)'))
        parser.add_argument('num_decimals', metavar='nd', type=int,
                            nargs='?', default=2,
                            help=('How many decimal places to round PCA ',
                                  'vectors, coefficients, and atomic ',
                                  'coordinates to. (default is 2)'))
        
        # Parse the arguments.
        args = parser.parse_args()

        # align_sel is selection if not specified.
        if args.align_sel == "NOT_SPECIFIED":
            align_sel = args.selection
        else:
            align_sel = args.align_sel

        # Ensure output_dir is formatted correctly.
        output_dir = args.output_dir
        if output_dir is not None and not output_dir.endswith(os.sep):
            output_dir = output_dir + os.sep

        # I want to convert it to a dictionary. Easier to use.
        params = vars(args)
        params['align_sel'] = align_sel
        params['output_dir'] = output_dir
    
    params['id'] = run_id
    
    # If we are told to find the files,
    # get first topology/coordinate/both file.
    if params['dir']:
        top_file, coor_file = find_top_and_coor_files()
        params['top_file'] = top_file
        params['coor_file'] = coor_file
    
    # Return the dictionary of parameters.
    return params    

