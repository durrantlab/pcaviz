import sys
import argparse
import random
import os
import json

from _Utils import find_top_and_coor_files

# Functions for anything to do with manipulating
# user parameters and input options.

def log_params(params):
    # Print out the parameters used as a list to write to json file.
    param_list = ["Parameters Used:"]
    for name in params:
        param_list.append(name)
        param_list.append(str(params[name]))
    return param_list
       
       
def get_params(function_type=None):

    print function_type
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
        parser.add_argument('--coor_file', metavar='c', type=str, nargs="*",
                            help=('The coordinate (trajectory)',
                                  ' filename (e.g., dcd)'))
        parser.add_argument('--selection', metavar='s', type=str, nargs="?",
                            const="name CA", default="name CA",
                            help=('The atom selection for the region of ',
                                  'interest (default "name CA"). ',
                                  'See https://goo.gl/kVeQuN for details on ',
                                  'how to construct a selection.'))
        parser.add_argument('--align_sel', metavar='as', type=str, nargs="?",
                            const="NOT_SPECIFIED", default="NOT_SPECIFIED",
                            help=('The atom selection for aligning (default, ',
                                  'whatever --selection is). You may wish to ',
                                  'align by one selection, but analyze a ',
                                  'different selection. For example, aligning ',
                                  'by all CA, but clustering on the CA of a ',
                                  'specific flexible loop.'))
        parser.add_argument('--output_dir', metavar='od', type=str, nargs="?",
                            const=None, default=None,
                            help=("The directory to store output files ",
                                  "(coordinate file's directory by default)."))
        parser.add_argument('--dont_randomize_output', action="store_true",
                            help=('Enable file-output overwrite. If present, ',
                                  'do not include a random string ',
                                  'in the output files.'))
        parser.add_argument('--dir', action="store_true",
                            help=('If present, no need to provide topology ',
                                  'and coordinate files, scripts will select ',
                                  'first trajectory in current directory.'))
        parser.add_argument('--ps_timestep', metavar='p', type=float,
                            nargs="?", const=0.002, default=0.002,
                            help=('The number of picoseconds per frame ',
                                  '(default 0.002 = 2 fs).'))

        # Then get parameters specific to the requested function.
        # Get all parameters needed for pruning and trimming.
        if function_type == 'Prune_Trim_Traj':
            parser.add_argument('--max_frames', metavar='m', type=int,
                                nargs="?", const=-1, default=-1,
                                help=('Max allowable number of frames ',
                                    '(default -1, do nothing). Stride the ',
                                    'trajectory as needed.'))
            parser.add_argument('--terminal_percent_to_keep', metavar='tpk',
                                type=float, nargs="?",
                                const=100., default=100.,
                                help=('Percentage of trajectory to keep, ',
                                    'trimming off frames from the beginning ',
                                    '(default 100).'))
            parser.add_argument('--save_pdbs', action="store_true",
                                help=('If present, pdbs of frame following ',
                                    ' pruning and trimming will be saved.'))
      
        # Get all parameters needed for performing PCA analysis.
        elif function_type == 'PCA':
            parser.add_argument('--cum_var', metavar='var', type=float,
                                nargs="?", const=0.95, 
                                help=('Cumulative variance type as decimal',
                                '(default 0.95)'))
            parser.add_argument('--axis_x_min', metavar='axmin', type=float,
                                nargs="?", const=None, default=None,
                                help=('Minimum value of the first (x) axis ',
                                    '(default None).'))
            parser.add_argument('--axis_x_max', metavar='axmax', type=float,
                                nargs="?", const=None, default=None,
                                help=('Maximum value of the first (x) axis ',
                                    '(default None).'))
            parser.add_argument('--axis_y_min', metavar='aymin', type=float,
                                nargs="?", const=None, default=None,
                                help=('Minimum value of the second (y) axis ',
                                    '(default None).'))
            parser.add_argument('--axis_y_max', metavar='aymax', type=float,
                                nargs="?", const=None, default=None,
                                help=('Maximum value of the second (y) axis ',
                                    '(default None).'))
            parser.add_argument('--clust_pca_method', metavar='cp', type=str,
                                nargs="?", const=None, default=None,
                                choices=["MeanShift", "kmeans"],
                                help=('The clustering method to use on PCA ',
                                    'projections. Can be "MeanShift" or ',
                                    '"kmeans", kmeans also requires ',
                                    '--num_clust be provided. (default does',
                                    ' not perform clustering).'))
            parser.add_argument('--num_clust', metavar='nc', type=int,
                                nargs="*", default=None,
                                help=('Number of clusters for kmeans to ',
                                    'identify. Requires one value for each',
                                    'coordinate file provided.'))
            parser.add_argument('--disp_frame', metavar='df', type=int,
                                nargs="*", default=-1,
                                help=('Displays all provided frames on all ',
                                    'PCA plots being created with labels ',
                                    '"F" + FRAME'))
            parser.add_argument('--on_server', action="store_true",
                                help=('If present, make the scripts work over',
                                    ' web server. Accepts to parameter.'))
        
        # Get all parameters needed for performing cluster analysis.
        elif function_type == 'cluster':
            parser.add_argument('--clust_method', metavar='cm', type=str,
                                nargs="?", const="AffinityPropagationNative",
                                default="AffinityPropagationNative",
                                choices=["AffinityPropagationNative","DBSCAN"],
                                help=('The clustering method to use. ',
                                    'Can be "AffinityPropagationNative" ',
                                    'or "DBSCAN" (default ',
                                    '"AffinityPropagationNative").'))
            parser.add_argument('--apn_pref', metavar='ap', type=float,
                                nargs="?", const=-1., default=-1.,
                                help=('Preference parameter used in the ',
                                    'Affinity Propagation algorithm for ',
                                    'clustering (default -1.0). A high ',
                                    'preference value results in many ',
                                    'clusters, a low preference will result',
                                    ' in fewer numbers of clusters.'))
        

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
    # get first topology/coordinate/both file and save the coor_file as a list.
    if params['dir']:
        top_file, coor_file = find_top_and_coor_files()
        params['top_file'] = top_file
        if coor_file is not None:
            params['coor_file'] = [coor_file]
    
    # Return the dictionary of parameters.
    return params    

