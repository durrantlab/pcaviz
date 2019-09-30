# PCAViz #

## Introduction ##

PCAViz is an open-source Python/JavaScript toolkit for sharing and visualizing
MD trajectories via a web browser. We release it under the terms of the
open-source GNU General Public License, version 2 (see `LICENSE.md`). This
package includes only the Compressor (Python) component of the toolkit. See
[http://git.durrantlab.com/jdurrant/pcaviz](http://git.durrantlab.com/jdurrant/pcaviz)
for the complete toolkit, including the JavaScript and WordPress components.

## PCAViz Compressor (Python) ##

The PCAViz Compressor is a Python script that uses lossy compression (i.e.,
principal component analysis and other strategies) to process trajectories.
The compressor saves this data as a minified JSON file that can itself be
further compressed using algorithms such as GZIP.

### Available User Parameters ###

```text
  -h, --help         show this help message and exit
  --top_file [t]     The topology filename (e.g., psf).
  --coor_file [c]    The coordinate (trajectory) filename (e.g., dcd).
  --selection [s]    Which atoms to save to the JSON output file (default:
                     "name CA C N O"). See https://goo.gl/kVeQuN to learn how
                     to construct an atom selection string.
  --output_dir [od]  The directory where output files should be saved. The
                     directory will be created if it does not exist. (default:
                     the directory where the coordinate file is located).
  --stride [ns]      How many frames to stride. For example, stride = 2 means
                     every other frame will be saved, and stride = 3 means
                     every third frame will be saved. (default: 1, no stride).
  --cum_var [var]    The target cumulative variance, as a float. PCAViz will
                     use the minimum number of principal components required
                     to capture at least this variance (default: 0.90).
  --precision [p]    The number of decimal places to retain when rounding PCA
                     vectors, coefficients, and atomic coordinates (default:
                     2, meaning round to the hundredths).
  --check_accuracy   Create a csv file containing the frame-to-frame RMSDs
                     between original- and decompressed-trajectory frames.
                     Useful for testing the impact of different settings on
                     atom-position accuracy.
  --test             Tests PCAViz to make sure all components are functioning.
```

### Examples of Use ###

1. Create a compressed JSON file from a topology (1J8K_example.psf) and a
   trajectory/coordinate (1J8K_example.dcd) file. <br>
   `python PCAViz.py --top_file examples/1J8K_example.psf --coor_file examples/1J8K_example.dcd`

2. PDB files can also contain multiple frames. In this case, the same file
   serves as the topology and trajectory file. <br>
   `python PCAViz.py --top_file examples/1J8K_example.pdb --coor_file examples/1J8K_example.pdb`

3. By default, PCAViz includes only the backbone atoms in the output. These
   atoms should be enough for cartoon-style visualization. But you can select
   your own atoms to include in the output. See https://goo.gl/kVeQuN to learn
   how to construct an atom-selection string. <br>
   `python PCAViz.py --top_file examples/1J8K_example.pdb --coor_file examples/1J8K_example.pdb --selection "name *"`

4. Striding the trajectory frames can reduce file sizes. PCAViz will
   interpolate between the remaining frames to fill in the frames that are
   missing. Here we keep only every other frame: <br>
   `python PCAViz.py --top_file examples/1J8K_example.pdb --coor_file examples/1J8K_example.pdb --stride 2`

5. PCAViz allows users to control the compression settings. Two settings are
   available. First, the user can specify how much of the cumulative variance
   should be explained by the principal components. PCAViz will pick the
   number of components required to meet that goal. Second, the user can
   specify how precisely PCAViz should represent numeric values (e.g.,
   principal-component coefficients). For example, here we tell PCAViz to
   produce a trajectory that accounts for 80% of the variance, and to round all
   numbers in the output JSON file to the nearest hundredth (two decimal
   places): <br>
   `python PCAViz.py --top_file examples/1J8K_example.pdb --coor_file examples/1J8K_example.pdb --cum_var 0.8 --precision 2`

6. To find the ideal --cum_var and --precision parameters, you may wish to
   check how closely the PCAViz-compressed trajectory matches the original
   trajectory. You can instruct PCAViz to output a CSV file that provides a
   frame-by-frame RMSD comparison between the two trajectories. This option
   also outputs an XYZ trajectory file that you can visually compare to the
   original. <br>
   `python PCAViz.py --top_file examples/1J8K_example.pdb --coor_file examples/1J8K_example.pdb --check_accuracy`

7. By default, PCAViz saves the compressed JSON file to the same directory
   where the coordinate file is located. You can specify a different output
   directory if needed. The directory will be created if it doesn't exist.
   <br>
   `python PCAViz.py --top_file examples/1J8K_example.pdb --coor_file examples/1J8K_example.pdb --output_dir "my_dir"`

8. For debugging purposes, PCAViz also includes an option to test whether the
   code is fully functional. <br>
   `python PCAViz.py --test`
