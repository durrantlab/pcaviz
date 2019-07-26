# PCAViz #

## Introduction ##

PCAViz is an open-source Python/JavaScript toolkit for sharing and visualizing
MD trajectories via a web browser.

## PCAViz Compressor (Python) ##

The PCAViz Compressor is a Python script that uses lossy compression (i.e.,
principal component analysis and other strategies) to process trajectories.
The compressor saves this data as a minified JSON file that can itself be
further compressed using algorithms such as GZIP.

### Location ###

This Git repository: `pcaviz-compressor-python/`

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

## PCAViz Interpreter (JavaScript) ##

The PCAViz Interpreter is a JavaScript library that runs in users' browsers.
It downloads a PCAViz JSON file and decompresses the trajectory so the atomic
coordinates can be fed to popular WebGL-based molecular-visualization
libraries (e.g., 3DMol.js, NGL Viewer, PV).

### Location ###

This Git repository: `pcaviz-interpreter-javascript/`

### Examples of Use ###

The below example shows how to incorporate PCAViz into an HTML page, together
with the molecular-visualization library
[3DMol.js](https://3dmol.csb.pitt.edu). The
`pcaviz-interpreter-javascript/examples/` directory includes additional HTML
examples showing how to use PCAViz with other molecular-visualization
libraries.

```html
<!DOCTYPE html>
<html>
    <head>
        <title>3DMol.js</title>
    </head>
    <body>
        <style>
            /* The below position the viewer div and pcaviz-controls div relative to
            each other. */
            #pcaviz-vis-and-controls {
                width: 60%;
                height: 400px;
                position: relative;
                display: flex;
                flex-direction: column;
                flex-grow: 1;
            }

            #pcaviz-viscontainer {
                width: 100%;
                height: 100%;
            }

            #pcaviz-controls {
                min-height: 40px;
                height: 40px;
            }

            /* You can also stylize the control buttons. */
            .pcaviz-button { padding: 7px; }
            .pcaviz-button img { height: 25px; }
            .pcaviz-slider-container { padding-left: 10px; }
        </style>

        <h1>3DMol.js</h1>

        <!-- Create containers for both 3DMol.js and the trajectory pcaviz-controls
        (e.g., play button) -->
        <div id="pcaviz-vis-and-controls">
            <div id="pcaviz-viscontainer"></div>
            <div id="pcaviz-controls"></div>
        </div>

        <!-- Load the javascript files -->
        <script src="https://code.jquery.com/jquery-3.3.1.min.js"></script>
        <script src="http://3Dmol.csb.pitt.edu/build/3Dmol-min.js"></script>
        <script src="PCAViz.min.js"></script>

        <script>
            // Setup 3DMol.js
            let element = jQuery('#pcaviz-viscontainer');
            let config = {backgroundColor: 'white'};
            let viewer = $3Dmol.createViewer(element, config);

            // Tell PCAViz that you want it to interface with the 3DMol.js
            // molecular-visualization library.
            viewerType = "3DMOLJS";

            // Create a PCAViz object with your desired parameters.
            let pcaViz = new PCAViz({
                // The 3DMol.js viewer object.
                viewer: viewer,

                // The type of viewer object (3DMOLJS, NGL, PV, or GENERIC).
                // 3DMOLJS in this case.
                viewerType: viewerType,

                // The 3DMol.js style scheme (i.e., how to represent the
                // model). Note that these styles are applied to all atoms of
                // the model. To apply different styles to different atom
                // selections, run PCAViz in GENERIC mode and specify a custom
                // render() function. See 3dmdoljs_generic.html and
                // https://3dmol.csb.pitt.edu/doc/$3Dmol.GLViewer.html#setStyle
                // for more information.
                visStyle: {cartoon:{}, stick:{radius:.5,colorscheme:'Jmol'}},

                // The number of milliseconds it should take to play the
                // simulation from start to finish.
                durationInMilliseconds: 10000,

                // How often to upate the atomic coordinates. Note that
                // updating every 16.67 ms gives 60 frames per second.
                updateFreqInMilliseconds: 16.67,

                // Whether to restart the simulation once it finishes (looping
                // playback).
                loop: true,

                // The HTML ID of the DOM element where PCAViz should add
                // playback pcaviz-controls (e.g., the play button).
                playerControlsID: "pcaviz-controls",

                // The number of frames PCAViz should consider when smoothing
                // the trajectory motions via a windowed average.
                windowAverageSize: 1,

                // The type of atom-coordinate caching to use. Acceptable
                // values are "none" (no caching), "continuous" (cache each
                // frame's coordinates after they are first calculated), and
                // "pre" (calculate and cache all frame coordinates before
                // starting the animation).
                caching: "none"
            });

            // Load in the PCAViz JSON file, created with the PCAViz
            // Compressor.
            pcaViz.io.loadJSON("sim.json", () => {
                // Now that the model/simulation is loaded, you might want to
                // center it using 3DMol.js functions.
                viewer.zoomTo();

                // Now start the animation. Note that any paramters passed to
                // the start function overwrite those used to create the
                // pcaViz object.
                pcaViz.player.start({
                    durationInMilliseconds: 10000,
                    updateFreqInMilliseconds: 16.67,  // 60 fps
                    loop: true,
                    windowAverageSize: 25
                });

                // Other available PCAViz functions:
                // pcaViz.player.stop();
                // pcaViz.player.toFrame(10);
            });
        </script>
    </body>
</html>
```

## PCAViz WordPress Plugin ##

To encourage use, an easy-to-install PCAViz-powered WordPress plugin enables
"plug-and-play" trajectory visualization. The PCAViz download also provides
examples showing how to use PCAViz in any webpage.

### Location ###

This Git repository: `pcaviz-wordpress/`

### Installation ###

You can install the PCAViz WordPress plugin via the official WordPress Plugin
Directory. Search for __XXXXX__. The
[wordpress.org](https://wordpress.org/support/article/managing-plugins/#automatic-plugin-installation)
site provides useful instructions detailing how to install plugins from their
directory.

You can also install the plugin using this repository's `pcaviz-wordpress/`
directory. Either copy `pcaviz-wordpress/` directly to your server's
`wp-content/plugins/` directory
([instructions](https://wordpress.org/support/article/managing-plugins/#manual-plugin-installation)),
or create a ZIP file of the `pcaviz-wordpress/` directory and upload it via
your WordPress admin screen
([instructions](https://wordpress.org/support/article/managing-plugins/#manual-upload-via-wordpress-admin)).

### Examples of Use ###

The PCAViz WordPress Plugin uses shortcodes that can be inserted into any
WordPress post or page.

#### Major Attributes ####

1. If you use the plugin without any attributes, it will allow your site
   visitors to select from any of the compressed JSON files you've uploaded to
   your WordPress media library. Note that these files must end in
   `.compressed.json` to be recognized. <br>
   `[pcaviz]`
2. If you want to show your site visitors a specific simulation present in
   your media library, simplify specify the name of the simulation. PCAViz
   will search through your media library and will choose the first JSON file
   ending in `.compressed.json` that includes your specified text in its
   title, URL, or file name. <br>
   `[pcaviz file="larp1"]`
3. You can also control the molecular styles. The WordPress plugin uses
   [3DMol.js](https://3dmol.csb.pitt.edu) to render molecules. Passing the
   plugin a [3DMol.js
   AtomStyleSpec](https://3dmol.csb.pitt.edu/doc/types.html#AtomStyleSpec)
   JSON string changes the way the molecules are displayed. <br>
   `[pcaviz visStyle='{"cartoon": {"style": "trace", "color": "grey", "opacity": 0.75}}']`
   Note that the [3DMol.js
   documentation](https://3dmol.csb.pitt.edu/doc/$3Dmol.GLViewer.html#setStyle)
   shows AtomStyleSpec examples as JavaScript objects. PCAViz accepts only
   JSON strings. The two look similar, but with important differences (e.g.,
   keys in JSON strings must always be quoted). We recommend using a [JSON
   validator](https://jsonformatter.curiousconcept.com) to check your
   AtomStyleSpec strings.

#### Minor Attributes ####

1. To explicitly specify the height and width of the viewer, in pixels:
   <br>
   `[pcaviz height=150 width=150]`
2. To align the viewer to the left, center, or right, and to optionally add a
   caption: <br>
   `[pcaviz align="right" caption="My molecule in motion!"]`
3. To hide the playback buttons: <br>
   `[pcaviz playback_buttons="false"]`
4. To prevent the plugin from looping the animation: <br>
   `[pcaviz loop="false"]`
5. The plugin starts playing the simulation automatically by default. To
   deactivate autoplay: <br>
   `[pcaviz autoplay="false"]`
6. To specify the duration of the animation and how frequently the atomic
   positions are updated: <br>
   `[pcaviz durationInMilliseconds=10000 updateFreqInMilliseconds=50]`
7. To smooth the animation by averaging the atomic coordinates over multiple
   frames: <br>
   `[pcaviz windowAverageSize=25]`
8. To control how the plugin caches calculated coordinates: <br>
   `[pcaviz caching="continuous"]` <br>
   Acceptable values are "none" (no caching), "continuous" (cache each frame's
   coordinates after they are first calculated), and "pre" (calculate and
   cache all frame coordinates before starting the animation).
