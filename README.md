PCAViz
======

Introduction
------------

PCAViz is an open-source Python/JavaScript toolkit for sharing and visualizing
MD trajectories via a web browser.

PCAViz Compressor (Python)
--------------------------

The PCAViz Compressor is a Python script that uses lossy compression (i.e.,
principal component analysis and other strategies) to process trajectories.
The compressor saves this data as a minified JSON file that can itself be
further compressed using algorithms such as GZIP.



PCAViz Interpreter (JavaScript)
-------------------------------

The PCAViz Interpreter is a JavaScript library that runs in users' browsers.
It downloads a PCAViz JSON file and decompresses the trajectory so the atomic
coordinates can be fed to popular WebGL-based molecular-visualization
libraries (e.g., 3DMol.js, NGL Viewer, PV).

The below example shows how to incorporate PCAViz into an HTML page, together
with the molecular-visualization library
[3DMol.js](https://3dmol.csb.pitt.edu). The `javascript/examples/` directory
includes additional HTML examples showing how to use PCAViz with other
molecular-visualization libraries.

```html
<!DOCTYPE html>
<html>
    <head>
        <title>3DMol.js</title>
    </head>
    <body>
        <style>
            /* The below position the viewer div and controls div relative to
            each other. */
            #vis-and-controls {
                width: 60%;
                height: 400px;
                position: relative;
                display: flex;
                flex-direction: column;
                flex-grow: 1;
            }

            #viscontainer {
                width: 100%;
                height: 100%;
            }

            #controls {
                min-height: 40px;
                height: 40px;
            }

            /* You can also stylize the control buttons. */
            .pcaviz-button { padding: 7px; }
            .pcaviz-button img { height: 25px; }
            .pcaviz-slider-container { padding-left: 10px; }
        </style>

        <h1>3DMol.js</h1>

        <!-- Create containers for both 3DMol.js and the trajectory controls
        (e.g., play button) -->
        <div id="vis-and-controls">
            <div id="viscontainer"></div>
            <div id="controls"></div>
        </div>

        <!-- Load the javascript files -->
        <script src="https://code.jquery.com/jquery-3.3.1.min.js"></script>
        <script src="http://3Dmol.csb.pitt.edu/build/3Dmol-min.js"></script>
        <script src="PCAViz.min.js"></script>

        <script>
            // Setup 3DMol.js
            let element = jQuery('#viscontainer');
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
                // model).
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
                // playback controls (e.g., the play button).
                playerControlsID: "controls",

                // The number of frames PCAViz should consider when smoothing
                // the trajectory motions via a windowed average.
                windowAverageSize: 1,
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

PCAViz WordPress Plugin
-----------------------

To encourage use, an easy-to-install PCAViz-powered WordPress plugin enables
"plug-and-play" trajectory visualization. The PCAViz download also provides
examples showing how to use PCAViz in any webpage.
