let element = jQuery('#pcaviz-viscontainer');
let config = {backgroundColor: 'white'};
let viewer = $3Dmol.createViewer(element, config);

// Run sims in the browsers
function makePCAViz(viewer, viewerType,  datafile) {
    element.show()
    let pcaViz = new PCAViz({
        viewer: viewer,
        viewerType: viewerType,
        visStyle: {cartoon:{}, stick:{radius:.5,colorscheme:'Jmol'}},
        durationInMilliseconds: 10000,
        updateFreqInMilliseconds: 16.67,  // 60 fps
        loop: true,
        playerControlsID: document.getElementById('pcaviz-controls') ? 'pcaviz-controls' : element,
        windowAverageSize: 1,
    });

    pcaViz.io.loadJSON(datafile, () => {
        // Zoom in on the model.
        viewer.zoomTo();

        // Start playing.
        pcaViz.player.start({
            durationInMilliseconds: 10000,
            updateFreqInMilliseconds: 16.67,  // 60 fps
            loop: true,
            windowAverageSize: 25
        });
    });
}
