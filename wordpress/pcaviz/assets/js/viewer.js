let element = jQuery('#viscontainer');
let config = {backgroundColor: 'white'};
let viewer = $3Dmol.createViewer(element, config);

// Run sims in the browsers
function makeBrowserSim(viewer, viewerType,  datafile) {
    element.show()
    let browserSim = new BrowserSim({
        viewer: viewer,
        viewerType: viewerType,
        visStyle: {cartoon:{}, stick:{radius:.5,colorscheme:'Jmol'}},
        durationInMilliseconds: 10000,
        updateFreqInMilliseconds: 16.67,  // 60 fps
        loop: true,
        playerControlsID: document.getElementById('controls') ? 'controls' : element,
        windowAverageSize: 1,
    });

    browserSim.io.loadJSON(datafile, () => {
        // Zoom in on the model.
        viewer.zoomTo();

        // Start playing.
        browserSim.player.start({
            durationInMilliseconds: 10000,
            updateFreqInMilliseconds: 16.67,  // 60 fps
            loop: true,
            windowAverageSize: 25
        });
    });
}
