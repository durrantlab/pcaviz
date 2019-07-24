let element = jQuery('#pcaviz-viscontainer');
let config = {backgroundColor: 'white'};
let viewer = $3Dmol.createViewer(element, config);

// Run sims in the browsers
function makePCAViz(viewer, viewerType, datafile, loop, autoplay, visStyle, durationInMilliseconds, updateFreqInMilliseconds, windowAverageSize, caching) {
    console.log(visStyle);

    // visStyle was protected, in case it contains quotes. Deprotect those values.
    // visStyle = visStyle.replace(/!SINGLEQUOTE!/g, "'");
    visStyle = visStyle.replace(/!QUOTE!/g, '"');

    console.log(visStyle);
    visStyle = JSON.parse(visStyle);
    console.log(visStyle);

    element.show()
    let pcaViz = new PCAViz({
        viewer: viewer,
        viewerType: viewerType,
        visStyle: visStyle, // {cartoon:{}, stick:{radius:.5,colorscheme:'Jmol'}},
        durationInMilliseconds: durationInMilliseconds,
        updateFreqInMilliseconds: updateFreqInMilliseconds,  // 60 fps
        loop: loop,
        playerControlsID: document.getElementById('pcaviz-controls') ? 'pcaviz-controls' : element,
        windowAverageSize: windowAverageSize,
        caching: caching
    });

    pcaViz.io.loadJSON(datafile, () => {
        // Zoom in on the model.
        viewer.zoomTo();

        if (autoplay) {
            // Start playing.
            pcaViz.player.start();
        }
    });
}
