// Setup 3DMol.js
let element = jQuery('#pcaviz-viscontainer');
let config = {backgroundColor: 'white'};
let viewer = $3Dmol.createViewer(element, config);

// Run sims in the browsers
function makePCAViz(viewer, viewerType) {
    let pcaViz = new PCAViz({
        viewer: viewer,
        viewerType: viewerType,
        visStyle: {cartoon:{}, stick:{radius:.5,colorscheme:'Jmol'}},
        durationInMilliseconds: 10000,
        updateFreqInMilliseconds: 16.67,  // 60 fps
        loop: true,
        windowAverageSize: 1,
    });

    pcaViz.io.loadJSON("/data.json", () => {
        pcaViz.player.start({
            durationInMilliseconds: 10000,
            updateFreqInMilliseconds: 16.67,  // 60 fps
            loop: true,
            windowAverageSize: 25
        });

        window.pcaViz = pcaViz;

        // setTimeout(() => {
        //    pcaViz.player.stop();
        //    setTimeout(() => {
        //        pcaViz.player.toFrame(159);
        //    }, 1000);
        // }, 5000);
    });
}
makePCAViz(viewer, '3DMOLJS');
