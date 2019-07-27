// Run sims in the browsers
function makePCAViz(uniqId, viewerType, datafile, loop, autoplay, visStyle, durationInMilliseconds, updateFreqInMilliseconds, windowAverageSize, caching) {
    // debugger;
    let element = jQuery('#pcaviz-viscontainer-' + uniqId);
    let config = {
        backgroundColor: 'white'
    };
    let viewer = $3Dmol.createViewer(element, config);

    // visStyle was protected, in case it contains quotes. Deprotect those values.
    visStyle = visStyle.replace(/!QUOTE!/g, '"');
    visStyle = JSON.parse(visStyle);

    // debugger;

    element.show();
    let pcaViz = new PCAViz({
        viewer: viewer,
        viewerType: viewerType,
        visStyle: visStyle, // {cartoon:{}, stick:{radius:.5,colorscheme:'Jmol'}},
        durationInMilliseconds: durationInMilliseconds,
        updateFreqInMilliseconds: updateFreqInMilliseconds, // 60 fps
        loop: loop,
        playerControlsID: document.getElementById('pcaviz-controls-' + uniqId) ? 'pcaviz-controls-' + uniqId : element,
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

function pcaVizSetCaption(uniqId) {
    let select = jQuery('#pca-file-input-' + uniqId + ' option:selected');
    let captionTxt = jQuery(select).data('caption');
    let captionDiv = jQuery('#pcaviz-figcaption-' + uniqId);
    captionDiv.html(captionTxt);
    jQuery('#pca-file-input-' + uniqId).hide();
}

function pcaVizToggleVisibility(uniqId) {
    jQuery('#pcaviz-collapsable-' + uniqId).show();
}
