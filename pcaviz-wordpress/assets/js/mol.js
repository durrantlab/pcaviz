            // Setup 3DMol.js
            let element = jQuery('#viscontainer');
            let config = {backgroundColor: 'white'};
            let viewer = $3Dmol.createViewer(element, config);

            // Run sims in the browsers
            function makeBrowserSim(viewer, viewerType) {
                let browserSim = new BrowserSim({
                    viewer: viewer,
                    viewerType: viewerType,
                    visStyle: {cartoon:{}, stick:{radius:.5,colorscheme:'Jmol'}},
                    durationInMilliseconds: 10000,
                    updateFreqInMilliseconds: 16.67,  // 60 fps
                    loop: true,
                    windowAverageSize: 1,
                });

                browserSim.io.loadJSON("/data.json", () => {
                    browserSim.player.start({
                        durationInMilliseconds: 10000,
                        updateFreqInMilliseconds: 16.67,  // 60 fps
                        loop: true,
                        windowAverageSize: 25
                    });

                    window.browserSim = browserSim;

                    // setTimeout(() => {
                    //    browserSim.player.stop();
                    //    setTimeout(() => {
                    //        browserSim.player.toFrame(159);
                    //    }, 1000);
                    // }, 5000);
                });
            }
            makeBrowserSim(viewer, '3DMOLJS');
