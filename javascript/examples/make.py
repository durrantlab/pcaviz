# This creates the demo files showing how to use the browser component.

import shutil

# Copy the js files.
shutil.copy("../BrowserSim.js", "./")
shutil.copy("../BrowserSim.min.js", "./")
shutil.copy("../data.json", "./")

# A function to return a general HTML template.
def HTML(title, body):
    return """<html>
    <head>
        <title>{title}</title>
    </head>
    <body>
{body}
    </body>
</html>""".format(
        title=title,
        body=body
    )

def demo_HTML(page):
    extra_browser_sim_params = "\n                    ".join(page.extra_browser_sim_params)
    if extra_browser_sim_params != "":
        extra_browser_sim_params = "\n                    " + extra_browser_sim_params

    new_body = """        <style>
            /* The below position the viewer div and controls div relative to each other. */
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
        <h1>""" + page.title + """</h1>
        <div id="vis-and-controls">
            <div id="viscontainer"></div>
            <div id="controls"></div>
        </div>

        <!-- Load the javascript files -->
        <script src="https://code.jquery.com/jquery-3.3.1.min.js"></script>
        """ + "\n        ".join(page.external_js) + """
        <script src="BrowserSim.min.js"></script>
        <script>
            // Setup """ + page.title + """
            """ + "\n            ".join(page.viewer_setup) + """

            // Run sims in the browsers
            function makeBrowserSim(viewer, viewerType) {
                let browserSim = new BrowserSim({
                    viewer: viewer,
                    viewerType: viewerType,
                    visStyle: """ + page.vis_style + """
                    durationInMilliseconds: 10000,
                    updateFreqInMilliseconds: 16.67,  // 60 fps
                    loop: true,
                    playerControlsID: "controls",
                    windowAverageSize: 1,""" + extra_browser_sim_params + """
                });

                browserSim.io.loadJSON("data.json", () => {
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
            makeBrowserSim(viewer, '""" + page.id + """');
        </script>
"""
    return HTML(page.title, new_body)

# Now get info about each page
class Page:
    def __init__(self, title, filename, viewer_setup, external_js, id, vis_style, extra_browser_sim_params=[]):
        self.title = title
        self.filename = filename
        self.viewer_setup = viewer_setup
        self.external_js = external_js
        self.id = id
        self.vis_style = vis_style
        self.extra_browser_sim_params = extra_browser_sim_params

pages = [
    Page("3DMol.js", "3dmoljs.html", [
            "let element = jQuery('#viscontainer');",
            "let config = {backgroundColor: 'white'};",
            "let viewer = $3Dmol.createViewer(element, config);"
        ], [
            '<script src="http://3Dmol.csb.pitt.edu/build/3Dmol-min.js"></script>'
        ],
        "3DMOLJS",
        "{cartoon:{}, stick:{radius:.5,colorscheme:'Jmol'}},"
        # "{sphere: {colorscheme: 'Jmol'}},"
    ),
    Page("3DMol.js Generic Mode", "3dmoljs_generic.html", [
            "let element = jQuery('#viscontainer');",
            "let config = {backgroundColor: 'white'};",
            "let viewer = $3Dmol.createViewer(element, config);"
        ], [
            '<script src="http://3Dmol.csb.pitt.edu/build/3Dmol-min.js"></script>'
        ],
        "GENERIC",
        "undefined,  // Not needed when in generic mode",
        [
            'loadPDBTxt: (pdbTxt, viewer, browserSim) => viewer.addModel(pdbTxt, "pdb"),',
            'updateAtomPositions: (newAtomCoors, model, viewer, browserSim) => {',
            '    var atoms = model.selectedAtoms({});',
            '    for (let atomIdx=0; atomIdx<atoms.length; atomIdx++) {',
            '        let coors = newAtomCoors[atomIdx];',
            '        atoms[atomIdx]["x"] = coors[0];',
            '        atoms[atomIdx]["y"] = coors[1];',
            '        atoms[atomIdx]["z"] = coors[2];',
            '    }',
            '},',
            'render: (model, viewer, browserSim) => {',
            '    model.setStyle(',
            '        {}, { sphere: { color: "green" }}',
            '    );',
            '    viewer.render();',
            '},'
        ]
    ),
    Page("NGLViewer", "ngl.html", [
            "var viewer = new NGL.Stage('viscontainer');",
            "viewer.viewer.setBackground( 0xFFFFFF );"
        ], [
            '<script src="https://unpkg.com/ngl@0.10.4/dist/ngl.js"></script>'
        ],
        "NGL",
        "undefined,  // No visStyle needed for NGLViewer"
    ),
    Page("NGLViewer Generic Mode", "ngl_generic.html", [
            "var viewer = new NGL.Stage('viscontainer');",
            "viewer.viewer.setBackground( 0xFFFFFF );"
        ], [
            '<script src="https://unpkg.com/ngl@0.10.4/dist/ngl.js"></script>'
        ],
        "GENERIC",
        "undefined,  // Not needed when in generic mode",
        [
            '// Normally loadPDBTxt() would return the model, but NGL Viewer',
            '// uses promises. So set the browserSim._model variable directly.',
            'loadPDBTxt: (pdbTxt, viewer, browserSim) => {',
            '    viewer.loadFile(',
            '        new Blob(',
            '            [pdbTxt],',
            '            {type: "text/plain"}',
            '        ),',
            '        {ext: "pdb", defaultRepresentation: true}',
            '    ).then((e) => {',
            '        browserSim._model = e;  // Cannot return, because it\'s a promise',
            '    });',
            '},',
            'updateAtomPositions: (newAtomCoors, model, viewer, browserSim) => {',
            '    var atomIdx = 0; window.m = model;',
            '    model.structure.eachAtom((atom, idx) => {',
            '        var coors = newAtomCoors[atomIdx];',
            '        atom["x"] = coors[0];',
            '        atom["y"] = coors[1];',
            '        atom["z"] = coors[2];',
            '        atomIdx++;',
            '    });',
            '},',
            'render: (model, viewer, browserSim) => model.rebuildRepresentations(),',
        ]
    ),
    Page("PV", "pv.html", [
            "var pvViewer = pv.Viewer(",
            "    document.getElementById('viscontainer'),",
            "    {quality:'medium', width:'auto', height: 'auto', antialias: true, outline: true}",
            ");",
            "// Not entirely encapsulated, so put both the library and the viewer in a single object.",
            "var viewer = {viewer: pvViewer, library: pv};"
        ], [
            '<script src="http://biasmv.github.io/pv/js/pv.min.js"></script>'
        ],
        "PV",
        "(structure) => { pvViewer.spheres('all', structure); },  // Callback to set visStyle."
    ),
    Page("PV Generic Mode", "pv_generic.html", [
            "var pvViewer = pv.Viewer(",
            "    document.getElementById('viscontainer'),",
            "    {quality:'medium', width:'auto', height: 'auto', antialias: true, outline: true}",
            ");",
            "// Not entirely encapsulated, so put both the library and the viewer in a single object.",
            "var viewer = {viewer: pvViewer, library: pv};"
        ], [
            '<script src="http://biasmv.github.io/pv/js/pv.min.js"></script>'
        ],
        "GENERIC",
        "undefined,  // Not needed when in generic mode",
        [
            'loadPDBTxt: (pdbTxt, viewer, browserSim) => viewer.library.io.pdb(pdbTxt),',
            'updateAtomPositions: (newAtomCoors, model, viewer, browserSim) => {',
            '    let atomIdx = 0;',
            '    model.eachAtom((atom, idx) => {',
            '        let coors = newAtomCoors[atomIdx];',
            '        atom._bV = coors;',
            '        atomIdx++;',
            '    });',
            '},',
            'render: (model, viewer, browserSim) => {',
            '    viewer.viewer.clear();  // Because encapsulated',
            '    viewer.viewer.spheres("all", model);',
            '},',
        ]
    ),
]

# Create the index.html file.
with open("index.html", "w") as f:
    linksTxt = "\n".join(
        ['        <p><a href="' + page.filename + '">' + page.title + '</a></p>' for page in pages]
    )
    f.write(HTML("BrowserSim Examples", linksTxt))

# Create the individual files.
for page in pages:
    with open(page.filename, "w") as f:
        f.write(
            demo_HTML(page)
        )
