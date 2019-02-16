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
            .mol-container {
                width: 60%;
                height: 400px;
                position: relative;
            }
        </style>
        <h1>""" + page.title + """</h1>
        """ + ('' if page.title == "JSMol" else '<div id="viscontainer" class="mol-container"></div>\n') + """
        <!-- Load the javascript files -->
        <script src="https://code.jquery.com/jquery-3.3.1.min.js"></script>
        """ + "\n        ".join(page.external_js) + """
        <script src="BrowserSim.min.js"></script>
        <script>
            // A PDB file to play with.
            var testPDB = "ATOM      0  X   XXX X 852       3.768   2.214   0.253  1.00  0.00           X\\nATOM      1  X   XXX X 852      -2.041  -2.608  -0.161  1.00  0.00           X\\nATOM      2  X   XXX X 852       0.079  -1.132   0.667  1.00  0.00           X\\nATOM      3  X   XXX X 852      -3.375  -2.411   2.668  1.00  0.00           X\\nATOM      4  X   XXX X 852      -0.392  -0.246  -0.161  1.00  0.00           X\\nATOM      5  X   XXX X 852       2.590   0.886  -0.253  1.00  0.00           X\\nATOM      6  X   XXX X 852      -0.235   2.804  -1.150  1.00  0.00           X\\nATOM      7  X   XXX X 852       3.454  -0.295  -1.219  1.00  0.00           X\\nATOM      8  X   XXX X 852      -0.392   0.000  -0.069  1.00  0.00           X\\nATOM      9  X   XXX X 852       0.942   0.836   0.322  1.00  0.00           X\\nATOM     10  X   XXX X 852      -0.863  -1.033  -1.265  1.00  0.00           X\\nATOM     11  X   XXX X 852       2.041  -1.328  -0.069  1.00  0.00           X\\nATOM     12  X   XXX X 852       3.375  -1.181  -1.449  1.00  0.00           X\\nATOM     13  X   XXX X 852       0.000   0.098   0.184  1.00  0.00           X\\nATOM     14  X   XXX X 852       1.334   0.344  -0.115  1.00  0.00           X\\nATOM     15  X   XXX X 852       2.276   0.541  -0.046  1.00  0.00           X\\nATOM     16  X   XXX X 852      -0.706   0.344   0.207  1.00  0.00           X\\nATOM     17  X   XXX X 852       0.942   0.640  -0.046  1.00  0.00           X\\nATOM     18  X   XXX X 852       2.041   0.148   0.023  1.00  0.00           X\\nATOM     19  X   XXX X 852      -1.413  -0.197  -0.092  1.00  0.00           X\\nATOM     20  X   XXX X 852       1.492   0.246   0.621  1.00  0.00           X\\nATOM     21  X   XXX X 852      -0.785   0.295   0.253  1.00  0.00           X\\nATOM     22  X   XXX X 852       0.628  -0.886   0.253  1.00  0.00           X\\nATOM     23  X   XXX X 852       0.863   0.148  -0.966  1.00  0.00           X\\nATOM     24  X   XXX X 852      -0.628  -0.049  -0.023  1.00  0.00           X\\nATOM     25  X   XXX X 852      -0.942   0.984   0.000  1.00  0.00           X\\nATOM     26  X   XXX X 852      -0.785   0.738  -0.115  1.00  0.00           X\\nATOM     27  X   XXX X 852      -0.078   0.541  -0.230  1.00  0.00           X\\nATOM     28  X   XXX X 852       0.392  -0.689  -0.230  1.00  0.00           X\\nATOM     29  X   XXX X 852      -0.471   0.492   0.138  1.00  0.00           X\\nATOM     30  X   XXX X 852       1.020  -0.246  -0.345  1.00  0.00           X\\nATOM     31  X   XXX X 852      -1.570   0.836  -0.115  1.00  0.00           X\\nATOM     32  X   XXX X 852       0.000  -0.738   0.207  1.00  0.00           X\\nATOM     33  X   XXX X 852       0.785  -0.344   0.069  1.00  0.00           X\\nATOM     34  X   XXX X 852      -1.727  -0.246   0.253  1.00  0.00           X\\nATOM     35  X   XXX X 852      -0.157  -0.197   0.000  1.00  0.00           X\\nATOM     36  X   XXX X 852       0.549   0.984  -0.391  1.00  0.00           X\\nATOM     37  X   XXX X 852       0.471   0.295   0.138  1.00  0.00           X\\nATOM     38  X   XXX X 852      -0.000   0.640  -0.069  1.00  0.00           X\\nATOM     39  X   XXX X 852      -0.942  -0.098  -0.023  1.00  0.00           X\\nATOM     40  X   XXX X 852      -1.256   0.197   0.092  1.00  0.00           X\\nATOM     41  X   XXX X 852      -0.314   0.098  -0.207  1.00  0.00           X\\nATOM     42  X   XXX X 852      -0.785  -0.836  -0.184  1.00  0.00           X\\nATOM     43  X   XXX X 852      -0.236  -0.787  -0.276  1.00  0.00           X\\nATOM     44  X   XXX X 852       1.020  -0.197  -0.253  1.00  0.00           X\\nATOM     45  X   XXX X 852       0.078   0.344   0.115  1.00  0.00           X\\nATOM     46  X   XXX X 852      -1.256  -0.148  -0.253  1.00  0.00           X\\nATOM     47  X   XXX X 852      -2.198  -0.443   0.000  1.00  0.00           X\\nATOM     48  X   XXX X 852       0.078  -0.443  -0.207  1.00  0.00           X\\nATOM     49  X   XXX X 852       0.314  -0.098  -0.345  1.00  0.00           X\\nATOM     50  X   XXX X 852      -0.471   0.148  -0.046  1.00  0.00           X\\nATOM     51  X   XXX X 852       0.942  -0.197   0.253  1.00  0.00           X\\nATOM     52  X   XXX X 852       0.785   0.246  -0.230  1.00  0.00           X\\nATOM     53  X   XXX X 852      -1.727   0.197   0.161  1.00  0.00           X\\nATOM     54  X   XXX X 852       1.021   0.689   0.069  1.00  0.00           X\\nATOM     55  X   XXX X 852      -0.314  -0.344  -0.207  1.00  0.00           X\\nATOM     56  X   XXX X 852      -0.157  -0.689  -0.092  1.00  0.00           X\\nATOM     57  X   XXX X 852       0.942   0.098   0.046  1.00  0.00           X\\nATOM     58  X   XXX X 852      -0.314  -0.590  -0.161  1.00  0.00           X\\nATOM     59  X   XXX X 852      -0.706   0.689   0.299  1.00  0.00           X\\nATOM     60  X   XXX X 852       0.314  -0.148  -0.023  1.00  0.00           X\\nATOM     61  X   XXX X 852      -0.078   0.394   0.092  1.00  0.00           X\\nATOM     62  X   XXX X 852       0.314   0.492  -0.299  1.00  0.00           X\\nATOM     63  X   XXX X 852       1.492   0.098   0.253  1.00  0.00           X\\nATOM     64  X   XXX X 852      -0.157  -0.344   0.115  1.00  0.00           X\\nATOM     65  X   XXX X 852       0.078  -0.344   0.276  1.00  0.00           X\\nATOM     66  X   XXX X 852       0.628  -0.148  -0.299  1.00  0.00           X\\nATOM     67  X   XXX X 852      -0.078  -1.132   0.092  1.00  0.00           X\\nATOM     68  X   XXX X 852       0.706  -0.492  -0.023  1.00  0.00           X\\nATOM     69  X   XXX X 852      -0.471  -0.295  -0.115  1.00  0.00           X\\nATOM     70  X   XXX X 852       0.236  -0.049   0.483  1.00  0.00           X\\nATOM     71  X   XXX X 852      -0.471  -0.148   0.299  1.00  0.00           X\\nATOM     72  X   XXX X 852      -0.863  -0.394   0.230  1.00  0.00           X\\nATOM     73  X   XXX X 852       0.785  -0.246   0.046  1.00  0.00           X\\nATOM     74  X   XXX X 852      -1.020  -0.246  -0.483  1.00  0.00           X\\nATOM     75  X   XXX X 852       0.549  -0.689   0.000  1.00  0.00           X\\nATOM     76  X   XXX X 852       0.235   0.197   0.230  1.00  0.00           X\\nATOM     77  X   XXX X 852       1.334   0.394  -0.069  1.00  0.00           X\\nATOM     78  X   XXX X 852      -0.236   0.344  -0.046  1.00  0.00           X\\nATOM     79  X   XXX X 852       0.549  -0.148   0.161  1.00  0.00           X\\nATOM     80  X   XXX X 852       0.392  -0.295   0.391  1.00  0.00           X\\nATOM     81  X   XXX X 852      -1.178   0.492  -0.667  1.00  0.00           X\\nATOM     82  X   XXX X 852      -0.078  -0.246  -0.230  1.00  0.00           X\\nATOM     83  X   XXX X 852      -0.628  -0.295   0.069  1.00  0.00           X\\nATOM     84  X   XXX X 852       0.157  -0.295  -0.207  1.00  0.00           X\\nATOM     85  X   XXX X 852       0.942  -0.492  -0.046  1.00  0.00           X\\nATOM     86  X   XXX X 852      -0.471  -0.640   0.161  1.00  0.00           X\\nATOM     87  X   XXX X 852      -0.942   0.541   0.460  1.00  0.00           X\\nATOM     88  X   XXX X 852       0.157   0.148   0.414  1.00  0.00           X\\nATOM     89  X   XXX X 852       0.706  -0.590   0.207  1.00  0.00           X\\nATOM     90  X   XXX X 852       1.492  -0.541   0.276  1.00  0.00           X\\nATOM     91  X   XXX X 852       0.471  -0.492   0.092  1.00  0.00           X\\nATOM     92  X   XXX X 852       0.628   1.279  -0.276  1.00  0.00           X\\nATOM     93  X   XXX X 852      -0.078   0.394  -0.253  1.00  0.00           X\\nATOM     94  X   XXX X 852      -1.020   0.148   0.207  1.00  0.00           X\\nATOM     95  X   XXX X 852       0.000  -0.049  -0.207  1.00  0.00           X\\nATOM     96  X   XXX X 852       0.157   0.197  -0.207  1.00  0.00           X\\nATOM     97  X   XXX X 852      -0.236  -0.295   0.023  1.00  0.00           X\\nATOM     98  X   XXX X 852       1.021   0.443  -0.092  1.00  0.00           X\\nATOM     99  X   XXX X 852       1.178   0.295   0.000  1.00  0.00           X\\nATOM    100  X   XXX X 852       0.628  -0.443   0.000  1.00  0.00           X\\nATOM    101  X   XXX X 852      -0.785  -0.689   0.207  1.00  0.00           X\\nATOM    102  X   XXX X 852       0.706  -0.541  -0.138  1.00  0.00           X\\nATOM    103  X   XXX X 852       0.393   0.246  -0.023  1.00  0.00           X\\nATOM    104  X   XXX X 852       1.099   0.394  -0.138  1.00  0.00           X\\nATOM    105  X   XXX X 852      -0.078   0.295   0.115  1.00  0.00           X\\nATOM    106  X   XXX X 852      -0.000  -0.049  -0.092  1.00  0.00           X\\nATOM    107  X   XXX X 852       0.942  -0.541  -0.023  1.00  0.00           X\\nATOM    108  X   XXX X 852      -0.079  -0.197   0.046  1.00  0.00           X\\nATOM    109  X   XXX X 852      -0.235  -0.098  -0.046  1.00  0.00           X\\nATOM    110  X   XXX X 852       0.079   0.541   0.023  1.00  0.00           X\\nATOM    111  X   XXX X 852      -0.157   0.098  -0.023  1.00  0.00           X\\nATOM    112  X   XXX X 852       0.628  -0.098   0.092  1.00  0.00           X\\nATOM    113  X   XXX X 852      -0.078   0.049  -0.138  1.00  0.00           X\\nATOM    114  X   XXX X 852       0.079  -0.098  -0.069  1.00  0.00           X\\nATOM    115  X   XXX X 852       0.157  -0.098   0.000  1.00  0.00           X\\nATOM    116  X   XXX X 852       0.078  -0.295   0.092  1.00  0.00           X\\nATOM    117  X   XXX X 852       0.157  -0.049   0.069  1.00  0.00           X\\nATOM    118  X   XXX X 852      -0.235   0.000   0.023  1.00  0.00           X\\nATOM    119  X   XXX X 852       0.157  -0.394   0.000  1.00  0.00           X\\nATOM    120  X   XXX X 852      -0.471  -0.197   0.000  1.00  0.00           X\\nATOM    121  X   XXX X 852       0.549  -0.098   0.069  1.00  0.00           X\\nATOM    122  X   XXX X 852      -0.078   0.246   0.046  1.00  0.00           X\\nATOM    123  X   XXX X 852       0.314   0.295  -0.092  1.00  0.00           X\\nATOM    124  X   XXX X 852      -0.078  -0.394  -0.138  1.00  0.00           X\\nATOM    125  X   XXX X 852       0.078  -0.098   0.069  1.00  0.00           X\\nATOM    126  X   XXX X 852      -0.157   0.000  -0.092  1.00  0.00           X\\nATOM    127  X   XXX X 852      -0.235  -0.148   0.046  1.00  0.00           X\\nATOM    128  X   XXX X 852      -0.314   0.049  -0.092  1.00  0.00           X\\nATOM    129  X   XXX X 852       0.235   0.000   0.092  1.00  0.00           X\\nATOM    130  X   XXX X 852      -0.471   0.148  -0.069  1.00  0.00           X\\nATOM    131  X   XXX X 852      -0.235   0.541   0.046  1.00  0.00           X\\nATOM    132  X   XXX X 852      -0.314  -0.295   0.046  1.00  0.00           X\\nATOM    133  X   XXX X 852       0.078   0.049   0.000  1.00  0.00           X\\nATOM    134  X   XXX X 852       0.549   0.098   0.046  1.00  0.00           X\\nATOM    135  X   XXX X 852      -0.078  -0.295   0.046  1.00  0.00           X\\nATOM    136  X   XXX X 852       0.078  -0.246   0.023  1.00  0.00           X\\nATOM    137  X   XXX X 852       0.078   0.295  -0.138  1.00  0.00           X\\nATOM    138  X   XXX X 852       0.314   0.049  -0.115  1.00  0.00           X\\nATOM    139  X   XXX X 852       0.235   0.197   0.046  1.00  0.00           X\\nATOM    140  X   XXX X 852       0.000   0.098  -0.046  1.00  0.00           X\\nATOM    141  X   XXX X 852       0.235  -0.049   0.046  1.00  0.00           X\\nATOM    142  X   XXX X 852      -0.314  -0.246   0.000  1.00  0.00           X\\nATOM    143  X   XXX X 852      -0.236  -0.098   0.115  1.00  0.00           X\\nATOM    144  X   XXX X 852       0.157  -0.344   0.000  1.00  0.00           X\\nATOM    145  X   XXX X 852      -0.314  -0.148   0.069  1.00  0.00           X\\nATOM    146  X   XXX X 852      -0.393   0.148  -0.069  1.00  0.00           X\\nATOM    147  X   XXX X 852       0.236  -0.197  -0.069  1.00  0.00           X\\nATOM    148  X   XXX X 852       0.078  -0.049   0.069  1.00  0.00           X\\nATOM    149  X   XXX X 852       0.000   0.197   0.069  1.00  0.00           X\\nATOM    150  X   XXX X 852       0.392  -0.394  -0.023  1.00  0.00           X\\nATOM    151  X   XXX X 852      -0.078  -0.049  -0.069  1.00  0.00           X\\nATOM    152  X   XXX X 852      -0.000  -0.098   0.138  1.00  0.00           X\\nATOM    153  X   XXX X 852      -2.512  -0.640  -0.621  1.00  0.00           X\\nATOM    154  X   XXX X 852       4.082  -0.295   0.851  1.00  0.00           X  ";

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
                    windowAverageSize: 1,""" + extra_browser_sim_params + """
                });

                browserSim.viewer.addPDBTxt(testPDB);

                browserSim.io.loadJSON("data.json", () => {
                    browserSim.player.start({
                        durationInMilliseconds: 10000,
                        updateFreqInMilliseconds: 16.67,  // 60 fps
                        loop: true,
                        windowAverageSize: 25
                    });

                    setTimeout(() => {
                        browserSim.player.stop();
                        setTimeout(() => {
                            browserSim.player.toFrame(159);
                        }, 1000);
                    }, 5000);
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
        "{sphere: {color: 'green'}},"
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