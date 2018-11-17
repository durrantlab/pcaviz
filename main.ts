declare var jQuery;

interface Params {
    "viewer"?: any;
    "viewerType"?: string;
    "visStyle"?: any;
    "durationInMilliseconds"?: number;
    "updateFreqInMilliseconds"?: number;
    "loop"?: boolean;
    "windowAverageSize"?: number;
    "halfWindowSize"?: number;  // Programmatically added.
    "parent"?: any;  // Optional because user doesn't need to specify. Always added programmatically.
}

class BrowserSim {
    // Note that the below are public so they can be accessed from other
    // classes. But they are not public to the user and so do not need to be
    // protected from closure compiler.
    public _numFrames: number = undefined;
    public _frameSize: number = undefined;
    public _frameData = [];

    public _numComponents: number = undefined;
    public _componentSize: number = undefined;
    public _componentData = [];

    public _averagePositions = undefined;
    public _params: Params = {};

    constructor(params: Params) {
        this["io"] = new _IO(this);  // This way so it survives closure compiler.

        this.updateParams(params);

        // Note that viewer provides a common interface for several
        // browser-based visualization libraries: (3DMol.js, ETC, ETC.) It is
        // optional. If included, this script will have greater control over
        // the visualization. But it can be ommited if the web developer wants
        // to take a more customized approach.
        if (params["viewer"] !== undefined) {
            this["viewer"] = new _Viewer(this);
            this["player"] = new _Player(this);
        }
    }

    public updateParams(updatedParams: any) {
        // Set the defaults and keeps any previous ones not now specified.
        let defaults: Params = {
            "viewer": undefined,
            "viewerType": undefined,
            "durationInMilliseconds": 10000,
            "updateFreqInMilliseconds": 10,
            "loop": true,
            "windowAverageSize": 1,
            "parent": this
        }

        this._params = jQuery.extend(defaults, this._params, updatedParams);

        // halfWindowSize is always derived from windowAverageSize
        this._params["halfWindowSize"] = Math.floor((this._params["windowAverageSize"] - 1) / 2);

        // Default visStyle depends on viewerType
        if (this._params["visStyle"] === undefined) {
            switch(this._params["viewerType"]) {
                case "3DMOLJS":
                    this._params["visStyle"] = {"line": {}};
                    break;
                case "NGL":
                    this._params["visStyle"] = undefined;
                    break;
                case "JSMOL":
                    this._params["visStyle"] = undefined;
                    break;
            }
        }
    }

    public getFrameCoors(frame: number): Float32Array[] {
        // Consider multiple frames if necessary.
        let framesToAvg = MathUtils._range(
            frame - this._params["halfWindowSize"],
            frame + this._params["halfWindowSize"] + 1
        );

        // Get the coefficients for each of those frames (an array of
        // Float32Array's)
        let framesCoefficients = framesToAvg.map((frameIdx): any => {
            // Make sure frame never out of bounds.
            if (frameIdx < 0) { frameIdx += this._numFrames; }
            frameIdx = frameIdx % this._numFrames;

            // Get the frame data (PCA coefficients).
            return this._frameData[frameIdx];
        });

        // Get the flattened coordinates for the frames (i.e., coordinates not
        // separated into triplets).
        let coorsFlattenedForFrames = framesCoefficients.map((frameCoefficients) => {
            // frameCoefficients is a Float32Array containing all the
            // coefficients for a given frame. Need to multiple those by the
            // corresponding components.

            let multipledComponents = this._componentData.map((component, componentIdx) => {
                return MathUtils.multiplyFloat32ArrayByScalar(
                    frameCoefficients[componentIdx], component
                );
            });

            let summedMultipledComponents = MathUtils.sumArrayOfFloat32Arrays(multipledComponents);

            return summedMultipledComponents;
        });

        // Now average the flattened coordinates over the frames.
        let summedCoorsOverFrames = coorsFlattenedForFrames.reduce((summedCoors, newCoors) => {
            return MathUtils.sumArrayOfFloat32Arrays([summedCoors, newCoors]);
        });

        let averageCoorsOverFrames = MathUtils.multiplyFloat32ArrayByScalar(
            1.0 / framesToAvg.length, summedCoorsOverFrames
        );

        // Reshape the averaged coordinates into a list of Float32Array triplets.
        let coors = MathUtils._range(0, this._componentSize / 3).map(i => {
            return new Float32Array([
                averageCoorsOverFrames[i],
                averageCoorsOverFrames[i + 1],
                averageCoorsOverFrames[i + 2]
            ]);
        })

        return coors;
    }
}

class _IO {
    private _parent = undefined;
    constructor(parent) {
        this._parent = parent;
    }

    public "loadJSON"(filename: string, callBack: any = () => {}) {
        jQuery.getJSON(filename, (data) => {
            // Setup the frames
            this._parent._numFrames = data["frames"].length;
            this._parent._frameSize = data["frames"][0].length;
            this._parent._frameData = [];

            // Set up the components
            this._parent._numComponents = data["vectors"].length;
            this._parent._componentSize = data["vectors"][0].length;
            this._parent._componentData = []

            // Here you will put the average positions, but not ready yet.
            this._parent._averagePositions = new Float32Array(
                this._parent._numComponents * this._parent._componentSize
            );

            // The length of the frame coefficients must equal the number of
            // components.
            if (this._parent._frameSize !== this._parent._numComponents) {
                throw new Error(`The number of coefficients per frame ` +
                                `(${this._parent._frameSize}) is not the same` +
                                `as the number of components (${this._parent._numComponents})`);
            }

            // Convert frames to array of typed arrays. It's faster.
            for (let idx1 in data["frames"]) {
                if (data["frames"].hasOwnProperty(idx1)) {
                    this._parent._frameData[idx1] = new Float32Array(data["frames"][idx1]);
                }
            }

            // Same with vectors.
            for (let idx1 in data["vectors"]) {
                if (data["vectors"].hasOwnProperty(idx1)) {
                    this._parent._componentData[idx1] = new Float32Array(data["vectors"][idx1]);
                }
            }
        }).done(() => {
            callBack();
        }).fail(() => {
            console.log( "error" );
        }); /* .always(function() {
            console.log( "complete" );
        }); */
    };

    public "makePDB"() {
        let frameIdx = 0;
        let numAtoms = this._parent._componentSize / 3;
        let pdbTxt = "";
        for (let frameIdx=0; frameIdx<this._parent._numFrames; frameIdx++) {
            pdbTxt += "MODEL " + frameIdx.toString() + "\n";
            for (let atomIdx=0; atomIdx<numAtoms; atomIdx++) {
                let coor = this._parent.getAtomCoors(frameIdx, atomIdx);
                pdbTxt += "ATOM  " + this._rjust(5, atomIdx.toString()) +
                          "  X   XXX X 852    " +
                          this._formatNumForPDB(coor[0]) +
                          this._formatNumForPDB(coor[1]) +
                          this._formatNumForPDB(coor[2]) +
                          "  1.00  0.00           X  \n";
            }
            pdbTxt += "ENDMDL\n";
        }

        return pdbTxt;
    };

    private _formatNumForPDB(val: number) {
        let valStr = val.toFixed(3);
        return this._rjust(8, valStr);
    };

    private _rjust(length: number, origString: string, padChar: string = " ") {
        while (origString.length < length) {
            origString = padChar + origString;
        }
        return origString;
    };
}

class _Viewer {
    private _model = undefined;
    private _updateLocked: boolean = false;  // Don't update if currently updating.
    private _updateAtomPosFun = undefined;
    private _render = undefined;
    private _parent = undefined;

    constructor(parent: any) {
        this._parent = parent;

        // Perform some checks.
        if (this._parent._params["viewer"] === undefined) {
            throw new Error("No viewer specified!");
        }

        if (this._parent._params["viewerType"] === undefined) {
            throw new Error("No viewer type specified!");
        }

        let validTypes = ["3DMOLJS", "NGL", "JSMOL"];
        if (validTypes.indexOf(this._parent._params["viewerType"]) === -1) {
            throw new Error("Specified viewer type, " + this._parent._params["viewerType"] + ", is invalid. Must be one of " + validTypes.join("/"));
        }

        // Rename some functions. This is to avoid having to check
        // this._viewerType every time.
        switch (this._parent._params["viewerType"]) {
            case "3DMOLJS":
                this["addPDBTxt"] = this._3DMolJS_AddPDBTxt;
                this._updateAtomPosFun = this._3DMolJS_UpdateAtomPos;
                this._render = this._3DMolJS_Render;
                break;
            case "NGL":
                this["addPDBTxt"] = this._NGL_AddPDBTxt;
                this._updateAtomPosFun = this._NGL_UpdateAtomPos;
                this._render = this._NGL_Render;
                break;
            case "JSMOL":
                this["addPDBTxt"] = this._JSMOL_AddPDBTxt;
                this._updateAtomPosFun = this._JSMOL_UpdateAtomPos;
                this._render = this._JSMOL_Render;
                break;
            default:
                // Should never be able to get here because of check above...
                throw new Error("Invalid viewer type specified: " + this._parent._params["viewerType"]);
        }
    }

    public updateAtomPos(frame: number) {
        if (this._updateLocked) { return; }
        this._updateLocked = true;
        this._updateAtomPosFun(frame);
        this._updateLocked = false;
        this._render();
    }

    private _3DMolJS_AddPDBTxt(pdbTxt: string) {
        this._model = this._parent._params["viewer"].addModel( pdbTxt, "pdb" );
        this._render();
    }

    private _3DMolJS_UpdateAtomPos(frame: number) {
        let newAtomCoors = this._parent.getFrameCoors(frame);

        let atoms = this._model.selectedAtoms({});
        for (let atomIdx=0; atomIdx<atoms.length; atomIdx++) {
            let coors = newAtomCoors[atomIdx];
            atoms[atomIdx]["x"] = coors[0];
            atoms[atomIdx]["y"] = coors[1];
            atoms[atomIdx]["z"] = coors[2];
        }
    }

    private _3DMolJS_Render() {
        // Must update styles to actually have atoms move. Annoying.
        this._model.setStyle({}, this._parent._params["visStyle"]);
        this._parent._params["viewer"].render();
    }

    private _NGL_AddPDBTxt(pdbTxt: string) {
        this._parent._params["viewer"].loadFile(
            new Blob(
                [pdbTxt],
                {type: 'text/plain'}
            ),
            {ext:'pdb', defaultRepresentation: true}
        ).then((e) => {
            this._model = e;
        });
    }

    private _NGL_UpdateAtomPos(frame: number) {
        let newAtomCoors = this._parent.getFrameCoors(frame);

        let atomIdx = 0;
        this._model["structure"]["eachAtom"]((atom, idx) => {
            let coors = newAtomCoors[atomIdx];
            atom["x"] = coors[0];
            atom["y"] = coors[1];
            atom["z"] = coors[2];
            atomIdx++;
        });
    }

    private _NGL_Render() {
        this._model["rebuildRepresentations"]();
    }

    private _JSMOL_AddPDBTxt(pdbTxt: string) {
        // let jsmolCmd = "data \"model browser_sim\"\n";
        // jsmolCmd += pdbTxt + "\n";
        // jsmolCmd += 'end "model browser_sim";show data "pdb"';

        // See
        // https://jmol-developers.narkive.com/UJZfiMkC/ie-specific-problem-with-jmolloadinline
        // for example.

        let test = `2
testing
C 1 1 1
O 2 2 2
`;

        // jQuery("body").append("<div id='test'>" + test + "</div>");
        setTimeout(() => {
            // window["t"] = test.split("\n");
            // let jsmolCmd = "jmolLoadInline('" + pdbTxt.replace(/\n/g, "\\n") + "');";

            let jsmolCmd = `load data "model example"
ATOM     30  X   XXX X 852       1.020  -0.246  -0.345  1.00  0.00           X
ATOM     31  X   XXX X 852      -1.570   0.836  -0.115  1.00  0.00           X
ATOM     32  X   XXX X 852       0.000  -0.738   0.207  1.00  0.00           X
ATOM     33  X   XXX X 852       0.785  -0.344   0.069  1.00  0.00           X
ATOM     34  X   XXX X 852      -1.727  -0.246   0.253  1.00  0.00           X
ATOM     35  X   XXX X 852      -0.157  -0.197   0.000  1.00  0.00           X
ATOM     36  X   XXX X 852       0.549   0.984  -0.391  1.00  0.00           X
ATOM     37  X   XXX X 852       0.471   0.295   0.138  1.00  0.00           X
end "model example";`
            console.log(jsmolCmd);
            // jsmolCmd = "background red";
            this._parent._params["viewer"]["library"].script(
                this._parent._params["viewer"]["applet"],
                jsmolCmd
            );
        }, 1000);

        // this._model = this._parent._params["viewer"].addModel( pdbTxt, "pdb" );
        // this._render();

        // viewerJSMol = {
        //     applet: jsmolApplet,
        //     library: Jmol
        // }

    }

    private _JSMOL_UpdateAtomPos(frame: number) {

    }

    private _JSMOL_Render() {

    }
}

class _Player {
    private _startTime = undefined;
    private _timer = undefined;
    private _parent = undefined;

    constructor(parent: any) {
        this._parent = parent;
    }

    public "start"(params: Params) {
        // loop: boolean = true, windowAverageSize: number = 1) {
        this._startTime = new Date().getTime();

        // Allow the user to change some of the parameters on start. So
        // duration doesn't always have to be fixed, for example.
        this._parent.updateParams(params);

        // Clear previous intervals.
        if (this._timer !== undefined) {
            clearInterval(this._timer);
        }

        this._timer = setInterval(() => {
            // How far along the animation are you?
            let curTime = new Date().getTime();
            let playRatio = (curTime - this._startTime) / this._parent._params["durationInMilliseconds"];

            // If you've gone over the end of the animation, start again from
            // the beginning.
            if ((!this._parent._params["loop"]) && (playRatio > 1.0)) {
                this["stop"]();
                return;
            }

            // Get the current frame.
            let curFrame = Math.floor(this._parent._numFrames * playRatio);

            this._parent["viewer"].updateAtomPos(curFrame);

        }, this._parent._params["updateFreqInMilliseconds"]);
    }

    public "stop"() {
        // Clear previous intervals.
        if (this._timer !== undefined) {
            clearInterval(this._timer);
        }
    }

    public "toFrame"(frame: number) {
        this._parent["viewer"].updateAtomPos(frame);
    }
}

module MathUtils {
    export function _range(a: number, b: number) {
        // Don't do this as a Float32Array, because map won't work as
        // expected. Map on typed array must return typed array.
        let rng = [];
        let i = a;
        while (i < b) {
            rng.push(i);
            i++
        }
        return rng;
    }

    export function sumArrayOfFloat32Arrays(arrayOfFloat32Arrays: Float32Array[]) {
        if (arrayOfFloat32Arrays.length === 1) {
            // Just one item in the array, so return that first item.
            return arrayOfFloat32Arrays[0];
        } else {
            // Multiple items. So need to sum them.
            return arrayOfFloat32Arrays.reduce((summedVals, newVals) => {
                return summedVals.map((v, i) => {
                    return v + newVals[i];
                })
            });
        }
    }

    export function multiplyFloat32ArrayByScalar(scalar: number, float32Array: Float32Array) {
        switch (scalar) {
            case 0.0:
                return new Float32Array(float32Array.length);
            case 1.0:
                return float32Array;
            default:
                return float32Array.map(v => scalar * v);
        }
    }
}

(<any>window)["BrowserSim"] = BrowserSim;  // To survive closure compiler.
