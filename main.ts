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
    public _frameData = undefined;

    public _numComponents: number = undefined;
    public _componentSize: number = undefined;
    public _componentData = undefined;

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
            "visStyle": {"line": {}},
            "durationInMilliseconds": 10000,
            "updateFreqInMilliseconds": 10,
            "loop": true,
            "windowAverageSize": 1,
            "parent": this
        }
        this._params = jQuery.extend(defaults, this._params, updatedParams);

        // halfWindowSize is always derived from windowAverageSize
        this._params["halfWindowSize"] = Math.floor((this._params["windowAverageSize"] - 1) / 2);

    }

    private _range(a: number, b: number) {
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

    public getFrameCoors(frame: number) {

        // Consider multiple frames if necessary.
        let framesToAvg = this._range(
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
            return this._frameData.slice(frameIdx * this._frameSize, (frameIdx + 1) * this._frameSize);
        });

        framesCoefficients.map((frameCoefficients) => {
            // frameCoefficients is a Float32Array containing all the
            // coefficients for a given frame. Need to multiple those by the
            // corresponding components.

            let multipledComponents = this._componentData.map((val, idx) => {
                let componentIdx = Math.floor(idx / this._componentSize);
                console.log(componentIdx);

            });


            // let component = this._componentData[idx];
            // debugger;


        });

        // Go through each framesToAvg and multiply components
        for (let idx in framesCoefficients) {
            if (framesCoefficients.hasOwnProperty(idx)) {
                let framesCoefficient = framesCoefficients[idx];

                [...Array(framesToAvg.length)].map((frameCoeffIdx) => {
                    let frameCoeff = framesCoefficient[frameCoeffIdx];
                    return
                });

            }
        }

        // Make an array for the XYZ of each atom. Put the average position in that.
        // let atomCoors = [...Array(this._componentSize / 3)].map((i) => new Float32Array([0, 0, 0]));


    }

    public getAtomCoors(frame: number, atomIdx: number) {
        // halfWindowSize

        // Make sure frame never out of bounds.
        frame = frame % this._numFrames;

        let thisFrameData = this._frameData.slice(frame * this._frameSize, (frame + 1) * this._frameSize);

        let x = 0;
        let y = 0;
        let z = 0;
        for (let componentIdx=0; componentIdx<this._numComponents; componentIdx++) {
            let idx = componentIdx * this._componentSize + 3 * atomIdx;
            x += thisFrameData[0] * this._componentData[idx];
            y += thisFrameData[1] * this._componentData[idx + 1];
            z += thisFrameData[2] * this._componentData[idx + 2];
        }

        // Remember to add in average position.
        x += this._averagePositions[3 * atomIdx]
        y += this._averagePositions[3 * atomIdx + 1]
        z += this._averagePositions[3 * atomIdx + 2]

        if (isNaN(x) || isNaN(y) || isNaN(z)) {
            // Note that the frame can never be invalid because it wraps
            // around with mod.

            let msg = "Warning: Invalid atom index: " + atomIdx.toString();
            console.log(msg);
        }

        return new Float32Array([x, y, z]);
    };
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
            this._parent._frameData = new Float32Array(
                this._parent._numFrames * this._parent._frameSize
            );

            // Set up the components
            this._parent._numComponents = data["vectors"].length;
            this._parent._componentSize = data["vectors"][0].length;
            this._parent._componentData = new Float32Array(
                this._parent._numComponents * this._parent._componentSize
            );

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
                    let idx1Int = parseInt(idx1, 10);
                    for (let idx2 in data["frames"][idx1]) {
                        if (data["frames"][idx1].hasOwnProperty(idx2)) {
                            let idx2Int = parseInt(idx2, 10);
                            this._parent._frameData[idx1Int * this._parent._frameSize + idx2Int] = data["frames"][idx1Int][idx2Int];
                        }
                    }
                }
            }

            // Same with vectors.
            for (let idx1 in data["vectors"]) {
                if (data["vectors"].hasOwnProperty(idx1)) {
                    let idx1Int = parseInt(idx1, 10);
                    for (let idx2 in data["vectors"][idx1]) {
                        if (data["vectors"][idx1].hasOwnProperty(idx2)) {
                            let idx2Int = parseInt(idx2, 10);
                            this._parent._componentData[idx1Int * this._parent._componentSize + idx2Int] = data["vectors"][idx1Int][idx2Int];
                        }
                    }
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

        let validTypes = ["3DMOLJS"];
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
            default:
                // Should never be able to get here because of check above...
                throw new Error("Invalid viewer type specified: " + this._parent._params["viewerType"]);
        }
    }

    public updateAtomPos(frame: number, halfWindowSize: number = 0) {
        if (this._updateLocked) { return; }
        this._updateLocked = true;
        this._updateAtomPosFun(frame);
        this._updateLocked = false;
        this._render();
    }

    private _3DMolJS_AddPDBTxt(pdbTxt: string) {
        this._model = this._parent._params["viewer"].addModel( pdbTxt, "pdb" );
        this._render();

        // this._model.setStyle({}, this._parent._params["visStyle"]); // {"sphere": {"color": 'spectrum'}});
        // this._parent._params["viewer"].render();
    }

    private _3DMolJS_UpdateAtomPos(frame: number, halfWindowSize: number = 0) {
        this._parent.getFrameCoors(frame);
        debugger;

        let atoms = this._model.selectedAtoms({});

        for (let atomIdx=0; atomIdx<atoms.length; atomIdx++) {
            let coors = this._parent._params["parent"].getAtomCoors(frame, atomIdx);
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
