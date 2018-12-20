declare var jQuery;
declare var THREE;

interface Params {
    "viewer"?: any;
    "viewerType"?: string;
    "visStyle"?: any;
    "durationInMilliseconds"?: number;
    "updateFreqInMilliseconds"?: number;
    "loop"?: boolean;
    "windowAverageSize"?: number;
    "halfWindowSize"?: number;    // Programmatically added.
    "parent"?: any;               // Optional because user doesn't need to
                                  // specify. Always added programmatically.
    "loadPDBTxt"?: any;           // Used if generic iterface
    "updateAtomPositions"?: any;  // Used if generic iterface
    "render"?: any;               // Used if generic iterface
}

class BrowserSim {
    // Note that the below are public so they can be accessed from other
    // classes. But they are not public to the user and so do not need to be
    // protected from closure compiler.

    public _numFrames: number = 0;
    public _frameSize: number = 0;

    // Keys are frame indecies. Values are Float32Arrays with the
    // coefficients.
    public _frameData = {};
    public _frameStride = 1;

    public _numComponents: number = undefined;
    public _componentSize: number = undefined;
    public _componentData = [];

    public _averagePositions = undefined;
    public _params: Params = {};

    /**
     * Constructor for the BrowserSim class.
     * @param  {Object<string,*>} params
     * @returns void
     */
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

    /**
     * Updates the parameters. Sets defaults, then overwrites those with
     * user-specified values (if given). Also calculates some derived values.
     * @param  {Object<string,*>} updatedParams
     * @returns void
     */
    public updateParams(updatedParams: any): void {
        // Set some common error text blurbs.
        let genericModeBut = 'You are running BrowserSim in GENERIC mode, but ' +
                             'you haven\'t specified ';
        let modelVar = 'The "model" variable is whatever your loadPDBTxt(...) ' +
                       'function returns.';
        let viewerVar = 'The "viewer" variable is the user-specified BrowserSim ' +
                        '"viewer" parameter.';
        let browserSimVar = 'The "browserSim" variable is the instantiated ' +
                               'BrowserSim object.';



        // Set the defaults and keeps any previous ones not now specified.
        let defaults: Params = {
            "viewer": undefined,
            "viewerType": undefined,
            "durationInMilliseconds": 10000,
            "updateFreqInMilliseconds": 16.67,  // 60 fps
            "loop": true,
            "windowAverageSize": 1,
            "parent": this,
            "loadPDBTxt": (pdbTxt: string, viewer?: any, browserSim?: any) => {
                throw new Error(
                    genericModeBut +
                    'a loadPDBTxt(pdbTxt, viewer, browserSim) function. This ' +
                    'function loads PDB text into the viewer and optionally ' +
                    'returns a model object. For your reference, the current ' +
                    'value of the "pdbTxt" variable (a string) starts with:\n\n' +
                    pdbTxt.toString().slice(0, 500)) + "\n\n" + viewerVar + "\n\n" +
                    browserSimVar + "\n\n";
            },
            "updateAtomPositions": (newAtomCoors: Float32Array[], model?: any, viewer?: any, browserSim?: any) => {
                let goodCoorRep = "[\n";
                for (let idx in newAtomCoors) {
                    if (newAtomCoors.hasOwnProperty(idx)) {
                        let idxVal = parseInt(idx, 10);
                        let coor = newAtomCoors[idx];
                        goodCoorRep += `  Float32Array [${coor[0]}, ${coor[1]}, ${coor[2]}],\n`;
                        if (idxVal > 10) {
                            break;
                        }
                    }
                }
                goodCoorRep += "  ...\n]";

                throw new Error(
                    genericModeBut + 'an updateAtomPositions(newAtomCoors, ' +
                    'model, viewer, browserSim) function. This function provides ' +
                    'the updated atom coordinates for the current frame. The ' +
                    'value of the "newAtomCoors" variable (a list of Float32Array ' +
                    'containing the new coordinates of the atoms) ' +
                    'looks like:\n\n' + goodCoorRep + '\n\n' + modelVar +
                    "\n\n" + viewerVar + "\n\n" + browserSimVar + "\n\n"
                );
            },
            "render": (model?: any, viewer?: any, browserSim?: any) => {
                throw new Error(
                    genericModeBut + 'a render(model, viewer, browserSim) function. ' +
                    'This function runs every time the atom coordinates ' +
                    'change, to update what the viewer renders. ' + modelVar +
                    ' ' + viewerVar + ' ' + browserSimVar + "\n\n"
                );
            }
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
                case "PV":
                    this._params["visStyle"] = (structure) => { this._params["viewer"]["viewer"]["spheres"]('all', structure); };
                    break;
                case "GENERIC":
                    this._params["visStyle"] = undefined;
                    break;
            }
        }
    }

    /**
     * Get the frame coordinates.
     * @param  {number} frame The frame number.
     * @returns Float32Array  A list of Float32Array containing the atom
     *                        coordinates, for each frame.
     */
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

        // Now average the flattened coordinates over the frame windowe.
        let summedCoorsOverFrames = coorsFlattenedForFrames.reduce((summedCoors, newCoors) => {
            return MathUtils.sumArrayOfFloat32Arrays([summedCoors, newCoors]);
        });

        let averageCoorsOverFrames = MathUtils.multiplyFloat32ArrayByScalar(
            1.0 / framesToAvg.length, summedCoorsOverFrames
        );

        // Add in the average coordinates from the JSON data.
        averageCoorsOverFrames = MathUtils.sumArrayOfFloat32Arrays([
            averageCoorsOverFrames, this._averagePositions
        ]);

        // Reshape the averaged coordinates into a list of Float32Array triplets.
        let coors = MathUtils._range(0, this._componentSize / 3).map(i => {
            return new Float32Array([
                averageCoorsOverFrames[3 * i],
                averageCoorsOverFrames[3 * i + 1],
                averageCoorsOverFrames[3 * i + 2]
            ]);
        })

        return coors;
    }
}

class _IO {
    private _parent = undefined;

    /**
     * The constructor of an _IO class.
     * @param  {*} parent The parent class (BrowserSim Object).
     */
    constructor(parent: BrowserSim) {
        this._parent = parent;
    }

    /**
     * Loads a JSON file written from the BrowserSim python script. Contains
     * all the information needed to visualize the simulation.
     * @param  {string}   path      The JSON path.
     * @param  {Function} callBack  The callback function once loaded.
     * @returns void
     */
    public "loadJSON"(path: string, callBack: any = () => {}): void {
        jQuery.getJSON(path, (data: any) => {
            // console.log(data);

            // Get some general info about the data.
            this._parent._frameSize = data["coeffs"][0].length;
            this._parent._numComponents = data["vecs"].length;

            // The length of the frame coefficients must equal the number of
            // components.
            if (this._parent._frameSize !== this._parent._numComponents) {
                throw new Error(`The number of coefficients per frame ` +
                                `(${this._parent._frameSize}) is not the same` +
                                `as the number of components (${this._parent._numComponents})`);
            }

            // Get the precision of the data (used to eliminate decimal points
            // from the JSON file).
            let precision = data["params"]["precision"];
            let precisionFactor = 10**-precision;

            // Set up the frames
            this._loadJSONSetupFrameCoeffs(data, precisionFactor);

            // Set up the components
            this._loadJSONSetupPCAComponents(data, precisionFactor);

            // Set up the average positions.
            this._loadJSONSetupAveragePos(data, precisionFactor);

            // Here you will put the average positions, but not ready yet.
            // this._parent._averagePositions = new Float32Array(
            //     this._parent._numComponents * this._parent._componentSize
            // );

        }).done(() => {
            callBack();
        }).fail(() => {
            console.log( "error" );
        }); /* .always(function() {
            console.log( "complete" );
        }); */
    };

    private _loadJSONSetupFrameCoeffs(data: any, precisionFactor: number) {
        // Setup the frames
        this._parent._numFrames = data["coeffs"].length;
        this._parent._frameData = {};  // Keys will be frame indexes,
                                        // values will be Float32Array
                                        // containing the PCA coefficients
                                        // for the corresponding frame.
        this._parent._frameStride = 5;  // TODO: SHOULD BE USER DEFINED.

        // Convert frames to array of typed arrays. It's faster.
        for (let idx1 in data["coeffs"]) {
            if (data["coeffs"].hasOwnProperty(idx1)) {
                let idx1Num = parseInt(idx1, 10);
                this._parent._frameData[idx1Num * this._parent._frameStride] = new Float32Array(
                    data["coeffs"][idx1Num]
                );
            }
        }

        // The coefficients need to be converted from int to floats.
        for (let idx1 in this._parent._frameData) {
            if (this._parent._frameData.hasOwnProperty(idx1)) {
                this._parent._frameData[idx1] = this._parent._frameData[idx1].map(v => precisionFactor * v);
            }
        }

        // Now go through and fill in the ones inbetween the explicitly
        // specified frames (with interpolation). Doing linear interpolation
        // for simplicity (rather than spline, for example).
        for (let frameIdx = 0;
             frameIdx < Object.keys(this._parent._frameData).length * this._parent._frameStride;
             frameIdx++) {
            if (this._parent._frameData[frameIdx] === undefined) {
                // This frame isn't defined.

                // Get the previous defined frame.
                let beforeFrame = this._parent._frameStride * Math.floor(
                    frameIdx / this._parent._frameStride
                );

                // Get the next defined frame.
                let afterFrame = beforeFrame + this._parent._frameStride;

                // If either doesn't exist, go to the next one.
                if ((this._parent._frameData[beforeFrame] === undefined) ||
                    (this._parent._frameData[afterFrame] === undefined)) {
                    continue;
                }

                // Get the ratio to interpolate (linear).
                let ratio = (frameIdx - beforeFrame) / this._parent._frameStride;

                // Get the difference between the two vectors.
                let delta = MathUtils.sumArrayOfFloat32Arrays([
                    this._parent._frameData[beforeFrame],
                    MathUtils.multiplyFloat32ArrayByScalar(
                        -1.0, this._parent._frameData[afterFrame]
                    )
                ]);

                // Calculate the coefficients of the undefined frame.
                let newFrame = MathUtils.sumArrayOfFloat32Arrays([
                    this._parent._frameData[beforeFrame],
                    MathUtils.multiplyFloat32ArrayByScalar(ratio, delta)
                ]);

                // Save it to the object.
                this._parent._frameData[frameIdx] = newFrame;
            }
        }
    }

    private _loadJSONSetupPCAComponents(data: any, precisionFactor: number) {
        // Set up the components
        this._parent._componentSize = data["vecs"][0].length;
        this._parent._componentData = []

        // Same with vectors.
        for (let idx1 in data["vecs"]) {
            if (data["vecs"].hasOwnProperty(idx1)) {
                let idx1Num = parseInt(idx1, 10);
                this._parent._componentData[idx1Num] = new Float32Array(data["vecs"][idx1Num]);
            }
        }

        // The vectors need to be converted from int to floats.
        for (let idx1 in this._parent._componentData) {
            if (this._parent._componentData.hasOwnProperty(idx1)) {
                let idxNum = parseInt(idx1, 10);
                this._parent._componentData[idxNum] = this._parent._componentData[idxNum].map(
                    v => precisionFactor * v
                );
            }
        }
    }

    private _loadJSONSetupAveragePos(data: any, precisionFactor: number) {
        // Make the average coordinates, converting them from int to float.
        this._parent._averagePositions = new Float32Array(
            data["coors"].map(v => precisionFactor * v)
        );
    }

    /**
     * Makes a multi-frame PDB file of the simulation. Good for debugging.
     * @returns string The text of the multi-frame PDB file.
     */
    public "makePDB"(): string {
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

    /**
     * Formats a number for use in a PDB line.
     * @param  {number} val The number
     * @returns string  The formatted number (rounded, padded, etc.).
     */
    private _formatNumForPDB(val: number): string {
        let valStr = val.toFixed(3);
        return this._rjust(8, valStr);
    };

    /**
     * Right justifies a string.
     * @param  {number}     length      The length of the string and padding.
     * @param  {string}     origString  The original string.
     * @param  {string=}  padChar     The character to use for padding.
     *                                  Space by default.
     * @returns string      The padded, right-justified string.
     */
    private _rjust(length: number, origString: string, padChar: string = " "): string {
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

    /**
     * The constructor of a _Viewer class.
     * @param  {*} parent The parent class (BrowserSim Object).
     */
    constructor(parent: any) {
        this._parent = parent;

        // Perform some checks.
        if (this._parent._params["viewer"] === undefined) {
            throw new Error("No viewer specified!");
        }

        if (this._parent._params["viewerType"] === undefined) {
            throw new Error("No viewer type specified!");
        }

        let validTypes = ["3DMOLJS", "NGL", "PV", "GENERIC"];
        if (validTypes.indexOf(this._parent._params["viewerType"]) === -1) {
            throw new Error("Specified viewer type, " +
                            this._parent._params["viewerType"] + ", is " +
                            "invalid. Must be one of " + validTypes.join("/"));
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
            case "PV":
                this["addPDBTxt"] = this._PV_AddPDBTxt;
                this._updateAtomPosFun = this._PV_UpdateAtomPos;
                this._render = this._PV_Render;
                break;
            case "GENERIC":
                this["addPDBTxt"] = this._GENERIC_AddPDBTxt;
                this._updateAtomPosFun = this._GENERIC_UpdateAtomPos;
                this._render = this._GENERIC_Render;
                break;
            default:
                // Should never be able to get here because of check above...
                throw new Error("Invalid viewer type specified: " + this._parent._params["viewerType"]);
        }
    }

    /**
     * Update the atom positions for a given frame.
     * @param  {number} frame The frame number.
     * @returns void
     */
    public updateAtomPos(frame: number): void {
        if (this._updateLocked) { return; }
        this._updateLocked = true;
        this._updateAtomPosFun(frame);
        this._updateLocked = false;
        this._render();
    }

    /**
     * The function to add PDB txt to the 3Dmol.js viewer.
     * @param  {string} pdbTxt The PDB text.
     * @returns void
     */
    private _3DMolJS_AddPDBTxt(pdbTxt: string): void {
        this._model = this._parent._params["viewer"].addModel( pdbTxt, "pdb" );
        this._render();
    }

    /**
     * Updates the atom positions in a 3Dmol.js viewer.
     * @param  {number} frame The frame number.
     * @returns void
     */
    private _3DMolJS_UpdateAtomPos(frame: number): void {
        let newAtomCoors = this._parent.getFrameCoors(frame);

        let atoms = this._model.selectedAtoms({});
        for (let atomIdx=0; atomIdx<atoms.length; atomIdx++) {
            let coors = newAtomCoors[atomIdx];
            atoms[atomIdx]["x"] = coors[0];
            atoms[atomIdx]["y"] = coors[1];
            atoms[atomIdx]["z"] = coors[2];
        }
    }

    /**
     * Renders the 3DMol.js viewer.
     * @returns void
     */
    private _3DMolJS_Render(): void {
        // Must update styles to actually have atoms move. Annoying.
        this._model.setStyle({}, this._parent._params["visStyle"]);
        this._parent._params["viewer"].render();
    }

    /**
     * The function to add PDB txt to the NGLViewer.
     * @param  {string} pdbTxt The PDB text.
     * @returns void
     */
    private _NGL_AddPDBTxt(pdbTxt: string): void {
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

    /**
     * Updates the atom positions in a NGLViewer.
     * @param  {number} frame The frame number.
     * @returns void
     */
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

    /**
     * Renders the NGLViewer.
     * @returns void
     */
    private _NGL_Render() {
        this._model["rebuildRepresentations"]();
    }

    /**
     * The function to add PDB txt to the PV Viewer.
     * @param  {string} pdbTxt The PDB text.
     * @returns void
     */
    private _PV_AddPDBTxt(pdbTxt: string): void {
        this._model = this._parent._params["viewer"]["library"]["io"]["pdb"](pdbTxt);
    }

    /**
     * Updates the atom positions in a PV Viewer.
     * @param  {number} frame The frame number.
     * @returns void
     */
    private _PV_UpdateAtomPos(frame: number) {
        let newAtomCoors = this._parent.getFrameCoors(frame);

        let atomIdx = 0;
        this._model["eachAtom"]((atom, idx) => {
            let coors = newAtomCoors[atomIdx];
            atom["_bV"] = coors;
            atomIdx++;
        });
    }

    /**
     * Renders the PV Viewer.
     * @returns void
     */
    private _PV_Render() {
        this._parent._params["viewer"]["viewer"].clear();

        // A callback function in the case of PV.
        this._parent._params["visStyle"](this._model);
    }

    /**
     * The function to add PDB txt if in GENERIC mode.
     * @param  {string} pdbTxt The PDB text.
     * @returns void
     */
    private _GENERIC_AddPDBTxt(pdbTxt: string): void {
        let candidateModel = this._parent._params["loadPDBTxt"](
            pdbTxt, this._parent._params["viewer"], this
        );

        if (candidateModel !== undefined) {
            this._model = candidateModel;
        }
    }

    /**
     * Updates the atom positions if in GENERIC mode.
     * @param  {number} frame The frame number.
     * @returns void
     */
    private _GENERIC_UpdateAtomPos(frame: number) {
        let newAtomCoors = this._parent.getFrameCoors(frame);
        this._parent._params["updateAtomPositions"](
            newAtomCoors, this._model, this._parent._params["viewer"], this
        );
    }

    /**
     * Renders if in GENERIC mode.
     * @returns void
     */
    private _GENERIC_Render() {
        this._parent._params["render"](
            this._model, this._parent._params["viewer"], this
        );
    }
}

class _Player {
    private _animationFrameID = undefined;
    private _parent = undefined;

    /**
     * The constructor of a _Player class.
     * @param  {*} parent The parent class (BrowserSim Object).
     */
    constructor(parent: any) {
        this._parent = parent;
    }

    /**
     * Starts playing the simulation animation.
     * @param  {Object<string,*>} params Parameter from the user.
     * @returns void
     */
    public "start"(params: Params): void {
        // Allow the user to change some of the parameters on start. So
        // duration doesn't always have to be fixed, for example.
        this._parent.updateParams(params);

        // Clear previous requestAnimationFrames.
        if (this._animationFrameID !== undefined) {
            cancelAnimationFrame(this._animationFrameID);
        }

        /**
         * The loop function.
         * @param  {number} timestamp The number of milliseconds since the
         * requestAnimationFrame started.
         */
        let timestampLastFire = 0;
        let loop = (timestamp: number) => {
            let msSinceLastFire = timestamp - timestampLastFire;

            // Enough time has passed that we should update.
            if (msSinceLastFire > this._parent._params["updateFreqInMilliseconds"]) {
                // Save the new timestampLastFire value.
                timestampLastFire = timestamp;

                // How far along the animation are you?
                let playRatio = (timestamp / this._parent._params["durationInMilliseconds"]) % 1.0;

                // If you've gone over the end of the animation, start again from
                // the beginning.
                if ((!this._parent._params["loop"]) && (playRatio > 1.0)) {
                    this["stop"]();
                    return;
                }

                // Get the current frame.
                let curFrame = Math.floor(this._parent._numFrames * playRatio);

                this._parent["viewer"].updateAtomPos(curFrame);
            }

            // Start the next iteration.
            requestAnimationFrame(loop);
        }

        this._animationFrameID = requestAnimationFrame(loop);
    }

    /**
     * Stops the simulation animation.
     * @returns void
     */
    public "stop"(): void {
        // Clear previous requestAnimationFrames.
        if (this._animationFrameID !== undefined) {
            cancelAnimationFrame(this._animationFrameID);
        }
    }

    /**
     * Switches the simulation animation to a specific frame.
     * @param  {number} frame The frame.
     * @returns void
     */
    public "toFrame"(frame: number): void {
        this._parent["viewer"].updateAtomPos(frame);
    }
}

module MathUtils {
    /**
     * Gets a list of numbers. Like Python's range function.
     * @param  {number} a     The first number.
     * @param  {number} b     The largest numer (plus one).
     * @returns Array<number> A list of the numbers.
     */
    export function _range(a: number, b: number): number[] {
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

    /**
     * Given a list of Float32Array, returns the sum of them (reduce).
     * @param  {Array<Float32Array>} arrayOfFloat32Arrays The list of Float32Array.
     * @returns Float32Array         The sum.
     */
    export function sumArrayOfFloat32Arrays(arrayOfFloat32Arrays: Float32Array[]): Float32Array {
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

    /**
     * Multiplies a Float32Array by a scalar.
     * @param  {number}       scalar        The scalar.
     * @param  {Float32Array} float32Array  The Float32Array.
     * @returns Float32Array  The multiplied Float32Array.
     */
    export function multiplyFloat32ArrayByScalar(scalar: number, float32Array: Float32Array): Float32Array {
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

// Polyfill requestAnimationFrame
// http://paulirish.com/2011/requestanimationframe-for-smart-animating/
// http://my.opera.com/emoller/blog/2011/12/20/requestanimationframe-for-smart-er-animating
// requestAnimationFrame polyfill by Erik Möller. fixes from Paul Irish and Tino Zijdel
// MIT license
(function() {
    var lastTime = 0;
    var vendors = ['ms', 'moz', 'webkit', 'o'];
    for(var x = 0; x < vendors.length && !window["requestAnimationFrame"]; ++x) {
        window["requestAnimationFrame"] = window[vendors[x]+'RequestAnimationFrame'];
        window["cancelAnimationFrame"] = window[vendors[x]+'CancelAnimationFrame']
                                         || window[vendors[x]+'CancelRequestAnimationFrame'];
    }

    if (!window["requestAnimationFrame"])
        window["requestAnimationFrame"] = function(callback, element) {
            var currTime = new Date().getTime();
            var timeToCall = Math.max(0, 16 - (currTime - lastTime));
            var id = window.setTimeout(function() {
                callback(currTime + timeToCall);
            }, timeToCall);
            lastTime = currTime + timeToCall;
            return id;
        };

    if (!window["cancelAnimationFrame"])
        window["cancelAnimationFrame"] = function(id) {
            clearTimeout(id);
        };
}());

(<any>window)["BrowserSim"] = BrowserSim;  // To survive closure compiler.
