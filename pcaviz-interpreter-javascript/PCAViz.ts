declare var THREE;
declare var process;
declare var ga;

interface Params {
    "viewer"?: any;
    "viewerType"?: string;
    "visStyle"?: any;
    "durationInMilliseconds"?: number;
    "updateFreqInMilliseconds"?: number;
    "loop"?: boolean;
    "windowAverageSize"?: number;
    "caching"?: string;  // none, continuous, or pre
    "cacheModeNum"?: number;  // A numerical representation of "caching", for
                              // speed.
    "halfWindowSize"?: number;    // Programmatically added.
    "parent"?: any;               // Optional because user doesn't need to
                                  // specify. Always added programmatically.
    "loadPDBTxt"?: any;           // Used if generic interface
    "updateAtomPositions"?: any;  // Used if generic interface
    "render"?: any;               // Used if generic interface
    "playerControlsID"?: string;  // Add player pcaviz-controls to div with this ID.
    "analytics"?: boolean;        // Whether to collect analytics data.
}

// Put it all in a big namespace to avoid having all the classes be global.
namespace PCAVizNameSpace {
    export class PCAViz {
        // Note that the below are public so they can be accessed from other
        // classes. But they are not public to the user and so do not need to be
        // protected from closure compiler.

        public _numFramesInJSON: number = 0;
        public _numFramesTotal: number = 0;  // _numFramesInJSON * _frameStride.
        public _frameSize: number = 0;

        // Keys are frame indecies. Values are Float32Arrays with the
        // coefficients.
        public _frameData = {};
        public _frameStride = 1;

        public _numComponents: number = undefined;
        public _componentSize: number = undefined;
        public _componentData = [];

        public _averagePositions = undefined;
        public _firstFramePositions = undefined;
        public _params: Params = {};

        public _res_info = undefined;

        public _cachedFrameCoors = {};

        private _playerControls;

        /**
         * Constructor for the PCAViz class.
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

            this._playerControls = new _PlayerControls(this);
        }

        /**
         * Updates the parameters. Sets defaults, then overwrites those with
         * user-specified values (if given). Also calculates some derived values.
         * @param  {Object<string,*>} updatedParams
         * @returns void
         */
        public updateParams(updatedParams: any): void {
            // Set some common error text blurbs.
            let genericModeBut = 'You are running PCAViz in GENERIC mode, but ' +
                                 'you haven\'t specified ';
            let modelVar = 'The "model" variable is whatever your loadPDBTxt(...) ' +
                           'function returns.';
            let viewerVar = 'The "viewer" variable is the user-specified PCAViz ' +
                            '"viewer" parameter.';
            let pcaVizVar = 'The "pcaViz" variable is the instantiated ' +
                                   'PCAViz object.';

            // Set the defaults and keeps any previous ones not now specified.
            let defaults: Params = {
                "viewer": undefined,
                "viewerType": undefined,
                "durationInMilliseconds": 10000,
                "updateFreqInMilliseconds": 16.67,  // 60 fps
                "loop": true,
                "caching": "none",
                "cacheModeNum": 0,  // corresponds to caching none.
                "windowAverageSize": 1,
                "playerControlsID": "",  // Player pcaviz-controls by default.
                "parent": this,
                "analytics": true,  // Collect analytics by default
                "loadPDBTxt": (pdbTxt: string, viewer?: any, pcaViz?: any) => {
                    throw new Error(
                        genericModeBut +
                        'a loadPDBTxt(pdbTxt, viewer, pcaViz) function. This ' +
                        'function loads PDB text into the viewer and optionally ' +
                        'returns a model object. For your reference, the current ' +
                        'value of the "pdbTxt" variable (a string) starts with:\n\n' +
                        pdbTxt.toString().slice(0, 500)) + "\n\n" + viewerVar + "\n\n" +
                        pcaVizVar + "\n\n";
                },
                "updateAtomPositions": (newAtomCoors: Float32Array[], model?: any, viewer?: any, pcaViz?: any) => {
                    let goodCoorRep = "[\n";
                    for (let idx in newAtomCoors) {
                        if (newAtomCoors.hasOwnProperty(idx)) {
                            let idxVal = +idx;
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
                        'model, viewer, pcaViz) function. This function provides ' +
                        'the updated atom coordinates for the current frame. The ' +
                        'value of the "newAtomCoors" variable (a list of Float32Array ' +
                        'containing the new coordinates of the atoms) ' +
                        'looks like:\n\n' + goodCoorRep + '\n\n' + modelVar +
                        "\n\n" + viewerVar + "\n\n" + pcaVizVar + "\n\n"
                    );
                },
                "render": (model?: any, viewer?: any, pcaViz?: any) => {
                    throw new Error(
                        genericModeBut + 'a render(model, viewer, pcaViz) function. ' +
                        'This function runs every time the atom coordinates ' +
                        'change, to update what the viewer renders. ' + modelVar +
                        ' ' + viewerVar + ' ' + pcaVizVar + "\n\n"
                    );
                }
            }

            this._params = jjQuery.extend(defaults, this._params, updatedParams);

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

            // Modify the updateFreqInMilliseconds parameter so FPS always <= 60
            // (browser limitation anyway).
            let fpsSixtyUpdateFreqInMilliseconds = 1000.0/60.0;
            if (this._params["updateFreqInMilliseconds"] < fpsSixtyUpdateFreqInMilliseconds) {
                this._params["updateFreqInMilliseconds"] = fpsSixtyUpdateFreqInMilliseconds;
            }

            // Slow playing animations may not be able to maintain even that FPS.
            // If the number of frames is available, further modify if necessary.
            if (this._numFramesInJSON > 0) {
                // Note that there are 0 frames by default. So if 0, there is no
                // frame data yet.
                let fpsPerSim = this._params["durationInMilliseconds"] / this._numFramesTotal;
                if (this._params["updateFreqInMilliseconds"] < fpsPerSim) {
                    this._params["updateFreqInMilliseconds"] = fpsPerSim;
                }
            }

            // Set the cacheModeNum based on the caching
            switch(this._params["caching"]) {
                case "none":
                    this._params["cacheModeNum"] = 0;
                    break;
                case "continuous":
                    this._params["cacheModeNum"] = 1;
                    break;
                case "pre":
                    this._params["cacheModeNum"] = 2;
                    break;
                default:
                    throw new Error(
                        `Invalid caching value: ${this._params["caching"]}. Valid ` +
                        `values are "none", "continuous", and "pre".`
                    );
            }

            // Pre-cache if necessary.
            this.cacheAllFrameCoorsIfNeeded();

            // Inform google analytics as appropriate
            if (this._params["analytics"]) {
                this.googleAnalytics();
            }
        }

        /**
         * Get the frame coordinates.
         * @param  {number} frame The frame number.
         * @returns Float32Array  A list of Float32Array containing the atom
         *                        coordinates, for each frame.
         */
        public getFrameCoors(frame: number): Float32Array[] {
            let coors = Array();

            // console.log(this._params["cacheModeNum"], this._cachedFrameCoors, frame, this._cachedFrameCoors[frame]);

            if ((this._params["cacheModeNum"] === 0) || (this._cachedFrameCoors[frame] === undefined)) {
                // So either cache is not turned on (none), or there's no cached
                // data for this frame.

                // Consider multiple frames if necessary.
                let framesToAvg = MathUtils._range(
                    frame - this._params["halfWindowSize"],
                    frame + this._params["halfWindowSize"] + 1
                );

                // Get the coefficients for each of those frames (output is an
                // array of Float32Array's). Note that map is slow, so using
                // standard for loop.
                let len = framesToAvg.length;
                let framesCoefficients = Array(len);  // An array of Float32Array (the coeff).
                for (let i = 0; i < len; i++) {
                    let frameIdx = framesToAvg[i];
                    // Make sure frame never out of bounds.
                    if (frameIdx < 0) { frameIdx += this._numFramesTotal; }
                    frameIdx = frameIdx % this._numFramesTotal;

                    // Get the frame data (PCA coefficients).
                    framesCoefficients[i] = this._frameData[frameIdx];
                }

                // Get the flattened coordinates for the frames (i.e., coordinates
                // not separated into triplets). Intentionally avoiding map()
                // functions here because they are slow.
                len = framesCoefficients.length;
                let len2 = this._componentData.length;
                let coorsFlattenedForFrames = Array(len);  // A list of Float32Array (flattened coors)
                for (let i = 0; i < len; i++) {
                    let frameCoefficients = framesCoefficients[i];
                    // frameCoefficients is a Float32Array containing all the
                    // coefficients for a given frame. Need to multiple those by the
                    // corresponding components.

                    let multipledComponents = Array(len2);
                    for (let componentIdx = 0; componentIdx < len2; componentIdx++) {
                        let component = this._componentData[componentIdx];
                        multipledComponents[componentIdx] = MathUtils.multiplyFloat32ArrayByScalar(
                            frameCoefficients[componentIdx], component
                        );
                    }

                    let summedMultipledComponents = MathUtils.sumArrayOfFloat32Arrays(multipledComponents);

                    coorsFlattenedForFrames[i] = summedMultipledComponents;
                }

                // Now average the flattened coordinates over the frame window.
                // Produces a long Float32array.
                let summedCoorsOverFrames = MathUtils.sumArrayOfFloat32Arrays(coorsFlattenedForFrames);

                let fac = 1.0 / framesToAvg.length;
                let averageCoorsOverFrames = MathUtils.multiplyFloat32ArrayByScalar(
                    fac, summedCoorsOverFrames
                );

                // Add in the average coordinates from the JSON data.
                averageCoorsOverFrames = MathUtils.sumArrayOfFloat32Arrays([
                    averageCoorsOverFrames, this._averagePositions
                ]);

                // Reshape the averaged coordinates into a list of Float32Array
                // triplets. Avoiding map() funcitons on purpose.
                let arange = MathUtils._range(0, this._componentSize / 3);
                len = arange.length;
                coors = Array(len);
                for (let i2 = 0; i2 < len; i2++) {
                    let i = arange[i2];
                    let three_i = 3 * i;
                    coors[i2] = new Float32Array([
                        averageCoorsOverFrames[three_i],
                        averageCoorsOverFrames[three_i + 1],
                        averageCoorsOverFrames[three_i + 2]
                    ]);
                }

                if (this._params["cacheModeNum"] > 0) {
                    this._cachedFrameCoors[frame] = coors;
                }
            } else {
                // Get the cached version.
                coors = this._cachedFrameCoors[frame];
                // console.log("cached");
            }
            return coors;
        }

        /**
         * Loads all frame data if the user has requested pre-caching.
         * @returns void
         */
        public cacheAllFrameCoorsIfNeeded(): void {
            if (this._params["cacheModeNum"] === 2) {
                this._cachedFrameCoors = {};  // Reset everything.
                for (let frameIdx = 0; frameIdx < this._numFramesTotal; frameIdx++) {
                    this.getFrameCoors(frameIdx);
                }
            }
        }

        private analyticsAlreadyCalled = false;  // To make sure called only once.

        /**
         * Record analytics if user allows.
         * @returns void
         */
        private googleAnalytics(): void {
            if (!this.analyticsAlreadyCalled) {
                this.analyticsAlreadyCalled = true;
                (function(i, s, o, g, r, a, m) {
                    i['GoogleAnalyticsObject'] = r;
                    i[r] = i[r] || function() {
                        (i[r].q = i[r].q || []).push(arguments)
                    }, i[r].l = 1 * new Date().getTime();
                    a = s.createElement(o),
                        m = s.getElementsByTagName(o)[0];
                    a.async = 1;
                    a.src = g;
                    m.parentNode.insertBefore(a, m)
                })(window, document, 'script', 'https://www.google-analytics.com/analytics.js', 'ga');
                ga('create', 'UA-144382730-1', {
                    'name': 'pcaviz'
                });
                // UA-144382730-1 reports to pcaviz
                ga('pcaviz.send', 'pageview');
            }
        }
    }

    class _IO {
        private _parent = undefined;

        /**
         * The constructor of an _IO class.
         * @param  {*} parent The parent class (PCAViz Object).
         */
        constructor(parent: PCAViz) {
            this._parent = parent;
        }

        /**
         * Loads a JSON file written from the PCAViz python script. Contains
         * all the information needed to visualize the simulation.
         * @param  {string}   path      The JSON path.
         * @param  {Function} callBack  The callback function once loaded.
         * @returns void
         */
        public "loadJSON"(path: string, callBack: any = () => {}): void {
            jjQuery.getJSON(path, (data: any) => {
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

                // Set up the first-frame positions.
                this._loadJSONSetupFirstFramePos(data, precisionFactor);

                // Create a PDB file from the JSON and load it.
                this._makePDBFromJSON(data);

                // Fire the callback.
                callBack();
            }); /*.done(() => {
                callBack();
            }).fail(() => {
                console.log( "error" );
            }); *//* .always(function() {
                console.log( "complete" );
            }); */
        };


        /**
         * Sets up the frame coefficients, including interpolating between them.
         * @param  {Object<string,*>}  data             The data from the JSON.
         * @param  {number}            precisionFactor  The precision to use
         *                                              (number of decimal
         *                                              points).
         */
        private _loadJSONSetupFrameCoeffs(data: any, precisionFactor: number) {
            // Setup the frames
            this._parent._numFramesInJSON = data["coeffs"].length;
            this._parent._frameData = {};  // Keys will be frame indexes,
                                            // values will be Float32Array
                                            // containing the PCA coefficients
                                            // for the corresponding frame.
            this._parent._frameStride = data["params"]["stride"];
            this._parent._numFramesTotal = this._parent._numFramesInJSON * this._parent._frameStride;

            // Convert frames to array of typed arrays. It's faster.
            let len = data["coeffs"].length;
            for (let idx1 = 0; idx1 < len; idx1++) {
                this._parent._frameData[idx1 * this._parent._frameStride] = new Float32Array(
                    data["coeffs"][idx1]
                );
            }

            // The coefficients need to be converted from int to floats.
            len = this._parent._frameData[0].length;
            let frameDataKeys = Object.keys(this._parent._frameData);
            let len2 = frameDataKeys.length;
            for (let i2 = 0; i2 < len2; i2++) {
                let idx1 = +frameDataKeys[i2];
                for (let i = 0; i < len; i++) {
                    this._parent._frameData[idx1][i] = precisionFactor * this._parent._frameData[idx1][i];
                }
            }

            // Now go through and fill in the ones inbetween the explicitly
            // specified frames (with interpolation). Doing linear interpolation
            // for simplicity (rather than spline, for example).
            let numExistingFrameData = Object.keys(this._parent._frameData).length;
            let numFrames = numExistingFrameData * this._parent._frameStride;
            let beforeFrameIdx = undefined;
            let afterFrameIdx = undefined;
            let deltaDefinedFrames = undefined;
            for (let frameIdx = 0; frameIdx < numFrames; frameIdx++) {
                if (this._parent._frameData[frameIdx] !== undefined) {
                    // This frame is defined.

                    // Save this as the previous defined frame.
                    beforeFrameIdx = frameIdx;

                    // Save the next defined frame.
                    afterFrameIdx = frameIdx + this._parent._frameStride;

                    // Get the difference between the two frame vectors.
                    if (this._parent._frameData[afterFrameIdx] !== undefined) {
                        // So if it isn't defined, keep previously caluclated
                        // deltas. Typically affects the end of the sim.
                        deltaDefinedFrames = MathUtils.sumArrayOfFloat32Arrays([
                            this._parent._frameData[afterFrameIdx],
                            MathUtils.multiplyFloat32ArrayByScalar(
                                -1.0, this._parent._frameData[beforeFrameIdx]
                            )
                        ]);
                    }
                } else {
                    // Get the ratio to interpolate (linear).
                    let ratio = (frameIdx - beforeFrameIdx) / this._parent._frameStride;

                    // Calculate the coefficients of the undefined frame.
                    let newFrame = MathUtils.sumArrayOfFloat32Arrays([
                        this._parent._frameData[beforeFrameIdx],
                        MathUtils.multiplyFloat32ArrayByScalar(ratio, deltaDefinedFrames)
                    ]);

                    // Save it to the object.
                    this._parent._frameData[frameIdx] = newFrame;
                }
            }
        }

        /**
         * Sets up the PCA vectors (components).
         * @param  {Object<string,*>}  data             The data from the JSON.
         * @param  {number}            precisionFactor  The precision to use
         *                                              (number of decimal
         *                                              points).
         */
        private _loadJSONSetupPCAComponents(data: any, precisionFactor: number) {
            // Set up the components
            this._parent._componentSize = data["vecs"][0].length;
            this._parent._componentData = []

            // Same with vectors.
            let len = data["vecs"].length;
            for (let idx1 = 0; idx1 < len; idx1++) {
                this._parent._componentData[idx1] = new Float32Array(data["vecs"][idx1]);
            }

            // The vectors need to be converted from int to floats.
            len = this._parent._componentData.length
            let len2 = this._parent._componentData[0].length;
            for (let idx1 = 0; idx1 < len; idx1++) {
                for (let i = 0; i < len2; i++) {
                    this._parent._componentData[idx1][i] = precisionFactor * this._parent._componentData[idx1][i];
                }
            }
        }

        /**
         * Sets up the average-position coordinates.
         * @param  {Object<string,*>}  data             The data from the JSON.
         * @param  {number}            precisionFactor  The precision to use
         *                                              (number of decimal
         *                                              points).
         */
        private _loadJSONSetupAveragePos(data: any, precisionFactor: number) {
            // Make the average coordinates, converting them from int to float.
            let len = data["coors"].length;
            this._parent._averagePositions = new Float32Array(len);
            for (let i = 0; i < len; i++) {
                this._parent._averagePositions[i] = precisionFactor * data["coors"][i];
            }
        }

        /**
         * Sets up the first-frame coordinates. This is necessary to get the bond
         * orders correct.
         * @param  {Object<string,*>}  data             The data from the JSON.
         * @param  {number}            precisionFactor  The precision to use
         *                                              (number of decimal
         *                                              points).
         */
        private _loadJSONSetupFirstFramePos(data: any, precisionFactor: number) {
            // Make the first-frame coordinates, converting them from int to float.
            let len = data["first_coors"].length;
            this._parent._firstFramePositions = new Float32Array(len);
            for (let i = 0; i < len; i++) {
                this._parent._firstFramePositions[i] = precisionFactor * data["first_coors"][i];
            }
        }

        /**
         * Makes a PDB file of the first-frame coordinates from the JSON file.
         * @param  {Object<string,*>}  data             The data from the JSON.
         * @returns void
         */
        private _makePDBFromJSON(data: any): void {
            let res_info = data["res_info"];
            this._parent._res_info = res_info;
            let curResID = "0";
            let curResName = "";
            let pdbTxt = "";

            // Reshape the averaged coordinates into a list of Float32Array
            // triplets.
            let arange = MathUtils._range(0, this._parent._componentSize / 3);
            let len = arange.length;
            let firstFrameCoors = Array(len);
            for (let i2 = 0; i2 < len; i2++) {
                let i = arange[i2];
                let three_i = 3 * i;
                firstFrameCoors[i2] = new Float32Array([
                    this._parent._firstFramePositions[three_i],
                    this._parent._firstFramePositions[three_i + 1],
                    this._parent._firstFramePositions[three_i + 2]
                ]);
            }

            len = res_info.length;
            for (let idx = 0; idx < len; idx++) {
                let v = res_info[idx];
                if (typeof(v) !== "string") {
                    curResID = v[0].toString();
                    curResName = v[1];
                    v = v[2];
                }

                let coor = firstFrameCoors[idx];

                pdbTxt += this._makePDBLine(
                    idx, v, curResName, "X", curResID, coor[0], coor[1], coor[2]
                );
            }

            this._parent["viewer"]["addPDBTxt"](pdbTxt);
        }

        /**
         * Makes a PDB line.
         * @param  {number} idx       The atom index.
         * @param  {string} atomName  The atom name.
         * @param  {string} resName   The residue name.
         * @param  {string} chain     The chain id.
         * @param  {string} resID     The residue id.
         * @param  {number} x         The X coordinate.
         * @param  {number} y         The Y coordinate.
         * @param  {number} z         The Z coordinate.
         * @returns string The PDB line.
         */
        private _makePDBLine(idx: number, atomName: string, resName: string,
                             chain: string, resID: string, x: number, y: number,
                             z: number): string {
            let element = atomName.substr(0, 2).toUpperCase();
            if (["CL", "BR", "ZN", "MG", "SE", "FE",
                 "AL", "MN", "CO", "NI", "CU"].indexOf(element) === -1) {
                element = " " + element.substr(0, 1);
            }

            return "ATOM  " + this._rjust(5, idx.toString()) +
                   this._rjust(5, atomName) + this._rjust(4, resName) +
                   " " + chain + this._rjust(4, resID) + "    " +
                   this._formatNumForPDB(x) +
                   this._formatNumForPDB(y) +
                   this._formatNumForPDB(z) +
                   "  1.00  0.00          " + element + "  \n";
        }

        /**
         * Makes a multi-frame PDB file of the simulation. Good for debugging.
         * @returns string The text of the multi-frame PDB file.
         */
        public "makePDB"(): string {
            let numAtoms = this._parent._componentSize / 3;
            let pdbTxt = "";
            for (let frameIdx = 0; frameIdx < this._parent._numFramesTotal; frameIdx++) {
                pdbTxt += "MODEL " + frameIdx.toString() + "\n";
                let coors = this._parent.getFrameCoors(frameIdx);
                let curResID = "";
                let curResName = "";
                for (let atomIdx = 0; atomIdx < numAtoms; atomIdx++) {
                    let v = this._parent._res_info[atomIdx];
                    if (typeof(v) !== "string") {
                        curResID = v[0].toString();
                        curResName = v[1];
                        v = v[2];
                    }

                    let coor = coors[atomIdx];
                    pdbTxt += this._makePDBLine(
                        atomIdx, v, curResName, "X", curResID, coor[0], coor[1], coor[2]
                    );
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
         * @param  {*} parent The parent class (PCAViz Object).
         */
        constructor(parent: any) {
            this._parent = parent;

            // NGLViewer specifically fulfills a promise when it has finished
            // loading a molecule, rather than returning the model
            // immediately. We need fake NGLViewer functions that can get
            // called until the actual _model is loaded (to avoid errors).
            this["_model"] = {
                "rebuildRepresentations": function() { return; },
                "structure":{
                    "eachAtom": function(a) { return; }
                }
            }

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
            // Note that by default it strips hydrogens! Unfortunate.
            this["_model"] = this._parent._params["viewer"]["addModel"](
                pdbTxt, "pdb", {"keepH": true}
            );
            this._render();
        }

        /**
         * Updates the atom positions in a 3Dmol.js viewer.
         * @param  {number} frame The frame number.
         * @returns void
         */
        private _3DMolJS_UpdateAtomPos(frame: number): void {
            let newAtomCoors = this._parent.getFrameCoors(frame);

            let atoms = this["_model"]["selectedAtoms"]({});
            let arrLen = atoms.length;
            for (let atomIdx=0; atomIdx<arrLen; atomIdx++) {
                atoms[atomIdx]["x"] = newAtomCoors[atomIdx][0];
                atoms[atomIdx]["y"] = newAtomCoors[atomIdx][1];
                atoms[atomIdx]["z"] = newAtomCoors[atomIdx][2];
            }
        }

        /**
         * Renders the 3DMol.js viewer.
         * @returns void
         */
        private _3DMolJS_Render(): void {
            // Must update styles to actually have atoms move. Annoying.
            this["_model"]["setStyle"]({}, this._parent._params["visStyle"]);
            this._parent._params["viewer"]["render"]();
        }

        /**
         * The function to add PDB txt to the NGLViewer.
         * @param  {string} pdbTxt The PDB text.
         * @returns void
         */
        private _NGL_AddPDBTxt(pdbTxt: string): void {
            this._parent._params["viewer"]["loadFile"](
                new Blob(
                    [pdbTxt],
                    {"type": 'text/plain'}
                ),
                {"ext":"pdb", "defaultRepresentation": true}
            ).then((e) => {
                // console.log(this);
                this["_model"] = e;
            });
        }

        /**
         * Updates the atom positions in a NGLViewer.
         * @param  {number} frame The frame number.
         * @returns void
         */
        private _NGL_UpdateAtomPos(frame: number) {
            // Unfortunately, this check is necessary because NGLViewer
            // returns a promise. If this fires before the promise is
            // fullfilled, you'll get an error.
            // if (this["_model"] !== undefined) {
                let newAtomCoors = this._parent.getFrameCoors(frame);

                let atomIdx = 0;
                this["_model"]["structure"]["eachAtom"]((atom, idx) => {
                    atom["x"] = newAtomCoors[atomIdx][0];
                    atom["y"] = newAtomCoors[atomIdx][1];
                    atom["z"] = newAtomCoors[atomIdx][2];
                    atomIdx++;
                });
            // }
        }

        /**
         * Renders the NGLViewer.
         * @returns void
         */
        private _NGL_Render() {
            // Unfortunately, this check is necessary because NGLViewer
            // returns a promise. If this fires before the promise is
            // fullfilled, you'll get an error.
            // if (this["_model"] !== undefined) {
                this["_model"]["rebuildRepresentations"]();
            // }
        }

        /**
         * The function to add PDB txt to the PV Viewer.
         * @param  {string} pdbTxt The PDB text.
         * @returns void
         */
        private _PV_AddPDBTxt(pdbTxt: string): void {
            this["_model"] = this._parent._params["viewer"]["library"]["io"]["pdb"](pdbTxt);
        }

        /**
         * Updates the atom positions in a PV Viewer.
         * @param  {number} frame The frame number.
         * @returns void
         */
        private _PV_UpdateAtomPos(frame: number) {
            let newAtomCoors = this._parent.getFrameCoors(frame);

            let atomIdx = 0;
            this["_model"]["eachAtom"]((atom, idx) => {
                atom["_bV"] = newAtomCoors[atomIdx];
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
            this._parent._params["visStyle"](this["_model"]);
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
                this["_model"] = candidateModel;
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
                newAtomCoors, this["_model"], this._parent._params["viewer"], this
            );
        }

        /**
         * Renders if in GENERIC mode.
         * @returns void
         */
        private _GENERIC_Render() {
            this._parent._params["render"](
                this["_model"], this._parent._params["viewer"], this
            );
        }
    }

    class _Player {
        private _animationFrameID = undefined;
        private _parent = undefined;

        /**
         * The constructor of a _Player class.
         * @param  {*} parent The parent class (PCAViz Object).
         */
        constructor(parent: any) {
            this._parent = parent;

            // You should automatically stop if the tab looses focus. Doing this
            // to prevent overheating laptops.
            // if (!runningUnderNodeJS) {
                // check to make sure in browser.
                window.onblur = () => {
                    this["stop"]();
                };
            // }
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

            let timestampLastFire: number = 0;
            let deltaStartRatio: number = this._parent._playerControls.getSliderValue();
            let loopTimestampOnStart = undefined;
            let msSinceLastFire: number = undefined;
            let playRatio: number = undefined;
            let curFrame: number = undefined;

            /**
             * The loop function. Note that this needs to be VERY optimized. Keep
             * if statements and such to a minimum.
             * @param  {number} timestamp The number of milliseconds since the
             *                            requestAnimationFrame started.
             */
            let loop = (timestamp: number) => {
                loopTimestampOnStart = (loopTimestampOnStart === undefined) ? timestamp : loopTimestampOnStart;

                // Redefine timestamp so it starts from 0 every time play is pressed.
                timestamp = timestamp - loopTimestampOnStart;

                msSinceLastFire = timestamp - timestampLastFire;  //  - loopTimestampOnStart;

                // Enough time has passed that we should update.
                if (msSinceLastFire > this._parent._params["updateFreqInMilliseconds"]) {
                    // Save the new timestampLastFire value.
                    timestampLastFire = timestamp;

                    // How far along the animation are you?
                    playRatio = (deltaStartRatio + timestamp / this._parent._params["durationInMilliseconds"]);

                    // If you've gone over the end of the animation and it's not
                    // set to loop, stop.
                    if (playRatio > 1.0) {
                        // You've gone over the end of the animation.
                        if (!this._parent._params["loop"]) {
                            // No loop requested, so stop.
                            this["stop"]();  // Sets this._animationFrameID to undefined.
                            playRatio = 1.0;
                            // this._parent._playerControls.setSlider(this._parent._numFramesTotal);
                            // return;
                        } else {
                            // Loop requested.
                            playRatio = playRatio % 1.0;
                        }
                    }

                    // Get the current frame.
                    curFrame = Math.floor(this._parent._numFramesTotal * playRatio);

                    this._parent["viewer"].updateAtomPos(curFrame);

                    // Also update UI.
                    this._parent._playerControls.setSlider(curFrame);
                }

                // Start the next iteration.
                if (this._animationFrameID !== undefined) {
                    requestAnimationFrame(loop);
                }
            }

            this._animationFrameID = requestAnimationFrame(loop);

            // Change ui if it is active.
            this._parent._playerControls.onPlayStart();
        }

        /**
         * Stops the simulation animation.
         * @returns void
         */
        public "stop"(): void {
            // Clear previous requestAnimationFrames.
            if (this._animationFrameID !== undefined) {
                cancelAnimationFrame(this._animationFrameID);
                this._animationFrameID = undefined;
            }

            this._parent._playerControls.onPlayStop();
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
            let len = arrayOfFloat32Arrays.length;
            if (len === 1) {
                // Just one item in the array, so return that first item.
                return arrayOfFloat32Arrays[0];
            } else {
                // Multiple items. So need to sum them.
                let len2 = arrayOfFloat32Arrays[0].length;
                let summed = new Float32Array(len2);
                for (let arrIdx = 0; arrIdx < len; arrIdx++) {
                    let arr = arrayOfFloat32Arrays[arrIdx];
                    for (let i = 0; i < len2; i++) {
                        let v = arr[i];
                        summed[i] = v + summed[i];
                    }
                }

                return summed;
            }
        }

        /**
         * Multiplies a Float32Array by a scalar.
         * @param  {number}       scalar        The scalar.
         * @param  {Float32Array} float32Array  The Float32Array.
         * @returns Float32Array  The multiplied Float32Array.
         */
        export function multiplyFloat32ArrayByScalar(scalar: number, float32Array: Float32Array): Float32Array {
            let len = float32Array.length;
            switch (scalar) {
                case 0.0:
                    return new Float32Array(len);
                case 1.0:
                    return float32Array;
                default:
                    let newFloat32Array = new Float32Array(len);
                    for (let i = 0; i < len; i++) {
                        newFloat32Array[i] = scalar * float32Array[i];
                    }
                    return newFloat32Array;
            }
        }
    }

    class _PlayerControls {
        private _parent = undefined;
        private _playDOM = undefined;
        private _pauseDOM = undefined;
        private _sliderDOM = undefined;
        private _uiActive = false;

        // SVGs of the play buttons.
        private _images = {
            play: `data:image/svg+xml;charset=utf-8;base64,PHN2ZyB2ZXJzaW9uPSIxLjIiIGJhc2VQcm9maWxlPSJ0aW55IiBpZD0iTGF5ZXJfMSIgeG1sbnM9Imh0dHA6Ly93d3cudzMub3JnLzIwMDAvc3ZnIiB4bWxuczp4bGluaz0iaHR0cDovL3d3dy53My5vcmcvMTk5OS94bGluayINCng9IjBweCIgeT0iMHB4IiB3aWR0aD0iNDBweCIgaGVpZ2h0PSI0MHB4IiB2aWV3Qm94PSIwIDAgNDAgNDAiIHhtbDpzcGFjZT0icHJlc2VydmUiPg0KPHBhdGggZmlsbD0iIzAxMDEwMSIgc3Ryb2tlPSIjMDEwMTAxIiBzdHJva2Utd2lkdGg9IjYuNjQ0OSIgc3Ryb2tlLWxpbmVqb2luPSJyb3VuZCIgZD0iTTcuMSwzNi40TDMyLjksMjBMNy4xLDMuNlYzNi40eiIvPg0KPC9zdmc+`,
            pause: `data:image/svg+xml;charset=utf-8;base64,PHN2ZyB2ZXJzaW9uPSIxLjIiIGJhc2VQcm9maWxlPSJ0aW55IiBpZD0iTGF5ZXJfMSIgeG1sbnM9Imh0dHA6Ly93d3cudzMub3JnLzIwMDAvc3ZnIiB4bWxuczp4bGluaz0iaHR0cDovL3d3dy53My5vcmcvMTk5OS94bGluayINCng9IjBweCIgeT0iMHB4IiB3aWR0aD0iNDBweCIgaGVpZ2h0PSI0MHB4IiB2aWV3Qm94PSIwIDAgNDAgNDAiIHhtbDpzcGFjZT0icHJlc2VydmUiPg0KPGc+DQo8cGF0aCBmaWxsPSIjMDEwMTAxIiBzdHJva2U9IiMwMTAxMDEiIHN0cm9rZS13aWR0aD0iNi40MTM0IiBzdHJva2UtbGluZWpvaW49InJvdW5kIiBkPSJNMy43LDMuNUgxNWMwLjEsMCwwLjEsMC4xLDAuMSwwLjN2MzIuMw0KYzAsMC4yLTAuMSwwLjMtMC4xLDAuM0gzLjdjLTAuMSwwLTAuMS0wLjEtMC4xLTAuM1YzLjlDMy41LDMuNywzLjYsMy41LDMuNywzLjUiLz4NCjxwYXRoIGZpbGw9IiMwMTAxMDEiIHN0cm9rZT0iIzAxMDEwMSIgc3Ryb2tlLXdpZHRoPSI2LjQxMzQiIHN0cm9rZS1saW5lam9pbj0icm91bmQiIGQ9Ik0yNSwzLjVoMTEuNGMwLjEsMCwwLjEsMC4xLDAuMSwwLjN2MzIuMw0KYzAsMC4yLTAuMSwwLjMtMC4xLDAuM0gyNWMtMC4xLDAtMC4xLTAuMS0wLjEtMC4zVjMuOUMyNC44LDMuNywyNC45LDMuNSwyNSwzLjUiLz4NCjwvZz4NCjwvc3ZnPg==`
        }

        /**
         * The constructor of an _IO class.
         * @param  {*} parent The parent class (PCAViz Object).
         */
        constructor(parent: PCAViz) {
            this._parent = parent;

            if (this._parent._params["playerControlsID"] !== "") {
                // Record that the UI is active.
                this._uiActive = true;

                // Get the control div.
                let controllerDiv = document.getElementById(
                    this._parent._params["playerControlsID"]
                );

                // Add the content to that div.
                let randomID = "id-" + Math.random().toString().replace(/\./g, "");
                controllerDiv.innerHTML = `
                    <div class="pcaviz-controls-container" style="display: flex; flex-direction: row;">
                        <div class="pcaviz-button-container" style="flex: none;">
                            <button id="${randomID}-play" class="pcaviz-play pcaviz-button" style="height: 100%;">
                                <img src="${this._images.play}" />
                            </button>
                            <button id="${randomID}-pause" class="pcaviz-pause pcaviz-button" style="height: 100%; display: none;">
                                <img src="${this._images.pause}" />
                            </button>
                        </div>
                        <div class="pcaviz-slider-container" style="flex: auto; display: inline-flex;">
                            <input id="${randomID}-slider" type="range" min="0" max="100" value="0" style="width: 100%;">
                        </div>
                    </div>
                `;

                // Keep track of the components added.
                this._playDOM = document.getElementById(randomID + "-play");
                this._pauseDOM = document.getElementById(randomID + "-pause");
                this._sliderDOM = document.getElementById(randomID + "-slider");

                // Make these elements clickable.
                this._playDOM.addEventListener("click", () => {
                    // console.log(this._parent._params);
                    this._parent["player"]["start"](this._parent._params);
                });
                this._pauseDOM.addEventListener("click", () => {
                    this._parent["player"]["stop"]();
                });
                this._sliderDOM.onchange = () => {
                    let ratio = this._sliderDOM.value / 100.0;
                    this._parent["player"]["toFrame"](
                        Math.round(this._parent._numFramesTotal * ratio)
                    );
                }
            }
        }

        /**
         * Fires when play starts. To update buttons.
         * @returns void
         */
        public onPlayStart(): void {
            if (this._uiActive === true) {
                // Make pause button visible.
                this._playDOM.style.display = "none";
                this._pauseDOM.style.display = "block";
            }
        }

        /**
         * Fires when play stops. To update buttons,
         * @returns void
         */
        public onPlayStop(): void {
            if (this._uiActive === true) {
                // Make pause button visible.
                this._playDOM.style.display = "block";
                this._pauseDOM.style.display = "none";
            }
        }

        /**
         * Updates the UI slider.
         * @param  {number} frame  The frame.
         * @returns void
         */
        public setSlider(frame: number): void {
            if (this._uiActive === true) {
                this._sliderDOM.value = Math.round(
                    100 * frame / this._parent._numFramesTotal
                );
            }
        }

        /**
         * Get's the current slider value, as a ratio.
         * @returns number  The ratio (0.0 to 1.0).
         */
        public getSliderValue(): number {
            if (this._uiActive === true) {
                return this._sliderDOM.value / 100.0;
            }
            return 0.0;
        }
    }

    // I don't want to require jQuery, but I need some of its functions.
    namespace jjQuery {

        /**
         * Extends an object with two additional objects. Inspired by
         * jQuery.extend.
         * @param  {Object<string, *>} a  The first object.
         * @param  {Object<string, *>} b  The second object.
         * @param  {Object<string, *>} c  The third object.
         * @returns The extended object.
         */
        export function extend(a: any, b: any, c: any): any {
            for (let k in b) {
                if (b.hasOwnProperty(k)) { a[k] = b[k]; }
            }
            for (let k in c) {
                if (c.hasOwnProperty(k)) { a[k] = c[k]; }
            }
            return a;
        }

        /**
         * Gets JSON from a remote file. Inspired by jQuery.getJSON.
         * @param  {string}       path          The remote URL.
         * @param  {function(*)}  callBackFunc  A callback function once done.
         * @returns void
         */
        export function getJSON(path: string, callBackFunc: any): void {
            var xhttp = new XMLHttpRequest();
            xhttp.onreadystatechange = function() {  // Shouldn't be arrow func.
                if (this.readyState == 4 && this.status == 200) {
                    callBackFunc(
                        JSON.parse(
                            this.responseText
                        )
                    );
                }
            };
            xhttp.open("GET", path, true);
            xhttp.send();
        }
    }
}

// Leave the below here in case you want to use from nodejs in the future...

// let runningUnderNodeJS = false;
// try {
//     window;
// } catch (err) {
//     // You must be running this using nodejs.
//     runningUnderNodeJS = true;
// }

// if (!runningUnderNodeJS) {
    // Polyfill requestAnimationFrame
    // http://paulirish.com/2011/requestanimationframe-for-smart-animating/
    // http://my.opera.com/emoller/blog/2011/12/20/requestanimationframe-for-smart-er-animating
    // requestAnimationFrame polyfill by Erik Möller. fixes from Paul Irish and Tino Zijdel
    // MIT license
    (function() {
        var lastTime = 0;
        var vendors = ['ms', 'moz', 'webkit', 'o'];
        for (var x = 0; x < vendors.length && !window["requestAnimationFrame"]; ++x) {
            window["requestAnimationFrame"] = window[vendors[x]+'RequestAnimationFrame'];
            window["cancelAnimationFrame"] = window[vendors[x]+'CancelAnimationFrame']
                                             || window[vendors[x]+'CancelRequestAnimationFrame'];
        }

        if (!window["requestAnimationFrame"])
            window["requestAnimationFrame"] = function(callback) {
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

    // To survive closure compiler. See fix_namespaces.py to see where
    // window["PCAVizNameSpace"] comes from.
    (<any>window)["PCAViz"] = window["PCAVizNameSpace"].PCAViz;
    // (<any>window)["PCAViz"] = PCAVizNameSpace.PCAViz;

// } else {
//     // It's running under node js. Note that you won't be using the closure
//     // compiled version of the code from nodejs, so no need to worry about
//     // that here.

//     // Get the json file from the command line and make the trajectory.
//     let jsonFile = process.argv[2];

//     // Make a fake jQuery.getJSON function. Overwriting the existing one to
//     // work with node.js.
//     let jsonData;
//     jjQuery.getJSON = function(path: string, callBackFunc: any) {
//         let fs = require("fs");

//         let content = fs.readFileSync(jsonFile);
//         let jsonData = JSON.parse(content);
//         callBackFunc(jsonData);

//         return {
//             done: function (func) {
//                 func();
//                 return {
//                     fail: function (func) {
//                         // Never allow fail.
//                         // func();
//                         return;
//                     }
//                 }
//             },
//         }
//     }

//     // Now create the
//     let pcaViz = new PCAViz({
//         viewer: {},
//         viewerType: 'GENERIC',
//         windowAverageSize: 1,
//         loadPDBTxt: (pdbTxt, viewer, pcaViz) => { return; },  // viewer.addModel(pdbTxt, "pdb"),
//         updateAtomPositions: (newAtomCoors, model, viewer, pcaViz) => { return ; },
//         render: (model, viewer, pcaViz) => {
//             return;
//         },
//     });

//     pcaViz["io"].loadJSON("data.json", () => {
//         let func = pcaViz["io"].makePDB.bind(pcaViz["io"]);
//         let pdbTxt = func(jsonData);
//         console.log(pdbTxt);
//     });
// }
