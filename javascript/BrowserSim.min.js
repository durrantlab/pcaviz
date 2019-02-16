var BrowserSim = /** @class */ (function () {
    /**
     * Constructor for the BrowserSim class.
     * @param  {Object<string,*>} params
     * @returns void
     */
    function BrowserSim(params) {
        // Note that the below are public so they can be accessed from other
        // classes. But they are not public to the user and so do not need to be
        // protected from closure compiler.
        this._numFrames = 0;
        this._frameSize = 0;
        // Keys are frame indecies. Values are Float32Arrays with the
        // coefficients.
        this._frameData = {};
        this._frameStride = 1;
        this._numComponents = undefined;
        this._componentSize = undefined;
        this._componentData = [];
        this._averagePositions = undefined;
        this._params = {};
        this["io"] = new _IO(this); // This way so it survives closure compiler.
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
    BrowserSim.prototype.updateParams = function (updatedParams) {
        var _this = this;
        // Set some common error text blurbs.
        var genericModeBut = 'You are running BrowserSim in GENERIC mode, but ' +
            'you haven\'t specified ';
        var modelVar = 'The "model" variable is whatever your loadPDBTxt(...) ' +
            'function returns.';
        var viewerVar = 'The "viewer" variable is the user-specified BrowserSim ' +
            '"viewer" parameter.';
        var browserSimVar = 'The "browserSim" variable is the instantiated ' +
            'BrowserSim object.';
        // Set the defaults and keeps any previous ones not now specified.
        var defaults = {
            "viewer": undefined,
            "viewerType": undefined,
            "durationInMilliseconds": 10000,
            "updateFreqInMilliseconds": 10,
            "loop": true,
            "windowAverageSize": 1,
            "parent": this,
            "loadPDBTxt": function (pdbTxt, viewer, browserSim) {
                throw new Error(genericModeBut +
                    'a loadPDBTxt(pdbTxt, viewer, browserSim) function. This ' +
                    'function loads PDB text into the viewer and optionally ' +
                    'returns a model object. For your reference, the current ' +
                    'value of the "pdbTxt" variable (a string) starts with:\n\n' +
                    pdbTxt.toString().slice(0, 500)) + "\n\n" + viewerVar + "\n\n" +
                    browserSimVar + "\n\n";
            },
            "updateAtomPositions": function (newAtomCoors, model, viewer, browserSim) {
                var goodCoorRep = "[\n";
                for (var idx in newAtomCoors) {
                    if (newAtomCoors.hasOwnProperty(idx)) {
                        var idxVal = parseInt(idx, 10);
                        var coor = newAtomCoors[idx];
                        goodCoorRep += "  Float32Array [" + coor[0] + ", " + coor[1] + ", " + coor[2] + "],\n";
                        if (idxVal > 10) {
                            break;
                        }
                    }
                }
                goodCoorRep += "  ...\n]";
                throw new Error(genericModeBut + 'an updateAtomPositions(newAtomCoors, ' +
                    'model, viewer, browserSim) function. This function provides ' +
                    'the updated atom coordinates for the current frame. The ' +
                    'value of the "newAtomCoors" variable (a list of Float32Array ' +
                    'containing the new coordinates of the atoms) ' +
                    'looks like:\n\n' + goodCoorRep + '\n\n' + modelVar +
                    "\n\n" + viewerVar + "\n\n" + browserSimVar + "\n\n");
            },
            "render": function (model, viewer, browserSim) {
                throw new Error(genericModeBut + 'a render(model, viewer, browserSim) function. ' +
                    'This function runs every time the atom coordinates ' +
                    'change, to update what the viewer renders. ' + modelVar +
                    ' ' + viewerVar + ' ' + browserSimVar + "\n\n");
            }
        };
        this._params = jQuery.extend(defaults, this._params, updatedParams);
        // halfWindowSize is always derived from windowAverageSize
        this._params["halfWindowSize"] = Math.floor((this._params["windowAverageSize"] - 1) / 2);
        // Default visStyle depends on viewerType
        if (this._params["visStyle"] === undefined) {
            switch (this._params["viewerType"]) {
                case "3DMOLJS":
                    this._params["visStyle"] = { "line": {} };
                    break;
                case "NGL":
                    this._params["visStyle"] = undefined;
                    break;
                case "PV":
                    this._params["visStyle"] = function (structure) { _this._params["viewer"]["viewer"]["spheres"]('all', structure); };
                    break;
                case "GENERIC":
                    this._params["visStyle"] = undefined;
                    break;
            }
        }
    };
    /**
     * Get the frame coordinates.
     * @param  {number} frame The frame number.
     * @returns Float32Array  A list of Float32Array containing the atom
     *                        coordinates, for each frame.
     */
    BrowserSim.prototype.getFrameCoors = function (frame) {
        var _this = this;
        // Consider multiple frames if necessary.
        var framesToAvg = MathUtils._range(frame - this._params["halfWindowSize"], frame + this._params["halfWindowSize"] + 1);
        // Get the coefficients for each of those frames (an array of
        // Float32Array's)
        var framesCoefficients = framesToAvg.map(function (frameIdx) {
            // Make sure frame never out of bounds.
            if (frameIdx < 0) {
                frameIdx += _this._numFrames;
            }
            frameIdx = frameIdx % _this._numFrames;
            // Get the frame data (PCA coefficients).
            return _this._frameData[frameIdx];
        });
        // Get the flattened coordinates for the frames (i.e., coordinates not
        // separated into triplets).
        var coorsFlattenedForFrames = framesCoefficients.map(function (frameCoefficients) {
            // frameCoefficients is a Float32Array containing all the
            // coefficients for a given frame. Need to multiple those by the
            // corresponding components.
            var multipledComponents = _this._componentData.map(function (component, componentIdx) {
                return MathUtils.multiplyFloat32ArrayByScalar(frameCoefficients[componentIdx], component);
            });
            var summedMultipledComponents = MathUtils.sumArrayOfFloat32Arrays(multipledComponents);
            return summedMultipledComponents;
        });
        // Now average the flattened coordinates over the frames.
        var summedCoorsOverFrames = coorsFlattenedForFrames.reduce(function (summedCoors, newCoors) {
            return MathUtils.sumArrayOfFloat32Arrays([summedCoors, newCoors]);
        });
        var averageCoorsOverFrames = MathUtils.multiplyFloat32ArrayByScalar(1.0 / framesToAvg.length, summedCoorsOverFrames);
        // Reshape the averaged coordinates into a list of Float32Array triplets.
        var coors = MathUtils._range(0, this._componentSize / 3).map(function (i) {
            return new Float32Array([
                averageCoorsOverFrames[i],
                averageCoorsOverFrames[i + 1],
                averageCoorsOverFrames[i + 2]
            ]);
        });
        return coors;
    };
    return BrowserSim;
}());
var _IO = /** @class */ (function () {
    /**
     * The constructor of an _IO class.
     * @param  {*} parent The parent class (BrowserSim Object).
     */
    function _IO(parent) {
        this._parent = undefined;
        this._parent = parent;
    }
    /**
     * Loads a JSON file written from the BrowserSim python script. Contains
     * all the information needed to visualize the simulation.
     * @param  {string}   path      The JSON path.
     * @param  {Function} callBack  The callback function once loaded.
     * @returns void
     */
    _IO.prototype["loadJSON"] = function (path, callBack) {
        var _this = this;
        if (callBack === void 0) { callBack = function () { }; }
        jQuery.getJSON(path, function (data) {
            // Setup the frames
            _this._parent._numFrames = data["coeffs"].length;
            _this._parent._frameSize = data["coeffs"][0].length;
            _this._parent._frameData = {};
            _this._parent._frameStride = 5; // TODO: SHOULD BE USER DEFINED.
            // Set up the components
            _this._parent._numComponents = data["vecs"].length;
            _this._parent._componentSize = data["vecs"][0].length;
            _this._parent._componentData = [];
            // Here you will put the average positions, but not ready yet.
            _this._parent._averagePositions = new Float32Array(_this._parent._numComponents * _this._parent._componentSize);
            // The length of the frame coefficients must equal the number of
            // components.
            if (_this._parent._frameSize !== _this._parent._numComponents) {
                throw new Error("The number of coefficients per frame " +
                    ("(" + _this._parent._frameSize + ") is not the same") +
                    ("as the number of components (" + _this._parent._numComponents + ")"));
            }
            // Convert frames to array of typed arrays. It's faster.
            var precision = data["params"]["precision"];
            for (var idx1 in data["coeffs"]) {
                if (data["coeffs"].hasOwnProperty(idx1)) {
                    var idx1Num = parseInt(idx1, 10);
                    _this._parent._frameData[idx1Num * _this._parent._frameStride] = new Float32Array(data["coeffs"][idx1Num].map(function (v) { return MathUtils.convertToDecimal(v, precision); }));
                }
            }
            // Now go through and fill in the ones inbetween the explicitly
            // specified frames (with interpolation). Doing linear
            // interpolation for simplicity (rather than spline, for example).
            for (var frameIdx = 0; frameIdx < data["coeffs"].length * _this._parent._frameStride; frameIdx++) {
                if (_this._parent._frameData[frameIdx] === undefined) {
                    // This frame isn't defined.
                    // Get the previous defined frame.
                    var beforeFrame = _this._parent._frameStride * Math.floor(frameIdx / _this._parent._frameStride);
                    // Get the next defined frame.
                    var afterFrame = beforeFrame + _this._parent._frameStride;
                    // If either doesn't exist, go to the next one.
                    if ((_this._parent._frameData[beforeFrame] === undefined) ||
                        (_this._parent._frameData[afterFrame] === undefined)) {
                        continue;
                    }
                    // Get the ratio to interpolate (linear).
                    var ratio = (frameIdx - beforeFrame) / _this._parent._frameStride;
                    // Get the difference between the two vectors.
                    var delta = MathUtils.sumArrayOfFloat32Arrays([
                        _this._parent._frameData[beforeFrame],
                        MathUtils.multiplyFloat32ArrayByScalar(-1.0, _this._parent._frameData[afterFrame])
                    ]);
                    // Calculate the coefficients of the undefined frame.
                    var newFrame = MathUtils.sumArrayOfFloat32Arrays([
                        _this._parent._frameData[beforeFrame],
                        MathUtils.multiplyFloat32ArrayByScalar(ratio, delta)
                    ]);
                    // Save it to the object.
                    _this._parent._frameData[frameIdx] = newFrame;
                }
            }
            // Same with vectors.
            for (var idx1 in data["vecs"]) {
                if (data["vecs"].hasOwnProperty(idx1)) {
                    var idx1Num = parseInt(idx1, 10);
                    _this._parent._componentData[idx1Num] = new Float32Array(data["vecs"][idx1Num].map(function (v) { return MathUtils.convertToDecimal(v, precision); }));
                }
            }
            console.log(data);
            debugger;
            _this.makePDBFromJSON(data);
        }).done(function () {
            callBack();
        }).fail(function () {
            console.log("error");
        }); /* .always(function() {
            console.log( "complete" );
        }); */
    };
    ;
    _IO.prototype.makePDBFromJSON = function (data) {
        console.log(data);
        console.log(this._parent.getFrameCoors(0));
        return "";
    };
    /**
     * Makes a multi-frame PDB file of the simulation. Good for debugging.
     * @returns string The text of the multi-frame PDB file.
     */
    _IO.prototype["makePDB"] = function () {
        var numAtoms = this._parent._componentSize / 3;
        var pdbTxt = "";
        for (var frameIdx = 0; frameIdx < this._parent._numFrames; frameIdx++) {
            pdbTxt += "MODEL " + frameIdx.toString() + "\n";
            for (var atomIdx = 0; atomIdx < numAtoms; atomIdx++) {
                var coor = this._parent.getAtomCoors(frameIdx, atomIdx);
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
    ;
    /**
     * Formats a number for use in a PDB line.
     * @param  {number} val The number
     * @returns string  The formatted number (rounded, padded, etc.).
     */
    _IO.prototype._formatNumForPDB = function (val) {
        var valStr = val.toFixed(3);
        return this._rjust(8, valStr);
    };
    ;
    /**
     * Right justifies a string.
     * @param  {number}     length      The length of the string and padding.
     * @param  {string}     origString  The original string.
     * @param  {string=}  padChar     The character to use for padding.
     *                                  Space by default.
     * @returns string      The padded, right-justified string.
     */
    _IO.prototype._rjust = function (length, origString, padChar) {
        if (padChar === void 0) { padChar = " "; }
        while (origString.length < length) {
            origString = padChar + origString;
        }
        return origString;
    };
    ;
    return _IO;
}());
var _Viewer = /** @class */ (function () {
    /**
     * The constructor of a _Viewer class.
     * @param  {*} parent The parent class (BrowserSim Object).
     */
    function _Viewer(parent) {
        this._model = undefined;
        this._updateLocked = false; // Don't update if currently updating.
        this._updateAtomPosFun = undefined;
        this._render = undefined;
        this._parent = undefined;
        this._parent = parent;
        // Perform some checks.
        if (this._parent._params["viewer"] === undefined) {
            throw new Error("No viewer specified!");
        }
        if (this._parent._params["viewerType"] === undefined) {
            throw new Error("No viewer type specified!");
        }
        var validTypes = ["3DMOLJS", "NGL", "PV", "GENERIC"];
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
    _Viewer.prototype.updateAtomPos = function (frame) {
        if (this._updateLocked) {
            return;
        }
        this._updateLocked = true;
        this._updateAtomPosFun(frame);
        this._updateLocked = false;
        this._render();
    };
    /**
     * The function to add PDB txt to the 3Dmol.js viewer.
     * @param  {string} pdbTxt The PDB text.
     * @returns void
     */
    _Viewer.prototype._3DMolJS_AddPDBTxt = function (pdbTxt) {
        this._model = this._parent._params["viewer"].addModel(pdbTxt, "pdb");
        this._render();
    };
    /**
     * Updates the atom positions in a 3Dmol.js viewer.
     * @param  {number} frame The frame number.
     * @returns void
     */
    _Viewer.prototype._3DMolJS_UpdateAtomPos = function (frame) {
        var newAtomCoors = this._parent.getFrameCoors(frame);
        var atoms = this._model.selectedAtoms({});
        for (var atomIdx = 0; atomIdx < atoms.length; atomIdx++) {
            var coors = newAtomCoors[atomIdx];
            atoms[atomIdx]["x"] = coors[0];
            atoms[atomIdx]["y"] = coors[1];
            atoms[atomIdx]["z"] = coors[2];
        }
    };
    /**
     * Renders the 3DMol.js viewer.
     * @returns void
     */
    _Viewer.prototype._3DMolJS_Render = function () {
        // Must update styles to actually have atoms move. Annoying.
        this._model.setStyle({}, this._parent._params["visStyle"]);
        this._parent._params["viewer"].render();
    };
    /**
     * The function to add PDB txt to the NGLViewer.
     * @param  {string} pdbTxt The PDB text.
     * @returns void
     */
    _Viewer.prototype._NGL_AddPDBTxt = function (pdbTxt) {
        var _this = this;
        this._parent._params["viewer"].loadFile(new Blob([pdbTxt], { type: 'text/plain' }), { ext: 'pdb', defaultRepresentation: true }).then(function (e) {
            _this._model = e;
        });
    };
    /**
     * Updates the atom positions in a NGLViewer.
     * @param  {number} frame The frame number.
     * @returns void
     */
    _Viewer.prototype._NGL_UpdateAtomPos = function (frame) {
        var newAtomCoors = this._parent.getFrameCoors(frame);
        var atomIdx = 0;
        this._model["structure"]["eachAtom"](function (atom, idx) {
            var coors = newAtomCoors[atomIdx];
            atom["x"] = coors[0];
            atom["y"] = coors[1];
            atom["z"] = coors[2];
            atomIdx++;
        });
    };
    /**
     * Renders the NGLViewer.
     * @returns void
     */
    _Viewer.prototype._NGL_Render = function () {
        this._model["rebuildRepresentations"]();
    };
    /**
     * The function to add PDB txt to the PV Viewer.
     * @param  {string} pdbTxt The PDB text.
     * @returns void
     */
    _Viewer.prototype._PV_AddPDBTxt = function (pdbTxt) {
        this._model = this._parent._params["viewer"]["library"]["io"]["pdb"](pdbTxt);
    };
    /**
     * Updates the atom positions in a PV Viewer.
     * @param  {number} frame The frame number.
     * @returns void
     */
    _Viewer.prototype._PV_UpdateAtomPos = function (frame) {
        var newAtomCoors = this._parent.getFrameCoors(frame);
        var atomIdx = 0;
        this._model["eachAtom"](function (atom, idx) {
            var coors = newAtomCoors[atomIdx];
            atom["_bV"] = coors;
            atomIdx++;
        });
    };
    /**
     * Renders the PV Viewer.
     * @returns void
     */
    _Viewer.prototype._PV_Render = function () {
        this._parent._params["viewer"]["viewer"].clear();
        // A callback function in the case of PV.
        this._parent._params["visStyle"](this._model);
    };
    /**
     * The function to add PDB txt if in GENERIC mode.
     * @param  {string} pdbTxt The PDB text.
     * @returns void
     */
    _Viewer.prototype._GENERIC_AddPDBTxt = function (pdbTxt) {
        var candidateModel = this._parent._params["loadPDBTxt"](pdbTxt, this._parent._params["viewer"], this);
        if (candidateModel !== undefined) {
            this._model = candidateModel;
        }
    };
    /**
     * Updates the atom positions if in GENERIC mode.
     * @param  {number} frame The frame number.
     * @returns void
     */
    _Viewer.prototype._GENERIC_UpdateAtomPos = function (frame) {
        var newAtomCoors = this._parent.getFrameCoors(frame);
        this._parent._params["updateAtomPositions"](newAtomCoors, this._model, this._parent._params["viewer"], this);
    };
    /**
     * Renders if in GENERIC mode.
     * @returns void
     */
    _Viewer.prototype._GENERIC_Render = function () {
        this._parent._params["render"](this._model, this._parent._params["viewer"], this);
    };
    return _Viewer;
}());
var _Player = /** @class */ (function () {
    /**
     * The constructor of a _Player class.
     * @param  {*} parent The parent class (BrowserSim Object).
     */
    function _Player(parent) {
        this._startTime = undefined;
        this._timer = undefined;
        this._parent = undefined;
        this._parent = parent;
    }
    /**
     * Starts playing the simulation animation.
     * @param  {Object<string,*>} params Parameter from the user.
     * @returns void
     */
    _Player.prototype["start"] = function (params) {
        var _this = this;
        // loop: boolean = true, windowAverageSize: number = 1) {
        this._startTime = new Date().getTime();
        // Allow the user to change some of the parameters on start. So
        // duration doesn't always have to be fixed, for example.
        this._parent.updateParams(params);
        // Clear previous intervals.
        if (this._timer !== undefined) {
            clearInterval(this._timer);
        }
        this._timer = setInterval(function () {
            // How far along the animation are you?
            var curTime = new Date().getTime();
            var playRatio = (curTime - _this._startTime) / _this._parent._params["durationInMilliseconds"];
            // If you've gone over the end of the animation, start again from
            // the beginning.
            if ((!_this._parent._params["loop"]) && (playRatio > 1.0)) {
                _this["stop"]();
                return;
            }
            // Get the current frame.
            var curFrame = Math.floor(_this._parent._numFrames * playRatio);
            _this._parent["viewer"].updateAtomPos(curFrame);
        }, this._parent._params["updateFreqInMilliseconds"]);
    };
    /**
     * Stops the simulation animation.
     * @returns void
     */
    _Player.prototype["stop"] = function () {
        // Clear previous intervals.
        if (this._timer !== undefined) {
            clearInterval(this._timer);
        }
    };
    /**
     * Switches the simulation animation to a specific frame.
     * @param  {number} frame The frame.
     * @returns void
     */
    _Player.prototype["toFrame"] = function (frame) {
        this._parent["viewer"].updateAtomPos(frame);
    };
    return _Player;
}());
var MathUtils;
(function (MathUtils) {
    /**
     * Gets a list of numbers. Like Python's range function.
     * @param  {number} a     The first number.
     * @param  {number} b     The largest numer (plus one).
     * @returns Array<number> A list of the numbers.
     */
    function _range(a, b) {
        // Don't do this as a Float32Array, because map won't work as
        // expected. Map on typed array must return typed array.
        var rng = [];
        var i = a;
        while (i < b) {
            rng.push(i);
            i++;
        }
        return rng;
    }
    MathUtils._range = _range;
    /**
     * Given a list of Float32Array, returns the sum of them (reduce).
     * @param  {Array<Float32Array>} arrayOfFloat32Arrays The list of Float32Array.
     * @returns Float32Array         The sum.
     */
    function sumArrayOfFloat32Arrays(arrayOfFloat32Arrays) {
        if (arrayOfFloat32Arrays.length === 1) {
            // Just one item in the array, so return that first item.
            return arrayOfFloat32Arrays[0];
        }
        else {
            // Multiple items. So need to sum them.
            return arrayOfFloat32Arrays.reduce(function (summedVals, newVals) {
                return summedVals.map(function (v, i) {
                    return v + newVals[i];
                });
            });
        }
    }
    MathUtils.sumArrayOfFloat32Arrays = sumArrayOfFloat32Arrays;
    /**
     * Multiplies a Float32Array by a scalar.
     * @param  {number}       scalar        The scalar.
     * @param  {Float32Array} float32Array  The Float32Array.
     * @returns Float32Array  The multiplied Float32Array.
     */
    function multiplyFloat32ArrayByScalar(scalar, float32Array) {
        switch (scalar) {
            case 0.0:
                return new Float32Array(float32Array.length);
            case 1.0:
                return float32Array;
            default:
                return float32Array.map(function (v) { return scalar * v; });
        }
    }
    MathUtils.multiplyFloat32ArrayByScalar = multiplyFloat32ArrayByScalar;
    /**
     * To safe on space, the JSON file removes decimals. Restore those here
     * for a given number.
     * @param  {number} val        The value without a decimal point.
     * @param  {number} precision  The precision.
     * @returns number            The value of the number with the decimal point.
     */
    function convertToDecimal(val, precision) {
        return val * Math.pow(10, -precision);
    }
    MathUtils.convertToDecimal = convertToDecimal;
})(MathUtils || (MathUtils = {}));
window["BrowserSim"] = BrowserSim; // To survive closure compiler.
