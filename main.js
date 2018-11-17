var BrowserSim = /** @class */ (function () {
    function BrowserSim(params) {
        // Note that the below are public so they can be accessed from other
        // classes. But they are not public to the user and so do not need to be
        // protected from closure compiler.
        this._numFrames = undefined;
        this._frameSize = undefined;
        this._frameData = [];
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
    BrowserSim.prototype.updateParams = function (updatedParams) {
        // Set the defaults and keeps any previous ones not now specified.
        var defaults = {
            "viewer": undefined,
            "viewerType": undefined,
            "durationInMilliseconds": 10000,
            "updateFreqInMilliseconds": 10,
            "loop": true,
            "windowAverageSize": 1,
            "parent": this
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
                case "JSMOL":
                    this._params["visStyle"] = undefined;
                    break;
            }
        }
    };
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
    function _IO(parent) {
        this._parent = undefined;
        this._parent = parent;
    }
    _IO.prototype["loadJSON"] = function (filename, callBack) {
        var _this = this;
        if (callBack === void 0) { callBack = function () { }; }
        jQuery.getJSON(filename, function (data) {
            // Setup the frames
            _this._parent._numFrames = data["frames"].length;
            _this._parent._frameSize = data["frames"][0].length;
            _this._parent._frameData = [];
            // Set up the components
            _this._parent._numComponents = data["vectors"].length;
            _this._parent._componentSize = data["vectors"][0].length;
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
            for (var idx1 in data["frames"]) {
                if (data["frames"].hasOwnProperty(idx1)) {
                    _this._parent._frameData[idx1] = new Float32Array(data["frames"][idx1]);
                }
            }
            // Same with vectors.
            for (var idx1 in data["vectors"]) {
                if (data["vectors"].hasOwnProperty(idx1)) {
                    _this._parent._componentData[idx1] = new Float32Array(data["vectors"][idx1]);
                }
            }
        }).done(function () {
            callBack();
        }).fail(function () {
            console.log("error");
        }); /* .always(function() {
            console.log( "complete" );
        }); */
    };
    ;
    _IO.prototype["makePDB"] = function () {
        var frameIdx = 0;
        var numAtoms = this._parent._componentSize / 3;
        var pdbTxt = "";
        for (var frameIdx_1 = 0; frameIdx_1 < this._parent._numFrames; frameIdx_1++) {
            pdbTxt += "MODEL " + frameIdx_1.toString() + "\n";
            for (var atomIdx = 0; atomIdx < numAtoms; atomIdx++) {
                var coor = this._parent.getAtomCoors(frameIdx_1, atomIdx);
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
    _IO.prototype._formatNumForPDB = function (val) {
        var valStr = val.toFixed(3);
        return this._rjust(8, valStr);
    };
    ;
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
        var validTypes = ["3DMOLJS", "NGL", "JSMOL"];
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
    _Viewer.prototype.updateAtomPos = function (frame) {
        if (this._updateLocked) {
            return;
        }
        this._updateLocked = true;
        this._updateAtomPosFun(frame);
        this._updateLocked = false;
        this._render();
    };
    _Viewer.prototype._3DMolJS_AddPDBTxt = function (pdbTxt) {
        this._model = this._parent._params["viewer"].addModel(pdbTxt, "pdb");
        this._render();
    };
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
    _Viewer.prototype._3DMolJS_Render = function () {
        // Must update styles to actually have atoms move. Annoying.
        this._model.setStyle({}, this._parent._params["visStyle"]);
        this._parent._params["viewer"].render();
    };
    _Viewer.prototype._NGL_AddPDBTxt = function (pdbTxt) {
        var _this = this;
        this._parent._params["viewer"].loadFile(new Blob([pdbTxt], { type: 'text/plain' }), { ext: 'pdb', defaultRepresentation: true }).then(function (e) {
            _this._model = e;
        });
    };
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
    _Viewer.prototype._NGL_Render = function () {
        this._model["rebuildRepresentations"]();
    };
    _Viewer.prototype._JSMOL_AddPDBTxt = function (pdbTxt) {
        // let jsmolCmd = "data \"model browser_sim\"\n";
        // jsmolCmd += pdbTxt + "\n";
        // jsmolCmd += 'end "model browser_sim";show data "pdb"';
        var _this = this;
        // See
        // https://jmol-developers.narkive.com/UJZfiMkC/ie-specific-problem-with-jmolloadinline
        // for example.
        var test = "2\ntesting\nC 1 1 1\nO 2 2 2\n";
        // jQuery("body").append("<div id='test'>" + test + "</div>");
        setTimeout(function () {
            // window["t"] = test.split("\n");
            // let jsmolCmd = "jmolLoadInline('" + pdbTxt.replace(/\n/g, "\\n") + "');";
            var jsmolCmd = "load data \"model example\"\nATOM     30  X   XXX X 852       1.020  -0.246  -0.345  1.00  0.00           X\nATOM     31  X   XXX X 852      -1.570   0.836  -0.115  1.00  0.00           X\nATOM     32  X   XXX X 852       0.000  -0.738   0.207  1.00  0.00           X\nATOM     33  X   XXX X 852       0.785  -0.344   0.069  1.00  0.00           X\nATOM     34  X   XXX X 852      -1.727  -0.246   0.253  1.00  0.00           X\nATOM     35  X   XXX X 852      -0.157  -0.197   0.000  1.00  0.00           X\nATOM     36  X   XXX X 852       0.549   0.984  -0.391  1.00  0.00           X\nATOM     37  X   XXX X 852       0.471   0.295   0.138  1.00  0.00           X\nend \"model example\";";
            console.log(jsmolCmd);
            // jsmolCmd = "background red";
            _this._parent._params["viewer"]["library"].script(_this._parent._params["viewer"]["applet"], jsmolCmd);
        }, 1000);
        // this._model = this._parent._params["viewer"].addModel( pdbTxt, "pdb" );
        // this._render();
        // viewerJSMol = {
        //     applet: jsmolApplet,
        //     library: Jmol
        // }
    };
    _Viewer.prototype._JSMOL_UpdateAtomPos = function (frame) {
    };
    _Viewer.prototype._JSMOL_Render = function () {
    };
    return _Viewer;
}());
var _Player = /** @class */ (function () {
    function _Player(parent) {
        this._startTime = undefined;
        this._timer = undefined;
        this._parent = undefined;
        this._parent = parent;
    }
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
    _Player.prototype["stop"] = function () {
        // Clear previous intervals.
        if (this._timer !== undefined) {
            clearInterval(this._timer);
        }
    };
    _Player.prototype["toFrame"] = function (frame) {
        this._parent["viewer"].updateAtomPos(frame);
    };
    return _Player;
}());
var MathUtils;
(function (MathUtils) {
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
})(MathUtils || (MathUtils = {}));
window["BrowserSim"] = BrowserSim; // To survive closure compiler.
