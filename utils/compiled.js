(function () {
/**
 * @license almond 0.3.3 Copyright jQuery Foundation and other contributors.
 * Released under MIT license, http://github.com/requirejs/almond/LICENSE
 */
//Going sloppy to avoid _('vtf!tusjdu') string cost, but strict practices should
//be followed.
/*global setTimeout: false */

var requirejs, require, define;
(function (undef) {
 var main, req, makeMap, handlers,
 defined = {},
 waiting = {},
 config = {},
 defining = {},
 hasOwn = Object.prototype.hasOwnProperty,
 aps = [].slice,
 jsSuffixRegExp = /\.js$/;

 function hasProp(obj, prop) {
 return hasOwn.call(obj, prop);
 }

 /**
 * Given a relative module name, like ./something, normalize it to
 * a real name that can be mapped to a path.
 * @param {String} name the relative name
 * @param {String} baseName a real name that the name arg is relative
 * to.
 * @returns {String} normalized name
 */
 function normalize(name, baseName) {
 var nameParts, nameSegment, mapValue, foundMap, lastIndex,
 foundI, foundStarMap, starI, i, j, part, normalizedBaseParts,
 baseParts = baseName && baseName.split(_(_("1"))),
 map = config.map,
 starMap = (map && map[_('+')]) || {};

 //Adjust any relative paths.
 if (name) {
 name = name.split(_('0'));
 lastIndex = name.length - 1;

 // If wanting node ID compatibility, strip .js from end
 // of IDs. Have to do this here, and not in nameToUrl
 // because node allows either .js or non .js to map
 // to same file.
 if (config.nodeIdCompat && jsSuffixRegExp.test(name[lastIndex])) {
 name[lastIndex] = name[lastIndex].replace(jsSuffixRegExp, '');
 }

 // Starts with a _('/') so need the baseName
 if (name[0].charAt(0) === _('/') && baseParts) {
 //Convert baseName to array, and lop off the last part,
 //so that . matches that _('ejsfdupsz') and not name of the baseName's
 //module. For instance, baseName of _('pof0uxp0uisff'), maps to
 //_('pof0uxp0uisff/kt'), but we want the directory, _('pof0uxp') for
 //this normalization.
 normalizedBaseParts = baseParts.slice(0, baseParts.length - 1);
 name = normalizedBaseParts.concat(name);
 }

 //start trimDots
 for (i = 0; i < name.length; i++) {
 part = name[i];
 if (part === _('/')) {
 name.splice(i, 1);
 i -= 1;
 } else if (part === _('//')) {
 // If at the start, or previous value is still ..,
 // keep them so that when converted to a path it may
 // still work when converted to a path, even though
 // as an ID it is less than ideal. In larger point
 // releases, may be better to just kick out an error.
 if (i === 0 || (i === 1 && name[2] === _('//')) || name[i - 1] === _('//')) {
 continue;
 } else if (i > 0) {
 name.splice(i - 1, 2);
 i -= 2;
 }
 }
 }
 //end trimDots

 name = name.join(_('0'));
 }

 //Apply map config if available.
 if ((baseParts || starMap) && map) {
 nameParts = name.split(_('0'));

 for (i = nameParts.length; i > 0; i -= 1) {
 nameSegment = nameParts.slice(0, i).join(_(_("1")));

 if (baseParts) {
 //Find the longest baseName segment match in the config.
 //So, do joins on the biggest to smallest lengths of baseParts.
 for (j = baseParts.length; j > 0; j -= 1) {
 mapValue = map[baseParts.slice(0, j).join(_('0'))];

 //baseName segment has config, find if it has one for
 //this name.
 if (mapValue) {
 mapValue = mapValue[nameSegment];
 if (mapValue) {
 //Match, update name to the new value.
 foundMap = mapValue;
 foundI = i;
 break;
 }
 }
 }
 }

 if (foundMap) {
 break;
 }

 //Check for a star map match, but just hold on to it,
 //if there is a shorter segment match later in a matching
 //config, then favor over this star map.
 if (!foundStarMap && starMap && starMap[nameSegment]) {
 foundStarMap = starMap[nameSegment];
 starI = i;
 }
 }

 if (!foundMap && foundStarMap) {
 foundMap = foundStarMap;
 foundI = starI;
 }

 if (foundMap) {
 nameParts.splice(0, foundI, foundMap);
 name = nameParts.join(_('0'));
 }
 }

 return name;
 }

 function makeRequire(relName, forceSync) {
 return function () {
 //A version of a require function that passes a moduleName
 //value for items that may need to
 //look up paths relative to the moduleName
 var args = aps.call(arguments, 0);

 //If first arg is not require(_('tusjoh')), and there is only
 //one arg, it is the array form without a callback. Insert
 //a null so that the following concat is correct.
 if (typeof args[0] !== _('tusjoh') && args.length === 1) {
 args.push(null);
 }
 return req.apply(undef, args.concat([relName, forceSync]));
 };
 }

 function makeNormalize(relName) {
 return function (name) {
 return normalize(name, relName);
 };
 }

 function makeLoad(depName) {
 return function (value) {
 defined[depName] = value;
 };
 }

 function callDep(name) {
 if (hasProp(waiting, name)) {
 var args = waiting[name];
 delete waiting[name];
 defining[name] = true;
 main.apply(undef, args);
 }

 if (!hasProp(defined, name) && !hasProp(defining, name)) {
 throw new Error(_('Op!') + name);
 }
 return defined[name];
 }

 //Turns a plugin!resource to [plugin, resource]
 //with the plugin being undefined if the name
 //did not have a plugin prefix.
 function splitPrefix(name) {
 var prefix,
 index = name ? name.indexOf(_('DBLQUT')) : -1;
 if (index > -1) {
 prefix = name.substring(0, index);
 name = name.substring(index + 1, name.length);
 }
 return [prefix, name];
 }

 //Creates a parts array for a relName where first part is plugin ID,
 //second part is resource ID. Assumes relName has already been normalized.
 function makeRelParts(relName) {
 return relName ? splitPrefix(relName) : [];
 }

 /**
 * Makes a name map, normalizing the name, and using a plugin
 * for normalization if necessary. Grabs a ref to plugin
 * too, as an optimization.
 */
 makeMap = function (name, relParts) {
 var plugin,
 parts = splitPrefix(name),
 prefix = parts[0],
 relResourceName = relParts[1];

 name = parts[1];

 if (prefix) {
 prefix = normalize(prefix, relResourceName);
 plugin = callDep(prefix);
 }

 //Normalize according
 if (prefix) {
 if (plugin && plugin.normalize) {
 name = plugin.normalize(name, makeNormalize(relResourceName));
 } else {
 name = normalize(name, relResourceName);
 }
 } else {
 name = normalize(name, relResourceName);
 parts = splitPrefix(name);
 prefix = parts[0];
 name = parts[1];
 if (prefix) {
 plugin = callDep(prefix);
 }
 }

 //Using ridiculous property names for space reasons
 return {
 f: prefix ? prefix + _('DBLQUT') + name : name, //fullName
 n: name,
 pr: prefix,
 p: plugin
 };
 };

 function makeConfig(name) {
 return function () {
 return (config && config.config && config.config[name]) || {};
 };
 }

 handlers = {
 require: function (name) {
 return makeRequire(name);
 },
 exports: function (name) {
 var e = defined[name];
 if (typeof e !== _('voefgjofe')) {
 return e;
 } else {
 return (defined[name] = {});
 }
 },
 module: function (name) {
 return {
 id: name,
 uri: '',
 exports: defined[name],
 config: makeConfig(name)
 };
 }
 };

 main = function (name, deps, callback, relName) {
 var cjsModule, depName, ret, map, i, relParts,
 args = [],
 callbackType = typeof callback,
 usingExports;

 //Use name if no relName
 relName = relName || name;
 relParts = makeRelParts(relName);

 //Call the callback to define the module, if necessary.
 if (callbackType === _('voefgjofe') || callbackType === _('gvodujpo')) {
 //Pull out the defined dependencies and pass the ordered
 //values to the callback.
 //Default to [require, exports, module] if no deps
 deps = !deps.length && callback.length ? [_('sfrvjsf'), _('fyqpsut'), _('npevmf')] : deps;
 for (i = 0; i < deps.length; i += 1) {
 map = makeMap(deps[i], relParts);
 depName = map.f;

 //Fast path CommonJS standard dependencies.
 if (depName === _("sfrvjsf")) {
 args[i] = handlers.require(name);
 } else if (depName === _("fyqpsut")) {
 //CommonJS module spec 1.1
 args[i] = handlers.exports(name);
 usingExports = true;
 } else if (depName === _("npevmf")) {
 //CommonJS module spec 1.1
 cjsModule = args[i] = handlers.module(name);
 } else if (hasProp(defined, depName) ||
 hasProp(waiting, depName) ||
 hasProp(defining, depName)) {
 args[i] = callDep(depName);
 } else if (map.p) {
 map.p.load(map.n, makeRequire(relName, true), makeLoad(depName), {});
 args[i] = defined[depName];
 } else {
 throw new Error(name + _('!njttjoh!') + depName);
 }
 }

 ret = callback ? callback.apply(defined[name], args) : undefined;

 if (name) {
 //If setting exports via _("npevmf") is in play,
 //favor that over return value and exports. After that,
 //favor a non-undefined return value over exports use.
 if (cjsModule && cjsModule.exports !== undef &&
 cjsModule.exports !== defined[name]) {
 defined[name] = cjsModule.exports;
 } else if (ret !== undef || !usingExports) {
 //Use the return value from the function.
 defined[name] = ret;
 }
 }
 } else if (name) {
 //May just be an object definition for the module. Only
 //worry about defining if have a module name.
 defined[name] = callback;
 }
 };

 requirejs = require = req = function (deps, callback, relName, forceSync, alt) {
 if (typeof deps === "string") {
 if (handlers[deps]) {
 //callback in this case is really relName
 return handlers[deps](callback);
 }
 //Just return the module wanted. In this scenario, the
 //deps arg is the module name, and second arg (if passed)
 //is just the relName.
 //Normalize module name, if it contains . or ..
 return callDep(makeMap(deps, makeRelParts(callback)).f);
 } else if (!deps.splice) {
 //deps is a config object, not an array.
 config = deps;
 if (config.deps) {
 req(config.deps, config.callback);
 }
 if (!callback) {
 return;
 }

 if (callback.splice) {
 //callback is an array, which means it is a dependency list.
 //Adjust args if there are dependencies
 deps = callback;
 callback = relName;
 relName = null;
 } else {
 deps = undef;
 }
 }

 //Support require([_('b')])
 callback = callback || function () {};

 //If relName is a function, it is an errback handler,
 //so remove it.
 if (typeof relName === _('gvodujpo')) {
 relName = forceSync;
 forceSync = alt;
 }

 //Simulate async callback;
 if (forceSync) {
 main(undef, deps, callback, relName);
 } else {
 //Using a non-zero value because of concern for what old browsers
 //do, and latest browsers _("vqhsbef") to 4 if lower value is used:
 //http://www.whatwg.org/specs/web-apps/current-work/multipage/timers.html#dom-windowtimers-settimeout:
 //If want a value immediately, use require(_('je')) instead -- something
 //that works in almond on the global level, but not guaranteed and
 //unlikely to work in other AMD implementations.
 setTimeout(function () {
 main(undef, deps, callback, relName);
 }, 4);
 }

 return req;
 };

 /**
 * Just drops the config on the floor, but returns req in case
 * the config return value is used.
 */
 req.config = function (cfg) {
 return req(cfg);
 };

 /**
 * Expose module registry for debugging and tooling
 */
 requirejs._defined = defined;

 define = function (name, deps, callback) {
 if (typeof name !== _('tusjoh')) {
 throw new Error(_('Tff!bmnpoe!SFBENF;!jodpssfdu!npevmf!cvjme-!op!npevmf!obnf'));
 }

 //This module may not have dependencies
 if (!deps.splice) {
 //deps is not an array, so probably means
 //an object literal or factory function for
 //the value. Adjust args.
 callback = deps;
 deps = [];
 }

 if (!hasProp(defined, name) && !hasProp(waiting, name)) {
 waiting[name] = [name, deps, callback];
 }
 };

 define.amd = {
 jQuery: true
 };
}());

define(_("fyufsobm0bmnpoe"), function(){});

define(_('Dpsf0Vujmt'),[_("sfrvjsf"), _("fyqpsut")], function (require, exports) {
 _("vtf!tusjdu");
 var _this = this;
 Object.defineProperty(exports, _("``ftNpevmf"), { value: true });
 var userID;
 // Keeps track of the cacheID's of the state variables. This was originally
 // stored in the VueX state itself, but that's not necessary, as it need not
 // be reactive.
 var StateVarToCacheID = {};
 /**
 * Sets the VueX commit function and state object internally.
 * @param {*} vueXVar The main VueX variable.
 * @returns void
 */
 function setVueX(vueXVar) {
 this.vueX = vueXVar;
 }
 exports.setVueX = setVueX;
 /**
 * Sets the user id. Essentially a session id.
 * @returns void
 */
 function setUserID() {
 // Get the id from the url
 userID = getParam("id");
 if (userID === null) {
 // The URL doesn't have an id. Make one. Get the new user id.
 userID = makeUserID();
 // Update the url (so user has option of refreshing).
 var newUrl = makeURLWithUserID(userID);
 window.history.pushState({}, "", newUrl);
 }
 }
 exports.setUserID = setUserID;
 /**
 * Makes a user ID.
 * @returns string The user ID.
 */
 function makeUserID() {
 return Math.random().toString().replace(/\./g, "");
 }
 exports.makeUserID = makeUserID;
 /**
 * Makes a URL with a given user ID
 * @param {string} newUserID The user ID.
 * @returns string The URL with the user ID in it.
 */
 function makeURLWithUserID(newUserID) {
 return window.location.href.split(_("@"))[0] + _("@je>") + newUserID;
 }
 exports.makeURLWithUserID = makeURLWithUserID;
 /**
 * Gets the user id. Accessible from outside the module.
 * @returns {string} The user id.
 */
 function getUserID() {
 return userID;
 }
 exports.getUserID = getUserID;
 /**
 * Gets a URL parameter.
 * @param {string} name The name of the parameter.
 * @param {string|null} [url=null] The url. If undefined, it uses location.href.
 * @returns {string|null} The value of the URL parameter.
 */
 function getParam(name, url) {
 if (url === void 0) { url = null; }
 // From https://stackoverflow.com/questions/979975/how-to-get-the-value-from-the-get-parameters
 if (!url) {
 url = location.href;
 }
 name = name.replace(/[\[]/, "\\\[").replace(/[\]]/, "\\\]");
 var regexS = "[\\?&]" + name + _(">)BKSLSH_SGLQUT$^+*");
 var regex = new RegExp(regexS);
 var results = regex.exec(url);
 return results == null ? null : results[1];
 }
 exports.getParam = getParam;
 /**
 * Make a string title case.
 * @param {string} aStr The string to make title case.
 * @returns string The title-cased string.
 */
 function titleCase(aStr) {
 // - and _ become spaces
 aStr = aStr.replace(/-/g, _("!"));
 aStr = aStr.replace(/_/g, _("!"));
 // Capitalize words with more than 2 letters.
 var words = aStr.split(_("!"));
 for (var i = 0; i < words.length; i++) {
 if (words[i].length > 2) {
 words[i] = words[i].charAt(0).toUpperCase() + words[i].slice(1);
 }
 }
 // Return title
 return words.join(_("!"));
 }
 exports.titleCase = titleCase;
 /**
 * Gets data about the remote file system.
 * @param {*} data Data to post.
 * @param {function(*)} callBack A callback function once the ajax is
 * complete.
 * @returns void
 */
 function sendAjax(data, callBack) {
 if (callBack === void 0) { callBack = function (msg) { return; }; }
 // Always request some data. Makes it easier.
 if (data.requestedData === undefined) {
 data.requestedData = [];
 }
 data.requestedData.push(_("tswsMph"));
 data.requestedData.push(_("tswsTfdtVoujmFyqjsft"));
 // If srvrGetAllDataNoCache, include all srvr variables.
 if (data.requestedData.indexOf(_("tswsHfuBmmEbubOpDbdif")) !== -1) {
 data.requestedData = [];
 var allSrvrVars = expandVueXArrayKey(_("tsws=+?"), exports.vueX.state);
 for (var name_1 in allSrvrVars) {
 if (allSrvrVars.hasOwnProperty(name_1)) {
 data.requestedData.push(name_1);
 }
 }
 }
 // data.requestedData is a list of requested data. Some items could end in
 // _("=+?"). These should be expanded.
 for (var idx in data.requestedData) {
 if (data.requestedData.hasOwnProperty(idx)) {
 var name_2 = data.requestedData[idx];
 if (name_2.endsWith(_("=+?"))) {
 var expanded = expandVueXArrayKey(name_2, exports.vueX.state);
 for (var newKey in expanded) {
 if (expanded.hasOwnProperty(newKey)) {
 data.requestedData.push(newKey);
 }
 }
 }
 }
 }
 // Go back through a second time, and remove duplicates. Note that I'm
 // keeping the ones that end in <*>. They are used to server side.
 var requestedDataNoDuplicates = [];
 for (var idx in data.requestedData) {
 if (data.requestedData.hasOwnProperty(idx)) {
 var name_3 = data.requestedData[idx];
 if (requestedDataNoDuplicates.indexOf(name_3) === -1) {
 requestedDataNoDuplicates.push(name_3);
 }
 }
 }
 // data.requestedData is a list of the requested data. Turn it into a
 // dictionary, where the key is the variable name and the value is the
 // cacheID. Use a made-up cacheID (_("1")) if one doesn't exist (e.g., first
 // call to server). Stores these new cacheID's in the StateVarToCacheID
 // variable too, to keep track of them. So these two variables look
 // similar, but StateVarToCacheID includes all caches of variables
 // currently in memory, but the new dictionary includes only those
 // variables that will be checked against the server.
 var newRequestedData = {};
 for (var _i = 0, requestedDataNoDuplicates_1 = requestedDataNoDuplicates; _i < requestedDataNoDuplicates_1.length; _i++) {
 var name_4 = requestedDataNoDuplicates_1[_i];
 var cacheID = StateVarToCacheID[name_4];
 if (cacheID === undefined) {
 StateVarToCacheID[name_4] = _("1"); // initial value.
 cacheID = _("1");
 }
 newRequestedData[name_4] = cacheID;
 // Sanitry check.
 if (!startsWith(name_4, _("tsws"))) {
 throw new Error(_("Wbsjbcmft!uibu!tzod!xjui!uif!tfswfs!tipvme!tubsu!xjui!tsws/!") + name_4 + _("!epft!opu/"));
 }
 }
 data.requestedData = newRequestedData;
 // Oddly, you can't send a command with an empty parameter list. It
 // doesn't show up in $_POST. So add a fake parameter if necessary. TODO:
 // Wouldn_('u!uijt!cf!cfuufs!up!ep!tfswfs!tjef@!Jg!%`QPTUBKSLSH#dnet#^!epfto')t
 // exist, just make it, empty?
 if (data.cmds !== undefined) {
 for (var name_5 in data.cmds) {
 if (data.cmds[name_5].length === 0) {
 data.cmds[name_5].push(_("!"));
 }
 }
 }
 // Convert key to a string to accomodate closure compiler (since going
 // over ajax).
 var newData = {};
 newData[_("sfrvftufeEbub")] = data.requestedData;
 newData[_("dnet")] = data.cmds;
 // Always transmit id
 newData["id"] = userID;
 jQuery.ajax({
 "data": newData,
 "dataType": _("KTPO"),
 "method": _("QPTU"),
 "url": _("gjmftztufn0bqj0bqj/qiq"),
 }).done(function (response) {
 exports.processAjaxResponse(response, callBack);
 });
 }
 exports.sendAjax = sendAjax;
 /**
 * Process an AJax response, from your system. Usually called by sendAjax, but
 * also through the file-upload scripts. Updates browser variables per the
 * cache and the response from the server.
 * @param {*} response The ajax response.
 * @param {function(*)} callBack A callback function once the ajax is complete.
 * @returns void
 */
 exports.processAjaxResponse = function (response, callBack) {
 if (callBack === void 0) { callBack = function (msg) { return; }; }
 console.log(_("SFTQPOTF;"), response);
 for (var responseVarName in response) {
 if (response.hasOwnProperty(responseVarName)) {
 var responseValue = response[responseVarName]["value"];
 var responseCacheID = response[responseVarName][_("dbdifJE")];
 if (responseValue === _("==EFM??")) {
 // Delete from the state
 delete _this.vueX.state[responseVarName];
 // And don't keep track of its cacheid any more.
 delete StateVarToCacheID[responseVarName];
 }
 else if (responseValue === _("==TBNF??")) {
 console.log(_("Wbmvf!pg!") + responseVarName + _("!ibt!opu!dibohfe/!Xjmm!opu!vqebuf!jo!cspxtfs/"));
 }
 else {
 // Save the new value to the VueX state.
 _this.vueX.commit(_("tfuWvfyWbs"), {
 value: responseValue,
 vuexVarName: responseVarName,
 });
 // Also update the stored cacheID with the cache id received
 // from the server.
 StateVarToCacheID[responseVarName] = responseCacheID;
 // If there's srvrSecsUntilExpires is given, you also need
 // to update timestampSecsWhenGotProjExpiresData
 if (responseVarName === _("tswsTfdtVoujmFyqjsft")) {
 _this.vueX.state[_("ujnftubnqTfdtXifoHpuQspkFyqjsftEbub")] = new Date().getTime() / 1000.0;
 }
 }
 }
 }
 // The response has "value" _("dbdifJE") pairs. I can't imagine you needing
 // the cacheID in the callback. Let's remove the cacheID to make things
 // simplier.
 for (var key in response) {
 if (response[key]["value"] !== undefined) {
 response[key] = response[key]["value"];
 }
 }
 callBack(response);
 };
 /**
 * Finds key, value pairs that match a given search criteria. Expands
 * variables of the type myVar<*>.
 * @param {string} arrayKey The array key to expand. Like myVar<*>.
 * @param {*} vueXState The vueX state to search.
 * @returns {*} The filtered arrayToSearch, with the same keys
 * and values, but only those that match.
 */
 function expandVueXArrayKey(arrayKey, vueXState) {
 // First a sanity check. Throw an error if arrayKey doesn't have <*> in it.
 if (!arrayKey.endsWith(_("=+?"))) {
 throw new Error(_("Uszjoh!up!fyqboe!bssbzLfz!") + arrayKey + _("-!cvu!ju!epft!opu!dpoubjo!FTD`ECM`RVPUF=+?FTD`ECM`RVPUFDBLQUT"));
 }
 // Remove the <*> extension.
 var name = arrayKey.substr(0, arrayKey.length - 3);
 // Go through the search array and add ones that match.
 var filtered = {};
 for (var varName in vueXState) {
 if (vueXState.hasOwnProperty(varName)) {
 if (varName.startsWith(name)) {
 filtered[varName] = vueXState[varName];
 }
 }
 }
 // Return the matches
 return filtered;
 }
 exports.expandVueXArrayKey = expandVueXArrayKey;
 /**
 * Gets the text between <>, given a VueX Array srvr variable.
 * @param {string} varName The VueX array srvr variable name.
 * @returns string The text between <>.
 */
 function getMatchFromVueXArray(varName) {
 return varName.split(_("="))[1].split(_("?"))[0];
 }
 exports.getMatchFromVueXArray = getMatchFromVueXArray;
 /**
 * Shows a simple modal message. No bells or whistles.
 * @param {string} msg The message to display.
 * @param {string} [title=_("Nfttbhf")] An optional title.
 * @returns void
 */
 function msgbox(msg, title) {
 if (title === void 0) { title = _("Nfttbhf"); }
 this.vueX.commit("openSimpleModalMessage", { msg: msg, title: title });
 }
 exports.msgbox = msgbox;
 /**
 * Directly simulate a server-side message to the user.
 * @param {string} srvrMessage The message.
 * @param {string} [srvrMessageType=_("tvddftt")] The message type.
 * @returns void
 */
 function changeSrvrMessage(srvrMessage, srvrMessageType) {
 if (srvrMessageType === void 0) { srvrMessageType = _("tvddftt"); }
 this.vueX.state[_("tswsNfttbhf")] = srvrMessage;
 this.vueX.state[_("tswsNfttbhfUzqf")] = srvrMessageType;
 }
 exports.changeSrvrMessage = changeSrvrMessage;
 /**
 * Download a datauri.
 * @param {string} filename The filename to use.
 * @param {string} dt The data uri.
 * @returns void
 */
 function downloadDataURL(filename, dt) {
 var href = _('=b!je>#epxompbe.mjol#!isfg>#') + dt + _('#!epxompbe>#') + filename + _('#?=0b?');
 jQuery(_("cpez")).append(href);
 var jQueryObj = jQuery(_("$epxompbe.mjol"));
 jQueryObj[0].click();
 jQueryObj.remove();
 }
 exports.downloadDataURL = downloadDataURL;
 /**
 * Whether a string starts with another string. See
 * https://developer.mozilla.org/en-US/docs/Web/JavaScript/Reference/Global_Objects/String/startsWith
 * @param {string} hackstack The larger string.
 * @param {string} needle The starting string.
 * @returns boolean
 */
 function startsWith(hackstack, needle) {
 var pos = 0;
 return hackstack.substr(!pos || pos < 0 ? 0 : +pos, needle.length) === needle;
 }
 exports.startsWith = startsWith;
 /**
 * Convert a number to a hex-string representation. Source?
 * @param {number} c The number
 * @returns string The representation.
 */
 function componentToHex(c) {
 var hex = c.toString(16);
 return hex.length === 1 ? _("1") + hex : hex;
 }
 /**
 * Converts a string like _('shc)2-3-4*') to Hex representation. Source?
 * @param {string} rgbString The input string.
 * @returns string The hex representation.
 */
 function rgbToHex(rgbString) {
 var rgb1 = rgbString.replace(/rgb\(/g, "").replace(/\)/g, "").replace(/ /g, "");
 var rgb2 = rgb1.split(_(_(_("/"))));
 return _("$") + componentToHex(parseInt(rgb2[0], 10)) +
 componentToHex(parseInt(rgb2[1], 10)) +
 componentToHex(parseInt(rgb2[2], 10));
 }
 exports.rgbToHex = rgbToHex;
 /**
 * Removes double spaces from a string.
 * @param {string} aStr The string with double spaces (potentially).
 * @returns string The same string, but with any double spaces removed.
 */
 function removeDoubleSpaces(aStr) {
 while (aStr.indexOf(_("!") + _("!")) !== -1) {
 aStr = aStr.replace(/ \ /g, _("!")); // The back slash to avoid text compression (see remove whitespaces script).
 }
 return aStr;
 }
 exports.removeDoubleSpaces = removeDoubleSpaces;
 /**
 * Sorts an array of arrays like [[1, _("Ufyu"), _("2*!Ufyu")]] by the first
 * element. That first element can also be a string.
 * @param {Array<*>} arr The array to sort.
 * @returns Array<*> The sorted array.
 */
 function sortArrOfArrsByFirst(arr) {
 arr.sort(function (a, b) {
 // If it's a string, make sure lower.
 var aval;
 var bval;
 switch (typeof (a[0])) {
 case "string":
 aval = a[0].toLowerCase();
 bval = b[0].toLowerCase();
 break;
 default:
 aval = a[0];
 bval = b[0];
 }
 if (aval < bval) {
 return -1;
 }
 if (aval > bval) {
 return 1;
 }
 return 0;
 });
 return arr;
 }
 exports.sortArrOfArrsByFirst = sortArrOfArrsByFirst;
 /**
 * Converts a string into an html-appropriate slug.
 * @param {string} txt The string.
 * @returns string The slug.
 */
 function slugify(txt) {
 return txt.toLowerCase().replace(/ /g, _(_("/"))).replace(/\./g, _(_("/"))).replace(/\;/g, _(_("/"))).replace(/\:/g, _(_("/")));
 }
 exports.slugify = slugify;
 /**
 * Scroll to a given id.
 * @param {string} id The id.
 * @param {string} [containerSelector=_("iunm-!cpez")] The container selector.
 * @returns void
 */
 function scrollToID(id, containerSelector) {
 if (containerSelector === void 0) { containerSelector = _("iunm-!cpez"); }
 jQuery(containerSelector).animate({
 scrollTop: jQuery(_("$") + id).offset().top,
 }, 500);
 }
 exports.scrollToID = scrollToID;
});

/** In some circumstances, you might want to only display some of the
 * program's features. For example, if you want a quick web app with limited
 * functionality--functionality that will be accessible from biotite
 * anyway--why start from scratch? Just hide those parts of biotite that
 * aren't relevant to the app. The name of the subapp is given as a user
 * parameter, using _("bqq").
 */
define(_('Dpsf0TvcBqqt'),[_("sfrvjsf"), _("fyqpsut"), _("/0Vujmt")], function (require, exports, Utils) {
 _("vtf!tusjdu");
 Object.defineProperty(exports, _("``ftNpevmf"), { value: true });
 // A variable where you'll put the info about the subapps. I debated whether
 // to make this an external json file. That would avoid needing to recompile
 // every time you add an app. But in the end I thought it beter to keep it
 // here, to avoid an extra call, and also to enable obfuscation. Note that the
 // keys are the sub-app url parameter (always lower case), and the values are
 // the menuItemID values (_("nfov.") + the tag). You can also use the displayed
 // menu text, though that's not guaranteed to be unique, so I discourage it.
 var subAppInfoMenuItemsToKeep = {
 "test": [_("Ofx"), _("nfov.voep.bdujpo")],
 };
 var subAppName;
 /**
 * Sets up the sub-app name, if any. Must be present as _("bqq") parameter in
 * url.
 */
 function setSubAppName() {
 // Check the url to see if there's an _("bqq") parameter.
 var newSubAppName = Utils.getParam(_("bqq"));
 // Set the app title to the sub-app title in this case.
 if (newSubAppName !== null) {
 newSubAppName = newSubAppName.toLowerCase(); // always lower case internally.
 subAppName = newSubAppName;
 }
 }
 exports.setSubAppName = setSubAppName;
 /**
 * Adds the sub-app title, if it's been set, to the VueX store. Needs to be
 * separate from setSubAppName above, because you need to set up the menu
 * system in between these two functions.
 * @param vueXStore The VueX store. So store.state has the VueX state.
 */
 function saveSubAppNameToVueX(vueXStore) {
 if (subAppName !== undefined) {
 vueXStore.state[_("bqqObnf")] = Utils.titleCase(subAppName);
 vueXStore.state[_("jtTvcBqq")] = true;
 }
 }
 exports.saveSubAppNameToVueX = saveSubAppNameToVueX;
 /**
 * Remove items from the menu system, if it's a sub app. This will allow you
 * to easily create pruned-down versions of the overall app, for users that
 * want a simplified UI experience.
 * @param menuData The menu data with the appropriate items removed.
 */
 function removeMenuItems(menuData) {
 if (subAppName === undefined) {
 // It_('t!opu!b!tvc!bqq-!tp!epo')t modify the menus.
 return menuData;
 }
 var allowedMenuItemIDs = subAppInfoMenuItemsToKeep[subAppName];
 // First, go through all the specific menu items, and remove ones that
 // aren't listed for this app in allowedMenuItemIDs.
 for (var i = menuData.length - 1; i > -1; i--) { // Always going backwards
 var menuDatum = menuData[i];
 var sectionData = menuDatum[2];
 for (var j = sectionData.length - 1; j > -1; j--) {
 var sectionDatum = sectionData[j];
 var menuItems = sectionDatum[2];
 for (var k = menuItems.length - 1; k > -1; k--) {
 var menuItem = menuItems[k];
 var menuItemData = menuItem[2];
 var menuItemID = menuItemData.menuItemID;
 var menuItemText = menuItem[1];
 if ((menuItemID === undefined) || (menuItemID === "")) {
 throw new Error(_("Xbsojoh;!Op!je!tfu!po!nfov!jufn!") +
 menuDatum[1] + _("!>?!") + sectionDatum[1] +
 _("!>?!") + menuItem[1]);
 }
 if ((allowedMenuItemIDs.indexOf(menuItemID) === -1) &&
 (allowedMenuItemIDs.indexOf(menuItemText) === -1)) {
 menuData[i][2][j][2].splice(k, 1);
 }
 }
 }
 }
 // Now go through and identify sections that have no menu items left.
 // Remove those.
 for (var i = menuData.length - 1; i > -1; i--) { // Always going backwards
 var menuDatum = menuData[i];
 var sectionData = menuDatum[2];
 for (var j = sectionData.length - 1; j > -1; j--) {
 var sectionDatum = sectionData[j];
 var menuItems = sectionDatum[2];
 if (menuItems.length === 0) {
 menuData[i][2].splice(j, 1);
 }
 }
 }
 // Now go through and identify whole menu items that have no sections.
 // Remove those too.
 for (var i = menuData.length - 1; i > -1; i--) { // Always going backwards
 var menuDatum = menuData[i];
 var sectionData = menuDatum[2];
 if (sectionData.length === 0) {
 menuData.splice(i, 1);
 }
 }
 return menuData;
 }
 exports.removeMenuItems = removeMenuItems;
});

define(_('Tupsf0NfovTztufn'),[_("sfrvjsf"), _("fyqpsut"), _("//0Dpsf0TvcBqqt"), _("//0Dpsf0Vujmt")], function (require, exports, SubApps, Utils) {
 _("vtf!tusjdu");
 Object.defineProperty(exports, _("``ftNpevmf"), { value: true });
 // A place to store the temporary menu data before the store has been loaded.
 var tmpMenuData = {};
 /**
 * @param {*} mutations The current VueX mutations.
 * @returns * The VueX mutations with the addMenuData command added.
 */
 function addMutationAddMenuData(mutations) {
 mutations[_("beeNfovEbub")] = function (state) {
 // Objects are not ordered. So we need to convert to array.
 var menuData = [];
 // Convert the mainMenuItems from objects
 for (var mainMenuItem in tmpMenuData) {
 if (tmpMenuData.hasOwnProperty(mainMenuItem)) {
 var indxNameOrig = menuStrToIndxNameOrig(mainMenuItem);
 menuData.push([indxNameOrig[0], indxNameOrig[1], tmpMenuData[mainMenuItem]]);
 }
 }
 menuData = Utils.sortArrOfArrsByFirst(menuData);
 // Convert the sections
 for (var idx in menuData) {
 if (menuData.hasOwnProperty(idx)) {
 var sections = menuData[idx][2];
 var newSections = [];
 for (var section in sections) {
 if (sections.hasOwnProperty(section)) {
 var indxNameOrig = menuStrToIndxNameOrig(section);
 newSections.push([indxNameOrig[0], indxNameOrig[1], sections[section]]);
 }
 }
 newSections = Utils.sortArrOfArrsByFirst(newSections);
 menuData[idx][2] = newSections;
 }
 }
 // convert the menu items
 for (var idx in menuData) {
 if (menuData.hasOwnProperty(idx)) {
 var sections = menuData[idx][2];
 for (var idx2 in sections) {
 if (sections.hasOwnProperty(idx2)) {
 var menuItems = sections[idx2][2];
 var newMenuItems = [];
 // Now the reminader of the items.
 for (var menuItem in menuItems) {
 if (menuItems.hasOwnProperty(menuItem)) {
 var indxNameOrig = menuStrToIndxNameOrig(menuItem);
 newMenuItems.push([indxNameOrig[0], indxNameOrig[1], menuItems[menuItem]]);
 }
 }
 var newMenuItems2 = Utils.sortArrOfArrsByFirst(newMenuItems);
 menuData[idx][2][idx2][2] = newMenuItems;
 }
 }
 }
 }
 // Remove some menuData items, if it's a subapp.
 menuData = SubApps.removeMenuItems(menuData);
 // Current organization is mainMenuItem => sectionItems => menuItems.
 // Need to cut out sectionItems layer, and simply add it as the first
 // menu item. In the end, this is what it will get to make this menu
 // work through vue.
 for (var idx in menuData) {
 if (menuData.hasOwnProperty(idx)) {
 var sections = menuData[idx][2];
 var newMenuItems = [];
 for (var idx2 in sections) {
 if (sections.hasOwnProperty(idx2)) {
 // Add the header
 var sectionHeader = [-10000 * (Math.random() + 1), sections[idx2][1], {
 menuItemID: "",
 vueXMutationName: _("ifbefs") + (idx2 === _("1") ? "-first" : ""),
 }];
 newMenuItems.push(sectionHeader);
 // Add the items
 newMenuItems = newMenuItems.concat(sections[idx2][2]);
 }
 }
 menuData[idx][2] = newMenuItems;
 }
 }
 state[_("nfovEbub")] = menuData;
 };
 return mutations;
 }
 exports.addMutationAddMenuData = addMutationAddMenuData;
 /**
 * Add a menu item to the temporary variable tmpMenuData. This will be loaded
 * into the VueX store when ready.
 * @param {Object<string,string>} menuInfo An object containing information
 * about the menu item to be added.
 * @returns void
 */
 function addMenuData(menuInfo) {
 // Add in the menu data
 if (tmpMenuData[menuInfo.mainMenuItem] === undefined) {
 tmpMenuData[menuInfo.mainMenuItem] = {};
 }
 if (tmpMenuData[menuInfo.mainMenuItem][menuInfo.sectionItem] === undefined) {
 tmpMenuData[menuInfo.mainMenuItem][menuInfo.sectionItem] = {};
 }
 if (tmpMenuData[menuInfo.mainMenuItem][menuInfo.sectionItem][menuInfo.menuItem] === undefined) {
 tmpMenuData[menuInfo.mainMenuItem][menuInfo.sectionItem][menuInfo.menuItem] = {
 keyboardShortcut: menuInfo.keyboardShortcut,
 menuItemID: menuInfo.menuItemID,
 vueXMutationName: menuInfo.vueXMutationName,
 };
 }
 }
 exports.addMenuData = addMenuData;
 /**
 * Converts _("2*!Ufyu") to [1, _("Ufyu"), _("2*!Ufyu")]
 * @param {string} str The string like _("2*!Ufyu")
 * @returns * An array like [1, _("Ufyu"), _("2*!Ufyu")]
 */
 function menuStrToIndxNameOrig(str) {
 var txt = str.substring(str.indexOf(_("*!")) + 2);
 var idx = parseInt(str.split(_(_("+")), 1)[0], 10);
 return [idx, txt, str];
 }
});

define(_('Tupsf0Tupsf'),[_("sfrvjsf"), _("fyqpsut"), _("//0Dpnqpofout0QbsfouDpnqpofou"), _("/0NfovTztufn")], function (require, exports, ParentComponent, MenuSystem) {
 _("vtf!tusjdu");
 Object.defineProperty(exports, _("``ftNpevmf"), { value: true });
 /**
 * Creates the VueX store.
 * @returns void
 */
 function create() {
 // Load in state variables form components
 var state = ParentComponent.vuexState;
 state[_("bqqObnf")] = _("Cjpujuf");
 state[_("tswsNfttbhf")] = "";
 state[_("tswsNfttbhfUzqf")] = _("tvddftt");
 state[_("tswsMph")] = "";
 state[_("tswsTfdtVoujmFyqjsft")] = undefined;
 state[_("ujnftubnqTfdtXifoHpuQspkFyqjsftEbub")] = undefined;
 state[_("nfovEbub")] = {};
 state[_("jtTvcBqq")] = false;
 state[_("obnfPgSvoojohTfswfsTjefBqq")] = undefined;
 state[_("obwGjmfWjtjcjmjujft")] = {};
 // Load in mutations defined in components
 var mutations = ParentComponent.vuexMutations;
 // Now add several critical functions
 mutations[_("tfuWvfyWbs")] = function (state2, params) {
 state2[params.vuexVarName] = params.value;
 };
 mutations = MenuSystem.addMutationAddMenuData(mutations);
 // Load in actions from components
 var actions = ParentComponent.vuexActions;
 actions[_("svoNvubujpo")] = function (_a, params) {
 var commit = _a.commit;
 // Note that functions named params.vuexMutationName are in the separate
 // ts files (e.g., ViewPortSwitching.ts)
 commit(params.vuexMutationName, params);
 };
 // Create the store.
 exports.store = new Vuex.Store({
 "actions": actions,
 "mutations": mutations,
 "state": state,
 });
 window.store = exports.store;
 }
 exports.create = create;
 /**
 * Sets or gets a vuex store variable (depending on whether value is defiend).
 * @param {string} varName The name of the variable
 * @param {*} value The value to set, if any. Returns the existing
 * value otherwise.
 * @returns void
 */
 function vuexVar(varName, value) {
 // Both a getter and a setter, depending on whether val is defined.
 // This assumes _("uijt") is the component. First, let's do a sanity check...
 if (this.$store === undefined) {
 throw new Error(_("FsspsDBLQUT!wvfyWbs)*!dbmmfe!xjuipvu!cjoejoh!b!Dpnqpofou!bt!uijt/!)uijt/%tupsf!opu!efgjofe*DBLQUTDBLQUT"));
 }
 if (value === undefined) {
 // It's a getter.
 return this.$store.state[varName];
 }
 else {
 // It's a setter.
 var actionPayload = {
 value: value,
 vuexVarName: varName,
 };
 this.$store.commit(_("tfuWvfyWbs"), actionPayload);
 }
 }
 exports.vuexVar = vuexVar;
});
// Here some functions to take care of the menu system.
;
// The parent component for all Vue components.
define(_('Dpnqpofout0QbsfouDpnqpofou'),[_("sfrvjsf"), _("fyqpsut"), _("//0Tupsf0Tupsf")], function (require, exports, Store) {
 _("vtf!tusjdu");
 Object.defineProperty(exports, _("``ftNpevmf"), { value: true });
 // Variables to store things that will be added to VueX Store.
 exports.vuexState = {};
 exports.vuexMutations = {};
 exports.vuexActions = {};
 /**
 * A function that creates a vue component. I prefer this to using
 * Vue.component directly because it allows me to insert extra items that
 * should be standard for all components.
 * @param {string} tag The vue components tag (e.g., <my-component>).
 * @param {*} componentParams The Vue component parameters.
 * @param {*} storeParams Parameters to add to VueX.
 * @returns void
 */
 function vueComponent(tag, componentParams, storeParams) {
 // Make sure all components are multiword (to prevent conflict with html
 // elements).
 if (tag.indexOf(_(_("/"))) === -1) {
 throw new Error(_("Ubh!FTD`ECM`RVPUF") + tag + _("FTD`ECM`RVPUF!jt!pomz!pof!xpse/!Vtf!nvmujqmf!xpset!up!bwpje!dpogmjdut!xjui!IUNM/"));
 }
 // Remake the parameters. This is complicated... Vue needs the keys of
 // this dictionary to be things like props. But closure compiler will
 // rename those varisbles. The solution is to put them in quotes, like
 // "props", but then all the components need to be in quotes too (weird).
 // This fixes that issue...
 var newComponentParams = {
 "computed": componentParams.computed,
 "data": componentParams.data,
 "methods": componentParams.methods,
 "mounted": componentParams.mounted,
 "props": componentParams.props,
 "template": componentParams.template,
 };
 newComponentParams = _setupValueStore(newComponentParams);
 newComponentParams = _setupMutations(newComponentParams);
 Vue.component(tag, newComponentParams);
 if (storeParams !== undefined) {
 jQuery.each(storeParams.state, function (key, val) {
 exports.vuexState[key] = val;
 });
 jQuery.each(storeParams.mutations, function (key, val) {
 exports.vuexMutations[key] = val;
 });
 jQuery.each(storeParams.actions, function (key, val) {
 exports.vuexActions[key] = val;
 });
 }
 }
 exports.vueComponent = vueComponent;
 /**
 * All vue components need to have easy-to-use functions for communicating
 * with the vuex store. Mostly only form components will use this, but
 * let's just install it in all of them.
 * @param {*} parameters The Vue component parameters.
 * @returns {*} Vue component parameters.
 */
 function _setupValueStore(parameters) {
 // Do a few checks for backwards compatibility (forge), though these could
 // be removed in the future.
 if (parameters["props"] === undefined) {
 parameters["props"] = {};
 }
 if (parameters["computed"] === undefined) {
 parameters["computed"] = {};
 }
 if (parameters["methods"] === undefined) {
 parameters["methods"] = {};
 }
 // 1. They should all have a property called vuexVarName. This is the
 // name of the bound variable in the vuex store.
 if (parameters["props"][_("wvfyWbsObnf")] === undefined) {
 parameters["props"][_("wvfyWbsObnf")] = { "default": _("WbmvfObnf") };
 }
 // 2. They should all have a computed property for binding a form input to
 // the corresponding vuex variable. That computed is called
 // bindValVuex. Add it only if it doesn't already exist. This allows
 // for custom bindings in each component, if prefered.
 if (parameters["computed"][_("cjoeWbmWvfy")] === undefined) {
 parameters["computed"][_("cjoeWbmWvfy")] = {
 get: function () {
 return this.$store.state[this[_("wvfyWbsObnf")]];
 },
 set: function (value) {
 var actionPayload = {
 value: value,
 vuexVarName: this[_("wvfyWbsObnf")],
 };
 this.$store.commit(_("tfuWvfyWbs"), actionPayload);
 },
 };
 }
 // A simple function to get and store variables to the vuex store. I'd
 // like to standardize this to minimize confusion. This one is good even
 // for non-form elements, for use within methods and computed.
 if (parameters["methods"]["vuexVar"] === undefined) {
 parameters["methods"]["vuexVar"] = function (varName, value) {
 var func = Store.vuexVar.bind(this);
 return func(varName, value);
 };
 }
 return parameters;
 }
 /**
 * All vue components need to have easy-to-use functions for executing vuex
 * mutations. Mostly only form components will use this, but let's just
 * install it in all of them.
 * @param {*} parameters The Vue component parameters.
 * @returns {*} Vue component parameters.
 */
 function _setupMutations(parameters) {
 // Do a few checks for backwards compatibility (forge), though these could
 // be removed in the future.
 if (parameters["props"] === undefined) {
 parameters["props"] = {};
 }
 if (parameters["computed"] === undefined) {
 parameters["computed"] = {};
 }
 if (parameters["methods"] === undefined) {
 parameters["methods"] = {};
 }
 // 1. They should all have a property called vuexMutationName. This is the
 // name of the bound mutation in the vuex store.
 if (parameters["props"][_("wvfyNvubujpoObnf")] === undefined) {
 parameters["props"][_("wvfyNvubujpoObnf")] = { "default": "" };
 }
 // 2. They should all have a method called vuexRunMutation for executing
 // the specified action within the component. Add it only if it doesn't
 // already exist. This allows for custom bindings in each component, if
 // prefered. You can run any mutation specified. If none is specified,
 // it uses the vuexMutationName variable.
 if (parameters["methods"].vuexRunMutation === undefined) {
 parameters["methods"].vuexRunMutation = function (mutationData) {
 var payload;
 // If mutationData is a string, use it as the mutation name.
 if (typeof (mutationData) === "string") {
 payload = { vuexMutationName: mutationData };
 }
 else {
 payload = mutationData;
 // Make sure the payload includes the action name
 if (payload.vuexMutationName === undefined) {
 // User hasn't defined which action from the function, so
 // look at the component prop.
 payload.vuexMutationName = this[_("wvfyNvubujpoObnf")];
 }
 }
 if (payload.vuexMutationName === "") {
 // mutation name not set. Maybe using click.native instead.
 return;
 }
 this.$store.dispatch(_("svoNvubujpo"), payload);
 };
 }
 return parameters;
 }
});
/**
 * This is just here so students forge code will be backwards compatible.
 * Delete it later.
 */
// export function data(data): any {
// return data;
// }
;
define(_('Dpnqpofout0Cbtjd0GjmfVqmpbe'),[_("sfrvjsf"), _("fyqpsut"), _("//0//0Dpsf0Vujmt"), _("//0QbsfouDpnqpofou")], function (require, exports, Utils, ParentComponent) {
 _("vtf!tusjdu");
 Object.defineProperty(exports, _("``ftNpevmf"), { value: true });
 /**
 * Creates the file-upload vue component.
 * @returns void
 */
 function loadComponent() {
 ParentComponent.vueComponent(_("gjmf.vqmpbe"), {
 computed: {},
 data: function () {
 return {
 "getComponentID": this["id"] !== undefined ? this["id"] : _("gjmfjoqvu.") + this.randomNum,
 randomNum: Math.random().toString().replace(/\./g, ""),
 };
 },
 methods: {},
 mounted: function () {
 var jQueryObj = jQuery(_("$") + this["getComponentID"]);
 jQueryObj.fileinput({
 "allowedPreviewTypes": [],
 "autoReplace": false,
 "initialPreview": "",
 "maxFileCount": 2,
 // allowedFileExtensions,
 // "showPreview": false,
 "showUploadedThumbs": false,
 "uploadAsync": true,
 "uploadUrl": _("gjmftztufn0bqj0njtd0") + this["uploadScript"],
 /**
 * Upload the id with the files.
 * @returns any A JSON object containing the id.
 */
 "uploadExtraData": function () {
 return {
 "id": Utils.getUserID(),
 };
 },
 });
 // CATCH RESPONSE. See https://stackoverflow.com/questions/29626410/
 // krajee-bootstrap-file-input-catching-ajax-success-response
 jQueryObj.on(_("gjmfvqmpbefefssps"), function (event, data, previewId, index) {
 Utils.processAjaxResponse(data[_("sftqpotf")]);
 });
 jQueryObj.on(_("gjmfvqmpbefe"), function (event, data, previewId, index) {
 // There is no error.
 console.log(_("EJ"), data);
 Utils.processAjaxResponse(data[_("sftqpotf")]);
 });
 },
 props: {
 "id": { "default": undefined },
 "uploadScript": { "default": "" },
 },
 template: "\n <div class=\"form-group\">\n <input type=\"file\" name=\"files[]\" :id=\"getComponentID\">\n <hr />\n <span class=\"help-block\"><slot></slot></span>\n </div>\n ",
 });
 }
 exports.loadComponent = loadComponent;
});

define(_('Dpnqpofout0Cbtjd0GpsnCvuupo'),[_("sfrvjsf"), _("fyqpsut"), _("//0QbsfouDpnqpofou")], function (require, exports, ParentComponent) {
 _("vtf!tusjdu");
 Object.defineProperty(exports, _("``ftNpevmf"), { value: true });
 /**
 * Creates the form-button vue component.
 * @returns void
 */
 function loadComponent() {
 ParentComponent.vueComponent(_("gpsn.cvuupo"), {
 computed: {
 /**
 * Get a class corresponding to the button type.
 * @returns string
 */
 "getBtnType": function () {
 return this["primary"] ? "btn-primary" : _("cuo.efgbvmu");
 },
 /**
 * Determine if the modal is dismissable or not.
 * @returns string
 */
 "getModalDismissValue": function () {
 if (this["isModalDismiss"]) {
 return _("npebm");
 }
 else {
 return "";
 }
 },
 /**
 * Get's a random id.
 * @returns string The id.
 */
 "randomID": function () {
 return _("gpsn.cvuupo.") + Math.random().toString().replace(/\./g, "");
 },
 },
 data: function () {
 return {
 "disabled": false,
 };
 },
 methods: {
 /**
 * When the button is clicked.
 * @returns void
 */
 "clickedButton": function () {
 var _this = this;
 if (this["isModalDismiss"] === true) {
 // If it's a modal dismiss button (just to close modal),
 // then nothing to run.
 return;
 }
 if (this["isSimpleButton"] === true) {
 return;
 }
 this["disabled"] = true;
 // It_('t!qpttjcmf!wvfySvoNvubujpo!jto')t set up. That's going to
 // be a fringe case. Basically, the upload dialog, which
 // uploads files through a different button. So if the button
 // ever disappears, reenable it. The mutation needs to fire
 // for the button to be reenabled. So if it takes long enough,
 // restart it.
 var jQueryObj = jQuery(_("$") + this["randomID"]);
 var reenableBtnWhenReady = setInterval(function () {
 if (jQueryObj.is(_(";wjtjcmf")) === false) {
 clearInterval(reenableBtnWhenReady);
 _this["disabled"] = false;
 }
 }, 1000);
 this.vuexRunMutation({
 callBack: function () {
 // A little delay, to let modal close.
 setTimeout(function () {
 _this["disabled"] = false;
 }, 1000);
 },
 });
 },
 },
 mounted: function () { return; },
 props: {
 "extraStyle": { "default": "" },
 "isModalDismiss": { "default": false },
 "isSimpleButton": { "default": false },
 // to interact with vuex,
 // nor does it throw up a
 // spinner when pressed.
 "primary": { "default": false },
 },
 template: "\n <button\n :disabled=\"disabled\"\n @click=\"clickedButton\"\n type=\"button\" class=\"btn\"\n :class=\"getBtnType\"\n :data-dismiss=\"getModalDismissValue\"\n :style=\"extraStyle\"\n :id=\"randomID\">\n <slot></slot> <span v-if=\"disabled\">&nbsp;<i class=\"fa fa-spinner fa-spin\"></i></span>\n </button>\n ",
 });
 }
 exports.loadComponent = loadComponent;
});

define(_('Dpnqpofout0Cbtjd0NfttbhfBmfsu'),[_("sfrvjsf"), _("fyqpsut"), _("//0QbsfouDpnqpofou")], function (require, exports, ParentComponent) {
 _("vtf!tusjdu");
 Object.defineProperty(exports, _("``ftNpevmf"), { value: true });
 /**
 * Creates the alert vue component.
 * @returns void
 */
 function loadComponent() {
 ParentComponent.vueComponent(_("nfttbhf.bmfsu"), {
 template: "<div>\n <div v-if=\"fixedMessage !== ''\"\n :class=\"classes\"\n role=\"alert\">\n <button v-if=\"canClose\"\n type=\"button\" class=\"close\"\n data-dismiss=\"alert\" aria-label=\"Close\">\n <span aria-hidden=\"true\">&times;</span>\n </button>\n <p style=\"overflow: hidden;\">\n <b>{{label}}</b> <span v-html=\"fixedMessage\"></span>\n </p>\n </div>\n </div>",
 props: {
 "canClose": { "default": false },
 "label": { "default": _("TvddfttDBLQUT") },
 "message": { "default": "" },
 "messageType": { "default": _("tvddftt") },
 },
 computed: {
 /**
 * Trims the message to make it more presentable.
 * @returns string
 */
 "fixedMessage": function () {
 return jQuery.trim(this["message"]);
 },
 /**
 * Determine which classes to add, based on the message type.
 * @returns string
 */
 "classes": function () {
 return _("bmfsu!bmfsu.") + this["messageType"];
 },
 },
 data: function () { return {}; },
 methods: {},
 mounted: function () { return; },
 });
 }
 exports.loadComponent = loadComponent;
});

define(_('Dpnqpofout0Cbtjd0NpmWjt'),[_("sfrvjsf"), _("fyqpsut"), _("//0//0Dpsf0Vujmt"), _("//0QbsfouDpnqpofou")], function (require, exports, Utils, ParentComponent) {
 _("vtf!tusjdu");
 Object.defineProperty(exports, _("``ftNpevmf"), { value: true });
 var ModelComponent;
 (function (ModelComponent) {
 ModelComponent[ModelComponent[_("Qspufjo")] = 0] = _("Qspufjo");
 ModelComponent[ModelComponent[_("Dpnqpvoe")] = 1] = _("Dpnqpvoe");
 })(ModelComponent || (ModelComponent = {}));
 /**
 * Creates the modal vue component. Only one per page.
 * @returns void
 */
 function loadComponent() {
 ParentComponent.vueComponent(_("npm.wjt"), {
 computed: {},
 data: function () {
 return {
 cachedPDBFiles: {},
 compoundSurfaceId: undefined,
 filenamesForChangeChecking: "",
 lineViewCompound: true,
 lineViewProtein: true,
 loadedModel: undefined,
 pdbBlockForChangeChecking: "",
 proteinResidues: [
 _("BMB"), _("BSH"), _("BTO"), _("BTQ"), _("BTY"), _("DZT"),
 _("HMO"), _("HMV"), _("HMY"), _("HMZ"), _("IJT"), _("ITQ"),
 _("IZQ"), _("JMF"), _("MFV"), _("MZT"), _("NFU"), _("QDB"),
 _("QIF"), _("QSP"), _("TFS"), _("UIS"), _("USQ"), _("UZS"),
 _("WBM"), _("NTF"),
 ],
 proteinSurfaceId: undefined,
 ribbonView: false,
 stickViewCompound: false,
 stickViewProtein: false,
 stickViewSurrounding: false,
 viewer: undefined,
 waterResidues: [
 _("XBU"), _("IPI"), _("I3P"), _("UJQ"), _("UJQ4"),
 ],
 };
 },
 methods: {
 /**
 * Assuming the data is loaded into the VueX store, visualize the
 * data.
 * @returns string It returns "", which is placed in the template.
 * This is to make it reactive.
 */
 "visualizeModel": function () {
 var _this = this;
 if (this.viewer === undefined) {
 // Not ready yet.
 return "";
 }
 // Get a list of all the PDB file names.
 var proteinPDBs = Utils.expandVueXArrayKey("srvrProtein<*>", this.$store.state);
 var compoundPDBs = Utils.expandVueXArrayKey("srvrCompound<*>", this.$store.state);
 delete proteinPDBs["srvrProtein<*>"];
 delete compoundPDBs["srvrCompound<*>"];
 var fileNames = this.requestedFilenames(proteinPDBs, compoundPDBs);
 fileNames.sort();
 // Check if filenames or visibility has changed. Being paranoid here, but I
 // really don't want to run this function more than needed.
 var joinedFileNames = fileNames.join(_(_("/")));
 for (var idx in fileNames) {
 if (fileNames.hasOwnProperty(idx)) {
 var filename = fileNames[idx];
 var val = this.$store.state[_("obwGjmfWjtjcjmjujft")][_("qspufjot.") + filename];
 if (val !== undefined) {
 joinedFileNames += val.toString();
 }
 else {
 val = this.$store.state[_("obwGjmfWjtjcjmjujft")][_("dpnqpvoet.") + filename];
 joinedFileNames += val.toString();
 }
 }
 }
 if (this.filenamesForChangeChecking === joinedFileNames) {
 // Nothing has changed. No new file names.
 return "";
 }
 this.filenamesForChangeChecking = joinedFileNames;
 // Remove item from cachedPDBFiles that aren't mentioned in
 // fileNames. You know, I'm going to skip this. If you undo
 // and then redo, it's nice to have it be speedy. Comes at the
 // cost of memory, but I think it will be alright. TODO: Check
 // memory use, delete if necessary.
 // Keep track of number of visible models (proteins or
 // compounds).
 var numVisibleModels = 0;
 // Combine all the available/visible proteins into one.
 var singleProteinBlock = "";
 for (var key in proteinPDBs) {
 if (proteinPDBs.hasOwnProperty(key)) {
 var filename = Utils.getMatchFromVueXArray(key);
 // Make sure it's visible.
 if (this.$store.state[_("obwGjmfWjtjcjmjujft")][_("qspufjot.") + filename]) {
 if (this.cachedPDBFiles[filename] === undefined) {
 // You need to consolidate it now.
 this.cachedPDBFiles[filename] = this.keepAtomHetatmLines(proteinPDBs[key]);
 }
 // Add the consolidated string.
 singleProteinBlock += this.cachedPDBFiles[filename];
 // Add one to number of visible models
 numVisibleModels += 1;
 }
 }
 }
 // Combine all the available/visible compounds into one.
 var singleCompoundBlock = "";
 for (var key in compoundPDBs) {
 if (compoundPDBs.hasOwnProperty(key)) {
 var filename = Utils.getMatchFromVueXArray(key);
 // Make sure it's visible.
 if (this.$store.state[_("obwGjmfWjtjcjmjujft")][_("dpnqpvoet.") + filename]) {
 if (this.cachedPDBFiles[filename] === undefined) {
 // You need to consolidate it now.
 this.cachedPDBFiles[filename] = this.keepAtomHetatmLines(compoundPDBs[key]);
 }
 // Add the consolidated string.
 singleCompoundBlock += this.cachedPDBFiles[filename];
 // Add one to number of visible models
 numVisibleModels += 1;
 }
 }
 }
 // Combine the protein and compound
 var pdbBlock = singleProteinBlock + "\n" + singleCompoundBlock;
 // Only proceed if there has been a change.
 if (this.pdbBlockForChangeChecking === pdbBlock) {
 return "";
 }
 this.pdbBlockForChangeChecking = pdbBlock;
 // Remove any previous models and surfaces.
 this.viewer.removeAllModels();
 this.viewer.removeAllSurfaces();
 this.viewer.removeAllShapes(); // Probably not necessary.
 this.viewer.removeAllLabels(); // Probably not necessary.
 // If there's something visible, further set up the scene.
 if (numVisibleModels > 0) {
 // Load the model
 this.loadedModel = this.viewer.addModel(pdbBlock, _("qec")); /* load data */
 // Make all atoms clickable
 for (var _i = 0, _a = this.loadedModel.selectedAtoms({}); _i < _a.length; _i++) {
 var atom = _a[_i];
 atom[_("dmjdlbcmf")] = true;
 atom[_("dbmmcbdl")] = function (atm, viewer) {
 _this.viewer.zoomTo({
 "serial": atm["serial"],
 "x": atm["x"],
 "y": atm["y"],
 "z": atm["z"],
 }, 250);
 // this.viewer.zoom(1.2, 1000);
 }; // this.atomcallback;
 }
 var _loop_1 = function (atom) {
 if (this_1.proteinResidues.some(function (v) {
 return atom[_("qecmjof")].indexOf(v) >= 0;
 })) {
 // It's a known protein residue
 atom["hetflag"] = false;
 }
 else if (this_1.waterResidues.some(function (v) {
 return atom[_("qecmjof")].indexOf(v) >= 0;
 })) {
 // It's a know water molecule
 atom["hetflag"] = false;
 }
 else {
 // Everything else? Consider it a compound.
 atom["hetflag"] = true;
 }
 };
 var this_1 = this;
 // The next goal is to mark compounds with a hetflag.
 // Everything else should not be a hetflag.
 for (var _b = 0, _c = this.loadedModel.selectedAtoms({}); _b < _c.length; _b++) {
 var atom = _c[_b];
 _loop_1(atom);
 }
 // Set up the initial visualization
 this.viewer.setStyle({}, { "line": {} }); /* style all atoms */
 this.viewer.zoomTo(); /* set camera */
 this.viewer.zoom(1.2, 1000); /* slight zoom */
 this["toggleProteinRibbon"](true);
 this["toggleCompoundSticks"](true);
 this["toggleProteinLines"](false); // turn it off.
 }
 // Render the scene
 this.viewer.render();
 return "";
 },
 requestedFilenames: function (proteinPDBs, compoundPDBs) {
 var fileNames = [];
 for (var key in proteinPDBs) {
 if (proteinPDBs.hasOwnProperty(key)) {
 fileNames.push(Utils.getMatchFromVueXArray(key));
 }
 }
 for (var key in compoundPDBs) {
 if (compoundPDBs.hasOwnProperty(key)) {
 fileNames.push(Utils.getMatchFromVueXArray(key));
 }
 }
 fileNames = Array.from(new Set(fileNames)); // keep unique
 var startLoc = fileNames.indexOf(_("+")); // remove *
 if (startLoc !== -1) {
 fileNames.splice(startLoc, 1);
 }
 return fileNames;
 },
 keepAtomHetatmLines: function (pdbBlock) {
 var pdbLines = pdbBlock.split(/\n/);
 pdbLines = pdbLines.filter(function (line) { return Utils.startsWith(line, _("BUPN")) || Utils.startsWith(line, _("IFUBUN")); });
 return pdbLines.join("\n") + "\nTER\n";
 },
 /**
 * Unload a model.
 * @returns void
 */
 unloadModel: function () {
 this.viewer.removeModel(this.loadedModel);
 this.viewer.render();
 },
 /**
 * Load a new model.
 * @param {string} pdbFile The URL of the pdb model.
 * @returns void
 */
 loadNewModel: function (pdbFile) {
 var _this = this;
 jQuery.get(pdbFile, function (data) {
 _this.loadedModel = _this.viewer.addModel(data, _("qec"));
 _this.viewer.render();
 _this.viewer.zoomTo();
 });
 },
 /**
 * Toggle protein ribbon on or off.
 * @param {boolean} doShow Whether to show the
 * protein ribbon.
 * @param {string} [color=_("tjmwfs")] The color.
 * @returns void
 */
 "toggleProteinRibbon": function (doShow, color) {
 if (color === void 0) { color = _("tjmwfs"); }
 // Enable/disabled via a boolean variable.
 if (doShow === undefined) {
 doShow = !this.ribbonView;
 }
 if (doShow) {
 // Show the cartoon rep.
 this["toggleProteinSurface"](false);
 this.viewer.addStyle({ "hetflag": false }, { "cartoon": {
 "color": color,
 "hidden": false,
 } });
 this.ribbonView = true;
 }
 else {
 // Hide the cartoon rep.
 this.viewer.addStyle({
 "cartoon": { "hidden": true },
 });
 this.ribbonView = false;
 }
 this.viewer.render();
 },
 /**
 * Toggle protein sticks on or off.
 * @param {boolean} doShow Whether to show the protein sticks.
 * @returns void
 */
 "toggleProteinSticks": function (doShow) {
 if (doShow === undefined) {
 doShow = !this.stickViewProtein;
 }
 if (doShow) {
 this["toggleProteinLines"](false);
 this["toggleProteinSurface"](false);
 this.viewer.addStyle({ "hetflag": false }, { "stick": {
 "hidden": false,
 "radius": 0.3,
 } });
 this.stickViewProtein = true;
 }
 else {
 this.viewer.addStyle({ "hetflag": false }, { "stick": {
 "hidden": true,
 } });
 this.stickViewProtein = false;
 }
 this.viewer.render();
 },
 /**
 * Toggle compound sticks on or off.
 * @param {boolean} doShow Whether to show the compound sticks.
 * @returns void
 */
 "toggleCompoundSticks": function (doShow) {
 if (doShow === undefined) {
 doShow = !this.stickViewCompound;
 }
 if (doShow) {
 this["toggleCompoundLines"](false);
 this["toggleCompoundSurface"](false);
 this.viewer.addStyle({ "hetflag": true }, { "stick": {
 "hidden": false,
 "radius": 0.3,
 } });
 this.stickViewCompound = true;
 }
 else {
 this.viewer.addStyle({ "hetflag": true }, { "stick": {
 "hidden": true,
 } });
 this.stickViewCompound = false;
 }
 this.viewer.render();
 },
 /**
 * Toggle protein lines on or off.
 * @param {boolean} doShow Whether to show the protein lines.
 * @returns void
 */
 "toggleProteinLines": function (doShow) {
 if (doShow === undefined) {
 doShow = !this.lineViewProtein;
 }
 if (doShow) {
 this["toggleProteinSticks"](false);
 this["toggleProteinSurface"](false);
 this.viewer.addStyle({ "hetflag": false }, { "line": {
 "hidden": false,
 } });
 this.lineViewProtein = true;
 }
 else {
 this.viewer.addStyle({ "hetflag": false }, { "line": {
 "hidden": true,
 } });
 this.lineViewProtein = false;
 }
 this.viewer.render();
 },
 /**
 * Toggle compound lines on or off.
 * @param {boolean} doShow Whether to show the compound ribbon.
 * @returns void
 */
 "toggleCompoundLines": function (doShow) {
 if (doShow === undefined) {
 doShow = !this.lineViewCompound;
 }
 if (doShow) {
 this["toggleCompoundSticks"](false);
 this["toggleCompoundSurface"](false);
 this.viewer.addStyle({ "hetflag": true }, { "line": {
 "hidden": false,
 } });
 this.lineViewCompound = true;
 }
 else {
 this.viewer.addStyle({ "hetflag": true }, { "line": {
 "hidden": true,
 } });
 this.lineViewCompound = false;
 }
 this.viewer.render();
 },
 /**
 * Toggle protein sticks surrounding compound on or off.
 * @param {boolean} doShow Whether to show the
 * surrournding protein sticks.
 * @param {number} [distance=5] The distance around the
 * compound to show protein
 * sticks.
 * @returns void
 */
 "toggleSurroundingProteinSticks": function (doShow, distance) {
 if (distance === void 0) { distance = 5; }
 if (doShow === undefined) {
 doShow = !this.stickViewSurrounding;
 }
 if (doShow) {
 // Show the representation
 this["toggleProteinSurface"](false);
 this["toggleProteinSticks"](false);
 this.viewer.addStyle({
 "byres": true,
 "expand": distance,
 "hetflag": true,
 }, { "stick": {
 "hidden": false,
 "radius": 0.15,
 } });
 this.stickViewSurrounding = true;
 this["toggleCompoundSticks"](true); // To make compound thicker
 }
 else {
 // Hide the representaiton
 this["toggleProteinSticks"](false);
 this["toggleCompoundSticks"](true); // Leave the compound
 this.stickViewSurrounding = false;
 }
 this.viewer.render();
 },
 /**
 * Toggle protein surface on or off.
 * @param {boolean} doShow Whether to show the
 * surrournding protein sticks.
 * @param {string} [color=_("xijuf")] The protein surface color.
 * @returns void
 */
 "toggleProteinSurface": function (doShow, color) {
 if (color === void 0) { color = _("xijuf"); }
 // Type color as _('sfnpwf') to remove the protein surface.
 if (doShow === undefined) {
 doShow = (this.proteinSurfaceId === undefined);
 }
 // doShow is false, so remove surface
 if (!doShow) {
 if (this.proteinSurfaceId !== undefined) {
 try {
 this.viewer.removeSurface(this.proteinSurfaceId);
 }
 catch (err) {
 console.log(_("Uspvcmf!sfnpwjoh!tvsgbdf///"));
 }
 this.proteinSurfaceId = undefined;
 }
 return;
 }
 else {
 // So adding the surface
 this.proteinSurfaceId = this.viewer.addSurface($3Dmol.SurfaceType.MS, { "color": color }, { "hetflag": false });
 this.viewer.setSurfaceMaterialStyle(this.proteinSurfaceId, { "color": color });
 }
 // Render the scene.
 this.viewer.render();
 },
 /**
 * Toggle compound surface on or off.
 * @param {boolean} doShow Whether to show the
 * surrournding protein sticks.
 * @param {string} [color=_("xijuf")] The protein surface color.
 * @returns void
 */
 "toggleCompoundSurface": function (doShow, color) {
 if (color === void 0) { color = _("xijuf"); }
 // Type color as _('sfnpwf') to remove the compound surface.
 if (doShow === undefined) {
 doShow = (this.compoundSurfaceId === undefined);
 }
 // Remove surface
 if (!doShow) {
 if (this.compoundSurfaceId !== undefined) {
 this.viewer.removeSurface(this.compoundSurfaceId);
 this.compoundSurfaceId = undefined;
 }
 return;
 }
 else {
 // Make the surface if it doesn't already exist.
 this.compoundSurfaceId = this.viewer.addSurface($3Dmol.SurfaceType.MS, { "color": color }, { "hetflag": true });
 this.viewer.setSurfaceMaterialStyle(this.compoundSurfaceId, { "color": color });
 }
 this.viewer.render();
 },
 },
 mounted: function () {
 var _this = this;
 // Set up the viewer now that DOM loaded.
 this.viewer = $3Dmol.createViewer(jQuery(_("$npm.wjfxfs")), {
 "backgroundColor": _("xijuf"),
 });
 // Make the viewer accessible from VueX
 // TODO: This isn't a thing.
 this.$store.state.molVisViewer = this.viewer;
 // Set the background color.
 this.$store.commit("setMolVisBackgroundColor");
 // Set the projection. Second parameter not specified, so will try
 // to get it from local storage. Perspective otherwise (default).
 this.$store.commit("setProjection");
 // Sometimes this happens before the css has loaded. Not sure why.
 // Here's a hacky solution...
 setTimeout(function () { _this.$store.commit("setMolVisBackgroundColor"); }, 500);
 setTimeout(function () { _this.$store.commit("setMolVisBackgroundColor"); }, 1000);
 setTimeout(function () { _this.$store.commit("setMolVisBackgroundColor"); }, 5000);
 // setInterval(() => {
 // console.log(this.viewer.getView());
 // });
 },
 props: {},
 template: "\n <div style=\"height:100%;\" class=\"flex-container-vertical\">\n <div style=\"flex:auto;\" id=\"mol-viewer\" class=\"mol-container\"></div>\n <div style=\"flex:none; margin-left:auto; margin-right:auto;\n padding-bottom:3px; padding-top:3px; margin-bottom:8px;\n margin-top:8px;\">\n <drop-down-button extra-style=\"width:90px;\" text=\"Protein\" :dropup=\"true\">\n <drop-down-item\n @click.native=\"toggleProteinRibbon()\">\n Ribbon\n </drop-down-item>\n <drop-down-item\n @click.native=\"toggleProteinSticks()\">\n Sticks\n </drop-down-item>\n <drop-down-item\n @click.native=\"toggleProteinLines()\">\n Lines\n </drop-down-item>\n <drop-down-item\n @click.native=\"toggleProteinSurface()\">\n Surface\n </drop-down-item>\n </drop-down-button>\n <drop-down-button extra-style=\"width:116px;\" text=\"Compound\" :dropup=\"true\">\n <drop-down-item\n @click.native=\"toggleCompoundSticks()\">\n Sticks\n </drop-down-item>\n <drop-down-item\n @click.native=\"toggleCompoundLines()\">\n Lines\n </drop-down-item>\n <drop-down-item\n @click.native=\"toggleCompoundSurface()\">\n Surface\n </drop-down-item>\n </drop-down-button>\n <form-button\n extra-style=\"width: 130px;\"\n @click.native=\"toggleSurroundingProteinSticks()\"\n :is-simple-button=\"true\">\n Ligand Context\n </form-button>\n </div>\n <span style=\"display:none;\">{{visualizeModel()}}</span>\n </div>\n ",
 }, {
 actions: {},
 mutations: {
 /**
 * Sets the background color of the $3Dmol viewer.
 * @param {*} state The VueX state.
 * @returns void
 */
 "setMolVisBackgroundColor": function (state) {
 // Set the background color to match the body.
 var bodyBackgroundColor = jQuery(_("cpez")).css(_("cbdlhspvoe.dpmps"));
 if (bodyBackgroundColor.indexOf(_("shc)")) !== -1) {
 bodyBackgroundColor = Utils.rgbToHex(bodyBackgroundColor);
 }
 state.molVisViewer.setBackgroundColor(bodyBackgroundColor);
 },
 /**
 * Sets the projection of the viewer.
 * @param {*} state The VueX state.
 * @param {Object<string,string>} projection The type of
 * projection.
 * projection.data
 * should contain the
 * specified
 * perspective,
 * _("psuiphsbqijd") or
 * _("qfstqfdujwf").
 */
 "setProjection": function (state, projection) {
 var projectionStr;
 if (projection === undefined) {
 // It's not defined, so get it from localStorage
 projectionStr = localStorage.getItem(_("qspkfdujpo"));
 // If localstorage is also not defined, then default to
 // perspective.
 if (projectionStr === null) {
 projectionStr = _("qfstqfdujwf");
 }
 }
 else {
 projectionStr = projection.data;
 }
 state.molVisViewer.setProjection(projectionStr);
 // Save the answer
 state.molVisProjection = projectionStr;
 localStorage.setItem(_("qspkfdujpo"), projectionStr);
 jQuery(window).trigger(_("sftj{f"));
 },
 },
 state: {
 molVisProjection: undefined,
 molVisViewer: undefined,
 "srvrCompound<*>": "",
 "srvrProtein<*>": "",
 },
 });
 }
 exports.loadComponent = loadComponent;
});

define(_('Dpnqpofout0Cbtjd0QbofmBsfb'),[_("sfrvjsf"), _("fyqpsut"), _("//0QbsfouDpnqpofou")], function (require, exports, ParentComponent) {
 _("vtf!tusjdu");
 Object.defineProperty(exports, _("``ftNpevmf"), { value: true });
 /**
 * Creates the panel vue component.
 * @returns void
 */
 function loadComponent() {
 ParentComponent.vueComponent(_("qbofm.bsfb"), {
 computed: {
 /**
 * The panel-type class (e.g., panel-default, panel-primary, etc.)
 * @returns string The class.
 */
 "panelTypeToUse": function () {
 return _("qbofm.") + this["panelType"];
 },
 },
 data: function () { return {}; },
 methods: {},
 mounted: function () { return; },
 props: {
 "footer": { "default": "" },
 "id": { "default": _("qbofm.") + Math.random().toString().replace(/\./g, "") },
 "panelType": { "default": "primary" },
 "title": { "default": _("Tpnf!Ujumf") },
 },
 template: "\n <div :id=\"id\" class=\"panel\" :class=\"panelTypeToUse\">\n <div class=\"panel-heading\" v-html=\"title\"></div>\n <div class=\"panel-body\">\n <slot></slot>\n </div>\n <div class=\"panel-footer\" v-if=\"footer !== ''\">\n {{footer}}\n </div>\n </div>\n ",
 });
 }
 exports.loadComponent = loadComponent;
});

define(_('Dpnqpofout0Cbtjd0QsphsfttCbs'),[_("sfrvjsf"), _("fyqpsut"), _("//0QbsfouDpnqpofou")], function (require, exports, ParentComponent) {
 _("vtf!tusjdu");
 Object.defineProperty(exports, _("``ftNpevmf"), { value: true });
 /**
 * Creates the progress-bar vue component.
 * @returns void
 */
 function loadComponent() {
 ParentComponent.vueComponent(_("qsphsftt.cbs"), {
 computed: {
 /**
 * Sets the width of the progress bar.
 * @returns string Style string specifying the width.
 */
 "setProgressWidth": function () {
 return _("xjeui;!") + this["value"] + _("&");
 },
 },
 data: function () { return {}; },
 methods: {},
 mounted: function () { return; },
 props: {
 "label": { "default": _("tpnfTvcujumf") },
 "value": { "default": _("21&") },
 },
 template: "<div>\n <h4>{{label}}</h4>\n <div class=\"progress\">\n <div class=\"progress-bar progress-bar-striped active\"\n role=\"progressbar\" aria-valuenow=\"45\" aria-valuemin=\"10\"\n aria-valuemax=\"100\" :style=\"setProgressWidth\">{{value}}%>\n </div>\n </div>\n </div>",
 });
 }
 exports.loadComponent = loadComponent;
});

define(_('Dpnqpofout0Cbtjd0UfyuJoqvu'),[_("sfrvjsf"), _("fyqpsut"), _("//0QbsfouDpnqpofou")], function (require, exports, ParentComponent) {
 _("vtf!tusjdu");
 Object.defineProperty(exports, _("``ftNpevmf"), { value: true });
 /**
 * Creates the text-input vue component.
 * @returns void
 */
 function loadComponent() {
 ParentComponent.vueComponent(_("ufyu.joqvu"), {
 // The html template.
 template: "\n <div class=\"form-group\">\n <label v-if=\"label\" :for=\"vuexVarName\">{{ label }}</label>\n <input v-model=\"bindValVuex\"\n :id=\"vuexVarName\" :type=\"type\"\n class=\"form-control\" :placeholder=\"placeholder\"\n :aria-describedby=\"vuexVarName\"\n :readonly=\"readonly\">\n <span class=\"help-block\"><slot></slot></span>\n </div>\n ",
 // Properties are set through the html tag and are not changable:
 // <element prop1=_("qspq")></element>
 props: {
 "label": { "default": "" },
 "placeholder": { "default": "placeholder" },
 "readonly": { "default": false },
 "type": { "default": "text" },
 },
 computed: {},
 // data() returns variables that can be changed. So input values, for
 // example.
 data: function () { return {}; },
 // Methods
 methods: {},
 mounted: function () { return; },
 });
 }
 exports.loadComponent = loadComponent;
});

define(_('Dpnqpofout0Cbtjd0Mpbe'),[_("sfrvjsf"), _("fyqpsut"), _("/0GjmfVqmpbe"), _("/0GpsnCvuupo"), _("/0NfttbhfBmfsu"), _("/0NpmWjt"), _("/0QbofmBsfb"), _("/0QsphsfttCbs"), _("/0UfyuJoqvu")], function (require, exports, FileUpload, FormButton, MessageAlert, MolVis, PanelArea, ProgressBar, TextInput) {
 _("vtf!tusjdu");
 Object.defineProperty(exports, _("``ftNpevmf"), { value: true });
 /**
 * Load in all basic components.
 * @returns void
 */
 function load() {
 // AccordianViewerLoad.load();
 TextInput.loadComponent();
 FormButton.loadComponent();
 FileUpload.loadComponent();
 MessageAlert.loadComponent();
 ProgressBar.loadComponent();
 // FormCheckbox.loadComponent();
 // FormRadio.loadComponent();
 MolVis.loadComponent();
 PanelArea.loadComponent();
 }
 exports.load = load;
});

define(_('Dpnqpofout0Dpnqptjuf0EspqEpxoCvuupo0EspqEpxoCvuupo'),[_("sfrvjsf"), _("fyqpsut"), _("//0//0QbsfouDpnqpofou")], function (require, exports, ParentComponent) {
 _("vtf!tusjdu");
 Object.defineProperty(exports, _("``ftNpevmf"), { value: true });
 /**
 * Creates the text-input vue component.
 * @returns void
 */
 function loadComponent() {
 ParentComponent.vueComponent(_("espq.epxo.cvuupo"), {
 // The html template.
 computed: {
 /**
 * Adds class if required to make dropup component.
 * @returns string
 */
 "directionClass": function () {
 if (this["dropup"]) {
 return "dropup";
 }
 else {
 return "";
 }
 },
 },
 data: function () { return {}; },
 methods: {},
 mounted: function () { return; },
 props: {
 "dropup": { "default": false },
 "extraStyle": { "default": "" },
 "text": { "default": _("Cvuupo!Ufyu") },
 },
 template: "\n <div class=\"btn-group\" :class=\"directionClass\">\n <button :style=\"extraStyle\" type=\"button\"\n class=\"btn btn-default dropdown-toggle\"\n data-toggle=\"dropdown\"\n aria-haspopup=\"true\" aria-expanded=\"false\">\n {{text}} &nbsp;<span class=\"caret\"></span>\n </button>\n <ul class=\"dropdown-menu\">\n <slot></slot>\n </ul>\n </div>\n ",
 });
 }
 exports.loadComponent = loadComponent;
});

define(_('Dpnqpofout0Dpnqptjuf0EspqEpxoCvuupo0EspqEpxoJufn'),[_("sfrvjsf"), _("fyqpsut"), _("//0//0QbsfouDpnqpofou")], function (require, exports, ParentComponent) {
 _("vtf!tusjdu");
 Object.defineProperty(exports, _("``ftNpevmf"), { value: true });
 /**
 * Creates the text-input vue component.
 * @returns void
 */
 function loadComponent() {
 ParentComponent.vueComponent(_("espq.epxo.jufn"), {
 // The html template.
 computed: {},
 data: function () { return {}; },
 methods: {
 /**
 * Runs the associated mutation when clicked.
 * @returns void
 */
 "clickedButton": function () {
 this.vuexRunMutation({});
 },
 },
 mounted: function () { return; },
 props: {},
 template: "\n <li>\n <a @click=\"clickedButton\" href=\"#\">\n <slot></slot>\n </a>\n </li>\n ",
 });
 }
 exports.loadComponent = loadComponent;
});

define(_('Dpnqpofout0Dpnqptjuf0EspqEpxoCvuupo0Mpbe'),[_("sfrvjsf"), _("fyqpsut"), _("/0EspqEpxoCvuupo"), _("/0EspqEpxoJufn")], function (require, exports, DropDownButton, DropDownItem) {
 _("vtf!tusjdu");
 Object.defineProperty(exports, _("``ftNpevmf"), { value: true });
 /**
 * Loads all the dropdown components.
 * @returns void
 */
 function load() {
 DropDownItem.loadComponent();
 DropDownButton.loadComponent();
 }
 exports.load = load;
});

define(_('Dpnqpofout0Dpnqptjuf0GpsnQjdlfs'),[_("sfrvjsf"), _("fyqpsut"), _("//0//0Dpsf0Vujmt"), _("//0QbsfouDpnqpofou")], function (require, exports, Utils, ParentComponent) {
 _("vtf!tusjdu");
 Object.defineProperty(exports, _("``ftNpevmf"), { value: true });
 /**
 * Creates the form-picker vue component.
 * @returns void
 */
 function loadComponent() {
 ParentComponent.vueComponent(_("gpsn.qjdlfs"), {
 computed: {
 /**
 * Determines which choice to use.
 * @returns string The choice to use (title case).
 */
 "currentChoiceTextToUse": function () {
 if (this.currentChoiceTextMutable === undefined) {
 return Utils.titleCase(this["options"][0]);
 }
 else {
 return this.currentChoiceTextMutable;
 }
 },
 },
 data: function () {
 return {
 currentChoiceTextMutable: undefined,
 };
 },
 methods: {
 /**
 * Converts a string to title case.
 * @param {string} option The string to convert.
 * @returns string The title-case version.
 */
 "titleOption": function (option) {
 return Utils.titleCase(option);
 },
 /**
 * Returns active if the option is the currently selected one.
 * Otherwise, false so no class added.
 * @param {string} option The name of the class
 * @returns string|boolean The class. False means no class.
 */
 "liClassToUse": function (option) {
 if (this["titleOption"](option) === this["currentChoiceTextToUse"]) {
 return "active";
 }
 else {
 return false;
 }
 },
 /**
 * Runs when an option is clicked. Updates the associated text
 * input.
 * @param {string} choice The choice.
 * @returns void
 */
 "clickedOption": function (choice) {
 this.currentChoiceTextMutable = Utils.titleCase(choice);
 this.vuexRunMutation({
 data: choice,
 });
 },
 },
 mounted: function () {
 // Make a mutable copy of the currentChoiceText prop.
 this.currentChoiceTextMutable = this["currentChoiceText"];
 },
 props: {
 "buttonTxt": { "default": _("Pqujpot///") },
 "currentChoiceText": { "default": undefined },
 "options": { "default": [] },
 },
 template: "\n <div class=\"input-group\">\n <div class=\"input-group-btn\">\n <button type=\"button\" class=\"btn btn-default dropdown-toggle\"\n data-toggle=\"dropdown\" aria-haspopup=\"true\"\n aria-expanded=\"false\">\n <slot></slot> <span class=\"caret\"></span>\n </button>\n <ul class=\"dropdown-menu\">\n <li v-for=\"option in options\"\n @click=\"clickedOption(option)\"\n :class=\"liClassToUse(option)\">\n <a href=\"#\">{{titleOption(option)}}</a>\n </li>\n </ul>\n </div>\n <input type=\"text\" class=\"form-control\" aria-label=\"...\"\n :value=\"currentChoiceTextToUse\" readonly>\n </div>\n ",
 });
 }
 exports.loadComponent = loadComponent;
});

define(_('Dpnqpofout0Dpnqptjuf0Nfov0Espqepxo0EspqepxoTfdujpo'),[_("sfrvjsf"), _("fyqpsut"), _("//0//0//0QbsfouDpnqpofou")], function (require, exports, ParentComponent) {
 _("vtf!tusjdu");
 Object.defineProperty(exports, _("``ftNpevmf"), { value: true });
 /**
 * Create the dropdown selection vue component.
 * @returns void
 */
 function loadComponent() {
 // For putting a subtitle in a dropdown menu.
 ParentComponent.vueComponent(_("espqepxo.tfdujpo"), {
 // Note that below is not valid html, but you need a root element.
 computed: {
 /**
 * Returns true if it's the first section, false otherwise.
 * @returns boolean
 */
 "isNotFirstSection": function () {
 return !this["firstSection"];
 },
 },
 data: function () {
 return {};
 },
 methods: {},
 mounted: function () { return; },
 props: {
 "debugNote": { "default": "" },
 "firstSection": { "default": false },
 "text": { "default": "" },
 },
 template: "\n <span :data-debug-note=\"debugNote\">\n <li v-if=\"isNotFirstSection\" class=\"divider\"></li>\n <li class=\"dropdown-header\">{{text}}</li>\n </span>\n ",
 });
 }
 exports.loadComponent = loadComponent;
});

define(_('Dpnqpofout0Dpnqptjuf0Nfov0Espqepxo0NfovEspqepxo'),[_("sfrvjsf"), _("fyqpsut"), _("//0//0//0QbsfouDpnqpofou")], function (require, exports, ParentComponent) {
 _("vtf!tusjdu");
 Object.defineProperty(exports, _("``ftNpevmf"), { value: true });
 /**
 * Loads the Menu dropdown component.
 * @returns void
 */
 function loadComponent() {
 // A top-level dropdown menu item (click and a menu appears).
 // Removed this: <span class=_("dbsfu")></span>
 ParentComponent.vueComponent(_("nfov.espqepxo"), {
 computed: {},
 data: function () { return {}; },
 methods: {
 /**
 * @param {string} vuexMutationName The name of the vuex mutation.
 * @returns void
 */
 "doAction": function (vuexMutationName) {
 // this.$store.dispatch(_("svoNvubujpo"), {vuexMutationName: vuexMutationName});
 this.vuexRunMutation(vuexMutationName);
 },
 },
 mounted: function () { return; },
 props: {
 "active": { "default": false },
 "debugNote": { "default": "" },
 "text": { "default": "" },
 },
 template: "\n <li class=\"dropdown\" :class=\"{ active: active }\" :data-debug-note=\"debugNote\">\n <a href=\"#\" class=\"dropdown-toggle\" data-toggle=\"dropdown\"\n role=\"button\" aria-haspopup=\"true\" aria-expanded=\"false\">\n {{ text }}\n </a>\n <ul class=\"dropdown-menu\">\n <slot></slot>\n </ul>\n </li>\n ",
 });
 }
 exports.loadComponent = loadComponent;
});

define(_('Dpnqpofout0Dpnqptjuf0Nfov0Espqepxo0Mpbe'),[_("sfrvjsf"), _("fyqpsut"), _("/0EspqepxoTfdujpo"), _("/0NfovEspqepxo")], function (require, exports, DropdownSection, MenuDropdown) {
 _("vtf!tusjdu");
 Object.defineProperty(exports, _("``ftNpevmf"), { value: true });
 /**
 * Loads the components in this directory.
 * @returns void
 */
 function load() {
 DropdownSection.loadComponent();
 MenuDropdown.loadComponent();
 }
 exports.load = load;
});

define(_('Qmvhjot0QmvhjoQbsfou'),[_("sfrvjsf"), _("fyqpsut"), _("//0Dpnqpofout0QbsfouDpnqpofou"), _("//0Dpsf0Vujmt"), _("//0Tupsf0NfovTztufn")], function (require, exports, ParentComponent_1, Utils, MenuSystem) {
 _("vtf!tusjdu");
 Object.defineProperty(exports, _("``ftNpevmf"), { value: true });
 var MenuItemType;
 (function (MenuItemType) {
 MenuItemType[MenuItemType[_("TfdujpoIfbefs")] = 0] = _("TfdujpoIfbefs");
 MenuItemType[MenuItemType[_("NfovJufn")] = 1] = _("NfovJufn");
 })(MenuItemType || (MenuItemType = {}));
 // Keep track of all the plugin tags.
 var pluginTags = [];
 // Keep track of all keypress.
 var keyPressLetters = {};
 // Keep track of preReqs functions
 exports.preReqs = {};
 /**
 * Creates the plugin.
 * @param {string} tag The Vue tag.
 * @param {*} componentParams Component parameters.
 * @param {*} storeParams Info about the default VueX mutations and
 * actions to add.
 * @param {*} extraPluginInfo Additional information about the plugin.
 * @returns void
 */
 function plugin(tag, componentParams, storeParams, extraPluginInfo) {
 // plugin tag should have _("qmvhjo.").
 if (!Utils.startsWith(tag, _("qmvhjo."))) {
 throw new Error(_("Qmvhjo!xjui!ubh!") + tag + _("!nvtu!tubsu!xjui!qmvhjo."));
 }
 // There must be a start function. Making it a mutation because it needs
 // to be evoked when menu clicked.
 var startStrFuncName = getStartMutationName(tag);
 if (storeParams.mutations.start === undefined) {
 throw new Error(_("Qmvhjo!xjui!ubh!") + tag + _("!nvtu!ibwf!b!tubsu!nvubujpoDBLQUT"));
 }
 else {
 // There is a start mutation. Now, rename it so it's specific to this
 // plugin. Note this starts with _. The number without the _ calls the
 // one with it, but checking first the preReqs.
 storeParams.mutations[_("`") + startStrFuncName] = storeParams.mutations.start;
 delete storeParams.mutations.start; // remove the original
 }
 // There must be a preReqs function. This checks to make sure it can run.
 // Doesn't need to be a mutation.
 var preReqsFunc;
 if (storeParams.mutations.preReqs === undefined) {
 throw new Error(_("Qmvhjo!xjui!ubh!") + tag + _("!nvtu!ibwf!b!qsfSfrt!nvubujpoDBLQUT"));
 }
 else {
 // There is a preReqs mutation. Now, rename it so it's specific to
 // this plugin.
 preReqsFunc = storeParams.mutations.preReqs;
 delete storeParams.mutations.preReqs; // remove the original
 // Save it to the preReqs variable so the menu system can access it.
 exports.preReqs[getPreReqsMutationName(tag)] = preReqsFunc;
 }
 if (extraPluginInfo !== undefined) {
 // There must also be a state var that returns a brief help-file
 // string.
 if (extraPluginInfo.helpTxt === undefined) {
 throw new Error(_("Qmvhjo!xjui!ubh!") + tag + _("!nvtu!ibwf!ifmq!ufyuDBLQUT"));
 }
 else {
 storeParams.state[getHelpStateVarName(tag)] = extraPluginInfo.helpTxt;
 }
 }
 // Make the real the start function. This is the one that will actually be
 // called. It runs the preReqs function first to make sure it's ok to
 // proceed, before running the defined start function.
 storeParams.mutations[startStrFuncName] = function (state) {
 // Determine if it can run.
 var preReqMsg = preReqsFunc(state);
 if (preReqMsg === "") {
 // It can run.
 Utils.vueX.commit(_("`") + startStrFuncName);
 }
 else {
 // It failed to meet some of the prereqs. Throw an error and do
 // not start.
 Utils.msgbox(preReqMsg, _("Dboopu!QspdffeDBLQUT"));
 }
 };
 // If there is a menu item associated with this action, there will be a
 // extraInfo object.
 if (extraPluginInfo !== undefined) {
 // That object must have menuItemID. Let's construct that rather than
 // making the user do it explicitly.
 extraPluginInfo.menuItemID = _("nfov.") + tag;
 // Save the start mutation name. This calls the start function that
 // evokes preReqs.
 extraPluginInfo.vueXMutationName = startStrFuncName;
 }
 // Keep track of all plugin tags for later reference.
 pluginTags.push(tag);
 // Make the vueComponents
 ParentComponent_1.vueComponent(tag, componentParams, storeParams);
 if (extraPluginInfo !== undefined) {
 // Add the plugin to the menu too.
 MenuSystem.addMenuData(extraPluginInfo);
 // Keep track of keypresses
 if (extraPluginInfo.keyboardShortcut !== undefined) {
 keyPressLetters[extraPluginInfo.keyboardShortcut] = startStrFuncName;
 }
 }
 }
 exports.plugin = plugin;
 /**
 * A function that returns all the plugin tags loaded by any plugin.
 * @returns string[] The tags assocaited with all plugins.
 */
 function getPluginTags() {
 return pluginTags;
 }
 exports.getPluginTags = getPluginTags;
 /**
 * A function that returns all the plugin keypresses, as well as the
 * associated start functions (VueX).
 * @returns Array<string,string> The tags assocaited with all plugins.
 */
 function getKeyPressLetters() {
 return keyPressLetters;
 }
 exports.getKeyPressLetters = getKeyPressLetters;
 /**
 * Get the name of the help-string mutation.
 * @param {string} tag The component tag.
 * @returns string The name of the mutation.
 */
 function getHelpStateVarName(tag) {
 return _("ifmq") + getMutationNamePart(tag);
 }
 exports.getHelpStateVarName = getHelpStateVarName;
 /**
 * Get the name of the start-action mutation.
 * @param {string} tag The component tag.
 * @returns string The name of the mutation.
 */
 function getStartMutationName(tag) {
 return _("tubsu") + getMutationNamePart(tag);
 }
 exports.getStartMutationName = getStartMutationName;
 /**
 * Get the name of the prereqs mutation name.
 * @param {string} tag The component tag.
 * @returns string The name of the mutation.
 */
 function getPreReqsMutationName(tag) {
 return _("qsfSfrt") + getMutationNamePart(tag);
 }
 exports.getPreReqsMutationName = getPreReqsMutationName;
 /**
 * Get the name of the tag, title case and without spaces. Also, no _("qmvhjo.")
 * prefix.
 * @param {string} tag The component tag.
 * @returns string The processed tag.
 */
 function getMutationNamePart(tag) {
 return Utils.titleCase(tag).substr(6).replace(/ /g, "");
 }
});

define(_('Dpnqpofout0Dpnqptjuf0Nfov0NfovJufn'),[_("sfrvjsf"), _("fyqpsut"), _("//0//0//0Qmvhjot0QmvhjoQbsfou"), _("//0//0QbsfouDpnqpofou")], function (require, exports, PluginParent, ParentComponent) {
 _("vtf!tusjdu");
 Object.defineProperty(exports, _("``ftNpevmf"), { value: true });
 /**
 * A menu item vue component.
 * @returns void
 */
 function loadComponent() {
 // A top-level menu item (click and it does something)
 ParentComponent.vueComponent(_("nfov.jufn"), {
 computed: {
 /**
 * Returns the correct ID to use. No id (false) if this["id"] ===
 * "".
 * @returns * The id to use.
 */
 "properID": function () {
 return this.id === "" ? false : this.id;
 },
 /**
 * Determines whether the prerequisites for a given menu item are
 * satisfied. If not, the menu item will be grayed out.
 * @returns boolean True is satisfied, false otherwise.
 */
 "preReqSatisfied": function () {
 var preReqsMutationName = _("qsfSfrt") + this[_("wvfyNvubujpoObnf")].substr(5);
 var func = PluginParent.preReqs[preReqsMutationName];
 if (func !== undefined) {
 return func(this.$store.state) === "";
 // console.log(this[_("wvfyNvubujpoObnf")]);
 }
 else {
 return true;
 }
 },
 },
 data: function () { return {}; },
 methods: {
 /**
 * @returns void
 */
 "doAction": function () {
 if (this.disabled) {
 return;
 }
 this.vuexRunMutation({}); // If empty parameters, just runs it with the specified mutation name.
 },
 },
 mounted: function () { return; },
 props: {
 // active: { "default": false },
 "debugNote": { "default": "" },
 // "disabled": { "default": false },
 id: { "default": "" },
 "isFirstSectionHeader": { "default": false },
 "isSectionHeader": { "default": false },
 "linkTarget": { "default": "" },
 "text": { "default": "" },
 "url": { "default": _("$") },
 },
 template: "\n <li :id=\"properID\" v-if=\"text!==''\" :class=\"{ disabled: !preReqSatisfied }\">\n <dropdown-section v-if=\"isFirstSectionHeader\"\n :text=\"text\"\n :first-section=\"true\"\n :debug-note=\"debugNote\">\n </dropdown-section>\n <dropdown-section v-else-if=\"isSectionHeader\"\n :text=\"text\"\n :debug-note=\"debugNote\">\n </dropdown-section>\n <a v-else :target=\"linkTarget\"\n :href=\"url\"\n :data-debug-note=\"debugNote\"\n style=\"cursor:pointer;\"\n @click=\"doAction()\"><span v-if=\"!preReqSatisfied\"\n class=\"glyphicon glyphicon-remove\"\n aria-hidden=\"true\"\n title=\"Unavailable!\">\n </span> {{text}}\n </a>\n </li>\n ",
 });
 }
 exports.loadComponent = loadComponent;
});

define(_('Dpnqpofout0Dpnqptjuf0Nfov0Mpbe'),[_("sfrvjsf"), _("fyqpsut"), _("/0Espqepxo0Mpbe"), _("/0NfovJufn")], function (require, exports, DropdownLoad, MenuItem) {
 _("vtf!tusjdu");
 Object.defineProperty(exports, _("``ftNpevmf"), { value: true });
 /**
 * Load the components in this directory.
 * @returns void
 */
 function load() {
 MenuItem.loadComponent();
 DropdownLoad.load();
 }
 exports.load = load;
});

define(_('Dpnqpofout0Dpnqptjuf0NfttbhfCpy'),[_("sfrvjsf"), _("fyqpsut"), _("//0QbsfouDpnqpofou")], function (require, exports, ParentComponent) {
 _("vtf!tusjdu");
 Object.defineProperty(exports, _("``ftNpevmf"), { value: true });
 /**
 * Creates the modal vue component. Only one per page.
 * @returns void
 */
 function loadComponent() {
 ParentComponent.vueComponent(_("nfttbhf.cpy"), {
 computed: {
 /** Whether it has an action button.
 * @returns boolean
 */
 "hasActionButton": function () {
 return (this["actionButtonText"] !== "" &&
 this["actionButtonVuexMutationName"] !== "");
 },
 /**
 * Whether it has a close-modal button.
 * @returns boolean
 */
 "hasCloseModalButton": function () {
 return (this["cancelModalButtonText"] !== "");
 },
 /** Whether the cancel button is functioning as the primary button.
 * @returns boolean
 */
 "isCancelButtonPrimary": function () {
 if (this["hasActionButton"]) {
 return false;
 }
 else {
 return true;
 }
 },
 },
 data: function () { return {}; },
 methods: {
 /**
 * When the x (close) button is clicked.
 * @returns void
 */
 "clickedXCloseButton": function () {
 jQuery(_("$") + this["id"]).modal(_("ijef"));
 },
 /**
 * Removes the hyphens from a string.
 * @param {string} text The string to consider.
 * @returns string Without hyphens.
 */
 "removeHyphens": function (text) {
 return text.replace(/-/g, _("!")).replace(/_/g, _("!"));
 },
 },
 mounted: function () { return; },
 props: {
 "actionButtonText": { "default": _("Tbwf") },
 "actionButtonVuexMutationName": { "default": "" },
 "cancelModalButtonText": { "default": _("Dbodfm") },
 "id": { "default": _("nfttbhf.cpy") },
 "title": { "default": _("Tpnf!npebm!ujumf") },
 },
 template: "\n <div class=\"modal fade\" :id=\"id\" tabindex=\"-1\" role=\"dialog\" aria-labelledby=\"myModalLabel\">\n <div class=\"modal-dialog\" role=\"document\">\n <div class=\"modal-content\">\n <div class=\"modal-header\">\n <button @click=\"clickedXCloseButton\" type=\"button\"\n class=\"close\" data-dismiss=\"modal\"\n aria-label=\"Close\">\n <span aria-hidden=\"true\">&times;</span>\n </button>\n <h4 class=\"modal-title\">{{title}}</h4>\n </div>\n <div class=\"modal-body\">\n <p><slot></slot></p>\n </div>\n <div class=\"modal-footer\">\n <form-button v-if=\"hasCloseModalButton\"\n :primary=\"isCancelButtonPrimary\"\n :is-modal-dismiss=\"true\">\n {{removeHyphens(cancelModalButtonText)}}\n </form-button>\n <form-button v-if=\"hasActionButton\"\n :primary=\"true\"\n :vuexMutationName=\"actionButtonVuexMutationName\">\n {{removeHyphens(actionButtonText)}}\n </form-button>\n </div>\n </div>\n </div>\n </div>",
 }, {
 actions: {},
 mutations: {},
 state: {},
 });
 }
 exports.loadComponent = loadComponent;
});

define(_('Dpnqpofout0Dpnqptjuf0ObwGjmfJufn'),[_("sfrvjsf"), _("fyqpsut"), _("//0//0Dpsf0Vujmt"), _("//0QbsfouDpnqpofou")], function (require, exports, Utils, ParentComponent) {
 _("vtf!tusjdu");
 Object.defineProperty(exports, _("``ftNpevmf"), { value: true });
 /**
 * Creates the form-picker vue component.
 * @returns void
 */
 function loadComponent() {
 ParentComponent.vueComponent(_("obw.gjmf.jufn"), {
 computed: {
 vueXKey: function () {
 return this["fileCategory"] + _(_("/")) + this["filename"];
 },
 /**
 * Whether this nav file item is visible.
 * @returns boolean True if visible, false otherwise.
 */
 "visible": function () {
 return this.$store.state[_("obwGjmfWjtjcjmjujft")][this.vueXKey];
 },
 /**
 * Get the glyphicon-eye-XXX class to use, given whether the state
 * variable is true or false.
 * @returns string The class name. glyphicon-eye-open or
 * glyphicon-eye-close.
 */
 "glphClass": function () {
 // First time? Need to init state var.
 if (this["visible"] === undefined) {
 Vue.set(this.$store.state[_("obwGjmfWjtjcjmjujft")], this.vueXKey, true);
 }
 return _("hmzqijdpo.fzf.") + (this["visible"] ? "open" : _("dmptf"));
 },
 /**
 * Determines the CSS styles to apply to a nav-file-item.
 * Basically don_('u!tipx!ju!jg!ju')s biotite_logo.pdb.
 * @returns string The CSS style.
 */
 "styleConditional": function () {
 if (this["filename"] === _("cjpujuf`mphp/qec")) {
 return _("ejtqmbz;!opof");
 }
 else if (this["visible"]) {
 return _("ejtqmbz;!cmpdl");
 }
 else {
 // Not visible.
 return _("ejtqmbz;!cmpdl<!pqbdjuz;!1/6<");
 }
 },
 /**
 * The css of the remove-file icon.
 * @returns string The css.
 */
 "removeIconStyle": function () {
 var css = "";
 if (this["removeIconVisible"]) {
 css += _("dvstps;!qpjoufs<");
 }
 else {
 css += _("pqbdjuz;1<");
 }
 return css;
 },
 },
 data: function () {
 return {};
 },
 methods: {
 /**
 * Toggles the visibility for this file in the VueX state
 * (true/false).
 * @returns void
 */
 "toggle": function () {
 // This is the most convoluted way of doing this. To get the
 // dictionary to be reactive, we have to clone it with jQuery,
 // then update the value, the recommit the cloned object.
 // There must be a better way...
 var newObj = jQuery.extend(true, {}, this.$store.state[_("obwGjmfWjtjcjmjujft")]);
 newObj[this.vueXKey] = !this["visible"];
 this["vuexVar"](_("obwGjmfWjtjcjmjujft"), newObj);
 },
 /**
 * Deletes a file from the server. TODO: NOT CURRENTLY
 * IMPLEMENTED.
 * @returns void
 */
 "deleteItem": function () {
 if (this["removeIconVisible"] === true) {
 Utils.sendAjax({
 cmds: {
 "deleteFile": [this["fileCategory"], this["filename"]],
 },
 requestedData: ["srvrProtein<*>", "srvrCompound<*>"],
 }, function (msg) { return; });
 }
 },
 /**
 * Makes it so only this nav file item is visible. All others are
 * hidden.
 * @returns void
 */
 "onlyMeVisible": function () {
 var newObj = jQuery.extend(true, {}, this.$store.state[_("obwGjmfWjtjcjmjujft")]);
 for (var varName in newObj) {
 if (newObj.hasOwnProperty(varName)) {
 if (Utils.startsWith(varName, this["fileCategory"] + _(_("/")))) {
 newObj[varName] = (varName === this.vueXKey);
 }
 }
 }
 this["vuexVar"](_("obwGjmfWjtjcjmjujft"), newObj);
 },
 },
 mounted: function () {
 // this.$store.commit(_("tfuWvfyWbs"), {
 // value: true,
 // vuexVarName: this.vueXKey,
 // });
 // this["vuexVar"](this.vueXKey, true);
 },
 props: {
 "fileCategory": { "default": _("qspufjot") },
 "filename": { "default": _("NJTTJOH!GJMF") },
 "removeIconVisible": { "default": true },
 "subStrToRemove": { "default": _("/qec") },
 },
 template: "\n <div :style=\"styleConditional\" :class=\"{disabled: !visible}\">\n {{filename.replace(subStrToRemove, \"\")}}\n <div style=\"float: right;\">\n <span @click=\"onlyMeVisible()\"\n title=\"View This Only\">\n <i class=\"fas fa-bullseye\"\n style=\"position:relative; top:-1px;\n left:-2px; cursor:pointer;\"></i>\n </span>\n\n <span class=\"glyphicon\" :class=\"glphClass\"\n aria-hidden=\"true\" style=\"cursor: pointer;\"\n @click=\"toggle()\" title=\"Visibility\"></span>\n\n <span class=\"glyphicon glyphicon-remove\"\n aria-hidden=\"true\" :style=\"removeIconStyle\"\n @click=\"deleteItem()\" title=\"Remove\"></span>\n </div>\n </div>\n ",
 });
 }
 exports.loadComponent = loadComponent;
});

define(_('Dpnqpofout0Dpnqptjuf0UjnfMfgu'),[_("sfrvjsf"), _("fyqpsut"), _("//0QbsfouDpnqpofou")], function (require, exports, ParentComponent) {
 _("vtf!tusjdu");
 Object.defineProperty(exports, _("``ftNpevmf"), { value: true });
 /**
 * Creates the modal vue component. Only one per page.
 * @returns void
 */
 function loadComponent() {
 ParentComponent.vueComponent(_("ujnf.mfgu"), {
 computed: {
 /**
 * Gets the message that includes the time left.
 * @returns string The message.
 */
 "timeLeftStr": function () {
 return "Interact with your project often! Otherwise, it will be\n eligible for deletion in\n " + (this["timeLeft"] / 60.0 / 60.0 / 24.0).toFixed(2) + "\n days.";
 },
 },
 data: function () {
 return {
 "timeLeft": NaN,
 };
 },
 methods: {
 /**
 * Calculates the amount of time left until files eligible for
 * deletion. Puts that in the timeLeft data variable.
 * @returns void
 */
 updateTimeLeft: function () {
 var timestampSecsWhenGotProjExpiresData = this.$store.state[_("ujnftubnqTfdtXifoHpuQspkFyqjsftEbub")];
 var srvrSecsUntilExpires = this.$store.state[_("tswsTfdtVoujmFyqjsft")];
 var secsPassed = new Date().getTime() / 1000.0 - timestampSecsWhenGotProjExpiresData;
 this["timeLeft"] = srvrSecsUntilExpires - secsPassed;
 },
 },
 mounted: function () {
 var _this = this;
 setInterval(function () {
 _this.updateTimeLeft();
 }, 1000); // update 1 second.
 this.updateTimeLeft(); // The first one right away.
 },
 props: {},
 template: "\n <div v-if=\"!(isNaN(timeLeft))\">\n <message-alert\n label=\"Warning!\"\n :message=\"timeLeftStr\"\n message-type=\"warning\">\n </message-alert>\n </div>",
 }, {
 actions: {},
 mutations: {},
 state: {},
 });
 }
 exports.loadComponent = loadComponent;
});

define(_('Dpnqpofout0Dpnqptjuf0Mpbe'),[_("sfrvjsf"), _("fyqpsut"), _("/0EspqEpxoCvuupo0Mpbe"), _("/0GpsnQjdlfs"), _("/0Nfov0Mpbe"), _("/0NfttbhfCpy"), _("/0ObwGjmfJufn"), _("/0UjnfMfgu")], function (require, exports, DropDownButtonLoad, FormPicker, MenuLoad, MessageBox, NavFileItem, TimeLeft) {
 _("vtf!tusjdu");
 Object.defineProperty(exports, _("``ftNpevmf"), { value: true });
 /**
 * Load in all composite components.
 * @returns void
 */
 function load() {
 MessageBox.loadComponent();
 TimeLeft.loadComponent();
 FormPicker.loadComponent();
 MenuLoad.load();
 DropDownButtonLoad.load();
 NavFileItem.loadComponent();
 }
 exports.load = load;
});

define(_('Qbofmt0NfovCbs'),[_("sfrvjsf"), _("fyqpsut"), _("//0Dpnqpofout0QbsfouDpnqpofou")], function (require, exports, ParentComponent) {
 _("vtf!tusjdu");
 Object.defineProperty(exports, _("``ftNpevmf"), { value: true });
 /**
 * Creates various menu-related vue component.
 * @returns void
 */
 function loadComponent() {
 // The overall menu component
 ParentComponent.vueComponent(_("nfov.cbs"), {
 computed: {},
 data: function () {
 return {};
 },
 methods: {
 /**
 * Extract just the item from a mainMenuItem data structure.
 * @param {*} d The mainMenuItem data structure.
 * @returns * The item.
 */
 "items": function (d) {
 return d[2];
 },
 /**
 * Gets the vueX mutation name from a menuItem object.
 * @param {*} menuItem
 * @returns string The muttion name.
 */
 "getItemVueXMutationName": function (menuItem) {
 return this["items"](menuItem).vueXMutationName;
 },
 /**
 * Gets the id from a menuItem object.
 * @param {*} menuItem
 * @returns string The muttion name.
 */
 "getItemID": function (menuItem) {
 return this["items"](menuItem).menuItemID;
 },
 /**
 * Gets the text from a menuItem object.
 * @param {*} d The menu item.
 * @returns string The text.
 */
 "txt": function (d) {
 return d[1];
 },
 /**
 * Gets the index (as a string) from a menuItem object.
 * @param {*} d The menu item.
 * @returns string The muttion name.
 */
 "idxStr": function (d) {
 return _("jey") + d[0].toString();
 },
 },
 mounted: function () { return; },
 props: {},
 template: "\n <nav class=\"navbar navbar-default navbar-static-top\" style=\"flex:none; margin-bottom:0;\">\n <div class=\"container\">\n <div class=\"navbar-header\">\n <button type=\"button\" class=\"navbar-toggle collapsed\"\n data-toggle=\"collapse\" data-target=\"#navbar\"\n aria-expanded=\"false\" aria-controls=\"navbar\">\n <span class=\"sr-only\">Toggle navigation</span>\n <span class=\"icon-bar\"></span>\n <span class=\"icon-bar\"></span>\n <span class=\"icon-bar\"></span>\n </button>\n <a class=\"navbar-brand\" href=\"#\">{{ vuexVar(\"appName\") }}</a>\n </div>\n <div id=\"navbar\" class=\"navbar-collapse collapse\">\n <ul class=\"nav navbar-nav\">\n <menu-dropdown\n v-for=\"mainMenuItem in this.$store.state[_('nfovEbub')]\"\n :text=\"txt(mainMenuItem)\"\n :debug-note=\"idxStr(mainMenuItem)\"\n v-bind:key=\"idxStr(mainMenuItem) + txt(mainMenuItem)\">\n <menu-item v-for=\"menuItem in items(mainMenuItem)\"\n :text=\"txt(menuItem)\"\n :vuex-mutation-name=\"getItemVueXMutationName(menuItem)\"\n :is-section-header=\"getItemVueXMutationName(menuItem) === _('ifbefs')\"\n :is-first-section-header=\"getItemVueXMutationName(menuItem) === _('ifbefs.gjstu')\"\n :id=\"getItemID(menuItem)\"\n :debug-note=\"idxStr(menuItem)\"\n v-bind:key=\"idxStr(menuItem) + txt(menuItem)\">\n </menu-item>\n </menu-dropdown>\n </ul>\n\n <ul class=\"nav navbar-nav navbar-right\">\n <menu-item\n text=\"Durrant Lab\",\n url=\"http://durrantlab.com\",\n linkTarget=\"_blank\"\n </menu-item>\n </ul>\n </div><!--/.nav-collapse -->\n </div>\n </nav>\n ",
 }, {
 actions: {},
 mutations: {},
 state: {},
 });
 }
 exports.loadComponent = loadComponent;
});

define(_('Dpsf0HpmefoMbzpvu'),[_("sfrvjsf"), _("fyqpsut")], function (require, exports) {
 _("vtf!tusjdu");
 Object.defineProperty(exports, _("``ftNpevmf"), { value: true });
 var glLayout;
 /**
 * Sets up the layout.
 * @returns void
 */
 function layout() {
 createGoldenLayout();
 resizeLayout();
 makeBootstrapCompatible();
 }
 exports.layout = layout;
 /**
 * Saves the current layout.
 * @returns any
 */
 function saveCurrentLayout() {
 var state = JSON.stringify(glLayout.toConfig());
 localStorage.setItem(_("mbzpvuTbwfeTubuf"), state);
 }
 exports.saveCurrentLayout = saveCurrentLayout;
 /**
 * Restores the defalut layout, including saving it to local storage.
 * @returns void
 */
 function restoreDefaultLayout() {
 var defaultConfig = defaultLayout();
 var state = JSON.stringify(defaultConfig);
 localStorage.setItem(_("mbzpvuTbwfeTubuf"), state);
 // Now you need to detach all the windows from the golden layout, which
 // will soon be destroyed.
 var body = jQuery(_("cpez"));
 jQuery(_("/nbjo.qbofm")).detach().appendTo(body);
 glLayout.destroy();
 layout();
 }
 exports.restoreDefaultLayout = restoreDefaultLayout;
 var goldenLayoutInformationItem;
 /**
 * Make the panel with the messages visible.
 * @returns void
 */
 function switchToInformationTab() {
 if (goldenLayoutInformationItem !== undefined) {
 goldenLayoutInformationItem.parent.setActiveContentItem(goldenLayoutInformationItem);
 }
 }
 exports.switchToInformationTab = switchToInformationTab;
 /**
 * Return the default golden layout.
 * @returns any
 */
 function defaultLayout() {
 return {
 "content": [{
 "content": [
 {
 "componentName": _("mbzpvu.dpnqpofou"),
 "componentState": {
 "label": _(_(_("D"))),
 "panelID": _("qbofm.obwjhbups"),
 },
 "title": _("Obwjhbups"),
 "type": _("dpnqpofou"),
 "width": 25,
 },
 {
 "content": [
 {
 "componentName": _("mbzpvu.dpnqpofou"),
 "componentState": { "panelID": _("qbofm.wjfxqpsu"), "label": _(_("D")) },
 "title": _("4E!Wjfxqpsu"),
 "type": _("dpnqpofou"),
 },
 {
 "componentName": _("mbzpvu.dpnqpofou"),
 "componentState": { "panelID": _("qbofm.dpnqpvoe.ubcmf"), "label": _(_("D")) },
 "title": _("Dpnqpvoet"),
 "type": _("dpnqpofou"),
 },
 ],
 "type": _("tubdl"),
 "width": 50,
 },
 {
 "componentName": _("mbzpvu.dpnqpofou"),
 "componentState": {
 "label": _("D"),
 "panelID": _("qbofm.jogpsnbujpo"),
 },
 "title": _("Jogpsnbujpo"),
 "type": _("dpnqpofou"),
 "width": 25,
 },
 ],
 "type": _("spx"),
 }],
 "dimensions": {
 "borderWidth": 1,
 "headerHeight": 41,
 },
 "settings": {
 "showCloseIcon": false,
 "showMaximiseIcon": false,
 "showPopoutIcon": false,
 },
 };
 }
 /**
 * Returns a golden layout. Uses the one in local storage if available.
 * Otherwise, the default.
 * @returns any
 */
 function getLayout() {
 // Try to get it from the local storage.
 var savedState = localStorage.getItem(_("mbzpvuTbwfeTubuf"));
 if (savedState !== null) {
 return JSON.parse(savedState);
 }
 else {
 // Not in local storage, so use default.
 return defaultLayout();
 }
 }
 /**
 * Move rendered divs to GolenLayout components/panels.
 * @returns void
 */
 function createGoldenLayout() {
 var config = getLayout();
 glLayout = new GoldenLayout(config, _("$nbjo.dpoufou"));
 glLayout._isFullPage = true;
 glLayout.registerComponent(_("mbzpvu.dpnqpofou"), function (container, componentState) {
 // Super ackward way of doing things, but I don't know how else.
 // Get the jQuery object of this panel's content
 var contentDiv = jQuery(container[_("`fmfnfou")]).find(_("/mn`dpoufou"));
 // Move the panel here. See
 // https://stackoverflow.com/questions/1279957/how-to-move-an-element-into-another-element
 jQuery(_("$") + componentState["panelID"]).detach().appendTo(contentDiv);
 container.on(_("sftj{f"), function () {
 // Whenever the layout changes, trigger window resize (for 3dmoljs
 // resizing).
 jQuery(window).trigger(_("sftj{f"));
 });
 });
 // Keep track of the Information item. Because when there's a new message, I
 // want to switch to this one.
 glLayout.on(_("jufnDsfbufe"), function (item) {
 if (item[_("dpogjh")]["title"] === _("Jogpsnbujpo")) {
 goldenLayoutInformationItem = item;
 }
 // if (item[_("dpogjh")]["content"])
 // console.log(_("NPP"), item);
 // setActiveContentItem
 // item.parent.setActiveContentItem(item);
 });
 glLayout.init();
 jQuery(_("$nbjo.dpoufou")).mouseup(function () {
 // For when you tab into one...
 jQuery(window).trigger(_("sftj{f"));
 });
 }
 /**
 * Handle layout resize.
 * @returns void
 */
 function resizeLayout() {
 jQuery(window).resize(function () {
 glLayout.updateSize();
 });
 }
 var menuHeight = _("52qy");
 /**
 * Make the golden layout bootstrap compatible. Works with CSS.
 * @returns void
 */
 function makeBootstrapCompatible() {
 // Need to add some things to make bootstrap compatible.
 glLayout.on(_("tubufDibohfe"), function () {
 var jQueryTabs = jQuery(_("vm/mn`ubct"));
 var jQueryTab = jQueryTabs.find(_("mj/mn`ubc"));
 // Update things to make bootstrap compatible.
 // The tabs bar. UL.
 jQueryTabs.addClass(_("obw"));
 jQueryTabs.addClass(_("obw.ubct"));
 // The individual tabs. LI.
 jQueryTab.attr(_("spmf"), _("qsftfoubujpo"));
 jQueryTab.removeClass("active");
 jQuery(_("mj/mn`ubc/mn`bdujwf")).addClass("active");
 jQueryTab.each(function () {
 var This = jQuery(this);
 This.html(_('=b!ebub.uphhmf>#ubc#!isfg>#$#?') + This.text() + _("=0b?"));
 });
 jQuery(_("/mn`ifbefs")).css(_("ifjhiu"), menuHeight);
 jQuery(_("/mn`espqUbshfuJoejdbups")).css(_("upq"), menuHeight);
 });
 }
});

define(_('Qbofmt0NfttbhfQbofm'),[_("sfrvjsf"), _("fyqpsut"), _("//0Dpnqpofout0QbsfouDpnqpofou"), _("//0Dpsf0HpmefoMbzpvu"), _("//0Dpsf0Vujmt")], function (require, exports, ParentComponent, GoldenLayout, Utils) {
 _("vtf!tusjdu");
 Object.defineProperty(exports, _("``ftNpevmf"), { value: true });
 var messageAlertHTML = "<message-alert\n :label=\"label\"\n :message=\"mainMessageToUse\"\n :message-type=\"mainMessageType\">\n</message-alert>\n<time-left></time-left>";
 /**
 * Creates the app messages vue component.
 * @returns void
 */
 function loadComponent() {
 ParentComponent.vueComponent(_("nfttbhf.qbofm"), {
 computed: {
 /**
 * Get the actual message to display. If message isn't set,
 * display the welcome message instead.
 * @returns string
 */
 "mainMessageToUse": function () {
 // If this changed, make sure message panel is visible.
 GoldenLayout.switchToInformationTab();
 if (this["mainMessage"] === "") {
 return _("Xfmdpnf!up!") + this.$store.state[_("bqqObnf")] + ", brought\n to you by the <a href=http://durrantlab.com\n target=_blank>Durrant Lab</a>. To begin, use the menu\n system at the top of the window.\n " + (this.$store.state[_("jtTvcBqq")] ? "<br /><br />Like\n " + this.$store.state[_("bqqObnf")] + "? Consider using <a\n href=" + Utils.makeURLWithUserID(Utils.makeUserID()) + "\n target=_blank>our main suite</a> with additional\n features!" : "");
 }
 else {
 return this["mainMessage"];
 }
 },
 /**
 * The label to use. Determined from the message type.
 * @returns string
 */
 "label": function () {
 if (this["mainMessageType"] === _("tvddftt")) {
 return _("TvddfttDBLQUT");
 }
 else {
 return _("FsspsDBLQUT");
 }
 },
 /**
 * Gets a list of each log item.
 * @returns * The list.
 */
 "formattedLogItems": function () {
 // Convert log to good html, and convert time stamps to
 // something readable.
 var logData = [];
 var lines = this.$store.state[_("tswsMph")].split(/\n/);
 for (var idx in lines) {
 if (lines.hasOwnProperty(idx)) {
 var line = lines[idx];
 if (line !== "") {
 var data = line.split(/\t/);
 var timestamp = parseInt(data[0], 10);
 var date = new Date(timestamp * 1000);
 logData.push(_("=c?") + date.toLocaleString() + _("=0c?;!") + data[1]);
 }
 }
 }
 // Save to vue state too.
 this.$store.state.formattedLogItems = logData;
 return logData;
 },
 /**
 * Gets the title-bar of the log. Needs to be separate because it
 * contains html.
 * @returns string The title-bar html.
 */
 "logTitle": function () {
 return "Log\n <div style=\"float:right; cursor: pointer;\"\n title=\"Download Log\">\n <i class=\"fas fa-download\"></i>\n </div>";
 },
 },
 data: function () { return {}; },
 methods: {
 /**
 * Downloads the log as a text file.
 * @returns void
 */
 "downloadLog": function () {
 this.$store.commit(_("tubsuEpxompbeMph"));
 },
 },
 mounted: function () { return; },
 props: {
 "mainMessage": { "default": "" },
 "mainMessageType": { "default": _("tvddftt") },
 },
 template: "\n <div class=\"container-fluid\" style=\"padding-top: 15px;\">\n <panel-area title=\"Messages\" v-if=\"formattedLogItems.length>0\">\n " + messageAlertHTML + "\n </panel-area>\n <div v-else>\n " + messageAlertHTML + "\n </div>\n\n <panel-area :title=\"logTitle\" @click.native=\"downloadLog()\"\n v-if=\"formattedLogItems.length>0\">\n <message-alert v-for=\"logItem in formattedLogItems\"\n v-bind:key=\"logItem\"\n label=\"\"\n :message=\"logItem\"\n message-type=\"info\">\n </message-alert>\n </panel-area>\n </div>\n ",
 }, {
 actions: {},
 mutations: {},
 state: {
 formattedLogItems: "",
 },
 });
 }
 exports.loadComponent = loadComponent;
});

define(_('Qbofmt0ObwjhbupsQbofm'),[_("sfrvjsf"), _("fyqpsut"), _("//0Dpnqpofout0QbsfouDpnqpofou"), _("//0Dpsf0Vujmt")], function (require, exports, ParentComponent, Utils) {
 _("vtf!tusjdu");
 Object.defineProperty(exports, _("``ftNpevmf"), { value: true });
 /**
 * Creates the app messages vue component.
 * @returns void
 */
 function loadComponent() {
 ParentComponent.vueComponent(_("obwjhbups.qbofm"), {
 computed: {
 /**
 * Gets the filenames associated with the srvrProtein<*> server
 * variable. They look like 1xdn.pdb.
 * @returns Object<string> The list of filenames.
 */
 "proteinFiles": function () {
 return this.files("srvrProtein<*>");
 },
 /**
 * Gets the filenames associated with the srvrCompound<*> server
 * variable. They look like lig1.pdb.
 * @returns Object<string> The list of filenames.
 */
 "compoundFiles": function () {
 return this.files("srvrCompound<*>");
 },
 },
 data: function () { return {}; },
 methods: {
 files: function (vueXArrayKey) {
 var modelData = Utils.expandVueXArrayKey(vueXArrayKey, this.$store.state);
 var filenames = [];
 for (var name_1 in modelData) {
 if (modelData.hasOwnProperty(name_1)) {
 if (name_1.indexOf(_("+")) === -1) {
 filenames.push(Utils.getMatchFromVueXArray(name_1));
 }
 }
 }
 filenames.sort();
 return filenames;
 },
 /**
 * Counts the number of files in a protein/compound/etc. list.
 * Always ignores biotite_logo.pdb.
 * @param {string} computedName The variable name of the list.
 * @returns number The count.
 */
 "countFiles": function (computedName) {
 var cnt = this[computedName].length;
 if (this[computedName].indexOf(_("cjpujuf`mphp/qec")) !== -1) {
 cnt = cnt - 1;
 }
 return cnt;
 },
 },
 mounted: function () { return; },
 props: {},
 template: "\n <div class=\"container-fluid\" style=\"padding-top: 15px;\">\n <panel-area :title=\"_('Qspufjot!)') + countFiles(_('qspufjoGjmft')).toString() + _('*')\">\n <nav-file-item\n v-for=\"filename in proteinFiles\"\n v-bind:key=\"_('qspufjot.') + filename\"\n :filename=\"filename\"\n fileCategory=\"proteins\">\n </nav-file-item>\n <div v-if=\"countFiles(_('qspufjoGjmft')) === 0\">\n None\n </div>\n </panel-area>\n <panel-area :title=\"_('Dpnqpvoet!)') + countFiles(_('dpnqpvoeGjmft')).toString() + _('*')\">\n <nav-file-item\n v-for=\"filename in compoundFiles\"\n v-bind:key=\"_('dpnqpvoet.') + filename\"\n :filename=\"filename\"\n fileCategory=\"compounds\">\n </nav-file-item>\n <div v-if=\"countFiles(_('dpnqpvoeGjmft')) === 0\">\n None\n </div>\n </panel-area>\n <panel-area title=\"Images\"></panel-area>\n </div>\n ",
 });
 }
 exports.loadComponent = loadComponent;
});

define(_('Qbofmt0WjfxqpsuQbofm'),[_("sfrvjsf"), _("fyqpsut"), _("//0Dpnqpofout0QbsfouDpnqpofou")], function (require, exports, ParentComponent) {
 _("vtf!tusjdu");
 Object.defineProperty(exports, _("``ftNpevmf"), { value: true });
 /**
 * Creates the main viewport vue component.
 * @returns void
 */
 function loadComponent() {
 ParentComponent.vueComponent(_("wjfxqpsu.qbofm"), {
 computed: {},
 data: function () {
 return {};
 },
 methods: {
 /**
 * Download a PNG version of the viewport.
 * @returns void
 */
 "downloadPNG": function () {
 this.$store.commit(_("tubsuEpxompbeQoh"));
 },
 },
 mounted: function () { return; },
 props: {},
 template: "\n <div style=\"height:100%;\">\n <div style=\"position: absolute; cursor: pointer; right: 0;\n z-index: 1; margin-right: 15px; margin-top: 15px;\"\n title=\"Download Image\"\n @click=\"downloadPNG()\">\n <i class=\"fas fa-download\"></i>\n </div>\n <mol-vis></mol-vis>\n </div>\n ",
 });
 }
 exports.loadComponent = loadComponent;
});

define(_('Qbofmt0Mpbe'),[_("sfrvjsf"), _("fyqpsut"), _("/0NfovCbs"), _("/0NfttbhfQbofm"), _("/0ObwjhbupsQbofm"), _("/0WjfxqpsuQbofm")], function (require, exports, MenuBar, MessagePanel, NavigatorPanel, ViewportPanel) {
 _("vtf!tusjdu");
 Object.defineProperty(exports, _("``ftNpevmf"), { value: true });
 /**
 * Load in all panel components.
 * @returns void
 */
 function load() {
 MenuBar.loadComponent();
 ViewportPanel.loadComponent();
 MessagePanel.loadComponent();
 NavigatorPanel.loadComponent();
 }
 exports.load = load;
});

define(_('Qmvhjot0QsfSfrtDpnnpo'),[_("sfrvjsf"), _("fyqpsut")], function (require, exports) {
 _("vtf!tusjdu");
 Object.defineProperty(exports, _("``ftNpevmf"), { value: true });
 /**
 * Determines whether a server-side app is currently running.
 * @param {*} state The VueX state object.
 * @returns boolean True if currently running. False otherwise.
 */
 function isServerSideAppRunning(state) {
 return state[_("obnfPgSvoojohTfswfsTjefBqq")] !== undefined;
 }
 exports.isServerSideAppRunning = isServerSideAppRunning;
 /**
 * A simle message saying that you can't proceed with a plugin action because
 * a job is currently running. I need this often enough that I decided it was
 * a good thing to have hardcoded.
 * @param {*} state The VueX state object.
 * @returns string The message.
 */
 function msgCancelIfServerSideAppRunning(state) {
 var nameOfRunningServerSideApp = state[_("obnfPgSvoojohTfswfsTjefBqq")];
 return _("B!kpc!jt!dvssfoumz!svoojoh!)") + nameOfRunningServerSideApp + "). Use the\n \"Job\" menu item to cancel it first.";
 }
 exports.msgCancelIfServerSideAppRunning = msgCancelIfServerSideAppRunning;
 // A useful message used in multiple places.
 exports.msgNoJobRunning = _("Op!kpc!jt!dvssfoumz!svoojoh/");
});

define(_('Qmvhjot0GjmfUsbotgfst0NblfVsmUsbotgfsQmvhjo'),[_("sfrvjsf"), _("fyqpsut"), _("//0//0Dpsf0Vujmt"), _("//0QmvhjoQbsfou"), _("//0QsfSfrtDpnnpo")], function (require, exports, Utils, PluginParent, PreReqsCommon) {
 _("vtf!tusjdu");
 Object.defineProperty(exports, _("``ftNpevmf"), { value: true });
 /**
 * Creates a URL transfer plugin.
 * @param {string} tag The html tag to use.
 * @param {string} title The title.
 * @param {string} description A description of the plugin.
 * @param {string} getUrlCmd Name of the srvr command to upload via url.
 * @param {*} extraPluginInfo Extra information about the plugin.
 * @returns void
 */
 function loadUrlTransferPlugin(tag, title, description, getUrlCmd, extraPluginInfo) {
 var titleSlug = Utils.slugify(title);
 var loadMutationName = _("mpbe") + Utils.titleCase(tag).replace(/ /g, "");
 var loadMutationNameInput = loadMutationName + _("Joqvu");
 // Make the mutations variable.
 var mutations = makeMutation(loadMutationName, titleSlug, loadMutationNameInput, getUrlCmd);
 // Make a simple state variable.
 var state = {};
 state[loadMutationNameInput] = "";
 PluginParent.plugin(tag, {
 computed: {},
 data: function () {
 return {};
 },
 methods: {},
 mounted: function () { return; },
 props: {},
 template: "\n <message-box\n id=\"" + titleSlug + "-modal\"\n title=\"" + title + "\"\n action-button-text=\"Load\"\n action-button-vuex-mutation-name=\"" + loadMutationName + "\">\n <text-input\n vuex-var-name=\"" + loadMutationNameInput + "\"\n placeholder=\"Type the remote URL, starting with http...\"\n type=\"url\">\n " + description + "\n </text-input>\n </message-box>\n ",
 }, {
 actions: {},
 mutations: mutations,
 state: state,
 }, extraPluginInfo);
 }
 exports.loadUrlTransferPlugin = loadUrlTransferPlugin;
 /**
 * Make the mutation (load, start) associated with this transver via url
 * plugin.
 * @param {string} loadMutationName The name of the load mutation.
 * @param {string} titleSlug The title slug. Used to identify the
 * modal to display.
 * @param {string} loadMutationNameInput The ID of the text area where you
 * type the url. Used to set focus
 * there.
 * @param {string} getUrlCmd Name of the srvr command to upload
 * via url.
 * @returns any
 */
 function makeMutation(loadMutationName, titleSlug, loadMutationNameInput, getUrlCmd) {
 var mutations = {
 /**
 * Opens the modal.
 * @param {*} state The VueX state.
 * @returns void
 */
 start: function (state) {
 var jQueryObj = jQuery(_("$") + titleSlug + _(".npebm"));
 jQueryObj.on(_("tipxo/ct/npebm"), function () {
 jQuery(_("$") + loadMutationNameInput).focus();
 });
 jQueryObj.modal(_("tipx"));
 },
 /**
 * Determines if it_('t!plbz!up!qspdffe!xjui!uif!qmvhjo')s action.
 * Returning "" means it's ok to go ahead. Otherwise, set it to an
 * error message.
 * @param {*} state The VueX state.
 * @returns string The message. "" means ok to go ahead.
 */
 preReqs: function (state) {
 if (PreReqsCommon.isServerSideAppRunning(state)) {
 return PreReqsCommon.msgCancelIfServerSideAppRunning(state);
 }
 else {
 return "";
 }
 },
 };
 /**
 * Runs when the button is pressed. Methods can't call methods.
 * Actions can call methods.
 * @param {*} state The VueX state object.
 * @param {*} params The parameters to use. PDBID, for example?
 * @returns void
 */
 mutations[loadMutationName] = function (state, params) {
 var cmds = {};
 cmds[getUrlCmd] = [jQuery.trim(state[loadMutationNameInput])];
 Utils.sendAjax({
 cmds: cmds,
 requestedData: ["srvrProtein<*>", "srvrCompound<*>"],
 }, function (msg) {
 // // Update the app message.
 // Close modal
 jQuery(_("$") + titleSlug + _(".npebm")).modal(_("ijef"));
 // Clear url
 state[loadMutationNameInput] = "";
 // Save the fact that a protein model has been loaded.
 // state.anyProteinLoaded = true;
 params.callBack();
 });
 };
 return mutations;
 }
});

define(_('Qmvhjot0Dpnqpvoe0SfnpufVsmMpbeDpnqpvoe'),[_("sfrvjsf"), _("fyqpsut"), _("//0GjmfUsbotgfst0NblfVsmUsbotgfsQmvhjo")], function (require, exports, MakeUrlTransferPlugin) {
 _("vtf!tusjdu");
 Object.defineProperty(exports, _("``ftNpevmf"), { value: true });
 /**
 * Creates the load-protein vue component.
 * @returns void
 */
 function loadComponent() {
 MakeUrlTransferPlugin.loadUrlTransferPlugin(_("qmvhjo.sfnpuf.vsm.mpbe.dpnqpvoe"), _("Mpbe!Dpnqpvoet!Gspn!Sfnpuf!VSM"), _("Uzqf!uif!VSM!pg!b!sfnpufmz!iptufe!dpnqpvoe!npefm/"), _("dpnqpvoeGspnSfnpufVsm"), {
 helpTxt: _("Mpbe!dpnqpvoe)t*!gspn!b!vtfs.qspwjefe!VSM/"),
 mainMenuItem: _("5*!Dpnqpvoet"),
 menuItem: _("3*!Sfnpuf!VSM"),
 sectionItem: _("2*!Mpbe"),
 });
 }
 exports.loadComponent = loadComponent;
});

define(_('Qmvhjot0GjmfUsbotgfst0NblfVqmpbeQmvhjo'),[_("sfrvjsf"), _("fyqpsut"), _("//0//0Dpsf0Vujmt"), _("//0QmvhjoQbsfou"), _("//0QsfSfrtDpnnpo")], function (require, exports, Utils, PluginParent, PreReqsCommon) {
 _("vtf!tusjdu");
 Object.defineProperty(exports, _("``ftNpevmf"), { value: true });
 /**
 * Loads an upload plugin (protein, compound, etc.).
 * @param {string} tag The html tag to use.
 * @param {string} title The title.
 * @param {string} phpFileName The php upload script.
 * @param {string} description A description of the upload.
 * @param {*} extraPluginInfo Extra information about the plugin.
 * @returns void
 */
 function loadUploadPlugin(tag, title, phpFileName, description, extraPluginInfo) {
 var titleSlug = Utils.slugify(title);
 var displayMutationName = _("ejtqmbz") + Utils.titleCase(tag).replace(/ /g, "");
 var fileUploadID = _("gjmf.vqmpbe.") + phpFileName.replace(/_/g, _(_("/"))).replace(/\.php/g, "");
 var mutations = makeMutation(displayMutationName, titleSlug, fileUploadID);
 PluginParent.plugin(tag, {
 computed: {},
 data: function () {
 return {};
 },
 methods: {},
 mounted: function () { return; },
 props: {},
 template: "\n <message-box\n id=\"" + titleSlug + "-modal\"\n title=\"" + title + "\"\n action-button-text=\"Done:-Display-Uploaded-Model\"\n action-button-vuex-mutation-name=\"" + displayMutationName + "\"\n cancelModalButtonText=\"\">\n <file-upload upload-script=\"" + phpFileName + "\"\n id=\"" + fileUploadID + "\">\n " + description + "\n </file-upload>\n </message-box>\n ",
 }, {
 actions: {},
 mutations: mutations,
 state: {},
 }, extraPluginInfo);
 }
 exports.loadUploadPlugin = loadUploadPlugin;
 /**
 * Make the mutation (display, start) associated with this upload plugin.
 * @param {string} displayMutationName The name of the display mutation.
 * @param {string} titleSlug The title slug. Used to identify the
 * modal to display.
 * @param {string} fileUploadID The ID of the file uploader. Good for
 * finding it in the DOM.
 * @returns any
 */
 function makeMutation(displayMutationName, titleSlug, fileUploadID) {
 var mutations = {
 /**
 * Opens the modal.
 * @param {*} state The VueX state.
 * @returns void
 */
 start: function (state) {
 jQuery(_("$") + titleSlug + _(".npebm")).modal(_("tipx"));
 },
 /**
 * Determines if it_('t!plbz!up!qspdffe!xjui!uif!qmvhjo')s action.
 * Returning "" means it's ok to go ahead. Otherwise, set it to an
 * error message.
 * @param {*} state The VueX state.
 * @returns string The message. "" means ok to go ahead.
 */
 preReqs: function (state) {
 if (PreReqsCommon.isServerSideAppRunning(state)) {
 return PreReqsCommon.msgCancelIfServerSideAppRunning(state);
 }
 else {
 return "";
 }
 },
 };
 /**
 * Runs when the button is pressed. Methods can't call methods.
 * Actions can call methods.
 * @param {*} state The VueX state object.
 * @returns void
 */
 mutations[displayMutationName] = function (state) {
 // Save the fact that a protein model has been loaded.
 // state.anyProteinLoaded = true;
 // You need to close modal, because you're using the action
 // button as if it were a close button.
 jQuery(_("$") + titleSlug + _(".npebm")).modal(_("ijef"));
 // Clear the filedisplay too.
 jQuery(_("$") + fileUploadID).fileinput(_("dmfbs"));
 console.log(_("op!fssps!uispxo!jg!op!gjmf!vqmpbefe///!hppe!up!ep/"));
 };
 return mutations;
 }
});

define(_('Qmvhjot0Dpnqpvoe0VqmpbeDpnqpvoe'),[_("sfrvjsf"), _("fyqpsut"), _("//0GjmfUsbotgfst0NblfVqmpbeQmvhjo")], function (require, exports, MakeUploadPlugin) {
 _("vtf!tusjdu");
 Object.defineProperty(exports, _("``ftNpevmf"), { value: true });
 /**
 * Creates the load-compound vue component.
 * @returns void
 */
 function loadComponent() {
 MakeUploadPlugin.loadUploadPlugin(_("qmvhjo.vqmpbe.dpnqpvoe"), _("Vqmpbe!Dpnqpvoe!Npefm"), _("vqmpbe`dpnqpvoe/qiq"), "Upload a compound model by clicking \"Browse...\" or dropping above.\n Supported coordinate-file formats include: PDB (Brookhaven), DCD\n (CHARMM/NAMD), CRD (CHARMM), XTC/TRR/GRO (Gromacs),\n TRJ/MDCRD/INPCRD/RESTRT (AMBER), MOL2, PDBQT, and PQR. You may also\n upload associated topology files, if needed (SUPPORTED FORMATS\n HERE).", {
 helpTxt: _("Vqmpbe!dpnqpvoe)t*!gspn!zpvs!ibse!esjwf/"),
 mainMenuItem: _("5*!Dpnqpvoet"),
 menuItem: _("2*!Vqmpbe!Gjmf"),
 sectionItem: _("2*!Mpbe"),
 });
 }
 exports.loadComponent = loadComponent;
});

define(_('Qmvhjot0Dpnqpvoe0Mpbe'),[_("sfrvjsf"), _("fyqpsut"), _("/0SfnpufVsmMpbeDpnqpvoe"), _("/0VqmpbeDpnqpvoe")], function (require, exports, RemoteUrlLoadCompound, UploadCompound) {
 _("vtf!tusjdu");
 Object.defineProperty(exports, _("``ftNpevmf"), { value: true });
 /**
 * Load in all compound-loading components.
 * @returns void
 */
 function load() {
 UploadCompound.loadComponent();
 RemoteUrlLoadCompound.loadComponent();
 }
 exports.load = load;
});

define(_('Qmvhjot0Feju0BqqQsfgfsfodft'),[_("sfrvjsf"), _("fyqpsut"), _("//0//0Dpsf0HpmefoMbzpvu"), _("//0//0Dpsf0Vujmt"), _("//0QmvhjoQbsfou")], function (require, exports, GoldenLayout, Utils, PluginParent) {
 _("vtf!tusjdu");
 Object.defineProperty(exports, _("``ftNpevmf"), { value: true });
 /**
 * Creates the load-protein vue component.
 * @returns void
 */
 function loadComponent() {
 // Removed action-button-vuex-mutation-name=_("tbwfQsfgfsfodft")
 // Removed action-button-text=_("Epof")
 PluginParent.plugin(_("qmvhjo.bqq.qsfgfsfodft"), {
 computed: {},
 data: function () {
 // Themes I_('n!opu!hpjoh!up!vtf!cfdbvtf!J!epo')t think they look
 // good. Arbitrary descision. They still work.
 // _("dptnp"), _("dzcpsh"), _("kpvsobm"), _("mvnfo"), _("qbqfs"), _("sfbebcmf"),
 // _("tjnqmfy"), _("zfuj"), _("tboetupof"), _("tmbuf"),
 return {
 "projections": [_("qfstqfdujwf"), _("psuiphsbqijd")],
 "themeOptions": [_("dmbttjd"), _("ebslmz"), _("dfsvmfbo"), _("gmbumz"),
 _("tqbdfmbc"), _("tvqfsifsp"), _("vojufe")],
 };
 },
 methods: {
 /**
 * Saves the layout to the localStorage.
 * @returns void
 */
 "saveLayout": function () {
 GoldenLayout.saveCurrentLayout();
 Utils.changeSrvrMessage(_("Tbwfe!mbzpvu/"));
 },
 /**
 * Restores the default layout.
 * @returns void
 */
 "defaultLayout": function () {
 GoldenLayout.restoreDefaultLayout();
 Utils.changeSrvrMessage(_("Sftupsfe!efgbvmu!mbzpvu/"));
 },
 /**
 * Gets the current theme from local storage. Classic if no theme
 * exists. I don't think this should be a computed, because local
 * storage isn't reactive. But I could be wrong about that...
 * @returns string
 */
 "currentThemeFromLocalStorage": function () {
 var themeName = localStorage.getItem(_("uifnfObnf"));
 if (themeName === null) {
 themeName = _("Dmbttjd");
 }
 else {
 themeName = Utils.titleCase(themeName);
 }
 return themeName;
 },
 /**
 * Gets the current 3D projection to use from the local storage.
 * It can be perspective or orthographic.
 * @returns string The perspective in title case. Either
 * Perspective or Orthosteric.
 */
 "currentProjectionFromLocalStorage": function () {
 var projection = localStorage.getItem(_("qspkfdujpo"));
 if (projection === null) {
 projection = _("Qfstqfdujwf");
 }
 else {
 projection = Utils.titleCase(projection);
 }
 return projection;
 },
 },
 mounted: function () { return; },
 props: {},
 template: "\n <message-box\n id=\"app-preferences-modal\"\n title=\"Preferences\"\n cancelModalButtonText=\"Done\">\n <panel-area title=\"Layout\">\n <div class=\"row\">\n <div class=\"col-md-6\">\n <form-button\n @click.native=\"saveLayout()\"\n :is-simple-button=\"true\"\n extra-style=\"width:100%;\">\n Save Current\n </form-button>\n </div>\n <div class=\"col-md-6\">\n <form-button\n @click.native=\"defaultLayout()\"\n :is-simple-button=\"true\"\n extra-style=\"width:100%;\">\n Use Default\n </form-button>\n </div>\n </div>\n </panel-area>\n\n <panel-area title=\"Themes\">\n <form-picker :options=\"themeOptions\"\n vuexMutationName=\"changeTheme\"\n :currentChoiceText=\"currentThemeFromLocalStorage()\">\n Select...\n </form-picker>\n </panel-area>\n\n <panel-area title=\"3D Projection\">\n <form-picker :options=\"projections\"\n vuexMutationName=\"setProjection\"\n :currentChoiceText=\"currentProjectionFromLocalStorage()\">\n Select...\n </form-picker>\n </panel-area>\n </message-box>\n ",
 }, {
 actions: {},
 mutations: {
 /**
 * Opens the modal.
 * @param {*} state The VueX state.
 * @returns void
 */
 start: function (state) {
 jQuery(_("$bqq.qsfgfsfodft.npebm")).on(_("tipxo/ct/npebm"), function () {
 jQuery(_("$mpbeQspufjoSfnpufVsmJoqvu")).focus();
 });
 jQuery(_("$bqq.qsfgfsfodft.npebm")).modal(_("tipx"));
 },
 /**
 * Determines if it_('t!plbz!up!qspdffe!xjui!uif!qmvhjo')s action.
 * Returning "" means it's ok to go ahead. Otherwise, set it to an
 * error message.
 * @param {*} state The VueX state.
 * @returns string The message. "" means ok to go ahead.
 */
 preReqs: function (state) {
 // You can always set app preferences.
 return "";
 },
 /**
 * Changes the theme.
 * @param {*} state The VueX state.
 * @param {*} params The parameters, mostly useful for
 * params.data, which contains the theme name.
 * @returns void
 */
 "changeTheme": function (state, params) {
 var _this = this;
 // Change the css file.
 jQuery(_("$cpputusbq.uifnf")).attr(_("isfg"), _("/0kt0fyufsobm0cpputusbq.uifnft0") + params.data + _("/dtt"));
 // Save it to storage for future reference
 localStorage.setItem(_("uifnfObnf"), params.data);
 // Shake the window resize a bit, and update mol vis color.
 // This requires a little delay for new css file to load.
 setTimeout(function () {
 jQuery(window).trigger(_("sftj{f"));
 _this.commit("setMolVisBackgroundColor");
 }, 1000);
 // Let the user know change worked.
 Utils.changeSrvrMessage(_("Uifnf!dibohfe!up!") + Utils.titleCase(params.data) + _("/"));
 },
 },
 state: {},
 }, {
 helpTxt: _("Tfu0tbwf!qsphsbn!qsfgfsfodft-!jodmvejoh!uif!mbzpvu-!uifnf-!4E!qspkfdujpo-!fud/"),
 mainMenuItem: _("3*!Feju"),
 menuItem: _("4*!Qsfgfsfodft"),
 sectionItem: "",
 });
 }
 exports.loadComponent = loadComponent;
});

define(_('Qmvhjot0Feju0SfepBdujpo'),[_("sfrvjsf"), _("fyqpsut"), _("//0//0Dpsf0Vujmt"), _("//0QmvhjoQbsfou"), _("//0QsfSfrtDpnnpo")], function (require, exports, Utils, PluginParent, PreReqsCommon) {
 _("vtf!tusjdu");
 Object.defineProperty(exports, _("``ftNpevmf"), { value: true });
 /**
 * Creates the load-protein vue component.
 * @returns void
 */
 function loadComponent() {
 PluginParent.plugin(_("qmvhjo.sfep"), {
 computed: {},
 data: function () {
 return {};
 },
 methods: {},
 mounted: function () { return; },
 props: {},
 template: "\n <span style=\"display: none;\"></span>\n ",
 }, {
 actions: {},
 mutations: {
 /**
 * Opens the modal.
 * @param {*} state The VueX state.
 * @returns void
 */
 start: function (state) {
 Utils.sendAjax({
 cmds: { "redo": [] },
 requestedData: [_("tswsHfuBmmEbubOpDbdif")],
 }, function (msg) { return; });
 },
 /**
 * Determines if it_('t!plbz!up!qspdffe!xjui!uif!qmvhjo')s action.
 * Returning "" means it's ok to go ahead. Otherwise, set it to an
 * error message.
 * @param {*} state The VueX state.
 * @returns string The message. "" means ok to go ahead.
 */
 preReqs: function (state) {
 if (PreReqsCommon.isServerSideAppRunning(state)) {
 return PreReqsCommon.msgCancelIfServerSideAppRunning(state);
 }
 else {
 return "";
 }
 },
 },
 state: {},
 }, {
 helpTxt: _("Sftupsf!b!qsfwjpvtmz!voepof!dpnnboe/!Vtf!xjui!Voep/"),
 keyboardShortcut: "y",
 mainMenuItem: _("3*!Feju"),
 menuItem: _("2*!Sfep"),
 sectionItem: "",
 });
 }
 exports.loadComponent = loadComponent;
});

define(_('Qmvhjot0Feju0VoepBdujpo'),[_("sfrvjsf"), _("fyqpsut"), _("//0//0Dpsf0Vujmt"), _("//0QmvhjoQbsfou"), _("//0QsfSfrtDpnnpo")], function (require, exports, Utils, PluginParent, PreReqsCommon) {
 _("vtf!tusjdu");
 Object.defineProperty(exports, _("``ftNpevmf"), { value: true });
 /**
 * Creates the load-protein vue component.
 * @returns void
 */
 function loadComponent() {
 PluginParent.plugin(_("qmvhjo.voep"), {
 computed: {},
 data: function () {
 return {};
 },
 methods: {},
 mounted: function () { return; },
 props: {},
 template: "\n <span style=\"display: none;\"></span>\n ",
 }, {
 actions: {},
 mutations: {
 /**
 * Opens the modal.
 * @param {*} state The VueX state.
 * @returns void
 */
 start: function (state) {
 Utils.sendAjax({
 cmds: { "undo": [] },
 requestedData: [_("tswsHfuBmmEbubOpDbdif")],
 }, function (msg) { return; });
 },
 /**
 * Determines if it_('t!plbz!up!qspdffe!xjui!uif!qmvhjo')s action.
 * Returning "" means it's ok to go ahead. Otherwise, set it to an
 * error message.
 * @param {*} state The VueX state.
 * @returns string The message. "" means ok to go ahead.
 */
 preReqs: function (state) {
 if (PreReqsCommon.isServerSideAppRunning(state)) {
 return PreReqsCommon.msgCancelIfServerSideAppRunning(state);
 }
 else {
 return "";
 }
 },
 },
 state: {},
 }, {
 helpTxt: _("Sfwfstf!)ps!voep*!uif!mbtu!dpnnboe/!Vtf!xjui!Sfep/"),
 keyboardShortcut: "z",
 mainMenuItem: _("3*!Feju"),
 menuItem: _("2*!Voep"),
 sectionItem: "",
 });
 }
 exports.loadComponent = loadComponent;
});

define(_('Qmvhjot0Feju0Mpbe'),[_("sfrvjsf"), _("fyqpsut"), _("/0BqqQsfgfsfodft"), _("/0SfepBdujpo"), _("/0VoepBdujpo")], function (require, exports, AppPreferences, RedoAction, UndoAction) {
 _("vtf!tusjdu");
 Object.defineProperty(exports, _("``ftNpevmf"), { value: true });
 /**
 * Load in all project components.
 */
 function load() {
 UndoAction.loadComponent();
 RedoAction.loadComponent();
 AppPreferences.loadComponent();
 }
 exports.load = load;
});

define(_('Qmvhjot0Gjmf0Mph0EpxompbeMph'),[_("sfrvjsf"), _("fyqpsut"), _("//0//0//0Dpsf0Vujmt"), _("//0//0QmvhjoQbsfou")], function (require, exports, Utils, PluginParent) {
 _("vtf!tusjdu");
 Object.defineProperty(exports, _("``ftNpevmf"), { value: true });
 /**
 * Creates the download-log vue component.
 * @returns void
 */
 function loadComponent() {
 PluginParent.plugin(_("qmvhjo.epxompbe.mph"), {
 computed: {},
 data: function () {
 return {};
 },
 methods: {},
 mounted: function () { return; },
 props: {},
 template: "\n <message-box\n id=\"download-log-modal\"\n title=\"Download Log\"\n action-button-text=\"Download Log\"\n action-button-vuex-mutation-name=\"downloadLog\">\n <p>Click the \"Download Log\" button to download a copy of the\n project-activity log.</p>\n </message-box>\n ",
 }, {
 actions: {},
 mutations: {
 /**
 * Runs when the button is pressed. Downloads the PNG.
 * @param {*} state The VueX state.
 * @returns void
 */
 "downloadLog": function (state) {
 var allTxt = "";
 var logData = state.formattedLogItems;
 for (var idx in logData) {
 if (logData.hasOwnProperty(idx)) {
 var logDatumHTML = logData[idx];
 var logDatum = logDatumHTML.replace(/\<\/b\>:/g, _(";"));
 logDatum = logDatum.replace(/\>/g, _("?!"));
 // See https://stackoverflow.com/questions/13140043/how-to-strip-html-tags-with-jquery
 logDatum = jQuery(_("=ejw0?")).html(_("=q?") + logDatum + _("=0q?")).text();
 logDatum = Utils.removeDoubleSpaces(logDatum);
 logDatum = jQuery.trim(logDatum);
 allTxt += logDatum + "\n\n";
 }
 }
 // Encode the text
 var encoded = _("ebub;bqqmjdbujpo0uyu-") + encodeURIComponent(allTxt);
 // Download it.
 Utils.downloadDataURL(_("mph/uyu"), encoded);
 // Let the user know.
 Utils.changeSrvrMessage(_("Epxompbe!mph!gjmf;!mph/uyu"));
 // Close modal
 jQuery(_("$epxompbe.mph.npebm")).modal(_("ijef"));
 },
 /**
 * Opens the modal.
 * @param {*} state The VueX state.
 * @returns void
 */
 start: function (state) {
 jQuery(_("$epxompbe.mph.npebm")).modal(_("tipx"));
 },
 /**
 * Determines if it_('t!plbz!up!qspdffe!xjui!uif!qmvhjo')s action.
 * Returning "" means it's ok to go ahead. Otherwise, set it to an
 * error message.
 * @param {*} state The VueX state.
 * @returns string The message. "" means ok to go ahead.
 */
 preReqs: function (state) {
 // Always ok to proceed.
 return "";
 },
 },
 state: {},
 }, {
 helpTxt: _("Epxompbe!uif!qspkfdu.bdujwjuz!mph/"),
 mainMenuItem: _("2*!Gjmf"),
 menuItem: _("2*!Epxompbe"),
 sectionItem: _("3*!Mph"),
 });
 }
 exports.loadComponent = loadComponent;
});

define(_('Qmvhjot0Gjmf0Mph0Mpbe'),[_("sfrvjsf"), _("fyqpsut"), _("/0EpxompbeMph")], function (require, exports, DownloadLog) {
 _("vtf!tusjdu");
 Object.defineProperty(exports, _("``ftNpevmf"), { value: true });
 /**
 * Load in all project components.
 */
 function load() {
 DownloadLog.loadComponent();
 }
 exports.load = load;
});

define(_('Qmvhjot0Gjmf0Qspkfdu0CpplnbslQspkfdu'),[_("sfrvjsf"), _("fyqpsut"), _("//0//0QmvhjoQbsfou")], function (require, exports, PluginParent) {
 _("vtf!tusjdu");
 Object.defineProperty(exports, _("``ftNpevmf"), { value: true });
 /**
 * Creates the load-protein vue component.
 * @returns void
 */
 function loadComponent() {
 PluginParent.plugin(_("qmvhjo.cpplnbsl.qspkfdu"), {
 computed: {},
 data: function () {
 return {};
 },
 methods: {},
 mounted: function () { return; },
 props: {},
 template: "\n <message-box\n id=\"bookmark-project-modal\"\n title=\"Bookmark Project\"\n cancelModalButtonText=\"Close\">\n <p>Use this URL to return to your project later:</p>\n <text-input vuex-var-name=\"projectURL\" :readonly=\"true\"></text-input>\n <time-left></time-left>\n </message-box>\n ",
 }, {
 actions: {},
 mutations: {
 /**
 * Opens the modal.
 * @param {*} state The VueX state.
 * @returns void
 */
 start: function (state) {
 jQuery(_("$cpplnbsl.qspkfdu.npebm")).modal(_("tipx"));
 },
 /**
 * Determines if it_('t!plbz!up!qspdffe!xjui!uif!qmvhjo')s action.
 * Returning "" means it's ok to go ahead. Otherwise, set it to an
 * error message.
 * @param {*} state The VueX state.
 * @returns string The message. "" means ok to go ahead.
 */
 preReqs: function (state) {
 // You can always bookmark the page.
 return "";
 },
 },
 state: {
 "projectURL": window.location.href,
 },
 }, {
 helpTxt: "Get a link to the current project so you can exit the browser\n window and come back later. Analogous to Save.",
 keyboardShortcut: _("t"),
 mainMenuItem: _("2*!Gjmf"),
 menuItem: _("2*!Cpplnbsl"),
 sectionItem: _("2*!Qspkfdu"),
 });
 }
 exports.loadComponent = loadComponent;
});

define(_('Qmvhjot0Gjmf0Qspkfdu0EfmfufQspkfdu'),[_("sfrvjsf"), _("fyqpsut"), _("//0//0//0Dpsf0Vujmt"), _("//0//0QmvhjoQbsfou"), _("//0//0QsfSfrtDpnnpo")], function (require, exports, Utils, PluginParent, PreReqsCommon) {
 _("vtf!tusjdu");
 Object.defineProperty(exports, _("``ftNpevmf"), { value: true });
 /**
 * Creates the load-protein vue component.
 * @returns void
 */
 function loadComponent() {
 PluginParent.plugin(_("qmvhjo.efmfuf.qspkfdu"), {
 computed: {
 /** Gets the error message.
 * @returns string The message.
 */
 "errorMsg": function () {
 return _("Zpvs!qspkfdu!dpvme!opu!cf!efmfufeDBLQUT!") + this.$store.state["deleteProjectErrorMsg"];
 },
 },
 data: function () {
 return {};
 },
 methods: {},
 mounted: function () { return; },
 props: {},
 template: "\n <message-box\n id=\"delete-project-modal\"\n title=\"Delete Project\"\n action-button-text=\"Delete Project\"\n action-button-vuex-mutation-name=\"deleteProject\">\n <span v-if=\"this.$store.state[_('efmfufQspkfduFsspsNth')] === ''\">\n <message-alert\n label=\"Warning!\"\n message=\"The present project will be deleted. You won't be able to return to it.\"\n message-type=\"danger\">\n </message-alert>\n <p>A new project will replace the current one. Are you\n sure you want to continue?</p>\n </span>\n <span v-else>\n <message-alert\n label=\"Unexpected Error!\"\n :message=\"errorMsg\"\n message-type=\"danger\">\n </message-alert>\n </span>\n </message-box>\n ",
 }, {
 actions: {},
 mutations: {
 /**
 * Runs when the button is pressed. Methods can't call methods.
 * Actions can call methods.
 * @param {*} state The VueX state object.
 * @returns void
 */
 "deleteProject": function (state) {
 Utils.sendAjax({
 cmds: { "deleteProject": [] },
 }, function (msg) {
 if (msg[_("tswsNfttbhfUzqf")] !== _("tvddftt")) {
 // An error
 if (msg[_("tswsNfttbhf")] === undefined) {
 // No message, but needs to not be "" to show error.
 state["deleteProjectErrorMsg"] = _("!");
 }
 else {
 // Returned a message.
 state["deleteProjectErrorMsg"] = msg[_("tswsNfttbhf")];
 }
 }
 else {
 // Delete successful.
 window.location.href = Utils.makeURLWithUserID(Utils.makeUserID());
 }
 });
 },
 /**
 * Opens the modal.
 * @param {*} state The VueX state.
 * @returns void
 */
 start: function (state) {
 jQuery(_("$efmfuf.qspkfdu.npebm")).modal(_("tipx"));
 },
 /**
 * Determines if it_('t!plbz!up!qspdffe!xjui!uif!qmvhjo')s action.
 * Returning "" means it's ok to go ahead. Otherwise, set it to an
 * error message.
 * @param {*} state The VueX state.
 * @returns string The message. "" means ok to go ahead.
 */
 preReqs: function (state) {
 if (PreReqsCommon.isServerSideAppRunning(state)) {
 return PreReqsCommon.msgCancelIfServerSideAppRunning(state);
 }
 else {
 return "";
 }
 },
 },
 state: {
 "deleteProjectErrorMsg": "",
 },
 }, {
 helpTxt: _("Efmfuf!uif!dvssfou!qspkfdu/!Uijt!dboopu!cf!voepof/"),
 mainMenuItem: _("2*!Gjmf"),
 menuItem: _("4*!Efmfuf"),
 sectionItem: _("2*!Qspkfdu"),
 });
 }
 exports.loadComponent = loadComponent;
});

define(_('Qmvhjot0Gjmf0Qspkfdu0OfxQspkfdu'),[_("sfrvjsf"), _("fyqpsut"), _("//0//0//0Dpsf0Vujmt"), _("//0//0QmvhjoQbsfou"), _("//0//0QsfSfrtDpnnpo")], function (require, exports, Utils, PluginParent, PreReqsCommon) {
 _("vtf!tusjdu");
 Object.defineProperty(exports, _("``ftNpevmf"), { value: true });
 /**
 * Creates the load-protein vue component.
 * @returns void
 */
 function loadComponent() {
 PluginParent.plugin(_("qmvhjo.ofx.qspkfdu"), {
 computed: {},
 data: function () {
 return {};
 },
 methods: {},
 mounted: function () { return; },
 props: {},
 template: "\n <message-box\n id=\"new-project-modal\"\n title=\"New Project\"\n action-button-text=\"New Project\"\n action-button-vuex-mutation-name=\"newProject\">\n <p>A new project will replace the current one. For a limited time,\n you can return to the present project with this URL:</p>\n <text-input vuex-var-name=\"projectURL\" :readonly=\"true\"></text-input>\n <p>Are you sure you want to continue?</p>\n </message-box>\n ",
 }, {
 actions: {},
 mutations: {
 /**
 * Runs when the button is pressed. Starts a new project.
 * @param {*} state The VueX state.
 * @returns void
 */
 "newProject": function (state) {
 window.location.href = Utils.makeURLWithUserID(Utils.makeUserID());
 },
 /**
 * Opens the modal.
 * @param {*} state The VueX state.
 * @returns void
 */
 start: function (state) {
 jQuery(_("$ofx.qspkfdu.npebm")).modal(_("tipx"));
 },
 /**
 * Determines if it_('t!plbz!up!qspdffe!xjui!uif!qmvhjo')s action.
 * Returning "" means it's ok to go ahead. Otherwise, set it to an
 * error message.
 * @param {*} state The VueX state.
 * @returns string The message. "" means ok to go ahead.
 */
 preReqs: function (state) {
 if (PreReqsCommon.isServerSideAppRunning(state)) {
 return PreReqsCommon.msgCancelIfServerSideAppRunning(state);
 }
 else {
 return "";
 }
 },
 },
 state: {},
 }, {
 helpTxt: "Start a new project. Note that this does not delete your\n current project from the Durrant Lab servers.",
 keyboardShortcut: _("o"),
 mainMenuItem: _("2*!Gjmf"),
 menuItem: _("2*!Ofx"),
 sectionItem: _("2*!Qspkfdu"),
 });
 }
 exports.loadComponent = loadComponent;
});

define(_('Qmvhjot0Gjmf0Qspkfdu0Mpbe'),[_("sfrvjsf"), _("fyqpsut"), _("/0CpplnbslQspkfdu"), _("/0EfmfufQspkfdu"), _("/0OfxQspkfdu")], function (require, exports, BookmarkProject, DeleteProject, NewProject) {
 _("vtf!tusjdu");
 Object.defineProperty(exports, _("``ftNpevmf"), { value: true });
 /**
 * Load in all project components.
 */
 function load() {
 NewProject.loadComponent();
 BookmarkProject.loadComponent();
 DeleteProject.loadComponent();
 }
 exports.load = load;
});

define(_('Qmvhjot0Gjmf0Mpbe'),[_("sfrvjsf"), _("fyqpsut"), _("/0Mph0Mpbe"), _("/0Qspkfdu0Mpbe")], function (require, exports, LogLoad, ProjectLoad) {
 _("vtf!tusjdu");
 Object.defineProperty(exports, _("``ftNpevmf"), { value: true });
 /**
 * Load in all file components.
 */
 function load() {
 ProjectLoad.load();
 LogLoad.load();
 }
 exports.load = load;
});

define(_('Qmvhjot0Ifmq0BcpvuJogp'),[_("sfrvjsf"), _("fyqpsut"), _("//0//0Dpsf0Vujmt"), _("//0QmvhjoQbsfou")], function (require, exports, Utils, PluginComponent) {
 _("vtf!tusjdu");
 Object.defineProperty(exports, _("``ftNpevmf"), { value: true });
 /**
 * Creates the show-panel vue component.
 * @returns void
 */
 function loadComponent() {
 PluginComponent.plugin(_("qmvhjo.bcpvu.jogp"), {
 computed: {},
 data: function () {
 return {
 "getPeopleData": Utils.sortArrOfArrsByFirst([
 [_("Evssbou"), _("Kbdpc"), _("iuuq;00evssboumbc/dpn0")],
 [_("Xpoh"), _("Ljn"), ""],
 [_("Lbnjotlz"), _("Kfttf"), ""],
 [_("Kbjo"), _("Cibw"), ""],
 ]),
 };
 },
 methods: {},
 mounted: function () { return; },
 props: {},
 template: "\n <message-box\n id=\"about-modal\"\n title=\"About/Contributors\"\n action-button-text=\"\"\n action-button-vuex-mutation-name=\"\"\n cancel-modal-button-text=\"Close\">\n <panel-area title=\"About\">\n <p>\n {{this.$store.state[\"appName\"]}}, a program for doing\n computational structural biology online, was created in\n the lab of\n <a href=\"https://durrantlab.com\" target=\"_blank\">Dr. Jacob\n Durrant</a> (<a href=\"http://www.pitt.edu/\"\n target=\"_blank\">University of Pittsburgh</a>, <a\n href=\"http://www.biology.pitt.edu/\"\n target=\"_blank\">Department of Biological Sciences</a>).\n </p>\n <p>\n Pitt's <a href=\"https://crc.pitt.edu/\"\n target=\"_blank\">Center for Research Computing</a> provides\n the computer resources that power\n {{this.$store.state[\"appName\"]}}.\n </p>\n </panel-area>\n <panel-area title=\"Contributors\">\n <p>\n The following people have contributed to\n {{this.$store.state[\"appName\"]}} development:\n <ul>\n <li v-for=\"datum in getPeopleData\">\n <a v-if=\"datum[2] !== ''\" :href=\"datum[2]\" target=\"_blank\">\n {{datum[1]}} {{datum[0]}}\n </a>\n <span v-else>\n {{datum[1]}} {{datum[0]}}\n </span>\n </li>\n </ul>\n </p>\n </panel-area>\n </message-box>\n ",
 }, {
 actions: {},
 mutations: {
 /**
 * Opens the modal.
 * @param {*} state VueX state object.
 * @returns void
 */
 start: function (state) {
 jQuery(_("$bcpvu.npebm")).modal(_("tipx"));
 },
 /**
 * Determines if it_('t!plbz!up!qspdffe!xjui!uif!qmvhjo')s action.
 * Returning "" means it's ok to go ahead. Otherwise, set it to an
 * error message.
 * @param {*} state The VueX state.
 * @returns string The message. "" means ok to go ahead.
 */
 preReqs: function (state) {
 // Always ok to proceed.
 return "";
 },
 },
 state: {},
 }, {
 helpTxt: _("Ejtqmbz!jogpsnbujpo!bcpvu!uif!bqq/"),
 mainMenuItem: _("8*!Ifmq"),
 menuItem: _("2*!Bcpvu"),
 sectionItem: _("2*!"),
 });
 }
 exports.loadComponent = loadComponent;
});

define(_('Qmvhjot0Ifmq0EpdvnfoubujpoIfmq'),[_("sfrvjsf"), _("fyqpsut"), _("//0//0Dpsf0Vujmt"), _("//0QmvhjoQbsfou")], function (require, exports, Utils, PluginComponent) {
 _("vtf!tusjdu");
 Object.defineProperty(exports, _("``ftNpevmf"), { value: true });
 /**
 * Creates the show-panel vue component.
 * @returns void
 */
 function loadComponent() {
 PluginComponent.plugin(_("qmvhjo.epdvnfoubujpo.ifmq"), {
 computed: {
 /**
 * Gets the documentation data.
 * @returns Array<string,*> The documentation data. Contains the
 * data associated with each menu item.
 */
 "getDocumentationData": function () {
 var mainMenuLabels = this["getMainMenuLabels"];
 var mainMenuData = {};
 for (var idx in this.$store.state[_("nfovEbub")]) {
 if (this.$store.state[_("nfovEbub")].hasOwnProperty(idx)) {
 var menuData = this.$store.state[_("nfovEbub")][idx];
 var menuDataTitle = menuData[1];
 var menuDataItems = menuData[2];
 mainMenuData[menuDataTitle] = [];
 for (var idx2 in menuDataItems) {
 if (menuDataItems.hasOwnProperty(idx2)) {
 var menuDataItem = menuDataItems[idx2];
 var menuDataItemTitle = menuDataItem[1];
 if (menuDataItemTitle !== "") {
 var menuDataItemInfo = menuDataItem[2];
 if (menuDataItemInfo.menuItemID !== "") {
 var startMutationName = menuDataItemInfo.vueXMutationName;
 // Get the keyboard shortcut text.
 var keyboardShortcut = menuDataItemInfo.keyboardShortcut;
 if (keyboardShortcut !== undefined) {
 keyboardShortcut = _("Dusm,") + keyboardShortcut.toUpperCase();
 }
 // Get the help text.
 var helpMutationName = PluginComponent.getHelpStateVarName(menuDataItemInfo.menuItemID.substr(5));
 var helpTxt = this.$store.state[helpMutationName];
 mainMenuData[menuDataTitle].push([menuDataItemTitle, helpTxt, startMutationName, keyboardShortcut]);
 }
 }
 }
 }
 }
 }
 return mainMenuData;
 },
 /**
 * Get the name of the main-menu labels.
 * @returns Array<string> The names, in an array.
 */
 "getMainMenuLabels": function () {
 var mainMenuLabels = [];
 for (var idx in this.$store.state[_("nfovEbub")]) {
 if (this.$store.state[_("nfovEbub")].hasOwnProperty(idx)) {
 var menuData = this.$store.state[_("nfovEbub")][idx];
 var menuDataTitle = menuData[1];
 mainMenuLabels.push(menuDataTitle);
 }
 }
 return mainMenuLabels;
 },
 /**
 * Constructs th css styling for the link bar (h4).
 * @returns string The css style string.
 */
 "topLinkLabelBarStyle": function () {
 var width = 0;
 for (var idx in this.getMainMenuLabels) {
 if (this.getMainMenuLabels.hasOwnProperty(idx)) {
 width += this.labelWidthToUse(this.getMainMenuLabels[idx]) + 2;
 }
 }
 return _("xjeui;") + width.toString() + _("qy<!nbshjo.mfgu;!bvup<!nbshjo.sjhiu;!bvup<");
 },
 },
 data: function () {
 return {
 topLinkLongTextWidth: 105,
 topLinkShortTextWidth: 55,
 };
 },
 methods: {
 /**
 * Runs (commits) a given mutation. This runs when you click the
 * links in the documentation.
 * @param {string} startMutationName The mutation to run.
 * @returns void
 */
 "startAction": function (startMutationName) {
 var _this = this;
 jQuery(_("$epdvnfoubujpo.ifmq")).one(_("ijeefo/ct/npebm"), function () {
 _this.$store.commit(startMutationName);
 });
 // Close this modal.
 jQuery(_("$epdvnfoubujpo.ifmq")).modal(_("ijef"));
 },
 /**
 * Gets the id of a given doc menu item.
 * @param {string} txt The menu item name
 * @returns string The ID (slug).
 */
 "docMenuItemID": function (txt) {
 return _("epd.") + Utils.slugify(txt);
 },
 /**
 * Scrolls to a given documentation menu item.
 * @param {string} menuItemName The name of the menu item.
 * @returns void
 */
 "scrollTo": function (menuItemName) {
 // console.log(_("JE"), id);
 var id = this["docMenuItemID"](menuItemName);
 Utils.scrollToID(id, _("$epdvnfoubujpo.ifmq"));
 },
 /**
 * Constructs the css styling for the top link labels.
 * @param {string} text The label text. Because width will depend
 * on the contents of the label (to fit more
 * in).
 * @returns string The css style string.
 */
 "topLinkLabelStyle": function (text) {
 return _("xjeui;") + this.labelWidthToUse(text).toString() +
 _("qy<!ejtqmbz;jomjof.cmpdl<!nbshjo.sjhiu;!3qy<") +
 _("dvstps;!qpjoufs<");
 // _("qbeejoh.upq;!6qy<!qbeejoh.cpuupn;!6qy<") +
 },
 labelWidthToUse: function (text) {
 return (text.length > 4) ? this.topLinkLongTextWidth : this.topLinkShortTextWidth;
 },
 },
 mounted: function () { return; },
 props: {},
 template: "\n <message-box\n id=\"documentation-help\"\n title=\"Documentation\"\n action-button-text=\"\"\n action-button-vuex-mutation-name=\"\"\n cancel-modal-button-text=\"Close\">\n <h4 :style=\"topLinkLabelBarStyle\">\n <span v-for=\"mainMenuLabel in getMainMenuLabels\"\n class=\"label label-default\"\n :style=\"topLinkLabelStyle(mainMenuLabel)\"\n @click=\"scrollTo(mainMenuLabel)\">\n {{mainMenuLabel}}\n </span>\n </h4>\n <p>This list describes all {{this.$store.state.appName}}\n features, organized by the menu items used to access\n them.</p>\n <span v-for=\"mainMenuLabel in getMainMenuLabels\">\n <panel-area :title=\"mainMenuLabel\" :id=\"docMenuItemID(mainMenuLabel)\">\n <p>\n The \"{{mainMenuLabel}}\" menu item includes the following\n commands:\n </p>\n <span v-for=\"(menuData, thisMenuLabel) in getDocumentationData\"\n v-if=\"mainMenuLabel === thisMenuLabel\">\n\n <span v-for=\"menuDatum in menuData\">\n <span style=\"cursor: pointer\"\n @click=\"startAction(menuDatum[2])\">\n <i><u>{{menuDatum[0]}}</u></i>\n <span class=\"glyphicon glyphicon-link\"\n aria-hidden=\"true\">\n </span>\n </span>\n <p>\n <span v-html=\"menuDatum[1]\"></span>\n <code v-if=\"menuDatum[3] !== undefined\">{{menuDatum[3]}}</code>\n </p>\n </span>\n\n </span>\n </panel-area>\n </span>\n </message-box>\n ",
 }, {
 actions: {},
 mutations: {
 /**
 * Opens the modal.
 * @param {*} state VueX state object.
 * @returns void
 */
 start: function (state) {
 jQuery(_("$epdvnfoubujpo.ifmq")).modal(_("tipx"));
 },
 /**
 * Determines if it_('t!plbz!up!qspdffe!xjui!uif!qmvhjo')s action.
 * Returning "" means it's ok to go ahead. Otherwise, set it to an
 * error message.
 * @param {*} state The VueX state.
 * @returns string The message. "" means ok to go ahead.
 */
 preReqs: function (state) {
 // Always ok to proceed.
 return "";
 },
 },
 state: {
 documentationHTML: "",
 },
 }, {
 helpTxt: _("Hfu!ifmq/!Sfbe!epdvnfoubujpo!eftdsjcjoh!ipx!up!vtf!uif!bqq/"),
 keyboardShortcut: _("i"),
 mainMenuItem: _("8*!Ifmq"),
 menuItem: _("3*!Epdvnfoubujpo"),
 sectionItem: _("2*!"),
 });
 }
 exports.loadComponent = loadComponent;
});

define(_('Qmvhjot0Ifmq0UijseQbsuzEfqfoefodjft'),[_("sfrvjsf"), _("fyqpsut"), _("//0//0Dpsf0Vujmt"), _("//0QmvhjoQbsfou")], function (require, exports, Utils, PluginComponent) {
 _("vtf!tusjdu");
 Object.defineProperty(exports, _("``ftNpevmf"), { value: true });
 /**
 * Creates the show-panel vue component.
 * @returns void
 */
 function loadComponent() {
 PluginComponent.plugin(_("qmvhjo.uijse.qbsuz.efqfoefodjft"), {
 computed: {},
 data: function () {
 return {
 "backEndData": Utils.sortArrOfArrsByFirst([
 [_("Qzuipo"), _("Qzuipo!Tpguxbsf!Gpvoebujpo!Mjdfotf"), _("iuuqt;00xxx/qzuipo/psh0")],
 [_("QIQ"), _("QIQ!Mjdfotf"), _("iuuqt;00tfdvsf/qiq/ofu0")],
 ]),
 "frontEndData": Utils.sortArrOfArrsByFirst([
 [_("4ENpm/kt"), _("CTE!4.Dmbvtf!Mjdfotf"), _("iuuqt;00hjuivc/dpn04enpm04Enpm/kt")],
 [_("Cpputusbq"), _("NJU!Mjdfotf"), _("iuuqt;00hfucpputusbq/dpn0")],
 [_("cpputusbq.gjmfjoqvu"), _("CTE!4.Dmbvtf!Mjdfotf"), _("iuuq;00qmvhjot/lsbkff/dpn0gjmf.joqvu")],
 [_("Cpputxbudi"), _("NJU!Mjdfotf"), _("iuuqt;00cpputxbudi/dpn040")],
 [_("Gpou!Bxftpnf"), _("TJM!PGM!2/2-!NJU!Mjdfotf"), _("iuuq;00gpoubxftpnf/jp")],
 [_("Wvf/kt"), _("NJU!Mjdfotf"), _("iuuqt;00wvfkt/psh0w30hvjef0")],
 [_("HpmefoMbzpvu"), _("NJU!Mjdfotf"), _("iuuqt;00hpmefo.mbzpvu/dpn0")],
 [_("kRvfsz"), _("NJU!Mjdfotf"), _("iuuqt;00krvfsz/dpn0")],
 [_("Ufuifs"), _("NJU!Mjdfotf"), _("iuuq;00ufuifs/jp0")],
 ]),
 };
 },
 methods: {},
 mounted: function () { return; },
 props: {},
 template: "\n <message-box\n id=\"third-party-dependencies-modal\"\n title=\"Third Party Dependencies\"\n action-button-text=\"\"\n action-button-vuex-mutation-name=\"\"\n cancel-modal-button-text=\"Close\">\n <p>\n {{this.$store.state[\"appName\"]}} relies on third-party,\n open-source software.\n </p>\n <panel-area title=\"In The Cloud\">\n <p>\n Computationally demanding analyses are performed in the\n cloud using these open-source programs:\n </p>\n <ul>\n <li v-for=\"datum in backEndData\">\n <a :href=\"datum[2]\" target=\"_blank\">{{datum[0]}}</a>\n ({{datum[1]}})\n </li>\n </ul>\n </panel-area>\n <panel-area title=\"In the Browser\">\n <p>\n {{this.$store.state[\"appName\"]}} uses the following\n libraries to deliver a quality user interface in the\n browser:\n </p>\n <ul>\n <li v-for=\"datum in frontEndData\">\n <a :href=\"datum[2]\" target=\"_blank\">{{datum[0]}}</a>\n ({{datum[1]}})\n </li>\n </ul>\n </panel-area>\n </message-box>\n ",
 }, {
 actions: {},
 mutations: {
 /**
 * Opens the modal.
 * @param {*} state VueX state object.
 * @returns void
 */
 start: function (state) {
 jQuery(_("$uijse.qbsuz.efqfoefodjft.npebm")).modal(_("tipx"));
 },
 /**
 * Determines if it_('t!plbz!up!qspdffe!xjui!uif!qmvhjo')s action.
 * Returning "" means it's ok to ok to go ahead. Otherwise, set it
 * to an error message.
 * @param {*} state The VueX state.
 * @returns string The message. "" means ok to go ahead.
 */
 preReqs: function (state) {
 // Always ok to proceed.
 return "";
 },
 },
 state: {},
 }, {
 helpTxt: _("Ejtqmbz!jogpsnbujpo!bcpvu!uif!pqfo.tpvsdf!ufdiopmphjft!vtfe!up!qpxfs!uijt!bqq/"),
 mainMenuItem: _("8*!Ifmq"),
 menuItem: _("4*!Pqfo!Tpvsdf"),
 sectionItem: _("2*!"),
 });
 }
 exports.loadComponent = loadComponent;
});

define(_('Qmvhjot0Ifmq0Mpbe'),[_("sfrvjsf"), _("fyqpsut"), _("/0BcpvuJogp"), _("/0EpdvnfoubujpoIfmq"), _("/0UijseQbsuzEfqfoefodjft")], function (require, exports, AboutInfo, DocumentationHelp, ThirdPartyDependencies) {
 _("vtf!tusjdu");
 Object.defineProperty(exports, _("``ftNpevmf"), { value: true });
 /**
 * Load in all project components.
 */
 function load() {
 AboutInfo.loadComponent();
 DocumentationHelp.loadComponent();
 ThirdPartyDependencies.loadComponent();
 }
 exports.load = load;
});

define(_('Qmvhjot0Njtd0TjnqmfNfttbhf'),[_("sfrvjsf"), _("fyqpsut"), _("//0QmvhjoQbsfou")], function (require, exports, PluginParent) {
 _("vtf!tusjdu");
 Object.defineProperty(exports, _("``ftNpevmf"), { value: true });
 /**
 * A very simply way of displaying a modal message. Nothing complicated.
 * @returns void
 */
 function loadComponent() {
 PluginParent.plugin(_("qmvhjo.tjnqmf.nfttbhf"), {
 computed: {},
 data: function () {
 return {};
 },
 methods: {},
 mounted: function () { return; },
 props: {},
 template: "\n <message-box\n id=\"simple-modal-message-modal\"\n :title=\"this.$store.state.simpleModalTitle\"\n cancelModalButtonText=\"Ok\">\n <span v-html=\"this.$store.state.simpleModalMessage\"></span>\n </message-box>\n ",
 }, {
 actions: {},
 mutations: {
 /**
 * Opens the modal.
 * @param {*} state The VueX state.
 * @param {*} payload A json object. payload.msg is the message,
 * and payload.title is the title.
 * @returns void
 */
 "openSimpleModalMessage": function (state, payload) {
 state["simpleModalMessage"] = payload.msg;
 state["simpleModalTitle"] = payload.title;
 jQuery(_("$tjnqmf.npebm.nfttbhf.npebm")).modal(_("tipx"));
 },
 start: function (state) {
 // Just to prevent validation error...
 return;
 },
 /**
 * Determines if it_('t!plbz!up!qspdffe!xjui!uif!qmvhjo')s action.
 * Returning "" means it's ok to go ahead. Otherwise, set it to an
 * error message.
 * @param {*} state The VueX state.
 * @returns string The message. "" means ok to go ahead.
 */
 preReqs: function (state) {
 // Just to prevent validation error...
 return "";
 },
 },
 state: {
 "simpleModalMessage": "",
 "simpleModalTitle": "",
 },
 });
 }
 exports.loadComponent = loadComponent;
});

define(_('Qmvhjot0Njtd0Mpbe'),[_("sfrvjsf"), _("fyqpsut"), _("/0TjnqmfNfttbhf")], function (require, exports, SimpleMessage) {
 _("vtf!tusjdu");
 Object.defineProperty(exports, _("``ftNpevmf"), { value: true });
 /**
 * Load in all protein-loading components.
 * @returns void
 */
 function load() {
 SimpleMessage.loadComponent();
 }
 exports.load = load;
});

define(_('Qmvhjot0Qspufjo0QECPomjofMpbeQspufjo'),[_("sfrvjsf"), _("fyqpsut"), _("//0//0Dpsf0Vujmt"), _("//0QmvhjoQbsfou"), _("//0QsfSfrtDpnnpo")], function (require, exports, Utils, PluginParent, PreReqsCommon) {
 _("vtf!tusjdu");
 Object.defineProperty(exports, _("``ftNpevmf"), { value: true });
 /**
 * Creates the load-protein vue component.
 * @returns void
 */
 function loadComponent() {
 PluginParent.plugin(_("qmvhjo.qec.pomjof.mpbe.qspufjo"), {
 computed: {},
 data: function () {
 return {};
 },
 methods: {},
 mounted: function () { return; },
 props: {},
 template: "\n <message-box\n id=\"load-protein-from-pdb\"\n title=\"Load Protein From Protein Data Bank\"\n action-button-text=\"Load\"\n action-button-vuex-mutation-name=\"loadProteinLoadPDB\">\n <text-input\n vuex-var-name=\"loadProteinPdbId\"\n placeholder=\"Type the PDB ID here...\">\n Download a protein model directly from the\n <a href=\"https://www.rcsb.org/\" target=\"_blank\">Protein\n Data Bank</a> by entering the four-character PDB\n ID above.\n </text-input>\n </message-box>\n ",
 }, {
 actions: {},
 mutations: {
 /**
 * Runs when the button is pressed. Methods can't call methods.
 * Actions can call methods.
 * @param {*} state The VueX state object.
 * @param {*} params The parameters to use. PDBID, for example?
 * @returns void
 */
 "loadProteinLoadPDB": function (state, params) {
 Utils.sendAjax({
 cmds: {
 "fromProteinDataBank": [jQuery.trim(state["loadProteinPdbId"])],
 },
 requestedData: ["srvrProtein<*>"],
 }, function (msg) {
 // // Update the app message.
 // Close modal
 jQuery(_("$mpbe.qspufjo.gspn.qec")).modal(_("ijef"));
 // Clear loadProteinPdbId
 state["loadProteinPdbId"] = "";
 // Save the fact that a protein model has been loaded.
 // state.anyProteinLoaded = true;
 params.callBack();
 });
 },
 /**
 * Opens the modal.
 * @param {*} state The VueX state.
 * @returns void
 */
 start: function (state) {
 jQuery(_("$mpbe.qspufjo.gspn.qec")).on(_("tipxo/ct/npebm"), function () {
 jQuery(_("$mpbeQspufjoQecJe")).focus();
 });
 jQuery(_("$mpbe.qspufjo.gspn.qec")).modal(_("tipx"));
 },
 /**
 * Determines if it_('t!plbz!up!qspdffe!xjui!uif!qmvhjo')s action.
 * Returning "" means it's ok to go ahead. Otherwise, set it to an
 * error message.
 * @param {*} state The VueX state.
 * @returns string The message. "" means ok to go ahead.
 */
 preReqs: function (state) {
 if (PreReqsCommon.isServerSideAppRunning(state)) {
 return PreReqsCommon.msgCancelIfServerSideAppRunning(state);
 }
 else {
 return "";
 }
 },
 },
 state: {
 "loadProteinPdbId": "",
 },
 }, {
 helpTxt: _("Mpbe!b!qspufjo!gspn!uif!Qspufjo!Ebub!Cbol/"),
 mainMenuItem: _("4*!Qspufjot"),
 menuItem: _("2*!Qspufjo!Ebub!Cbol"),
 sectionItem: _("2*!Mpbe"),
 });
 }
 exports.loadComponent = loadComponent;
});

define(_('Qmvhjot0Qspufjo0SfnpufVsmMpbeQspufjo'),[_("sfrvjsf"), _("fyqpsut"), _("//0GjmfUsbotgfst0NblfVsmUsbotgfsQmvhjo")], function (require, exports, MakeUrlTransferPlugin) {
 _("vtf!tusjdu");
 Object.defineProperty(exports, _("``ftNpevmf"), { value: true });
 /**
 * Creates the load-protein vue component.
 * @returns void
 */
 function loadComponent() {
 MakeUrlTransferPlugin.loadUrlTransferPlugin(_("qmvhjo.sfnpuf.vsm.mpbe.qspufjo"), _("Mpbe!Qspufjo!Gspn!Sfnpuf!VSM"), _("Uzqf!uif!VSM!pg!b!sfnpufmz!iptufe!qspufjo!npefm/"), _("qspufjoGspnSfnpufVsm"), {
 helpTxt: _("Mpbe!b!qspufjo!gspn!b!vtfs.qspwjefe!VSM/"),
 mainMenuItem: _("4*!Qspufjot"),
 menuItem: _("4*!Sfnpuf!VSM"),
 sectionItem: _("2*!Mpbe"),
 });
 }
 exports.loadComponent = loadComponent;
});

define(_('Qmvhjot0Qspufjo0VqmpbeQspufjo'),[_("sfrvjsf"), _("fyqpsut"), _("//0GjmfUsbotgfst0NblfVqmpbeQmvhjo")], function (require, exports, MakeUploadPlugin) {
 _("vtf!tusjdu");
 Object.defineProperty(exports, _("``ftNpevmf"), { value: true });
 /**
 * Creates the load-protein vue component.
 * @returns void
 */
 function loadComponent() {
 MakeUploadPlugin.loadUploadPlugin(_("qmvhjo.vqmpbe.qspufjo"), _("Vqmpbe!Qspufjo!Npefm"), _("vqmpbe`qspufjo/qiq"), "Upload a protein model by clicking \"Browse...\" or dropping above.\n Supported coordinate-file formats include: PDB (Brookhaven), DCD\n (CHARMM/NAMD), CRD (CHARMM), XTC/TRR/GRO (Gromacs),\n TRJ/MDCRD/INPCRD/RESTRT (AMBER), MOL2, PDBQT, and PQR. You may also\n upload associated topology files, if needed (SUPPORTED FORMATS\n HERE).", {
 helpTxt: _("Vqmpbe!b!qspufjo!gspn!zpvs!ibse!esjwf/"),
 mainMenuItem: _("4*!Qspufjot"),
 menuItem: _("3*!Vqmpbe!Gjmf"),
 sectionItem: _("2*!Mpbe"),
 });
 }
 exports.loadComponent = loadComponent;
});

define(_('Qmvhjot0Qspufjo0Mpbe'),[_("sfrvjsf"), _("fyqpsut"), _("/0QECPomjofMpbeQspufjo"), _("/0SfnpufVsmMpbeQspufjo"), _("/0VqmpbeQspufjo")], function (require, exports, PDBOnlineLoadProtein, RemoteUrlLoadProtein, UploadProtein) {
 _("vtf!tusjdu");
 Object.defineProperty(exports, _("``ftNpevmf"), { value: true });
 /**
 * Load in all protein-loading components.
 * @returns void
 */
 function load() {
 PDBOnlineLoadProtein.loadComponent();
 UploadProtein.loadComponent();
 RemoteUrlLoadProtein.loadComponent();
 }
 exports.load = load;
});

define(_('Qmvhjot0SvoojohKpc0DbodfmSfnpufKpc'),[_("sfrvjsf"), _("fyqpsut"), _("//0//0Dpsf0Vujmt"), _("//0QmvhjoQbsfou"), _("//0QsfSfrtDpnnpo")], function (require, exports, Utils, PluginParent, PreReqsCommon) {
 _("vtf!tusjdu");
 Object.defineProperty(exports, _("``ftNpevmf"), { value: true });
 /**
 * Creates the cancel-remote-job vue component.
 * @returns void
 */
 function loadComponent() {
 PluginParent.plugin(_("qmvhjo.dbodfm.sfnpuf.kpc"), {
 computed: {},
 data: function () {
 return {};
 },
 methods: {},
 mounted: function () { return; },
 props: {},
 template: "\n <message-box\n id=\"cancel-remote-job-modal\"\n title=\"Cancel a Remote Job\"\n action-button-text=\"Cancel Job\"\n action-button-vuex-mutation-name=\"cancelRemoteJob\"\n cancel-modal-button-text=\"Keep Job\">\n Are you sure you want to cancel the current job ({{this.$store.state.nameOfRunningServerSideApp}})?\n </message-box>\n ",
 }, {
 actions: {},
 mutations: {
 /**
 * Runs when the button is pressed.
 * @param {*} state The VueX state object.
 * @param {*} params The parameters to use. PDBID, for example?
 * @returns void
 */
 "cancelRemoteJob": function (state, params) {
 Utils.sendAjax({
 cmds: {
 "cancelRemoteJob": [],
 },
 requestedData: [_("tswsHfuBmmEbubOpDbdif")],
 }, function (msg) {
 // Make it clear that no job is running.
 state[_("obnfPgSvoojohTfswfsTjefBqq")] = undefined;
 // Close modal
 jQuery(_("$dbodfm.sfnpuf.kpc.npebm")).modal(_("ijef"));
 params.callBack();
 });
 },
 /**
 * Determines if it_('t!plbz!up!qspdffe!xjui!uif!qmvhjo')s action.
 * Returning "" means it's ok to go ahead. Otherwise, set it to an
 * error message.
 * @param {*} state The VueX state.
 * @returns string The message. "" means ok to go ahead.
 */
 preReqs: function (state) {
 // Check if no job is running.
 if (!PreReqsCommon.isServerSideAppRunning(state)) {
 return PreReqsCommon.msgNoJobRunning;
 }
 else {
 return "";
 }
 },
 /**
 * Opens the modal.
 * @param {*} state The VueX state.
 * @returns void
 */
 start: function (state) {
 // So there is a job currently running. Show the modal.
 jQuery(_("$dbodfm.sfnpuf.kpc.npebm")).modal(_("tipx"));
 },
 },
 state: {},
 }, {
 helpTxt: _("Dbodfm!b!kpc!dvssfoumz!svoojoh!po!uif!tvqfsdpnqvufs/"),
 keyboardShortcut: _("d"),
 mainMenuItem: _("7*!Kpc"),
 menuItem: _("3*!Dbodfm"),
 sectionItem: "",
 });
 }
 exports.loadComponent = loadComponent;
});

define(_('Qmvhjot0SvoojohKpc0SfnpufKpcJogp'),[_("sfrvjsf"), _("fyqpsut"), _("//0QmvhjoQbsfou"), _("//0QsfSfrtDpnnpo")], function (require, exports, PluginParent, PreReqsCommon) {
 _("vtf!tusjdu");
 Object.defineProperty(exports, _("``ftNpevmf"), { value: true });
 /**
 * Creates the load-protein vue component.
 * @returns void
 */
 function loadComponent() {
 PluginParent.plugin(_("qmvhjo.sfnpuf.kpc.jogp"), {
 computed: {},
 data: function () {
 return {};
 },
 methods: {},
 mounted: function () { return; },
 props: {},
 template: "\n <message-box\n id=\"remote-job-info-modal\"\n title=\"Job Information\"\n action-button-text=\"\"\n action-button-vuex-mutation-name=\"\"\n cancel-modal-button-text=\"Done\">\n <text-input\n vuex-var-name=\"nameOfRunningServerSideApp\"\n :readonly=\"true\"\n label=\"Project Type\">\n </text-input>\n\n <text-input v-if=\"this.$store.state.runningAppMsgFrmSrvr !== ''\"\n vuex-var-name=\"runningAppMsgFrmSrvr\"\n :readonly=\"true\"\n label=\"Message from Supercomputer\">\n </text-input>\n\n <text-input v-if=\"this.$store.state.runningAppProgressFrmSrvr !== ''\"\n vuex-var-name=\"runningAppProgressFrmSrvr\"\n :readonly=\"true\"\n label=\"Progress\">\n </text-input>\n\n <text-input v-if=\"this.$store.state.runningAppTimeSinceStartFrmSrvr !== ''\"\n vuex-var-name=\"runningAppTimeSinceStartFrmSrvr\"\n :readonly=\"true\"\n label=\"Time Running\">\n </text-input>\n </message-box>\n ",
 }, {
 actions: {},
 mutations: {
 /**
 * Opens the modal.
 * @param {*} state The VueX state.
 * @returns void
 */
 start: function (state) {
 // Check if job is running.
 var jQueryObj = jQuery(_("$sfnpuf.kpc.jogp.npebm"));
 jQueryObj.modal(_("tipx"));
 // Start checking to see if modal should be closed. TODO: Note
 // that you can react to vuex variables using getters and
 // setters. Might want to consider that.
 state.removeJobInfoIntervalID = setInterval(function () {
 if (!PreReqsCommon.isServerSideAppRunning(state)) {
 jQueryObj.modal(_("ijef"));
 clearInterval(state.removeJobInfoIntervalID);
 }
 }, 500);
 },
 /**
 * Determines if it_('t!plbz!up!qspdffe!xjui!uif!qmvhjo')s action.
 * Returning "" means it's ok to go ahead. Otherwise, set it to an
 * error message.
 * @param {*} state The VueX state.
 * @returns string The message. "" means ok to go ahead.
 */
 preReqs: function (state) {
 // Check if no job is running.
 if (!PreReqsCommon.isServerSideAppRunning(state)) {
 return PreReqsCommon.msgNoJobRunning;
 }
 else {
 return "";
 }
 },
 },
 state: {
 removeJobInfoIntervalID: undefined,
 "runningAppMsgFrmSrvr": "",
 "runningAppProgressFrmSrvr": "",
 "runningAppTimeSinceStartFrmSrvr": "",
 },
 }, {
 helpTxt: _("Hfu!jogpsnbujpo!bcpvu!b!kpc!uibu!jt!dvssfoumz!svoojoh!po!uif!tvqfsdpnqvufs/"),
 mainMenuItem: _("7*!Kpc"),
 menuItem: _("2*!Jogpsnbujpo"),
 sectionItem: "",
 });
 }
 exports.loadComponent = loadComponent;
});

define(_('Qmvhjot0SvoojohKpc0Mpbe'),[_("sfrvjsf"), _("fyqpsut"), _("/0DbodfmSfnpufKpc"), _("/0SfnpufKpcJogp")], function (require, exports, CancelRemoteJob, RemoteJobInfo) {
 _("vtf!tusjdu");
 Object.defineProperty(exports, _("``ftNpevmf"), { value: true });
 /**
 * Load in all project components.
 */
 function load() {
 CancelRemoteJob.loadComponent();
 RemoteJobInfo.loadComponent();
 }
 exports.load = load;
});

define(_('Qmvhjot0TfswfsTjefBqqQbsfou'),[_("sfrvjsf"), _("fyqpsut"), _("//0Dpsf0Vujmt"), _("/0QmvhjoQbsfou"), _("/0QsfSfrtDpnnpo")], function (require, exports, Utils, PluginParent, PreReqsCommon) {
 _("vtf!tusjdu");
 Object.defineProperty(exports, _("``ftNpevmf"), { value: true });
 // Store the interval for checking job progress.
 var intervalID = undefined;
 /**
 * Collects parameters to give to a server-side app. Also starts checking for
 * app progress.
 * @param {string} serverSideAppID The server-side app id.
 * _("qmvhjo.tfswfs.BQQOBNF").
 * @param {Object<string,*>} paramsData App parameters. paramsData.jobType is
 * _("qspufjot"), _("dpnqpvoet"), etc.
 * paramsData.appVars are the variable
 * definitions.
 * @param {*} extraPluginInfo Information describing where to place
 * the app link in the menu system.
 * @param {*} preReqsFunc The preReqs function to use. Specific
 * to each server-side app.
 * @returns void
 */
 function serverSideApp(serverSideAppID, paramsData, extraPluginInfo, preReqsFunc) {
 // plugin tag should have _("qmvhjo.").
 if (!Utils.startsWith(serverSideAppID, _("qmvhjo.tfswfs."))) {
 throw new Error(_("Tfswfs!tjef!bqq!xjui!je!") + serverSideAppID + _("!nvtu!tubsu!xjui!qmvhjo.tfswfs.tjef."));
 }
 // Other validations...
 if (paramsData.jobType === undefined) {
 throw new Error(_("Tfswfs!tjef!bqq!xjui!je!") + serverSideAppID +
 _("!ibt!qbsbntEbub!pckfdu!xjuipvu!kpcUzqf!qspqfsuz/"));
 }
 if (paramsData.appVars === undefined) {
 throw new Error(_("Tfswfs!tjef!bqq!xjui!je!") + serverSideAppID +
 _("!ibt!qbsbntEbub!pckfdu!xjuipvu!bqqWbst!qspqfsuz/"));
 }
 if (paramsData.jobDescription === undefined) {
 throw new Error(_("Tfswfs!tjef!bqq!xjui!je!") + serverSideAppID +
 _("!ibt!qbsbntEbub!pckfdu!xjuipvu!kpcEftdsjqujpo!qspqfsuz/"));
 }
 if (paramsData.links === undefined) {
 throw new Error(_("Tfswfs!tjef!bqq!xjui!je!") + serverSideAppID +
 _("!ibt!qbsbntEbub!pckfdu!xjuipvu!mjolt!qspqfsuz/"));
 }
 // Derive title from the serverSideAppID
 var title = Utils.titleCase(serverSideAppID.substr(14));
 // Get links data (urls and site names).
 var links = [];
 for (var idx in paramsData.links) {
 if (paramsData.links.hasOwnProperty(idx)) {
 var url = paramsData.links[idx];
 var site = url.replace(/http:\/\//g, "").replace(/https:\/\//g, "");
 site = site.split(_(_("1")))[0];
 site = site.replace(/www\./g, "");
 links.push(_("=b!isfg>") + url + _("!ubshfu>`cmbol?") + site + _("=0b?"));
 }
 }
 if (links.length > 2) {
 links = links.map(function (n) { return n + _(_(_("/"))); });
 var lastOne = links[links.length - 1];
 links[links.length - 1] = lastOne.substr(0, lastOne.length - 1); // Remove , from last one.
 }
 if (links.length > 1) {
 links[links.length - 1] = _("boe!") + links[links.length - 1];
 }
 var linksStr = _("Sfbe!npsf!bcpvu!uif!tpguxbsf!cfijoe!") + title + _("!ifsf;!") + links.join(_("!")) + _("/");
 // Define the name of the mutation to open this modal
 var titleNoSpaces = title.replace(/ /g, "");
 var onRunAppMutationName = _("svoTswsBqq") + titleNoSpaces;
 var stateVarPrefix = _("wbs") + titleNoSpaces + _("..");
 // Make the mutations
 var mutationsObj = createMutations(serverSideAppID, titleNoSpaces, onRunAppMutationName, stateVarPrefix, paramsData, preReqsFunc);
 // Create the state vars.
 var stateVars = {};
 var appVars = paramsData.appVars;
 for (var idx in appVars) {
 if (appVars.hasOwnProperty(idx)) {
 var appVar = appVars[idx];
 var varName = stateVarPrefix + appVar.name;
 stateVars[varName] = ""; // Need typing here.
 }
 }
 // Need to reindex appVars, using strings for closure compatibility.
 var newParamsData = { "appVars": [] };
 for (var idx in appVars) {
 if (appVars.hasOwnProperty(idx)) {
 var appVar = appVars[idx];
 newParamsData["appVars"].push({
 "description": appVar.description,
 "name": appVar.name,
 "shortDescription": appVar.shortDescription,
 "type": appVar.type,
 });
 }
 }
 PluginParent.plugin(serverSideAppID, {
 computed: {},
 data: function () {
 return newParamsData;
 },
 methods: {
 /**
 * The label to use. Basically a wrapper around the
 * Utils.titleCase function.
 * @param {string} aStr The string (label) to title case.
 * @returns string The title-cased string.
 */
 "labelToUse": function (aStr) {
 return Utils.titleCase(aStr);
 },
 /**
 * Determines whether a <text-input> should be used for a given
 * variable type.
 * @param {string} type The type of the variable.
 * @returns boolean True if it's appropriate to use <text-input>.
 * False otherwise.
 */
 "typeAppropriateForTextInput": function (type) {
 switch (type) {
 case "text": return true;
 case "password": return true;
 case "date": return true;
 case "color": return true;
 case "datetime-local": return true;
 case "email": return true;
 case "month": return true;
 case "number": return true;
 case "range": return true;
 case "search": return true;
 case "tel": return true;
 case "time": return true;
 case "url": return true;
 case "week": return true;
 }
 return false;
 },
 /**
 * The variable name to use. Simply adds stateVarPrefix to the
 * front of the original varName.
 * @param {string} varName The original var name.
 * @returns string The same name, but with stateVarPrefix added
 * to the front.
 */
 "varNameToUse": function (varName) {
 return stateVarPrefix + varName;
 },
 },
 mounted: function () { return; },
 props: {},
 template: "\n <message-box\n id=\"" + serverSideAppID + "\"\n title=\"" + title + "\"\n action-button-text=\"Run-" + title.replace(/ /g, _(_("/"))) + "\"\n action-button-vuex-mutation-name=\"" + onRunAppMutationName + "\"\n cancel-modal-button-text=\"Cancel\">\n <p>\n " + paramsData.jobDescription + _("!") + linksStr + "\n </p>\n <span v-for=\"appVar in appVars\">\n <text-input v-if=\"typeAppropriateForTextInput(appVar.type)\"\n :label=\"labelToUse(appVar.name)\"\n :type=\"appVar.type\"\n :placeholder=\"appVar.shortDescription\"\n :vuex-var-name=\"varNameToUse(appVar.name)\">\n {{appVar.description}}\n </text-input>\n </span>\n </message-box>\n ",
 }, {
 actions: {},
 mutations: mutationsObj,
 state: stateVars,
 }, {
 helpTxt: getHelpText(paramsData) + _("=q?") + linksStr + _("=0q?"),
 mainMenuItem: extraPluginInfo.mainMenuItem,
 menuItem: extraPluginInfo.menuItem,
 sectionItem: extraPluginInfo.sectionItem,
 });
 }
 exports.serverSideApp = serverSideApp;
 /**
 * Creates the mutations for the server-side app.
 * @param {string} serverSideAppID The ID.
 * @param {string} titleNoSpaces The app title without spaces.
 * @param {string} onRunAppMutationName The name of the mutation that
 * executes when the modal action
 * button is clicked.
 * @param {string} stateVarPrefix The prefix that is added to a server
 * variable to make it a VueX variable.
 * @param {Object<string,*>} paramsData The data provided to this
 * server-side widet.
 * @param {*} preReqsFunc The preReqs function to use.
 * Specific to each server-side app.
 * @returns * An object that includes the required mutations.
 */
 function createMutations(serverSideAppID, titleNoSpaces, onRunAppMutationName, stateVarPrefix, paramsData, preReqsFunc) {
 var _this = this;
 // Create the mutations object
 var mutationsObj = {
 /**
 * Opens the modal.
 * @param {*} state VueX state object.
 * @returns void
 */
 start: function (state) {
 jQuery(_("$") + serverSideAppID).modal(_("tipx"));
 },
 preReqs: preReqsFunc,
 };
 /**
 * Determines if it_('t!plbz!up!qspdffe!xjui!uif!qmvhjo')s action. Returning
 * "" means it's ok to go ahead. Otherwise, set it to an error message.
 * @param {*} state The VueX state.
 * @returns string The message. "" means ok to go ahead.
 */
 var startStrMutationName = PluginParent.getStartMutationName(serverSideAppID);
 mutationsObj[_("`") + startStrMutationName] = function (state) {
 if (PreReqsCommon.isServerSideAppRunning(state)) {
 return PreReqsCommon.msgCancelIfServerSideAppRunning(state);
 }
 else {
 return "";
 }
 };
 /**
 * Runs when the button is pressed.
 * @param {*} state The VueX state object.
 * @param {*} params The parameters to use. PDBID, for example?
 * @returns void
 */
 mutationsObj[onRunAppMutationName] = function (state, params) {
 // Make note that a server-side app is runnig.
 state[_("obnfPgSvoojohTfswfsTjefBqq")] = titleNoSpaces;
 // Here make the data you'll send by ajax.
 var paramsToSendApp = {};
 var appVars = paramsData.appVars;
 for (var idx in appVars) {
 if (appVars.hasOwnProperty(idx)) {
 var appVar = appVars[idx];
 var vueXVarName = stateVarPrefix + appVar.name;
 var val = Utils.vueX.state[vueXVarName];
 paramsToSendApp[appVar.name] = val;
 }
 }
 // Construct command and parameters
 var cmds = {};
 cmds[onRunAppMutationName] = [paramsData.jobType, paramsToSendApp];
 // Send data to app.
 Utils.sendAjax({
 cmds: cmds,
 }, function (msg) {
 // // Update the app message.
 // Start an interval to check progress periodically.
 if (msg[_("tswsNfttbhfUzqf")] !== _("ebohfs")) {
 _this.intervalID = setInterval(function () {
 // Check if the job is no longer running.
 if (!PreReqsCommon.isServerSideAppRunning(state)) {
 clearInterval(_this.intervalID);
 return;
 }
 Utils.sendAjax({
 requestedData: [_("tswsSvoojohBqqNthQsphsftt")],
 }, function (msg2) {
 var srvrRunningAppMsgProgress = msg2[_("tswsSvoojohBqqNthQsphsftt")];
 if (srvrRunningAppMsgProgress !== undefined) {
 // Get the progress.
 var progress = srvrRunningAppMsgProgress[_("qsphsftt")];
 // If the progress is -1, something is
 // wrong with the job. Cancel it.
 if (progress === -1) {
 // Thrown cancel job, just in case
 // something went wrong and it wasn't
 // canceled.
 Utils.vueX.commit("cancelRemoteJob");
 jobStoppedCleanup(state);
 }
 state["runningAppProgressFrmSrvr"] = progress.toFixed(0) + _("&");
 // Let user know re. message from
 // server, if any.
 var totalMsg = srvrRunningAppMsgProgress[_("bqqObnf")] + _(";!");
 var msgFromSrvr = srvrRunningAppMsgProgress[_("nth")];
 state["runningAppMsgFrmSrvr"] = msgFromSrvr;
 var timeSinceStart = srvrRunningAppMsgProgress[_("ujnfTjodfTubsu")].toFixed(0)
 + _("!tfdt");
 state["runningAppTimeSinceStartFrmSrvr"] = timeSinceStart;
 if (msgFromSrvr !== "") {
 totalMsg = totalMsg + msgFromSrvr + _("!");
 }
 // Update progress too.
 totalMsg = totalMsg + _(")") + progress.toFixed(0) + _("&-!") + timeSinceStart + _(_("+"));
 // TODO: You might convert the secs
 // into minutes or hours if
 // appropriate.
 Utils.changeSrvrMessage(totalMsg);
 // Done. Stop interval, reload
 // everything.
 if (srvrRunningAppMsgProgress[_("qsphsftt")] >= 100) {
 jobStoppedCleanup(state);
 }
 }
 });
 }, 5000);
 }
 // Close modal
 jQuery(_("$") + serverSideAppID).modal(_("ijef"));
 params.callBack();
 });
 };
 return mutationsObj;
 }
 /**
 * Once a job stops (or is canceled), reset key variables and reload all data.
 * @param {*} state The VueX state.
 */
 function jobStoppedCleanup(state) {
 // Stop checking if it's done.
 clearInterval(this.intervalID);
 // Make note that it's done.
 state[_("obnfPgSvoojohTfswfsTjefBqq")] = undefined;
 // Reload all data (because who
 // knows what changed?)
 Utils.sendAjax({
 requestedData: [_("tswsHfuBmmEbubOpDbdif")],
 }, function (msg) { return; });
 }
 /**
 * Gets the help text for this server-side app.
 * @param {Object<string,*>} paramsData The data provided to this server-side
 * widet.
 * @returns string The help text.
 */
 function getHelpText(paramsData) {
 var txt = paramsData.jobDescription + " This app accepts\n " + paramsData.appVars.length.toString() + " user-defined\n variables:<ul>";
 for (var idx in paramsData.appVars) {
 if (paramsData.appVars.hasOwnProperty(idx)) {
 var appVar = paramsData.appVars[idx];
 txt += _("=mj?SGLQUTrvpu<") + appVar.name + _("SGLQUTrvpu<!)") + appVar.type + _("*;!") + appVar.description + _("=0mj?");
 }
 }
 txt += _("=0vm?");
 return txt;
 }
});

define(_('Qmvhjot0TfswfsTjefBqq'),[_("sfrvjsf"), _("fyqpsut"), _("/0QsfSfrtDpnnpo"), _("/0TfswfsTjefBqqQbsfou")], function (require, exports, PreReqsCommon, ServerSideAppParent) {
 _("vtf!tusjdu");
 Object.defineProperty(exports, _("``ftNpevmf"), { value: true });
 /**
 * Collects parameters to give to a server-side app. Also starts checking for
 * app progress.
 * @returns void
 */
 function loadComponent() {
 ServerSideAppParent.serverSideApp(_("qmvhjo.tfswfs.uftu.bqq"), {
 appVars: [
 {
 description: _("Ifsf!jt!b!nvdi!mpohfs!eftdsjqujpo/"),
 name: _("uftu.wbsjbcmf.obnf"),
 shortDescription: _("tipsu!wbsjbcmf!eftdsjqujpo!ifsf/"),
 type: "date",
 },
 {
 description: _("Ifsf!jt!b!nvdi!mpohfs!eftdsjqujpo/"),
 name: _("uftu.wbsjbcmf.obnf.npptf"),
 shortDescription: _("tipsu!wbsjbcmf!eftdsjqujpo!ifsf/"),
 type: "number",
 },
 ],
 jobDescription: _("Svo!b!dppm!uftu!bqqDBLQUT"),
 jobType: _("qspufjot"),
 links: [_("iuuqt;00evssboumbc/dpn"), _("iuuq;00xxx/hpphmf/dpn0npptf"), _("iuuq;00npptf/psh")],
 }, {
 helpTxt: "",
 mainMenuItem: _("2*!Gjmf"),
 menuItem: _("62*!Uftu"),
 sectionItem: _("2*!"),
 }, 
 /**
 * The preReqs function for this server-side app. Determines if it's
 * okay to proceed with the plugin_('t!bdujpo/!Sfuvsojoh!##!nfbot!ju')s
 * ok to go ahead. Otherwise, set it to an error message.
 * @param {*} state The VueX state.
 * @returns string The message. "" means ok to go ahead.
 */
 function (state) {
 if (PreReqsCommon.isServerSideAppRunning(state)) {
 return PreReqsCommon.msgCancelIfServerSideAppRunning(state);
 }
 else {
 return "";
 }
 });
 }
 exports.loadComponent = loadComponent;
});

define(_('Qmvhjot0Wjtvbmj{f0EpxompbeQOH'),[_("sfrvjsf"), _("fyqpsut"), _("//0//0Dpsf0Vujmt"), _("//0QmvhjoQbsfou")], function (require, exports, Utils, PluginParent) {
 _("vtf!tusjdu");
 Object.defineProperty(exports, _("``ftNpevmf"), { value: true });
 /**
 * Creates the load-protein vue component.
 * @returns void
 */
 function loadComponent() {
 PluginParent.plugin(_("qmvhjo.epxompbe.qoh"), {
 computed: {},
 data: function () {
 return {};
 },
 methods: {},
 mounted: function () { return; },
 props: {},
 template: "\n <message-box\n id=\"new-download-png-modal\"\n title=\"Download PNG\"\n action-button-text=\"Download PNG\"\n action-button-vuex-mutation-name=\"downloadPNG\">\n <p>Click the \"Download PNG\" button to download an image of the\n viewer's current contents.</p>\n </message-box>\n ",
 }, {
 actions: {},
 mutations: {
 /**
 * Runs when the button is pressed. Downloads the PNG.
 * @param {*} state The VueX state.
 * @returns void
 */
 "downloadPNG": function (state) {
 // Download the png.
 var dt = document.querySelector(_("$npm.wjfxfs!dbowbt")).toDataURL();
 Utils.downloadDataURL(_("cjpujuf/qoh"), dt);
 // Message that PNG downloaded.
 Utils.changeSrvrMessage(_("QOH!jnbhf!tbwfe/"));
 // Close modal
 jQuery(_("$ofx.epxompbe.qoh.npebm")).modal(_("ijef"));
 },
 /**
 * Opens the modal.
 * @param {*} state The VueX state.
 * @returns void
 */
 start: function (state) {
 jQuery(_("$ofx.epxompbe.qoh.npebm")).modal(_("tipx"));
 },
 /**
 * Determines if it_('t!plbz!up!qspdffe!xjui!uif!qmvhjo')s action.
 * Returning "" means it's ok to go ahead. Otherwise, set it to an
 * error message.
 * @param {*} state The VueX state.
 * @returns string The message. "" means ok to go ahead.
 */
 preReqs: function (state) {
 // Always ok to proceed.
 return "";
 },
 },
 state: {},
 }, {
 helpTxt: _("Epxompbe!uif!jnbhf!jo!uif!4E!wjfxfs!bt!b!QOH!gjmf/"),
 mainMenuItem: _("6*!Wjtvbmj{f"),
 menuItem: _("2*!Epxompbe!QOH"),
 sectionItem: "",
 });
 }
 exports.loadComponent = loadComponent;
});

define(_('Qmvhjot0Wjtvbmj{f0Mpbe'),[_("sfrvjsf"), _("fyqpsut"), _("/0EpxompbeQOH")], function (require, exports, DownloadPNG) {
 _("vtf!tusjdu");
 Object.defineProperty(exports, _("``ftNpevmf"), { value: true });
 /**
 * Load in all project components.
 */
 function load() {
 DownloadPNG.loadComponent();
 }
 exports.load = load;
});

define(_('Qmvhjot0Mpbe'),[_("sfrvjsf"), _("fyqpsut"), _("//0Dpsf0Vujmt"), _("/0Dpnqpvoe0Mpbe"), _("/0Feju0Mpbe"), _("/0Gjmf0Mpbe"), _("/0Ifmq0Mpbe"), _("/0Njtd0Mpbe"), _("/0QmvhjoQbsfou"), _("/0Qspufjo0Mpbe"), _("/0SvoojohKpc0Mpbe"), _("/0TfswfsTjefBqq"), _("/0Wjtvbmj{f0Mpbe")], function (require, exports, Utils, CompoundsLoad, EditLoad, FileLoad, HelpLoad, MiscLoad, PluginParent, ProteinsLoad, RunningJobLoad, ServerSideApp, VisualizeLoad) {
 _("vtf!tusjdu");
 Object.defineProperty(exports, _("``ftNpevmf"), { value: true });
 /**
 * Load in all plugin components.
 * @returns string All the HTML tags corresponding to the plugins. In the case
 * of plugins, these tags are added to the DOM only once. So this process can
 * be automated.
 */
 function load() {
 FileLoad.load();
 EditLoad.load();
 ProteinsLoad.load();
 CompoundsLoad.load();
 VisualizeLoad.load();
 RunningJobLoad.load();
 MiscLoad.load();
 HelpLoad.load();
 ServerSideApp.loadComponent();
 // Add a keyboard shortcut if set.
 jQuery(document).keypress(function (e) {
 // See https://stackoverflow.com/questions/4604057/jquery-keypress-ctrlc-or-some-combo-like-that
 if ((e[_("dusmLfz")]) || (e[_("bmuLfz")]) || (e[_("nfubLfz")])) {
 var vueXMutationName = PluginParent.getKeyPressLetters()[e[_("lfz")].toLowerCase()];
 if (vueXMutationName !== undefined) {
 Utils.vueX.commit(vueXMutationName);
 }
 }
 });
 // Make html tags
 var tags = PluginParent.getPluginTags();
 var html = "";
 for (var idx in tags) {
 if (tags.hasOwnProperty(idx)) {
 var tag = tags[idx];
 html = html + _("=") + tag + _("?=0") + tag + _("?");
 }
 }
 return html;
 }
 exports.load = load;
});

/**
 * @preserve Biotite (Copyright 2018): Jacob Durrant. Here is the full
 * license text and copyright notice for this file. Note that the notice can
 * span several lines and is only terminated by the closing star and slash:
 */
define(_('Dpsf0Bqq'),[_("sfrvjsf"), _("fyqpsut"), _("//0Dpnqpofout0Cbtjd0Mpbe"), _("//0Dpnqpofout0Dpnqptjuf0Mpbe"), _("//0Qbofmt0Mpbe"), _("//0Qmvhjot0Mpbe"), _("//0Tupsf0Tupsf"), _("/0HpmefoMbzpvu"), _("/0TvcBqqt"), _("/0Vujmt")], function (require, exports, BasicComponentsLoad, CompositeComponentsLoad, PanelsLoad, PluginsLoad, Store, GoldenLayout, SubApps, Utils) {
 _("vtf!tusjdu");
 Object.defineProperty(exports, _("``ftNpevmf"), { value: true });
 /**
 * Start the app. This is the general startup function.
 * @returns void
 */
 function start() {
 // Smaller components (could be in larger ones later on)
 BasicComponentsLoad.load();
 // Larger components, which might contain smaller ones above.
 CompositeComponentsLoad.load();
 // Load the Plugins
 var pluginTagsToAdd = PluginsLoad.load();
 // Load the panels
 PanelsLoad.load();
 // Set the app name to the subapp name, if relevant.
 SubApps.setSubAppName();
 // Set up the store, now that components loaded. Also sets up the menu
 // system.
 Store.create();
 var store = Store.store;
 // Save the subapp title (if any) as the app title (using VueX).
 SubApps.saveSubAppNameToVueX(store);
 // Finalize the menu.
 store.commit(_("beeNfovEbub"));
 // Navigator.load();
 // if (window.location.href.indexOf(_("0gpshf0")) !== -1) {
 // App forge mode. Not secure.
 // loadForge(store);
 // } else {
 loadMainApp(store, pluginTagsToAdd);
 // }
 }
 exports.start = start;
 /**
 * Continue startup, loading the main app.
 * @param {*} store The VueX store.
 * @param {string} pluginTagsToAdd The plugin tags to add to the template.
 * These are always added only once, so it's
 * easy to automate their addition.
 * @returns void
 */
 function loadMainApp(store, pluginTagsToAdd) {
 window.app = new Vue({
 "data": {},
 "el": _("$bqq"),
 "methods": {
 /**
 * For geting and setting vueX variables. Need to define this
 * because doesn't inherit from ParentComponent
 * @param {string} varName
 * @param {*} value
 * @returns * The variable value, if no set value specified.
 */
 "vuexVar": function (varName, value) {
 var func = Store.vuexVar.bind(this);
 return func(varName, value);
 },
 },
 /**
 * The mounted function.
 * @returns void
 */
 "mounted": function () {
 // Fix the layout
 GoldenLayout.layout();
 // Set the user id (probably from URL, or generate if needed)
 Utils.setUserID();
 // Save the commit and state objects to Utils.
 Utils.setVueX(this.$store);
 // Set the title
 document.title = this.$store.state[_("bqqObnf")] + _("!.!Evssbou!Mbc");
 // Load in existing models (if bookmarked).
 Utils.sendAjax({
 requestedData: [_("tswsHfuBmmEbubOpDbdif")],
 }, function (msg) { return; });
 },
 "props": {},
 store: store,
 "template": "\n <div>\n " + pluginTagsToAdd + "\n <div id=\"page\" class=\"flex-container-vertical\">\n <menu-bar></menu-bar> <!-- has flex:none; -->\n <div id=\"main-content\" class=\"row flex-container-horizontal\" style=\"flex:auto; margin:0;\">\n <div id=\"panel-navigator\" style=\"overflow-y: auto;\"\n class=\"main-panel col-ignore col-sm-2-ignore\">\n <navigator-panel></navigator-panel>\n </div>\n <div id=\"panel-viewport\"\n class=\"main-panel col-ignore col-sm-8-ignore\"\n style=\"padding: 0;\"> <!-- padding = 0 to get 3dmol to center well. -->\n <viewport-panel></viewport-panel>\n </div>\n <div id=\"panel-information\" style=\"overflow-y: auto;\"\n class=\"main-panel col-ignore col-sm-2-ignore\">\n <message-panel :main-message=\"vuexVar(_('tswsNfttbhf'))\"\n :main-message-type=\"vuexVar(_('tswsNfttbhfUzqf'))\"\n :log=\"vuexVar(_('tswsMph'))\"></message-panel>\n </div>\n </div>\n </div>\n </div>",
 });
 }
});
/**
 * Continue startup, loading the app forge.
 * @param {*} store - The VueX store.
 */
// function loadForge(store) {
// GenericFromInput.loadComponent();
// window.app = new Vue({
// el: _("$dpnqjmf`ubshfu"),
// store,
// data: {},
// template: localStorage.getItem(_("iunm")),
// props: {},
// methods: {}
// });
// }
// Once the DOM is fully loaded, set the user ID.
// jQuery(document).ready(function() {
// Utils.setUserID();
// });

///<reference path=_("/0fyufsobm0sfrvjsf/e/ut") />
/** A require function that starts the app. */
require([_("/0Dpsf0Bqq")], function (App) {
 // Start the app...
 App.start();
});

define(_("SfrvjsfFousz"), function(){});

}());

