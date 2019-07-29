=== PCAViz ===
Contributors: durrantlab
Donate link: http://durrantlab.com
Tags: PCAViz, molecular dynamics simulations
Requires at least: 5.2.2
Tested up to: 5.2.2
Stable tag: trunk
Requires PHP: 7.3.3
License: GPLv2
License URI: https://www.gnu.org/licenses/gpl-2.0.html

Display molecular dynamics simulations in the browser!

== Description ==

To encourage use, an easy-to-install PCAViz-powered WordPress plugin enables
"plug-and-play" trajectory visualization. The PCAViz download also provides
examples showing how to use PCAViz in any webpage.

== Installation ==

You can install the PCAViz WordPress plugin via the official WordPress Plugin
Directory. Search for __XXXXX__. The
[wordpress.org](https://wordpress.org/support/article/managing-plugins/#automatic-plugin-installation)
site provides useful instructions detailing how to install plugins from their
directory.

You can also install the plugin using this repository's `pcaviz-wordpress/`
directory. Either copy `pcaviz-wordpress/` directly to your server's
`wp-content/plugins/` directory
([instructions](https://wordpress.org/support/article/managing-plugins/#manual-plugin-installation)),
or create a ZIP file of the `pcaviz-wordpress/` directory and upload it via
your WordPress admin screen
([instructions](https://wordpress.org/support/article/managing-plugins/#manual-upload-via-wordpress-admin)).

== Frequently Asked Questions ==

= How can I compress my MD trajectories to use with this plugin? =

Download the [PCAViz
Compressor](http://git.durrantlab.com/jdurrant/pcaviz) to compress your
trajectories.

== Screenshots ==

1. This screen shot description corresponds to screenshot-1.(png|jpg|jpeg|gif). Note that the screenshot is taken from
the /assets directory or the directory that contains the stable readme.txt (tags or trunk). Screenshots in the /assets
directory take precedence. For example, `/assets/screenshot-1.png` would win over `/tags/4.3/screenshot-1.png`
(or jpg, jpeg, gif).
2. This is the second screen shot

== Changelog ==

= 1.0 =
* The original version.

== Upgrade Notice ==

= 1.0 =
The original version.

== Examples of Use ==

The PCAViz WordPress Plugin uses shortcodes that can be inserted into any
WordPress post or page.

### Major Attributes ###

1. If you use the plugin without any attributes, it will allow your site
   visitors to select from any of the compressed JSON files you've uploaded to
   your WordPress media library. Note that these files must end in
   `.compressed.json` to be recognized. <br>
   `[pcaviz]`
2. If you want to show your site visitors a specific simulation present in
   your media library, simplify specify the name of the simulation. PCAViz
   will search through your media library and will choose the first JSON file
   ending in `.compressed.json` that includes your specified text in its
   title, URL, or file name. <br>
   `[pcaviz file="larp1"]`
3. You can also control the molecular styles. The WordPress plugin uses
   [3DMol.js](https://3dmol.csb.pitt.edu) to render molecules. Passing the
   plugin a [3DMol.js
   AtomStyleSpec](https://3dmol.csb.pitt.edu/doc/types.html#AtomStyleSpec)
   JSON string changes the way the molecules are displayed. <br>
   `[pcaviz visStyle='{"cartoon": {"style": "trace", "color": "grey", "opacity": 0.75}}']` <br>
   Note that the [3DMol.js
   documentation](https://3dmol.csb.pitt.edu/doc/$3Dmol.GLViewer.html#setStyle)
   shows AtomStyleSpec examples as JavaScript objects. PCAViz accepts only
   JSON strings. The two look similar, but with important differences (e.g.,
   keys in JSON strings must always be quoted). We recommend using a [JSON
   validator](https://jsonformatter.curiousconcept.com) to check your
   AtomStyleSpec strings.

### Minor Attributes ###

1. To explicitly specify the height and width of the viewer, in pixels:
   <br>
   `[pcaviz height=150 width=150]`
2. To align the viewer to the left, center, or right, and to optionally add a
   caption: <br>
   `[pcaviz align="right" caption="My molecule in motion!"]`
3. To instead use the caption associated with the `.compressed.json` file in
   your media library: <br>
   `[pcaviz mediaLibraryCaption="true"]` <br>
   (Note that this option is automatically set to true when no `file`
   attribute is specified.)
4. To hide the playback buttons: <br>
   `[pcaviz playback_buttons="false"]`
5. To prevent the plugin from looping the animation: <br>
   `[pcaviz loop="false"]`
6. The plugin starts playing the simulation automatically by default. To
   deactivate autoplay: <br>
   `[pcaviz autoplay="false"]`
7. To specify the duration of the animation and how frequently the atomic
   positions are updated: <br>
   `[pcaviz durationInMilliseconds=10000 updateFreqInMilliseconds=50]`
8. To smooth the animation by averaging the atomic coordinates over multiple
   frames: <br>
   `[pcaviz windowAverageSize=25]`
9. To control how the plugin caches calculated coordinates: <br>
   `[pcaviz caching="continuous"]` <br>
   Acceptable values are "none" (no caching), "continuous" (cache each frame's
   coordinates after they are first calculated), and "pre" (calculate and
   cache all frame coordinates before starting the animation).
