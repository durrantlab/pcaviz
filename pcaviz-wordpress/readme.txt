=== PCAViz ===
Contributors: durrantlab
Donate link: http://durrantlab.com/pcaviz/
Tags: PCAViz, molecular dynamics simulations
Requires at least: 5.2.2
Tested up to: 5.2.2
Stable tag: trunk
Requires PHP: 7.3.3
License: GPLv2
License URI: https://www.gnu.org/licenses/gpl-2.0.html

Display molecular dynamics simulations in the browser!

== Description ==

To encourage use, an easy-to-install PCAViz-powered WordPress plugin enables "plug-and-play" trajectory visualization. The PCAViz download also provides examples showing how to use PCAViz in any webpage.

== Installation ==

You can install the PCAViz WordPress plugin via the official WordPress Plugin Directory. Search for __XXXXX__. The [wordpress.org](https://wordpress.org/support/article/managing-plugins/#automatic-plugin-installation) site provides useful instructions detailing how to install plugins from their directory.

You can also install the plugin using this repository's `pcaviz-wordpress/` directory. Either copy `pcaviz-wordpress/` directly to your server's `wp-content/plugins/` directory ([instructions](https://wordpress.org/support/article/managing-plugins/#manual-plugin-installation)), or create a ZIP file of the `pcaviz-wordpress/` directory and upload it via your WordPress admin screen ([instructions](https://wordpress.org/support/article/managing-plugins/#manual-upload-via-wordpress-admin)).

== Frequently Asked Questions ==

= How can I compress my MD trajectories to use with this plugin? =

Download the [PCAViz Compressor](http://git.durrantlab.com/jdurrant/pcaviz) to compress your trajectories.

== Screenshots ==

1. An example of the PCAViz WordPress plugin.

== Changelog ==

= 1.0 =
* The original version.

== Upgrade Notice ==

= 1.0 =
The original version.

== Examples of Use ==

The PCAViz WordPress Plugin uses shortcodes that can be inserted into any WordPress post or page.

### Major Attributes ###

1. If you use the plugin without any attributes, it will allow your site visitors to select from any of the compressed JSON files you've uploaded to your WordPress media library. Note that these files must end in `.compressed.json` to be recognized.
   `[pcaviz]`
2. If you want to show your site visitors a specific simulation present in your media library, simplify specify the name of the simulation. PCAViz will search through your media library and will choose the first JSON file ending in `.compressed.json` that includes your specified text in its title, URL, or file name.
   `[pcaviz file="larp1"]`
3. You can also control the molecular styles. The WordPress plugin uses [3DMol.js](https://3dmol.csb.pitt.edu) to render molecules. Passing the plugin a [3DMol.js AtomStyleSpec](https://3dmol.csb.pitt.edu/doc/types.html#AtomStyleSpec) JSON string changes the way the molecules are displayed.
   `[pcaviz visStyle='{"cartoon": {"style": "trace", "color": "grey", "opacity": 0.75}}']`
   Note that the [3DMol.js documentation](https://3dmol.csb.pitt.edu/doc/$3Dmol.GLViewer.html#setStyle) shows AtomStyleSpec examples as JavaScript objects. PCAViz accepts only JSON strings. The two look similar, but with important differences (e.g., keys in JSON strings must always be quoted). We recommend using a [JSON validator](https://jsonformatter.curiousconcept.com) to check your AtomStyleSpec strings.

### Minor Attributes ###

1. To explicitly specify the height and width of the viewer, in pixels:
   `[pcaviz height=150 width=150]`
2. To align the viewer to the left, center, or right, and to optionally add a caption:
   `[pcaviz align="right" caption="My molecule in motion!"]`
3. To instead use the caption associated with the `.compressed.json` file in your media library:
   `[pcaviz mediaLibraryCaption="true"]`
   (Note that this option is automatically set to true when no `file` attribute is specified.)
4. To hide the playback buttons:
   `[pcaviz playback_buttons="false"]`
5. To prevent the plugin from looping the animation:
   `[pcaviz loop="false"]`
6. The plugin starts playing the simulation automatically by default. To deactivate autoplay:
   `[pcaviz autoplay="false"]`
7. To specify the duration of the animation and how frequently the atomic positions are updated:
   `[pcaviz durationInMilliseconds=10000 updateFreqInMilliseconds=50]`
8. To smooth the animation by averaging the atomic coordinates over multiple frames:
   `[pcaviz windowAverageSize=25]`
9. To control how the plugin caches calculated coordinates:
   `[pcaviz caching="continuous"]`
   Acceptable values are "none" (no caching), "continuous" (cache each frame's coordinates after they are first calculated), and "pre" (calculate and cache all frame coordinates before starting the animation).
