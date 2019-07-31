<?php
/**
 * Plugin Name: PCAViz
 *
 * @package     PCAViz
 * @author      Yuri Kochnev, Jacob Durrant
 * @copyright   2019 Durrant Lab
 * @license     Apache-2.0
 *
 * @wordpress-plugin
 * Plugin Name: PCAViz
 * Plugin URI:  http://durrantlab.com/pcaviz/
 * Description: Display molecular dynamics simulations in the browser!
 * Version:     1.0.0
 * Author:      Yuri Kochnev, Jacob Durrant
 * Author URI:  http://durrantlab.com
 * Text Domain: pcaviz
 * License:     Apache-2.0
 * License URI: https://www.apache.org/licenses/LICENSE-2.0
 */

/**
 * Ask the user to go to the plugin's settings page and opt-in to the
 * google-analytics option.
 *
 * @return void
 */
function pcaviz_plugin_activation() {
    $notices = get_option('pcaviz_plugin_deferred_admin_notices', array());
    $notices[] = "<b>Please visit the PCAViz <a href='/wp-admin/options-general.php?page=pcaviz'>".
    "settings page</a> to enable Google-analytics tracking.</b> Allowing us "
    ."to collect and report minimal data on PCAViz usage will improve our "
    ."chances of getting grants to fund our ongoing work. Not to worry. We "
    ."respect user privacy!";
    update_option('pcaviz_plugin_deferred_admin_notices', $notices);

    // Copy default trajectories to the media library
    pcaviz_add_default_files();
}
register_activation_hook( __FILE__, 'pcaviz_plugin_activation' );

/**
 * Adds sample JSON trajectory files to the media library.
 *
 * @return void
 */
function pcaviz_add_default_files() {
    // See https://wordpress.stackexchange.com/questions/256830/programmatically-adding-images-to-media-library

    // Get info re. the upload directory.
    // $upload_dir = wp_upload_dir();

    // Need to require these files. See
    // https://codex.wordpress.org/Function_Reference/media_handle_sideload
	if ( !function_exists('media_handle_upload') ) {
		require_once(ABSPATH."wp-admin".'/includes/image.php');
		require_once(ABSPATH."wp-admin".'/includes/file.php');
		require_once(ABSPATH."wp-admin".'/includes/media.php');
    }

    // Get plugin sim url.
    $sim_url = plugins_url('assets/sims/', __FILE__);

    // Get a list of all these previously uploaded json files, so you don't
    // reuploaded.
    $uploaded_jsons = pcaviz_get_all_json_files();
    $uploaded_json_filenames = Array();
    foreach ($uploaded_jsons as $uploaded_json) {
        $uploaded_json_filenames[] = basename($uploaded_json["url"]);
    }

    // ob_start();  // For debugging

    // Get the source file in the assets directory.
    $all_flnms = plugin_dir_path(__FILE__).'assets/sims/*.compressed.json';
    foreach(glob($all_flnms) as $src_flnm) {
        // Basename
        $bsnm = basename($src_flnm);

        if (in_array($bsnm, $uploaded_json_filenames)) {
            // Already uploaded;
            continue;
        }

        // Get the URL of the file.
        $url = $sim_url.$bsnm;

        // Download it to a temporary file.
        $tmp = download_url($url);

        $file_array = array(
            'name' => basename($url),
            'tmp_name' => $tmp,
            'type' => 'text/plain'
        );

        // If there are download errors, unlink the file.
        if (is_wp_error($tmp)) {
            @unlink($file_array['tmp_name']);
            return $tmp;
        }

        // Now get title and description.
        $sim_inf = file($src_flnm.'.inf');
        $title = $sim_inf[0];
        $caption = $sim_inf[1];

        // Collect all the data about the simulation.
        $attachment = array(
            'post_mime_type' => 'text/plain',
            'post_title' => $title,
            'post_content' => '',
            'post_status' => 'inherit',
            'post_excerpt' => $caption
        );

        // "Upload" the file, without atttaching it to any post.
        $post_id = '0';
        $id = media_handle_sideload($file_array, $post_id, $title, $attachment);

        // Delete the temporary file if there were errors moving it into the
        // media library.
        if (is_wp_error($id)) {
            @unlink($file_array['tmp_name']);
            return $id;
        }

        // You can get the url of the attachment like this, in case you ever
        // need it...
        // $value = wp_get_attachment_url($id);
    }

    // return ob_get_clean();  // For debugging
}

/**
 * Displays any notices (e.g., the output of display_notice()).
 *
 * @return void
 */
function pcaviz_plugin_admin_notices() {
    if ($notices = get_option('pcaviz_plugin_deferred_admin_notices')) {
        foreach ($notices as $notice) {
            echo "<div class='updated'><p>$notice</p></div>";
        }
        delete_option('pcaviz_plugin_deferred_admin_notices');
    }
}
add_action('admin_notices', 'pcaviz_plugin_admin_notices');

/**
 * Notify wordpress filters that we want to allow JSON files.
 *
 * @param array $existing_mimes  An array containing mime times that are
 *                               already allowed.
 * @return array  The same array, with json information added.
 */
function pcaviz_myme_types($existing_mimes=array()) {
    $existing_mimes['json'] = 'text/plain';
    return $existing_mimes;
}
add_filter( 'upload_mimes', 'pcaviz_myme_types' );

// Register the necessary scripts.
wp_register_script('3Dmol-min', plugin_dir_url(__FILE__) . 'assets/js/3Dmol-nojquery-min.js', array('jquery'));
wp_enqueue_script('3Dmol-min');

wp_register_script('PCAViz-min', plugin_dir_url(__FILE__) . 'assets/js/PCAViz.min.js');  // Should not require jquery
wp_enqueue_script('PCAViz-min');

wp_register_script('PCAViz-support', plugin_dir_url(__FILE__) . 'assets/js/viewer.js', array('PCAViz-min'));
wp_enqueue_script('PCAViz-support');

wp_enqueue_script('jquery');

/**
 * The main function for displaying PCAViz in the browser.
 *
 * @param array  $atts      The attributes.
 * @param string $content   Not actually used.
 * @param string $tag       Not actually used.
 * @return string The HTML required to render PCAViz.
 */
function pcaviz_main($atts = [], $content = null, $tag = '') {
    // a shortcode should never produce any output but return text! see
    // https://stackoverflow.com/a/40500555/5907621 for details
    ob_start();

    // Normalize attribute keys, lowercase.
    $atts = array_change_key_case((array)$atts, CASE_LOWER);

    // Get a random ID so you can have multiple pcaviz viewers on a page.
    $uniqid = uniqid();

    // Set some default values.

    // We need a default width for 3dmoljs.
    if (is_null($atts['width'])) {
        $atts['width'] = 333;
    }

    // Good to set default height too.
    if (is_null($atts['height'])) {
        // Golden ratio
        $atts['height'] = $atts['width'] / 1.61803398875;
    }

    // Loop
    if (is_null($atts["loop"])) {
        $atts["loop"] = "true";
    }

    // Autoplay
    if (is_null($atts["autoplay"])) {
        $atts["autoplay"] = "true";
    }

    // VisStyle
    if (is_null($atts["visstyle"])) {
        $atts["visstyle"] = "{'cartoon':{}, 'stick':{'radius':0.5,'colorscheme':'Jmol'}}";
    }
    $atts["visstyle"] = pcaviz_protect_quotes($atts["visstyle"]);

    // durationinmilliseconds
    if (is_null($atts["durationinmilliseconds"])) {
        $atts["durationinmilliseconds"] = 10000;
    }

    // updatefreqinmilliseconds
    if (is_null($atts["updatefreqinmilliseconds"])) {
        $atts["updatefreqinmilliseconds"] = 16.67;
    }

    // windowaveragesize
    if (is_null($atts["windowaveragesize"])) {
        $atts["windowaveragesize"] = 1;
    }

    // caching
    if (is_null($atts["caching"])) {
        $atts["caching"] = "none";  // none, continuous, or pre
    }
    $atts["caching"] = trim($atts["caching"]);
    if (!in_array($atts["caching"], Array("none", "continuous", "pre"))) {
        echo "<span style='color:red;'><b>Error! The caching attribute must ".
        "be either \"none\", \"continuous\", or \"pre\".</b></span>";
    }

    // medialibrarycaption
    if (is_null($atts["medialibrarycaption"])) {
        $atts["medialibrarycaption"] = "false";
    }
    if (is_null($atts["file"])) {
        // No file specified? They use media library captions. Because it will
        // depend on the one selected.
        $atts["medialibrarycaption"] = "true";
    }

    // The content is wrapped in a .pcaviz-container div.
    echo "<div style='width:$atts[width]px;' id='pcaviz-container-".$uniqid."' class='pcaviz-container wp-block-image";
    if (!is_null($atts['align'])) {
        echo " align$atts[align]";
    }
    echo "'>";

    // If the user agreed to let us collect statistics...
    if(get_option('pcaviz_option_name') === 'checked'){
        $use_analytics = "true";
    } else {
        $use_analytics = "false";
    }

    // Add the required CSS styles.
    wp_register_style('pcaviz-viscontainer', plugin_dir_url(__FILE__) . 'assets/css/viscontainer.css');
    wp_enqueue_style('pcaviz-viscontainer');

    // If the height attribute was set by user, we overwrite the css setting.
    if(!is_null($atts['height'])) {
        echo "<style>#pcaviz-vis-and-controls-".$uniqid." {height: $atts[height]px !important;}</style>";
    }

    // If the width attribute was set by user, we overwrite the css setting.
    if(!is_null($atts['width'])) {
        echo "<style>#pcaviz-vis-and-controls-".$uniqid." {width: $atts[width]px !important;}</style>";
    }

    // The canvas is wrapped into "pcaviz-collapsable" div. There's a script bellow
    // that toggles its visibility.
    echo "    <div id='pcaviz-collapsable-".$uniqid."' style='display: none;'>";

    // Setup the playback pcaviz-controls.
    echo "    <div class='pcaviz-vis-and-controls' id='pcaviz-vis-and-controls-".$uniqid."'>";
    echo "        <div class='pcaviz-viscontainer' id='pcaviz-viscontainer-".$uniqid."'></div>";
    echo "        <div class='pcaviz-controls' id='pcaviz-controls-".$uniqid."'";

    // Only display playback buttons if the user requested them in the shortcode.
    if ((!is_null($atts['playback_buttons'])) && ($atts['playback_buttons'] === 'false')) {
        echo " style='display: none'";
    }
    echo ">";

    echo "</div>";  // #pcaviz-controls
    echo "</div>";  // #pcaviz-vis-and-controls

    $jsons = pcaviz_get_all_json_files();

    // If the user didn't include a 'file' attribute in the shortcode, display
    // a dropdown listing all JSON files from the media library.
    $html = "";
    $medialibrarycaption = "";
    if (is_null($atts['file'])) {
        $html .= "<select style='width:$atts[width]px' id='pca-file-input-".$uniqid."' onchange='pcaVizToggleVisibility(\"".$uniqid."\");pcaVizSetCaption(\"".$uniqid."\");makePCAViz(\"".$uniqid."\", \"3DMOLJS\", this.options[this.selectedIndex].value, $atts[loop], $atts[autoplay], \"$atts[visstyle]\", $atts[durationinmilliseconds], $atts[updatefreqinmilliseconds], $atts[windowaveragesize], \"$atts[caching]\", $use_analytics)'>
                  <option selected>Select file</option>";
        foreach ($jsons as $json){
            $html .= "<option data-caption='$json[caption]' value=$json[url]>$json[title] | $json[date_time]</option>";
        }
        $html .= '</select>';
    } else {
        // The user provided a 'file' attribute in the shortcode. Search if
        // the user-provided file name is a substring of any file in the media
        // library. Output javascript code to display that file once it is
        // found.
        $file_found = FALSE;
        foreach ($jsons as $json){
            // if (strpos($json["file_name"], preg_replace('/\\.[^.\\s]{3,4}$/', '', $atts['file'])) !== false){
            $file_to_test = $atts['file'];
            if (
                    (strpos(strtolower($json["file_name"]), preg_replace('/\\.[^.\\s]{3,4}$/', '', $file_to_test)) !== false) |
                    (strpos(strtolower($json["url"]), preg_replace('/\\.[^.\\s]{3,4}$/', '', $file_to_test)) !== false) |
                    (strpos(strtolower($json["title"]), preg_replace('/\\.[^.\\s]{3,4}$/', '', $file_to_test)) !== false)
                )
                {
                // Keep trying to show the visibility. Doing it this way
                // because theme might change where javascript files are
                // loaded (footer).
                $html .= "<script>
                          let pcaVizFileSpecSetup = () => {
                              if (pcaVizToggleVisibility !== undefined) {
                                  pcaVizToggleVisibility(\"".$uniqid."\");
                                  makePCAViz(\"".$uniqid."\", \"3DMOLJS\", \"$json[url]\", $atts[loop], $atts[autoplay], \"$atts[visstyle]\", $atts[durationinmilliseconds], $atts[updatefreqinmilliseconds], $atts[windowaveragesize], \"$atts[caching]\", $use_analytics);
                              } else {
                                  setTimeout(pcaVizFileSpecSetup, 100);
                              }
                          };
                          setTimeout(pcaVizFileSpecSetup, 100);
                      </script>";
                $file_found = TRUE;
                $medialibrarycaption = $json["caption"];
                break;
            }
        }
        if (!$file_found) {
            // File not found.
            $html .= "<span style='color:red;'><b>Error! No trajectory file "
                  ."(*.compressed.json) in the media library has a file name, "
                  ."url, or title that contains the string \""
                  .$file_to_test."\".</b></span>";
        }
    }

    // If it's medialibrarycaption, use that caption.
    if ($atts["medialibrarycaption"] == "true") {
        $atts['caption'] = $medialibrarycaption;
    }

    // Add a caption if it's set.
    if(!is_null($atts['caption'])) {
        // Caption must not have quote
        $atts['caption'] = str_replace('"', '', $atts['caption']);

        echo "<figcaption class='pcaviz-figcaption' id='pcaviz-figcaption-".$uniqid."'>$atts[caption]</figcaption>";
    }

    echo "</div>";  // #pcaviz-collapsable

    // Load the PCAViz javascript code.
    echo $html;

    echo "</div>";  // .pcaviz-container

    // Return all the above code as a string.
    return ob_get_clean();
}

/**
 * Gets information about the uploaded .compressed.json files.
 *
 * @return array The information about the .compressed.json files.
 */
function pcaviz_get_all_json_files() {
    // Get a list of all the JSON files that have been uploaded to the media
    // library. Put that list in a variable called $jsons.
    $query_json_args = array(
        'post_type'      => 'attachment',
        'post_mime_type' => 'text/plain',
        'post_status'    => 'inherit',
        'posts_per_page' => - 1,
    );
    $query_jsons = new WP_Query( $query_json_args  );
    $jsons = array();
    foreach ($query_jsons->posts as $json) {
        $url = wp_get_attachment_url($json->ID);
        if (substr($url, -strlen(".compressed.json")) === ".compressed.json") {
            // Filename ends with ".compressed.json"

            // Caption must not have quote
            $caption = str_replace('"', '', $json->post_excerpt);

            $jsons[] = array(
                "url" => $url,
                "file_name" => $json->post_name,
                "date_time" => $json->post_date,
                "title" => $json->post_title,
                "caption" => $caption
            );
        }
    }

    return $jsons;
}


/**
 * Replaces ' with ". Then replaces " with !QUOTE!. viewer.js will replace
 * !QUOTE! with " to restore the quote.
 *
 * @param $str $str
 * @return string The string with the quotes protected.
 */
function pcaviz_protect_quotes($str) {
    $str = str_replace("'", '"', $str);
    $str = str_replace('"', '!QUOTE!', $str);
    return $str;
}

/**
 * This method draws the plugin's settings page.
 *
 * @return void
 */
function pcaviz_options_page() { ?>
    <div>
        <?php screen_icon(); ?>
        <h2>PCAViz Plugin Settings</h2>
        <form method="post" action="options.php">
            <?php settings_fields( 'pcaviz_options_group' ); ?>
            <h3>Please allow us to collect limited data on PCAViz usage.</h3>
            <p>
                As an academic group, the <a href="http://durrantlab.com"
                target="_blank">Durrant Lab</a>'s research depends on
                government grants. Being able to cite usage statistics in
                reports and articles will help us convince reviewers that our
                work is worth funding. <b>We understand and respect your
                privacy and that of your visitors.</b> Your personal data will
                never be made public. We will only use summarized, aggregate
                data in our reports.
            </p>
            <p>
                Thank you! This helps us a lot!
            </p>
            <p>
                <input type="checkbox" id="pcaviz_option_name" name="pcaviz_option_name" value="checked" <?php echo get_option('pcaviz_option_name'); ?>/>
                Allow Analytics
            </p>
            <p>
                (If you are using a caching plugin, you may need to
                clear your cache for any changes to take effect.)
            </p>
            <?php  submit_button(); ?>
            <hr />
            <h3 id="examplesofuse">Examples of Use</h3>

            <p>The PCAViz WordPress Plugin uses shortcodes that can be inserted into any WordPress post or page.</p>

            <h4 id="majorattributes">Major Attributes</h4>

            <ol>
                <li>If you use the plugin without any attributes, it will allow your site visitors to select from any of the compressed JSON files you've uploaded to your WordPress media library. Note that these files must end in
                    <code>.compressed.json</code> to be recognized.
                    <br>
                    <code>[pcaviz]</code></li>

                <li>If you want to show your site visitors a specific simulation present in your media library, simplify specify the name of the simulation. PCAViz will search through your media library and will choose the first JSON file ending in <code>.compressed.json</code> that includes your specified text in its title, URL, or file name.
                    <br>
                    <code>[pcaviz file="larp1"]</code></li>

                <li>You can also control the molecular styles. The WordPress plugin uses
                    <a href="https://3dmol.csb.pitt.edu">3DMol.js</a> to render molecules. Passing the plugin a <a href="https://3dmol.csb.pitt.edu/doc/types.html#AtomStyleSpec">3DMol.js
            AtomStyleSpec</a> JSON string changes the way the molecules are displayed.
                    <br>
                    <code>[pcaviz visStyle='{"cartoon": {"style": "trace", "color": "grey", "opacity": 0.75}}']</code>
                    <br> Note that the <a href="https://3dmol.csb.pitt.edu/doc/$3Dmol.GLViewer.html#setStyle">3DMol.js
            documentation</a> shows AtomStyleSpec examples as JavaScript objects. PCAViz accepts only JSON strings. The two look similar, but with important differences (e.g., keys in JSON strings must always be quoted). We recommend using a <a href="https://jsonformatter.curiousconcept.com">JSON
            validator</a> to check your AtomStyleSpec strings.</li>
            </ol>

            <h4 id="minorattributes">Minor Attributes</h4>

            <ol>
                <li>To explicitly specify the height and width of the viewer, in pixels:
                    <br>
                    <code>[pcaviz height=150 width=150]</code></li>

                <li>To align the viewer to the left, center, or right, and to optionally add a caption:
                    <br>
                    <code>[pcaviz align="right" caption="My molecule in motion!"]</code></li>

                <li>To instead use the caption associated with the <code>.compressed.json</code> file in your media library:
                    <br>
                    <code>[pcaviz mediaLibraryCaption="true"]</code>
                    <br>(Note that this option is automatically set to true when no <code>file</code> attribute is specified.)</li>

                <li>To hide the playback buttons:
                    <br>
                    <code>[pcaviz playback_buttons="false"]</code></li>

                <li>To prevent the plugin from looping the animation:
                    <br>
                    <code>[pcaviz loop="false"]</code></li>

                <li>The plugin starts playing the simulation automatically by default. To deactivate autoplay:
                    <br>
                    <code>[pcaviz autoplay="false"]</code></li>

                <li>To specify the duration of the animation and how frequently the atomic positions are updated:
                    <br>
                    <code>[pcaviz durationInMilliseconds=10000 updateFreqInMilliseconds=50]</code></li>

                <li>To smooth the animation by averaging the atomic coordinates over multiple frames:
                    <br>
                    <code>[pcaviz windowAverageSize=25]</code></li>

                <li>To control how the plugin caches calculated coordinates:
                    <br>
                    <code>[pcaviz caching="continuous"]</code>
                    <br> Acceptable values are "none" (no caching), "continuous" (cache each frame's coordinates after they are first calculated), and "pre" (calculate and cache all frame coordinates before starting the animation).</li>
            </ol>
        </form>
    </div>
<?php }

/**
 * Register the settings.
 *
 * @return void
 */
function pcaviz_register_settings() {
   add_option( 'pcaviz_option_name', 'true');
   register_setting( 'pcaviz_options_group', 'pcaviz_option_name', 'pcaviz_callback' );
}
add_action( 'admin_init', 'pcaviz_register_settings' );

/**
 * Register the options page.
 *
 * @return void
 */
function pcaviz_register_options_page() {
  add_options_page('PCAViz Settings', 'PCAViz', 'manage_options', 'pcaviz', 'pcaviz_options_page');
}
add_action('admin_menu', 'pcaviz_register_options_page');

/**
 * Let wordpress know about our shortcode.
 *
 * @return void
 */
function pcaviz_shortcodes_init() {
    add_shortcode( 'pcaviz', 'pcaviz_main' );
}
add_action('init', 'pcaviz_shortcodes_init');

?>
