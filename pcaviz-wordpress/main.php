<?php
/**
 * Plugin Name: PCAviz
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
    $txt = pcaviz_add_default_files();
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
    $upload_dir = wp_upload_dir();

    // ob_start();

    // Get the source file in the assets directory.
    $all_flnms = plugin_dir_path(__FILE__).'assets/sims/*.compressed.json';
    foreach(glob($all_flnms) as $src_flnm) {
        // Get the contents of that file.
        $data = file_get_contents($src_flnm);

        // And the base filename.
        $filename = basename($src_flnm);

        // Make the pcaviz directory in uploads if necessary.
        $example_json_outdir = $upload_dir['basedir'].'/pca-sample-trajs';
        wp_mkdir_p($example_json_outdir);

        // Pick the file.
        $file = $example_json_outdir.'/'.$filename;

        // If the file already exists, skip copying it.
        if (!file_exists($file)) {
            // Save the file to that location.
            file_put_contents($file, $data);

            // Get information about the file
            $sim_inf = file($src_flnm.'.inf');

            $wp_filetype = wp_check_filetype($filename, null);

            // Save the attachment to the media library.
            $attachment = array(
                'post_mime_type' => $wp_filetype['type'],
                'post_title' => $sim_inf[0],
                'post_content' => '',
                'post_status' => 'inherit',
                'post_excerpt' => $sim_inf[1]
            );

            $attach_id = wp_insert_attachment($attachment, $file);
        }
    }

    // return ob_get_clean();
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

wp_register_script('PCAViz.min', plugin_dir_url(__FILE__) . 'assets/js/PCAViz.min.js');
wp_enqueue_script('PCAViz.min');
wp_enqueue_script('jquery');
wp_register_script('GoogleAnalytics', plugin_dir_url( __FILE__ ) . 'assets/js/googleanalytics.js');
wp_enqueue_script('GoogleAnalytics');

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

    if (is_null($atts["loop"])) {
        $atts["loop"] = "true";
    }
    $atts["loop"] = strtolower($atts["loop"]);

    if (is_null($atts["autoplay"])) {
        $atts["autoplay"] = "true";
    }
    $atts["autoplay"] = strtolower($atts["autoplay"]);

    if (is_null($atts["visstyle"])) {
        $atts["visstyle"] = "{'cartoon':{}, 'stick':{'radius':0.5,'colorscheme':'Jmol'}}";
    }
    $atts["visstyle"] = pcaviz_protect_quotes($atts["visstyle"]);

    if (is_null($atts["durationinmilliseconds"])) {
        $atts["durationinmilliseconds"] = 10000;
    }

    if (is_null($atts["updatefreqinmilliseconds"])) {
        $atts["updatefreqinmilliseconds"] = 16.67;
    }

    if (is_null($atts["windowaveragesize"])) {
        $atts["windowaveragesize"] = 1;
    }

    if (is_null($atts["caching"])) {
        $atts["caching"] = "none";  // none, continuous, or pre
    }
    $atts["caching"] = strtolower($atts["caching"]);
    $atts["caching"] = trim($atts["caching"]);
    if (!in_array($atts["caching"], Array("none", "continuous", "pre"))) {
        echo "<span style='color:red;'><b>Error! The caching attribute must ".
        "be either \"none\", \"continuous\", or \"pre\".</b></span>";
    }

    // The content is wrapped in a .pcaviz-container div.
    echo "<div style='width:$atts[width]px;' class='pcaviz-container wp-block-image";
    if (!is_null($atts['align'])) {
        echo " align$atts[align]";
    }
    echo "'>";

    // If the user agreed to let us collect statistics...
    if(get_option('pcaviz_option_name') === 'checked'){
        echo "<script>ga('pcaviz.send', 'pageview');</script>";
    }

    // Add the required CSS styles.
    wp_register_style('pcaviz-viscontainer', plugin_dir_url(__FILE__) . 'assets/css/viscontainer.css');
    wp_enqueue_style('pcaviz-viscontainer');

    // If the height attribute was set by user, we overwrite the css setting.
    if(!is_null($atts['height'])) {
        echo "<style>#pcaviz-vis-and-controls {height: $atts[height]px !important;}</style>";
    }

    // If the width attribute was set by user, we overwrite the css setting.
    if(!is_null($atts['width'])) {
        echo "<style>#pcaviz-vis-and-controls {width: $atts[width]px !important;}</style>";
    }

    // The canvas is wrapped into "pcaviz-collapsable" div. There's a script bellow
    // that toggles its visibility.
    echo "    <div id='pcaviz-collapsable'>";

    // Setup the playback pcaviz-controls.
    echo "    <div id='pcaviz-vis-and-controls'>";
    echo "        <div id='pcaviz-viscontainer'></div>";
    echo "        <div id='pcaviz-controls'";

    // Only display playback buttons if the user requested them in the shortcode.
    if ((!is_null($atts['playback_buttons'])) && ($atts['playback_buttons'] === 'false')) {
        echo " style='display: none'";
    }
    echo ">";

    echo "</div>";  // #pcaviz-controls
    echo "</div>";  // #pcaviz-vis-and-controls

    // Add a caption if it's set.
    if(!is_null($atts['caption'])) {
        echo "<figcaption>$atts[caption]</figcaption>";
    }

    echo "</div>";  // #pcaviz-collapsable

    // Load the PCAViz javascript code.
    echo "<script type='text/javascript' src='"
          .plugin_dir_url(__FILE__)."assets/js/viewer.js'"
          ."></script>";

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
            $jsons[] = array(
                "url" => $url,
                "file_name" => $json->post_name,
                "date_time" => $json->post_date,
                "title" => $json->post_title
            );
        }
    }

    // Add javascript for toggling the "pcaviz-collapsable" div.
    echo "<script>
              jQuery('#pcaviz-collapsable').hide();
              function togglevisibility() {
                  jQuery('#pcaviz-collapsable').show();
              }
          </script>";

    // If the user didn't include a 'file' attribute in the shortcode, display
    // a dropdown listing all JSON files from the media library.
    if (is_null($atts['file'])) {
        echo "<select style='width:$atts[width]px' id='pca-file-input' onchange='togglevisibility();makePCAViz(viewer, \"3DMOLJS\", this.options[this.selectedIndex].value, $atts[loop], $atts[autoplay], \"$atts[visstyle]\", $atts[durationinmilliseconds], $atts[updatefreqinmilliseconds], $atts[windowaveragesize], \"$atts[caching]\")'>
                  <option selected>Select file</option>";
        foreach ($jsons as $json){
            echo "<option value=$json[url]>$json[title] | $json[date_time]</option>";
        }
        echo '</select>';
    } else {
        // The user provided a 'file' attribute in the shortcode. Search if
        // the user-provided file name is a substring of any file in the media
        // library. Output javascript code to display that file once it is
        // found.
        $file_found = FALSE;
        foreach ($jsons as $json){
            // if (strpos($json["file_name"], preg_replace('/\\.[^.\\s]{3,4}$/', '', $atts['file'])) !== false){
            $file_to_test = strtolower($atts['file']);
            if (
                    (strpos(strtolower($json["file_name"]), preg_replace('/\\.[^.\\s]{3,4}$/', '', $file_to_test)) !== false) |
                    (strpos(strtolower($json["url"]), preg_replace('/\\.[^.\\s]{3,4}$/', '', $file_to_test)) !== false) |
                    (strpos(strtolower($json["title"]), preg_replace('/\\.[^.\\s]{3,4}$/', '', $file_to_test)) !== false)
                )
                {
                echo "<script>
                          togglevisibility();
                          makePCAViz(viewer, \"3DMOLJS\", \"$json[url]\", $atts[loop], $atts[autoplay], \"$atts[visstyle]\", $atts[durationinmilliseconds], $atts[updatefreqinmilliseconds], $atts[windowaveragesize], \"$atts[caching]\");
                      </script>";
                $file_found = TRUE;
                break;
            }
        }
        if (!$file_found) {
            // File not found.
            echo "<span style='color:red;'><b>Error! No trajectory file "
                 ."(*.compressed.json) in the media library has a file name, "
                 ."url, or title that contains the string \""
                 .$file_to_test."\".</b></span>";
        }
    }

    echo "</div>";  // .pcaviz-container

    // Return all the above code as a string.
    return ob_get_clean();
}

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
        <h2>PCAviz Plugin Settings</h2>
        <form method="post" action="options.php">
            <?php settings_fields( 'pcaviz_options_group' ); ?>
            <h3>Please allow us to collect limited data on PCAViz usage.</h3>
            <p>
                As an academic group, our research depends on government grants.
                Being able to cite usage statistics in reports and articles will
                help us convince reviewers that our work is worth funding. <b>We
                understand and respect your privacy and that of your
                visitors.</b> Your personal data will never be made public. We
                will only use summarized, aggregate data in our reports.
            </p>
            <p>
                Thank you! This helps us a lot!
            </p>
            <p>
                <input type="checkbox" id="pcaviz_option_name" name="pcaviz_option_name" value="checked" <?php echo get_option('pcaviz_option_name'); ?>/>
                Allow Analytics
            </p>
            <?php  submit_button(); ?>
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
  add_options_page('Page Title', 'PCAviz', 'manage_options', 'pcaviz', 'pcaviz_options_page');
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
