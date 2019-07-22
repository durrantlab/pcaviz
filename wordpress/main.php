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
    $notices[] = "<b>Please visit the PCAViz <a href='/wp-admin/options-general.php?page=durrant'>".
    "settings page</a> to enable Google-analytics tracking.</b> Allowing us "
    ."to collect and report minimal data on PCAViz usage will improve our "
    ."chances of getting grants to fund our ongoing work. Not to worry. We "
    ."respect user privacy!";
    update_option('pcaviz_plugin_deferred_admin_notices', $notices);
}
register_activation_hook( __FILE__, 'pcaviz_plugin_activation' );

/**
 * Displays any notices (e.g., the output of display_notice()).
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
function my_myme_types($existing_mimes=array()) {
    $existing_mimes['json'] = 'text/plain';
    return $existing_mimes;
}
add_filter( 'upload_mimes', 'my_myme_types' );

// Register the necessary scripts.
// wp_register_script('3Dmol-min', plugin_dir_url(__FILE__) . 'assets/js/3Dmol-min.js');
wp_register_script('3Dmol-min', plugin_dir_url(__FILE__) . 'assets/js/3Dmol-nojquery-min.js', array('jquery'));
wp_enqueue_script('3Dmol-min');

wp_register_script('BrowserSim.min', plugin_dir_url(__FILE__) . 'assets/js/BrowserSim.min.js');
wp_enqueue_script('BrowserSim.min');
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

    // If the user agreed to let us collect statistics...
    if(get_option('durrant_option_name') === 'checked'){
        echo "<script>ga('pcaviz.send', 'pageview');</script>";
    }

    // Add the required CSS styles.
    wp_register_style('viscontainer', plugin_dir_url(__FILE__) . 'assets/css/viscontainer.css');
    wp_enqueue_style('viscontainer');

    // If the height attribute was set by user, we overwrite the css setting.
    if(!is_null($atts['height'])) {
        echo "<style>#vis-and-controls {height: $atts[height] !important;}</style>";
    }

    // If the width attribute was set by user ,we overwrite the css setting.
    if(!is_null($atts['width'])) {
        echo "<style>#vis-and-controls {width: $atts[width] !important;}</style>";
    }

    // The canvas is wrapped into "collapsable" div. There's a script bellow
    // that toggles its visibility.
    echo "<div id='collapsable' class='wp-block-image'";
    if (!is_null($atts['align'])) {
        echo " align=$atts[align]";
    }
    echo ">";

    // Setup the playback controls.
    echo "<div id='vis-and-controls'>
              <div id='viscontainer'></div>
              <div id='controls'";

    // Only display playback buttons if the user requested them in the shortcode.
    if ((!is_null($atts['playback_buttons'])) && ($atts['playback_buttons'] === 'false')) {
        echo " style='display: none'";
    }
    echo ">";

    echo "</div>";  // #controls
    echo "</div>";  // #vis-and-controls

    // Add a caption if it's set.
    if(!is_null($atts['caption'])) {
        echo "<figcaption>$atts[caption]</figcaption>";
    }

    echo "</div>";  // #collapsable

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
        $jsons[] = array(
            "url" => wp_get_attachment_url($json->ID),
            "file_name" => $json->post_name,
            "date_time" => $json->post_date
        );
    }

    // Add javascript for toggling the "collapsable" div.
    echo "<script>
              jQuery('#collapsable').hide();
              function togglevisibility() {
                  jQuery('#collapsable').show();
              }
          </script>";

    // If the user didn't include a 'file' attribute in the shortcode, display
    // a dropdown listing all JSON files from the media library.
    if (is_null($atts['file'])) {
        echo "<select id='file-input' onchange='togglevisibility();makeBrowserSim(viewer, \"3DMOLJS\", this.options[this.selectedIndex].value)'>
                  <option selected>Select file</option>";
        foreach ($jsons as $json){
            echo "<option value=$json[url]>$json[file_name] | $json[date_time]</option>";
        }
        echo '</select>';
    } else {
        // The user provided a 'file' attribute in the shortcode. Search if
        // the user-provided file name is a substring of any file in the media
        // library. Output javascript code to display that file once it is
        // found.
        foreach ($jsons as $json){
            if (strpos($json[file_name], preg_replace('/\\.[^.\\s]{3,4}$/', '', $atts['file'])) !== false){
                echo "<script>
                          togglevisibility();
                          makeBrowserSim(viewer, \"3DMOLJS\", \"$json[url]\");
                      </script>";
                break;
            }
        }
    }

    // Return all the above code as a string.
    return ob_get_clean();
}

/**
 * This method draws the plugin's settings page.
 *
 * @return void
 */
function durrant_options_page() { ?>
    <div>
        <?php screen_icon(); ?>
        <h2>PCAviz Plugin Settings</h2>
        <form method="post" action="options.php">
            <?php settings_fields( 'durrant_options_group' ); ?>
            <h3>Please allow us to collect minimal data on PCAViz usage.</h3>
            <p>
                As an academic group, our research depends on government grants.
                Being able to cite usage statistics in reports and articles will
                help us convince reviewers that our work is worth funding. <b>We
                understand and respect your privacy and that of your
                visitors.</b> Your detailed data will never be made public. We
                will only use summarized, aggregate data in our reports.
            </p>
            <p>
                Thank you! This helps us a lot!
            </p>
            <table>
                <tr valign="top">
                    <th scope="row">
                        <label for="durrant_option_name">Allow Analytics</label>
                    </th>
                    <td>
                        <input type="checkbox" id="durrant_option_name" name="durrant_option_name" value="checked" <?php echo get_option('durrant_option_name'); ?>/>
                    </td>
                </tr>
            </table>
            <?php  submit_button(); ?>
        </form>
    </div>
<?php }

// the settings are to be registered

/**
 * Register the settings.
 *
 * @return void
 */
function durrant_register_settings() {
   add_option( 'durrant_option_name', 'true');
   register_setting( 'durrant_options_group', 'durrant_option_name', 'durrant_callback' );
}
add_action( 'admin_init', 'durrant_register_settings' );

/**
 * Register the options page.
 *
 * @return void
 */
function durrant_register_options_page() {
  add_options_page('Page Title', 'PCAviz', 'manage_options', 'durrant', 'durrant_options_page');
}
add_action('admin_menu', 'durrant_register_options_page');

/**
 * Let wordpress know about our shortcode.
 *
 * @return void
 */
function durrant_shortcodes_init() {
    add_shortcode( 'pcaviz', 'pcaviz_main' );
}
add_action('init', 'durrant_shortcodes_init');

?>
