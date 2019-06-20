<?php
/**
* Plugin Name: PCA viz Durran browser sim
*/

// we need to ask user in notification area to go to plugin' settings page and
// activate opt-in google-analytics option
register_activation_hook( __FILE__, 'pcaviz_plugin_activation' );
function pcaviz_plugin_activation() {
    $notices = get_option('pcaviz_plugin_deferred_admin_notices', array());
    $notices[] = "Please visit <a href='/wp-admin/options-general.php?page=durrant'>settings page</a>.";
    update_option('pcaviz_plugin_deferred_admin_notices', $notices);
}

add_action('admin_notices', 'pcaviz_plugin_admin_notices');
function pcaviz_plugin_admin_notices() {
    if($notices = get_option('pcaviz_plugin_deferred_admin_notices')) {
        foreach ($notices as $notice) {
            echo "<div class='updated'><p>$notice</p></div>";
        }
        delete_option('pcaviz_plugin_deferred_admin_notices');
    }
}

function display_notice( ) {?>

	<div class="notice notice-success is-dismissible">
		<p><?php _e('Please visit <a href="/wp-admin/options-general.php?page=durrant">plugin settings page</a> and activate statistics usage analysis!', 'shapeSpace'); ?></p>
	</div>

<?php
}

// the following is to notify wordpress' filters that we want to allow JSON
// format
function my_myme_types( $existing_mimes=array() ) {
    $existing_mimes['json'] = 'text/plain';
    return $existing_mimes;
}
add_filter( 'upload_mimes', 'my_myme_types' );

// registering necessary scripts
wp_register_script('3Dmol-min', plugin_dir_url(__FILE__) . 'assets/js/3Dmol-min.js');
wp_enqueue_script('3Dmol-min');
wp_register_script('BrowserSim.min', plugin_dir_url(__FILE__) . 'assets/js/BrowserSim.min.js');
wp_enqueue_script('BrowserSim.min');
wp_enqueue_script('jquery');
wp_register_script('GoogleAnalytics', plugin_dir_url( __FILE__ ) . 'assets/js/googleanalytics.js');
wp_enqueue_script('GoogleAnalytics');


// this method is doing most of the job
function browser_sim( $atts = [], $content = null, $tag = '' ){
    // a shortcode should never produce any output but return text! see https://stackoverflow.com/a/40500555/5907621 for details
    ob_start();
    // normalize attribute keys, lowercase
    $atts = array_change_key_case((array)$atts, CASE_LOWER);
    // if user agreed to share statistics with us...
    if(get_option('durrant_option_name') === 'checked'){
        echo "
            <script>
                ga('pcaviz.send', 'pageview');
            </script>
            ";
    }

    wp_register_style( 'viscontainer', plugin_dir_url(__FILE__) . 'assets/css/viscontainer.css');
    wp_enqueue_style('viscontainer');

    // if height attribute was set by user we overwrite the css setting
    if(!is_null($atts['height']))
        echo "
            <style>
            #vis-and-controls {
            height: $atts[height] !important;
            }
            </style>
                ";
    // if width attribute was set by user we overwrite the css default with it
    if(!is_null($atts['width']))
        echo "
            <style>
            #vis-and-controls {
                width: $atts[width] !important;
            }
            </style>
        ";

    // the canvas is wrapped into "collapsable" div; there's a script bellow
    // that toggles its visibility
    echo "
    <div id='collapsable' class='wp-block-image'";
        if(!is_null($atts['align']))
            echo "align=$atts[align]";
    echo ">
        <div id='vis-and-controls'>
             <div id='viscontainer'></div>
            <div id='controls'
            ";
    // we only display playback buttons if user mentioned them in shortcode
    if(!is_null($atts['playback_buttons']))
        if($atts['playback_buttons'] === 'false')
            echo "
                style='display: none'
                ";
    echo "
        ></div></div>";
    if(!is_null($atts['caption']))
        echo "<figcaption>$atts[caption]</figcaption>";
    echo "
        </div>
    <script type='text/javascript' src=";
echo plugin_dir_url(__FILE__);
echo "assets/js/viewer.js></script>";
    // collecting information on available JSON files in media library
    $query_json_args = array(
        'post_type'      => 'attachment',
        'post_mime_type' => 'text/plain',
        'post_status'    => 'inherit',
        'posts_per_page' => - 1,
    );
    $query_jsons = new WP_Query( $query_json_args  );
    $jsons = array();
    foreach ($query_jsons->posts as $json){
        $jsons[] = array("url" => wp_get_attachment_url($json->ID), "file_name" => $json->post_name, "date_time" => $json->post_date);
    } // now all JSON urls are in $jsons array
    // the script bellow toggles "collapsable" div visibility
    echo "
        <script>
            $('#collapsable').hide();
            function togglevisibility(){
                $('#collapsable').show();
            }
        </script>
    ";
    // if 'file' attribute wasn't provided with the shortcode we display
    // a dropdown listing all JSON files from medialibrary
    if (is_null($atts['file'])) {
        echo "
            <select id='file-input' onchange='togglevisibility();makeBrowserSim(viewer, \"3DMOLJS\", this.options[this.selectedIndex].value)'>
                <option selected>Select file</option>
        ";
        foreach ($jsons as $json){
            echo "
                <option value=$json[url]>$json[file_name] | $json[date_time]</option>
                ";
        }
        echo '</select>';
    } else {
        // searching if the user provided file name is a substring of any file
        // in medialibrary
        foreach ($jsons as $json){
            if(strpos($json[file_name], preg_replace('/\\.[^.\\s]{3,4}$/', '', $atts['file'])) !== false){
                echo "
                    <script>
                        togglevisibility();
                        makeBrowserSim(viewer, \"3DMOLJS\", \"$json[url]\");
                    </script>
                ";
                break;

            }
        }
    }
    return ob_get_clean();
}

// this method draws plugin' specific settings page
function durrant_options_page()
{
?>
  <div>
  <?php screen_icon(); ?>
  <h2>PCAviz Browser Sim Plugin</h2>
  <form method="post" action="options.php">
  <?php settings_fields( 'durrant_options_group' ); ?>
  <h3>Please allow us to collect some minimal data on PCA-VIZ plugin usage.</h3>
<p>We use analytics to mention it in futhure papers. <b>We don't spy on you!</b>
The collected data will never be made public. We understand and respect your privacy and that of your visitors.
The summarised information will be used in reports&articles. Thank you!</p>
  <p>This helps us a lot!</p>
  <table>
  <tr valign="top">
  <th scope="row"><label for="durrant_option_name">Allow Analytics</label></th>
  <td><input type="checkbox" id="durrant_option_name" name="durrant_option_name" value="checked" <?php echo get_option('durrant_option_name'); ?>/></td>
  </tr>
  </table>
  <?php  submit_button(); ?>
  </form>
  </div>
<?php
}

// the settings are to be registered
function durrant_register_settings() {
   add_option( 'durrant_option_name', 'true');
   register_setting( 'durrant_options_group', 'durrant_option_name', 'durrant_callback' );
}
add_action( 'admin_init', 'durrant_register_settings' );

function durrant_register_options_page() {
  add_options_page('Page Title', 'PCAviz', 'manage_options', 'durrant', 'durrant_options_page');
}
add_action('admin_menu', 'durrant_register_options_page');

// this funtion notifies wordpress of our shortcode
function durrant_shortcodes_init()
{
    add_shortcode( 'pcaviz', 'browser_sim' );
}
// the shortcode is to be initialized
add_action('init', 'durrant_shortcodes_init');

?>
