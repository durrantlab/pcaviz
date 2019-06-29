echo USAGE:   ./compress_decompress.sh {cum_var} {precision}
echo EXAMPLE: ./compress_decompress.sh 0.8 2
echo

export cum_var=${1}
export precision=${2}

# Make the JSON file
python ../python/PCA.py --top_file ../python/test_files/LARP1_first_frame_top.pdb \
    --coor_file ../python/test_files/LARP1_test_traj.pdb --cum_var ${cum_var} \
    --selection "not name H*" --precision ${precision} --check_accuracy

# Move the JSON file to the output directory, renaming it.
mkdir -p output
export new_json_filename="./output/LARP1_test_traj__${cum_var}__${precision}.json"
mv ../python/test_files/LARP1_test_traj.compressed.json ${new_json_filename}
echo JSON moved to ${new_json_filename}

# Report size of JSON file.
echo -n "JSON file size: "
ls -lth ${new_json_filename} | awk '{print $5}'

# Now run the JSON through the javascript. To run javascript from the
# commandline, you need to install node: https://nodejs.org/en/download/
node ../javascript/BrowserSim.js ${new_json_filename} > ${new_json_filename}.pdb
echo 3D coordinates from JSON file saved to ${new_json_filename}.pdb

# Make a quick VMD visualization, just so you can visually inspect the
# awesomeness.
cat vis.vmd.template | sed "s|NEWPDB|${new_json_filename}.pdb|g" > vis.vmd
echo To compare the two PDBs visually, use vis.vmd.
