# This BASH script compiles all components of PCAViz.

############ Compile the JavaScript ############
echo "Compiling JavaScript Interpreter"
cd ../pcaviz-interpreter-javascript/

echo "    Remove old files."
rm PCAViz.js PCAViz.min.js closure.errs.txt

echo "    Compile TypeScript."
tsc --target ES5 PCAViz.ts

echo "    Compile JavaScript with closure compiler."
# export formatting="--formatting=PRETTY_PRINT"
export formatting=""
java -jar ../dev_tools/utils/closure-compiler-v20180506.jar $formatting \
   --compilation_level=ADVANCED_OPTIMIZATIONS \
   --externs='../dev_tools/utils/custom_extern.js' --js_output_file='PCAViz.min.js' 'PCAViz.js' \
   2> closure.errs.txt

echo "    Display and delete closure errors."
cat closure.errs.txt
rm closure.errs.txt

echo "    Preserve PCAVizNameSpace global namespace."
python ../dev_tools/utils/fix_namespaces.py

echo "    Copy js file to examples directory."
cp PCAViz.min.js examples/

echo "    Make ZIP of JavaScript directory."
cd ..
rm pcaviz-interpreter-javascript.zip
zip -r pcaviz-interpreter-javascript.zip pcaviz-interpreter-javascript/
cd dev_tools

############ Compile WordPress Plugin ############
echo "Compiling WordPress Plugin."
cd ../pcaviz-wordpress/

echo "    Copy compiled PCAViz.min.js to the wordpress directory."
cp ../pcaviz-interpreter-javascript/PCAViz.min.js ./assets/js/

# Remake wordpress plugin ZIP file.
echo "    Make ZIP file of wordpress directory."
cd ../
rm pcaviz-wordpress.zip
zip -r pcaviz-wordpress.zip pcaviz-wordpress
cd dev_tools

############ Compile Pytthon Compressor ############
echo "Compiling Python Compressor."
cd ../pcaviz-compressor-python/

echo "    Delete any cruft test files."
ls examples/* | grep -v "\.dcd" | grep -v "\.pdb" | grep -v "\.psf" | awk '{print "rm " $1}' | bash

echo "    Make ZIP file of Python Compressor."
cd ../
rm pcaviz-compressor-python.zip
zip -r pcaviz-compressor-python.zip pcaviz-compressor-python
cd -
