# Clean old files
rm PCAViz.js PCAViz.min.js closure.errs.txt

# Typescript compile
echo "Compile typescript."
tsc --target ES5 PCAViz.ts

# Closure compile
echo "Closure compile."
# export formatting="--formatting=PRETTY_PRINT"
export formatting=""

# --externs='utils/jquery-1.9.js'
java -jar utils/closure-compiler-v20180506.jar $formatting \
   --compilation_level=ADVANCED_OPTIMIZATIONS \
   --externs='utils/custom_extern.js' --js_output_file='PCAViz.min.js' 'PCAViz.js' \
   2> closure.errs.txt

#cp PCAViz.js PCAViz.min.js

# Fix a namespace issue.
python utils/fix_namespaces.py

# Make example files
#echo "Make example files."
#cd examples
#python make.py

cp PCAViz.min.js examples/

echo "Start server."
python -m SimpleHTTPServer 8000
