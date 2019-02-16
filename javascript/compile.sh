# Typescript compile
echo "Compile typescript."
tsc --target ES5 BrowserSim.ts

# Closure compile
echo "Closure compile."
# export formatting="--formatting=PRETTY_PRINT"
export formatting=""

#java -jar utils/closure-compiler-v20180506.jar $formatting \
#   --compilation_level=ADVANCED_OPTIMIZATIONS --externs='utils/jquery-1.9.js' \
#   --externs='utils/custom_extern.js' --js_output_file='BrowserSim.min.js' 'BrowserSim.js' \
#   2> closure.errs.txt

# If you don't closure compile...
cp BrowserSim.js BrowserSim.min.js

# Make example files
cd examples
python make.py

#echo "Start server."
python -m SimpleHTTPServer 8000
