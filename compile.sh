# Typescript compile
echo "Compile typescript."
tsc --target ES5 main.ts

# Closure compile
echo "Closure compile."
# export formatting="--formatting=PRETTY_PRINT"
export formatting=""

java -jar utils/closure-compiler-v20180506.jar $formatting --compilation_level=ADVANCED_OPTIMIZATIONS --externs='utils/jquery-1.9.js' --externs='utils/twitter-bootstrap-2.1.1-externs.js' --externs='utils/custom_extern.js' --js_output_file='main.min.js' 'main.js' 2> closure.errs.txt

#echo "Start server."
#python -m SimpleHTTPServer 8000
