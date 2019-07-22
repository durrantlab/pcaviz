cd ../javascript/
./compile.sh
cd -
cp ../javascript/BrowserSim.min.js ./pcaviz/assets/js/BrowserSim.min.js
rsync -rv pcaviz durrantj@durrantlab.pitt.edu:/var/www/html/wp-content/plugins/

