# Closure compiler takes the top-level namespace and converts it into a
# single-letter variable. This variable can potentially clash with other
# libraries. So we need to process the closure output to fix this problem.
# Could find a way to assign a namespace to the window like
# window["namespace"] to survive closure. This is a hackish but functional
# solution.

import re

js = open("PCAViz.min.js", "r").read()

# First, get the variable name that closure chose.
var_name = re.findall(r"^var (.);", js)[0]

# Now replace it with PCAVizNameSpace
js = js.replace("var " + var_name + ";", "var PCAVizNameSpace;")
js = js.replace(
    "(" + var_name + "||(" + var_name + "={}))",
    "(PCAVizNameSpace||(PCAVizNameSpace={}))"
)

open("PCAViz.min.js", "w").write(js)
