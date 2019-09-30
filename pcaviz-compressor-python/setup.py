import setuptools
# import glob

with open("PCAViz/README.md", "r") as fh:
    long_description = fh.read()

# extra_files = glob.glob("*.md") + glob.glob("examples/*")

setuptools.setup(
    name="pcaviz-durrantlab",
    version="1.7",
    author="Jacob Durrant",
    author_email="durrantj@pitt.edu",
    description="A script to compress molecular dynamics simulations so they can be visualized with PCAViz in a web browser.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="http://git.durrantlab.com/jdurrant/pcaviz",
    packages=setuptools.find_packages(),
    package_data={'PCAViz': ['*.md', 'examples/*']},
    classifiers=[
        "License :: OSI Approved :: GNU General Public License v2 (GPLv2)",
        "Operating System :: OS Independent",
    ],
    scripts=['bin/pcaviz']
)
