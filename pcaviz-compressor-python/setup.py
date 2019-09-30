import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="pcaviz-durrantlab",
    version="1.4",
    author="Jacob Durrant",
    author_email="durrantj@pitt.edu",
    description="A script to compress molecular dynamics simulations so they can be visualized with PCAViz in a web browser.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="http://git.durrantlab.com/jdurrant/pcaviz",
    packages=setuptools.find_packages(),
    classifiers=[
        "License :: OSI Approved :: GNU General Public License v2 (GPLv2)",
        "Operating System :: OS Independent",
    ],
    scripts=['pcaviz']
)
