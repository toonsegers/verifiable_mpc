import pathlib
from setuptools import setup

# The directory containing this file
HERE = pathlib.Path(__file__).parent

# The text of the README file
README = (HERE / "README.md").read_text()

# This call to setup() does all the work
setup(
    name="verifiable-mpc",
    version="0.0.1",
    description="Implements the verfiable MPC scheme.",
    long_description=README,
    long_description_content_type="text/markdown",
    url="https://github.com/toonsegers/verifiable_mpc",
    author="Toon Segers",
    author_email="a.j.m.segers@tue.nl",
    license="MIT",
    classifiers=[
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
    ],
    packages=["verifiable_mpc", "verifiable_mpc.tools", "verifiable_mpc.ac20", "verifiable_mpc.trinocchio"],
    include_package_data=True,
    # install_requires=["mpyc", "secgroups"],
    python_requires='>=3.6',
)
