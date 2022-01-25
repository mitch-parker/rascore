import setuptools
from pathlib import Path

this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

setuptools.setup(
    name="rascore",
    version="0.1.19",
    author="mitch-parker",
    author_email="mitch.isaac.parker@gmail.com",
    description=" tool for analyzing the conformations of RAS structures",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/mitch-parker/rascore",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    package_dir={"": "src"},
    packages=setuptools.find_packages(where="src"),
    python_requires=">=3.7.1,<3.10",
)