import setuptools

with open("readme.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="superfit", # Replace with your own username
    version="0.0.1",
    author="Samantha Goldwasser, Ido Irani",
    author_email="Samantha.goldwasser@weizmann.ac.il",
    description="Classify supernova spectra according to Howell et al. 2004",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/samanthagoldwasser25/superfit",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)