import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="EnrichR",
    version="0.0.1",
    author="Sergei RYBALKO",
    author_email="benja1972@gmail.com",
    description="Enrichment analysis with clustering",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/Benja1972/bioinformatics",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
