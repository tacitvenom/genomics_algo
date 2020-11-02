import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="genomics_algo",
    version="0.0.1",
    author="Vasu Sharma",
    author_email="vasu.sharma314@gmail.com",
    description="A python package with algorithms relevant for DNA sequencing",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/tacitvenom/genomics_algo",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.6",
)
