import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="pyfb",
    version="0.0.1",
    author="Daniel G. Trabada",
    description="python scripts for fireball",
    long_description=long_description,
    url="https://github.com/fireball-QMD/progs/pyfb",
    long_description_content_type="text/markdown",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: GPLv3",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3',
)
