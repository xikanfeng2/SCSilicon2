import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="SCSilicon2",
    version="1.0.1",
    author="Xikang Feng",
    author_email="fxk@nwpu.edu.cn",
    description="SCSilicon2: a single-cell genomics simulator which cost-effectively simulates single-cell genomics reads with haplotype-specific copy number annotation.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/xikanfeng2/SCSilicon2",
    project_urls={
        "Bug Tracker": "https://github.com/xikanfeng2/SCSilicon2/issues",
    },
    classifiers = [
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    package_dir={"": "src"},
    packages=setuptools.find_packages(where="src"),
    python_requires=">=3.6",
    install_requires=[
        'pandas>=0.23.4',
        'matplotlib>=3.0.2',
        'networkx>=3.2.1'
    ],
)
