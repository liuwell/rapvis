import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="rapvis", # Replace with your own username
    version="0.0.1",
    author="Wei Liu",
    author_email="liuwei_hnnd@163.com",
	keywords=['RNAseq', 'Visualization', 'different Expressed','GO', 'Enrichment',
          'Bioinformatics', 'Computational Biology',],
    description="A tool for RNAseq processing and visualization",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/liuwell",
    packages=setuptools.find_packages(),
	install_requires=[
                        'numpy',
                        'scipy',
                        'pandas',
                        'matplotlib',
                        'seaborn',
                        'gseapy',
                        ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
