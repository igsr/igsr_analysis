import setuptools

# read the contents of your README file
from os import path
this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, 'README4pypi.md'), encoding='utf-8') as f:
    long_description = f.read()


setuptools.setup(
    name='igsr_analysis',
    version='1.1.1',
    description='Code that is relevant for the analysis (Mapping, BAM qc, Variant Calling, Filtering etc...) of IGSR data',
    license="Apache License 2.0",
    author='Ernesto Lowy',
    author_email='ernestolowy@gmail.com',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url="https://github.com/igsr/igsr_analysis",
    packages=setuptools.find_packages(),
    install_requires=['pandas','PyMySQL','scipy','numpy','matplotlib','matplotlib-venn'], #external packages as dependencies
)
