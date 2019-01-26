from distutils.core import setup
from setuptools import find_packages
from cpo_pipeline import __version__

setup(
    name='cpo_pipeline',
    version=__version__,
    description='An analysis pipeline for the purpose of investigating Carbapenemase-Producing Organisms.',
    url='http://github.com/Public-Health-Bioinformatics/cpo-pipeline',
    author='Dan Fornika, Justin Jia, Matthew Croxen, William Hsiao',
    author_email='dan.fornika@bccdc.ca, bja20@sfu.ca, matthew.croxen@bccdc.ca, william.hsiao@bccdc.ca',
    license='MIT',
    packages=['cpo_pipeline'],
    include_package_data=True,
    keywords = "molecular epidemiology",
    project_urls = {
        "Bug Reports": "https://github.com/Public-Health-Bioinformatics/cpo-pipeline/issues",
        "Change Log":  "https://github.com/Public-Health-Bioinformatics/cpo-pipeline/CHANGES.md",
        "Source":      "https://github.com/Public-Health-Bioinformatics/cpo-pipeline",
    },
    scripts=[
        'bin/cpo-pipeline',
    ],
    zip_safe=False,
    python_requires = '>=3.4',
    install_requires = [
        "drmaa >=0.7.9, ==0.7.*",
    ],
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",  
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.4",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
        
    ],
)
