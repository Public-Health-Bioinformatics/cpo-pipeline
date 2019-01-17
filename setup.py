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
    scripts=[
        'bin/cpo-pipeline',
        'bin/cpo-multi',
    ],
    zip_safe=False,
)
