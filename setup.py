from setuptools import setup
import subprocess


with open('README.md', 'r') as f:
    long_description = f.read()


setup(
   name='sahara_work',
   version='0.1.0',
   description='Criticality implementation + animal databases.',
   license="",
   long_description=long_description,
   package_dir={'sahara_work': '.'},
   author='\
           (Hengen Lab Washington University in St. Louis)',
   author_email='',
   maintainer='Kiran Bhaskaran-Nair,\
           (Hengen Lab Washington University in St. Louis)',
   maintainer_email='',
   url="https://github.com/hengenlab/sahara_work",
   download_url="https://github.com/hengenlab/sahara_work",
   packages=['sahara_work'],
   install_requires=[],
   classifiers=[
        'Development Status :: 1 - Pre-Alpha',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
    ],
   scripts=[]
)
