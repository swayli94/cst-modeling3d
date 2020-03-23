from setuptools import setup, find_packages
from cst_modeling import __version__

with open('Readme.md') as f:
      long_description = f.read()

setup(name='cst_modeling3d',
      version=__version__,
      description='This is the module of surface/airfoil modeling',
      long_description=long_description,
      keywords='CST modeling',
      download_url='https://github.com/swayli94/cst-modeling3d/',
      license='MIT',
      author='Runze LI',
      author_email='swayli94@gmail.com',
      packages=find_packages(exclude=['example']),
      install_requires=['numpy', 'scipy', 'matplotlib'],
      classifiers=[
            'Programming Language :: Python :: 3'
      ]
)

