from setuptools import setup

# What you enter as the name variable of the setup 
# is what you will use to import your environment
# (for eg. here, import samo_database).

setup(name='cst_modeling',
      version='0.1',
      description='This is the module of surface/airfoil modeling',
      author='Runze LI',
      author_email='swayli94@gmail.com',
      install_requires=['copy','numpy','scipy','matplotlib']
)
