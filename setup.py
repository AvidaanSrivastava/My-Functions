# my_utils/setup.py

from setuptools import setup, find_packages

setup(
    name='common_functions',
    version='0.1',
    packages=find_packages(),
    install_requires=[
        'numpy',
    ],
    author='Avidaan Srivastava',
    description='Astrophysics utility functions (e.g., semi-major axis, equilibrium temperature)',
)
