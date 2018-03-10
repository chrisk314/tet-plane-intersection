from distutils.core import setup
from setuptools import find_packages

setup(
    name='tet_plane_intersect',
    version='0.1',
    description='Calculates intersections between tetrahedra and planes.',
    packages=find_packages(),
    install_requires=[
        'numpy>=1.14.0',
    ],
    author='Chris Knight',
    author_email='chrisk314@gmail.com',
)
