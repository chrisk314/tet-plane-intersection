
from setuptools import setup, Extension, find_packages

module_tetinter = Extension(
    'tet_plane_intersect.src.libwrap_tetinter',
    include_dirs=['./tet_plane_intersect/src'],
    libraries=[],
    library_dirs=[],
    extra_compile_args=['-std=c11', '-fopenmp'],
    extra_link_args=['-lgomp'],
    language='c',
    sources=['./tet_plane_intersect/src/plane_normal_tetrahedron_intersect.c']
)

setup(
    name='tet_plane_intersect',
    version='0.2',
    description='Calculates intersections between tetrahedra and planes.',
    packages=find_packages(),
    install_requires=[
        'numpy>=1.14.0',
    ],
    ext_modules=[module_tetinter],
    author='Chris Knight',
    author_email='chrisk314@gmail.com',
)
