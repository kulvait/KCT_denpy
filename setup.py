from setuptools import setup

#Following https://uoftcoders.github.io/studyGroup/lessons/python/packages/lesson/

#Look here https://stackoverflow.com/questions/458550/standard-way-to-embed-version-into-python-package
exec(open('denpy/version.py').read())

petra_requires = ['h5py', 'pandas']

setup(
    # Needed to silence warnings (and to be a worthwhile package)
    name='denpy',
    url='https://github.com/kulvait/KCT_denpy',
    author='Vojtech Kulvait',
    author_email='vojtech.kulvait@hereon.de',
    # Needed to actually package something
    packages=['denpy', 'denpy.scripts'],
    # Needed for dependencies
    install_requires=['numpy','pydicom'] + petra_requires,
    # TIFF stacks to zarr / DEN (denpy.TIFF, convertTiffStack): pip install denpy[zarr]
    extras_require={'zarr': ['zarr>=3', 'tifffile', 'imagecodecs', 'numcodecs']},
    entry_points={'console_scripts': ['convertTiffStack=denpy.scripts.convertTiffStack:main']},
    # *strongly* suggested for sharing
    version=__version__,
    # The license can be anything you like
    license='GPL3',
    description='Python package for processing various tomographic data and DEN format.',
)
