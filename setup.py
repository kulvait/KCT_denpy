from setuptools import setup

#Following https://uoftcoders.github.io/studyGroup/lessons/python/packages/lesson/

#Look here https://stackoverflow.com/questions/458550/standard-way-to-embed-version-into-python-package
exec(open('denpy/version.py').read())

setup(
    # Needed to silence warnings (and to be a worthwhile package)
    name='denpy',
    url='https://github.com/kulvait/KCT_denpy',
    author='Vojtech Kulvait',
    author_email='vojtech.kulvait@ovgu.de',
    # Needed to actually package something
    packages=['denpy'],
    # Needed for dependencies
    install_requires=['numpy','pydicom'],
    # *strongly* suggested for sharing
    version=__version__,
    # The license can be anything you like
    license='GPL3',
    description='Python package for DEN and DICOM processing for perfusion vizualization.',
)
