"""
BAP: Bead-based scATAC-seq data Processing
"""
from setuptools import find_packages, setup
from distutils.core import setup, Extension

dependencies = ['click', 'pytest', 'snakemake', 'optparse-pretty', 'regex', 'pysam', 'ruamel.yaml']

setup(
    name='bap',
    version='0.3.6',
    url='https://github.com/caleblareau/bap',
    license='MIT',
    author='Caleb Lareau',
    author_email='clareau@broadinstitute.org',
    description='Bead-based scATAC-seq data Processing.',
    long_description=__doc__,
    packages=find_packages(exclude=['tests']),
    include_package_data=True,
    zip_safe=False,
    platforms='any',
    install_requires=dependencies,
    entry_points={
        'console_scripts': [
            'bap = bap.cli:main',
            'annotateAlleles = bap.snpAnnotate.cli_anno:main'
        ],
    },
    classifiers=[
        # As from http://pypi.python.org/pypi?%3Aaction=list_classifiers
        # 'Development Status :: 1 - Planning',
        # 'Development Status :: 2 - Pre-Alpha',
         'Development Status :: 3 - Alpha',
        # 'Development Status :: 4 - Beta',
        # 'Development Status :: 5 - Production/Stable',
        # 'Development Status :: 6 - Mature',
        # 'Development Status :: 7 - Inactive',
        'Environment :: Console',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Operating System :: POSIX',
        'Operating System :: MacOS',
        'Operating System :: Unix',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Topic :: Software Development :: Libraries :: Python Modules',
    ]
)
