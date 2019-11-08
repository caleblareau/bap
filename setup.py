"""
bap: Bead-based scATAC-seq data Processing
"""
from setuptools import find_packages, setup
from distutils.core import setup, Extension

dependencies = ['python-Levenshtein','biopython','fuzzysearch','fuzzywuzzy','click', 'pytest', 'snakemake', 'optparse-pretty', 'regex', 'pysam', 'ruamel.yaml']

setup(
    name='bap-atac',
    version='0.6.4a',
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
            'bap = bap.cli_old_dontEdit:main',
            'bap-barcode = bap.barcode.cli_barcode:main',
            'bap2 = bap.cli_bap2:main',
            'bap-frag = bap.cli_bap_frag:main',
            'bap-bulk-frag = bap.cli_bap_bulk_frag:main',
            'bap-reanno = bap.cli_reanno:main'
        ],
    },
    classifiers=[
        # As from http://pypi.python.org/pypi?%3Aaction=list_classifiers
        # 'Development Status :: 1 - Planning',
        # 'Development Status :: 2 - Pre-Alpha',
        # 'Development Status :: 3 - Alpha',
         'Development Status :: 4 - Beta',
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
