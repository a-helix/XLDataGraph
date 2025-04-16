from setuptools import setup, find_packages

setup(
    name='xldg',
    version='0.3.0',
    description='"XLDataGraph is a library for crosslinking data analysis and visualization.',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    author='a-helix',
    author_email='sorokin.biochemistry@gmail.com',
    url='https://github.com/a-helix/XLDataGraph',
    packages=find_packages(where='xldg', exclude=['tests*']),
    package_dir={'': 'xldg/src'},
    install_requires=[
        'requests>=2.32.3',
        'pyCirclize>=1.6.0',
        'matplotlib_venn>=1.1.2',
    ],
    classifiers=[
        'Programming Language :: Python :: 3',
        'Operating System :: OS Independent',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
    ],
    python_requires='>=3.7',
    license='GPL-3.0',
    test_suite='tests',
)
