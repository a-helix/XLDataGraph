from setuptools import setup, find_packages

setup(
    name='xldg',
    version='0.2.2',
    description='A Python library for protein-protein interaction analysis',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    author='a-helix',
    author_email='your.email@example.com',
    url='https://github.com/a-helix/XLDataGraph',
    packages=find_packages(where='xldg', exclude=['tests*']),
    package_dir={'': 'xldg/src'},
    install_requires=[
        'requests>=2.32.3',
        'pyCirclize>=1.6.0',
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
