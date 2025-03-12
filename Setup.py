from setuptools import setup, find_packages

setup(
    name='xldg',
    version='0.2.2',
    description='A Python library for protein-protein interaction analysis',
    author='a-helix',
    author_email='your.email@example.com',
    url='https://github.com/a-helix/XLDataGraph',
    packages=find_packages(where='xldg', exclude=['tests*']),
    package_dir={'xldg': 'xldg/src/xldg'},
    install_requires=[
        'requests>=2.32.3',
        'pyCirclize>=1.9.0',
    ],
    classifiers=[
        'Programming Language :: Python :: 3',
        'Operating System :: OS Independent',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
    ],
    python_requires='>=3.7',
    license='GPL-3.0',
)