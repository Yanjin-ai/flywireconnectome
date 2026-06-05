from setuptools import setup, find_packages

setup(
    name='mcis_connectome',
    version='1.0.0',
    description='Maximum Common Induced Subgraph across Drosophila connectomes',
    author='Yanjin Li',
    author_email='geraldineliyanjin@gmail.com',
    url='https://github.com/Yanjin-ai/flywireconnectome',
    packages=find_packages('src'),
    package_dir={'': 'src'},
    python_requires='>=3.9',
    install_requires=['pandas>=1.5', 'numpy>=1.23', 'networkx>=3.0', 'pyarrow>=10.0'],
    entry_points={'console_scripts': ['mcis-connectome=mcis_connectome.cli:main']},
)
