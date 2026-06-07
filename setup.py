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
    # loose analysis modules (run as `python src/<name>.py`, also importable)
    py_modules=['mcis_paths', 'run_analysis', 'incremental_mcis', 'mcis_watch',
                'exact_ilp', 'conservation_track', 'spectral_mcis'],
    python_requires='>=3.10',
    install_requires=[
        'pandas>=1.5', 'numpy>=1.23', 'networkx>=3.0', 'pyarrow>=10.0',
        'scipy>=1.9', 'matplotlib>=3.5', 'seaborn>=0.12', 'pulp>=2.7',
    ],
    extras_require={'test': ['pytest>=7.0']},
    entry_points={'console_scripts': [
        'mcis-connectome=mcis_connectome.cli:main',
        'mcis-watch=mcis_watch:main',
    ]},
)
