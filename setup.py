from setuptools import find_packages, setup

import versioneer

setup(
    name='infresnel',
    packages=find_packages(),
    install_requires=[
        'colorcet',
        'ipympl',
        'ipywidgets',
        'jupyterlab',
        'matplotlib',
        'numba',
        'pip',
        'pygmt',
        'rioxarray',
        'tqdm',
    ],
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
)
