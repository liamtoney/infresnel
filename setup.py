from setuptools import find_packages, setup

setup(
    name='infresnel', packages=find_packages(), install_requires=['pygmt', 'rioxarray']
)
