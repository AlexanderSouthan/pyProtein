from setuptools import setup, find_packages

setup(
    name='pyProtein',
    version='0.0.4',
    packages=find_packages(where='src'),
    install_requires=['numpy', 'pandas', 'scipy']
)
