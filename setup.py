"""setup file for pyinharmc"""
from setuptools import setup, find_packages

setup(
    name="pyinharmc",
    version="0.0.1",
    author="Koutoupis Stefanos",
    author_email="stefanoskoutoupis@gmail.com",
    description="A python package for inharmonicity computation",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires=[],
)
