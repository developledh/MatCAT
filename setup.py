from setuptools import setup, find_packages

setup(
    name="MatCat",
    version="1.0",
    packages=find_packages(exclude=["__pycache__"]),  
    install_requires=[],
    author="developledh",
    author_email="dongho0724@naver.com",
    description="MatCat: A module for automating VASP calculations",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/developledh/MatCAT/",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.8",
)
