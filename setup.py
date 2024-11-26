from setuptools import setup

setup(
    name="neurodamus-sonata-template",
    version="0.1.0",
    description="A simple SONATA project",
    author="Kerem Kurban",
    author_email="kerem.kurban@epfl.ch",
    url="https://github.com/KeremKurban/neurodamus-sonata-template",
    keywords=["sonata", "neurodamus", "simulation"],
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
    ],
    install_requires=[
        "h5py",
    ],
)
