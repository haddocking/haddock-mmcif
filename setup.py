from setuptools import setup, find_packages

setup(
    name="haddock2mmcif",
    license="Apache License 2.0",
    version="0.1.0",
    author="Rodrigo Honorato",
    description="Encode a HADDOCK run to mmCIF",
    author_email="",
    include_package_data=True,
    packages=find_packages("src"),
    package_dir={"": "src"},
    classifiers=[],
    python_requires=">=3.9, <4",
    install_requires=["ihm"],
    entry_points={
        "console_scripts": [
            "haddock2mmcif=haddock2mmcif.cli:main",
        ],
    },
)
