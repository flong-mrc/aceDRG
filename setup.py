from setuptools import setup, find_packages

setup(
    name="acedrg",
    version="315",
    packages=find_packages(),
    install_requires=["rdkit", "servalcat"],  # Add dependencies if needed
    entry_points={
        "console_scripts": [
            "acedrg=acedrg.acedrg:main",
        ],
    },
)