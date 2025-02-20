from setuptools import setup, find_packages

setup(
    name="acedrg",
    version="315",
    packages=find_packages(),
    install_requires=["rdkit", "servalcat",
                      "aceDRG-tables @ git+https://github.com/flong-mrc/aceDRG-tables.git"],  # Add dependencies if needed
    entry_points={
        "console_scripts": [
            "acedrg=acedrg.acedrg:main",
        ],
    },
)