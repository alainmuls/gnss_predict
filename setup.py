from setuptools import setup, find_packages

setup(
    name="gnss_predict",
    version="1.0",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    package_data={
        "gnss_predict": ["GNSS/*", "tle/*", "utils/*"]
    },
    include_package_data=True,
    install_requires=[
        "ephem",
        "numpy",
        "matplotlib",
        "cartopy",
        "rich",
        "pandas",
    ],
    entry_points={
        "console_scripts": [
            "gnss_predict=gnss_predict.main:main",
        ],
    },
)
