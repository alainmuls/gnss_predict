from setuptools import setup, find_packages

setup(
    name="gnss_tracking",
    version="1.0",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    package_data={
        "gnss_tracking": ["GNSS/*", "tle/*", "utils/*"]
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
            "gnss_tracking=gnss_tracking.main:main",
        ],
    },
)
