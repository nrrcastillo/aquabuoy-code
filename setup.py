from setuptools import setup, find_packages

setup(
    name="aquabuoy_code",
    version="1.0",
    packages=find_packages(where="aquabuoy"),
    package_dir={"": "aquabuoy"},
    install_requires=["numpy", "scipy", "matplotlib", "pandas"],
    description="Python port for the original MATLAB code used in the Aquabuoy project",
    url="https://github.com/nrrcastillo/aquabuoy-code",
    author="Nicolas Castillo",
    author_email="ncastillo@seattleu.edu"
)