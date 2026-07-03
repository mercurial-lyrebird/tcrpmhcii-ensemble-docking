
from setuptools import setup, find_packages


with open("README.md", "r") as f:
    description = f.read()


if __name__ == "__main__":

    setup(
        name="pmhclib",
        description="TODO",
        long_description=description,
        long_description_content_type="text/markdown",
        author="TODO",
        version="0.0.0",
        packages=find_packages()
    )
