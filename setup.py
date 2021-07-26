from setuptools import setup

# Doing it as suggested here:
# https://packaging.python.org/guides/single-sourcing-package-version/
# (number 3)

version = {}
with open("malada/version.py") as fp:
    exec(fp.read(), version)

setup(
    name="malada",
    version=version["__version__"],
    description="MALA Data acquisition",
    url="https://github.com/mala-project/malada",
    author="Lenz Fiedler",
    author_email="l.fiedler@hzdr.de",
    license="MIT",
    packages=["malada"],
    zip_safe=False,
    install_requires=open('requirements.txt').read().splitlines(),
    python_requires='<3.9',
)
