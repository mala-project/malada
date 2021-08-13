# Installation of MALADA

## Python packages

In order to run MALA you have to have the following packages installed:

* ase
* numpy
* abc
* scipy
* xml

See also the `requirements.txt` file.
You can install each Python package with

```sh
$ pip install packagename
```

or all with

```sh
$ pip install -r requirements.txt
```

or just install the package along with dependencies

```sh
$ pip install -e .
```

(note: the `-e` is absolutely crucial, so that changes in the code will be
reflected system wide)
