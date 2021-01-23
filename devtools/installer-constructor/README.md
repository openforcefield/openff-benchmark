## Summary

This repository contains the build script for the Open Force Field partner benchmarking conda deployment.

_Special thanks to @loriab and the Psi4 community, as this is basically a copy of [their build system](https://github.com/psi4/psi4meta/tree/master/conda-recipes/constructor-cutter-unified)_


## Usage


* run in an environment that has `cookiecutter` and `constructor` installed.

* edit `cookiecutter/cookiecutter.json` for control. edit `cookiecutter/{{.../construct.yaml` for templating

* dir `build/` is regenerated each time.

* may want to clear out ~/.conda/constructor

* `python run.py`

* watch out for `py_` in buildstring as this means a noarch, which can be troublesome
