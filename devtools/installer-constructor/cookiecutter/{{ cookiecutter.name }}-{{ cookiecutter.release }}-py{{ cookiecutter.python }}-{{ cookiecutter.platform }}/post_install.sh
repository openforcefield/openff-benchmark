#!/bin/sh

# Explicitly move noarch packages into `lib/python?.?/site-packages` as a
# workaround to [this issue][i86] with lack of `constructor` support for
# `noarch` packages.
#
# [i86]: https://github.com/conda/constructor/issues/86#issuecomment-330863531
if [[ -e site-packages ]]; then
    cp -r site-packages/* $PREFIX/lib/python{{ cookiecutter.python }}/site-packages
    rm -r site-packages
fi
