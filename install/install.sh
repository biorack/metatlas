#!/usr/bin/env bash


DIRNAME=$(python -c "import imp; print(imp.find_module('pymzml')[1])")

cp psi-ms-1.2.0.obo $DIRNAME/obo

cp obo.patch $DIRNAME
pushd $DIRNAME
patch obo.py < obo.patch
rm obo.patch
popd
