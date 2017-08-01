#!/usr/bin/env bash

autoreconf --install && \
cd build && \
../configure && \
make
