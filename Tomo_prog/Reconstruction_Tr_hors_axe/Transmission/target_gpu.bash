#!/bin/bash

make cleanobj
rm -f Makefile
ln -s Makefile-cuda Makefile
