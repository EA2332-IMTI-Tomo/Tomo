#!/bin/bash

cnet
rm *.sparcf
cmucl -eval "(declaim-debug)

             (compile-file \"Tools.lisp\")
             (mapcar #'compile-file (directory \"./\*.lisp\"))
             (quit)"


