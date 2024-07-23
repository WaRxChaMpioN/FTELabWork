#!/bin/bash
. $WM_PROJECT_DIR/etc/bashrc
(cd solid && FOAM_ABORT=true pisoFoam -world solid -parallel)
