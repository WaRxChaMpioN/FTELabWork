#!/bin/bash
. $WM_PROJECT_DIR/etc/bashrc
(cd fluid && FOAM_ABORT=true interFoam -world fluid -parallel)
