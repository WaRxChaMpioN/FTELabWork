#!/bin/bash
. $WM_PROJECT_DIR/etc/bashrc
(cd fluid && FOAM_ABORT=true pisoFoam -world fluid >& log.pisoFoam1)
