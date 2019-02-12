#!/usr/bin/env bash

# This puts everything we want into the zip file
# Should be run from the parent directory of TFOCS

#PACKAGE='TFOCS_v1.0a.zip'
PACKAGE='TFOCS_v1.0aa.zip'

zip $PACKAGE -r TFOCS/* --exclude \*TFOCS/hide/* \*TFOCS/experiments/* \
    \*conflicted* \*old_examples/* \*.sh \*.m~ \*.zip \*TFOCS/userguide/* \
    \*solver_sNuclerBP_packSVD.m \*prox_nuclearP.m \*TFOCS/@packSVD/*

zip -u $PACKAGE --junk-paths TFOCS/userguide/tfocs_userguide.pdf
