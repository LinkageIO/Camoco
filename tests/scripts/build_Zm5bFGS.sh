#!/bin/sh

camoco build-refgen \
    ../raw/RefGen/ZmB73_5b_FGS.gff.gz \
    Zm5bFGS \
    'Maize RefGen' \
    5b \
    Zeamays \
    --chrom-feature chromosome \
    --ID-attr ID \
    --gene-feature gene
