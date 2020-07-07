#!/bin/bash

rsync -rav "Snakefile" test/.;
rsync -rav "rules" test/.;
#rsync -rav "config" test/.;
rsync -rav "resources" test/.;

#rsync -rav --delete --exclude=".*" test/* mkreuzer@binfservms01.unibe.ch:/data/projects/p532_BX-rhizosphere-metagenomics/soil_metagenomics/test/. ;
rsync -rav --delete --exclude=".*" test/* mkreuzer@binfservms01.unibe.ch:/data/projects/p532_BX-rhizosphere-metagenomics/soil_metagenomics/test/. ;

rsync -rav "Snakefile" mkreuzer@binfservms01.unibe.ch:/data/projects/p532_BX-rhizosphere-metagenomics/soil_metagenomics/. ;
rsync -rav "rules" mkreuzer@binfservms01.unibe.ch:/data/projects/p532_BX-rhizosphere-metagenomics/soil_metagenomics/. ;
rsync -rav "config" mkreuzer@binfservms01.unibe.ch:/data/projects/p532_BX-rhizosphere-metagenomics/soil_metagenomics/. ;
rsync -rav "resources" mkreuzer@binfservms01.unibe.ch:/data/projects/p532_BX-rhizosphere-metagenomics/soil_metagenomics/. ;
