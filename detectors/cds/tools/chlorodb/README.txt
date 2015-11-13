#
# Make chloro DB's for use by CDS prediction program
#

mkdir myDB
mkdir myDB/download
cp GBKorEMBLEntries myDB/download
[optional
 cp ORG_HOME/data/cds/chlorodb/parameters.sh myDB
 <edit> myDB/parameters.sh
]

go_chlorodb.sh myDB

and... wait... (typical time : 5h for 500 Entries)

then (after checking) 

  mv ORG_HOME/data/cds/chlorodb ORG_HOME/data/cds/chlorodb.old
  mv myDB ORG_HOME/data/cds/chlorodb

#
# notes
#

this requires an installed R > 3.0.1

with the following packages

    igraph          # <- mandatory

    grid            # <- the following are not needed
    gridExtra       #    by scripts, but just to
    vcd             #    produce graphics for models
    plotrix         #    by lib/plot.models.r

you can check and install them by running :

 lib/install.rpackages.r 
