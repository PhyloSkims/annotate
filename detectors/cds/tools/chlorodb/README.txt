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

calculation of models currently requires R
(without any specific package) this will be
replaced in the future...

in addition, optional graphics output (plot.models.r) 
requires the following graphic packages :

    grid           
    gridExtra      
    vcd            
    plotrix        

