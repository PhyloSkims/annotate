#
# Compare CDS prediction to Reality
#

go_compare.sh reference predicted

  reference: Gbk/Embl file
  predicted: Gbk/Embl file

  output on stdout 

#
# Summarize comparisons
#
  
go_summarize.sh *.cmp

  *.cmp comparison files (output of go_compare.sh)
  
  output :
    compare.txt : space delimited format summary
    compare.pdf : graphics summary
    
