#
# summarize comparison
#
# run after go_compare.sh
#
# usage: go_summarize.sh *.cmp
#
#
unsetenv ORG_SOURCED

setenv ORG_HOME `dirname $0`/../../../..
source $ORG_HOME/scripts/csh_init.sh

NeedArg 1

egrep '^#|^MATCH' $* | awk -f $LIB_DIR/summary.cmp.awk > compare.txt

$LIB_DIR/summarize_cmp.r 



exit 0


