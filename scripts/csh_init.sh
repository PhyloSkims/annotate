#
# csh setup
# usage: should normaly be sourced from other csh scripts as :
# usage: source `dirname $0`/csh_setup.csh
#

if ($?ORG_SOURCED == 0) then

  setenv ORG_SOURCED 1
  setenv ORG_VERSION 1
  
  # --------------------------------------
  # adapt one of these if needed
  # --------------------------------------
  #
  # o AwkCmd  : name of awk command (on some system should be 'gawk')
  #
  # o Verbose : default verbosity (may be changed by -v option)
  #
  
  setenv AwkCmd        "gawk"
  setenv Verbose       0
  
  # --------------------------------------
  # don't change hereafter (normally)
  # --------------------------------------

  # ORG_HOME : the absolute path to the ORG.Annot home directory
  # note: should be set before usually, the following will work
  #       only for scripts located in /scripts
  
  if ($?ORG_HOME == 0) setenv ORG_HOME `dirname $0`/..
  setenv ORG_HOME `cd $ORG_HOME && pwd -P`
  
  setenv ORG_PORTNAME `$ORG_HOME/config/guess_port`	 # The architecture running
                                                     # the ORG.Annot instance

  setenv BIN_DIR "$ORG_HOME/ports/$ORG_PORTNAME/bin" # Directory containing 
                                                     # binaries for this port

  setenv SCRIPT_DIR "$ORG_HOME/scripts"              # Directory containing
                                                     # scripts utilities

  setenv PROG_DIR `dirname $0`		                 # Directory containing
  setenv PROG_DIR `cd $PROG_DIR && pwd -P`           # the main script file

  setenv LIB_DIR "$PROG_DIR/../lib"                  # Directory containing
  setenv LIB_DIR `cd $LIB_DIR && pwd -P`             # the main script libraries

  setenv CALL_DIR `pwd -P`                           # Directory from where the
                                                     # main script is called

  setenv DATA_DIR "$ORG_HOME/data"       # Directory containing reference data
                                         # for the annotation

  setenv IR_DATA_DIR "$DATA_DIR/ir"  	 # Directory containing data related to
										 # IRs detection

  setenv TRNA_DATA_DIR "$DATA_DIR/trna"  # Directory containing data related to
										 # tRNAs detection

  setenv RRNA_DATA_DIR "$DATA_DIR/rrna"  # Directory containing data related to
										 # rRNAs detection

  setenv CDS_DATA_DIR "$DATA_DIR/cds"    # Directory containing data related to
										 # CDSs detection

  # setup the LC_ALL environment variable (for Linux mostly)
  # so the GNU tools (like sort) will work properly
  
  setenv LANG C
  setenv LC_ALL C

endif

# --------------------------------------
# path should be set each time
# --------------------------------------

set path = ($SCRIPT_DIR $BIN_DIR $path)

# --------------------------------------
# alias should be sourced each time
# --------------------------------------

alias Debug 'if ($Verbose) echo "# "\!:* >> /dev/stderr' 

alias Notify 'echo "# "\!:* >> /dev/stderr' 

alias Error 'echo "[0;31m# Error "\!:2-*"[m" >> /dev/stderr; Exit \!:1' 

alias Exit 'set Stat = \!:1; Debug "<--- $0 [$Stat]"; exit \!:1'

alias Stop 'Exit $Stat'

alias GetStatus 'set Stat = $status; Debug "status: $Stat"'

alias Exec 'Debug "execute: \!:*"; \!:* ; GetStatus'

alias OnError 'if ($Stat != 0)'

alias OnSuccess 'if ($Stat == 0)'

alias Usage 'egrep "^#  *usage:" $0 | sed -e "s/^# *//1"'

alias ExitUsage 'Usage; Exit 1'

alias Abort 'Error \!:1 "\!:2-* [abort]"'

alias CheckAbort 'GetStatus; OnError eval Abort \!:*'

alias NeedArg 'if ($#Argv < \!:1) eval ExitUsage'

alias Shift 'shift Argv'

alias NeedFile 'if (! -e \!:1) eval Error 5 "\!:1 : file not found"'

alias NeedDir 'if (! -d \!:1) eval Error 5 "\!:1 : directory not found"'

alias Cat 'awk '"'"'{print "# " $0}'"'"' \!:*'

alias AssignUndef 'if ($?\!:1 == 0) set \!:1=\!:2-*'

# --------------------------------------
# reset Stat each time
# --------------------------------------

set Stat = 0

# --------------------------------------
# process arguments to catch -v and -h
# and setup Argv, Argn (leave argv, argn untouched)
# Opts contains any standard option used (-v)
# --------------------------------------

set Argv = ()
set Opts = ()
@ Argn = 0
foreach arg ($argv) 
  switch ("$arg")
    case "-v":
      setenv Verbose 1
      set Opts = ($Opts $arg)
      breaksw
    case "-h":
      egrep "^#  *usage:" $0 | sed -e "s/^#  *//1" >> /dev/stderr
      kill $$  # exit will not work from sourced file
      breaksw
    default:
        set Argv = ($Argv "$arg")
        @ Argn ++
  endsw
end

# --------------------------------------
# notification for debug
# --------------------------------------

Debug "---> $0"
