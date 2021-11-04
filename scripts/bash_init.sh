#
# Bash file to be sourced at the begining of each bash script 
# for setting up basic variables and functions
#

export AwkCmd="gawk"

# setup the LC_ALL environment variable (for Linux mostly)
# so the GNU tools (like sort) will work properly
  
export LANG=C
export LC_ALL=C


########################
#
# General usage functions
#
#

function getAbsolutePath {
    [[ -d $1 ]] && { cd "$1"; echo "$(pwd -P)"; } || 
    { cd "$(dirname "$1")"; echo "$(pwd -P)/$(basename "$1")"; }
}

# Manage temp directory

function mktempdir {
	local TMP_DIR
	TMP_DIR=$(mktemp -d -t "$1_proc_$$_XXXXXX")
	logdebug "Creating temp directory $TMP_DIR"
	echo $TMP_DIR
}

function pushTmpDir {
	local TMP_DIR
	TMP_DIR=$(mktempdir "$1")
	pushd $TMP_DIR >& /dev/null
	TMP_DIR_STACK="$TMP_DIR $TMP_DIR_STACK"
	logdebug "Pushing temp directory $TMP_DIR"
	logdebug "Stack : ${TMP_DIR_STACK}"
}

function popTmpDir {
	local TMP_DIR
	TMP_DIR=$($AwkCmd '{print $1}' <<< $TMP_DIR_STACK)
	TMP_DIR_STACK=$($AwkCmd '{$1="";print $0}' <<< $TMP_DIR_STACK)
	popd  >& /dev/null
	rm -rf $TMP_DIR >& /dev/null
	logdebug "Poping temp directory $TMP_DIR"
	logdebug "Stack : ${TMP_DIR_STACK}"
}

# Logging functions
 

function errcho {
	>&2 echo $*
}

function openLogFile {
	ORG_LOGFILE=$1
	export ORG_LOGFILE
	touch ${ORG_LOGFILE}
}


function loginfo {
    errcho `date +'%Y-%m-%d %H:%M:%S'` "[OA INFO   ] $$ -- $1"
    if [[ ! -z ${ORG_LOGFILE} ]]; then
          echo `date +'%Y-%m-%d %H:%M:%S'` "[OA INFO   ] $$ -- $1" >> ${ORG_LOGFILE}
    fi
}

function logerror {
    errcho `date +'%Y-%m-%d %H:%M:%S'` "[OA ERROR  ] $$ -- $1"
    if [[ ! -z ${ORG_LOGFILE} ]]; then
          echo `date +'%Y-%m-%d %H:%M:%S'` "[OA ERROR   ] $$ -- $1" >> ${ORG_LOGFILE}
    fi
}

function logwarning {
    errcho `date +'%Y-%m-%d %H:%M:%S'` "[OA WARNING] $$ -- $1"
    if [[ ! -z ${ORG_LOGFILE} ]]; then
          echo `date +'%Y-%m-%d %H:%M:%S'` "[OA WARNING] $$ -- $1" >> ${ORG_LOGFILE}
    fi
}

function logdebug {
	if [[ ! -z ${ORG_DEBUG} ]]; then
	    errcho `date +'%Y-%m-%d %H:%M:%S'` "[OA DEBUG  ] $$ -- $1"
	    if [[ ! -z ${ORG_LOGFILE} ]]; then
	          echo `date +'%Y-%m-%d %H:%M:%S'` "[OA DEBUG  ] $$ -- $1" >> ${ORG_LOGFILE}
	    fi
	fi
}

# 
# Asserts that the number of arguments passed to the script
# is at least equal to the first argument of the function.
#
# needarg 3
#
# requires that the script is called at least with 3 arguments
#
__ORG_ANNOT_ARGS_SCRIPT_COUNT__=$#
function needarg {
	if (( $__ORG_ANNOT_ARGS_SCRIPT_COUNT__ < $1 )) ; then
	    logerror "not enougth arguments provided"
		exit 1
	fi
}

function needfile {
	if [[ ! -e $1 ]] ; then 
		logerror "File $1 doesn't exist"
		exit 1
	fi
}

function needdir {
	if [[ ! -d $1 ]] ; then 
		logerror "Directory $1 doesn't exist"
		exit 1
	fi
}
function assignundef {
    local value=$(eval echo \${$1+x})
    if [[ -z "$value" ]] ; then 
       eval $1=$2
    fi
}

# Sequence related functions	
	
# Counts how many sequences are stored in a fasta file
# 	- $1 : The fasta file to count
function fastaCount {
	grep '^>' $1 | wc -l
}


# compute the sequence length from a fasta sequence
# 	- $1 : The fasta file to cut
function seqlength {
  cat $1 | \
  wc |\
  $AwkCmd -v t="`head -1 $1 | wc -c`" '{print $3 - t - $1 + 1}'
}

function readfirstfastaseq {
	awk '(/^>/ && first) {on=0}
	     (on) {seq=seq $1}
	     (/^>/ && ! first) { first = 1
	                         on = 1
	                       } 
	     END {print seq}
	    ' $*
}

# extract a subseq from a fasta sequence
# 	- $1 : The fasta file to cut
#   - $2 : First position of the subsequence (first position is numered 1), 
#   - $3 : End of the subsequence (included in the subsequence)
function cutseq {
	$AwkCmd -v from=$2 -v end=$3 'function printfasta(seq) {            \
									seqlen=length(seq);             \
									for (i=1; i <= seqlen; i+=60)    \
									   print substr(seq,i,60);      \
								}                                   \
																	\
								/^>/   {print $0}                   \
								! /^>/ {seq=seq$0}                  \
								END {printfasta(substr(seq,from,end-from+1))}' $1
}

# Joins a set of sequences stored in a fasta file into 
# a single sequence
# 	- $1 : The fasta file containing the sequences to join
function joinfasta {
	$AwkCmd '(NR==1 && /^>/) {print $0}                                 \
	     ! /^>/          {print $0}' "${1}" |                           \
		 formatfasta
}

function fasta1line {
	$AwkCmd  '(/^>/ && seq !="") {print seq}          \
	          /^>/               {print $0;seq=""}    \
	          !/^>/              {seq=seq $0}         \
	          END                {print seq}' "${1}"
}

function formatfasta {
	$AwkCmd  'function printfasta(seq) {                            \
									seqlen=length(seq);             \
									for (i=1; i <= seqlen; i+=60)   \
									   print substr(seq,i,60);      \
								   }                                \
								(seq && /^>/) { printfasta(seq);    \
                                                seq=""}             \
								/^>/   { print $0 }                 \
								! /^>/ { seq=seq $0 }               \
								END    { printfasta(seq)}' "${1}"
}

# Reverse complement a DNA string 
#    - $1 : The DNA string to reverse complement
function reversecomp {
	echo $1 \
	| tr 'Aa' '@!' | tr 'Tt' 'Aa' | tr '@!' 'Tt' \
	| tr 'Cc' '@!' | tr 'Gg' 'Cc' | tr '@!' 'Gc' \
	| tr 'Mm' '@!' | tr 'Kk' 'Mm' | tr '@!' 'Kk' \
	| tr 'Rr' '@!' | tr 'Yy' 'Rr' | tr '@!' 'Yy' \
	| tr 'Ww' '@!' | tr 'Ss' 'Ww' | tr '@!' 'Ss' \
	| tr 'Bb' '@!' | tr 'Vv' 'Bb' | tr '@!' 'Vv' \
	| tr 'Dd' '@!' | tr 'Hh' 'Dd' | tr '@!' 'Hh' \
	| rev
}

#
# Process management related function
#

function timeoutcmd() {
  local seconde=$1
  shift

  $* &

  local mainpid=$!
  sleep $seconde &
  local sleeppid=$!

  local nproc=$(ps $mainpid $sleeppid | tail -n +2 | wc -l)

  while (( nproc > 1 )) ; do
    sleep 1
    nproc=$(ps $mainpid $sleeppid | tail -n +2 | wc -l)
  done

  local timealive=$(ps $sleeppid | tail -n +2 | wc -l)

  if (( timealive > 0 )) ; then
    kill -9 $sleeppid
  else
        if (( $(ps $mainpid | tail -n +2 | wc -l) > 0 )) ; then
                kill -9 $mainpid
        logwarning "Timeout after ${seconde}s on command : $*"
                return 1
        fi
  fi
}


#
#
########################



########################
#
# Local variable definitions
#
#

# The absolute path to the ORG.Annot home directory
ORG_HOME=`getAbsolutePath $(dirname ${BASH_SOURCE[0]})/..`

ORG_PORTNAME=`${ORG_HOME}/config/guess_port`	# The architecture running the ORG.Annot instance

BIN_DIR="${ORG_HOME}/ports/${ORG_PORTNAME}/bin" # Directory containing binaries for this port

SCRIPT_DIR="${ORG_HOME}/scripts"                # Directory containing scripts utilities

PROG_DIR="$(getAbsolutePath $(dirname $0))"		# Directory containing the main script file

LIB_DIR="$(getAbsolutePath ${PROG_DIR}/../lib)"	# Directory containing the main script libraries

CALL_DIR="$(getAbsolutePath $(pwd))"            # Directory from where the main script is called


DATA_DIR="${ORG_HOME}/data"                     # Directory containing reference data for the annotation

IR_DATA_DIR="${DATA_DIR}/ir"  				  	# Directory containing data related to
												# IRs detection
												
TRNA_DATA_DIR="${DATA_DIR}/trna"  				# Directory containing data related to
												# tRNAs detection
												
CDS_DATA_DIR="${DATA_DIR}/cds"  				# Directory containing data related to
												# CDSs detection
								
RRNA_DATA_DIR="${DATA_DIR}/rrna"  				# Directory containing data related to 
												# rRNAs detection
								
NUCRRNA_DATA_DIR="${DATA_DIR}/nucrrna"  				# Directory containing data related to 
												# rRNAs detection
								
ITS_DATA_DIR="${DATA_DIR}/its"  				# Directory containing data related to 
												# rRNAs detection
								

#
#
########################



########################
#
# Altering the environment
#
#

# We alter the path to include the bin dir corresponding to the port
PATH="${SCRIPT_DIR}:${BIN_DIR}:${PATH}"
export PATH

# Force to basic international setting for a correction behaviour of AWK on mac with float
export LANG=C
export LC_ALL=C

