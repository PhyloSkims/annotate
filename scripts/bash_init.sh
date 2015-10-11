#
# Bash file to be sourced at the begining of each bash script 
# for setting up basic variables and functions
#

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

function pushTmpDir {
	TMP_DIR=$(mktemp -d -t "$1_proc_$$_")
	pushd $TMP_DIR >& /dev/null
	TMP_DIR_STACK="$TMP_DIR $TMP_DIR_STACK"
	logdebug "Pushing temp directory $TMP_DIR"
	logdebug "Stack : ${TMP_DIR_STACK}"
}

function popTmpDir {
	TMP_DIR=$(echo $TMP_DIR_STACK | awk '{print $1}')
	TMP_DIR_STACK=$(echo $TMP_DIR_STACK | awk '{$1="";print $0}')
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
  awk -v t="`head -1 $1 | wc -c`" '{print $3 - t - $1 + 1}'
}

# extract a subseq from a fasta sequence
# 	- $1 : The fasta file to cut
#   - $2 : First position of the subsequence (first position is numered 1), 
#   - $3 : End of the subsequence (included in the subsequence)
function cutseq {
	awk -v from=$2 -v end=$3 'function printfasta(seq) {            \
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
	awk '(NR==1 && /^>/) {print $0}                                 \
	     ! /^>/          {print $0}' $1 |                           \
		 formatfasta
}

function formatfasta {
	awk  'function printfasta(seq) {                                \
									seqlen=length(seq);             \
									for (i=1; i <= seqlen; i+=60)   \
									   print substr(seq,i,60);      \
								   }                                \
								/^>/   { print $0 }                 \
								! /^>/ { seq=seq $0 }               \
								END    { printfasta(seq)}' $1
}


#
#
########################



########################
#
# Local variable definitions
#
#

# The absolute path to the ORG.Annote home direcotory
ORG_HOME=`getAbsolutePath $(dirname ${BASH_SOURCE[0]})/..`

 
ORG_PORTNAME=`${ORG_HOME}/config/guess_port`	# The architecture running the ORG.Annnote instance

PROG_DIR="$(getAbsolutePath $(dirname $0))"		# Directory containing the main script file

DATA_DIR="${ORG_HOME}/data"                     # Directory containing reference data for the annotation

CALL_DIR="$(getAbsolutePath $(pwd))"            # Directory from where the main script is called


IR_DATA_DIR="${DATA_DIR}/ir"  				  	# Directory containing data related to the 
												# Inverted repeat strucuture
								
TRNA_DATA_DIR="${DATA_DIR}/trna"  				# Directory containing data related to the 
												# Inverted repeat strucuture
								

#
#
########################



########################
#
# Altering the environment
#
#

# We alter the path to include the bin dir corresponding to the port
PATH="${ORG_HOME}/ports/${ORG_PORTNAME}/bin:${PATH}"
export PATH

# Force to basic international setting for a correction behaviour of AWK on mac with float
export LANG=C
export LC_ALL=C


