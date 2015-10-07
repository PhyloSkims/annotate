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
}

function popTmpDir {
	popd  >& /dev/null
	rm -rf $TMP_DIR >& /dev/null
}

# Logging functions
 

function errcho {
	>&2 echo $*
}

function openLogFile {
	LOGFILE=$1
	touch ${LOGFILE}
}


function loginfo {
    errcho `date +'%Y-%m-%d %H:%M:%S'` "[OA INFO   ] $$ -- $1"
    if [[ ! -z ${LOGFILE} ]]; then
          echo `date +'%Y-%m-%d %H:%M:%S'` "[OA INFO   ] $$ -- $1" >> ${LOGFILE}
    fi
}

function logerror {
    errcho `date +'%Y-%m-%d %H:%M:%S'` "[OA ERROR  ] $$ -- $1"
    if [[ ! -z ${LOGFILE} ]]; then
          echo `date +'%Y-%m-%d %H:%M:%S'` "[OA ERROR   ] $$ -- $1" >> ${LOGFILE}
    fi
}

function logwarning {
    errcho `date +'%Y-%m-%d %H:%M:%S'` "[OA WARNING] $$ -- $1"
    if [[ ! -z ${LOGFILE} ]]; then
          echo `date +'%Y-%m-%d %H:%M:%S'` "[OA WARNING] $$ -- $1" >> ${LOGFILE}
    fi
}

# Sequence related functions	
	
function fastaCount {
	grep '^>' $1 | wc -l
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

