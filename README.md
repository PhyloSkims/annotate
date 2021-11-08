# ORG.Annot - Organelle Annotator

**version 1.0.0 - Sept 2015**

## -1- Installation

# First check if binaries have been already compiled for your port

```bash
scripts/check_port.sh
```

if it prints :

```
checking port compilation OK
```

then everything is fine

if it prints :

```
port not yet compiled
```

then, you should first compile binaries for your port
please consult: src/README.txt for details

if it prints something like

```
XXX version A.B.C (should be >= X.Y.Z)
please consider installing XXX from src/_unix_tools_
```

then, you should first compile some unix tools binaries for your port
please consult: `src/README.txt` for details

## -2- Distribution organisation


- config    :  config file for (re)compiling tools
- data      :  internal data files (e.g. ChloroDB)
- detectors :  features specific detector scripts
- ports     :  binaries for various ports
- scripts   :  main scripts directory
- src       :  sources for (re)compiling tools

## -3- Usage


> note: all scripts are located in scripts, so you may add
>       this directory in your path
>       e.g. using csh 
>       set path = ($path <ROOTDIR>/scripts)

@TODO : scripts description

test
