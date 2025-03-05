FROM ubuntu:jammy AS builder
LABEL dockerfile.version="1"
LABEL software="phyloskims"

##
#
# Install org-annot
#
##

WORKDIR /
RUN apt-get update  --fix-missing && apt upgrade -y
RUN apt-get install -y build-essential 
RUN apt-get install -y git tcsh bash gawk parallel gettext zlib1g-dev libglib2.0-0
RUN git clone https://github.com/PhyloSkims/annotate.git org-annotate

#
# Compile tools
#

WORKDIR /org-annotate/src
RUN make


##
#
# Install Organelle Assembler
#
##

# WORKDIR "/" 
# RUN apt install -y python3-pip python3-dev python3-venv
# RUN git clone https://git.metabarcoding.org/org-asm/org-asm.git
# RUN python3 -m venv org-assembler 
# WORKDIR "/org-asm" 
# RUN bash -c '../org-assembler/bin/pip3 install -r requirements.txt'
# RUN bash -c '../org-assembler/bin/python3 setup.py install --no-serenity' 
# WORKDIR "/" 

##
#
# Install the final image grouping previously compiled softwares
#
##

FROM ubuntu:jammy AS phyloskims
WORKDIR /
RUN apt update  --fix-missing && apt upgrade -y
RUN apt-get install -y tcsh bash gawk parallel zlib1g libglib2.0-0 clustalo muscle cd-hit hmmer
COPY --from=builder /org-annotate /org-annotate
RUN rm -rf /org-annotate/src
# COPY --from=builder /org-assembler /org-assembler
# COPY --from=builder /org-assembler/bin/oa /bin
COPY organnot /bin
RUN  chmod +x /bin/organnot
#
# Recompile HMMER models
#
WORKDIR /org-annotate/data/its/ITSx_db/HMMs
RUN rm *.h3*
RUN for f in *.hmm ; do hmmpress $f ; done
WORKDIR /

RUN mkdir -p /data
WORKDIR /data