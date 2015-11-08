#!/usr/bin/env gawk -f
function genomeid() {
  if (gid=="") {
    gid="XXXXXXX";
  }

  return gid;
}

function trnalib() {
  "echo $CAUTRNADB" | getline ref;
  return ref
}

function awkPID() {
  "bash -c 'echo $PPID'" | getline pid
  return pid
}

function tmpName(prefix,suffix) {
  return prefix awkPID() suffix
}

function rm(filename) {
  system("rm -f " filename);
}

function printfasta(id,seq,filename) {
  if (filename=="")
    filename = tmpName(id "_",".fasta");

	seqlen=length(seq);

  print ">" id > filename;
	for (i=1; i <= seqlen; i+=60)
			print substr(seq,i,60) >> filename;

  close(filename);
  return filename;

	}

function epissage(intron,seq) {
  if (intron != "") {
     l=length(intron);
     intron=substr(intron,2,l-2);
     split(intron,intronpos,",");
     lseq=length(seq);
     lintron2=lseq - intronpos[2] + 1;
     seq = substr(seq,1,intronpos[1]) substr(seq,intronpos[2],lintron2);
    }

  return seq;
}

function patchtRNA(anticodon,trna,seq) {
  if (anticodon == "cat") {
    file=printfasta(trna "_" anticodon,seq,"");

    command= "sumatra -t 0.9 -x -n " file " " trnalib() " 2>> /dev/null";
    while ((command | getline output) > 0) {
      split(output,field," ");
      match(field[2],"trn.M?");
      trna=substr(field[2],RSTART,RLENGTH);
      n[trna]+=field[5];
    }
    close(command)

    nmax=0;
    for (i in n) {
      dist=n[i];
      if (n[i] > nmax) {
        nmax=n[i];
        trna=i;
      }
    }

  }
  else {
    trna="trn" AA1[substr(trna,6,3)];
  }

  return trna;
}

function gene2product(gene) {
  return "tRNA-" AA3[substr(gene,4,1)];
}

function emblTRNA(geneid,trna,loc,anti,intron,seq) {
	if (loc ~ /^c/) {
	     sub("c\\[","complement(",loc); 
	     sub("\\]",")",loc);
	     sub(",","..",loc)}
	else {
		sub("\\[","",loc);
	    sub("\\]","",loc);
	    sub(",","..",loc)}
	    
	anti=toupper(anti);
    gsub("T","U",anti);
    product=gene2product(trna);
	
	if (intron!="") {
	   l=length(intron);
       intron=substr(intron,2,l-2);
       split(intron,intronpos,",");
	   ib=intronpos[0];
	   ie=intronpos[1];
	   match(loc,"[0-9][0-9]*");
	   gb=substr(loc,RSTART,RLENGTH);
	   sub("\\.\\.",".." (gb + ib -2) "," (gb + ie) "..",loc); \
	   sub("complement","complement(join",loc);\
	   if (substr(loc,1,1) ~ /[0-9]/) {
	      loc="join("loc} 
	      loc=loc")"; 
	      }
	      
	   print "FT   tRNA            " loc;
	   print "FT                   /gene=\""trna"\"";
	   print "FT                   /anticodon=\""anti"\"";	   
	   print "FT                   /product=\""product"("anti")\"";
	    
}

function gffTRNA(geneid,trna,loc,anti,intron,seq) {
  if (loc ~ /^c/) {
    complement="-";
    loc=substr(loc,3,length(loc)-3);
  }
  else {
    complement="+";
    loc=substr(loc,2,length(loc)-2);
  }

  split(loc,pos,",");
  anti=toupper(anti);
  gsub("T","U",anti);
  product=gene2product(trna);

  printf("%s\taragorn\tgene\t%d\t%d\t.\t%s\t.\tID=genetrn%d;gbkey=Gene;gene=%s\n",
         genomeid(),pos[1],pos[2],complement,geneid,trna);
  printf("%s\taragorn\ttRNA\t%d\t%d\t.\t%s\t.\tID=trn%d;parent=genetrn%d;gbkey=tRNA;gene=%s;Note=anticodon: %s;product=%s\n",
                genomeid(),pos[1],pos[2],complement,geneid,geneid,trna,anti,product);
  if (intron=="") {
    printf("%s\taragorn\texon\t%d\t%d\t.\t%s\t.\tID=exontrn%d;parent=trn%d;gbkey=tRNA;gene=%s;Note=anticodon: %s;product=%s\n",
                  genomeid(),pos[1],pos[2],complement,geneid,geneid,trna,anti,product);

  }
  else {
    l=length(intron);
    intron=substr(intron,2,l-2);
    split(intron,intronpos,",");
    lseq=pos[2]-pos[1]+1;
    lintron2=lseq - intronpos[2] + 1;
    seq = substr(seq,1,intronpos[1]) substr(seq,intronpos[2],lintron2);

    printf("%s\taragorn\texon\t%d\t%d\t.\t%s\t.\tID=exontrn%da;parent=trn%d;gbkey=tRNA;gene=%s;Note=anticodon: %s;product=%s\n",
                  genomeid(),pos[1],pos[1]+intronpos[1]-2,complement,geneid,geneid,trna,anti,product);

    printf("%s\taragorn\tintron\t%d\t%d\t.\t%s\t.\tID=exontrn%di;parent=trn%d;gbkey=intron;gene=%s;Note=anticodon: %s;product=%s\n",
                  genomeid(),pos[1]+intronpos[1]-1,pos[1]+intronpos[2]-2,complement,geneid,geneid,trna,anti,product);

    printf("%s\taragorn\texon\t%d\t%d\t.\t%s\t.\tID=exontrn%db;parent=trn%d;gbkey=tRNA;gene=%s;Note=anticodon: %s;product=%s\n",
                  genomeid(),pos[1]+intronpos[2]-1,pos[2],complement,geneid,geneid,trna,anti,product);

  }
}

BEGIN {
    AA1["Ala"]="A";
    AA1["Cys"]="C";
    AA1["Asp"]="D";
    AA1["Glu"]="E";
    AA1["Phe"]="F";
    AA1["Gly"]="G";
    AA1["His"]="H";
    AA1["Ile"]="I";
    AA1["Lys"]="K";
    AA1["Leu"]="L";
    AA1["Met"]="M";
    AA1["Asn"]="N";
    AA1["Pyl"]="O";
    AA1["Pro"]="P";
    AA1["Gln"]="Q";
    AA1["Arg"]="R";
    AA1["Ser"]="S";
    AA1["Thr"]="T";
    AA1["Sec"]="U";
    AA1["Val"]="V";
    AA1["Trp"]="W";
    AA1["Tyr"]="Y";

    AA3["A"]="Ala";
    AA3["C"]="Cys";
    AA3["D"]="Asp";
    AA3["E"]="Glu";
    AA3["F"]="Phe";
    AA3["G"]="Gly";
    AA3["H"]="His";
    AA3["I"]="Ile";
    AA3["K"]="Lys";
    AA3["L"]="Leu";
    AA3["M"]="Met";
    AA3["N"]="Asn";
    AA3["O"]="Pyl";
    AA3["P"]="Pro";
    AA3["Q"]="Gln";
    AA3["R"]="Arg";
    AA3["S"]="Ser";
    AA3["T"]="Thr";
    AA3["U"]="Sec";
    AA3["V"]="Val";
    AA3["W"]="Trp";
    AA3["Y"]="Tyr";
    AA3["f"]="fMet";

  }


/^>/ { id = substr($1,2);
     }

/^[0-9]+ +genes found/ \
     { nbgene=$1
     }

 ((geneid != "") && /^[0-9]+/ && ! /genes found/) \
     { seq=epissage(intron,seq);
       trna=patchtRNA(anti,trna,seq);
#       print geneid,trna,loc,anti,"'"intron"'",seq;
       emblTRNA(geneid,trna,loc,anti,intron,seq);
       seq=""
     }


/^[0-9]+/ && ! /genes found/ \
     { geneid=$1;
       trna  =$2;
       loc   =$3;
       lseq  =$4;
       x=$5;
       split($5,intron_desc,"i");
       anti  =substr(intron_desc[1],2,3);
       intron=intron_desc[2];
     }

/^[^>0-9]/ \
     { seq=seq $1

     }

END  { seq=epissage(intron,seq);
       trna=patchtRNA(anti,trna,seq);
#       print geneid,trna,loc,anti,"'"intron"'",seq;
       emblTRNA(geneid,trna,loc,anti,intron,seq);
     }
