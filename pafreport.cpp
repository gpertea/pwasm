#include <stdio.h>
#include <ctype.h>
#include <math.h>
#include "GArgs.h"
#include "GStr.h"
#include "GHash.hh"
#include "GList.hh"
#include "GFaSeqGet.h"
#include "GapAssem.h"
#include "codons.h"

#define USAGE "Usage:\n\
 pafreport <paf_with_cg_cs> -r <refseq.fa> [-s <summary.txt>]\n\
    [-o <diff_report.dfa>][-w <outfile.mfa>] [-G|-F|-C|-N]\n\
    \n\
   <paf_with_cg_cs> is the input PAF file with high quality query sequence(s)\n\
      aligned to many target sequences using minimap2 --cs\n\
   -r provide the fasta file with query sequence(s) (required)\n\
   -o write difference data for each alignment into <diff_report.dfa>\n\
   -s write event summary counts into <summary.txt>\n\
   -w write MSA as multifasta into <outfile.mfa>\n\
   -G gene CDS analysis mode (default for query<100K; assumes -C)\n\
   -F full genome alignment mode (default for query>100Kb; assumes -N)\n\
   -R reverse mapping mode: input PAF has mappings of reads/contigs to\n\
      a single reference sequence\n\
   -A like -R mode, but assemble all read mappings into a consensus and\n\
      report differences on the assembly instead of individual mappings\n\
   -C perform codon impact analysis\n\
   -N skip codon impact analysis\n"

#define LOG_MSG_CLIPMAX "Overlap between %s and target %s rejected due to clipmax=%4.2f constraint.\n"
#define LOG_MSG_OVLCLIP "Overlap between %s and %s invalidated by the %dnt clipping of %s at %d' end.\n"
//-------- global variables
bool debugMode=false;
bool removeConsGaps=false;
//PAF report options:
bool fullgenomeAlns=false; //full genome vs genome alignments, all query-target alignments are retained
                     //otherwise only first query-target alignment is processed
bool skipCodAn=false; //skip codon impact assessments (-N/-C)?

bool verbose=false;
bool revMapping=false; // -R : reverse mappings of multiple reads/contigs vs a reference sequence
bool makeConsensus=false; // -A : assemble reverse mappings into a consensus before analyzing differences

//methylation motifs:
//TODO: these should be loaded from an external file
const char* const metmot[] = { "CCTGG", "CCAGG", "GATC", "GTAC", NULL };
//look for these or a homopolymer around the indel/substitution event

void printDebugAln(FILE* f, GSeqAlign* a);

float clipmax=0;
//--------------------------------
struct GapData {
	int pos; //position of gap on sequence
	int len; //length of gap (0 means 1)
	GapData(int p=0, int l=1):pos(p),len(l) { }
};

struct AlnInfo {
  char reverse;
  const char* r_id;
  int r_len;
  int r_alnstart;
  int r_alnend;
  const char* t_id;
  int t_len; //raw length directly from the PAF field
  int t_alnstart;
  int t_alnend;
  AlnInfo(char strand=0):reverse(2), r_id(NULL), r_len(0), r_alnstart(0), r_alnend(0),
		  t_id(NULL), t_len(0), t_alnstart(0), t_alnend(0) {
		if (strand=='+') reverse=0;
		 else if (strand=='-') reverse=1;
  }
  AlnInfo(char strand, const char* rid, const char* rlen, const char* rstart, const char* rend,
		  const char* tid, const char* tlen, const char* tstart, const char* tend):r_id(rid), t_id(tid) {
	reverse = (strand=='-') ? 1 : 0;
	r_len = atoi(rlen);
	r_alnstart=atoi(rstart);
	r_alnend=atoi(rend);
	t_len = atoi(tlen);
	t_alnstart=atoi(tstart);
	t_alnend=atoi(tend);
  }
  void init(const char* rid, const char* rlen, const char* rstart, const char* rend,
		    const char* tid, const char* tlen, const char* tstart, const char* tend,
		    char strand=0) {
	r_id=rid;
	t_id=tid;
	r_len = atoi(rlen);
	r_alnstart=atoi(rstart);
	r_alnend=atoi(rend);
	t_len = atoi(tlen);
	t_alnstart=atoi(tstart);
	t_alnend=atoi(tend);
	if (strand) {
		if (strand=='+') reverse=0;
		 else if (strand=='-') reverse=1;
	}
  }
};

struct TDiffInfo {
	char evt; //event code: I=insertion, D=deletion, S=substitution
	int evtlen; //event length (bases; 0 for deletion)
	GStr evtbases; //bases inserted, deleted or newly substituted (by these bases)
	GStr evtsub;  //for substitutions only: the original base(s)
	//potential diff cause: homopolymer or methylation motif
	int dcoffset; //rloc - (motif/hpoly start location)
	char hpoly; //if homopolymer found, the repeated base
	int hpolylen; //length of homopolymer, if found
	GStr motif; //motif found

	int rloc; //location of the event on theref query sequence
	//---
	//TODO: tloc must be converted to GLOBAL coordinates, on the full target sequence
	int tloc; //location of the event within the aligned target region, on the *aligned strand*
	GStr tctx; //target context: event+5 bases before and after the event (at least 10 bases)
	char flags; //indicates if a methylation motif or a homopolymer was found in the context
	TDiffInfo(char e=0):evt(e), evtlen(0), evtbases("",8), evtsub("",8),
			dcoffset(0),hpoly(0),hpolylen(0),motif(),
			rloc(0),tloc(0), tctx("",18), flags(0) { }
	void init(char e, int len, int rpos, int tpos) {
		evt=e;
		evtlen=len;
		rloc=rpos;
		tloc=tpos;
		evtbases.clear(8);
		evtsub.clear(8);
		tctx.clear(18);
		flags=0;
	}
	void setTContext(GStr& tseq) {
		 int evt_len=evtlen;
		 if (evt=='D') evt_len=0;
		 int thi=(evt_len+10)>>1;
		 int tc_start=tloc-thi;
		 if (tc_start<0) tc_start=0;
		 int tc_end=tloc+thi;
		 if (tc_end>=tseq.length()) tc_end=tseq.length()-1;
		 tctx=tseq.substr(tc_start, tc_end-tc_start);
	}
    bool operator<(TDiffInfo o) {
    	return (rloc<o.rloc);
    }
};

class PAFAlignment {
 public:
  AlnInfo alninfo;
  GVec<GapData> rgaps; //reference query gaps
  GVec<GapData> tgaps; //target gaps
  GVec<TDiffInfo> tdiffs; //indel/substitution events on this target sequence
  //char* seqname; //aligned target sequence (read, or genome mapping)
  char* cs; //cs tag value = difference string (short)
  int edist; //number of 1bp edits (edit distance) = NM tag value
  int alnscore; // AS tag value
  char* cigar;

  int seqlen; //length of this tseq (target sequence, unedited)
  int offset; //offset of this mapping relative to reference
  int clip5; //amount to clip on the left end (0 for PAF on tseq)
  int clip3; //amount to clip on the right end (0 for PAF on tseq)
  char reverse; //0, or 1 if this mapping is reverse complemented
  GStr rlabel;
  GStr tlabel;
  void parseErr(int fldno, const char* line);
  void printDiffInfo(FILE* f, GASeq& refseq);
  PAFAlignment(GDynArray<char*>& t, AlnInfo& alni, GASeq& refseq, GStr& tseq, const char* line);
   //this also rebuilds tseq with the target sequence, by transforming refseq according to cs string
  ~PAFAlignment() { GFREE(cs); GFREE(cigar); }
  //int nextRefGap(int& pos);
  //int nextSeqGap(int& pos);
 };

//returns the end of a space delimited token
void skipSp(char*& p) {
 while (*p==' ' || *p=='\t' || *p=='\n') p++;
}
char* endSpToken(char* str) {
 if (str==NULL || *str==0) return NULL;
 char* p=str;
 for (;*p!=0;p++)
  if (*p==' ' || *p=='\t' || *p=='\n') return p;
 return p; // *p is '\0'
}

//========================================================
//====================     main      =====================
//========================================================
int main(int argc, char * const argv[]) {
 //GArgs args(argc, argv, "DGvd:o:c:");
 GArgs args(argc, argv, "DGFRACNvd:p:r:o:m:w:c:s:");
 int e;
 if ((e=args.isError())>0) {
    GMessage("%s\nInvalid argument: %s\n", USAGE, argv[e]);
    exit(1);
 }
 if (args.getOpt('h')!=NULL) GError("%s\n", USAGE);
 debugMode=(args.getOpt('D')!=NULL);
 //removeConsGaps=(args.getOpt('G')==NULL);
 fullgenomeAlns=(args.getOpt('F')!=NULL);
 bool geneCDSalns=(args.getOpt('G')!=NULL);
 if (fullgenomeAlns && geneCDSalns) {
	 GMessage("%s Error: cannot use both -G and -F!\n",USAGE);
	 exit(1);
 }
 revMapping=(args.getOpt('R')!=NULL);
 makeConsensus=(args.getOpt('A')!=NULL);
 if (makeConsensus) revMapping=true;
 bool forceCoding=(args.getOpt('C')!=NULL);
 bool forceNonCoding=(args.getOpt('N')!=NULL);
 if (forceCoding && forceNonCoding) {
	 GMessage("%s Error: cannot use both -N and -C!\n",USAGE);
	 exit(1);
 }
 verbose=(args.getOpt('v')!=NULL);
 if (debugMode) {
	  verbose=true;
	  //outf=stdout;
 }
 MSAColumns::removeConsGaps=removeConsGaps;
 GStr infile;
 if (args.startNonOpt()) {
   infile=args.nextNonOpt();
 }
 //==
 FILE* inf;
 if (!infile.is_empty()) {
    inf=fopen(infile, "r");
    if (inf==NULL)
       GError("Cannot open input file %s!\n",infile.chars());
    }
  else
   inf=stdin;
  GStr s=args.getOpt('c');
  if (!s.is_empty()) {
      bool ispercent=(s[-1]=='%');
      if (ispercent) s.trimR('%');
      int c=s.asInt();
      if (c<=0) GError("Error: invalid -c <clipmax> (%d) option provided (must be "
                             "a positive integer)!\n",c);
      if (ispercent && c>99)
                GError("Error: invalid percent value (%d) for -c option "
                       " (must be an integer between 1 and 99)!\n",c);
      if (ispercent) {
            clipmax = float(c)/100;
            int maxp=iround(100*clipmax);
            if (verbose) fprintf(stderr,
                  "Percentual max clipping set to %d%%\n",
                                               maxp);
            }
         else {
          clipmax=c;
          if (verbose) fprintf(stderr, "Max clipping set to %d bases\n",
                                               c);
          }

      } //clipmax option

  GStr outfile=args.getOpt('o');
  FILE* freport=stdout;
  if (!outfile.is_empty()) {
     freport=fopen(outfile, "w");
     if (freport==NULL)
        GError("Cannot open file %s for writing!\n",outfile.chars());
  }
  //else freport=stdout;
  //************************************************
  s=args.getOpt('r');
  if (s.is_empty()) GError("Error: query FASTA file (-r) is required!\n");
  int64 fsize=fileSize(s.chars());
  if (fsize<=0) GError("Error: invalid FASTA file %s !\n", s.chars());
  GFastaDb qfasta(s.chars(), false);
  if (!fullgenomeAlns && !geneCDSalns && fsize>120000) {
		  fullgenomeAlns=true;
  }
  skipCodAn=(fullgenomeAlns || forceNonCoding);
  if (!skipCodAn && !forceCoding && fsize>120000) {
	  skipCodAn=true;
  }
  GStr msafile=args.getOpt('w');
  FILE* fmsa=NULL;
  if (!msafile.is_empty()) {
	 if (fullgenomeAlns) {
		 GMessage("%s Error: can only generate MSA for -G mode!\n",USAGE);
		 exit(1);
	 }
     fmsa=fopen(msafile, "w");
     if (fmsa==NULL)
        GError("Cannot open file %s for writing!\n",msafile.chars());
  }
  s=args.getOpt('s');
  GHash<int> alnpairs(true); //key is q_name+'~'+t_name ; only used if !fullgenomeAlns

  GHash<GASeq> rseqs(false); //query seqs to use as reference for MSAs
  GList<GSeqAlign> msalns(false,   true, false);
                      // unsorted, free, not unique
  //msalns.setSorted(compareOrdnum);
  GLineReader* linebuf=new GLineReader(inf);
  char* line;
  GSeqAlign *refMSA=NULL; //this will be the final MSA
  int numalns=0;
  GDynArray<char*> t(24);
  GASeq *refseq=NULL;
  GASeq *refseq_rc=NULL;
  while ((line=linebuf->getLine())!=NULL) {
   if (line[0]=='#') continue;
   bool firstRefAln=false;
   GStr lstr(line);
   int numt=strsplit(line, t, '\t');
   if (numt<15)
	   GError("Error: invalid PAF fline (num. fields=%d):\n%s\n", numt,lstr.chars());
   //if (strcmp(t[0], refseq->getId())!=0)
   //   GError("Error: expected reference name (%s) in PAF line, but found %s instead!\n",
   //    		  refseq->getId(),t[0]);
   AlnInfo al(*t[4]);
   if (revMapping)
	    al.init(t[5], t[6], t[7], t[8],
	    		t[0], t[1], t[2], t[3]);
   else al.init(t[0], t[1], t[2], t[3],
	            t[5], t[6], t[7], t[8]);
   if (strcmp(al.r_id, al.t_id)==0) {
	  if (verbose) GMessage("Skipping alignment of qry seq to itself.\n");
      continue; //skip redundant inclusion of reference as its own child
   }
   GStr qtpair(al.r_id);
   if (!fullgenomeAlns) { //gene CDS mode
	   qtpair.append('~');
	   qtpair.append(al.t_id);
	   int* vcount=alnpairs.Find(qtpair.chars());
	   if (vcount==NULL) {
		   alnpairs.Add(qtpair.chars(), new int(0));
	   }
	   else {
         ++(*vcount);
	     if (*vcount==1)
	       	GMessage("Warning: alignment %s to %s already seen, ignoring \n",
	                            al.r_id, al.t_id);
	     continue;
	   }
   }
   //-- load the current alignment
   ++numalns;
   //retrieve current refseq
   if (refseq==NULL || strcmp(refseq->name(), al.r_id)!=0) {
      delete refseq_rc;
      //check if this refseq was seen before
      GASeq* vrseq=rseqs.Find(al.r_id);
      if (vrseq!=NULL) {
        refseq=vrseq;
      }
      else {
        GFaSeqGet* seq=qfasta.fetch(al.r_id);
        if (seq==NULL) GError("Error: could not retrieve sequence for %s !\n",al.r_id);
        refseq=new GASeq(al.r_id, NULL, seq->seq(), seq->getseqlen());
        refseq->setFlag(GA_flag_IS_REF);
        refseq->allupper();
      }
      if (!revMapping) {
        refseq_rc=new GASeq(*refseq);
        refseq_rc->reverseComplement();
      }
   }

   if (al.r_len!=refseq->getSeqLen())
	    GError("Error: ref seq len in this PAF line (%d) differs from loaded sequence length(%d)!\n%s\n",
		     al.r_len, refseq->getSeqLen(),lstr.chars());

   GASeq *r_seq = refseq;
   if (al.reverse && !revMapping)
      r_seq=refseq_rc;
   GStr tseq("", al.t_alnend-al.t_alnstart+2);
   PAFAlignment* aln=new PAFAlignment(t, al, *r_seq, tseq, lstr.chars());
   //this also fills tseq and sets aln->offset to the tseq mapping offset on refseq
   //GStr tlabel(al.t_id, 26);
   //GStr rlabel(al.r_id, 26);
   if (fullgenomeAlns) {
	   aln->rlabel+=':';
	   aln->rlabel+=(al.r_alnstart+1);
	   aln->rlabel+='-';aln->rlabel+=al.r_alnend;
   }
   aln->tlabel+=":";aln->tlabel+=(al.t_alnstart+1);
   aln->tlabel+='-';aln->tlabel+=al.t_alnend;
   if (al.reverse) aln->tlabel+='-';
   else aln->tlabel+='+';
   //if (freport) {
   if (qfasta.faIdx->getCount()==1 && !fullgenomeAlns)
   	aln->rlabel="";
   if (!makeConsensus) {
        aln->printDiffInfo(freport, *refseq);
   }
   GASeq* taseq=new GASeq(aln->tlabel.chars(), "", tseq.chars(), tseq.length(), al.r_alnstart);
   if (!revMapping) taseq->revcompl=aln->reverse;
   GASeq* rseq=refseq; //only for the first ref alignment
   if (refseq->msa!=NULL) { //already have a MSA
     //need to merge this alignment into existing MSA refseq->msa
	//rseq is just an instance of refseq in this alignment
	 rseq=new GASeq(refseq->getId(), 0, al.r_len);
			 // al.r_alnstart, al.r_len-al.r_alnend, 0);
   }
   else {
	   firstRefAln=true;
   }
   if (fmsa!=NULL) {
	   //once a gap, always a gap
	   //propagate gaps in ref seq from the current alignment
	   for (int g=0;g<aln->rgaps.Count();++g) {
		   GapData &gd=aln->rgaps[g];
		   rseq->setGap(gd.pos, gd.len);
	   }
	   //propagate gaps in target sequence as well
	   for (int g=0;g<aln->tgaps.Count();++g) {
		   GapData &gd=aln->tgaps[g];
		   taseq->setGap(gd.pos, gd.len);
	   }
	   GSeqAlign *newmsa=new GSeqAlign(rseq, taseq);
	   //this also sets rseq->msa and taseq->msa to newmsa
	   if (firstRefAln) { //first alignment with refseq so rseq==refseq
			 //newmsa->incOrd();
			 newmsa->ordnum=numalns;
			 msalns.Add(newmsa);
			 refMSA=newmsa;
	   }
	   else { // rseq != refseq, but an instance of refseq in this current aln
			//incrementally add this alignment, as a MSA, to existing MSA,
			//propagating gaps through rseq instance
		   refseq->msa->addAlign(refseq, newmsa, rseq);
		   refMSA=refseq->msa;
		   delete newmsa; //incorporated, no longer needed
	   }
   } //MSA constructions
//   if (!fullgenomeAlns) {
//	   rseqs.Add(al.r_id, )
//   }
   // new pairwise alignment added to MSA
   // debug print the progressive alignment
   /*
   if (debugMode) {
    //for (int a=0;a<alns.Count();a++) {
    	printDebugAln(stderr, refMSA);
     // }
    }
    */
   delete aln; //no longer needed
   if (linebuf->isEof()) break;
 } //----> PAF line parsing loop

  //*************** now build MSAs and write them
  delete linebuf;
  //print all MSAs found
  //for (int i=0;i<alns.Count();i++) {
   //GSeqAlign* a=alns.Get(i);
   if (debugMode) { //write plain text alignment file
     fprintf(stderr,">MSA (%d)\n",refMSA->Count());
     refMSA->print(stderr, 'v');
   }
   if (fmsa) {
     refMSA->writeMSA(fmsa);
     fclose(fmsa);
   }
 //} // for each MSA cluster
 // --------- D O N E --------
  msalns.Clear();
  delete refseq_rc;
  //delete refMSA;
  rseqs.Clear();
  //fflush(outf);
  if (freport!=stdout) fclose(freport);
  if (inf!=stdin) fclose(inf);
}

//---------------------- RefAlign class
void PAFAlignment::parseErr(int fldno, const char* line) {
 fprintf(stderr, "Error parsing input line (field %d):\n%s\n",
         fldno, line);
 exit(3);
}

void revCompl(GStr& s) {
	s.tr(IUPAC_DEFS, IUPAC_COMP);
	s.reverse();
}
const char* CIGAR_ERROR="Error parsing cigar string from line: %s (cigar position: %s)\n";
const char* CS_ERROR="Error parsing cs string from line: %s (cs position: %s)\n";

//this also rebuilds tseq with the target sequence, by transforming refseq according to cs string
PAFAlignment::PAFAlignment(GDynArray<char*>& t, AlnInfo& al, GASeq& refseq, GStr& tseq,
		const char* line):rgaps(), tgaps(), tdiffs(), cs(NULL), edist(-1), alnscore(0),
		cigar(NULL), seqlen(0), offset(0), clip5(0), clip3(0), reverse(0),
		rlabel(), tlabel() {
  //toffs must be the offset of the alignment start from the beginning of the full, original refseq
  // so it's r_clip5 (or r_clip3 if it's revcomp)
  //seqname=Gstrdup(t[5]);
  alninfo=al; //do we still need to keep it? *_id pointers will be obsolete
  rlabel=al.r_id;
  tlabel=al.t_id;
  //if (!revMapping) {
	reverse=al.reverse;
  //}
  //clip5=0; //always cutting out the aligned region from the target sequence
  //clip3=0;
  offset=al.r_alnstart;
  //int base_ofs=offset; //mismatch base offset -- different for
  if (reverse && !revMapping) { //offset is on the reverse complement ref string
	   offset=al.r_len-al.r_alnend;
  }
  seqlen=al.t_alnend-al.t_alnstart; //check this against the cigar string, and the cs string
  int tnum=t.Count();
  char got=0;
  char gotall=1+2+4+8;
  for (int i=12;i<tnum;++i) {
	  if (startsWith(t[i], "NM:i:")) {
		  edist=atoi(t[i]+5);
		  got|=1;
		  if (got==gotall) break;
		  continue;
	  }
	  if (startsWith(t[i], "AS:i:")) {
		  alnscore=atoi(t[i]+5);
		  got|=2;
		  if (got==gotall) break;
		  continue;
	  }
	  if (startsWith(t[i], "cg:Z:")) {
          cigar=Gstrdup(t[i]+5);
          got|=4;
          if (got==gotall) break;
          continue;
	  }
	  if (startsWith(t[i], "cs:Z:")) {
          cs=Gstrdup(t[i]+5);
          got|=8;
          if (got==gotall) break;
          continue;
	  }
  }
  if (cigar==NULL || cigar[0]==0) GError(CIGAR_ERROR, line, 0);
  int qpos = 0;  //should match length of aligned region of reference qseq, minus offset
  int tpos = 0; //will end up with the length of target genome region being aligned
  int eff_t_len=al.t_alnend-al.t_alnstart; //effective target sequence length
  //target sequence is to be trimmed to the aligned region
  char *p=cs;
  TDiffInfo dif;
  //interpret cs string to fill tseq
  tseq="";
  qpos=0; //ref query location in the alignment
  tpos=0;
  char tch=0,qch=0;
  int q_pos=0; //real ref query location (adjusted for offset and strand)
  int s_pos=0; //similar for qry
  int e_len=0;
  while (*p != 0) {
    char op=*p;
    int cl=0;
    p++;
    if (revMapping) {
    	if (op=='-') op='+';
    	else if (op=='+') op='-';
    }
    switch (op) {
      case ':':
    	if (!parseInt(p,cl)) GError(CS_ERROR, line, p);
    	tseq.appendmem(refseq.getSeq()+offset+qpos, cl);
    	qpos+=cl;
    	tpos+=cl;
    	break;
      case '*': //substitution
    	if (revMapping) {
   	       qch=toupper(*p);++p;
  	       tch=toupper(*p);++p;
    	} else {
	       tch=toupper(*p);++p;
	       qch=toupper(*p);++p;
    	}
	    q_pos=offset+qpos;
	    if (qch!=refseq.getSeq()[q_pos]) {
	    	GError("Error: base mismatch %c != qstr[%d] (%c) at line\n%s\n",
	    			qch, q_pos, refseq.getSeq()[q_pos], line);
	    }
	    //merge adjacent substitutions into a single event
	    if (tdiffs.Count()>0 && tdiffs.Last().evt=='S' &&
	    		tdiffs.Last().rloc==q_pos-tdiffs.Last().evtbases.length()) {
	    	tdiffs.Last().evtbases.append((char)toupper(tch));
	    	tdiffs.Last().evtsub.append((char)toupper(qch));
	    }
	    else {
			dif.init('S', 1, q_pos, tpos);
			dif.evtbases.append((char)toupper(tch));
			dif.evtsub.append((char)toupper(qch));
			//keep substitution location on reverse to simplify merging
			/*
			if (reverse) {
				revCompl(dif.evtbases);
				dif.rloc=al.r_len-q_pos;
			}
			*/
			tdiffs.Add(dif);
	    }
	    //
	    tseq.append((char)tolower(tch));
	    ++qpos;
	    ++tpos;
   	    break;
      case '-': //gap in ref qseq (insertion in tseq)
    	//this gap is at position (offset+qpos) in ref query if reverse==0
    	// BUT at rlen-(offset+qpos) if reverse==1
    	//s_pos= reverse ? eff_t_len-tpos : tpos;
    	s_pos=tpos;
    	//these bases are only in tseq (insert)
    	while (isalpha(tch=toupper(*p))) {
    		++p;
    		++tpos;
    		tseq.append((char)tolower(tch));
    	}
		//insert in tseq
    	e_len=tpos-s_pos;
    	q_pos=offset+qpos;
		dif.init('I', e_len, q_pos, s_pos);
		dif.evtbases.append(tseq.substr(-e_len));
		if (reverse && !revMapping) {
			revCompl(dif.evtbases);
			dif.rloc=al.r_len-q_pos;
		}
		this->tdiffs.Add(dif);
    	//qpos unchanged, offset+qpos is the location of the gap in qseq
    	break;
      case '+': //gap in tseq
    	s_pos=qpos;
    	while (isalpha(qch=toupper(*p))) {
    		++p;
    		++qpos;
    	}
    	e_len=qpos-s_pos;
    	q_pos=s_pos+offset;
    	//these bases are missing in tseq (deletion)
		dif.init('D', e_len, q_pos, tpos);
		dif.evtbases.append(refseq.getSeq()+q_pos, e_len);
		if (reverse && !revMapping) {
			revCompl(dif.evtbases);
			dif.rloc=al.r_len-q_pos-e_len;
		}
		this->tdiffs.Add(dif);
    	//tpos unchanged, it is the location of the gap in tseq
    	break;
      case '~':
    	GError("Error: spliced alignments not supported! at line:\n%s\n", line);
    	break;
      default:
        GError("Error: unhandled event at %s in cs, line:\n%s\n",p,line);
    }
  } //interpret CS string;
  //fill in context for differences:
  for (int d=0;d<tdiffs.Count();d++) {
	  tdiffs[d].setTContext(tseq);
	  if (reverse) {
		if (revMapping) {
			if (tdiffs[d].evt=='S' || tdiffs[d].evt=='I')
				 tdiffs[d].tloc+=tdiffs[d].evtbases.length();
		    tdiffs[d].tloc=tseq.length()-tdiffs[d].tloc;
		//this location is relative to the beginning of the alignment
		} else { // !revMapping
			revCompl(tdiffs[d].tctx);
			//also reverse the location so it shows as if tseq was reversed
			if (tdiffs[d].evt=='S') {
			  //substitutions were kept on reverse to simplify merging
			  //so now it's time to adjust that
				revCompl(tdiffs[d].evtbases);
				revCompl(tdiffs[d].evtsub);
				tdiffs[d].rloc=al.r_len-tdiffs[d].rloc-tdiffs[d].evtbases.length();
			}
		} //reverse && !revMapping
	  } //if reverse
  } //for each tdiff
  if (reverse & !revMapping) tdiffs.Reverse();
  //parse cigar string to get the gaps
  p=cigar;
  int mbases = 0; //count "aligned" bases (includes mismatches)
  qpos=0;
  tpos=0;
  GapData gap;
  // gpos = current genomic position (will end up as right coordinate on the genome)
  // rpos = read position (will end up as the length of the read)
  // cop = CIGAR operation, cl = operation length
  while (*p != '\0') {
	  int cl=0;
	  if (!parseInt(p,cl)) GError(CIGAR_ERROR, line, p);
	  char cop=*p;
	  if (cop=='\0') GError(CIGAR_ERROR, line, p);
	  if (revMapping) {
		  if (cop=='D') cop='I';
		  else if (cop=='I') cop='D';
	  }
	  switch (cop) {
	   case 'X':
	   case 'M':
	   case '=':
		 //printf("[%d-%d]", gpos, gpos + cl - 1);
		 tpos+=cl;
		 qpos+=cl;
		 ++mbases;
		 break;
	   case 'P': //padding;silent deletion from padded reference
		 // printf("[%d-%d]", pos, pos + cl - 1); // Spans positions, No Coverage
		 // neither gpos nor rpos are advanced by this operation
		 break;
	   case 'H':
		 // printf("[%d]", pos);  // No coverage
		 // neither gpos nor rpos are advanced by this operation
		 break;
	   case 'S':
		 //soft clipped bases on query
		 qpos+=cl;
		 GMessage("Warning: soft clipping shouldn't be found in this application!\n%s\n", line);
		 break;
	   case 'I':
		 //(gap in genomic (target) seq)
		 // tpos is not advanced by this operation
		 gap.pos= reverse? eff_t_len-tpos  : tpos;
		 gap.len=cl;
		 tgaps.Add(gap);
		 qpos+=cl;
		 break;
	   case 'D':
		 //deletion in target sequence relative to the query (gap in query ref seq)
		 gap.pos=offset+qpos;
		 if (reverse && !revMapping) //actual location on qry ref for a reverse complement match:
			  gap.pos=al.r_len-gap.pos;
		 gap.len=cl;
		 rgaps.Add(gap);
		 //insertion in tseq
		 tpos += cl;
		 break;
	   case 'N':
		 // intron
		 //special skip operation, not contributing to "edit distance",
		 //   so num_mismatches should not update
		 tpos+=cl;
		 //shouldn't really happen in genomic PAF
		 gap.pos=offset+qpos;
		 if (reverse && !revMapping) //actual location on qry ref for a reverse complement match:
			  gap.pos=al.r_len-gap.pos;
		 gap.len=cl;
		 rgaps.Add(gap);
		 GMessage("Warning: introns are not supported by this application!\n%s\n", line);
		 break;
	   default:
		 GError("Error: unhandled cigar_op %c (len %d) in %s\n", cop, cl, line);
	  } //switch cigar op
	++p;
  } // interpret_CIGAR string
  if (eff_t_len!=tpos || tseq.length()!=tpos)
   	GError("Error: tseq alignment length mismatch (%d vs %d(%d-%d)) at line:%s\n",tpos, eff_t_len, al.t_alnend, al.t_alnstart, line);
  if (al.r_alnend-al.r_alnstart!=qpos)
   	GError("Error: ref alignment length mismatch (%d vs %d-%d) at line:%s\n",qpos, al.r_alnend, al.r_alnstart, line);
}

int getRefContext(GASeq& refseq, int rloc, GStr& rctx) {
	//rctx must be empty already!
	int ctxstart=rloc-4;
	int evtloc=4; //local position of event start in rctx
	if (ctxstart<0) { evtloc+=ctxstart; ctxstart=0; }
	else if (ctxstart+8>=refseq.getSeqLen()) {
		evtloc+=refseq.getSeqLen()-ctxstart-9;
		ctxstart=refseq.getSeqLen()-9;
	}
    rctx.append(refseq.getSeq()+ctxstart, 9);
	rctx.upper();
	return evtloc;
}

bool hpolyCheck(TDiffInfo& d, GStr& rctx, int rctxloc) {
	if (d.evtbases.length()>1) {
		char c=d.evtbases[0];
		for(int i=1;i<d.evtbases.length();i++)
			if (c!=d.evtbases[i]) return false;
	}
	//reference context: 4 bases around rloc (total 9 bases)
	char ch=d.evtbases[0];
	GStr cseed(ch);
	cseed.append(ch);cseed.append(ch);cseed.append(ch);
	int l=rctx.index(cseed);
	if (l>=0 && l<=rctxloc && l+cseed.length()>=rctxloc) return true;
	return false;
}

//bool mmotifCheck(GStr& stat, GStr& rctx, int rctxloc) {
int mmotifCheck(GStr& stat, GStr& rctx) {
	int m=0;
	while (metmot[m]!=NULL) {
	  int mpos=rctx.index(metmot[m]);
  	  if (mpos>=0) { //should we check if it's actually including rctxloc?
  		  stat="motif ";
  		  stat.append(metmot[m]);
  		  return m+1;
  	  }
  	  ++m;
    }
	return 0;
}

/*
bool mmotifCheck(GASeq& refseq, TDiffInfo& d, GStr& stat) {
	//for deletions, also search for methylation motifs in evtbases
	if (d.evt=='D') {
	  GStr r_ev(d.evtbases);
	  revCompl(r_ev);
	  int m=0;
      while (metmot[m]!=NULL) {
    	  if (d.evtbases.index(metmot[m])>=0 ||
    		r_ev.index(metmot[m])>=0) {
    		stat="motif ";
    		stat.append(metmot[m]);
    		return true;
    	  }
    	++m;
	  }
	}
//now similarly for the context -- does it even make sense?
	GStr ctx(d.context);
	ctx.upper();
	GStr r_c(ctx);
	revCompl(r_c);
	int m=0;
	while (metmot[m]!=NULL) {
  	  if (ctx.index(metmot[m])>=0 ||
  			r_c.index(metmot[m])>=0) {
  		  stat="motif ";
  		  stat.append(metmot[m]);
  		  stat.append(" ctx");
  		  return true;
  	  }
  	++m;
    }
	return false;
}
*/
void predictImpact(GStr& txt, TDiffInfo& di, GStr& r_trseq, int r_offset) {
  GStr modseq(r_trseq);
  //apply the modification
  if (di.evt=='S') {
	//aminoacid substitution or stop codon introduced here?
	//affected aminoacid locations
	int aaofs=-1;
    GStr orig;
    GStr mod;
    GVec<int> aamods;
    for (int i=0;i<di.evtbases.length();++i) {
    	if (toupper(modseq[di.rloc-r_offset+i])!=toupper(di.evtsub[i]))
    		GError("Error: modseq[%d] not matching di.evtsub[%d] !\n", di.rloc-r_offset+i, i);
    	modseq[di.rloc-r_offset+i]=di.evtbases[i];
        int ao=(di.rloc-r_offset+i)/3;
        if (ao!=aaofs) { //aa affected
        	aaofs=ao;
        	aamods.Add(ao);
        }
    }
    //compare affected aminoacids
    int cdiff=0;
    for (int i=0;i<aamods.Count();++i) {
    	char aa=translateCodon(r_trseq.chars()+ (aamods[i] * 3));
    	char maa=translateCodon(modseq.chars()+ (aamods[i] * 3));
    	if (aa!=maa) {//not a synonymous codon
    		if (cdiff) txt.append(", ");
    		txt.append("AA");
    		int aapos=aamods[i]+di.rloc/3;
    		txt.append(aapos);
    		txt.append('|');
    		txt.append(aa);txt.append(':');txt.append(maa);
    		++cdiff;
    		if (maa=='.') {
    			txt.append("|premature stop at AA");
    			txt.append(aapos);
    		}
    	}
    }
    if (txt.is_empty()) txt.append("synonymous");
    return;
  }
  if (di.evt=='I') {
    //frame shift introducing?
    modseq.insert(di.evtbases, di.rloc-r_offset); //check if it's done right!
  } else if (di.evt=='D') {
	modseq.cut(di.rloc-r_offset, di.evtlen); //check if it's done right!
  }
  else GError("Error: unrecognized editing event (%c)!\n", di.evt);
  //for I/D look for premature stop codons down the road
  int aamodc=0;
  GStr aa4("",4);
  GStr maa4("",4);
  for (int i=0;i+2<modseq.length();i+=3) {
	  char aamod=translateCodon(modseq.chars()+i);
	  if (aamod=='.') {
		  txt.append("premature stop at AA");
		  int aapos=1+(i+r_offset)/3;
		  txt.append(aapos);
		  break;
	  }
	  if (i>0 && aamodc<4) {
		  ++aamodc;
		  if (i+2<r_trseq.length()) {
			  char aa=translateCodon(r_trseq.chars()+i);
			  aa4.append(aa);
		  }
		  maa4.append(aamod);
	  }
  }
  if (txt.is_empty()) {
	  if (!aa4.is_empty() && !maa4.is_empty()) {
		  txt.append("frame shift ");
		  //char a0=translateCodon(r_trseq.chars()+3);
		  //char a1=translateCodon(modseq.chars()+3);
		  //txt.append(a0);txt.append(':');txt.append(a1);
		  //txt.append("; ");
		  txt.append(aa4);txt.append('+');
		  txt.append(':');
		  txt.append(maa4);txt.append('+');
	  }
  }
}

void PAFAlignment::printDiffInfo(FILE* f, GASeq& refseq) {
  double cov=((alninfo.r_alnend-alninfo.r_alnstart)*100.00)/alninfo.r_len;
  if (rlabel.is_empty())
	  fprintf(f, ">%s coverage:%.2f score=%d edit_distance=%d\n",
			  tlabel.chars(), cov, alnscore, edist);
  else
	  fprintf(f, ">%s--%s coverage:%.2f score=%d edit_distance=%d\n",rlabel.chars(),
		  tlabel.chars(), cov, alnscore, edist);
  for (int i=0;i<tdiffs.Count();++i) {
	TDiffInfo& di=tdiffs[i];
	di.evtbases.upper();
	//GStr r_c(di.context);
	int aapos=(int)(di.rloc/3);
	char aa=translateCodon(refseq.getSeq()+3*aapos);
	++aapos;
	GStr status;
	//check for homopolymers at the location
	GStr rctx("",9);
	int rctxloc=getRefContext(refseq, di.rloc, rctx);
	if (hpolyCheck(di, rctx, rctxloc)) status="homopolymer";
	//impact if this edit were applied
	int r_trloc=3*(aapos-2); //start editing before
	if (r_trloc<0) r_trloc=0;
	GStr r_trseq(refseq.getSeq()+r_trloc, di.evtlen);
	//int mmatchidx=0;
	if (status.is_empty())
		//mmatchidx=
		mmotifCheck(status, rctx); // rctxloc);
	GStr impact;
	//if (status.is_empty() || (mmatchidx>0 && strlen(metmot[mmatchidx-1])<5))
	if (!skipCodAn)
	    predictImpact(impact, di, r_trseq, r_trloc);
	if (status.is_empty()) status="[unknown]";
	GStr tcontext(di.tctx);
	const int MAX_EVLEN=12; //maximum event length to display
	if (tcontext.length()>10+MAX_EVLEN) {
		int dlen=tcontext.length()-10;
		tcontext=di.tctx.substr(0,5);
		tcontext.append('[');
		tcontext.append(dlen);
		tcontext.append(']');
		tcontext.append(di.tctx.substr(-5));
	}
	GStr evtbases(di.evtbases);
	if (evtbases.length()>MAX_EVLEN) {
		int dlen=evtbases.length();
		evtbases="[";
		evtbases.append(dlen);
		evtbases.append(']');
	}
	GStr evtsub(di.evtsub);
	if (evtsub.length()>MAX_EVLEN) {
		int dlen=evtsub.length();
		evtsub="[";
		evtsub.append(dlen);
		evtsub.append(']');
	}
	int rloc=di.rloc+1;
	int tloc=alninfo.t_alnstart+di.tloc+1;
	if (di.evt=='S')
    	fprintf(f, "%c\t%d\t%d(%c)\t%s:%s\t%d\t%s\t%s\t%s\t%s\n", di.evt, rloc, aapos, aa,
    			evtsub.chars(),evtbases.chars(),
    			tloc, tcontext.chars(), rctx.chars(), status.chars(), impact.chars());
	else {
    	GStr fmt("%c\t%d\t%d(%c)\t");
    	if (di.evt=='I') fmt.append(":%s\t%d\t%s\t%s\t%s\t%s\n");
    			    else fmt.append("%s:\t%d\t%s\t%s\t%s\t%s\n");
    	fprintf(f, fmt.chars(), di.evt, di.rloc+1, aapos, aa, evtbases.chars(), tloc,
    			tcontext.chars(),
    			rctx.chars(), status.chars(), impact.chars());
	}
  }
}

void printDebugAln(FILE* f, GSeqAlign* a) {
  fprintf(f,">DebugAlign (%d)\n", a->Count());
  a->print(f,'=');
}

