#include <stdio.h>
#include <ctype.h>
#include <math.h>
#include "GArgs.h"
#include "GStr.h"
#include "GHash.hh"
#include "GList.hh"
#include "GapAssem.h"
#define USAGE "Usage:\n\
 pfa2msa <pfa_with_cg_cs> -r <refseq.fa>\n\
    [-o <outfile_maf>]\n\
    \n\
   <pfa_with_cg_cs> is the input PAF file where a single query was aligned\n\
      to many (larger) target sequences using minimap2 --cs -P\n\
   -r a fasta file with the query sequence to use as reference (required)\n\
   -o write output to file <outfile_maf> instead of stdout\n"

#define LOG_MSG_CLIPMAX "Overlap between %s and target %s rejected due to clipmax=%4.2f constraint.\n"
#define LOG_MSG_OVLCLIP "Overlap between %s and %s invalidated by the %dnt clipping of %s at %d' end.\n"
//-------- global variables
bool debugMode=false;
bool removeConsGaps=false;
//char* ref_prefix=NULL;
bool verbose=false;
//int rlineno=0;

FILE* outf=NULL;

GHash<GASeq> seqs(false);

//list of collected MSAs
GList<GSeqAlign> alns(true, true,false);
                 // sorted, free element, not unique

void printDebugAln(FILE* f, GList<GSeqAlign>& aa, int a);

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
  AlnInfo():reverse(2), r_id(NULL), r_len(0), r_alnstart(0), r_alnend(0),
		  t_id(NULL), t_len(0), t_alnstart(0), t_alnend(0) {}
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
  void set(char strand, const char* rid, const char* rlen, const char* rstart, const char* rend,
		  const char* tid, const char* tlen, const char* tstart, const char* tend) {
	r_id=rid;
	t_id=tid;
	reverse = (strand=='-') ? 1 : 0;
	r_len = atoi(rlen);
	r_alnstart=atoi(rstart);
	r_alnend=atoi(rend);
	t_len = atoi(tlen);
	t_alnstart=atoi(tstart);
	t_alnend=atoi(tend);
  }
};

class PAFAlignment {
  //char* linecpy;
  //int lineno; //input file line#
	/*
  char* gpos;
  char* gpos_onref;
  char* gaps;//char string with gaps in this read
  char* gaps_onref;//char string with gaps in reference
  */
 public:
  AlnInfo alninfo;
  GVec<GapData> rgaps; //reference query gaps
  GVec<GapData> tgaps; //target gaps
  char* seqname; //aligned target sequence (read, or genome mapping)
  char* cs; //cs tag value = difference string (short)
  int nm; //number of 1bp edits (edit distance)
  char* cigar;

  int seqlen; //length of this tseq (target sequence, unedited)
  int offset; //offset of this mapping relative to reference
  int clip5; //amount to clip on the left end (0 for PAF on tseq)
  int clip3; //amount to clip on the right end (0 for PAF on tseq)
  char reverse; //0, or 1 if this mapping is reverse complemented
  void parseErr(int fldno, const char* line);
  PAFAlignment(GDynArray<char*>& t, AlnInfo& alni, GASeq& refseq, GStr& tseq, const char* line);
   //this also rebuilds tseq with the target sequence, by transforming refseq according to cs string
  ~PAFAlignment() { GFREE(seqname); GFREE(cs); GFREE(cigar); }
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
 GArgs args(argc, argv, "DGvd:p:r:o:c:");
 int e;
 if ((e=args.isError())>0)
    GError("%s\nInvalid argument: %s\n", USAGE, argv[e]);
 if (args.getOpt('h')!=NULL) GError("%s\n", USAGE);
 debugMode=(args.getOpt('D')!=NULL);
 removeConsGaps=(args.getOpt('G')==NULL);
 verbose=(args.getOpt('v')!=NULL);
 if (debugMode) verbose=true;
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
  if (!outfile.is_empty()) {
     outf=fopen(outfile, "w");
     if (outf==NULL)
        GError("Cannot open file %s for writing!\n",outfile.chars());
     }
   else outf=stdout;
  //************************************************
  s=args.getOpt('r');
  if (s.is_empty()) GError("Error: query sequence file (-r) is required!\n");
  GFastaFile rfa(s.chars());
  FastaSeq faseq;
  FastaSeq* r=rfa.getFastaSeq(&faseq);
  if (r==NULL || faseq.getSeqLen()==0)
	  GError("Error loading FASTA sequence from file %s !\n",s.chars());
  GASeq *refseq=new GASeq(faseq, true); //take over faseq data
  refseq->lowercase();
  GASeq refseq_rc(*refseq);
  refseq_rc.reverseComplement();
  GLineReader* linebuf=new GLineReader(inf);
  char* line;
  alns.setSorted(compareOrdnum);
  refseq->setFlag(GA_FLAG_IS_REF);
  seqs.Add(refseq->getId(),refseq);

  int numalns=0;
  GDynArray<char*> t(24);
  while ((line=linebuf->getLine())!=NULL) {
   if (line[0]=='#') continue;
   bool firstRefAln=false;
   GStr lstr(line);
   int numt=strsplit(line, t, '\t');
   if (numt<15)
	   GError("Error: invalid PAF fline (num. fields=%d):\n%s\n", numt,lstr.chars());
   if (strcmp(t[0], refseq->getId())!=0)
	   GError("Error: expected reference name (%s) in PAF line, but found %s instead!\n",
	    		  refseq->getId(),t[0]);
   AlnInfo al(*t[4], t[0], t[1], t[2], t[3], t[5], t[6], t[7], t[8]);
   if (strcmp(al.r_id, al.t_id)==0) {
	  if (verbose) GMessage("Skipping alignment of qry seq to itself.\n");
      continue; //skip redundant inclusion of reference as its own child
   }

   if (al.r_len!=refseq->getSeqLen())
	    GError("Error: ref seq len in this PAF line (%d) differs from loaded sequence length(%d)!\n%s\n",
		     al.r_len, refseq->getSeqLen(),lstr.chars());
   if (seqs.Find(al.t_id)) {
        GMessage("Warning: target sequence %s alignment already seen before, ignoring other alignments\n",
                             al.t_id);
        continue;
   }
   //-- load the current alignment
   ++numalns;
   GASeq *r_seq = refseq;
   if (al.reverse) r_seq=&refseq_rc;

   GStr tseq("",al.t_alnend-al.t_alnstart+2);
   PAFAlignment* aln=new PAFAlignment(t, al, *r_seq, tseq, lstr.chars());
   //this also fills tseq and sets aln->offset to the tseq mapping offset on refseq
   GStr tlabel(al.t_id, 21);
   tlabel+=":";tlabel+=al.t_alnstart;tlabel+='-';
   tlabel+=al.t_alnend;
   if (al.reverse) tlabel+='-';
   else tlabel+='+';
   //char* t_seq=Gstrdup(tseq.chars());
   GASeq* taseq=new GASeq(tlabel.chars(), "", tseq.chars(), tseq.length(), al.r_alnstart);
   taseq->revcompl=aln->reverse;
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
   /*
   if (taseq->revcompl==1) {
         taseq->reverseGaps();
   }
   */
   GSeqAlign *newmsa=new GSeqAlign(rseq, taseq);
   //this also sets rseq->msa and taseq->msa to newmsa
   if (firstRefAln) { //first alignment with refseq so rseq==refseq
	     //newmsa->incOrd();
	     newmsa->ordnum=numalns;
	     alns.Add(newmsa);
   }
   else { // rseq != refseq, but an instance of refseq in this current aln
	    //incrementally add this alignment, as a MSA, to existing MSA,
	    //propagating gaps through rseq instance
	   refseq->msa->addAlign(refseq, newmsa, rseq);
	   delete newmsa; //incorporated, no longer needed
	   //seqs.Add(taseq->name(),taseq); //this should be added for firstRefAln too, right?
   }

   seqs.Add(taseq->name(),taseq);

   // new pairwise alignment added to MSA
   // debug print the progressive alignment
   if (debugMode) {
    for (int a=0;a<alns.Count();a++) {
    	printDebugAln(outf, alns, a);
      }
    }
   delete aln; //no longer needed
   if (linebuf->isEof()) break;
 } //----> PAF line parsing loop

  //*************** now build MSAs and write them
  delete linebuf;
  if (verbose) {
     //fprintf(stderr, "\n%d input lines processed.\n",rlineno);
     fprintf(stderr, "Refining and printing %d MSA(s)..\n", alns.Count());
     fflush(stderr);
     }

  //print all MSAs found
  for (int i=0;i<alns.Count();i++) {
   GSeqAlign* a=alns.Get(i);
   if (debugMode) { //write plain text alignment file
     fprintf(outf,">Alignment%d (%d)\n",i+1, a->Count());
     a->print(outf, 'v');
     }
   else {//write a real ACE file
     //a->buildMSA();
     GStr ctgname;
     ctgname.format("AorContig%d",i+1);
     a->writeACE(outf, ctgname.chars());
     a->freeMSA(); //free MSA and seq memory
     }
   } // for each PMSA cluster
  // oooooooooo D O N E oooooooooooo
  alns.Clear();
  seqs.Clear();

  fflush(outf);
  if (outf!=stdout) fclose(outf);
  if (inf!=stdin) fclose(inf);
}

//---------------------- RefAlign class
void PAFAlignment::parseErr(int fldno, const char* line) {
 fprintf(stderr, "Error parsing input line (field %d):\n%s\n",
         fldno, line);
 exit(3);
}

const char* CIGAR_ERROR="Error parsing cigar string from line: %s (cigar position: %s)\n";
const char* CS_ERROR="Error parsing cs string from line: %s (cs position: %s)\n";

//this also rebuilds tseq with the target sequence, by transforming refseq according to cs string
PAFAlignment::PAFAlignment(GDynArray<char*>& t, AlnInfo& al, GASeq& refseq, GStr& tseq, const char* line):rgaps(),tgaps(),
		cs(NULL), nm(-1), cigar(NULL) {
  //toffs must be the offset of the alignment start from the beginning of the full, original refseq
  // so it's r_clip5 (or r_clip3 if it's revcomp)
  seqname=Gstrdup(t[5]);
  alninfo=al; //do we still need to keep it? *_id pointers will be obsolete
  reverse=al.reverse;
  clip5=0; //always cutting out the aligned region from the target sequence
  clip3=0;
  offset=al.r_alnstart;
  //int base_ofs=offset; //mismatch base offset -- different for
  if (al.reverse) { //offset is on the reverse complement ref string
	   offset=al.r_len-al.r_alnend;
	   //al.r_alnend=al.r_len-al.r_alnstart;
	   //al.r_alnstart=offset;
  }
  seqlen=al.t_alnend-al.t_alnstart; //check this against the cigar string, and the cs string
  int tnum=t.Count();
  //parse cigar to get the gaps
  for (int i=12;i<tnum;++i) {
	  if (startsWith(t[i], "NM:i:")) {
		  nm=atoi(t[i]+5);
		  if (cigar && cs) break;
		  continue;
	  }
	  if (startsWith(t[i], "cg:Z:")) {
          cigar=Gstrdup(t[i]+5);
          if (cs && nm>=0) break;
          continue;
	  }
	  if (startsWith(t[i], "cs:Z:")) {
          cs=Gstrdup(t[i]+5);
          if (cigar && nm>=0) break;
          continue;
	  }
  }
  if (cigar==NULL || cigar[0]==0) GError(CIGAR_ERROR, line, 0);
  // gpos = current genomic position (will end up as right coordinate on the genome)
  // rpos = read position (will end up as the length of the read)
  // cop = CIGAR operation, cl = operation length
  int mbases = 0; //count "aligned" bases (includes mismatches)
  int qpos = 0;  //should match length of aligned region of reference qseq, minus offset
  int tpos = 0; //will end up with the length of target genome region being aligned
  int eff_t_len=al.t_alnend-al.t_alnstart; //effective target sequence length
  //target sequence is to be trimmed to the aligned region
  GapData gap;
  char *p=cigar;
  while (*p != '\0') {
	  int cl=0;
	  if (!parseInt(p,cl)) GError(CIGAR_ERROR, line, p);
	  char cop=*p;
	  if (cop=='\0') GError(CIGAR_ERROR, line, p);
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
		 qpos+=cl;
		 if (reverse) gap.pos=eff_t_len-tpos;
		 else gap.pos=tpos;
		 gap.len=cl;
		 tgaps.Add(gap);
		 break;
	   case 'D':
		 //deletion in reference sequence relative to the query (gap in query ref seq)
		 tpos += cl;
		 gap.pos=offset+qpos;
		 if (reverse) //actual location on qry ref for a reverse complement match:
			  gap.pos=al.r_len-gap.pos;
		 gap.len=cl;
		 rgaps.Add(gap);
		 break;
	   case 'N':
		 // intron
		 //special skip operation, not contributing to "edit distance",
		 //   so num_mismatches should not update
		 tpos+=cl;
		 //shouldn't really happen in genomic PAF
		 gap.pos=offset+qpos;
		 if (reverse) //actual location on qry ref for a reverse complement match:
			  gap.pos=al.r_len-gap.pos;
		 gap.len=cl;
		 rgaps.Add(gap);
		 break;
	   default:
		 GError("Error: unhandled cigar_op %c (len %d) in %s\n", cop, cl, line);
	  } //switch cigar op
	++p;
  } // interpret_CIGAR string
  if (eff_t_len!=tpos)
   	GError("Error: tseq alignment length mismatch (%d vs %d(%d-%d)) at line:%s\n",tpos, eff_t_len, al.t_alnend, al.t_alnstart, line);
  if (al.r_alnend-al.r_alnstart!=qpos)
   	GError("Error: ref alignment length mismatch (%d vs %d-%d) at line:%s\n",qpos, al.r_alnend, al.r_alnstart, line);

  //interpret cs string to fill tseq
  tseq="";
  p=cs;
  qpos=0; //ref query location in the alignment
  tpos=0;
  char tch=0,qch=0;
  int q_pos=0;
  while (*p != 0) {
    char op=*p;
    int cl=0;
    p++;
    switch (op) {
      case ':':
    	if (!parseInt(p,cl)) GError(CS_ERROR, line, p);
    	tseq.appendmem(refseq.getSeq()+offset+qpos, cl);
    	qpos+=cl;
    	tpos+=cl;
    	break;
      case '*': //substitution
	    tch=tolower(*p);++p;
	    qch=tolower(*p);++p;
	    q_pos=offset+qpos;
	    if (qch!=refseq.getSeq()[q_pos]) {
	    	GError("Error: base mismatch %c != qstr[%d] (%c) at line\n%s\n",
	    			qch, q_pos, refseq.getSeq()[q_pos], line);
	    }
	    tseq.append((char)toupper(tch));
	    ++qpos;
	    ++tpos;
   	    break;
      case '-': //gap in qseq (insertion in q, deletion in tseq)
    	//this insertion is at position (offset+qpos) in ref query if reverse==0
        // BUT at rlen-(offset+qpos) if reverse==1
    	while (isalpha(tch=tolower(*p))) {
    		  ++p;
    		  ++tpos;
    		  tseq.append((char)toupper(tch));
    	}
    	//qpos unchanged, offset+qpos is the location of the gap in qseq
    	break;
      case '+': //gap in tseq (insertion in tseq, deletion in qseq)
    	while (isalpha(qch=tolower(*p))) {
    		  ++p;
    		  ++qpos;
    	}
    	//tpos unchanged, it is the location of the gap in tseq
    	break;
      case '~':
    	GError("Error: spliced alignments not supported! at line:\n%s\n", line);
    	break;
      default:
        GError("Error: unhandled event at %s in cs, line:\n%s\n",p,line);
    }
  } //interpret CS string;
 }

void printDebugAln(FILE* f, GList<GSeqAlign>& aa, int a) {
  ++a;
  fprintf(f,">[%d]DebugAlign (%d)\n", a, aa[a-1]->Count());
  aa[a-1]->print(f,'=');
}


