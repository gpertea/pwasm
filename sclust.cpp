//un/comment this for tracing options
//#define DBGTRACE 1
#include <stdio.h>
#include <ctype.h>
#include "GBase.h"
#include "GArgs.h"
#include "GStr.h"
#include "GHash.hh"
#include "GList.hh"

#ifdef DBGTRACE
#include <gcl/GShMem.h>
#endif
/* TO DO:
 -w  : sub-cluster bridge detection threshold: when a cluster has at least \n\
       2 \"seed\" sequences having both at least the specified\n\
       weight <c_value>, write about it in <warning_file> (default: 10)\n\
*/

#define usage "Performs weighted seeded clustering by filtering tabulated hits.\n\
The hits are expected to be sorted in descending order, by score\n\
Usage:\n\
 sclust [<hits_file>] [-H] [-o <out_file>] [-L] [-C] [-x <excludefile>]\n\
 [-s <only_these>] [SEQFLT={ET|EST|ET2EST}] [HEAVY=xxxx] [ETBONUS=xxxx]\n\
 [SCOV=xx] [SCORE=xx] [PID=xx] [MAXCL=nnnnn]\n\
Parameters:\n\
HEAVY=XXXX mininum score of a seed to be considered \"heavy\" so its \n\
          cluster be not allowed to merge with other \"heavy\" clusters. \n\
          Default value is 5000.\n\
ETBONUS=XXXX default: 3000\n\
 -H  : disable the fasta-style header for output clusters\n\
 -o  : write output in <out_file> instead of stdout\n\
 -L  : redirect stderr to file log_<outfile>\n\
 -C  : do not compact/adopt clusters before displaying results \n\
\n\
Optional pairing filters: \n\
 -s     :  only consider pairs/hits between sequences in <only_these>\n\
 -x     :  discard any pairs/hits related to sequences listed in <excludefile>\n\
SEQFLT=     filter sequence names for specific sequence type clustering: \n\
            EST  - EST-only clustering (use only ESTvsEST hits) \n\
            ET   - ETs/NPs only clustering\n\
            ET2EST - ignore EST-EST pairs, uses only EST-ET or ET-ET pairs\n\
SCOV=xx     is the minimum overlap length from the shorter sequence (percentage)\n\
            (minimum overlap length as percentage of the shorter sequence)\n\
LCOV=xx     is the minimum overlap length from the longer sequence (percentage)\n\
            (minimum overlap length as percentage of the longer sequence)\n\
OVHANG=xx   is the maximum overhang (in nucleotides) to be allowed\n\
            for any overlap(default 5000)\n\
PID=xx      is the minimum percent of identity\n\
SCORE=xx    is the minimum score for a match to be considered\n\
MAXCL=xxxxx do not allow clusters to grow past this number of sequences\n\
            (default 40000)\n\
\n\
If <hits_file> is not given, hits are expected at stdin\n\
Each input line must have these tabulated fields:\n\n\
q_name q_len q_n5 q_n3 hit_name hit_len hit_n5 hit_n3 pid score score2 strand\n\
"

int maxcl=60000;

int heavy_score=5000;
int ET_bonus=3000;
bool sep_seeds=true;
bool compact_clusters=true;

unsigned char isET(char* name) {
 return (startsWith(name, "np|") ||
        startsWith(name, "et|") || 
        startsWith(name, "egad|") ||
        startsWith(name, "preegad|")
        ) ? 1 : 0; 
 }

int comparePtr(void* p1, void* p2) {
  return (p1>p2) ? -1 : ((p1<p2)?1:0);
  }

class CNode;

CNode* tracenode=NULL;

struct CLinkData {
 char* id1;
 char* id2;
 CNode* n1;
 CNode* n2;
 int len1;
 int len2;
 int ovl1;
 int ovl2;
 int cov1;
 int cov2; 
 int scov;
 int lcov;
 int score;
 int pid;
 char ext1;
 char ext2;
 bool minus;
 CLinkData(char* seq1, char* seq2, int l1, int l2, int o1, int o2, int sc, int p_id, char x1, char x2, bool neg) {
   id1=seq1;id2=seq2;len1=l1;len2=l2;ovl1=o1;ovl2=o2;score=sc; pid=p_id;
   ext1=x1;
   ext2=x2;
   minus=neg;
   n1=NULL;n2=NULL;
   cov1=(ovl1*100)/len1;
   cov2=(ovl2*100)/len2;
   scov=(len1>len2)?cov2:cov1;
   lcov=(len1<len2)?cov2:cov1;
   }
 };

struct BLink {
 CNode* n1;
 CNode* n2;
 int score;
 int scov;
 int pid;
 BLink( CNode* nn1, CNode* nn2, int s, int v, int p) {
  n1=nn1;n2=nn2;scov=v, score=s; pid=p;
  }
 BLink(CLinkData &l) {
  n1=l.n1;n2=l.n2;scov=l.scov, score=l.score; pid=l.pid;
  }
};


class GCluster :public GList<CNode> {
 public:
  bool isHeavy;
  bool stored; //stored in tclusters 
  int weight; //it's only the weight of the heaviest sequence in the cluster
  int total; //total number of sequences (accounts for embedded ones too!)
  BLink* bestLink; //link to a node from another cluster
  CNode* maxnode;
  bool operator==(GCluster& d){
     return (this==&d);
     }
  bool operator>(GCluster& d){
     return (this>&d);
     }
  bool operator<(GCluster& d){
     return (this<&d);
     }
  GCluster():GList<CNode>(true,false,false) {
   //default is: sorted by Node pointer, don't free nodes, unique
    isHeavy=false;
    weight=0;
    total=0;
    maxnode=NULL;
    stored=false;
    bestLink=NULL;
    }
  GCluster(bool sorted, bool free_elements=false, bool beUnique=false):
     GList<CNode>(sorted,free_elements,beUnique) {
     isHeavy=false;
     weight=0;
     total=0;
     maxnode=NULL;
     stored=false;
     bestLink=NULL;
     }
  void addSeq(CNode* n, CLinkData* pl=NULL);
  void setBestLink(CLinkData& l);
  bool betterLink(CLinkData& l);
};



class CNode {
 public:
   char* id; //actual node ID
   int len; //sequence length
   int weight;
   bool isFull;
   //unsigned char isFull;
   //unsigned char cFlag;
   unsigned char is_ET;
   /*
   char cl_ext; // -1=backward tail, +1=forward tail, 0=undefined, 2 = in the middle!
   bool cl_minus; //(true if the overlap is anti-sense with maxnode)
   */
   //links to high weight seeds are kept here, and not added to the t-cluster!
   
   GCluster* cluster; /* pointer to a list of nodes (the t-cluster
                         containing this sequence) */
   CNode(char* name, int slen) {
     is_ET=isET(name);
     len=slen;
     weight=len/10 + is_ET ? ET_bonus : 0;
     isFull=(startsWith(name, "np|") || startsWith(name, "et|"));
     //cFlag=0;
     id=Gstrdup(name);    

     cluster=new GCluster();
     cluster->addSeq(this);
    }

   //add adirectly in another sequence's cluster
   CNode(char* name, CNode* pairnode, int slen) {
     is_ET=isET(name);
     isFull=(startsWith(name, "np|") || startsWith(name, "et|"));
     //cFlag=0;
     len=slen;
     //cl_ext=0;     
     //cl_minus=false;
     weight=len/10 + is_ET ? ET_bonus :0;
     id=Gstrdup(name);
     cluster=pairnode->cluster; //same cluster
     cluster->addSeq(this); //add self to the cluster
    }

   ~CNode() {
     }
  //comparison: just pointer based:
    bool operator==(CNode& d){
     return (this==&d);
     }
    bool operator>(CNode& d){
     return (this>&d);
     }
    bool operator<(CNode& d){
     return (this<&d);
     }
 };

//-- global hash containing all the nodes with structured links info
// holds all the nodes
// associated with CNode (its links and its clusters)


GHash<CNode> all;

 #ifdef DBGTRACE
   GShMem mem("tclust", 4096);
   GShMem memsum("tclustsum", 1024);
   GStr slog;
 #endif


bool asfasta;
bool flt_ET_only=false;
bool flt_EST_only=false;
bool flt_EST2ET=false;
bool ETcheck=false;
bool flt_seqOnly=false;

FILE* outf=NULL;
FILE* fwarng=NULL;
FILE* flinks=NULL;

GHash<int> xcludeList(false);
GHash<int> seqonlyList(false);

//all clusters are stored in a unique sorted list
GHash<CNode> seedList(false); //do not try to free the nodes
GList<GCluster> tclusters(true,false,false);

GList<GCluster> wasted(false,true,false);

//=========== functions used:
bool addPair(CLinkData& pair);
int compareCounts(void* p1, void* p2);
int compareWeights(void* p1, void* p2);
int compareID(void* p1, void* p2); //compare node IDs
int getNextValue(char*& p, char* line);
char* getNextWord(char*& p, char* line);
bool seq_filter(char* s1, char* s2); //returns true if the pair can pass through the filter
int readNames(FILE* f, GHash<int>& xhash);

void showCl(GCluster* C) {
  GMessage("{%s",C->Get(0)->id);
  for (int i=1;i<C->Count();i++) 
     GMessage(" %s",C->Get(i)->id);
  GMessage("}");
}


//========================================================
//====================     main      =====================
//========================================================
#define BUF_LEN 1024 //maximum input line length 
#define ERR_INVALID_PAIR "Invalid input line encountered:\n%s\n"
#define ERR_SEED_PARSE "Error parsing seed file at line:\n%s\n"

int main(int argc, char * const argv[]) {
 char inbuf[BUF_LEN]; // incoming buffer for sequence lines.
 GArgs args(argc, argv, 
     "htHLCSs:o:l:x:w:HEAVY=ETBONUS=SEQFLT=PID=SCOV=LCOV=OVHANG=SCORE=MAXCL=");
 int e;
 if ((e=args.isError())>0)
    GError("%s\nInvalid argument: %s\n", usage, argv[e]);
 if (args.getOpt('h')!=NULL) GError("%s\n", usage);

 GStr infile;
 if (args.startNonOpt()) {
        infile=args.nextNonOpt();
        GMessage("Given file: %s\n",infile.chars());
        }
 if (args.getOpt('L')!=NULL && args.getOpt('o')!=NULL) {
   GStr flog(args.getOpt('o'));
   flog="log_"+flog;
   //stderr=fopen(flog.chars(),"w");
   }
 GMessage("Command line was:\n ");
 for (int i=0;i<argc-1;i++) 
   GMessage("%s ", argv[i]);
 GMessage("%s\n", argv[argc-1]);
                
 asfasta=(args.getOpt('H')==NULL);
 sep_seeds=(args.getOpt('S')!=NULL);
 compact_clusters=(args.getOpt('C')==NULL);
 //tabflt=(args.getOpt('t')!=NULL); 
 int minpid=0, minscov=0, minlcov=0, minscore=0, maxovhang=1000; 
 
 
 GStr s=args.getOpt('x');
 FILE* fxclude=NULL;
 int c;
 if (!s.is_empty()) {
   if ((fxclude=fopen(s, "r"))==NULL)
      GError("Cannot open exclusion file '%s'!\n", s.chars());
   c=readNames(fxclude, xcludeList);
   GMessage("Loaded %d sequences to be excluded.\n", c);
   }

 /*FILE* fseeds=NULL;
 s=args.getOpt('r');
  if (s.is_empty())
   GError("%s\n Seed sequence file was not provided.\n", usage);   */
   
 /*if (!s.is_empty()) {
   if ((fseeds=fopen(s, "r"))==NULL)
      GError("Cannot open seed file '%s'!\n", s.chars());
        
   c=readSeeds(fseeds);
   GMessage("Loaded %d seed sequences.\n", c);
   }*/
   
     
 FILE* fseqonly=NULL;
 s=args.getOpt('s');
 if (!s.is_empty()) {
   if ((fseqonly=fopen(s, "r"))==NULL)
      GError("Cannot open exclusive pairs file '%s'!\n", s.chars());
   int c=readNames(fseqonly, seqonlyList);
   GMessage("Loaded %d sequence names - only consider pairs including them.\n", c);
   flt_seqOnly=(c>0);
   }
  
   
   
 GStr wfile=args.getOpt('w'); //write warnings here
 if (wfile.is_empty())
    wfile="sclust_warnings";
 /*if ((fwarng=fopen(wfile, "w"))==NULL)
      GError("Cannot write warnings file '%s'!\n", 
          wfile.chars());
*/
  
 s=args.getOpt("SEQFLT");
 if (!s.is_empty()) {
   if (s=="ET") {
     flt_ET_only=true;
     s=" ETs/NPs only (discard ESTs)";
     }
   else if (s=="EST") {
     flt_EST_only=true;
     s=" EST only (discard ETs/NPs)";
     }
   else if (s=="EST2ET" || s=="ET2EST")  {
     flt_EST2ET=true;
     s=" links to ETs/NPs only (discard EST-EST links)";
     }
   else GError("%s\nIncorrect SEQFLT option ('%s' not recognized)", s.chars());
   GMessage("Using sequence filter: %s\n", s.chars());
   ETcheck=true;
   }
 s=args.getOpt("HEAVY"); 
 if (!s.is_empty()) {
    //tabflt=true;
    heavy_score = s.asInt();
   }
 GMessage("HEAVY score=%d\n",  heavy_score);
 s=args.getOpt("ETBONUS"); 
 if (!s.is_empty()) {
    //tabflt=true;
    ET_bonus = s.asInt();
   }
 GMessage("ET BONUS=%d\n",  ET_bonus);
 s=args.getOpt("PID"); 
 if (!s.is_empty()) {
    //tabflt=true;
    minpid = s.asInt();
    if (minpid>0) GMessage("PID=%d\n",  minpid);
   }
 s=args.getOpt("SCOV");   
 if (!s.is_empty()) {
    //tabflt=true;
    minscov = s.asInt();
    if (minscov>0) GMessage("SCOV=%d\n",  minscov);
   }
 s=args.getOpt("MAXCL");
 if (!s.is_empty()) {
    //tabflt=true;
    maxcl = s.asInt();
    if (maxcl>0) GMessage("Maximum cluster size: %d\n",  maxcl);
    }
 s=args.getOpt("LCOV");   
 if (!s.is_empty()) {
    //tabflt=true;
    minlcov = s.asInt();
    if (minlcov>0) GMessage("LCOV=%d\n",  minlcov);
   }
 s=args.getOpt("OVHANG");
 if (!s.is_empty()) {
    maxovhang = s.asInt();
   }
 if (maxovhang<1000) GMessage("OVHANG=%d\n",  maxovhang);
 s=args.getOpt("SCORE");   
 if (!s.is_empty()) {
    //tabflt=true;
    minscore = s.asInt();
    if (minscore>0) GMessage("SCORE=%d\n",  minscore);
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
   
 //======== main program loop
 char* line;
 //if (tabflt) 
   {
    while ((line=fgets(inbuf, BUF_LEN-1,inf))!=NULL) {
        int l=strlen(line);
        if (line[l-1]=='\n') line[l-1]='\0';
        GStr linecpy(line);
        if (strlen(line)<=1) continue;
        char* tabpos=line;
        //find the 1st tab
        tabpos=strchr(line, '\t');
        if (tabpos==NULL || tabpos==line)
            GError(ERR_INVALID_PAIR, line);
        *tabpos='\0'; //so line would be the first node name

        if (xcludeList.hasKey(line)) continue;
        tabpos++; //tabpos is now on the first char of the second field (q_len)
        int score, score2, scov, lcov, ovh_l, ovh_r, pid;
        //skip 3 other tabs delimited
        //read the query length:
        int qlen=getNextValue(tabpos, line);
        tabpos++;
        int q5=getNextValue(tabpos,line);
        tabpos++;
        int q3=getNextValue(tabpos,line);
        tabpos++;
        if (q5==q3) GError(ERR_INVALID_PAIR, line);
        bool minus=false;
        if (q5>q3) {
          minus=true;
          Gswap(q5,q3);
          }
        //now we should be on the first char of the hitname field
        /*
        while (isspace(*tabpos)) tabpos++; //skip any spaces in this second node name field
        if (*tabpos=='\0') GError(ERR_INVALID_PAIR, line);
        //add a string termination after this field
        char* p=tabpos; while (!isspace(*p) && *p!='\0') p++;
        *p='\0';
        */
        char* p=tabpos;
        tabpos=getNextWord(p,(char*)linecpy.chars());
        //now tabpos contains the exact second sequence string
        //if (strstr(line, TRACESEQ)!=NULL || strstr(tabpos, TRACESEQ)!=NULL) {

        if (strcmp(line, tabpos)==0) {
          GMessage("Warning: self pairing found for node %s\n",line);
          continue;
          }
         
        if (xcludeList.hasKey(tabpos)) continue;        

        if (!seq_filter(line, tabpos)) continue;

        p++; //move on the first char of the hitlen 
        int hitlen=getNextValue(p,line);
        p++;
        int h5=getNextValue(p,line);
        p++;
        int h3=getNextValue(p,line);
        p++;
        if (h5==h3) GError(ERR_INVALID_PAIR, line);
        if (h5>h3) {
           Gswap(h5,h3);
           minus=!minus;
           }
        pid = getNextValue(p,line); p++;
        score = getNextValue(p,line);p++;
        score2 = getNextValue(p,line);
        //compute coverages:
        char q_ext=0;
        char h_ext=0;        
        if (minus) {
          ovh_r=GMIN(q5-1, hitlen-h3);
          ovh_l=GMIN(h5-1, qlen-q3);
          /*setting q_ext (how hit extends it: backward=-1, forward=+1)          
          if (hitlen-h3-q5+1>100 && h5-1<100)
                q_ext=-1;
              else if (h5-1-qlen+q3>100 && hitlen-h3<100) 
                    q_ext=1;
                 else if (hitlen-h3-q5+1>100 && h5-1-qlen+q3>100)
                        q_ext=2;
          if (qlen-q3-h5+1>100 && q5-1<100)
                h_ext=-1;
              else if (q5-1-hitlen+h3>100 && qlen-q3<100) h_ext=1;
                else if (qlen-q3-h5+1>100 && q5-1-hitlen+h3>100)
                       h_ext=2;*/
          }
         else {
          ovh_r=GMIN(hitlen-h3, qlen-q3);
          ovh_l=GMIN(h5-1, q5-1);
          /*
          if (h5-1-q5+1>100 && hitlen-h3<100)
              q_ext=-1;
            else if (hitlen-h3-qlen+q3>100 && h5-1<100)
                       q_ext=1;
            else if (h5-1-q5+1>100 && hitlen-h3-qlen+q3>100)
                   q_ext=2;           
          if (q5-1-h5+1>100 && qlen-q3<100)
              h_ext=-1;
            else  
             if (qlen-q3-hitlen+h3>100 && q5-1<100)
                       h_ext=1;
            else if (q5-1-h5+1>100 && qlen-q3-hitlen+h3>100)
                        h_ext=2;
                        */
          }
        /*ovh_r=minus ?(GMIN(q5-1, hitlen-h3)) :(GMIN(hitlen-h3, qlen-q3));
         ovh_l=minus ?(GMIN(h5-1, qlen-q3)) :(GMIN(h5-1, q5-1));*/
        if (hitlen>qlen) { //query is shorter
          scov = (int) ((double)(q3-q5+1)*100)/qlen;
          lcov = (int) ((double)(h3-h5+1)*100)/hitlen;
          }
        else { //hit sequence is shorter
          lcov = (int) ((double)(q3-q5+1)*100)/qlen;
          scov = (int) ((double)(h3-h5+1)*100)/hitlen;
          }
          
        if (scov>=minscov && lcov>=minlcov && pid>=minpid 
             && score>=minscore && ovh_r <= maxovhang && ovh_l <= maxovhang) {
           //GMessage("%s(%d) %d-%d  | %s(%d) %d-%d, pid=%d%%, score=%d, scov=%d%%, lcov=%d%%\n",
           //          line,qlen,q5,q3,tabpos,hitlen, h5,h3,pid, score, scov, lcov);
           //if (fpairs!=NULL) fprintf(fpairs, "%s\n", linecpy.chars());
           //if (strstr(line, TRACESEQ)!=NULL || strstr(tabpos, TRACESEQ)!=NULL)
           // GMessage("Adding pair '%s' - '%s'  \n",line, tabpos);

           CLinkData pair(line, tabpos, qlen, hitlen, q3-q5+1, h3-h5+1, score2, pid, q_ext, h_ext, minus);
           if (addPair(pair)) {
               if (flinks!=NULL)
                   fprintf(flinks, "%s\n", linecpy.chars());
               }
           }
      } //while
     } //tabulated hits case
  if (inf!=stdin) fclose(inf);
  //if (fpairs!=NULL && fpairs!=stdout) fclose(fpairs);
  //fclose(fwarng);
  //==========
  //
  GStr outfile=args.getOpt('o');
  if (!outfile.is_empty()) {
     outf=fopen(outfile, "w");
     if (outf==NULL)
        GError("Cannot open file %s for writing!\n",outfile.chars());
     }
   else outf=stdout;
  //scan tclusters; for small clusters 
  //link them to their best link found (adopt singletons)
  //tclusters.setSorted(&compareWeights); //heavier clusters first
  if  (compact_clusters) {
    tclusters.setSorted(false);    
    GMessage("Compacting small clusters...\n");
    bool was_merge;
    do {
      was_merge=false;
      for (int i=0;i<tclusters.Count(); i++) {
	GCluster* cl=tclusters[i];
	GCluster* b=NULL;
	/* if (cl!=NULL && cl->maxnode==NULL) GError("maxnode is NULL!\n"); */
	if (cl==NULL) continue;

	if (cl->total==1) {
               tclusters.Put(i,NULL);
               }
	if (cl->weight>ET_bonus || (cl->maxnode!=NULL && cl->maxnode->is_ET==1)) {
           continue; //no adoption for clusters with ETs or already big enough
           }
	if (cl->bestLink==NULL || cl->bestLink->n2->cluster==cl || cl->bestLink->n2->cluster==NULL) {
        	 //if (cl->total==1) tclusters.Put(i,NULL);
        	 continue; //no merging info..
        	 }
	b=cl->bestLink->n2->cluster; // cluster to merge to
	was_merge=true;
	for (int j=0;j<cl->Count();j++)
        	 b->addSeq(cl->Get(j));
	tclusters.Put(i,NULL); //mark the current cluster as no good for merging
	cl->Clear();
	cl->total=0;
	wasted.Add(cl);
	}
      tclusters.Pack();
      } while (was_merge);
    }  
  tclusters.setSorted(&compareCounts); //largest clusters first
  GMessage("Total clusters: %d \n", tclusters.Count());
  if (tclusters.Count()>0) {
     GMessage("Largest cluster has %d nodes\n", tclusters[0]->total);
     for (int i=0; i<tclusters.Count() && tclusters[i]->total>1; i++) {
        GCluster* c=tclusters[i];
        c->setSorted(&compareID);
        if (asfasta) {
          if (sep_seeds)
             fprintf(outf,">SCL%d\t%d\t%d\n", i+1, c->total, c->Count());
            else
             fprintf(outf,">SCL%d\t%d\n", i+1, c->total);
          }
        if (sep_seeds) {
          for (int j=0; j<c->Count(); j++) {
            CNode* n=c->Get(j);
            fprintf(outf, "%s\n", n->id);
            }
          }
         else { //all the sequence names on one line
          fprintf(outf,"%s",c->Get(0)->id);
          for (int j=1; j<c->Count(); j++) {
            fprintf(outf," %s",c->Get(j)->id);
            }
          fprintf(outf,"\n");
          }
        }
     fflush(outf);
     } //for each cluster
  //GMessage("all.Clear:\n");   
  all.Clear(); //this frees the nodes
  //GMessage("tclusters.Clear:\n");
  tclusters.setFreeItem(true);   
  tclusters.Clear(); /* this frees the clusters, 
     but should not try to free the nodes inside */
  if (outf!=stdout) fclose(outf);
  GMessage("*** all done ***\n");
  //getc(stdin);
}

//
int getNextValue(char*& p, char* line) {
 //expects p to be on the first char of a numeric field
 char* start, *end;
 while (isspace(*p) && *p!='\0') p++; //skip any spaces;
 if (*p=='\0') GError(ERR_INVALID_PAIR, line);
 start=p;
 while (!isspace(*p) && *p!='\0') p++;
 //now p should be on the next space or eos;
 end=p-1;
 //while (isspace(*end)) end--;
 *(end+1)='\0';
 double d=atof(start);
 return (int)d;
 }
 
//
char* getNextWord(char*& p, char* line) {
 //expects p to be on the first char of a numeric field
 char* start, *end;
 while (isspace(*p) && *p!='\0') p++; //skip any spaces;
 //if (*p=='\0') GError(ERR_INVALID_PAIR, line);
 if (*p=='\0') return p;
 start=p;
 while (!isspace(*p) && *p!='\0') p++;
 //now p should be on the next space or eos;
 end=p-1; 
 //while (isspace(*end)) end--;
 //char c=*(end+1);
 *(end+1)='\0';
 return start;
 //backup
 }
 
//=====
int compareCounts(void* p1, void* p2) {
 int c1=((GCluster*)p1)->total;
 int c2=((GCluster*)p2)->total;
 return (c1>c2)?-1:((c1<c2)?1:0);
}

int compareWeights(void* p1, void* p2) {
 int c1=((GCluster*)p1)->weight;
 int c2=((GCluster*)p2)->weight;
 return (c1>c2)?1:((c1<c2)?-1:0);
}


//=====
int compareID(void* p1, void* p2) {
 return strcmp(((CNode*)p1)->id,((CNode*)p2)->id);
}

void storeCluster(GCluster* cl) {
 tclusters.Add(cl);
 cl->stored=true;
}

bool badLink(CLinkData& l) {
 if (l.n1->cluster==l.n2->cluster && l.n1->cluster!=NULL) return true;
 //don't directly link to a full NP unless there is at least 95%, 80% coverage 
 //of the shorter read
 if ((l.n1->isFull || l.n2->isFull) && l.pid<94 && l.scov<75) return true;
 //don't directly link two seeds
 //unless there is at least 75% coverage for the shorter sequence 
 //with pid>=98
 if (l.n2->weight>=ET_bonus+l.len2 && l.n1->weight>=ET_bonus+l.len1 && 
    (l.pid<97 || l.scov < 80)) return true;
      
 //unconditionally attach a small EST cluster to anything...
 if ((l.n1->weight<ET_bonus && l.n1->cluster->total<100 && l.n1->cluster->weight<ET_bonus) ||
     (l.n2->weight<ET_bonus && l.n2->cluster->total<100 && l.n2->cluster->weight<ET_bonus)) {
        //if (traceLink) GMessage("    Link to small EST cluster ACCEPTED !\n");
        return false;
        }
 //don't link an ordinary seq, already in a good cluster, to a seed, unless more than 60% of
 //the incoming seq overlaps with more than 96% pid 
 if (l.n1->weight>=ET_bonus+l.len1 && l.n2->weight<ET_bonus && l.n2->cluster->weight>ET_bonus &&
      (l.cov2 < 65 || l.pid<96)) return true;
 if (l.n2->weight>=ET_bonus+l.len2 && l.n1->weight<ET_bonus && l.n1->cluster->weight>ET_bonus &&
      (l.cov1 < 65 || l.pid<96)) return true; 


 //don't link two ordinary seqs if any of them is in an cluster dominated by a full NP
 //unless the link is at least 60% coverage with PID>96
 if (!l.n1->isFull && !l.n2->isFull &&
     (l.n1->cluster->maxnode->isFull || l.n2->cluster->maxnode->isFull) 
        //&& (l.cov1<65 || l.cov2<65 || l.pid<96)) 
        ) return true;

 
 //don't link two ordinary seqs if any of them is in an ET cluster 
 // unless the link is at least 60% coverage with PID>96
 if (l.n1->is_ET==0 && l.n2->is_ET==0 &&
        //(l.n1->cluster->weight>=ET_bonus+l.n1->cluster->maxnode->len  || 
        // l.n2->cluster->weight>=ET_bonus+l.n2->cluster->maxnode->len) &&
       (l.n1->cluster->maxnode->is_ET || l.n2->cluster->maxnode->is_ET) 
        && (l.cov1<65 || l.cov2<65 || l.pid<96))
        return true;

 //don't link two ordinary seqs if at least one them is in a heavy cluster!
 if (l.n1->weight<ET_bonus && l.n2->weight<ET_bonus &&
      (l.n1->cluster->isHeavy || l.n2->cluster->isHeavy)) return true;

 //don't link too clusters if the resulting cluster is over a max count threshold
 if (l.n1->cluster->total+l.n2->cluster->total>maxcl) {
   GMessage("Clusters containing %s and %s were kept apart because of MAXCL(%d) size limit reached.\n",
      l.n1->cluster->maxnode->id, l.n2->cluster->maxnode->id, maxcl);
   return true;
   }

 //don't link two heavy clusters, unless the link is performed by the best seqs in that cluster
 if (l.n1->cluster->isHeavy && l.n2->cluster->isHeavy && 
   (l.n1->cluster->maxnode!=l.n1 || l.n2->cluster->maxnode!=l.n2 || l.scov<75 || l.pid<98))
   return true;
 return false;
}

//create a link between two existing nodes
//(unless they are both 2 full-length seeds)
bool linkNodes(CLinkData& p) {
 bool badlink=badLink(p);
 if (badlink) {
   //for each of these clusters, keep this link for later (small cluster adoption) 
   if (p.n1->cluster->betterLink(p))
      p.n1->cluster->setBestLink(p);
   if (p.n2->cluster->betterLink(p))
      p.n2->cluster->setBestLink(p);
   return false; //ignore the link if it doesn't look good
   }

//we have to merge the two clusters
//add to the bigger cluster the sequences from the smaller one
//(big fish swallows the smaller one)
GCluster* src, *dest;
if (p.n1->cluster->total>p.n2->cluster->total) {
  // n1 is the bigger fish
  dest = p.n1->cluster;
  src  = p.n2->cluster;
  }
 else {
  dest = p.n2->cluster;
  src  = p.n1->cluster;
  }
//merge the clusters (in the bigger one)
//and remove the smaller cluster  
for (register int i=0; i<src->Count(); i++) {
     CNode* n= src->Get(i);
     dest->addSeq(n);  /* check duplicates? they should never be there!
          if (dest->Add(n)<0)  {
                   delete n;
                   src->Put(i,NULL);
                   } */
     //each node from smaller cluster should now point to the bigger (current)cluster
     n->cluster=dest;
     }
 dest->isHeavy = dest->isHeavy || src->isHeavy;
 dest->weight = GMAX(dest->weight, src->weight);
 if (!dest->isHeavy && dest->weight + dest->total/2 > heavy_score)
        dest->isHeavy=true;
 int idx;
 if (src->stored && (idx=tclusters.IndexOf(src))>=0)
        tclusters.Delete(idx);
 src->Clear();
 wasted.Add(src);
 if (!dest->stored)
      storeCluster(dest); //this will not add if already there, but it's better to check
 //also set the maxnode extension sense if provided
 return true;
}

//-----
//both nodes are new - create global entries for them
//and a new shared t-cluster containing them
// so they are  NOT seeds and they passed the filters
void newPair(CLinkData& l) {
 l.n1=new CNode(l.id1, l.len1); // this creates a new cluster with only n1 in it
 l.n2=new CNode(l.id2, l.len2); // 
 //p.n1->bestHit=p.n2;
 //p.n1->bestHitScore=p.score;
 all.shkAdd(l.n1->id,l.n1);
 all.shkAdd(l.n2->id,l.n2);
 storeCluster(l.n1->cluster);
 storeCluster(l.n2->cluster);
 linkNodes(l);
}



//create a new node with name id1, paired with an existing node n2
bool addNode(CLinkData& p) {
 //use the same cluster
 bool linked;
 if (p.n1==NULL) { //p.n1 is the new node name
   p.n1=new CNode(p.id1, p.len1);

   //p.n1->bestHit=p.n2;
   //p.n1->bestHitScore=p.score;
   all.shkAdd(p.n1->id,p.n1);
   storeCluster(p.n1->cluster);
   linked=linkNodes(p);
   }
  else { //p.n2 is the new node's name
   p.n2=new CNode(p.id2, p.len2);
   //p.n2->bestHit=p.n1;
   //p.n2->bestHitScore=p.score;
   all.shkAdd(p.n2->id,p.n2);
   storeCluster(p.n2->cluster);
   linked=linkNodes(p);
   }
   
 return linked;
}

/* addPair() should slowly build the t-clusters and update
   the links for both n1 and n2
*/

bool seq_filter(char* s1, char* s2) {
 //also checks if any of the sequences are "seed"
 /* if (seedList.Find(s1)==NULL 
          && seedList.Find(s2)==NULL)
     return false; */
 bool isET1, isET2;
 if (ETcheck) {
   isET1=isET(s1);
   isET2=isET(s2);
   if (flt_ET_only) 
     return (isET1 && isET2); //both are et
   else if (flt_EST_only)
     return (!isET1 && !isET2);
   else if (flt_EST2ET)
     return (isET1 || isET2);
   }
 if (flt_seqOnly) {
   return (seqonlyList.hasKey(s1) && seqonlyList.hasKey(s2));
   }
 return true; //passed
 }
 
bool addPair(CLinkData& p) {
 CNode* n1=all.Find(p.id1);
 CNode* n2=all.Find(p.id2); 
 if (n1==NULL) {
   //###### n1 is new
   if (n2==NULL) //n2 is also new 
        newPair(p);
          //this should set d1 and d2 to each new CNode*
      else { //n1 is new, but n2 already has a t-cluster
        p.n2=n2;
        return addNode(p);
        }
   }
  else {
   //###### n1 has a t-cluster
   if (n2==NULL) {//n2 is new
        p.n1=n1;
        return addNode(p);
        }
      else {//n2 also has a t-cluster, merge them
        p.n1=n1;
        p.n2=n2;
        return linkNodes(p);
        }
   }
 return true;
}


#ifdef DBGTRACE
void memlog(int level, const char* msg1, 
         GCluster* C=NULL, 
         const char* msg2=NULL) {
 GStr s;
 if (level>0) {
   s.format("%d>", level);
   s.padR(level+2);
   }
 if (msg1!=NULL) s+=msg1;
 if (C!=NULL) {
     GStr s2;
     s2.format("(%d){%s",C->Count(), C->Get(0)->id);
     s+=s2;
     for (int i=1;i<C->Count();i++)  {
        s2.format(" %s",C->Get(i)->id);
        s+=s2;
        }
     s+="}";
     }
 if (msg2!=NULL) s+=msg2;
 mem.log(s.chars());
 }
#endif  


int readNames(FILE* f, GHash<int>& xhash) {
  int c;
  int count=0;
  char name[256]; 
  int len=0;
  for(;;) {
    c=getc(f);
    if (isspace(c) || c==EOF) {
      if (len>0) {
        name[len]='\0';      
        xhash.Add(name, NULL);
        count++;
        len=0;
        }
      if (c==EOF) break;  
      continue;
      }
    //a non-space
    name[len]=(char) c;
    len++;
    if (len>255) {
      name[len-1]='\0';
      GError("Error reading file: name too long ('%s') !\n",name);      
      }
    }
 return count;   
}

/* record format: 

>NRCL# <numseqs> <refseq> <refseqlen> <flags>
<seq1> <seq2> ... 

Computing seq weight:
 weight = <numseqs>/2 + <seq_en> + ET_Bonus 
*/

void GCluster::addSeq(CNode* n, CLinkData* pl) {
 Add(n);

 if (n->weight==0) {
   n->weight = n->len/10 + (n->is_ET ? ET_bonus : 0);
   //update the sequence weight if not set yet
   }
 if (n->weight>weight) { //switch to a better cluster representative!
   weight=n->weight; // cluster's weight is set to this better sequence
   if (n->isFull) { maxnode=n; isHeavy=true; }
   }
 if (maxnode==NULL) maxnode=n;
 total++;
 if  (!isHeavy && (weight + total/2 > heavy_score)) { 
        /* GMessage("heavy after adding sequence '%s' (weight=%d), cl weight=%d, "
             "total=%d\n", n->id, n->weight, weight, total); */
        isHeavy = true;
        }
}

/*
int readSeeds(FILE* f) {
  char *p, *id, *flags;
  char buf[256];
  int c;
  int len, count=0, scount=0;
  int numseqs;
  //assumes: first line is always a defline 
  //         deflines are always shorter than 255
  //         sequence names are always shorter than 255 chars
  CNode* seed;
  while (!feof(f)) {
   if (fgets(buf,255,f)==NULL) break;
   len=strlen(buf);
   if (buf[len-1]!='\n') GError(ERR_SEED_PARSE, buf); //very long defline ?!?
   buf[len-1]='\0';
   if (buf[0]!='>') GError(ERR_SEED_PARSE,buf);
   p=buf+1;
   p=strpbrk(p,"\t "); //skip NRCL#
   if (p==NULL) GError(ERR_SEED_PARSE,buf);
   while (isspace(*p) && *p!='\0') p++;
   //s to the first nonspace
   if (*p=='\0') GError(ERR_SEED_PARSE,buf);
   numseqs=0;
   numseqs=getNextValue(p, buf); //should return number of seqs within this seed
   if (numseqs<=0) GError(ERR_SEED_PARSE,buf);   
   p++; //p is probably on the first char of the refseq
   id=getNextWord(p,buf);//s must now be the seed seq name
   p++;
   len=getNextValue(p,buf); //length of the seed
   p++;
   flags=getNextWord(p,buf);
   //if (strlen(flags)>4) GError(ERR_SEED_PARSE, buf);
   count++;
   int w=len/10+numseqs/2;
   if (strchr(flags,'G')!=NULL || strchr(flags,'F')!=NULL) w+=ET_bonus; //has an NP:
   seed=new CNode(id, w, len);
   //if (strchr(flags,'F')) seed->isFull=1;
    seedList.shkAdd(seed->id, seed);
   //if (numseqs>1) //it's sure it won't become a singleton
   storeCluster(seed->cluster);
   all.shkAdd(seed->id,seed); //add it to the global hash too
   //--now read all component sequences:
   len=0;
   scount=0;
   while ((c=getc(f))!=EOF) {
    if (isspace(c)) { //word separator found
      if (len>0) {
        buf[len]='\0';
        scount++;
        if (strcmp(seed->id, buf)!=0) { //if it's not the seed itself 
           CNode* n; 
           //GMessage("Adding '%s' in seed '%s'\n", buf, seed->id);
           n=new CNode(seed, buf,200); //add directly as a contained sequence
           all.shkAdd(n->id, n); //add it to the global hash too, just in case
           //anyway, exclude it from further analysis!
           xcludeList.shkAdd(n->id, new int(1));
           }
        len=0;
        }
      if (c=='\n') break;
      continue;
      } //isspace
    //a non-space: append to the current name
    buf[len]=(char) c;
    len++;
    if (len>255) {
      buf[len-1]='\0';
      GError("Error reading seqname: too long ('%s') !\n",buf);
      }
    }
   //GMessage("numseqs=%d, scount=%d\n", numseqs, scount);
   if (scount!=numseqs) GError(ERR_SEED_PARSE, seed->id);
   }
 fclose(f);
 return count;
}

*/

void GCluster::setBestLink(CLinkData& l) {
   if (bestLink==NULL) {
     bestLink = new BLink(l);
     }
    else {
     if (l.n1->cluster==this) { //n1 is always from this cluster
        bestLink->n1=l.n1;
        bestLink->n2=l.n2;
	}else { 
        bestLink->n1=l.n2;
        bestLink->n2=l.n1;
	}
	
     bestLink->score=l.score;
     bestLink->pid=l.pid;
     bestLink->scov=l.scov;
     }
   }
   
bool GCluster::betterLink(CLinkData &l) {
    if (bestLink==NULL) return true;
    if (l.n1->cluster==l.n2->cluster) return false;
    if (bestLink->n1->cluster!=bestLink->n2->cluster) {
      //should we make the move? the score is worse or equal
      if (l.score>bestLink->score) return true; //should not happen, if sorted
      if (l.score==bestLink->score) {
	if (l.pid>bestLink->pid) return true;
        if (l.scov>bestLink->scov) return true;
	int cmp1=strcmp(l.n1->id,bestLink->n1->id);
	int cmp2=strcmp(l.n2->id,bestLink->n2->id);
        if (cmp1>0) return true; //stupid last resort comparison
        if (cmp1==0 && cmp2>0) return true; //stupid last resort comparison
        }
      }
     return false; 
   }

