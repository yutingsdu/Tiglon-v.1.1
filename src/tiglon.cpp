// Some of this program was borrowed from StringTie and modified by Ting Yu.

//#define GFF_DEBUG 1 //debugging guides loading
#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <map>
#include "rlink.h"
#include "tmerge.h"
#ifndef NOTHREADS
#include "GThreads.h"
#endif

//#define GMEMTRACE 1  //debugging mem allocation

#ifdef GMEMTRACE
#include "proc_mem.h"
#endif

#define VERSION "2.1.1"



//#define DEBUGPRINT 1

#ifdef DEBUGPRINT
#define DBGPRINT(x) GMessage(x)
#define DBGPRINT2(a,b) GMessage(a,b)
#define DBGPRINT3(a,b,c) GMessage(a,b,c)
#define DBGPRINT4(a,b,c,d) GMessage(a,b,c,d)
#define DBGPRINT5(a,b,c,d,e) GMessage(a,b,c,d,e)
#else
#define DBGPRINT(x)
#define DBGPRINT2(a,b)
#define DBGPRINT3(a,b,c)
#define DBGPRINT4(a,b,c,d)
#define DBGPRINT5(a,b,c,d,e)
#endif
#define USAGE "estimate_expressiong_level " " usage:\n\
 estimate_expressiong_level input.bam -G *.gtf -e \n\
"

FILE* f_out=NULL;
FILE* c_out=NULL;
//#define B_DEBUG 1
#ifdef B_DEBUG
 FILE* dbg_out=NULL;
#endif

GStr outfname;
GStr out_dir;
GStr tmp_path;
GStr tmpfname;
GStr genefname;
bool guided=false;
bool trim=true;
bool eonly=false; // parameter -e ; for mergeMode includes estimated coverage sum in the merged transcripts 
bool longreads=false;
bool rawreads=false;
bool nomulti=false;
bool enableNames=false;
bool includecov=false;
bool fr_strand=false;
bool rf_strand=false;
//bool complete=true; // set by parameter -i the reference annotation contains partial transcripts
bool retained_intron=false; // set by parameter -i for merge option
bool geneabundance=false;
//bool partialcov=false;
int num_cpus=1;
int mintranscriptlen=200; // minimum length for a transcript to be printed
//int sensitivitylevel=1;
uint junctionsupport=10; // anchor length for junction to be considered well supported <- consider shorter??
int junctionthr=1; // number of reads needed to support a particular junction
float readthr=1;     // read coverage per bundle bp to accept it; // paper uses 3
float COV_YU=1;//YU
bool Yuflag = true;//SFlag
float singlethr=4.75;
uint bundledist=50;  // reads at what distance should be considered part of separate bundles
uint runoffdist=200;
float mcov=1; // fraction of bundle allowed to be covered by multi-hit reads paper uses 1
int allowed_nodes=1000;
//bool adaptive=true; // adaptive read coverage -> depends on the overall gene coverage

int no_xs=0; // number of records without the xs tag

float fpkm_thr=1;
float tpm_thr=1;

// different options of implementation reflected with the next three options
bool includesource=true;
//bool EM=false;
//bool weight=false;

float isofrac=0.01;
bool isunitig=true;
GStr label("TIGLON");
GStr ballgown_dir;

GFastaDb* gfasta=NULL;

GStr guidegff;

bool debugMode=false;
bool verbose=false;
bool ballgown=false;

//int maxReadCov=1000000; //max local read coverage (changed with -s option)
//no more reads will be considered for a bundle if the local coverage exceeds this value
//(each exon is checked for this)

bool forceBAM = false; //useful for stdin (piping alignments into StringTie)

bool mergeMode = false; //--merge option
bool keepTempFiles = false; //--keeptmp

int GeneNo=0; //-- global "gene" counter
double Num_Fragments=0; //global fragment counter (aligned pairs)
double Frag_Len=0;
double Cov_Sum=0;
//bool firstPrint=true; //just for writing the GFF header before the first transcript is printed

GffNames* gseqNames=NULL; //used as a dictionary for genomic sequence names and ids

int refseqCount=0;


#ifdef GMEMTRACE
 double maxMemRS=0;
 double maxMemVM=0;
 GStr maxMemBundle;
#endif


#ifndef NOTHREADS
//single producer, multiple consumers
//main thread/program is always loading the producer
GMutex dataMutex; //manage availability of data records ready to be loaded by main thread
GVec<int> dataClear; //indexes of data bundles cleared for loading by main thread (clear data pool)
GConditionVar haveBundles; //will notify a thread that a bundle was loaded in the ready queue
                           //(or that no more bundles are coming)
int bundleWork=1; // bit 0 set if bundles are still being prepared (BAM file not exhausted yet)
                  // bit 1 set if there are Bundles ready in the queue

//GFastMutex waitMutex;
GMutex waitMutex; // controls threadsWaiting (idle threads counter)

int threadsWaiting; // idle worker threads
GConditionVar haveThreads; //will notify the bundle loader when a thread
                          //is available to process the currently loaded bundle

GConditionVar haveClear; //will notify when bundle buf space available

GMutex queueMutex; //controls bundleQueue and bundles access

GFastMutex printMutex; //for writing the output to file

GFastMutex logMutex; //only when verbose - to avoid mangling the log output

GFastMutex bamReadingMutex;

GFastMutex countMutex;

#endif

GHash<int> excludeGseqs; //hash of chromosomes/contigs to exclude (e.g. chrM)

bool NoMoreBundles=false;
bool moreBundles(); //thread-safe retrieves NoMoreBundles
void noMoreBundles(); //sets NoMoreBundles to true
//--
void processOptions(GArgs& args);
char* sprintTime();

void processBundle(BundleData* bundle);
//void processBundle1stPass(BundleData* bundle); //two-pass testing

void writeUnbundledGuides(GVec<GRefData>& refdata, FILE* fout, FILE* gout=NULL);

#ifndef NOTHREADS

bool noThreadsWaiting();

void workerThread(GThreadData& td); // Thread function

//prepare the next free bundle for loading
int waitForData(BundleData* bundles);
#endif


TInputFiles bamreader;
map<string,vector<pair<int,int> > > ChrTransBorder_Map;
bool check_overlap(int l, int r, vector<pair<int,int> > V)
{
    double R=-1.0;
    double L = 1.0*(r-l+1);
    for(size_t i=0;i<V.size();i++)
    {
	int A = V[i].first, B = V[i].second;
	if(l<=A &&( A<=r && r<=B))
	    R = double(1.0*(r-A+1)/L);
	else if((A<=l && l<=B) &&(A<=r && r<=B))
	    return true;
	else if(l<=A && B<=r)
	    R = double(1.0*(B-A+1)/L);
	else if((A<=l&&l<=B) && B<=r)
	    R = double(1.0*(B-l+1)/L);

	if(R>=0.85) return true;
	else if(R>0) return false;

    }
    return false;
}
int main(int argc, char* argv[]) {

 // == Process arguments.
 GArgs args(argc, argv,
   "debug;help;version;conservative;keeptmp;rseq=;bam;fr;rf;merge;exclude=zEihvteuLRx:n:j:s:D:G:C:S:l:m:o:a:j:c:f:p:g:P:M:Bb:A:F:T:");
 args.printError(USAGE, true);

 processOptions(args);

 GVec<GRefData> refguides; // plain vector with transcripts for each chromosome

 //table indexes for Ballgown Raw Counts data (-B/-b option)
 GPVec<RC_TData> guides_RC_tdata(true); //raw count data or other info for all guide transcripts
 GPVec<RC_Feature> guides_RC_exons(true); //raw count data for all guide exons
 GPVec<RC_Feature> guides_RC_introns(true);//raw count data for all guide introns

 GVec<int> alncounts(30,0); //keep track of the number of read alignments per chromosome [gseq_id]

 int bamcount=bamreader.start(); //setup and open input files
#ifndef GFF_DEBUG
 if (bamcount<1) {
	 GError("%sError: no input files provided!\n",USAGE);
 }
#endif

#ifdef DEBUGPRINT
  verbose=true;
#endif

const char* ERR_BAM_SORT="\nError: the input alignment file is not sorted!\n";

 if(guided) { // read guiding transcripts from input gff file
	 if (verbose) {
		 printTime(stderr);
		 GMessage(" Loading reference annotation (guides)..\n");
	 }
   FILE* f=fopen(guidegff.chars(),"r");
   if (f==NULL) GError("Error: could not open reference annotation file (%s)!\n",
       guidegff.chars());
   //                transcripts_only    sort by location?
   GffReader gffr(f,       true,             true); //loading only recognizable transcript features
   gffr.setRefAlphaSorted(); //alphabetical sorting of refseq IDs
   gffr.showWarnings(verbose);
   //        keepAttrs    mergeCloseExons   noExonAttrs
   gffr.readAll(false,          true,        true);
   //the list of GffObj is in gffr.gflst, sorted by chromosome and start-end coordinates
   //collect them in other data structures, if it's kept for later call gffobj->isUsed(true)
   // (otherwise it'll be deallocated when gffr is destroyed due to going out of scope)
   refseqCount=gffr.gseqtable.Count();
   if (refseqCount==0 || gffr.gflst.Count()==0) {
	   GError("Error: could not any valid reference transcripts in %s (invalid GTF/GFF file?)\n",
			   guidegff.chars());
   }
   refguides.setCount(refseqCount); //maximum gseqid
   uint c_tid=0;
   uint c_exon_id=0;
   uint c_intron_id=0;
   GList<RC_Feature> uexons(true, false, true); //sorted, free items, unique
   GList<RC_Feature> uintrons(true, false, true);
   //assign unique transcript IDs based on the sorted order
   int last_refid=-1;
   bool skipGseq=false;
   for (int i=0;i<gffr.gflst.Count();i++) {
	   GffObj* m=gffr.gflst[i];
	   if (last_refid!=m->gseq_id) {
		   //chromosome switch
		   if (ballgown) { //prepare memory storage/tables for all guides on this chromosome
			   uexons.Clear();
			   uintrons.Clear();
		   }
		   last_refid=m->gseq_id;
		   skipGseq=excludeGseqs.hasKey(m->getGSeqName());
	   }
	   //sanity check: make sure there are no exonless "genes" or other
	   if (skipGseq) continue;
	   if (m->exons.Count()==0) {
		    if (verbose)
		    	GMessage("Warning: exonless GFF %s feature with ID %s found, added implicit exon %d-%d.\n",
		    			m->getFeatureName(), m->getID(), m->start, m->end);
		    m->addExon(m->start, m->end); //should never happen!
	   }
	   //DONE: always keep a RC_TData pointer around, with additional info about guides
	   RC_TData* tdata=new RC_TData(*m, ++c_tid);
	   m->uptr=tdata;
	   guides_RC_tdata.Add(tdata);
	   if (ballgown) { //already gather exon & intron info for all ref transcripts
		   tdata->rc_addFeatures(c_exon_id, uexons, guides_RC_exons,
		          c_intron_id, uintrons, guides_RC_introns);
	   }
	   GRefData& grefdata = refguides[m->gseq_id];
	   grefdata.add(&gffr, m); //transcripts already sorted by location
   }
	 if (verbose) {
		 printTime(stderr);
		 GMessage(" %d reference transcripts loaded.\n", gffr.gflst.Count());
	 }
 }

#ifdef GFF_DEBUG
  for (int r=0;r<refguides.Count();++r) {
	  GRefData& grefdata = refguides[r];
      for (int k=0;k<grefdata.rnas.Count();++k) {
    	  GMessage("#transcript #%d : %s (%d exons)\n", k, grefdata.rnas[k]->getID(), grefdata.rnas[k]->exons.Count());
    	  grefdata.rnas[k]->printGff(stderr);
      }
  }
  GMessage("GFF Debug mode, exiting...\n");
  exit(0);
#endif

 // --- here we do the input processing
 gseqNames=GffObj::names; //might have been populated already by gff data
 gffnames_ref(gseqNames);  //initialize the names collection if not guided


 GHash<int> hashread;      //read_name:pos:hit_index => readlist index

 GList<GffObj>* guides=NULL; //list of transcripts on a specific chromosome

 int currentstart=0, currentend=0;
 int ng_start=0;
 int ng_end=-1;
 int ng=0;
 GStr lastref;
 bool no_ref_used=true;
 int lastref_id=-1; //last seen gseq_id
 // int ncluster=0; used it for debug purposes only

 //Ballgown files
 FILE* f_tdata=NULL;
 FILE* f_edata=NULL;
 FILE* f_idata=NULL;
 FILE* f_e2t=NULL;
 FILE* f_i2t=NULL;
if (ballgown)
 Ballgown_setupFiles(f_tdata, f_edata, f_idata, f_e2t, f_i2t);
#ifndef NOTHREADS
//model: one producer, multiple consumers
#define DEF_TSTACK_SIZE 8388608
 int tstackSize=GThread::defaultStackSize();
 size_t defStackSize=0;
if (tstackSize<DEF_TSTACK_SIZE) defStackSize=DEF_TSTACK_SIZE;
 if (verbose) {
	 if (defStackSize>0){
		 int ssize=defStackSize;
		 GMessage("Default stack size for threads: %d (increased to %d)\n", tstackSize, ssize);
	 }
	 else GMessage("Default stack size for threads: %d\n", tstackSize);
 }
 GThread* threads=new GThread[num_cpus]; //bundle processing threads

 GPVec<BundleData> bundleQueue(false); //queue of loaded bundles
 //the consumers take (pop) bundles out of this queue for processing
 //the producer populates this queue with bundles built from reading the BAM input

 BundleData* bundles=new BundleData[num_cpus+1];
 //bundles[0..num_cpus-1] are processed by threads, loading bundles[num_cpus] first

 dataClear.setCapacity(num_cpus+1);
 for (int b=0;b<num_cpus;b++) {
	 threads[b].kickStart(workerThread, (void*) &bundleQueue, defStackSize);
	 bundles[b+1].idx=b+1;
	 dataClear.Push(b);
   }
 BundleData* bundle = &(bundles[num_cpus]);
#else
 BundleData bundles[1];
 BundleData* bundle = &(bundles[0]);
#endif
 GBamRecord* brec=NULL;
 bool more_alns=true;
 TAlnInfo* tinfo=NULL; // for --merge
 int prev_pos=0;
 bool skipGseq=false;
 while (more_alns) {
	 bool chr_changed=false;
	 int pos=0;
	 const char* refseqName=NULL;
	 char xstrand=0;
	 int nh=1;
	 int hi=0;
	 int gseq_id=lastref_id;  //current chr id
	 bool new_bundle=false;
	 //delete brec;
	 if ((brec=bamreader.next())!=NULL) 
	 {

	 }
	 else { //no more alignments
		 more_alns=false;
		 new_bundle=true; //fake a new start (end of last bundle)
	 }

	 if (new_bundle || chr_changed) {

		 if (chr_changed) {
			 lastref=refseqName;
			 lastref_id=gseq_id;
			 currentend=0;
		 }

		 if (!more_alns) {
			 noMoreBundles();
			 break;
		 }
         }
 } //for each read alignment

 //cleaning up
 delete brec;
 //bamreader.bclose();
 bamreader.stop(); //close all BAM files


 delete gfasta;

#ifndef NOTHREADS
 for (int t=0;t<num_cpus;t++)
	 threads[t].join();
 if (verbose) {
   printTime(stderr);
   GMessage(" All threads finished.\n");
 }
 delete[] threads;
 delete[] bundles;
#else
 if (verbose) {
    printTime(stderr);
    GMessage(" Done.\n");
 }
#endif

#ifdef B_DEBUG
 fclose(dbg_out);
#endif
 if (mergeMode && guided )
	 writeUnbundledGuides(refguides, f_out);

 fclose(f_out);
 if (c_out && c_out!=stdout) fclose(c_out);

 if(verbose && no_xs>0)
	 GMessage("Number spliced alignments missing the XS tag (skipped): %d\n",no_xs);



 if (!keepTempFiles) {
   tmp_path.chomp('/');
   remove(tmp_path);
 }


 gffnames_unref(gseqNames); //deallocate names collection


#ifdef GMEMTRACE
 if(verbose) GMessage(" Max bundle memory: %6.1fMB for bundle %s\n", maxMemRS/1024, maxMemBundle.chars());
#endif
} // -- END main

//----------------------------------------
char* sprintTime() {
	static char sbuf[32];
	time_t ltime; /* calendar time */
	ltime=time(NULL);
	struct tm *t=localtime(&ltime);
	sprintf(sbuf, "%02d_%02d_%02d:%02d:%02d",t->tm_mon+1, t->tm_mday,
			t->tm_hour, t->tm_min, t->tm_sec);
	return(sbuf);
}

void processOptions(GArgs& args) {


	if (args.getOpt('h') || args.getOpt("help")) {
		fprintf(stdout,"%s",USAGE);
	    exit(0);
	}
	if (args.getOpt("version")) {
	   fprintf(stdout,"%s\n",VERSION);
	   exit(0);
	}


	 longreads=(args.getOpt('L')!=NULL);
	 if(longreads) {
		 bundledist=0;
		 singlethr=1.5;
	 }


	if (args.getOpt("conservative")) {
	  isofrac=0.05;
	  singlethr=4.75;
	  readthr=1.5;
	  trim=false;
	}

	if (args.getOpt('t')) {
		trim=false;
	}

	if (args.getOpt("fr")) {
		fr_strand=true;
	}
	if (args.getOpt("rf")) {
		rf_strand=true;
		if(fr_strand) GError("Error: --fr and --rf options are incompatible.\n");
	}

	 debugMode=(args.getOpt("debug")!=NULL || args.getOpt('D')!=NULL);
	 forceBAM=(args.getOpt("bam")!=NULL); //assume the stdin stream is BAM instead of text SAM
	 mergeMode=(args.getOpt("merge")!=NULL);
	 keepTempFiles=(args.getOpt("keeptmp")!=NULL);
	 //adaptive=!(args.getOpt('d')!=NULL);
	 verbose=(args.getOpt('v')!=NULL);
	 if (verbose) {
	     fprintf(stderr, "Running StringTie " VERSION ". Command line:\n");
	     args.printCmdLine(stderr);
	 }
	 //complete=!(args.getOpt('i')!=NULL);
	 // trim=!(args.getOpt('t')!=NULL);
	 includesource=!(args.getOpt('z')!=NULL);
	 //EM=(args.getOpt('y')!=NULL);
	 //weight=(args.getOpt('w')!=NULL);

	 GStr s=args.getOpt('m');
	 if (!s.is_empty()) {
	   mintranscriptlen=s.asInt();
	   if (!mergeMode) {
		   if (mintranscriptlen<30)
			   GError("Error: invalid -m value, must be >=30)\n");
	   }
	   else if (mintranscriptlen<0) GError("Error: invalid -m value, must be >=0)\n");
	 }
	 else if(mergeMode) mintranscriptlen=50;

	 /*
	 if (args.getOpt('S')) {
		 // sensitivitylevel=2; no longer supported from version 1.0.3
		 sensitivitylevel=1;
	 }
	*/

	 s=args.getOpt("rseq");
	 if (s.is_empty())
		 s=args.getOpt('S');
	 if (!s.is_empty()) {
		 gfasta=new GFastaDb(s.chars());
	 }

     s=args.getOpt('x');
     if (!s.is_empty()) {
    	 //split by comma and populate excludeGSeqs
    	 s.startTokenize(" ,\t");
    	 GStr chrname;
    	 while (s.nextToken(chrname)) {
    		 excludeGseqs.Add(chrname.chars(),new int(0));
    	 }
     }

     /*
	 s=args.getOpt('n');
	 if (!s.is_empty()) {
		 sensitivitylevel=s.asInt();
		 if(sensitivitylevel<0) {
			 sensitivitylevel=0;
			 GMessage("sensitivity level out of range: setting sensitivity level at 0\n");
		 }
		 if(sensitivitylevel>3) {
			 sensitivitylevel=3;
			 GMessage("sensitivity level out of range: setting sensitivity level at 2\n");
		 }
	 }
	*/


	 s=args.getOpt('p');
	 if (!s.is_empty()) {
		   num_cpus=s.asInt();
		   if (num_cpus<=0) num_cpus=1;
	 }

	 s=args.getOpt('a');
	 if (!s.is_empty()) {
		 junctionsupport=(uint)s.asInt();
	 }

	 s=args.getOpt('j');
	 if (!s.is_empty()) junctionthr=s.asInt();


	 rawreads=(args.getOpt('R')!=NULL);
	 if(rawreads) {
		 if(!longreads) {
			 if(verbose) GMessage("Enable longreads processing\n");
			 longreads=true;
			 bundledist=0;
		 }
		 readthr=0;
	 }

	 s=args.getOpt('c');
	 if (!s.is_empty()) {
		 COV_YU = (float)s.asDouble();//YU
		 //readthr=(float)s.asDouble();
		 if (readthr<0.001 && !mergeMode) {
			 GError("Error: invalid -c value, must be >=0.001)\n");
		 }
	 }
	 else if(mergeMode) readthr=0;


	 s=args.getOpt('g');
	 if (!s.is_empty()) {
		 bundledist=s.asInt();
		 if(bundledist>runoffdist) runoffdist=bundledist;
	 }
	 else if(mergeMode) bundledist=250; // should figure out here a reasonable parameter for merge

	 s=args.getOpt('F');
	 if (!s.is_empty()) {
		 fpkm_thr=(float)s.asDouble();
	 }
	 //else if(mergeMode) fpkm_thr=0;

	 s=args.getOpt('T');
	 if (!s.is_empty()) {
		 tpm_thr=(float)s.asDouble();
	 }
	 //else if(mergeMode) tpm_thr=0;

	 s=args.getOpt('l');
	 if (!s.is_empty()) label=s;
	 else if(mergeMode) label="TIGLON";//"MSTRG";

	 s=args.getOpt('f');
	 if (!s.is_empty()) {
		 float yuisofrac = (float)s.asDouble();
		 //isofrac=(float)s.asDouble();
		 if(yuisofrac == 1)  Yuflag = true;
		 else Yuflag = false; //YU spk
		 if(isofrac>=1) GError("Miminum isoform fraction (-f coefficient: %f) needs to be less than 1\n",isofrac);
	 }
	 else if(mergeMode) isofrac=0.01;
	 s=args.getOpt('M');
	 if (!s.is_empty()) {
		 mcov=(float)s.asDouble();
	 }

	 genefname=args.getOpt('A');
	 if(!genefname.is_empty()) {
		 geneabundance=true;
	 }

	 tmpfname=args.getOpt('o');
	 //YU
	 /* gingko
	 if(fileExists(tmpfname.chars()))
	 {
	    //cout<<"HH"<<endl;
	 }
	 else Gmkdir(tmpfname.chars());
	 */
	 // coverage saturation no longer used after version 1.0.4; left here for compatibility with previous versions
	 s=args.getOpt('s');
	 if (!s.is_empty()) {
		 singlethr=(float)s.asDouble();
		 if (readthr<0.001 && !mergeMode) {
			 GError("Error: invalid -s value, must be >=0.001)\n");
		 }
	 }

	 if (args.getOpt('G')) {
	   guidegff=args.getOpt('G');
	   if (fileExists(guidegff.chars())>1)
	        guided=true;
	   else GError("Error: reference annotation file (%s) not found.\n",
	             guidegff.chars());
	 }

	 enableNames=(args.getOpt('E')!=NULL);

	 retained_intron=(args.getOpt('i')!=NULL);

	 nomulti=(args.getOpt('u')!=NULL);

	 //isunitig=(args.getOpt('U')!=NULL);

	 eonly=(args.getOpt('e')!=NULL);
	 if(eonly && rawreads) {
		 if(verbose) GMessage("Error: can not use -e and -R at the same time; parameter -e will be ignored\n");
	 }
	 else if(eonly && mergeMode) {
		 eonly=false;
		 includecov=true;
	 }
	 else if(eonly && !guided)
		 GError("Error: invalid -e usage, GFF reference not given (-G option required).\n");


	 ballgown_dir=args.getOpt('b');
	 ballgown=(args.getOpt('B')!=NULL);
	 if (ballgown && !ballgown_dir.is_empty()) {
		 GError("Error: please use either -B or -b <path> options, not both.");
	 }
	 if ((ballgown || !ballgown_dir.is_empty()) && !guided)
		 GError("Error: invalid -B/-b usage, GFF reference not given (-G option required).\n");

	 /* s=args->getOpt('P');
	 if (!s.is_empty()) {
		 if(!guided) GError("Error: option -G with reference annotation file has to be specified.\n");
		 c_out=fopen(s.chars(), "w");
		 if (c_out==NULL) GError("Error creating output file %s\n", s.chars());
		 partialcov=true;
	 }
	 else { */
		 s=args.getOpt('C');
		 if (!s.is_empty()) {
			 if(!guided) GError("Error: invalid -C usage, GFF reference not given (-G option required).\n");
			 c_out=fopen(s.chars(), "w");
			 if (c_out==NULL) GError("Error creating output file %s\n", s.chars());
		 }
	 //}
	int numbam=args.startNonOpt();
#ifndef GFF_DEBUG
	if (numbam < 1 ) {
	 	 GMessage("%s\nError: no input file provided!\n",USAGE);
	 	 exit(1);
	}
#endif
	const char* ifn=NULL;
	while ( (ifn=args.nextNonOpt())!=NULL) {
		//input alignment files
		bamreader.Add(ifn);
	}
	//deferred creation of output path
	outfname="stdout";
	out_dir="./";
	 if (!tmpfname.is_empty() && tmpfname!="-") {
		 if (tmpfname[0]=='.' && tmpfname[1]=='/')
			 tmpfname.cut(0,2);
		 outfname=tmpfname;
		 int pidx=outfname.rindex('/');
		 if (pidx>=0) {//path given
			 out_dir=outfname.substr(0,pidx+1);
			 tmpfname=outfname.substr(pidx+1);
		 }
	 }
	 else { // stdout
		tmpfname=outfname;
		char *stime=sprintTime();
		tmpfname.tr(":","-");
		tmpfname+='.';
		tmpfname+=stime;
	 }
	 if (out_dir!="./") {
		 if (fileExists(out_dir.chars())==0) {
			//directory does not exist, create it
			//YU
			if (Gmkdir(out_dir.chars()) && !fileExists(out_dir.chars())) {
				GError("Error: cannot create directory %s!\n", out_dir.chars());
			}
		 }
	 }
	 if (!genefname.is_empty()) {
		 if (genefname[0]=='.' && genefname[1]=='/')
		 			 genefname.cut(0,2);
	 //attempt to create the gene abundance path
		 GStr genefdir("./");
		 int pidx=genefname.rindex('/');
		 if (pidx>=0) { //get the path part
				 genefdir=genefname.substr(0,pidx+1);
				 //genefname=genefname.substr(pidx+1);
		 }
		 if (genefdir!="./") {
			 if (fileExists(genefdir.chars())==0) {
				//directory does not exist, create it
				if (Gmkdir(genefdir.chars()) || !fileExists(genefdir.chars())) {
					GError("Error: cannot create directory %s!\n", genefdir.chars());
				}
			 }
		 }

	 }
	 { //prepare temp path
		 GStr stempl(out_dir);
		 stempl.chomp('/');
		 stempl+="/tmp.XXXXXXXX";
		 char* ctempl=Gstrdup(stempl.chars());
	     //Gmktempdir(ctempl);//YU //
	     tmp_path=ctempl;
	     tmp_path+='/';
	     GFREE(ctempl);
	 }
	
	 //tmpfname=tmp_path+tmpfname;//YU
	 if (ballgown) ballgown_dir=out_dir;
	   else if (!ballgown_dir.is_empty()) {
			ballgown=true;
			ballgown_dir.chomp('/');ballgown_dir+='/';
			if (fileExists(ballgown_dir.chars())==0) {
				//directory does not exist, create it
				if (Gmkdir(ballgown_dir.chars()) && !fileExists(ballgown_dir.chars())) {
					GError("Error: cannot create directory %s!\n", ballgown_dir.chars());
				}

			}
	   }
#ifdef B_DEBUG
	 GStr dbgfname(tmpfname);
	 dbgfname+=".dbg";
	 dbg_out=fopen(dbgfname.chars(), "w");
	 if (dbg_out==NULL) GError("Error creating debug output file %s\n", dbgfname.chars());
#endif

	 if(mergeMode) {
/*//YU
		 f_out=stdout;
		 if(outfname!="stdout") {
			 f_out=fopen(outfname.chars(), "w");
			 if (f_out==NULL) GError("Error creating output file %s\n", outfname.chars());
		 }
		 fprintf(f_out,"# ");
		 args.printCmdLine(f_out);
		 fprintf(f_out,"# StringTie version %s\n",VERSION);
*/
	 }
	 else {
		 //tmpfname = "transborrow" + tmpfname;
		 tmpfname=outfname;//gingko
		 //cout<<"HH2: "<<outfname<<" "<<tmpfname<<endl;
		 outfname += ".gtf";
		 char *stime=sprintTime();
                 //tmpfname.tr(":","-");
                 tmpfname+=".";
                 tmpfname+=stime;
		 tmpfname+=".tmp.gtf";
		 //cout<<tmpfname<<endl;
		// tmpfname = "yu.gtf";//YU
		 f_out=fopen(tmpfname.chars(), "w");
		 if (f_out==NULL) GError("Error creating output file %s\n", tmpfname.chars());
	 }
}

//---------------
bool moreBundles() { //getter (interogation)
	bool v=true;
#ifndef NOTHREADS
  GLockGuard<GFastMutex> lock(bamReadingMutex);
#endif
  v = ! NoMoreBundles;
  return v;
}

void noMoreBundles() {
#ifndef NOTHREADS
		bamReadingMutex.lock();
		NoMoreBundles=true;
		bamReadingMutex.unlock();
		queueMutex.lock();
		bundleWork &= ~(int)0x01; //clear bit 0;
		queueMutex.unlock();
		bool areThreadsWaiting=true;
		do {
		  waitMutex.lock();
		   areThreadsWaiting=(threadsWaiting>0);
		  waitMutex.unlock();
		  if (areThreadsWaiting) {
		    DBGPRINT("##> NOTIFY ALL workers: no more data!\n");
		    haveBundles.notify_all();
		    current_thread::sleep_for(10);
		    waitMutex.lock();
		     areThreadsWaiting=(threadsWaiting>0);
		    waitMutex.unlock();
		    current_thread::sleep_for(10);
		  }
		} while (areThreadsWaiting); //paranoid check that all threads stopped waiting
#else
	  NoMoreBundles=true;
#endif
}

void processBundle(BundleData* bundle) {
	if (verbose) {
	#ifndef NOTHREADS
		GLockGuard<GFastMutex> lock(logMutex);
	#endif
		printTime(stderr);
		GMessage(">bundle %s:%d-%d [%lu alignments (%d distinct), %d junctions, %d guides] begins processing...\n",
				bundle->refseq.chars(), bundle->start, bundle->end, bundle->numreads, bundle->readlist.Count(), bundle->junction.Count(),
                bundle->keepguides.Count());
	#ifdef GMEMTRACE
			double vm,rsm;
			get_mem_usage(vm, rsm);
			GMessage("\t\tstart memory usage: %6.1fMB\n",rsm/1024);
			if (rsm>maxMemRS) {
				maxMemRS=rsm;
				maxMemVM=vm;
				maxMemBundle.format("%s:%d-%d(%d)", bundle->refseq.chars(), bundle->start, bundle->end, bundle->readlist.Count());
			}
	#endif
	}
#ifdef B_DEBUG
	for (int i=0;i<bundle->keepguides.Count();++i) {
		GffObj& t=*(bundle->keepguides[i]);
		RC_TData* tdata=(RC_TData*)(t.uptr);
		fprintf(dbg_out, ">%s (t_id=%d) %s%c %d %d\n", t.getID(), tdata->t_id, t.getGSeqName(), t.strand, t.start, t.end );
		for (int fe=0;fe < tdata->t_exons.Count(); ++fe) {
			RC_Feature& exoninfo = *(tdata->t_exons[fe]);
			fprintf(dbg_out, "%d\texon\t%d\t%d\t%c\t%d\t%d\n", exoninfo.id, exoninfo.l, exoninfo.r,
					    exoninfo.strand, exoninfo.rcount, exoninfo.ucount);
			if (! (exoninfo==*(bundle->rc_data->guides_RC_exons->Get(exoninfo.id-1))))
				 GError("exoninfo with id (%d) not matching!\n", exoninfo.id);
		}
		for (int fi=0;fi < tdata->t_introns.Count(); ++fi) {
			RC_Feature& introninfo = *(tdata->t_introns[fi]);
			fprintf(dbg_out, "%d\tintron\t%d\t%d\t%c\t%d\t%d\n", introninfo.id, introninfo.l, introninfo.r,
					introninfo.strand, introninfo.rcount, introninfo.ucount);
			if (! (introninfo==*(bundle->rc_data->guides_RC_introns->Get(introninfo.id-1))))
				 GError("introninfo with id (%d) not matching!\n", introninfo.id);
		}
		//check that IDs are properly assigned
		if (tdata->t_id!=(uint)t.udata) GError("tdata->t_id(%d) not matching t.udata(%d)!\n",tdata->t_id, t.udata);
		if (tdata->t_id!=bundle->rc_data->guides_RC_tdata->Get(tdata->t_id-1)->t_id)
			 GError("tdata->t_id(%d) not matching rc_data[t_id-1]->t_id (%d)\n", tdata->t_id, bundle->rc_data->g_tdata[tdata->t_id-1]->t_id);

	}
#endif
	infer_transcripts(bundle);

	if (ballgown && bundle->rc_data) {
		rc_update_exons(*(bundle->rc_data));
	}
	if (bundle->pred.Count()>0 || ((eonly || geneabundance) && bundle->keepguides.Count()>0)) {
#ifndef NOTHREADS
		GLockGuard<GFastMutex> lock(printMutex);
#endif
		if(mergeMode) GeneNo=printMergeResults(bundle, GeneNo,bundle->refseq);
		else GeneNo=printResults(bundle, GeneNo, bundle->refseq);
	}

	if (bundle->num_fragments) {
		#ifndef NOTHREADS
				GLockGuard<GFastMutex> lock(countMutex);
		#endif
		Num_Fragments+=bundle->num_fragments;
		Frag_Len+=bundle->frag_len;
		Cov_Sum+=bundle->sum_cov;
	}

	if (verbose) {
		#ifndef NOTHREADS
				GLockGuard<GFastMutex> lock(logMutex);
		#endif
	  /*
	  SumReads+=bundle->sumreads;
	  SumFrag+=bundle->sumfrag;
	  NumCov+=bundle->num_cov;
	  NumReads+=bundle->num_reads;
	  NumFrag+=bundle->num_frag;
	  NumFrag3+=bundle->num_fragments3;
	  SumFrag3+=bundle->sum_fragments3;
	  fprintf(stderr,"Number of fragments in bundle: %g with length %g\n",bundle->num_fragments,bundle->frag_len);
	  fprintf(stderr,"Number of fragments in bundle: %g with sum %g\n",bundle->num_fragments,bundle->frag_len);
	  */
	  printTime(stderr);
	  GMessage("^bundle %s:%d-%d done (%d processed potential transcripts).\n",bundle->refseq.chars(),
	  		bundle->start, bundle->end, bundle->pred.Count());
	#ifdef GMEMTRACE
		    double vm,rsm;
		    get_mem_usage(vm, rsm);
		    GMessage("\t\tfinal memory usage: %6.1fMB\n",rsm/1024);
		    if (rsm>maxMemRS) {
			    maxMemRS=rsm;
			    maxMemVM=vm;
			    maxMemBundle.format("%s:%d-%d(%d)", bundle->refseq.chars(), bundle->start, bundle->end, bundle->readlist.Count());
		    }
	#endif
	    }
	bundle->Clear();
}

#ifndef NOTHREADS

bool noThreadsWaiting() {
	waitMutex.lock();
	int v=threadsWaiting;
	waitMutex.unlock();
	return (v<1);
}

void workerThread(GThreadData& td) {
	GPVec<BundleData>* bundleQueue = (GPVec<BundleData>*)td.udata;
	//wait for a ready bundle in the queue, until there is no hope for incoming bundles
	DBGPRINT2("---->> Thread%d starting..\n",td.thread->get_id());
	DBGPRINT2("---->> Thread%d locking queueMutex..\n",td.thread->get_id());
	queueMutex.lock(); //enter wait-for-notification loop
	while (bundleWork) {
		DBGPRINT3("---->> Thread%d: waiting.. (queue len=%d)\n",td.thread->get_id(), bundleQueue->Count());
		waitMutex.lock();
		 threadsWaiting++;
		queueMutex.unlock();
		waitMutex.unlock();
		haveThreads.notify_one(); //in case main thread is waiting
		current_thread::yield();
		queueMutex.lock();
		while (bundleWork && bundleQueue->Count()==0) {
		    haveBundles.wait(queueMutex);//unlocks queueMutex and wait until notified
		               //when notified, locks queueMutex and resume
		}
		waitMutex.lock();
		if (threadsWaiting>0) threadsWaiting--;
		waitMutex.unlock();
		DBGPRINT3("---->> Thread%d: awakened! (queue len=%d)\n",td.thread->get_id(),bundleQueue->Count());
		BundleData* readyBundle=NULL;
		if ((bundleWork & 0x02)!=0 && (readyBundle=bundleQueue->Pop())!=NULL) { //is bit 1 set?
				if (bundleQueue->Count()==0)
					 bundleWork &= ~(int)0x02; //clear bit 1 (queue is empty)
				//Num_Fragments+=readyBundle->num_fragments;
				//Frag_Len+=readyBundle->frag_len;
				queueMutex.unlock();
				processBundle(readyBundle);
				DBGPRINT3("---->> Thread%d processed bundle #%d, now locking back dataMutex and queueMutex\n",
						td.thread->get_id(), readyBundle->idx);
				dataMutex.lock();
				dataClear.Push(readyBundle->idx);
				DBGPRINT3("---->> Thread%d pushed bundle #%d into dataClear",
										td.thread->get_id(), readyBundle->idx);
				dataMutex.unlock();
				DBGPRINT2("---->> Thread%d informing main thread and yielding", td.thread->get_id());
				haveClear.notify_one(); //inform main thread
				current_thread::yield();
				DBGPRINT2("---->> Thread%d processed bundle, now locking back queueMutex\n", td.thread->get_id());
				queueMutex.lock();
				DBGPRINT2("---->> Thread%d locked back queueMutex\n", td.thread->get_id());

		}
		//haveThreads.notify_one();
	} //while there is reason to live
	queueMutex.unlock();
	DBGPRINT2("---->> Thread%d DONE.\n", td.thread->get_id());
}

//prepare the next available bundle slot for loading
int waitForData(BundleData* bundles) {
	int bidx=-1;
	dataMutex.lock();
	DBGPRINT("  #waitForData: locking dataMutex");
	while (dataClear.Count()==0) {
		DBGPRINT("  #waitForData: dataClear.Count is 0, waiting for dataMutex");
		haveClear.wait(dataMutex);
	}
	bidx=dataClear.Pop();
	if (bidx>=0) {
	  bundles[bidx].status=BUNDLE_STATUS_LOADING;
	}

	DBGPRINT("  #waitForData: unlocking dataMutex");
	dataMutex.unlock();
	return bidx;
}

#endif

void writeUnbundledGenes(GHash<CGene>& geneabs, const char* refseq, FILE* gout) {
				 //write unbundled genes from this chromosome
	geneabs.startIterate();
	while (CGene* g=geneabs.NextData()) {
		const char* geneID=g->geneID;
		const char* geneName=g->geneName;
		if (geneID==NULL) geneID=".";
		if (geneName==NULL) geneName=".";
	    fprintf(gout, "%s\t%s\t%s\t%c\t%d\t%d\t0.0\t0.0\t0.0\n",
	    		geneID, geneName, refseq,
	    		g->strand, g->start, g->end);
	}
	geneabs.Clear();
}

void writeUnbundledGuides(GVec<GRefData>& refdata, FILE* fout, FILE* gout) {
 for (int g=0;g<refdata.Count();++g) {
	 GRefData& crefd=refdata[g];
	 if (crefd.rnas.Count()==0) continue;
	 GHash<CGene> geneabs;
	 //gene_id abundances (0), accumulating coords
	 for (int m=0;m<crefd.rnas.Count();++m) {
		 GffObj &t = *crefd.rnas[m];
		 RC_TData &td = *(RC_TData*) (t.uptr);
		 if (td.in_bundle) {
			 if (gout && m==crefd.rnas.Count()-1)
			 		writeUnbundledGenes(geneabs, crefd.gseq_name, gout);
			 continue;
		 }
		 //write these guides to output
		 //for --merge and -e
		 if (mergeMode || eonly) {
			  fprintf(fout, "%s\t%s\ttranscript\t%d\t%d\t.\t%c\t.\t",
					  crefd.gseq_name, t.getTrackName(), t.start, t.end, t.strand);
			  if (t.getGeneID())
				  fprintf(fout, "gene_id \"%s\"; ", t.getGeneID());
			  fprintf(fout, "transcript_id \"%s\";",t.getID());
			  if (eonly) {
				if (t.getGeneName())
					  fprintf(fout, " gene_name \"%s\";", t.getGeneName());
			    fprintf(fout, " cov \"0.0\"; FPKM \"0.0\"; TPM \"0.0\";");
			  }
			  else { //merge_mode
				  if (t.getGeneName())
					  fprintf(fout, " gene_name \"%s\";", t.getGeneName());
				  if (t.getGeneID())
					  fprintf(fout, " ref_gene_id \"%s\";", t.getGeneID());
			  }
			  fprintf(fout, "\n");
			  for (int e=0;e<t.exons.Count();++e) {
				  fprintf(fout,"%s\t%s\texon\t%d\t%d\t.\t%c\t.\t",
						  crefd.gseq_name, t.getTrackName(), t.exons[e]->start, t.exons[e]->end, t.strand);
				  if (t.getGeneID())
					  fprintf(fout, "gene_id \"%s\"; ",  t.getGeneID());
				  fprintf(fout,"transcript_id \"%s\"; exon_number \"%d\";",
						  t.getID(), e+1);

				  if (t.getGeneName())
						  fprintf(fout, " gene_name \"%s\";", t.getGeneName());
				  if (eonly) {
					  fprintf(fout, " cov \"0.0\";");
				  }
				  fprintf(fout, "\n");
			  }
		 }
		 if (gout!=NULL) {
			 //gather coords for this gene_id
			 char* geneid=t.getGeneID();
			 if (geneid==NULL) geneid=t.getGeneName();
			 if (geneid!=NULL) {
				 CGene* gloc=geneabs.Find(geneid);
				 if (gloc) {
					 if (gloc->strand!=t.strand)
						 GMessage("Warning: gene \"%s\" (on %s) has reference transcripts on both strands?\n",
								  geneid, crefd.gseq_name);
					 if (t.start<gloc->start) gloc->start=t.start;
					 if (t.end>gloc->end) gloc->end=t.end;
				 } else {
					 //add new geneid locus
					 geneabs.Add(geneid, new CGene(t.start, t.end, t.strand, t.getGeneID(), t.getGeneName()));
				 }
			 }
			 if (m==crefd.rnas.Count()-1)
				  writeUnbundledGenes(geneabs, crefd.gseq_name, gout);
		 } //if geneabundance
	 }
 }
}




