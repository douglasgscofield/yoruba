// yoruba_seda.cpp  (c) Douglas G. Scofield, Dept. Plant Physiology, Umeå University
//
// Seda finds duplicated read pairs in BAM files.
//
// Seda is (*I think*) the Yoruba (Nigeria) verb for 'to copy'.
//
// Uses BamTools C++ API for reading BAM files


// CHANGELOG
//
//
//
// TODO
// --- implement --as-single-end
// --- deal with ordering issue and pairs
// --- double-check the unseen-mates removed issue, make sure it is consistent
// --- make dupMap a class
// --- compare against picard MarkDuplicates and samtools rmdup
// xxx add sorted check, abort if detected not coordinate sorted
// xxx implement the --{single,paired}-end-only options
// xxx incorporate read group into duplicate check
// xxx double-check the better mapping quality solution


#include "yoruba_seda.h"

using namespace std;
using namespace BamTools;
using namespace yoruba;

static string       input_file;         // set from command line
static string       output_file;        // defaults to stdout, set with -o FILE
enum detect_t { DETECT_as_single, DETECT_paired_only, DETECT_single_only, DETECT_all };
static detect_t     opt_detect = DETECT_all;
static bool         opt_remove;         // set with --remove
static bool         opt_duplicatefile;  // set with --duplicate-file FILE
static string       duplicate_file;     // set with --duplicate-file FILE, holds FILE
#ifdef _WITH_DEBUG
static int32_t      opt_debug = 1;
static int64_t      opt_reads = -1;
static int64_t      opt_progress = 100000; // 100000;
static int64_t      last_n_reads_mod = 0;  // helps with progress output during pass1
#endif
static const string delim = "'";
static const string sep = "\t";
static const string endline = "\n";


//-------------------------------------


#ifdef _STANDALONE
int 
main(int argc, char* argv[]) {
    return main_seda(argc, argv);
}
#endif


//-------------------------------------


static int
usage()
{
    cerr << endl;
    cerr << "Usage:   " << YORUBA_NAME << " duplicate [options] <in.bam>" << endl;
    cerr << "         " << YORUBA_NAME << " seda      [options] <in.bam>" << endl;
    cerr << endl;
    cerr << "Either command invokes this function." << endl;
    cerr << endl;
    cerr << "\
Determines duplicate reads in a BAM file, marks them as duplicates, and removes\n\
them on option.\n\
\n\
NOTE: THIS COMMAND IS INCOMPLETE, DO NOT USE\n\
\n\
Options: --as-single-end                     all reads treated as single-end, ignore pairing\n\
         --single-end-only                   only look for duplicates in single-end reads\n\
         --paired-end-only                   only look for duplicates in paired-end reads\n\
         --remove                            remove reads from the output BAM\n\
         --duplicate-file FILE               write duplicate reads to BAM file FILE,\n\
                                             note this does not currently imply --remove\n\
         --o FILE | -o FILE | --output FILE  output file name [default is stdout]\n\
         --? | -? | --help                   longer help\n\
\n";
#ifdef _WITH_DEBUG
    cerr << "\
         --debug INT     debug info level INT [" << opt_debug << "]\n\
         --reads INT     only process INT reads [" << opt_reads << "]\n\
         --progress INT  print reads processed mod INT [" << opt_progress << "]\n\
\n";
#endif
    cerr << "Seda is the Yoruba (Nigeria) verb for 'to copy'." << endl;
    cerr << endl;

    return EXIT_FAILURE;
}


//-------------------------------------


typedef list<BamAlignment>            alignmentList;
typedef alignmentList::iterator       alignmentListI;
typedef alignmentList::const_iterator alignmentListCI;

enum dup_t { dupMap_singleend = -1, dupMap_UNSET = 0, 
    dupMap_paired_one = 1, dupMap_paired_both = 2 };
//typedef map<string, dup_t>            dupMap;
typedef std::tr1::unordered_map<string, dup_t>  dupMap;
typedef dupMap::iterator              dupMapI;
typedef dupMap::const_iterator        dupMapCI;

// local functions, will soon create class out of dupMap
static void listAlignments(const alignmentList& al_set);
static bool isDuplicate(const BamAlignment& al_i, const BamAlignment& al_j);
static void diagnoseDuplicate(const BamAlignment& al_i, const BamAlignment& al_j);
static void determineDuplicates(alignmentList& al_set, alignmentList& al_dups);
static void dump_dupMap(const dupMap& this_dm);
static void update_dupMap(alignmentList& al_set, dupMap& this_dm);
static void query_dupMap(const dupMap& this_dm);
static void clear_dupMap(dupMap& this_dm);
static void clear_dupMap(dupMap& this_dm, dup_t val);
// a member for the class, to get stats for last update of the dupMap
// static void querylastupdate_dupMap(const dupMap& this_dm);

//-------------------------------------


int 
yoruba::main_seda(int argc, char* argv[])
{
    SamProgram new_program;

    new_program.ID = YORUBA_NAME;
    new_program.ID = new_program.ID + " " + argv[0];
    new_program.Name = YORUBA_NAME;
    new_program.Version = YORUBA_VERSION;
    new_program.CommandLine = YORUBA_NAME;
    for (int i = 0; i < argc; ++i)
        new_program.CommandLine = new_program.CommandLine + " " + argv[i];

    //----------------- Command-line options

	if( argc < 2 ) {
		return usage();
	}
    
    enum { OPT_output, OPT_as_single, OPT_single_only, OPT_paired_only,
        OPT_remove, OPT_duplicatefile,
#ifdef _WITH_DEBUG
        OPT_debug, OPT_reads, OPT_progress,
#endif
        OPT_help };

    CSimpleOpt::SOption seda_options[] = {
        { OPT_as_single,       "--as-single-end",   SO_NONE },
        { OPT_single_only,     "--single-end-only", SO_NONE },
        { OPT_paired_only,     "--paired-end-only", SO_NONE },
        { OPT_remove,          "--remove",          SO_NONE },
        { OPT_duplicatefile,   "--duplicate-file",  SO_REQ_SEP },
        { OPT_help,            "--?",               SO_NONE }, 
        { OPT_help,            "-?",                SO_NONE }, 
        { OPT_help,            "--help",            SO_NONE },
        { OPT_output,          "--o",               SO_REQ_SEP },
        { OPT_output,          "-o",                SO_REQ_SEP },
        { OPT_output,          "--output",          SO_REQ_SEP },
#ifdef _WITH_DEBUG
        { OPT_debug,           "--debug",           SO_REQ_SEP },
        { OPT_reads,           "--reads",           SO_REQ_SEP },
        { OPT_progress,        "--progress",        SO_REQ_SEP },
#endif
        SO_END_OF_OPTIONS
    };

    CSimpleOpt args(argc, argv, seda_options);

    while (args.Next()) {
        if (args.LastError() != SO_SUCCESS) {
            cerr << NAME << " invalid argument '" << args.OptionText() << "'" << endl;
            return usage();
        }
        if (args.OptionId() == OPT_help) {
            return usage();
        } else if (args.OptionId() == OPT_as_single) {
            opt_detect = DETECT_as_single;
        } else if (args.OptionId() == OPT_single_only) {
            opt_detect = DETECT_single_only;
        } else if (args.OptionId() == OPT_paired_only) {
            opt_detect = DETECT_paired_only;
        } else if (args.OptionId() == OPT_remove) {
            opt_remove = true;
        } else if (args.OptionId() == OPT_duplicatefile) {
            opt_duplicatefile = true; duplicate_file = args.OptionArg();
        } else if (args.OptionId() == OPT_output) {
            output_file = args.OptionArg();
#ifdef _WITH_DEBUG
        } else if (args.OptionId() == OPT_debug) {
            opt_debug = args.OptionArg() ? atoi(args.OptionArg()) : opt_debug;
        } else if (args.OptionId() == OPT_reads) {
            opt_reads = strtoll(args.OptionArg(), NULL, 10);
        } else if (args.OptionId() == OPT_progress) {
            opt_progress = args.OptionArg() ? strtoll(args.OptionArg(), NULL, 10) : opt_progress;
#endif  // end debug options
        } else {
            cerr << NAME << " unprocessed argument '" << args.OptionText() << "'" << endl;
            return EXIT_FAILURE;
        }
    }

    if (args.FileCount() > 1) {
        cerr << NAME << " requires at most one BAM file specified as input" << endl;
        return usage();
    } else if (args.FileCount() == 1) {
        input_file = args.File(0);
    } else if (input_file.empty()) {
        cerr << NAME << " can't currently read from stdin, ask Doug about it" << endl;
        return EXIT_FAILURE;
        input_file = "/dev/stdin";
    }

    // set up output; if file not specified, use stdout or its equivalent
    if (output_file.empty())
        output_file = "/dev/stdout";


    // in this map, key is a read name, value is true if read paired, false if
    // single-end 
    //
    // on pass 1, if a read name is in this map, then it is either a known
    // duplicate (single-end with val = false) or it is an unmated read of a
    // pair that might be a duplicate (val = true) or is a mated read that is
    // an established duplicate.  if it was an unmated read of a pair its name
    // is removed if the mates are not themselves duplicates.
    //
    // on pass 2, if a read name is in this map, then it is a known duplicate,
    // and is either single-end (val = false) or paired-end (val = true)
    //
    // this map might get very, very large... maybe i'm better off reading in
    // all reads for say 10kbp, and then slurp in steps of 5kbp keep the
    // this-contig reads out of the map entirely
    //
    // also check whether string.c_str() is a faster key than string itself

    //----------------- Open files, start reading data

	BamReader reader;

	if (! reader.Open(input_file)) {
        cerr << NAME << "could not open BAM input" << endl;
        return EXIT_FAILURE;
    }

#ifdef _BAMTOOLS_EXTENSION
    const SamHeader& header = reader.GetConstSamHeader();
#else
    SamHeader header = reader.GetHeader();
#endif

    BamWriter writer;
    BamWriter writer_dups;

    if (! writer.Open(output_file, header, reader.GetReferenceData())) {
        cerr << NAME << " could not open output " << output_file << endl;
        return EXIT_FAILURE;
    }

    if (opt_duplicatefile && ! writer_dups.Open(duplicate_file, header, reader.GetReferenceData())) {
        cerr << NAME << " could not open duplicate output file  " << duplicate_file << endl;
        return EXIT_FAILURE;
    }


    //----------------- Pass 1: Determine which reads are duplicates


    dupMap dup_map;

    int64_t n_reads = 0;  // number of reads processed
    int64_t n_reads_pass1 = 0;  // number of reads processed
    int64_t n_reads_written_to_output = 0;  // number of reads written to output BAM
    int64_t n_reads_written_to_dups = 0;  // number of reads written to duplicates BAM
    int64_t n_reads_removed = 0;  // number of reads removed from output BAM

	BamAlignment al;  // holds the current read from the BAM file

    // it might be more efficient to manage my own list, std::list<> is fine for now
    alignmentList al_set;

    int32_t last_RefID = -2;
    int32_t last_Position = -1;

    if (reader.GetNextAlignment(al)) {
        al_set.push_back(al);
        last_RefID = al.RefID;
        last_Position = al.Position;
        ++n_reads;
        IF_DEBUG(3) cerr << "beginning with " << al_set.size() << " alignments, al.RefID = " 
            << al.RefID << " al.Position = " << al.Position << endl;
    }

	while (! al_set.empty() && (opt_reads < 0 || n_reads < opt_reads)) {

        IF_DEBUG(3) cerr << al_set.size() << " alignments at start of alignment-reading loop" << endl;

        bool al_remaining;

        while ((al_remaining = reader.GetNextAlignment(al)) 
                && al.RefID == last_RefID 
                && al.Position == last_Position ) {
            al_set.push_back(al);
            IF_DEBUG(3) cerr << al_set.size() << " alignments, al.RefID = " << al.RefID 
                    << " al.Position = " << al.Position << endl;
            ++n_reads;
        }

        if (al_remaining && ! isCoordinateSorted(al.RefID, al.Position, last_RefID, last_Position)) {
            cerr << NAME << " input is not coordinate-sorted, " << al.Name << " out of position" << endl;
            return EXIT_FAILURE;
        }

        // all alignments in al_set share RefID and Position

        IF_DEBUG(2) cerr << "read " << al_set.size() << " alignments at Ref = " << last_RefID 
                << " Pos = " << last_Position << endl;

        if (al_set.size() > 1) {
            alignmentList al_dups;  // holds duplicates detected

            IF_DEBUG(2) listAlignments(al_set);
            determineDuplicates(al_set, al_dups);  // which reads here are potential duplicates?
            assert(al_set.empty());  // still true?
            update_dupMap(al_dups, dup_map);  // add duplicates to set for pass 2
            assert(al_dups.empty());

        } else {

            al_set.clear();  // just one read here, no duplicates

        }

        if (al_remaining) {
            al_set.push_back(al);
            last_RefID = al.RefID;
            last_Position = al.Position;
            ++n_reads;
        }
        
        // because we eat reads in chunks, we rarely hit n_reads % opt_progress == 0
        if ((opt_progress || DEBUG(1)) && (n_reads % opt_progress <= last_n_reads_mod))
            cerr << NAME << "[pass1] " << n_reads << " reads examined"
                << ", last at Ref = " << last_RefID << " Pos = " << last_Position
                << ", size of dupMap = " << dup_map.size()
                << endl;
        last_n_reads_mod = n_reads % opt_progress;
	}

    if (opt_progress || DEBUG(1)) {
        cerr << NAME << "[pass1] " << n_reads << " reads examined"
            << ", last at Ref = " << last_RefID << " Pos = " << last_Position
            << ", size of dupMap = " << dup_map.size()
            << endl;
    }

    { // clean the map: remove PE reads with unseen mates
        size_t initial_size = dup_map.size();
        clear_dupMap(dup_map, dupMap_paired_one);
        size_t n_removed = initial_size - dup_map.size();
        if (n_removed || DEBUG(1))
            cerr << NAME << "[pass1] map size was " << initial_size 
                << ", removed " << n_removed << " PE reads with unseen mates, size now is " 
                << dup_map.size() << endl;
    }

    n_reads_pass1 = n_reads;


    //----------------- Pass 2: dup_map holds names of duplicate reads


    IF_DEBUG(1) {
        cerr << NAME << "[pass2] ";
        query_dupMap(dup_map);
    }

    n_reads = 0;
    int64_t n_dupMap_entries_decremented = 0;
    int64_t n_dupMap_entries_erased_SE = 0;
    int64_t n_dupMap_entries_erased_PE = 0;

    reader.Rewind();

    // we nead the read name, so GetNextAlignmentCore() is insufficient
	while (reader.GetNextAlignment(al) && (opt_reads < 0 || n_reads < opt_reads)) {

        ++n_reads;

        dupMapI dupI = dup_map.find(al.Name);

        if (dupI == dup_map.end()) {  // we did not find this read name in dup_map
            
            al.SetIsDuplicate(false);

            writer.SaveAlignment(al);
            ++n_reads_written_to_output;

        } else {  // read name found in dup_map

            al.SetIsDuplicate(true);

            if (opt_duplicatefile) {
                writer_dups.SaveAlignment(al);
                ++n_reads_written_to_dups;
            }

            if (opt_remove) {
                ++n_reads_removed;
            } else {
                writer.SaveAlignment(al);
                ++n_reads_written_to_output;
            }

            if (dupI->second == dupMap_singleend) {
                dup_map.erase(dupI);
                ++n_dupMap_entries_erased_SE;
            } else if (dupI->second == dupMap_paired_one) {  // second of pair
                dup_map.erase(dupI);
                ++n_dupMap_entries_erased_PE;
            } else if (dupI->second == dupMap_paired_both) {
                dupI->second = dupMap_paired_one;
                ++n_dupMap_entries_decremented;
            } else {
                cerr << NAME << " unknown dupMap value for '" << dupI->first << "': " 
                    << dupI->second << endl;
                return EXIT_FAILURE;
            }
        }

        if ((opt_progress && DEBUG(1)) && n_reads % opt_progress == 0) 
            cerr << NAME << "[pass2] dupMap operations:"
                << " erased " << n_dupMap_entries_erased_SE << " SE, "
                << " erased " << n_dupMap_entries_erased_PE << " PE, "
                << " decremented " << n_dupMap_entries_decremented << " PE halves" << endl;
        if ((opt_progress || DEBUG(1)) && n_reads % opt_progress == 0) 
            cerr << NAME << "[pass2] "
                << n_reads << " reads seen, last at RefID = " << al.RefID << " Pos = " << al.Position << ", "
                << n_reads_written_to_output << " written to " << output_file << ", "
                << n_reads_written_to_dups << " written to " << duplicate_file << ", "
                << n_reads_removed << " removed" << endl;
	}

    if (opt_progress && DEBUG(1))
        cerr << NAME << "[pass2] dupMap operations: "
            << " erased " << n_dupMap_entries_erased_SE << " SE, "
            << " erased " << n_dupMap_entries_erased_PE << " PE, "
            << " decremented " << n_dupMap_entries_decremented << " PE halves" << endl;
    if ((opt_progress || DEBUG(1)) && n_reads % opt_progress == 0) 
        cerr << NAME << "[pass2] "
            << n_reads << " reads seen, "
            << n_reads_written_to_output << " written to " << output_file << ", "
            << n_reads_written_to_dups << " written to " << duplicate_file << ", "
            << n_reads_removed << " removed" << endl;

    IF_DEBUG(2) {
        cerr << n_reads_pass1 << " reads in pass 1" << endl;
        cerr << n_reads << " reads in pass 2" << endl;
    }

	reader.Close();
	writer.Close();
    if (opt_duplicatefile)
        writer_dups.Close();

	return EXIT_SUCCESS;
}


//-------------------------------------
//-------------------------------------  local functions
//-------------------------------------


static void
listAlignments(const alignmentList& al_set)
{
    for (alignmentListCI i = al_set.begin(); i != al_set.end(); ++i) {
        printAlignmentInfo(cerr, (*i), 99);
    }
}


//-------------------------------------


static void
determineDuplicates(alignmentList& al_set, alignmentList& al_dups)
{
    const string HERE = "determineDuplicates():";
    size_t initial_size = al_set.size();
    IF_DEBUG(2) cerr << HERE << " received " << initial_size << " reads" << endl;

    alignmentListI al_i, al_j;

    // pass 0, exclude easy cases first

    int n0_paired_single_only = 0;
    int n0_single_paired_only = 0;
    int n0_unmapped = 0;
    int n0_mate_unmapped = 0;

    al_i = al_set.begin();
    while (al_i != al_set.end()) {
        if (opt_detect == DETECT_single_only && al_i->IsPaired()) {
            IF_DEBUG(3) cerr << HERE << " " << al_i->Name << " is paired and --single-only, excluded" << endl;
            al_i = al_set.erase(al_i); 
            ++n0_paired_single_only;
        } else if (opt_detect == DETECT_paired_only && ! al_i->IsPaired()) {
            IF_DEBUG(3) cerr << HERE << " " << al_i->Name << " is single and --paired-only, excluded" << endl;
            al_i = al_set.erase(al_i); 
            ++n0_single_paired_only;
        } else if (! al_i->IsMapped()) { // no dup if not mapped
            IF_DEBUG(3) cerr << HERE << " " << al_i->Name << " is not mapped, excluded" << endl;
            al_i = al_set.erase(al_i); 
            ++n0_unmapped;
        } else if (opt_detect != DETECT_as_single // ignore mate if --as-single-end
                   && al_i->IsPaired() && ! al_i->IsMateMapped()) { // no dup if mate not mapped
            IF_DEBUG(3) cerr << HERE << " " << al_i->Name << " has a mate that is not mapped, excluded" << endl;
            al_i = al_set.erase(al_i); 
            ++n0_mate_unmapped;
        } else {
            ++al_i;
        }
    }

    IF_DEBUG(2) {
        if (al_set.size() < initial_size || DEBUG(3)) {
            cerr << HERE << " done with pass 0, " << al_set.size() << " reads left";
            if (n0_paired_single_only) cerr << ", paired w/ single-only = " << n0_paired_single_only;
            if (n0_single_paired_only) cerr << ", single w/ paired-only = " << n0_single_paired_only;
            if (n0_unmapped) cerr << ", unmapped = " << n0_unmapped;
            if (n0_mate_unmapped) cerr << ", mate unmapped = " << n0_mate_unmapped;
            cerr << endl;
        }
    }

    if (al_set.empty()) {
        return;
    }

    // Starting with first read, check later reads as dups against it, keeping
    // the best.
    // 
    // If dups are found, we add then to al_dups.  The best read remains in al_set
    // until the end of its pass, then it is deleted from al_set and *not* added
    // to al_dups, so the best read will not be claimed as a dup.
    //
    // There is an issue here for pairs... my keeping one arbitrarily means that
    // its mate is indicated as a dup prior to my having evaluated it.  Maybe I
    // should consider the reads at a site in read-name order, or impose some other
    // ordering on the way in which I consider the reads to be dups.  I hope I do
    // not have to keep more info...
    //
    // al_set: the set of reads we are scanning for dups; we remove each read
    // as we consider it and move to the next.  at the end, this should be
    // empty
    //
    // al_dups: contains the duplicates as determined from all
    // post-single-end-cases reads, so the presence of an alignment in al_dups
    // means that it is a duplicate

    int cycle = 0;
    al_i = al_j = al_set.begin();
    do {
        ++al_j;  // compare to the next one...
        ++cycle;

        bool found_a_match = false;

        IF_DEBUG(2) cerr << HERE << " starting cycle, looking for duplicates of " << al_i->Name << endl;

        while (al_j != al_set.end()) {
            // al_i is the "best" read that we're comparing al_j against

            if (isDuplicate((*al_i), (*al_j))) {

                found_a_match = true;
                if (al_j->MapQuality <= al_i->MapQuality) {
                    al_dups.push_back((*al_j));  // add second read to dups list
                    al_j = al_set.erase(al_j);   // erase second read from al_set, continue with next read
                } else {
                    IF_DEBUG(2) cerr << HERE << " second read has better map quality" << endl;
                    al_dups.push_back((*al_i));  // add first read to dups list
                    al_set.erase(al_i);          // erase first read from al_set
                    al_i = al_j;                 // reset best read for these dups
                    ++al_j;                      // continue with next read
                }

            } else {
                ++al_j;
            }

        }

        // so al_i either had no duplicates (found_a_match == false) or is the
        // best of the duplicates seen in al_set, and the other duplicates in
        // al_set have been removed from al_set and placed in al_dups.  in
        // either case, i want to remove al_i from al_set now, and restart
        // looking for dups at the beginning of al_set, which has an alignment
        // that has not yet been determined to be a dup.

        if (! found_a_match)
            IF_DEBUG(2) cerr << HERE << " " << al_i->Name << " had no duplicates" << endl;

        al_set.erase(al_i);
        al_i = al_j = al_set.begin();

        IF_DEBUG(2) cerr << HERE << " end of cycle " << cycle << ", " << al_set.size() 
                << " reads left to consider" << endl;

    } while (! al_set.empty());

    IF_DEBUG(2) {
        if (al_dups.size() > 0 || DEBUG(2))
            cerr << HERE << " *** received " << initial_size << " reads, returning " 
                << al_set.size() << " non-dup reads, " << al_dups.size() << " dup reads" << endl;
        IF_DEBUG(3) listAlignments(al_dups);
    }
}


//-------------------------------------


static bool
isDuplicate(const BamAlignment& al_i, const BamAlignment& al_j)
{
    const string HERE = "isDuplicate():";
    string i_tag, j_tag;
    // we already know that these alignments are mapped, and to the same reference at the same position
    if (   al_j.RefID                 == al_i.RefID        // same reference
        && al_j.Position              == al_i.Position     // same position
        && al_j.IsReverseStrand()     == al_i.IsReverseStrand()   // same orientation
        && al_j.GetTag("RG", j_tag)   == al_i.GetTag("RG", i_tag)  // has a RG tag?
        && j_tag                      == i_tag             // RG tag is the same?
        && (opt_detect == DETECT_as_single  // ignore pair-dependent qualities with --as-single-end
           || (   al_j.IsPaired()            == al_i.IsPaired()     // same pairing
               && al_j.MateRefID             == al_i.MateRefID      // mates mapped to same sequence
               && al_j.MatePosition          == al_i.MatePosition // mates mapped to same position
               && al_j.IsMateReverseStrand() == al_i.IsMateReverseStrand())) // mates same orientation
        && al_j.QueryBases.length()   == al_i.QueryBases.length() // same read length
        && al_j.AlignedBases.length() == al_i.AlignedBases.length() // same alignment length
        // need to include some notion of optical distance?
        ) {
        IF_DEBUG(2) cerr << HERE << " " << al_j.Name << " is a duplicate of " << al_i.Name << endl;
        return true;
    }
    return false;
}


//-------------------------------------


static void
diagnoseDuplicate(const BamAlignment& al_i, const BamAlignment& al_j)
{
    const string HERE = "diagnoseDuplicate():";
    cerr << HERE << " " << al_i.Name << " vs " << al_j.Name << endl;
    printAlignmentInfo(cerr, al_i, 99);
    printAlignmentInfo(cerr, al_j, 99);
    string i_tag, j_tag;
    // we already know that these alignments are mapped, and to the same reference at the same position
    if (al_j.RefID != al_i.RefID)
        { cerr << HERE << " mismatch RefID" << endl; return; }
    if (al_j.Position != al_i.Position)
        { cerr << HERE << " mismatch Position" << endl; return; }
    if (al_j.IsReverseStrand() != al_i.IsReverseStrand())
        { cerr << HERE << " mismatch IsReverseStrand()" << endl; return; }
    if (al_j.GetTag("RG", j_tag) != al_i.GetTag("RG", i_tag))
        { cerr << HERE << " mismatch GetTag(\"RG\")" << endl; return; }
    if (j_tag != i_tag)
        { cerr << HERE << " mismatch RG: tag value" << endl; return; }
    if (opt_detect == DETECT_as_single)
        { cerr << HERE << " skipped pair stuff, DETECT_as_single" << endl; return; }
    else {
        if (al_j.IsPaired() != al_i.IsPaired())
            { cerr << HERE << " mismatch IsPaired()" << endl; return; }
        if (al_j.MateRefID != al_i.MateRefID)
            { cerr << HERE << " mismatch MateRefID" << endl; return; }
        if (al_j.MatePosition != al_i.MatePosition)
            { cerr << HERE << " mismatch MatePosition" << endl; return; }
        if (al_j.IsMateReverseStrand() != al_i.IsMateReverseStrand())
            { cerr << HERE << " mismatch IsMateReverseStrand()" << endl; return; }
    }
    if (al_j.QueryBases.length() != al_i.QueryBases.length())
        { cerr << HERE << " mismatch QueryBases.length()" << endl; return; }
    if (al_j.AlignedBases.length() != al_i.AlignedBases.length())
        { cerr << HERE << " mismatch AlignedBases.length()" << endl; return; }
    cerr << HERE << " " << al_i.Name << " and " << al_j.Name << " are duplicates" << endl;
}


//-------------------------------------


static void
dump_dupMap(const dupMap& this_dm)
{
    cerr << "* dupMap contains " << this_dm.size() << " elements" << endl;
    size_t i = 1;
    for (dupMapCI dmI = this_dm.begin(); dmI != this_dm.end(); ++dmI, ++i) {
        cerr << "* " << i << " " << dmI->first << " = " << dmI->second << endl;
    }
    cerr << "* dupMap contains " << this_dm.size() << " elements" << endl;
}


//-------------------------------------


static void
query_dupMap(const dupMap& this_dm) 
{
    const string HERE = "query_dupMap():";
    size_t this_size = this_dm.size();
    int64_t n_SE = 0, n_PE2 = 0, n_PE1 = 0, n_UNSET = 0, n_other = 0;
    for (dupMapCI dupI = this_dm.begin(); dupI != this_dm.end(); ++dupI) {
        switch(dupI->second) {
            case dupMap_singleend:   ++n_SE; break;
            case dupMap_UNSET:       ++n_UNSET; break;
            case dupMap_paired_one:  ++n_PE1; break;
            case dupMap_paired_both: ++n_PE2; break;
            default: ++n_other; break;
        }
    }
    cerr << HERE << " size: " << this_size << ", values: SE " << n_SE << ", PE2 " << n_PE2 
        << ", PE1 " << n_PE1 << ", PE0 " << n_UNSET << ", other " << n_other << endl;
}


//-------------------------------------


static void
update_dupMap(alignmentList& al_set, dupMap& this_dm)
{
    const string HERE = "update_dupMap():";
    IF_DEBUG(2) cerr << HERE << " received " << al_set.size() 
        << " duplicate alignments" << endl;

    if (al_set.empty())
        return;

    IF_DEBUG(2) cerr << "*********************************************" << endl;

    int n_reads_received = al_set.size();
    int n_reads_found_in_map = 0;
    int n_SE_found_in_map = 0;
    int n_SE_added = 0;
    int n_PE_first_added = 0;
    int n_PE_second_added = 0;
    int n_PE_mate_upstream = 0;

    alignmentListI aLI_i;

    for (aLI_i = al_set.begin(); aLI_i != al_set.end(); ++aLI_i) {

        dupMapI dupI = this_dm.find(aLI_i->Name);

        if (dupI != this_dm.end()) {
            ++n_reads_found_in_map;
            IF_DEBUG(2) cerr << HERE << " " << aLI_i->Name 
                << " in dupMap, val = " << dupI->second << endl;
        }

        if (! aLI_i->IsPaired()) {  // single-end

            if (dupI != this_dm.end()) {
                ++n_SE_found_in_map;
                cerr << HERE << " ERROR, SE read name already seen for '"
                        << aLI_i->Name << "', is this a duplicate read name??" << endl;
            }

            this_dm[aLI_i->Name] = dupMap_singleend;  // add to map as SE
            ++n_SE_added;
            IF_DEBUG(3) cerr << HERE << " " << aLI_i->Name
                << " SE, set dupMap = -1" << endl;

        } else {  // paired-end

            if (dupI == this_dm.end()) {  // not in map

                if (aLI_i->MateRefID >= 0 && isMateUpstream((*aLI_i))) { 

                    // if mate is upstream and not in the dupMap, it wasn't a dup
                    ++n_PE_mate_upstream;
                    IF_DEBUG(2) cerr << HERE << " " << aLI_i->Name 
                        << " PE, dupMap no mate found" << ", mate UPSTREAM, NOT DUP" << endl;

                } else {

                    this_dm[aLI_i->Name] = dupMap_paired_one;  // add to map as first read of PE
                    ++n_PE_first_added;
                    IF_DEBUG(2) cerr << HERE << " " << aLI_i->Name 
                        << " PE, dupMap no mate found" << ", set dupMap = 1" << endl;

                }

            } else {

                if (dupI->second == dupMap_paired_both) {
                    cerr << HERE << " ERROR, two PE reads already seen for '"
                        << dupI->first << "', is this a duplicate read name??" << endl;
                }
                if (dupI->second != dupMap_paired_one) {
                    cerr << HERE << " ERROR, PE read name already seen for '"
                        << dupI->first << "' and entry is inconsistent" << endl;
                }
                dupI->second = dupMap_paired_both;
                ++n_PE_second_added;
                IF_DEBUG(2) cerr << HERE << " " << aLI_i->Name 
                    << " PE, update dupMap = " << (*dupI).second << endl;

            }
        }
    }

    assert(aLI_i == al_set.end());  // should be none left
    al_set.clear();

    IF_DEBUG(2) {
        cerr << HERE << " received " << n_reads_received;
        cerr << ", found " << n_reads_found_in_map << " in map";
        if (n_SE_found_in_map)
            cerr << ", *** found " << n_SE_found_in_map << " SE reads already in map ***";
        cerr << ", added " << n_SE_added << " SE, " 
            << n_PE_first_added << " PE first, " 
            << n_PE_second_added << " PE second";
        if (n_PE_mate_upstream)
            cerr << ", discarded " << n_PE_mate_upstream << " PE with non-dup mate upstream";
        cerr << endl;
        IF_DEBUG(3) dump_dupMap(this_dm);
        IF_DEBUG(2) cerr << "*********************************************" << endl;
    }
}


//-------------------------------------


static void
clear_dupMap(dupMap& this_dm) 
{
    const string HERE = "clear_dupMap(1 arg):";
    IF_DEBUG(2) cerr << HERE << " size: " << this_dm.size() << endl;
    this_dm.clear();
}


//-------------------------------------


static void
clear_dupMap(dupMap& this_dm, dup_t val) 
{
    const string HERE = "clear_dupMap(2 args):";
    size_t initial_size = this_dm.size();
    dupMapI dupI = this_dm.begin();
    while (dupI != this_dm.end()) {
        if (dupI->second == val) {
            dupMapI erase_it = dupI;
            ++dupI;
            this_dm.erase(erase_it);
        } else {
            ++dupI;
        }
    }
    IF_DEBUG(2) cerr << HERE << " size at entrance: " << initial_size 
        << " size at exit: " << this_dm.size() << endl;
}


