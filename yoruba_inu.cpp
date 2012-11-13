// yoruba_inu.cpp  (c) Douglas G. Scofield, Dept. Plant Physiology, Umeå University
//
// Inu (English command is inside) summarizes the contents of a BAM file.
//
// Inu reads the BAM file structure and summarizes the header, references and read
// contents.  It can also check the validity of the header.
//
// Inu is the Yoruba (Nigeria) noun for 'inside'.
//
// Uses BamTools C++ API for handling BAM files

// CHANGELOG
//
// --- Add option to dump tags
// --- Add options to dump sequence, aligned sequence, qualities?

#include "yoruba_inu.h"

using namespace std;
using namespace BamTools;
using namespace yoruba;

static string       input_file;  // defaults to stdin, set from command line
static int64_t      opt_reads_to_report = 10;
static bool         opt_continue = false;
static bool         opt_validate = false;
static int32_t      opt_refs_to_report = 10;
#ifdef _WITH_DEBUG
static int32_t      opt_debug = 0;
static int32_t      debug_progress = 100000;
static int64_t      opt_reads = -1;
static int64_t      opt_progress = 0; // 100000;
#endif
static const string delim = "'";
static const string sep = "\t";
static const string endline = "\n";


//-------------------------------------


#ifdef _STANDALONE
int 
main(int argc, char* argv[]) {
    return main_inu(argc, argv);
}
#endif


//-------------------------------------


static int
usage()
{
    cerr << endl;
    cerr << "\
Usage:   " << YORUBA_NAME << " inside [options] <in.bam>\n\
         " << YORUBA_NAME << " inu [options] <in.bam>\n\
\n\
Summarizes the contents of the BAM file <in.bam>.  Either command\n\
invokes this function.\n\
\n\
Output includes:\n\
   (1) header lines exclusive of reference sequences\n\
   (2) the first " << opt_refs_to_report << " reference sequences\n\
   (3) mapping characteristics of the first " << opt_reads_to_report << " reads\n\
\n\
Options: --reads-to-report INT   print this many reads [" << opt_reads_to_report << "]\n\
         --refs-to-report INT    print this many references [" << opt_refs_to_report << "]\n\
         --continue              continue counting reads until the end of the BAM\n\
         --validate              check validity using BamTools API; very strict\n\
         -? | --help             longer help\n\
\n";
#ifdef _WITH_DEBUG
    cerr << "\
         --debug INT      debug info level INT [" << opt_debug << "]\n\
         --reads INT      only process INT reads [" << opt_reads << "]\n\
         --progress INT   print reads processed mod INT [" << opt_progress << "]\n\
\n";
#endif
    cerr << "Inu is the Yoruba (Nigeria) noun for 'inside'." << endl;
    cerr << endl;

    return EXIT_FAILURE;
}


//-------------------------------------


int 
yoruba::main_inu(int argc, char* argv[])
{
    //----------------- Command-line options

	if( argc < 2 ) {
		return usage();
	}

    enum { OPT_reads_to_report, OPT_refs_to_report, OPT_continue, OPT_validate, 
#ifdef _WITH_DEBUG
        OPT_debug, OPT_reads, OPT_progress,
#endif
        OPT_help };

    CSimpleOpt::SOption inu_options[] = {
        { OPT_refs_to_report,  "--refs-to-report",  SO_REQ_SEP },
        { OPT_reads_to_report, "--reads-to-report", SO_REQ_SEP },
        { OPT_continue,        "--continue",        SO_NONE },
        { OPT_validate,        "--validate",        SO_NONE },
        { OPT_help,            "--help",            SO_NONE },
        { OPT_help,            "-?",                SO_NONE }, 
#ifdef _WITH_DEBUG
        { OPT_debug,           "--debug",           SO_REQ_SEP },
        { OPT_reads,           "--reads",           SO_REQ_SEP },
        { OPT_progress,        "--progress",        SO_REQ_SEP },
#endif
        SO_END_OF_OPTIONS
    };

    CSimpleOpt args(argc, argv, inu_options);

    while (args.Next()) {
        if (args.LastError() != SO_SUCCESS) {
            cerr << NAME << " invalid argument '" << args.OptionText() << "'" << endl;
            return usage();
        }
        if (args.OptionId() == OPT_help)       return usage();
        else if (args.OptionId() == OPT_reads_to_report) 
            opt_reads_to_report = strtoll(args.OptionArg(), NULL, 10);
        else if (args.OptionId() == OPT_refs_to_report) 
            opt_refs_to_report = strtol(args.OptionArg(), NULL, 10);
        else if (args.OptionId() == OPT_continue)  opt_continue = true;
        else if (args.OptionId() == OPT_validate) opt_validate = true;
#ifdef _WITH_DEBUG
        else if (args.OptionId() == OPT_debug) 
            opt_debug = args.OptionArg() ? atoi(args.OptionArg()) : opt_debug;
        else if (args.OptionId() == OPT_reads) 
            opt_reads = strtoll(args.OptionArg(), NULL, 10);
        else if (args.OptionId() == OPT_progress) 
            opt_progress = args.OptionArg() ? strtoll(args.OptionArg(), NULL, 10) : opt_progress;
#endif
        else {
            cerr << NAME << " unprocessed argument '" << args.OptionText() << "'" << endl;
            return EXIT_FAILURE;
        }
    }

    if (DEBUG(1) && ! opt_progress)
        opt_progress = debug_progress;

    if (args.FileCount() > 1) {
        cerr << NAME << " requires at most one BAM file specified as input" << endl;
        return usage();
    } else if (args.FileCount() == 1) {
        input_file = args.File(0);
    } else if (input_file.empty()) {
        input_file = "/dev/stdin";
    }

    //----------------- Open file, start reading data

	BamReader reader;

	if (! reader.Open(input_file)) {
        cerr << NAME << " could not open BAM input" << endl;
        return EXIT_FAILURE;
    }

#ifdef _BAMTOOLS_EXTENSION
    const SamHeader& header = reader.GetConstSamHeader();
#else
    SamHeader header = reader.GetHeader();
#endif

    if (opt_validate) {
        if (! header.IsValid(true)) { // this check is very strict
            cout << NAME << " header not well-formed, errors are:" << endl;
            cout << header.GetErrorString() << endl;
        }
    }

    //----------------- Header metadata

    if (header.HasVersion() || header.HasSortOrder() || header.HasGroupOrder()) {
        cout << NAME << "[headerline]";
        if (header.HasVersion()) 
            cout << sep << "VN:" << delim << header.Version << delim;
        if (header.HasSortOrder()) 
            cout << sep << "SO:" << delim << header.SortOrder << delim;
        if (header.HasGroupOrder()) 
            cout << sep << "GO:" << delim << header.GroupOrder << delim;
        cout << endl;
    } else cout << NAME << "[headerline] no header line found" << endl;

    //----------------- Reference sequences

    const RefVector& refs = reader.GetReferenceData();

    if (header.HasSequences()) {
        int32_t ref_count = reader.GetReferenceCount();
        if (ref_count > opt_refs_to_report)
            cout << NAME << "[ref] displaying the first " << opt_refs_to_report 
                << " reference sequences" << endl;
        for (int32_t i = 0; i < ref_count && i < opt_refs_to_report; ++i) {
            cout << NAME << "[ref] " << i << " ";
            cout << "@SQ";
            // these tags must exist for a reference sequence
            cout << sep << "NM:" << delim << refs[i].RefName << delim
                << sep << "LN:" << refs[i].RefLength
                << endline;
        }
        cout << NAME << "[ref] " << ref_count << " reference sequences found" << endl;
    } else cout << NAME << "[ref] no reference sequences found" << endl;

    //----------------- Read groups

    if (header.HasReadGroups()) {
        const string prefix = NAME "[readgroup] ";
        printReadGroupDictionary(cout, header.ReadGroups, prefix, "@RG", "\t", "'", "\n");
    } else cout << NAME << "[readgroup] no read group dictionary found" << endl;

    //----------------- Programs

    if (header.HasPrograms()) {
        for (SamProgramConstIterator pcI = header.Programs.ConstBegin();
                pcI != header.Programs.ConstEnd(); ++pcI) {
            cout << NAME << "[program] ";
            cout << "@PG";
            if ((*pcI).HasID())
                cout << sep << "ID:" << delim << (*pcI).ID << delim;
            if ((*pcI).HasName())
                cout << sep << "PN:" << delim << (*pcI).Name << delim;
            if ((*pcI).HasCommandLine())
                cout << sep << "CL:" << delim << (*pcI).CommandLine << delim;
            if ((*pcI).HasPreviousProgramID())
                cout << sep << "PP:" << delim << (*pcI).PreviousProgramID << delim;
            if ((*pcI).HasVersion())
                cout << sep << "VN:" << delim << (*pcI).Version << delim;
            cout << endline;
        }
    } else cout << NAME << "[program] no program information found" << endl;


    //----------------- Comments

    if (! header.Comments.empty()) {
        for (vector<string>::const_iterator vI = header.Comments.begin();
                vI < header.Comments.end(); ++vI) {
            cout << NAME << "[program] ";
            cout << "@CO";
            cout << sep << delim << (*vI) << delim << endline;
        }
    } else cout << NAME << "[comment] no comment lines found" << endl;

    //----------------- Reads

	BamAlignment al;  // holds the current read from the BAM file

    int64_t n_reads = 0;  // number of reads processed

    if (opt_reads_to_report) {
        cout << NAME << "[read] printing the first " << opt_reads_to_report << " reads" << endl;
    }

	while (reader.GetNextAlignmentCore(al) && (opt_reads < 0 || n_reads < opt_reads)) {

        ++n_reads;

        if (n_reads <= opt_reads_to_report) {
            al.BuildCharData();
            cout << NAME << "[read] ";
            printAlignmentInfo(cout, al, refs, 99);
        }

        if (opt_progress && n_reads % opt_progress == 0)
            cerr << NAME << "[read] " << n_reads << " reads processed..." << endl;

        if (! opt_continue && n_reads == opt_reads_to_report)
            break;
	}

    cout << NAME << "[read] " << n_reads << " reads examined from the BAM file" << endl;

	reader.Close();

	return EXIT_SUCCESS;
}

