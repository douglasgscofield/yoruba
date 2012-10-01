// yoruba_inu.cpp  (c) Douglas G. Scofield, Dept. Plant Physiology, Umeå University
//
// Gbagbe (English command is forget) dynamically reassociates reads with reference
// sequences in a BAM file.
//
// Gbagbe dynamically reassociating reads with reference sequences in a BAM
// file.  In automatic mode, gbagbe removes reference sequences from the BAM
// header that are not actually referred to by reads.  In directed mode, gbagbe
// will replace the BAM header with one that is provided by the user, and
// reassociate reads with reference sequences in the order they appear in the
// new header.
//
// Inu reads the BAM file structure and summarizes the header, references and read
// contents.  It can also check the validity of the header.
//
// Gbagbe is the Yoruba (Nigeria) verb for 'to forget'.
//
// Uses BamTools C++ API for handling BAM files


// CHANGELOG
//
//
//
// TODO
//
// Update README.md
// Update with options
// Rebuild the BAM file index on option
//

#include "yoruba_gbagbe.h"

using namespace std;
using namespace BamTools;
using namespace yoruba;

static string       input_file;  // will eventually default to stdin, set from command line
static string       output_file;  // defaults to stdout, set with -o FILE
#ifdef _WITH_DEBUG
static int32_t      opt_debug = 0;
static int64_t      opt_reads = -1;
static int64_t      opt_progress = 100000; // 100000;
#endif
static const string delim = "'";
static const string sep = "\t";
static const string endline = "\n";


//-------------------------------------


#ifdef _STANDALONE
int 
main(int argc, char* argv[]) {
    return main_gbagbe(argc, argv);
}
#endif


//-------------------------------------


static int
usage()
{
    cerr << endl;
    cerr << "Usage:   " << YORUBA_NAME << " forget [options] <in.bam>" << endl;
    cerr << "         " << YORUBA_NAME << " gbagbe [options] <in.bam>" << endl;
    cerr << endl;
    cerr << "Either command invokes this function." << endl;
    cerr << endl;
    cerr << "\
Dynamically reduces the number of reference sequences in <in.bam>.\n\
\n\
Options: --o FILE | -o FILE | --output FILE  output file name [default is stdout]\n\
         --? | -? | --help      longer help\n\
\n";
#ifdef _WITH_DEBUG
    cerr << "\
         --debug INT     debug info level INT [" << opt_debug << "]\n\
         --reads INT     only process INT reads [" << opt_reads << "]\n\
         --progress INT  print reads processed mod INT [" << opt_progress << "]\n\
\n";
#endif
    cerr << "Gbagbe is the Yoruba (Nigeria) verb for 'to forget'." << endl;
    cerr << endl;

    return 1;
}


//-------------------------------------


int 
yoruba::main_gbagbe(int argc, char* argv[])
{
    SamHeader new_header;  // the new header we are creating
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
    
    enum { OPT_output,
#ifdef _WITH_DEBUG
        OPT_debug, OPT_reads, OPT_progress,
#endif
        OPT_help };

    CSimpleOpt::SOption gbagbe_options[] = {
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

    CSimpleOpt args(argc, argv, gbagbe_options);

    while (args.Next()) {
        if (args.LastError() != SO_SUCCESS) {
            cerr << NAME << " invalid argument '" << args.OptionText() << "'" << endl;
            return usage();
        }
        if (args.OptionId() == OPT_help) {
            return usage();
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
            return 1;
        }
    }

    if (args.FileCount() > 1) {
        cerr << NAME << " requires at most one BAM file specified as input" << endl;
        return usage();
    } else if (args.FileCount() == 1) {
        input_file = args.File(0);
    } else if (input_file.empty()) {
        cerr << NAME << " can't currently read from stdin, ask Doug about it" << endl;
        return 1;
        input_file = "/dev/stdin";
    }

    // set up output; if file not specified, use stdout or its equivalent
    if (output_file.empty()) {
        output_file = "/dev/stdout";
    }

    // gbagba workflow
    //
    // create new SamHeader object
    // read header, get const SamHeader& ref
    // if refs.empty()
    //   exit 1, no references to forget
    // to the new header object, copy over 
    //   metadata
    //   read group dictionary
    //   programs
    //     update the program list with yoruba forget
    //   comments
    // allocate refd[], vector of uint64_t length reader.GetReferenceCount()
    // while (reader.GetNextAlignmentCore(al) {
    //   if (al.IsMapped())
    //     ++refd[al.RefID];
    // }
    // Go through refd[], for each nonzero entry add its value in the sequence of nonzero entries
    // Create references in new header object 
    // rewind BAM file to beginning of reads
    //

    //----------------- Open file, start reading data

	BamReader reader;

	if (! reader.Open(input_file)) {
        cerr << NAME << "could not open BAM input" << endl;
        return 1;
    }

    if (reader.GetReferenceCount() == 0) {
        cerr << NAME << "no reference sequences found in BAM header" << endl;
        reader.Close();
        return 1;
    }

#ifdef _BAMTOOLS_EXTENSION
    const SamHeader& header = reader.GetConstSamHeader();
#else
    SamHeader header = reader.GetHeader();
#endif

    //----------------- Header metadata

    new_header.Version = header.Version;
    new_header.SortOrder = header.SortOrder;
    new_header.GroupOrder = header.GroupOrder;

    //----------------- Read groups

    new_header.ReadGroups = header.ReadGroups;

    //----------------- Programs

    for (SamProgramConstIterator pcI = header.Programs.ConstBegin(); 
            pcI != header.Programs.End(); ++pcI) {
        if (pcI->ID != new_program.ID) {
            // I would prefer to use duplicate program IDs in the SAM header,
            // if the same program worked over the file twice or more, but the
            // SAM spec says no, so move this one to last if it's already there
            SamProgram this_program = *pcI;  // workaround for weird overloading issue
            new_header.Programs.Add(this_program);
        }
    }
    new_header.Programs.Add(new_program);

    //----------------- Comments

    new_header.Comments = header.Comments;

    //----------------- Pass 1: Determine which references are used

    if (opt_progress)
        cerr << NAME << "[pass1] " << reader.GetReferenceCount() 
            << " references in the old BAM" << endl;
    else
        _DEBUG(1) cerr << NAME << "[pass1] " << reader.GetReferenceCount() 
            << " references in the old BAM" << endl;
    vector<int64_t> refs_used( reader.GetReferenceCount() );
    int64_t n_reads = 0;  // number of reads processed
	BamAlignment al;  // holds the current read from the BAM file

	while (reader.GetNextAlignmentCore(al) && (opt_reads < 0 || n_reads < opt_reads)) {

        ++n_reads;
        if (al.IsMapped()) {
            assert(al.RefID >= 0);
            ++refs_used[al.RefID];
        }
        if (opt_progress && n_reads % opt_progress == 0)
            cerr << NAME << "[pass1] " << n_reads << " reads examined..." << endl;
 
	}
    if (opt_progress)
        cerr << NAME << "[pass1] " << n_reads << " reads examined" << endl;
    else
        _DEBUG(1) cerr << NAME << "[pass1] " << n_reads << " reads examined" << endl;

    //----------------- Pass 2: create new reference set

    const RefVector& old_refs = reader.GetReferenceData();
    RefVector        new_refs;

    assert(new_header.Sequences.IsEmpty());  // new_refs contains the new @SQ info
    int32_t new_RefID = 0;
    for (size_t i = 0; i < refs_used.size(); ++i) {
        if (refs_used[i] > 0) {
            new_refs.push_back(old_refs[i]);
            // new_header.Sequences.Add(old_refs[i].RefName, old_refs[i].RefLength + 100000);
            //new_header.Sequences.Add(old_refs[i].RefName, 100000);
            refs_used[i] = new_RefID++;  // entry now contains new reference ID
        } else refs_used[i] = -1;
    }
    assert(new_refs.size() == (size_t)new_RefID);
    //assert(! new_header.Sequences.IsEmpty());  // new_refs contains the new @SQ info
    if (opt_progress)
        cerr << NAME << "[pass2] " << new_RefID << " references in the new BAM" << endl;
    else
        _DEBUG(1) cerr << NAME << "[pass2] " << new_RefID << " references in the new BAM" << endl;
    _DEBUG(1) {
        for (size_t i = 0; i < new_refs.size(); ++i) {
            cerr << NAME << "[pass2] " << i << " new_refs RefName='" << new_refs[i].RefName
                << "' RefLength=" << new_refs[i].RefLength << endl;
        }
        for (SamSequenceConstIterator sscI = header.Sequences.ConstBegin(); 
                sscI != header.Sequences.ConstEnd(); ++sscI) {
            cerr << NAME << "[pass2] " << " Sequences RefName='" << sscI->Name
                << "' RefLength=" << sscI->Length << endl;
        }
    }

    //----------------- Pass 2: second pass through reads, write new BAM file

    BamWriter writer;

    if (! writer.Open(output_file, new_header, new_refs)) {
        cerr << NAME << " could not open output " << output_file << endl;
        return 1 ;
    }

    int64_t n_reads_pass1 = n_reads;
    n_reads = 0;
    int64_t n_reads_rerefd = 0;  // number of reads given re-references

    reader.Rewind();

	while (reader.GetNextAlignmentCore(al) && (opt_reads < 0 || n_reads < opt_reads)) {

        ++n_reads;

        if (al.IsMapped()) {
            assert(al.RefID >= 0);  // it was valid before...
            if (al.RefID != refs_used[al.RefID]) {  // the reference ID is different
                ++n_reads_rerefd;
                al.RefID = refs_used[al.RefID];
                assert(al.RefID >= 0);  // and it is valid after
            }
        }

        writer.SaveAlignment(al);

        if (opt_progress && n_reads % opt_progress == 0)
            cerr << NAME << "[pass2] " << n_reads << " reads processed, " 
                << n_reads_rerefd << " rereferenced" << endl;
 
	}
    if (opt_progress)
        cerr << NAME << "[pass2] " << n_reads << " reads examined from the BAM file" << endl;
    else
        _DEBUG(1) cerr << NAME << "[pass2] " << n_reads << " reads examined from the BAM file" << endl;
    assert(n_reads == n_reads_pass1);
    if (opt_progress)
        cerr << NAME << "[pass2] " << n_reads_rerefd << " reads rereferenced" << endl;
    else 
        _DEBUG(1) cerr << NAME << "[pass2] " << n_reads_rerefd << " reads rereferenced" << endl;

	reader.Close();
	writer.Close();

	return 0;
}

