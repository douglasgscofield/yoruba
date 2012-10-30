// yoruba_gbagbe.cpp  (c) Douglas G. Scofield, Dept. Plant Physiology, Umeå University
//
// Gbagbe (English command is forget) dynamically reassociates reads with reference
// sequences in a BAM file.
//
// Gbagbe dynamically reassociating reads with reference sequences in a BAM
// file, removing reference sequences from the BAM header that are not actually
// referred to by reads.  Mates that reference sequences that the aligned reads
// do not reference will also have their reference sequences kept in the header,
// unless --no-mate is specified, then the mates' reference sequences will be 
// changed to -1.
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

#include "yoruba_gbagbe.h"

using namespace std;
using namespace BamTools;
using namespace yoruba;

static string       input_file;
static string       output_file;  // defaults to stdout, set with -o FILE
static bool         opt_mate = true; // with --no-mate, forget references for mates too
#ifdef _WITH_DEBUG
static int32_t      opt_debug = 0;
static int64_t      opt_reads = -1;
static int64_t      opt_progress = 0; // 1000000;
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
usage(bool longer = false)
{
    cerr << endl;
    cerr << "Usage:   " << YORUBA_NAME << " forget [options] <in.bam>" << endl;
    cerr << "         " << YORUBA_NAME << " gbagbe [options] <in.bam>" << endl;
    cerr << endl;
    cerr << "Either command invokes this function." << endl;
    cerr << endl;
    cerr << "\
Dynamically reduces the number of reference sequences in <in.bam>.\n\
\n";
    if (longer) cerr << "\
Alignments within a BAM file mention reference sequences.  The BAM header may\n\
contain many reference sequence descriptions (@SQ lines) to which no alignments \n\
refer.  This function removes these unmentioned reference sequence descriptions,\n\
and can be especially helpful for reducing space and loading time when the BAM \n\
file contains alignments extracted from a region of a larger BAM containing\n\
alignments against many, many reference sequences.  If the original BAM file \n\
mentioned a few hundred references, it is probably not much of a problem. \n\
10 million reference descriptions, on the other hand, take a while to load... \n\
\n\
Each paired-end read with an aligned mate mentions the reference sequence of\n\
its mate.  These reference sequence descriptions will also be kept in the\n\
output BAM file unless the --no-mate option is given.  With this option, such\n\
mates will have their reference sequence ID set to -1, which indicates a missing \n\
reference sequence description.\n\
\n";
    cerr << "\
Options: --no-mate                             also forget references for paired-end mates\n\
         --o FILE | -o FILE | --output FILE    output file name [default is stdout]\n\
         --? | -? | --help                     longer help\n\
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

    return EXIT_FAILURE;
}


//-------------------------------------


int 
yoruba::main_gbagbe(int argc, char* argv[])
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
    
    enum { OPT_output, OPT_nomate, OPT_mate,
#ifdef _WITH_DEBUG
        OPT_debug, OPT_reads, OPT_progress,
#endif
        OPT_help };

    CSimpleOpt::SOption gbagbe_options[] = {
        { OPT_nomate,          "--no-mate",         SO_NONE }, 
        //{ OPT_mate,            "--mate",            SO_NONE }, 
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
            return usage(true);
        //} else if (args.OptionId() == OPT_mate) {
        //    opt_mate = true;
        } else if (args.OptionId() == OPT_nomate) {
            opt_mate = false;
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
    if (output_file.empty()) {
        output_file = "/dev/stdout";
    }


    //----------------- Open input BAM, create header for output BAM


	BamReader reader;

	if (! reader.Open(input_file)) {
        cerr << NAME << "could not open BAM input" << endl;
        return EXIT_FAILURE;
    }

    if (reader.GetReferenceCount() == 0) {
        cerr << NAME << "no reference sequences found in BAM header" << endl;
        reader.Close();
        return EXIT_FAILURE;
    }

#ifdef _BAMTOOLS_EXTENSION
    const SamHeader& header = reader.GetConstSamHeader();
#else
    SamHeader header = reader.GetHeader();
#endif

    //----------------- Create the header for the new BAM

    SamHeader new_header;

    new_header.Version = header.Version;
    new_header.SortOrder = header.SortOrder;
    new_header.GroupOrder = header.GroupOrder;

    new_header.ReadGroups = header.ReadGroups;

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

    new_header.Comments = header.Comments;


    //----------------- Pass 1: Determine which references are used


    if (opt_progress || DEBUG(1))
        cerr << NAME << "[pass1] " << reader.GetReferenceCount() 
            << " references in the input BAM" << endl;

    vector<int64_t> refs_used( reader.GetReferenceCount() );

    int64_t n_reads = 0;  // number of reads processed
	BamAlignment al;  // holds the current read from the BAM file

	while (reader.GetNextAlignmentCore(al) && (opt_reads < 0 || n_reads < opt_reads)) {

        ++n_reads;
        if (al.IsMapped()) {
            if (al.RefID < 0) {
                cerr << NAME << "[pass1] missing reference sequence from input bam" << endl;
                return EXIT_FAILURE;
            }
            ++refs_used[al.RefID];
        }
        if (opt_mate && al.IsPaired() && al.IsMateMapped() && al.MateRefID >= 0)
            // if a reference is missing for a mapped mate then MateRefID == -1,
            // an unmapped mate has our RefID and Position, so not a reference "use"
            ++refs_used[al.MateRefID];

        if ((opt_progress || DEBUG(1)) && n_reads % opt_progress == 0)
            cerr << NAME << "[pass1] " << n_reads << " reads examined..." << endl;
 
	}
    if (opt_progress || DEBUG(1))
        cerr << NAME << "[pass1] " << n_reads << " reads examined" << endl;


    //----------------- Pass 2: Create new reference set


    const RefVector& old_refs = reader.GetReferenceData();
    RefVector        new_refs;
    RefVector        blank_refs;

    assert(new_header.Sequences.IsEmpty());  // new_refs contains the new @SQ info
    int32_t new_RefID = 0;
    for (size_t i = 0; i < refs_used.size(); ++i) {

        if (refs_used[i] > 0) {
            new_refs.push_back(old_refs[i]);
#ifdef _BAMTOOLS_EXTENSION
            SamSequenceConstIterator refI = header.Sequences.ConstFind(old_refs[i].RefName);
            if (refI == header.Sequences.ConstEnd()) {
                cerr << NAME << "[pass2] internal header inconsistency, " << old_refs[i].RefName 
                    << " not found in the input header" << endl;
                return EXIT_FAILURE;
            }
            const SamSequence& existing_ref = (*refI);
#else
            const SamSequence& existing_ref = header.Sequences[old_refs[i].RefName];
#endif
            new_header.Sequences.Add(existing_ref);
            refs_used[i] = new_RefID;  // entry now contains new reference ID
            ++new_RefID;

        } else {
            refs_used[i] = -1;
        }
    }
    assert(new_refs.size() == (size_t)new_RefID);

    if (opt_progress || DEBUG(1))
        cerr << NAME << "[pass2] " << new_RefID << " reference"
            << (new_RefID == 1 ? "" : "s") << " in the output BAM" << endl;

    IF_DEBUG(1) {
        for (size_t i = 0; i < new_refs.size(); ++i) {
            cerr << NAME << "[pass2] " << i << "] SN:" << new_refs[i].RefName
                << "  LN:" << new_refs[i].RefLength << endl;
        }
        IF_DEBUG(2) {
            SamSequenceConstIterator sscI = new_header.Sequences.ConstBegin();
            if (sscI == new_header.Sequences.ConstEnd()) {
                cerr << NAME "[pass2] no entries in new_header.Sequences" << endl;
            } else {
                for (; sscI != new_header.Sequences.ConstEnd(); ++sscI)
                    cerr << NAME << "[pass2] new_header " << " Sequences RefName=" << sscI->Name
                        << "  RefLength=" << sscI->Length << endl;
            }
        }
    }


    //----------------- Pass 2: Second pass through reads, write new BAM file


    BamWriter writer;

    IF_DEBUG(2) {
        cerr << "********* BEGIN new_header.ToString()" << endl;
        cerr << new_header.ToString();
        cerr << "********* END   new_header.ToString()" << endl;
    }


    if (! writer.Open(output_file, new_header, new_refs)) {
    //if (! writer.Open(output_file, new_header, blank_refs)) {
        cerr << NAME << " could not open output " << output_file << endl;
        return EXIT_FAILURE;
    }
    if (false) {
        reader.Close();
        writer.Close();
        return EXIT_FAILURE;
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
                if (al.IsPaired()) {
                    if (al.MateRefID == al.RefID
                        || (opt_mate && al.MateRefID >= 0)) {
                        // update the mate RefID
                        al.MateRefID = refs_used[al.MateRefID];
                    } else if (al.IsMateMapped()) {
                        // mate ref is now unavailable
                        al.MateRefID = -1;
                    }
                }
                al.RefID = refs_used[al.RefID];
                assert(al.RefID >= 0);  // and it is valid after
            }
        }

        writer.SaveAlignment(al);

        if (opt_progress && n_reads % opt_progress == 0)
            cerr << NAME << "[pass2] " << n_reads << " reads examined, " 
                << n_reads_rerefd << " rereferenced in the output BAM" << endl;
 
	}
    if (opt_progress || DEBUG(1))
        cerr << NAME << "[pass2] " << n_reads << " reads examined, " 
            << n_reads_rerefd << " reads rereferenced in the output BAM" << endl;
    assert(n_reads == n_reads_pass1);

	reader.Close();
	writer.Close();

	return 0;
}

