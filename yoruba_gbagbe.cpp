// yoruba_gbagbe.cpp  (c) Douglas G. Scofield, douglasgscofield@gmail.com
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
// --- deal with references mentioned in multiple-mapping and other tags
// --- much more robust --list file reading (eg find/create a class for reading delimited files)
// xxx add --usage-file
// xxx more usefully handle a missing reference sequence in the input (al.RefID == -1)
// xxx implement mention-by-name
// xxx clean up mention-by-mates (currently doing message output)
// xxx add --usage
//

#include "yoruba_gbagbe.h"

using namespace std;
using namespace BamTools;
using namespace yoruba;

static string       input_file;
static string       output_file;
static bool         opt_usageonly = false;
static string       usage_file;
static bool         opt_mate = true;
static string       list_file;
#ifdef _WITH_DEBUG
static int32_t      opt_debug = 0;
static int32_t      debug_progress = 100000;
static int64_t      opt_reads = -1;
static int64_t      opt_progress = 0; // 1000000;
#endif
static const string sep = "\t";
class refStats {  // holds statistics for reference sequences
    public:
        string  ref; 
        int32_t old_id;
        int64_t m_read, m_mate; 
        bool    m_name, no_mate; 
        int32_t new_id;
};


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
    cerr << "\n\
Dynamically reduce the number of reference sequences in <in.bam>.\n\
Either command invokes this function.\n\
\n\
NOTE: This command does not adjust reference sequence mentioned within tags.\n\
      There are some de facto standards for these mentions, for example bwa\n\
      with multiply-mapped reads, and this command will handle these as I\n\
      learn of them and have time to implement them.\n\
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
\n\
A list of reference sequences to keep regardless of whether they are referred\n\
to can be provided with the --list option.  The file provided can be in BED\n\
format or contains whitespace-separated fields with the reference sequence name\n\
as the first field.\n\
\n";
    cerr << "\
Options: --no-mate                 also forget references for paired-end mates\n\
         --usage-only              analyze reference usage, do not produce output BAM\n\
         --usage-file FILE         write per-reference usage details to FILE\n\
         -L FILE | --list FILE     file containing names of reference sequences to keep\n\
         -o FILE | --output FILE   output file name [default is stdout]\n\
         -? | --help               longer help\n\
\n";
#ifdef _WITH_DEBUG
    cerr << "\
         --debug INT      debug info level INT [" << opt_debug << "]\n\
         --reads INT      only process INT reads [" << opt_reads << "]\n\
         --progress INT   print reads processed mod INT [" << opt_progress << "]\n\
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
    
    enum { OPT_output, OPT_nomate, OPT_usageonly, OPT_usagefile, OPT_list,
#ifdef _WITH_DEBUG
        OPT_debug, OPT_reads, OPT_progress,
#endif
        OPT_help };

    CSimpleOpt::SOption gbagbe_options[] = {
        { OPT_nomate,          "--no-mate",         SO_NONE }, 
        { OPT_usageonly,       "--usage-only",      SO_NONE }, 
        { OPT_usagefile,       "--usage-file",      SO_REQ_SEP }, 
        { OPT_help,            "--help",            SO_NONE },
        { OPT_help,            "-?",                SO_NONE }, 
        { OPT_list,            "--list",            SO_REQ_SEP },
        { OPT_list,            "-L",                SO_REQ_SEP },
        { OPT_output,          "--output",          SO_REQ_SEP },
        { OPT_output,          "-o",                SO_REQ_SEP },
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
        } else if (args.OptionId() == OPT_nomate) {
            opt_mate = false;
        } else if (args.OptionId() == OPT_usageonly) {
            opt_usageonly = true;
        } else if (args.OptionId() == OPT_usagefile) {
            usage_file = args.OptionArg();
        } else if (args.OptionId() == OPT_list) {
            list_file = args.OptionArg();
        } else if (args.OptionId() == OPT_output) {
            output_file = args.OptionArg();
#ifdef _WITH_DEBUG
        } else if (args.OptionId() == OPT_debug) {
            opt_debug = args.OptionArg() ? atoi(args.OptionArg()) : opt_debug;
        } else if (args.OptionId() == OPT_reads) {
            opt_reads = strtoll(args.OptionArg(), NULL, 10);
        } else if (args.OptionId() == OPT_progress) {
            opt_progress = args.OptionArg() ? strtoll(args.OptionArg(), NULL, 10) : opt_progress;
#endif
        } else {
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
        cerr << NAME << " can't currently read from stdin, ask Doug about it" << endl;
        return EXIT_FAILURE;
        input_file = "/dev/stdin";
    }

    // set up output; if file not specified, use stdout or its equivalent
    if (output_file.empty()) {
        output_file = "/dev/stdout";
    }


    //----------------- If --list option used, open file and read in list of references.

    //typedef std::tr1::unordered_map<string, int32_t> bedMap;
    typedef std::tr1::unordered_map<string, bool> nameMap;
    nameMap name_map;
    // I can do much better than this...
    if (! list_file.empty()) {
        if (opt_progress || DEBUG(1))
            cerr << NAME << "[pass1] reading reference sequence names from "
                << list_file << endl;
        ifstream list_stream(list_file.c_str());
        vector<string> fields;
        char buf[10001];
        int32_t nl = 0;
        while (list_stream.getline(buf, 10000)) {
            ++nl;
            IF_DEBUG(1) cerr << "list line " << nl << " " << buf << endl;
            // ignore lines beginning with "#"
            if (buf[0] == '#') continue;
            stringstream line_stream(buf);
            string ref = "";
            line_stream >> noskipws >> ref;
            IF_DEBUG(1) cerr << "list line " << nl << " ref :" << ref << ":" << endl;
            name_map[ref] = true;

        }
        list_stream.close();
    }


    //----------------- Open input BAM, create header for output BAM


	BamReader reader;

    if (opt_progress || DEBUG(1))
        cerr << NAME << "[pass1] opening input BAM and reading references..." << endl;

	if (! reader.Open(input_file)) {
        cerr << NAME << "[pass1] could not open BAM input" << endl;
        return EXIT_FAILURE;
    }

    if (reader.GetReferenceCount() == 0) {
        cerr << NAME << "[pass1] no reference sequences found in BAM header" << endl;
        reader.Close();
        return EXIT_FAILURE;
    }

#ifdef _IF_BAMTOOLS_IS_BROKEN
    SamHeader header = reader.GetHeader();
#else
    const SamHeader& header = reader.GetConstSamHeader();
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


    if (true || opt_progress || DEBUG(1))
        cerr << NAME << "[pass1] " << reader.GetReferenceCount() 
            << " references in the input BAM" << endl;

    vector<int64_t> refs_mentioned( reader.GetReferenceCount() );
    vector<int64_t> refs_mentioned_mate( reader.GetReferenceCount() );
    int64_t n_unref_mentioned = 0;
    int64_t n_unref_mentioned_mate = 0;

    int64_t n_reads = 0;  // number of reads processed
	BamAlignment al;  // holds the current read from the BAM file

	while (reader.GetNextAlignmentCore(al) && (opt_reads < 0 || n_reads < opt_reads)) {

        ++n_reads;
        if (al.IsMapped()) {
            if (al.RefID < 0) {
                ++n_unref_mentioned;
                // cerr << NAME << "[pass1] missing reference sequence from input bam" << endl;
                // return EXIT_FAILURE;
            } else {
                ++refs_mentioned[al.RefID];
            }
        }
        if (al.IsPaired() && al.IsMateMapped() && al.MateRefID >= 0) {
            // an unmapped mate has our RefID and Position, so not a reference "use"
            ++refs_mentioned_mate[al.MateRefID];
        } else if (al.IsPaired() && al.IsMateMapped() && al.MateRefID < 0) {
            // if a reference is missing for a mapped mate then MateRefID == -1,
            ++n_unref_mentioned_mate;
        }
        // FIXME handle at least a subset of reference mentions within tags

        if ((opt_progress || DEBUG(1)) && n_reads % opt_progress == 0)
            cerr << NAME << "[pass1] " << n_reads << " reads examined..." << endl;
 
	}
    if (opt_progress || DEBUG(1))
        cerr << NAME << "[pass1] " << n_reads << " reads examined" << endl;


    //----------------- Pass 2: Create new reference set


    const RefVector& old_refs = reader.GetReferenceData();
    RefVector        new_refs;
    int32_t          n_refs_mention = 0;
    int32_t          n_refs_mate = 0;
    int32_t          n_refs_mate_not_kept = 0;
    int32_t          n_refs_name = 0;

    vector<refStats> refs_stats; // don't allocate vector if not needed
    if (! usage_file.empty()) {
        refs_stats.resize(reader.GetReferenceCount() + 1);
        size_t i_unref = refs_stats.size() - 1;  // last entry, for mentions of ref -1
        for (size_t i = 0; i < refs_mentioned.size(); ++i) {
            refs_stats[i].ref = old_refs[i].RefName;
            refs_stats[i].old_id = i;
            refs_stats[i].m_read = refs_mentioned[i];
            refs_stats[i].m_mate = refs_mentioned_mate[i];
            refs_stats[i].m_name = ! name_map.empty() && name_map.find(refs_stats[i].ref) != name_map.end();
            refs_stats[i].no_mate = false;
            refs_stats[i].new_id = -1;
        }
        refs_stats[i_unref].ref = "*";
        refs_stats[i_unref].old_id = -1;
        refs_stats[i_unref].m_read = n_unref_mentioned;
        refs_stats[i_unref].m_mate = n_unref_mentioned_mate;
        refs_stats[i_unref].m_name = false;
        refs_stats[i_unref].no_mate = false;
        refs_stats[i_unref].new_id = -1;
    }

    assert(new_header.Sequences.IsEmpty());  // new_refs contains the new @SQ info
    int32_t new_RefID = 0;
    for (size_t i = 0; i < refs_mentioned.size(); ++i) {

        if (   (refs_mentioned[i] > 0)  // if any of the reasons for keeping it are true
            || (opt_mate && refs_mentioned_mate[i] > 0)
            || (! name_map.empty() && name_map.find(old_refs[i].RefName) != name_map.end())) {

            if (refs_mentioned[i] > 0) {
                ++n_refs_mention;
            } else if (opt_mate && refs_mentioned_mate[i] > 0) {
                ++n_refs_mate;
            } else if (! name_map.empty() && name_map.find(old_refs[i].RefName) != name_map.end()) {
                ++n_refs_name;
            }

            new_refs.push_back(old_refs[i]);
#if defined(_BAMTOOLS_EXTENSION)  && ! defined(_IF_BAMTOOLS_IS_BROKEN)
            SamSequenceConstIterator refI = header.Sequences.ConstFind(old_refs[i].RefName);
            if (refI == header.Sequences.ConstEnd()) {
                cerr << NAME << "[pass2] internal error, " << old_refs[i].RefName 
                    << " not found in the input header" << endl;
                return EXIT_FAILURE;
            }
            const SamSequence& existing_ref = (*refI);
#else
            const SamSequence& existing_ref = header.Sequences[old_refs[i].RefName];
#endif
            new_header.Sequences.Add(existing_ref);
            refs_mentioned[i] = new_RefID;  // entry now contains new reference ID
            ++new_RefID;

        } else {
            if (refs_mentioned_mate[i]) {
                ++n_refs_mate_not_kept;
                if (! usage_file.empty()) 
                    refs_stats[i].no_mate = true;
            }
            refs_mentioned[i] = -1;
        }
        if (! usage_file.empty())
            refs_stats[i].new_id = refs_mentioned[i];
    }
    assert(new_refs.size() == (size_t)new_RefID);

    if (true || opt_progress || DEBUG(1)) {
        cerr << NAME << "[pass2] " << new_RefID 
            << " references kept in the output BAM" << endl;
        if (n_refs_mention)
            cerr << NAME << "[pass2] " << n_refs_mention
                << " references were mentioned by mapped reads" << endl;
        if (n_refs_mate)
            cerr << NAME << "[pass2] " << n_refs_mate
                << " additional references were mentioned by mapped mates" << endl;
        if (! list_file.empty() || n_refs_name) {
            cerr << NAME << "[pass2] " << n_refs_name 
                << " additional references were kept by name";
            if (! list_file.empty())
                cerr << " (out of " << name_map.size() << " named in " << list_file << ")";
            cerr << endl;
        }
        if (n_refs_mate_not_kept)
            cerr << NAME << "[pass2] " << n_refs_mate_not_kept
                << " references were mentioned by mapped mates but not kept (--no-mate)" << endl;
    }

    IF_DEBUG(2) {
        for (size_t i = 0; i < new_refs.size(); ++i) {
            cerr << NAME << "[pass2] " << i << "] SN:" << new_refs[i].RefName
                << "  LN:" << new_refs[i].RefLength << endl;
        }
        IF_DEBUG(3) {
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

    if (! usage_file.empty()) {
        ofstream usage_stream(usage_file.c_str());
        usage_stream << "ref" << sep << "input_id" << sep << "m_read" << sep << "m_mate" 
            << sep << "m_name" << sep << "no_mate" << sep << "output_id" << endl;
        for (size_t i = 0; i < refs_stats.size(); ++i) {
            usage_stream << refs_stats[i].ref;
            usage_stream << sep << refs_stats[i].old_id;
            usage_stream << sep << refs_stats[i].m_read;
            usage_stream << sep << refs_stats[i].m_mate;
            usage_stream << sep << refs_stats[i].m_name;
            usage_stream << sep << refs_stats[i].no_mate;
            usage_stream << sep << refs_stats[i].new_id;
            usage_stream << endl;
        }
        usage_stream.close();
        cerr << NAME << " per-reference usage in " << usage_file << " (--usage-file)" << endl;
    }

    if (opt_usageonly) {
	    reader.Close();
        cerr << NAME << " no output BAM produced (--usage-only)" << endl;
	    return EXIT_SUCCESS;
    }


    //----------------- Pass 2: Second pass through reads, write new BAM file


    BamWriter writer;

    IF_DEBUG(2) {
        cerr << "********* BEGIN new_header.ToString()" << endl;
        cerr << new_header.ToString();
        cerr << "********* END   new_header.ToString()" << endl;
    }


    if (! writer.Open(output_file, new_header, new_refs)) {
        cerr << NAME << " could not open output " << output_file << endl;
        return EXIT_FAILURE;
    }

    int64_t n_reads_pass1 = n_reads;
    n_reads = 0;
    int64_t n_reads_rerefd = 0;  // number of reads given re-references
    int64_t n_mates_derefd = 0;  // number of mates that have references removed

    reader.Rewind();

	while (reader.GetNextAlignmentCore(al) && (opt_reads < 0 || n_reads < opt_reads)) {

        ++n_reads;

        if (al.IsMapped()) {
            assert(al.RefID >= 0);  // it was valid before...
            if (al.RefID != refs_mentioned[al.RefID]) {  // the reference ID is different
                ++n_reads_rerefd;  // strictly rereferenced
                if (al.IsPaired()) {
                    if (al.MateRefID == al.RefID
                        || refs_mentioned[al.MateRefID] >= 0) {
                        // update the mate RefID
                        al.MateRefID = refs_mentioned[al.MateRefID];
                    } else if (al.IsMateMapped()) {
                        // mate ref is now unavailable
                        al.MateRefID = -1;
                        ++n_mates_derefd;
                    }
                }
                al.RefID = refs_mentioned[al.RefID];
                assert(al.RefID >= 0);  // and it is valid after
            }
        }

        writer.SaveAlignment(al);

        if (opt_progress && n_reads % opt_progress == 0) {
            cerr << NAME << "[pass2] " << n_reads << " reads rereferenced";
            if (! opt_mate)
                cerr << ", "<< n_mates_derefd << " mates dereferenced";
            cerr << "..." << endl;
        }
 
	}
    if (true || opt_progress || DEBUG(1)) {
        cerr << NAME << "[pass2] " << n_reads << " reads rereferenced";
        if (! opt_mate)
            cerr << ", "<< n_mates_derefd << " mates dereferenced";
        cerr << endl;
    }
    assert(n_reads == n_reads_pass1);

	reader.Close();
	writer.Close();

	return EXIT_SUCCESS;
}

