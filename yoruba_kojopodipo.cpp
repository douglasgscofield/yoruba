// yoruba_kojopodipo.cpp  (c) Douglas G. Scofield, Dept. Plant Physiology, Umeå University
//
// Kojopodipo (English command is readgroup) adds or modifies the read group in a BAM file.
//
// Kojopodipo rewrites the header and then annotates each read.  I've started this because
// picard's AddOrChangeReadGroups.jar uses a lot of memory and seems to take a long time for
// what should be a simple task.  We'll see if those opinions hold up :-)
//
// Kojopodipo is the Yoruba (Nigeria) verb for 'to group'.
//
// Uses BamTools C++ API for handling BAM files


// CHANGELOG
//
//
//
// TODO
//
// Write usage
// Write README.md
//

#include "yoruba_kojopodipo.h"

using namespace std;
using namespace BamTools;
using namespace yoruba;

// options
static string       input_file;  // defaults to stdin, set from command line
static string       output_file;  // defaults to stdout, set with -o FILE
static bool         opt_noreplace = false;
static bool         opt_onlyreplace = false;
static bool         opt_clear = false;
// leave debug options in
#ifdef _WITH_DEBUG
#endif
static int32_t      opt_debug = 0;
static int64_t      opt_reads = -1;
static int64_t      opt_progress = 100000;
static int64_t      debug_reads_to_report = 1;


//-------------------------------------


#ifdef _STANDALONE
int 
main(int argc, char* argv[]) {
    return main_kojopodipo(argc, argv);
}
#endif


//-------------------------------------


static int
usage(bool long_help = false)
{
    cerr << endl;
    cerr << "Usage:   " << YORUBA_NAME << " readgroup [options] <in.bam>" << endl;
    cerr << "         " << YORUBA_NAME << " kojopodipo [options] <in.bam>" << endl;
    cerr << endl;
    cerr << "Either command invokes this function." << endl;
    cerr << endl;
    if (long_help) {
    cerr << "\
Add or replace read group information in the BAM file <in.bam>.\n\
\n\
Read group information appears in two places in the BAM file:\n\
   (1) the read group dictionary, found in the header, which contains the \n\
       read group ID and any other information associated with the ID, such \n\
       as library, sample name, etc., and thus defines the read group;\n\
   (2) the RG tag on each read, which specifies one of the IDs that appear\n\
       in the read group dictionary, and thus declares the read to be part\n\
       of the given read group.\n\
\n\
By default, all reads in the BAM file will be given the supplied read group.\n\
If the dictionary already defines a read group with the same ID, its definition\n\
will be replaced with the supplied information.  If the dictionary contains\n\
other read groups, their definitions will remain in the BAM file header but\n\
all reads will be given the supplied read group.\n\
\n\
Other behaviour can be specified using --no-replace, --only-replace and \n\
--clear.  See table below.\n\
\n\
The only argument required to specify a valid read group is --ID or --id.\n\
\n";
    }
    cerr << "Options: --ID STR | --id STR                 read group identifier" << endl;
    cerr << "         --LB STR | --library STR            library" << endl;
    cerr << "         --SM STR | --sample-name STR        sample name" << endl;
    cerr << "         --DS STR | --description STR        description" << endl;
    cerr << "         --DT STR | --date STR               date" << endl;
    cerr << "         --PG STR | --programs STR           programs used" << endl;
    cerr << "         --PL STR | --platform STR           sequencing platform" << endl;
    cerr << "         --PU STR | --platform-unit STR      platform unit" << endl;
    cerr << "         --PI STR | --predicted-insert STR   predicted median insert size" << endl;
    cerr << "         --FO STR | --flow-order STR         flow order" << endl;
    cerr << "         --KS STR | --key-sequence STR       key sequence" << endl;
    cerr << "         --CN STR | --sequencing-center STR  sequencing center" << endl;
    cerr << "         --o FILE | -o FILE | --output FILE  output file name [default is stdout]" << endl;
    cerr << "         --no-replace                        abort if the read group exists" << endl;
    cerr << "         --only-replace                      replace just this read group" << endl;
    cerr << "         --clear                             clear all read group information" << endl;
    cerr << "         --? | -? | --help                   longer help" << endl;
    cerr << endl;
#ifdef _WITH_DEBUG
    cerr << "         --debug INT     debug info level INT [" << opt_debug << "]" << endl;
    cerr << "         --reads INT     process at most this many reads [" << opt_reads << "]" << endl;
    cerr << "         --progress INT  print reads processed mod INT [" << opt_progress << "]" << endl;
    cerr << endl;
#endif
    if (long_help) {
    cerr << "\
No formatting restrictions are imposed on any of the read group elements. It\n\
is the user's responsibility to ensure that they conform to the SAM definitions\n\
<http://samtools.sourceforge.net/SAM1.pdf> or to any other tool requirements.\n\
\n\
If the output file is not specified, then output is written to stdout.\n\
\n\
The --no-replace option will abort if the given read group ID is found in the\n\
dictionary, and will only add read group information to reads that don't\n\
don't already have it.\n\
\n\
The --only-replace option modifies information for only those reads in the\n\
supplied read group (same ID). Read group information for other reads,\n\
including those without any other read group information, is unchanged.\n\
\n\
The --clear option removes all read group information from all reads.  If\n\
specified with options defining a read group, then the read group dictionary\n\
will be cleared prior to defining the new read group.\n\
\n\
Only one of these may be supplied at a time.  To summarizing the effects of these options:\n\
\n\
                      Read read group (RG) tag status                          \n\
                ---------------------------------------------                  \n\
                    no RG    |    RG matches   |  RG does not                  \n\
Option                       |       --ID      |  match --ID    RG dictionary  \n\
--------------  ---------------------------------------------  ----------------\n\
                                                                               \n\
only --ID etc.             new RG set for all reads             RG added       \n\
                                                                               \n\
--no-replace     new RG set         abort         no change     RG added; abort\n\
                 from --ID                                      if present     \n\
                                                                               \n\
--only-replace   no change        no change       no change     RG updated     \n\
                                                                from options   \n\
                                                                               \n\
--clear          no change        RG removed      RG removed    cleared        \n\
  no --ID                                                                      \n\
                                                                               \n\
--clear                    new RG set for all reads             cleared, then  \n\
  with --ID                                                     RG added       \n\
\n\
\n";
    cerr << "Kojopodipo is the Yoruba (Nigeria) verb for 'to group'." << endl;
    cerr << endl;
    }

    return 1;
}


//-------------------------------------


int 
yoruba::main_kojopodipo(int argc, char* argv[])
{
    SamReadGroup new_readgroup;  // the read group we are creating
    SamProgram   new_program;    // the program info for yoruba, added to the header

    // first process options

	if( argc < 2 ) {
		return usage();
	}

    enum { OPT_ID, OPT_LB, OPT_SM, OPT_DS, OPT_DT, OPT_PG, OPT_PL, OPT_PU, OPT_PI, OPT_FO,
        OPT_KS, OPT_CN, OPT_output, OPT_noreplace, OPT_onlyreplace, OPT_clear,
#ifdef _WITH_DEBUG
        OPT_debug, OPT_reads, OPT_progress,
#endif
        OPT_help };

    CSimpleOpt::SOption kojopodipo_options[] = {
        { OPT_ID, "--ID", SO_REQ_SEP }, { OPT_ID, "--id", SO_REQ_SEP },
        { OPT_LB, "--LB", SO_REQ_SEP }, { OPT_LB, "--library", SO_REQ_SEP },
        { OPT_SM, "--SM", SO_REQ_SEP }, { OPT_SM, "--sample-name", SO_REQ_SEP },
        { OPT_DS, "--DS", SO_REQ_SEP }, { OPT_DS, "--description", SO_REQ_SEP },
        { OPT_DT, "--DT", SO_REQ_SEP }, { OPT_DT, "--date", SO_REQ_SEP },
        { OPT_PG, "--PG", SO_REQ_SEP }, { OPT_PG, "--programs", SO_REQ_SEP },
        { OPT_PL, "--PL", SO_REQ_SEP }, { OPT_PL, "--platform", SO_REQ_SEP },
        { OPT_PU, "--PU", SO_REQ_SEP }, { OPT_PU, "--platform-unit", SO_REQ_SEP },
        { OPT_PI, "--PI", SO_REQ_SEP }, { OPT_PI, "--predicted-insert", SO_REQ_SEP },
        { OPT_FO, "--FO", SO_REQ_SEP }, { OPT_FO, "--flow-order", SO_REQ_SEP },
        { OPT_KS, "--KS", SO_REQ_SEP }, { OPT_KS, "--key-sequence", SO_REQ_SEP },
        { OPT_CN, "--CN", SO_REQ_SEP }, { OPT_CN, "--sequencing-center", SO_REQ_SEP },
        { OPT_output, "--o", SO_REQ_SEP }, 
        { OPT_output, "-o", SO_REQ_SEP }, 
        { OPT_output, "--output", SO_REQ_SEP },
        { OPT_noreplace, "--no-replace", SO_NONE },
        { OPT_onlyreplace, "--only-replace", SO_NONE },
        { OPT_clear, "--clear", SO_NONE },
        { OPT_help, "--?", SO_NONE }, 
        { OPT_help, "-?", SO_NONE }, 
        { OPT_help, "--help", SO_NONE },
// leave debug options in
#ifdef _WITH_DEBUG
#endif
        { OPT_debug, "--debug", SO_REQ_SEP },
        { OPT_reads, "--reads", SO_REQ_SEP },
        { OPT_progress, "--progress", SO_REQ_SEP },
// end debug options
        SO_END_OF_OPTIONS
    };

    CSimpleOpt args(argc, argv, kojopodipo_options);

    while (args.Next()) {
        if (args.LastError() != SO_SUCCESS) {
            cerr << NAME << " invalid argument '" << args.OptionText() << "'" << endl;
            return usage();
        }
        if (args.OptionId() == OPT_help) return usage(true);
        else if (args.OptionId() == OPT_ID) new_readgroup.ID = args.OptionArg();
        else if (args.OptionId() == OPT_LB) new_readgroup.Library = args.OptionArg();
        else if (args.OptionId() == OPT_SM) new_readgroup.Sample = args.OptionArg();
        else if (args.OptionId() == OPT_DS) new_readgroup.Description = args.OptionArg();
        else if (args.OptionId() == OPT_DT) new_readgroup.ProductionDate = args.OptionArg();
        else if (args.OptionId() == OPT_PG) new_readgroup.Program = args.OptionArg();
        else if (args.OptionId() == OPT_PL) new_readgroup.SequencingTechnology = args.OptionArg();
        else if (args.OptionId() == OPT_PU) new_readgroup.PlatformUnit = args.OptionArg();
        else if (args.OptionId() == OPT_PI) new_readgroup.PredictedInsertSize = args.OptionArg();
        else if (args.OptionId() == OPT_FO) new_readgroup.FlowOrder = args.OptionArg();
        else if (args.OptionId() == OPT_KS) new_readgroup.KeySequence = args.OptionArg();
        else if (args.OptionId() == OPT_CN) new_readgroup.SequencingCenter = args.OptionArg();
        else if (args.OptionId() == OPT_output) output_file = args.OptionArg();
        else if (args.OptionId() == OPT_noreplace) opt_noreplace = true;
        else if (args.OptionId() == OPT_onlyreplace) opt_onlyreplace = true;
        else if (args.OptionId() == OPT_clear) opt_clear = true;
// leave debug options in 
#ifdef _WITH_DEBUG
#endif
        else if (args.OptionId() == OPT_debug) 
            opt_debug = args.OptionArg() ? atoi(args.OptionArg()) : opt_debug;
        else if (args.OptionId() == OPT_reads) 
            opt_reads = strtoll(args.OptionArg(), NULL, 10);
        else if (args.OptionId() == OPT_progress) 
            opt_progress = args.OptionArg() ? strtoll(args.OptionArg(), NULL, 10) : opt_progress;
// end debug options
        else {
            cerr << NAME << " unprocessed argument '" << args.OptionText() << "'" << endl;
            return 1;
        }
    }

    // set up input location; if file not specified, use /dev/stdin
    _DEBUG(0) {
        for (int i = 0; i < args.FileCount(); ++i)
            cerr << NAME << " file argument " << i << ": " << args.File(i) << endl;
    }
    if (args.FileCount() > 1) {
        cerr << NAME << " requires at most one BAM file specified as input" << endl;
        return usage();
    } else if (args.FileCount() == 1) {
        input_file = args.File(0);
    } else if (input_file.empty()) {  // don't replace if not empty, a defauult is set
        input_file = "/dev/stdin";
    }

    // set up output; if file not specified, use /dev/stdout
    if (output_file.empty()) {
        output_file = "/dev/stdout";
    }

    // check option semantics
    if (! opt_clear && new_readgroup.ID.empty()) {
        cerr << NAME << " must define a read group using --ID or --id" << endl;
        return usage();
    }
    if (opt_noreplace + opt_onlyreplace + opt_clear > 1) {
        cerr << NAME << " use only one of --no-replace, --only-replace or --clear" << endl;
        return usage(true);
    }

	BamReader reader;

	if (! reader.Open(input_file)) {
        cerr << NAME << "could not open BAM input" << endl;
        return 1;
    }

    SamHeader header = reader.GetHeader();

    _DEBUG(0) cerr << NAME << " finished reading header from " << input_file << endl;

    new_program.ID = YORUBA_NAME;
    new_program.ID = new_program.ID + " " + argv[0];
    new_program.Name = YORUBA_NAME;
    new_program.Version = YORUBA_VERSION;
    new_program.CommandLine = YORUBA_NAME;
    for (int i = 0; i < argc; ++i)
        new_program.CommandLine = new_program.CommandLine + " " + argv[i];

    _DEBUG(0) { 
        if (opt_reads >= 0) 
            cerr << NAME << " modifying up to " << opt_reads << " reads" << endl; 
        else
            cerr << NAME << " processing all reads" << endl; 
        cerr << NAME << " opt_noreplace = " << opt_noreplace << endl;
        cerr << NAME << " opt_onlyreplace = " << opt_onlyreplace << endl;
        cerr << NAME << " opt_clear = " << opt_clear << endl;
        cerr << NAME << " new_readgroup.ID = " << new_readgroup.ID << endl;
        cerr << NAME << " new_readgroup.Library = " << new_readgroup.Library << endl;
        cerr << NAME << " new_readgroup.Sample = " << new_readgroup.Sample << endl;
        cerr << NAME << " new_readgroup.Description = " << new_readgroup.Description << endl;
        cerr << NAME << " new_readgroup.ProductionDate = " << new_readgroup.ProductionDate << endl;
        cerr << NAME << " new_readgroup.Program = " << new_readgroup.Program << endl;
        cerr << NAME << " new_readgroup.SequencingTechnology = " << new_readgroup.SequencingTechnology << endl;
        cerr << NAME << " new_readgroup.PlatformUnit = " << new_readgroup.PlatformUnit << endl;
        cerr << NAME << " new_readgroup.PredictedInsertSize = " << new_readgroup.PredictedInsertSize << endl;
        cerr << NAME << " new_readgroup.FlowOrder = " << new_readgroup.FlowOrder << endl;
        cerr << NAME << " new_readgroup.KeySequence = " << new_readgroup.KeySequence << endl;
        cerr << NAME << " new_readgroup.SequencingCenter = " << new_readgroup.SequencingCenter << endl;
        cerr << NAME << " new_program.ID = " << new_program.ID << endl;
        cerr << NAME << " new_program.Name = " << new_program.Name << endl;
        cerr << NAME << " new_program.Version = " << new_program.Version << endl;
        cerr << NAME << " new_program.CommandLine = " << new_program.CommandLine << endl;
    }

    if (header.HasReadGroups()) {

        _DEBUG(0) printReadGroupDictionary(cerr, header.ReadGroups);

        if (opt_clear)
            header.ReadGroups.Clear();

    }

    if (header.ReadGroups.Contains(new_readgroup.ID)) {
        if (opt_noreplace) {
            cerr << NAME << " BAM already contains read group '" << new_readgroup.ID << "', aborting" << endl;
            return 1;
        }
        header.ReadGroups.Remove(new_readgroup.ID);
    }

    header.ReadGroups.Add(new_readgroup);

    _DEBUG(0) printReadGroupDictionary(cerr, header.ReadGroups);

    header.Programs.Add(new_program);
	
    BamWriter writer;

    if (! writer.Open(output_file, header, reader.GetReferenceData())) {
        cerr << NAME << " could not open output " << output_file << endl;
        return 1 ;
    }
    _DEBUG(0) cerr << NAME << " writing output to " << output_file << endl;

	BamAlignment al;  // holds the current read from the BAM file

    int64_t n_reads = 0;  // number of reads processed

	while (reader.GetNextAlignmentCore(al) && (opt_reads < 0 || n_reads < opt_reads)) {

        // Nicely, BamTools supports 'lazy' evaluation of character data in alignments.
        // Since we don't need character data from each read with the --only-replace
        // option, use GetNextAlignmentCore() and allow character data to be built
        // only when it's required.  This saves on the order of 18% of execution time
        // on my 64-bit iMac when --only-replace is used.

        ++n_reads;

        _DEBUG(1) printAlignmentInfo(cerr, al);

        string RG_tag;

        if (! opt_onlyreplace) { 

            if (al.GetTag("RG", RG_tag) && ! opt_noreplace) {

                _DEBUG(0) if (n_reads <= debug_reads_to_report) {
                    cerr << NAME << " " << n_reads << " read before processing: ";
                    printAlignmentInfo(cerr, al);
                }
                _DEBUG(1) cerr << NAME << " " << al.Name << " has tag: RG:Z:'" << RG_tag << "'" << endl;

                if (! al.EditTag("RG", "Z", new_readgroup.ID)) {
                    cerr << NAME << " could not edit tag for read " << al.Name << endl;
                    return 1;
                }

                _DEBUG(0) if (n_reads <= debug_reads_to_report) {
                    cerr << NAME << " " << n_reads << " read after processing: ";
                    printAlignmentInfo(cerr, al);
                }

            } else {

                _DEBUG(0) if (n_reads <= debug_reads_to_report) {
                    cerr << NAME << " " << n_reads << " read before processing: ";
                    printAlignmentInfo(cerr, al);
                }

                if (! al.AddTag("RG", "Z", new_readgroup.ID)) {
                    cerr << NAME << " could not add tag to read " << al.Name << endl;
                    return 1;
                }

                _DEBUG(0) if (n_reads <= debug_reads_to_report) {
                    cerr << NAME << " " << n_reads << " read after processing: ";
                    printAlignmentInfo(cerr, al);
                }

            }

        }

        _DEBUG(1) { printAlignmentInfo(cerr, al); cerr << endl; }

        writer.SaveAlignment(al);

        if (opt_progress && n_reads % opt_progress == 0)
            cerr << NAME << " " << n_reads << " reads processed..." << endl;
	}

    _DEBUG(0) cerr << NAME << " " << n_reads << " reads processed" << endl;

	reader.Close();
	writer.Close();

	return 0;
}

