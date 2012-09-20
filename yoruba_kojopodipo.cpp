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
static bool         other_rg_opts = false;  // read group options other than --ID were given
static bool         opt_dictionary; 
static string       dictionary_string; 
static bool         opt_replace;
static string       replace_string;
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
other read groups, their definitions will remain in the BAM file header (if\n\
present) but all reads will be given the supplied read group.\n\
\n\
Other behaviour can be specified using --replace and --clear.  See table below.\n\
\n\
The only argument required to specify a valid read group is --ID or --id.\n\
\n";
    }
    cerr << "Options: --ID STR | --id STR                 read group identifier" << endl;
    cerr << "         --LB STR | --library STR            read group library" << endl;
    cerr << "         --SM STR | --sample-name STR        read group sample name" << endl;
    cerr << "         --DS STR | --description STR        read group description" << endl;
    cerr << "         --DT STR | --date STR               read group date" << endl;
    cerr << "         --PG STR | --programs STR           read group programs used" << endl;
    cerr << "         --PL STR | --platform STR           read group sequencing platform" << endl;
    cerr << "         --PU STR | --platform-unit STR      read group platform unit" << endl;
    cerr << "         --PI STR | --predicted-insert STR   read group predicted median insert size" << endl;
    cerr << "         --FO STR | --flow-order STR         read group flow order" << endl;
    cerr << "         --KS STR | --key-sequence STR       read group key sequence" << endl;
    cerr << "         --CN STR | --sequencing-center STR  read group sequencing center" << endl;
    cerr << endl;
    cerr << "         --o FILE | -o FILE | --output FILE  output file name [default is stdout]" << endl;
    cerr << "         --replace STR                       replace read group STR with --ID" << endl;
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
The --replace option will replace the identified read group to have the name\n\
provided in --ID, in both its dictionary entry and on reads.  If only --ID\n\
is provided, then the read group is simply renamed.  If any other read group\n\
options are given, then the read group is redefined as well.\n\
\n\
The --clear option removes all read group information from all reads.  If\n\
specified with options defining a read group, then the read group dictionary\n\
will be cleared prior to defining the new read group.\n\
\n\
Only one of these may be supplied at a time.  To summarizing the effects of these options:\n\
\n\
                      Read read group (RG) tag status                        \n\
                --------------------------------------------                 \n\
                  no RG    |    RG matches    |  RG does not                 \n\
Option                     |       STR        |   match STR   RG dictionary  \n\
--------------  --------------------------------------------  ---------------\n\
                                                                             \n\
only --ID etc.           new RG set for all reads             RG added       \n\
                                                                             \n\
--replace STR    no change  RG changed to --ID   no change    RG STR updated with\n\
                                  to --ID                     --ID; entry replaced\n\
                                                              if any other RG\n\
                                                              options\n\
                                                                             \n\
--clear          no change       RG removed      RG removed   cleared        \n\
  no --ID                                                                    \n\
                                                                             \n\
--clear                  new RG set for all reads             cleared, then  \n\
  with --ID                                                   RG added       \n\
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
    SamReadGroup new_rg;  // the read group we are creating
    SamProgram   new_program;    // the program info for yoruba, added to the header

    // first process options

	if( argc < 2 ) {
		return usage();
	}

    new_program.ID = YORUBA_NAME;
    new_program.ID = new_program.ID + " " + argv[0];
    new_program.Name = YORUBA_NAME;
    new_program.Version = YORUBA_VERSION;
    new_program.CommandLine = YORUBA_NAME;
    for (int i = 0; i < argc; ++i)
        new_program.CommandLine = new_program.CommandLine + " " + argv[i];

    enum { OPT_ID, OPT_LB, OPT_SM, OPT_DS, OPT_DT, OPT_PG, OPT_PL, OPT_PU, OPT_PI, OPT_FO,
        OPT_KS, OPT_CN, OPT_dictionary, OPT_output, OPT_replace, OPT_clear,
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
        { OPT_output,      "--o", SO_REQ_SEP }, 
        { OPT_output,      "-o", SO_REQ_SEP }, 
        { OPT_output,      "--output", SO_REQ_SEP },
        { OPT_dictionary,  "--dictionary", SO_REQ_SEP },
        { OPT_replace,     "--replace", SO_REQ_SEP },
        { OPT_clear,       "--clear", SO_NONE },
        { OPT_help,        "--?", SO_NONE }, 
        { OPT_help,        "-?", SO_NONE }, 
        { OPT_help,        "--help", SO_NONE },
// leave debug options in
#ifdef _WITH_DEBUG
#endif
        { OPT_debug,       "--debug", SO_REQ_SEP },
        { OPT_reads,       "--reads", SO_REQ_SEP },
        { OPT_progress,    "--progress", SO_REQ_SEP },
// end debug options
        SO_END_OF_OPTIONS
    };

    CSimpleOpt args(argc, argv, kojopodipo_options);

    while (args.Next()) {
        if (args.LastError() != SO_SUCCESS) {
            cerr << NAME << " invalid argument '" << args.OptionText() << "'" << endl;
            return usage();
        }
        if (args.OptionId() == OPT_help) {
            return usage(true);
        } else if (args.OptionId() == OPT_ID) {
            new_rg.ID = args.OptionArg();
        } else if (args.OptionId() == OPT_LB) {
            new_rg.Library = args.OptionArg(); other_rg_opts = true;
        } else if (args.OptionId() == OPT_SM) {
            new_rg.Sample = args.OptionArg(); other_rg_opts = true;
        } else if (args.OptionId() == OPT_DS) {
            new_rg.Description = args.OptionArg(); other_rg_opts = true;
        } else if (args.OptionId() == OPT_DT) {
            new_rg.ProductionDate = args.OptionArg(); other_rg_opts = true;
        } else if (args.OptionId() == OPT_PG) {
            new_rg.Program = args.OptionArg(); other_rg_opts = true;
        } else if (args.OptionId() == OPT_PL) {
            new_rg.SequencingTechnology = args.OptionArg(); other_rg_opts = true;
        } else if (args.OptionId() == OPT_PU) {
            new_rg.PlatformUnit = args.OptionArg(); other_rg_opts = true;
        } else if (args.OptionId() == OPT_PI) {
            new_rg.PredictedInsertSize = args.OptionArg(); other_rg_opts = true;
        } else if (args.OptionId() == OPT_FO) {
            new_rg.FlowOrder = args.OptionArg(); other_rg_opts = true;
        } else if (args.OptionId() == OPT_KS) {
            new_rg.KeySequence = args.OptionArg(); other_rg_opts = true;
        } else if (args.OptionId() == OPT_CN) {
            new_rg.SequencingCenter = args.OptionArg(); other_rg_opts = true;
        } else if (args.OptionId() == OPT_output) {
            output_file = args.OptionArg();
        } else if (args.OptionId() == OPT_dictionary) {
            opt_dictionary = true; dictionary_string = args.OptionArg();
        } else if (args.OptionId() == OPT_replace) {
            opt_replace = true; replace_string = args.OptionArg();
        } else if (args.OptionId() == OPT_clear) {
            opt_clear = true;
// leave debug options in 
#ifdef _WITH_DEBUG
#endif
        } else if (args.OptionId() == OPT_debug) {
            opt_debug = args.OptionArg() ? atoi(args.OptionArg()) : opt_debug;
        } else if (args.OptionId() == OPT_reads) {
            opt_reads = strtoll(args.OptionArg(), NULL, 10);
        } else if (args.OptionId() == OPT_progress) {
            opt_progress = args.OptionArg() ? strtoll(args.OptionArg(), NULL, 10) : opt_progress;
// end debug options
        } else {
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
    if (! opt_clear && ! opt_dictionary && new_rg.ID.empty()) {
        cerr << NAME << " must define a read group using --ID or --id" << endl;
        return usage();
    }
    if (opt_replace + opt_clear > 1) {
        cerr << NAME << " use only one of --replace or --clear" << endl;
        return usage(true);
    }

	BamReader reader;

	if (! reader.Open(input_file)) {
        cerr << NAME << "could not open BAM input" << endl;
        return 1;
    }

    SamHeader header = reader.GetHeader();

    _DEBUG(0) { 
        if (opt_reads >= 0) 
            cerr << NAME << " modifying up to " << opt_reads << " reads" << endl; 
        else
            cerr << NAME << " processing all reads" << endl; 
        cerr << NAME << " opt_dictionary = " << opt_dictionary << endl;
        cerr << NAME << " dictionary_string = " << dictionary_string << endl;
        cerr << NAME << " opt_replace = " << opt_replace << endl;
        cerr << NAME << " replace_string = " << replace_string << endl;
        cerr << NAME << " opt_clear = " << opt_clear << endl;
        cerr << NAME << " new_rg.ID = " << new_rg.ID << endl;
        cerr << NAME << " new_rg.Library = " << new_rg.Library << endl;
        cerr << NAME << " new_rg.Sample = " << new_rg.Sample << endl;
        cerr << NAME << " new_rg.Description = " << new_rg.Description << endl;
        cerr << NAME << " new_rg.ProductionDate = " << new_rg.ProductionDate << endl;
        cerr << NAME << " new_rg.Program = " << new_rg.Program << endl;
        cerr << NAME << " new_rg.SequencingTechnology = " << new_rg.SequencingTechnology << endl;
        cerr << NAME << " new_rg.PlatformUnit = " << new_rg.PlatformUnit << endl;
        cerr << NAME << " new_rg.PredictedInsertSize = " << new_rg.PredictedInsertSize << endl;
        cerr << NAME << " new_rg.FlowOrder = " << new_rg.FlowOrder << endl;
        cerr << NAME << " new_rg.KeySequence = " << new_rg.KeySequence << endl;
        cerr << NAME << " new_rg.SequencingCenter = " << new_rg.SequencingCenter << endl;
        cerr << NAME << " new_program.ID = " << new_program.ID << endl;
        cerr << NAME << " new_program.Name = " << new_program.Name << endl;
        cerr << NAME << " new_program.Version = " << new_program.Version << endl;
        cerr << NAME << " new_program.CommandLine = " << new_program.CommandLine << endl;
    }

    //-------------------------------------  @RG: read group dictionary
    // the read group dictionary may not exist though there are RG tags on reads

    _DEBUG(0) cerr << NAME << " read group dictionary before modifying it" << endl;
    _DEBUG(0) printReadGroupDictionary(cerr, header.ReadGroups);

    if (opt_clear) {
        if (header.HasReadGroups())
            header.ReadGroups.Clear();
    }

    if (opt_dictionary) {
        if (header.HasReadGroups())
            header.ReadGroups.Clear();
        SamReadGroupDictionary rgd = parseReadGroupDictionaryString(dictionary_string);
        if (rgd.IsEmpty()) {
            cerr << NAME << " error parsing read group dictionary" << endl;
            return 1;
        }
        header.ReadGroups.Add(rgd);
        _DEBUG(0) cerr << NAME << " dictionary after adding parsed string '" << dictionary_string << "'" << endl;
        _DEBUG(0) printReadGroupDictionary(cerr, header.ReadGroups);
    }

    if (opt_replace) {  // --replace was given
        if (header.ReadGroups.Contains(replace_string)) {
            if (other_rg_opts) {  // more than --ID was given, replace entry and add new one
                header.ReadGroups.Remove(replace_string);
                if (header.ReadGroups.Contains(new_rg.ID))  // remove entry for new name, if it exists
                    header.ReadGroups.Remove(new_rg.ID);
                header.ReadGroups.Add(new_rg);
            } else {  // only --ID given, so simply rename
                header.ReadGroups[replace_string].ID = new_rg.ID;
            }
        } else {
            header.ReadGroups.Add(new_rg.ID);
        }
    }  else {
        if (header.ReadGroups.Contains(new_rg.ID))
            header.ReadGroups.Remove(new_rg.ID);
        header.ReadGroups.Add(new_rg.ID);
    }

    _DEBUG(0) cerr << NAME << " read group dictionary after modifying it" << endl;
    _DEBUG(0) printReadGroupDictionary(cerr, header.ReadGroups);

    //-------------------------------------  @PG: programs

    if (header.Programs.Contains(new_program.ID)) {
        // I would prefer to use duplicate program IDs in the SAM header,
        // if the same program worked over the file twice or more
        SamProgram& prog = header.Programs[new_program.ID];
        prog.Name = new_program.Name;
        prog.Version = new_program.Version;
        prog.CommandLine = new_program.CommandLine;
    } else {
        header.Programs.Add(new_program);
    }
	
    //-------------------------------------  open output

    BamWriter writer;

    if (! writer.Open(output_file, header, reader.GetReferenceData())) {
        cerr << NAME << " could not open output " << output_file << endl;
        return 1 ;
    }

	BamAlignment al;  // holds the current read from the BAM file
    int64_t n_reads = 0;  // number of reads processed

    //-------------------------------------  loop through reads in BAM file

	while (reader.GetNextAlignment(al) && (opt_reads < 0 || n_reads < opt_reads)) {

        ++n_reads;

        _DEBUG(0) if (n_reads <= debug_reads_to_report) {
            cerr << NAME << " " << n_reads << " read before processing: ";
            printAlignmentInfo(cerr, al);
        }

        string RG_tag;

        if (opt_clear) {
            al.RemoveTag("RG");
        }

        if (opt_replace) {

            // only modify reads with an RG tag matching replace_string
            if (al.GetTag("RG", RG_tag) && RG_tag == replace_string) {
                if (! al.EditTag("RG", "Z", new_rg.ID)) {
                    cerr << NAME << " could not edit tag for read " << al.Name << endl;
                    return 1;
                }
            }

        } else if (! new_rg.ID.empty()) {

            al.AddTag("RG", "Z", new_rg.ID);

        }

        _DEBUG(0) if (n_reads <= debug_reads_to_report) {
            cerr << NAME << " " << n_reads << " read after processing: ";
            printAlignmentInfo(cerr, al);
        }

        writer.SaveAlignment(al);

        if (opt_progress && n_reads % opt_progress == 0)
            cerr << NAME << " " << n_reads << " reads processed..." << endl;
	}

    _DEBUG(0) cerr << NAME << " " << n_reads << " reads processed" << endl;

	reader.Close();
	writer.Close();

	return 0;
}


//-------------------------------------


const SamReadGroupDictionary
yoruba::parseReadGroupDictionaryString(const string& in)
{
    SamReadGroupDictionary rgd, empty_rgd;
    SamReadGroup rg;
    uint32_t pos, prev_pos;
    string this_tag, this_val;
    if (in.empty() || in.substr(0, 5) != "@RG\\t") {
        cerr << NAME << " dictionary string does not begin with '@RG\\t'" << endl;
        return empty_rgd;
    }
    pos = prev_pos = 5;
    this_tag = this_val = "";
    while (pos < in.length()) {
        switch (in[pos]) {
            case ':': 
                if (prev_pos < pos) {
                    this_tag = in.substr(prev_pos, (pos - prev_pos));
                    _DEBUG(0) cerr << "after ':', this_tag = " << this_tag << endl;
                } else return empty_rgd;
                prev_pos = pos + 1;
                break;
            case '\\': 
                if (pos + 1 >= in.length() 
                    || (in[pos + 1] != 't' && in[pos + 1] != 'n')) {
                    cerr << NAME << " only '\\t' and '\\n' escape chars allowed in dictionary string" << endl;
                    return empty_rgd;
                }

                if (prev_pos < pos) {
                    this_val = in.substr(prev_pos, (pos - prev_pos));
                    _DEBUG(0) cerr << "after '\\', this_val = " << this_val << endl;
                    if (this_tag == "ID") rg.ID = this_val;
                    else if (this_tag == "LB") rg.Library = this_val;
                    else if (this_tag == "SM") rg.Sample = this_val;
                    else if (this_tag == "DS") rg.Description = this_val;
                    else if (this_tag == "DT") rg.ProductionDate = this_val;
                    else if (this_tag == "PG") rg.Program = this_val;
                    else if (this_tag == "PL") rg.SequencingTechnology = this_val;
                    else if (this_tag == "PU") rg.PlatformUnit = this_val;
                    else if (this_tag == "PI") rg.PredictedInsertSize = this_val;
                    else if (this_tag == "FO") rg.FlowOrder = this_val;
                    else if (this_tag == "KS") rg.KeySequence = this_val;
                    else if (this_tag == "CN") rg.SequencingCenter = this_val;
                    else return empty_rgd;
                    if (in[pos + 1] == 'n') {
                        if (pos + 2 < in.length() && in.substr((pos + 2), 5) != "@RG\\t") {
                            cerr << NAME << " dictionary string does not continue with '@RG\\t'" << endl;
                            return empty_rgd;
                        }
                        rgd.Add(rg);
                        rg.Clear();
                        pos = pos + 5;
                    }
                }
                this_tag = this_val = "";
                pos = pos + 1;
                prev_pos = pos + 1;
                break;
            default:
                break;
        }
        ++pos;
    }
    if (! this_tag.empty() || ! this_val.empty()) {
        cerr << NAME << " malformed read group string '" << in << "'" << endl;
        return empty_rgd;
    }
    return rgd;
}

