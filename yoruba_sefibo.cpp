// yoruba_sefibo.cpp  (c) Douglas G. Scofield, Dept. Plant Physiology, Umeå University
//
// Sefibo calculates insert sizes based on read mappings within BAM files.
//
// Sefibo only works with reads mapping to the same chromosome/contig.  Beyond that 
// restriction, many options for insert size calculation are provided.
//
// Sefibo is the Yoruba (Nigeria) noun for 'insert' (insertion point, etc.).
//
// Uses BamTools C++ API for reading BAM files


// CHANGELOG
//
//
//
// TODO
//
// process name-sorted BAM files with --name/-n
// regular expression handling for matching strings
// command line opotion processing via SimpleOpt.h
// debugging options
// lightweight alignment class to reduce memory usage?
//
// Command line options
//
//   --read-orientation/-r FR|FF|RR|RF
//   
//   --insert-type/-t outer(5'L to 5'R)|inner (3'L to 3'R)|left(5'L to 3'R|right(3'R to 5'L)
//
//   --quantiles/-q <comma-separated list of quantiles to produce, 0-1>
//
//   --bam-insert-size/-b  calculate statistics using insert size recorded in the BAM file
//
//   --better-estimate  Adjust insert size estimate using complement of Kristoffer's gap size 
//                      estimator.  In short, the set of reads that both map to the same contig
//                      represents a biased sample of reads in terms of insert size.  Read
//                      pairs that fail to be 'properly paired' also provide insert size
//                      information, albeit partial, as a result of the tails that hang off of
//                      their 3' ends.
//
// Consider implementing options below
//
//   --chrom <string>   Restrict examination to pairs involving chromosome <string>, may be 
//   -c <string>        specified multiple times for multiple chromosomes.  <string> may be
//                      a regular expression
//
//   --chrom-file <string> file containing chromosome names/regular expressions, as above
//   -C <file>
//
//   --chrom-exclude    Only consider pairs for which at least one read maps to a chromosomes
//   -x                 that DOES NOT MATCH the --chrom/--chrom-file name(s)
//
//   --chrom-and        Only consider pairs for which both reads map to chromosomes that
//   -a                 match (don't match, -x) the --chrom/--chrom-file name(s)



// First, some definitions

// tail
//
// The distance between the start of a mapped read and the appropriate end of
// the reference sequence.  For a forward-oriented read, this is the distance
// between the leftmost mapped position and the 3' end of the reference, and
// has positive sign.  For a reverse-oriented read, this is the distance
// between the rightmost mapped position and the 5' end of the reference, and
// has negative sign.  The absolute value of the tail is the total number of
// bases.

// total tail
//
// The sum of the absolute values of the tails for two link pair candidate
// reads.

#include "yoruba_sefibo.h"

using namespace std;
using namespace BamTools;
using namespace yoruba;

//options
static string  output_bam_filename = "test.bam";
static int64_t pairs_to_process = 20;
static int64_t max_read_length = 101;
static int64_t link_pair_total_tail = 1000;
static int64_t link_pair_crit_tail = 999999;
static bool    link_pair_diff_chrom = true;
static int64_t mate_tail_est_crit = link_pair_total_tail + max_read_length;
static bool    debug_ref_mate = false;


//-------------------------------------


#ifdef _STANDALONE
int 
main(int argc, char* argv[]) {
    return main_sefibo(argc, argv);
}
#endif


//-------------------------------------


static int
usage(bool long_help = false)
{
    cerr << endl;
    cerr << "Usage:   " << YORUBA_NAME << " insertsize [options] <in.bam>" << endl;
    cerr << "         " << YORUBA_NAME << " sefibo [options] <in.bam>" << endl;
    cerr << endl;
    cerr << "Either command invokes this function." << endl;
    cerr << endl;
    cerr << "\
Calculate the insert size distribution among alignments in <in.bam>.\n\
\n\
Options: --read-orientation | -r  FR|FF|RR|RF          expected orientation of reads\n\
         --insert-type | -t  outer|inner|left|right    insert type to calculate\n\
         --quantiles | -q LIST                         list of quantiles to report for distribution\n\
         --better-estimate | -b                        improve insertion distribution calculation\n\
\n";
    if (long_help) {
        cerr << "\
By default, all reads in the BAM file will be given the supplied read group.\n\
If the dictionary already defines a read group with the same ID, its definition\n\
will be replaced with the supplied information.  If the dictionary contains\n\
other read groups, their definitions will remain in the BAM file header but\n\
all reads will be given the supplied read group.\n\
\n\
The --insert-type argument specifies the manner in which the insert size should\n\
be calculated.  Assuming orientations relative to the forward strand, insert size\n\
is calculates as follows: \n\
\n\
    outer   [default] maximal 5' extent of the 5'-most read to maximal 3' extent \n\
            of the 3'-most read, inclusive;\n\
    inner   maximal 3' extent of the 5'-most read to maximal 5' extent \n\
            of the 3'-most read, exclusive;\n\
    left    maximal 5' extent of the 5'-most read to maximal 5' extent \n\
            of the 3'-most read, inclusive of the left, exclusive of the right;\n\
    right   maximal 3' extent of the 5'-most read to maximal 3' extent \n\
            of the 3'-most read, exclusive of the left, inclusive of the right;\n\
\n\
\n\
\n\
\n\
\n";
    }
    cerr << "         --? | -? | --help                   longer help" << endl;
    cerr << endl;
#ifdef _WITH_DEBUG
    cerr << "         --debug INT     debug info level INT [" << opt_debug << "]" << endl;
    cerr << "         --reads INT     process at most this many reads [" << opt_reads << "]" << endl;
    cerr << "         --progress INT  print reads processed mod INT [" << opt_progress << "]" << endl;
    cerr << endl;
#endif
    cerr << "Sefibo is the Yoruba (Nigeria) noun for 'insert'." << endl;
    cerr << endl;
    }

    return EXIT_FAILURE;
}


//-------------------------------------


int 
yoruba::main_sefibo(int argc, char* argv[]) {

	// first process options
	if( argc < 2 ) {
		return usage();
	}

	string filename = argv[1];
	//cerr << "Printing alignments from file: " << filename << endl;
	
	BamReader reader;
	if (!reader.Open(filename)) {
        cerr << "could not open filename " << filename << ", exiting" << endl;
        return EXIT_FAILURE;
    }

    // Header can't be used to accurately determine sort order because samtools never
    // changes it; instead, check after loading each read as is done with "samtools index"

    const SamHeader header = reader.GetHeader();

    const RefVector refs = reader.GetReferenceData();

	
    BamWriter writer;
    if (! output_bam_filename.empty()) {
        if (! writer.Open(output_bam_filename, header, refs)) {
            cerr << "Could not open BAM output file " << output_bam_filename << endl;
            exit(1);
        }
        cerr << filename << ": Done opening BAM output file " << output_bam_filename << endl;
    }

    alignmentMap read1Map;  // a single map, for all reads awaiting their mate
    typedef map<string,int64_t> stringMap;
    typedef stringMap::iterator stringMapI;
    stringMap ref_mates;
    // alignmentMap read1Map, read2Map;

	BamAlignment full_al;
    int64_t count = 0;
    int64_t max_reads_in_map = 0;
    int64_t n_reads_skipped_unmapped = 0;
    int64_t n_reads_skipped_mate_unmapped = 0;
    int64_t n_reads_skipped_wont_see_mate = 0;
    int64_t n_reads_skipped_mate_tail_est = 0;
    int64_t n_reads_skipped_ref_mate = 0;
    int64_t n_reads = 0;
    int64_t n_singleton_reads = 0;
    int64_t last_RefID = -1;
    int64_t last_Position = -1;

    cerr << filename << ": Looking for up to " << pairs_to_process << " link pairs,"
        << " total tail = " << link_pair_total_tail 
        << " critical tail = " << link_pair_crit_tail 
        << ", must be on diff chromosome = " << link_pair_diff_chrom << endl;

	while (reader.GetNextAlignment(full_al) 
           && (! pairs_to_process || count < pairs_to_process)) {

        BamAlignment al = full_al;

        //printAlignmentInfo(al, refs);
        //++count;
        ++n_reads;

        if (last_RefID < 0) last_RefID = al.RefID;
        if (last_Position < 0) last_Position = al.Position;
        if (al.RefID > last_RefID) {
            // We've moved to the next reference sequence
            // Clean up reads with mates expected here that haven't been seen
            if (debug_ref_mate) {
                cerr << "MISSED " << ref_mates.size() << " ref_mates on this reference "
                    << last_RefID << " " << refs[last_RefID].RefName << endl;
            }
            for (stringMapI rmI = ref_mates.begin(); rmI != ref_mates.end(); ++rmI) {
                ++n_reads_skipped_ref_mate;
                read1Map.erase(read1Map.find(rmI->first));
                ref_mates.erase(ref_mates.find(rmI->first));
            }
            last_RefID = al.RefID;
            last_Position = al.Position;
        } else if (al.RefID < last_RefID) {
            cerr << filename << " does not appear to be sorted, chromosome out of order: "
                 << last_RefID << " (" << refs[last_RefID].RefName << ") "
                 << al.RefID << " (" << refs[al.RefID].RefName << ") " << endl;
            exit(1);
        } else if (al.Position < last_Position) {
            cerr << filename << " does not appear to be sorted, reads out of order: "
                 << last_Position << " (last read) "
                 << al.Position << " (" << al.Name << ")" << endl;
            exit(1);
        }

        if (! al.IsMapped()) { ++n_reads_skipped_unmapped; continue; }

        if (! al.IsMateMapped()) { ++n_reads_skipped_mate_unmapped; continue; }

        alignmentMapI mI = read1Map.find(al.Name);

        if (mI == read1Map.end()) {
            // the read name has not been seen before

            if (al.MateRefID < al.RefID
                || (al.MateRefID == al.RefID && al.MatePosition < al.Position)) {
                // we should have seen its mate earlier, so skip it
                ++n_reads_skipped_wont_see_mate;
                continue;
            }

            // If the mate likely to also be a link pair candidate, add the read
            int64_t mate_tail_est = readTailS(al.IsMateMapped(), al.IsMateReverseStrand(),
                            al.MatePosition, refs[al.MateRefID].RefLength, max_read_length);
            if (mate_tail_est <= mate_tail_est_crit) {
                // the mate tail estimate suggests it might be a link pair candidate
                read1Map[al.Name] = al;  // add the read to the map
            } else {
                // the mate tail estimate appears too long for the mate to be a candidate
                ++n_reads_skipped_mate_tail_est;
                continue;
            }

            if (read1Map.size() > max_reads_in_map) max_reads_in_map = read1Map.size();
            if (al.MateRefID == al.RefID && al.MatePosition >= al.Position) {
                // the mate is expected later on this contig
                ref_mates[al.Name] = al.MateRefID;
            }

        } else {
            // get the mate's alignment, and process the pair

            const BamAlignment& al_mate = mI->second;

            if (processReadPair(al, al_mate, refs, link_pair_total_tail, 
                                link_pair_crit_tail, link_pair_diff_chrom)) {
                ++count;

                // write to the new BAM file, if the string is not empty
                if (! output_bam_filename.empty()) {
                    writer.SaveAlignment(al_mate);  // the first one seen
                    writer.SaveAlignment(al);  // the second one seen
                }
            }

            read1Map.erase(mI);

            if (al.MateRefID == al.RefID) {
                stringMapI rmI = ref_mates.find(al.Name);
                if (rmI == ref_mates.end()) {
                    cerr << "expected a ref_mate, couldn't find its name: " << al.Name << endl;
                    exit(1);
                }
                ref_mates.erase(rmI);
            }

        }

	}

	cerr << "===============================" << endl;
    cerr << read1Map.size() << " alignments left in read1Map" << endl;
    cerr << max_reads_in_map << " maximum number of reads in read1Map" << endl;
    cerr << count << " pairs processed" << endl;
	cerr << "===============================" << endl;
    cerr << n_reads << " total reads" << endl;
    cerr << n_singleton_reads << " singleton reads" << endl;
    cerr << n_reads_skipped_unmapped << " reads skipped because unmapped" << endl;
    cerr << n_reads_skipped_mate_unmapped << " reads skipped because mate unmapped" << endl;
    cerr << n_reads_skipped_wont_see_mate << " reads skipped because mate won't be seen" << endl;
    cerr << n_reads_skipped_mate_tail_est << " reads skipped because mate tail appears too long" << endl;
    cerr << n_reads_skipped_ref_mate << " reads skipped because mate not on reference" << endl;

	reader.Close();
    if (! output_bam_filename.empty()) {
	    writer.Close();
    }
	return EXIT_SUCCESS;
}

