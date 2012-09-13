// gbagbe.cpp  (c) Douglas G. Scofield, Dept. Plant Physiology, Umeå University
//
// Gbagbe dynamically reassociating reads with reference sequences in a BAM
// file.  In automatic mode, gbagbe removes reference sequences from the BAM
// header that are not actually referred to by reads.  In directed mode, gbagbe
// will replace the BAM header with one that is provided by the user, and
// reassociate reads with reference sequences in the new header.
//
// Gbagbe is the Yoruba (Nigeria) word for 'to forget'.
//
// ## When to use gbagbe
//
// Gbagbe can be particularly helpful when the set of reference sequences to
// which reads are mapped is quite large, as for a genome assembly containing
// many contigs, and a BAM file is created from a mapping that represents a set
// of reads mapped to a portion of the reference sequences, for example the set
// of reads overlapping a region of a contig in which a particular locus is
// found.
//
// This is because the BAM file resulting from a 'samtools view' command with a
// region specified contains the full header of the original BAM file, although
// the extracted region may refer to only one or a few of the reference
// sequences of the original BAM.  When the original BAM header contains a few
// reference sequences, this is not a problem.  When the BAM header contains
// references to millions of contigs, the header can take up much more space in
// the final BAM file than the reads mapped to the region of interest.  This
// wastes disk space and increases the time required for downstream tools to
// read the BAM file.
//
// ## Alternatives to gbagbe
//
// 1. This same task can be accomplished with samtools by converting the BAM
// file to SAM format without the header (leaving off -h), and then
// reconverting back to BAM:
//
//    samtools view this.bam | samtools view -Sb > this_reheader.bam
//
// 2. The 'samtools reheader' very quickly replaces the header of a BAM file
// with a new header, but does not adjust references to @SQ lines in reads.
// This is useful for adding/subtracting other '@' lines, but much caution
// should be used if @SQ lines are different in the new header, specifically
// all reads must refer to reference sequences that are adjustment to the
// reference sequence references in reads.
//
// Gbagbe has two modes of operation, automatic and directed.
//
//    1) Automatic mode, in which the header of the BAM file is reduced to
//    contain only those reference sequences to which reads actually refer.
//    This requires two passes through the BAM file, the first to determine
//    which reference sequences to keep, and the second to reheader the file.
//
//    2) Directed mode (options -f, -h), in which a header is provided via
//    command-line option.  The BAM file is read and reads are adjusted to
//    refer to reference sequences in the new header.  If a read refers to a
//    reference sequence whose name does not appear in the new header, or has a
//    position that is inconsistent with the new header (mapped beyond the end
//    of the stated length if its reference), gbagbe exits with an error.
//    Gbagbe is only able to read from stdin in directed mode.
//
//       gbagbe -f newheader.sam in.bam > out.bam
//
//       gbagbe -h header-source.bam in.bam > out.bam
//
//       samtools 
//
// Gbagbe uses BamTools C++ API for reading BAM files


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

// link pair candidate
//
// A single read mapped such that it is mapped near the end of the contig and
// is oriented so that its mate may lie off the contig, taking paired-end
// library insert size into account.

//            rd1>
//    |===========>    or     |==========>   or etc.
//                              <rd2

// link pair
//
// A pair of reads in which each read is a link pair candidate and each is mapped
// to a different contig.

//            rd1>
//    |===========>  and  |==========>   or other such compatible configurations
//                          <rd2

// link pair orphan
//
// A pair of reads in which one read is a link pair candidate and the other is
// unmapped

//            rd1>                   with rd1 unmapped
//    |===========>          or      |==========>        or etc.
//    with rd2 unmapped                <rd2

// broken link pair
//
// A pair of reads in which one read is a link pair candidate and the other is
// mapped but is not a link pair candidate

//            rd1>                   rd1> 
//    |===========>   or      |==========>    |==========>  or etc.
//      rd2>                                          <rd2

//
//
// One read is mapped near the end of a contig, its mate is not mapped to the
// contig oriented in such a way that its mate is likely to be off the contig: 
//
// Find reads in BAM files that are mapped near the end of a contig and
// oriented so they point off the contig

// Std C/C++ includes
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <string>
#include <list>
using namespace std;

// BamTools includes
#include "api/BamMultiReader.h"
#include "api/BamWriter.h"
#include "api/BamAlignment.h"
using namespace BamTools;

// SimpleOpt includes: http://code.jellycan.com/simpleopt, http://code.google.com/p/simpleopt/
#include "SimpleOpt.h"

// Ibeji includes
#include "ibejiAlignment.h"
#include "processReadPair.h"
#include "utils.h"
using namespace ibeji;

string  output_bam_filename = "test.bam";
int32_t pairs_to_process = 20;
int32_t max_read_length = 101;
int32_t link_pair_total_tail = 1000;
int32_t link_pair_crit_tail = 999999;
bool    link_pair_diff_chrom = true;

int32_t mate_tail_est_crit = link_pair_total_tail + max_read_length;

bool    debug_ref_mate = false;

int 
main(int argc, char* argv[]) {

	// validate argument count
	if( argc != 2 ) {
		cerr << "USAGE: " << argv[0] << " <input BAM file> " << endl;
		exit(1);
	}

	string filename = argv[1];
	//cerr << "Printing alignments from file: " << filename << endl;
	
	BamReader reader;
	if (!reader.Open(filename)) {
        cerr << "could not open filename " << filename << ", exiting" << endl;
        return 1;
    }
    cerr << filename << ": Done opening" << endl;

    // Header can't be used to accurately determine sort order because samtools never
    // changes it; instead, check after loading each read as is done with "samtools index"

    // We don't need to load an index (right?)
	// if (!reader.LocateIndex()) {
    //     const string index_filename = filename + ".bai";
	//     if (!reader.OpenIndex(index_filename)) {
    //         cerr << "could not open index" << endl;
    //     }
    // }

    // gbagbe workflow
    //
    // create newHeaderSQ, newHeaderOther data structures
    // *** directed mode
    // if (-h <file.bam>)
    //     extract header from <file.bam>, save in newHeader
    //     splitHeader(newHeader, "@SQ", newHeaderSQ, newHeaderOther)
    // else if (-f <file.sam>)
    //     translate <file.sam>, save in newHeader (is there an API available?)
    //     die if any non-header lines
    //     splitHeader(newHeader, "@SQ", newHeaderSQ, newHeaderOther)
    // else
    // *** automatic mode
    //     read header from in.bam, save in oldHeader
    //     splitHeader(oldHeader, "@SQ", oldHeaderSQ, oldHeaderOther)
    //     int32_t n_oldRefID;
    //     for (n_oldRefID = 1; ; ++n_oldRefID) {
    //         refmap[n_oldRefID] = struct { std::string ref_name; uint64 ref_length; uint64 ref_nreads; }
    //     }
    //     while alignments in in.bam
    //         ++refmap[al.RefID]->ref_nreads;
    //     end while
    //     int32_t newRefID;
    //     foreach key in refmap[]  // preserve original order
    //         if refmap[key]->ref_nreads > 0
    //             pushd newHeader, { std::string ref_name; uint64 ref_length; }
    //         fi
    //     end foreach
    // fi

            


    const SamHeader header = reader.GetHeader();
    cerr << filename << ": Done getting header" << endl;
    const RefVector refs = reader.GetReferenceData();
    cerr << filename << ": Done getting reference data" << endl;
	
    BamWriter writer;
    if (! output_bam_filename.empty()) {
        if (! writer.Open(output_bam_filename, header, refs)) {
            cerr << "Could not open BAM output file " << output_bam_filename << endl;
            exit(1);
        }
        cerr << filename << ": Done opening BAM output file " << output_bam_filename << endl;
    }

    alignmentMap read1Map;  // a single map, for all reads awaiting their mate
    typedef map<string,int32_t> stringMap;
    typedef stringMap::iterator stringMapI;
    stringMap ref_mates;
    // alignmentMap read1Map, read2Map;

	BamAlignment full_al;
    int32_t count = 0;
    uint32_t max_reads_in_map = 0;
    int32_t n_reads_skipped_unmapped = 0;
    int32_t n_reads_skipped_mate_unmapped = 0;
    int32_t n_reads_skipped_wont_see_mate = 0;
    int32_t n_reads_skipped_mate_tail_est = 0;
    int32_t n_reads_skipped_ref_mate = 0;
    int32_t n_reads = 0;
    int32_t n_singleton_reads = 0;
    int32_t last_RefID = -1;
    int32_t last_Position = -1;

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
            int32_t mate_tail_est = readTailS(al.IsMateMapped(), al.IsMateReverseStrand(),
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
	return 0;
}

