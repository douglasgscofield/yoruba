// kojopodipo.cpp  (c) Douglas G. Scofield, Dept. Plant Physiology, Umeå University
//
// Kojopodipo (English command is readgroup) adds or modifies the read group in a BAM file.
//
// Kojopodipo rewrites the header and then annotates each read.  I've started this because
// picard's AddOrChangeReadGroups.jar uses a lot of memory and seems to take a long time for
// what should be a simple task.  We'll see if those opinions hold up :-)
//
// Kojopodipo is the Yoruba (Nigeria) word for group.
//
// Uses BamTools C++ API for reading BAM files


// CHANGELOG
//
//
//
// TODO
//
// Process options
// Open BAM file
// Read BAM header
// Attach read group information to header
// Create new BAM file
// Write newe header to new BAM file
// While alignments left in BAM file
//   Read alignment
//   Attach read group information to alignment
//   Write alignment to new BAM file
// Close BAM file
//
// Command line options
//
//   --id                / --RGID  <string>
//   --sequencing-center / --RGCN  <string>
//   --description       / --RGCN  <string>
//   --date              / --RGDT  <string>
//   --flow-order        / --RGFO  <string>
//   --key-sequence      / --RGKS  <string>
//   --library           / --RGLB  <string>
//   --programs          / --RGPG  <string>
//   --predicted-insert  / --RGPI  <string>
//   --platform          / --RGPL  <string>
//   --platform-unit     / --RGPU  <string>
//   --sample-name       / --RGSM  <string>
//
// No formatting restrictions are imposed on any of these elements, and none
// are added by default by kojopodipo.
//



// Std C/C++ includes
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <string>
#include <list>
using namespace std;

// BamTools includes
#include "api/BamReader.h"
#include "api/BamWriter.h"
#include "api/BamAlignment.h"
#include "api/SamHeader.h"
#include "api/SamReadGroup.h"
#include "api/SamReadGroupDictionary.h"
using namespace BamTools;

// SimpleOpt includes: http://code.jellycan.com/simpleopt, http://code.google.com/p/simpleopt/
#include "SimpleOpt.h"

// Ibeji includes
#include "utils.h"
using namespace yoruba;

// Finally, declarations of functions defined in this file, eventually to be moved to header file
//
void printReadGroupDictionary(std::ostream& os, const SamReadGroupDictionary& rgd);
void printReadGroup(std::ostream& os, const SamReadGroup& rg);

const int32_t reads_to_process = 2000000000;

const bool    debug_ref_mate = false;

int 
main(int argc, char* argv[]) {

	// validate argument count
	if( argc != 2 ) {
		cerr << "USAGE: " << argv[0] << " <input BAM file> " << endl;
		exit(1);
	}

    // Initial set of read group information
    SamReadGroup new_readgroup;
    new_readgroup.ID = "rg_ID";
    new_readgroup.Library = "rg_LB";
    new_readgroup.Sample = "rg_SM";
    new_readgroup.SequencingTechnology = "rg_PL";

    const string  output_bam_filename = "testout.bam";

	string filename = argv[1];
	cerr << "Modifying read groups in file: " << filename << endl;
	
	BamReader reader;
	if (!reader.Open(filename)) {
        cerr << "could not open filename " << filename << ", exiting" << endl;
        return 1;
    }
    cerr << filename << ": Done opening" << endl;

    SamHeader header = reader.GetHeader();
    cerr << filename << ": Done getting header" << endl;
    // No reason to explicitly load reference sequences, we never use them
    // const RefVector refs = reader.GetReferenceData();
    // cerr << filename << ": Done getting reference data" << endl;

    // now for the all-important read groups

    if (header.HasReadGroups()) {
        cerr << filename << ": Already has read groups in its dictionary" << endl;
        printReadGroupDictionary(cerr, header.ReadGroups);
        cerr << filename << ": Clearing read group dictionary" << endl;
        header.ReadGroups.Clear();
    } else {
        cerr << filename << ": Read group dictionary is empty" << endl;
    }
    header.ReadGroups.Add(new_readgroup);
    cerr << filename << ": Added read group" << endl;
    printReadGroupDictionary(cerr, header.ReadGroups);
    cerr << endl;
	
    BamWriter writer;
    if (! output_bam_filename.empty()) {
        if (! writer.Open(output_bam_filename, header, reader.GetReferenceData())) {
            cerr << "Could not open BAM output file " << output_bam_filename << endl;
            exit(1);
        }
        cerr << filename << ": Done opening BAM output file " << output_bam_filename << endl;
    }

	BamAlignment al;
    int32_t count = 0;
    int32_t n_reads = 0;

    cerr << filename << ": Modifying up to " << reads_to_process << " reads" << endl;

	while (reader.GetNextAlignment(al) 
           && (! reads_to_process || count < reads_to_process)) {

        ++count;
        ++n_reads;

        // printAlignmentInfo(al, refs);

        string RG_tag;
        if (al.GetTag("RG", RG_tag)) {
            // cout << "alignment " << al.Name << " already has read group tag: RG:Z:'" 
            //    << RG_tag << "'" << endl;
            if (! al.EditTag("RG", "Z", new_readgroup.ID)) {
                cerr << "could not edit tag for read " << al.Name << endl;
                return 1;
            }
        } else {
            if (! al.AddTag("RG", "Z", new_readgroup.ID)) {
                cerr << "could not add tag to read " << al.Name << endl;
                return 1;
            }
        }

        // printAlignmentInfo(al, refs);
        // cout << endl;

        writer.SaveAlignment(al);

	}

	cerr << "===============================" << endl;
    cerr << count << " reads processed" << endl;
    cerr << n_reads << " total reads" << endl;

	reader.Close();
    if (! output_bam_filename.empty()) {
	    writer.Close();
    }
	return 0;
}


//-------------------------------------


void
printReadGroupDictionary(std::ostream& os, const SamReadGroupDictionary& rgd)
{
    os << "SamReadGroupDictionary" << std::endl;
    if (rgd.IsEmpty()) {
        os << "** empty **" << std::endl;
        return;
    }
    int i = 0;
    for (SamReadGroupConstIterator rgdI = rgd.ConstBegin(); 
            rgdI != rgd.ConstEnd(); ++rgdI, ++i) {
        os << i << "-th entry" << std::endl;
        printReadGroup(os, *rgdI);
    }
}


//-------------------------------------


void
printReadGroup(std::ostream& os, const SamReadGroup& rg)
{
    os << "SamReadGroup" << std::endl;
    if (rg.HasID())
        os << "@RG ID:'" << rg.ID << "'" << std::endl;
    if (rg.HasSequencingCenter()) 
        os << "@RG CN:'" << rg.SequencingCenter << "'" << std::endl;
    if (rg.HasDescription())
        os << "@RG DS:'" << rg.Description << "'" << std::endl;
    if (rg.HasProductionDate())   
        os << "@RG DT:'" << rg.ProductionDate << "'" << std::endl;
    if (rg.HasFlowOrder()) 
        os << "@RG FO:'" << rg.FlowOrder << "'" << std::endl;
    if (rg.HasKeySequence()) 
        os << "@RG KS:'" << rg.KeySequence << "'" << std::endl;
    if (rg.HasLibrary()) 
        os << "@RG LB:'" << rg.Library << "'" << std::endl;
    if (rg.HasProgram()) 
        os << "@RG PG:'" << rg.Program << "'" << std::endl;
    if (rg.HasPredictedInsertSize()) 
        os << "@RG PI:'" << rg.PredictedInsertSize << "'" << std::endl;
    if (rg.HasSequencingTechnology()) 
        os << "@RG PL:'" << rg.SequencingTechnology << "'" << std::endl;
    if (rg.HasPlatformUnit()) 
        os << "@RG PU:'" << rg.PlatformUnit << "'" << std::endl;
    if (rg.HasSample()) 
        os << "@RG SM:'" << rg.Sample << "'" << std::endl;
}

