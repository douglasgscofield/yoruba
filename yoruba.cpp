// yoruba.cpp  (c) Douglas G. Scofield, Dept. Plant Physiology, Umeå University
//
// Yoruba queries and manipulates BAM files.  The goal is to provide either 
// additional functionality or more efficient and/or flexible versions of 
// operations provided by samtools and Picard tools.
//
// Yoruba is a major language spoken around Nigeria.
//
// Uses BamTools C++ API for handling BAM files


// CHANGELOG
//
//
//
// TODO
// --- add typedefs for reference id (int32_t), number of reads (int64_t), etc. 
//     consistent with real world and SAM specification... or does BamTools 
//     provide these?


#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <string>
#include <ctime>

#undef _STANDALONE
#undef _IMPLEMENTED
#define _YORUBA_MAIN

#include "yoruba.h"
#include "yoruba_gbagbe.h"
#include "yoruba_inu.h"
#include "yoruba_kojopodipo.h"
#include "yoruba_seda.h"
#include "yoruba_util.h"
#ifdef _IMPLEMENTED
#include "yoruba_sefibo.h"
#include "yoruba_ibeji.h"
#endif

using namespace std;
using namespace yoruba;

// Yep, this looks a bit like samtools' bamtk.c.  I like its brevity and the simplicity
// of the "tool <command> [options]" framework.

static int
usage()
{
    cerr << endl;
    cerr << "Program: " << YORUBA_NAME << " -- query and manipulate BAM files" << endl;
    cerr << "Version: " << YORUBA_VERSION << endl << endl;
    cerr << "Usage:   " << YORUBA_NAME << " <command> [options]" << endl << endl;
    cerr << "Each command has two names, one in English and one in Yoruba" << endl << endl;
    cerr << "Command: forget     | gbagbe       remove unused reference sequences" << endl;
    cerr << "         inside     | inu          display summary of BAM file contents" << endl;
    cerr << "         readgroup  | kojopodipo   add or modify read group information" << endl;
    cerr << "         duplicate  | seda         mark (and remove) duplicate reads" << endl;
#ifdef _IMPLEMENTED
    cerr << "         insertsize | sefibo       calculates insert sizes" << endl;
    cerr << "         twinreads  | ibeji        find reads paired in various ways" << endl;
#endif
    cerr << endl;

    return EXIT_FAILURE;
}

int
main (int argc, char* argv[])
{
    clock_t time_start = clock();
    if (argc < 2) return usage();
    string cmd = argv[1];
    int retval = EXIT_SUCCESS;
    if (cmd == "forget" || cmd == "gbagbe") 
        retval = main_gbagbe(argc-1, argv+1);
    else if (cmd == "inside" || cmd == "inu") 
        retval = main_inu(argc-1, argv+1);
    else if (cmd == "readgroup" || cmd == "kojopodipo") 
        retval = main_kojopodipo(argc-1, argv+1);
    else if (cmd == "duplicate" || cmd == "seda") 
        retval = main_seda(argc-1, argv+1);
#ifdef _IMPLEMENTED
    else if (cmd == "insert" || cmd == "sefibo") 
        retval = main_sefibo(argc-1, argv+1);
    else if (cmd == "twinreads" || cmd == "ibeji") 
        retval = main_ibeji(argc-1, argv+1);
#endif
    else {
        cerr << "Unrecognized command '" << argv[1] << "'" << endl;
        retval = EXIT_FAILURE;
    }
    float runtime = ((float)(clock() - time_start)) / CLOCKS_PER_SEC;
    cerr << NAME << " runtime " << fixed << setprecision(3) << runtime << " sec" << endl;
    return retval;
}

