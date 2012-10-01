
#include <cstdlib>
#include <iostream>
#include <string>

#undef _STANDALONE
#undef _IMPLEMENTED
#define _YORUBA_MAIN

#include "yoruba.h"
#include "yoruba_gbagbe.h"
#include "yoruba_inu.h"
#include "yoruba_kojopodipo.h"
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
    cerr << "         Each command has two names, one in English and one in Yoruba" << endl << endl;
    cerr << "Command: forget     | gbagbe       remove unused reference sequences" << endl;
    cerr << "         inside     | inu          display summary of bam file contents" << endl;
    cerr << "         readgroup  | kojopodipo   add or modify read group information" << endl;
#ifdef _IMPLEMENTED
    cerr << "         insertsize | sefibo       calculates insert sizes" << endl;
    cerr << "         twinreads  | ibeji        find reads paired in various ways" << endl;
#endif
    cerr << endl;

    return 1;
}

int
main (int argc, char* argv[])
{
    if (argc < 2) return usage();
    string cmd = argv[1];
    if (cmd == "forget" || cmd == "gbagbe") 
        return main_gbagbe(argc-1, argv+1);
    else if (cmd == "inside" || cmd == "inu") 
        return main_inu(argc-1, argv+1);
    else if (cmd == "readgroup" || cmd == "kojopodipo") 
        return main_kojopodipo(argc-1, argv+1);
#ifdef _IMPLEMENTED
    else if (cmd == "insert" || cmd == "sefibo") 
        return main_sefibo(argc-1, argv+1);
    else if (cmd == "twinreads" || cmd == "ibeji") 
        return main_ibeji(argc-1, argv+1);
#endif
    else {
        cerr << "Unrecognized command '" << argv[1] << "'" << endl;
        return 1;
    }
    return 0;
}

