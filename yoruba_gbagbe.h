// yoruba_gbagbe.h  (c) Douglas G. Scofield, Dept. Plant Physiology, Umeå University
//
// Header file for yoruba_gbagbe.cpp
//
// Gbagbe is the Yoruba (Nigeria) verb for 'to forget'.
//
// Uses BamTools C++ API for reading BAM files

#ifndef _YORUBA_GBAGBE_H_
#define _YORUBA_GBAGBE_H_


// Std C/C++ includes
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <string>
#include <list>
#include <fstream>
#include <sstream>
#include <map>
#include <tr1/unordered_map>

// BamTools includes: my own fork of https://github.com/pezmaster31/bamtools
#include "api/BamAux.h"
#include "api/BamReader.h"
#include "api/BamAlignment.h"
#include "api/SamHeader.h"
#include "api/SamReadGroup.h"
#include "api/SamReadGroupDictionary.h"
#include "api/SamProgram.h"
#include "api/SamProgramChain.h"
#include "api/BamWriter.h"

// SimpleOpt includes: http://code.jellycan.com/simpleopt, http://code.google.com/p/simpleopt/
#include "SimpleOpt.h"

// Yoruba includes
#include "yoruba.h"
#include "yoruba_util.h"

#ifndef _YORUBA_MAIN
#define NAME "[yoruba_forget]"
#endif

// Functions defined in yoruba_gbagbe.cpp
//
namespace yoruba {

int  main_gbagbe(int argc, char* argv[]);

}  // namespace yoruba

#endif // _YORUBA_GBAGBE_H_
