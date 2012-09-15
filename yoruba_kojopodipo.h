// kojopodipo.h  (c) Douglas G. Scofield, Dept. Plant Physiology, Umeå University
//
// Header file for Kojopodipo
//
// Kojopodipo (English command is readgroup) adds or modifies the read group in a BAM file.
//
// Kojopodipo is the Yoruba (Nigeria) word for group.
//
// Uses BamTools C++ API for reading BAM files

#ifndef _YORUBA_KOJOPODIPO_H_
#define _YORUBA_KOJOPODIPO_H_


// Std C/C++ includes
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <string>
#include <list>

// BamTools includes: https://github.com/pezmaster31/bamtools
#include "api/BamReader.h"
#include "api/BamWriter.h"
#include "api/BamAlignment.h"
#include "api/SamHeader.h"
#include "api/SamReadGroup.h"
#include "api/SamReadGroupDictionary.h"

// SimpleOpt includes: http://code.jellycan.com/simpleopt, http://code.google.com/p/simpleopt/
#include "SimpleOpt.h"

// Yoruba includes
#include "yoruba.h"
#include "yoruba_utils.h"

#define NAME "[yoruba_readgroup]:"

// Functions defined in yoruba_kojopodipo.cpp
//
namespace yoruba {

int  main_kojopodipo(int argc, char* argv[]);
void printReadGroupDictionary(std::ostream& os, const BamTools::SamReadGroupDictionary& rgd);
void printReadGroup(std::ostream& os, const BamTools::SamReadGroup& rg);

}  // namespace yoruba

#endif // _YORUBA_KOJOPODIPO_H_
