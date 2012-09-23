// yoruba_inu.h  (c) Douglas G. Scofield, Dept. Plant Physiology, Umeå University
//
// Header file for yoruba_inu.cpp
//
// Inu is the Yoruba (Nigeria) noun for 'inside'.
//
// Uses BamTools C++ API for reading BAM files

#ifndef _YORUBA_INU_H_
#define _YORUBA_INU_H_


// Std C/C++ includes
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <string>
#include <list>

// BamTools includes: https://github.com/pezmaster31/bamtools
#include "api/BamReader.h"
#include "api/BamAlignment.h"
#include "api/SamHeader.h"
#include "api/SamReadGroup.h"
#include "api/SamReadGroupDictionary.h"
#include "api/SamProgram.h"
#include "api/SamProgramChain.h"

// SimpleOpt includes: http://code.jellycan.com/simpleopt, http://code.google.com/p/simpleopt/
#include "SimpleOpt.h"

// Yoruba includes
#include "yoruba.h"
#include "yoruba_util.h"

#ifndef _YORUBA_MAIN
#define NAME "[yoruba_inside]"
#endif

// Functions defined in yoruba_inu.cpp
//
namespace yoruba {

int  main_inu(int argc, char* argv[]);

}  // namespace yoruba

#endif // _YORUBA_INU_H_
