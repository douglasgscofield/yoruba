// yoruba_sefibo.h  (c) Douglas G. Scofield, Dept. Plant Physiology, Umeå University

#ifndef _YORUBA_SEFIBO_H_
#define _YORUBA_SEFIBO_H_

// Std C/C++ includes
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <string>
#include <list>

// BamTools includes: https://github.com/pezmaster31/bamtools
#include "api/BamMultiReader.h"
#include "api/BamReader.h"
#include "api/BamWriter.h"
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
// #include "ibejiAlignment.h"  lightweight alignment class not completed
#include "processReadPair.h"  // needed for now, probably not in future

#ifndef _YORUBA_MAIN
#define NAME "[yoruba_insertsize]"
#endif

// Functions defined in yoruba_sefibo.cpp
//
namespace yoruba {

int  main_sefibo(int argc, char* argv[]);

}  // namespace yoruba

#endif // _YORUBA_SEFIBO_H_
