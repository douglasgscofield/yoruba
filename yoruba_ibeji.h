// yoruba_kojopodipo.h  (c) Douglas G. Scofield, douglasgscofield@gmail.com

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
#include "api/SamProgram.h"
#include "api/SamProgramChain.h"

// SimpleOpt includes: http://code.jellycan.com/simpleopt, http://code.google.com/p/simpleopt/
#include "SimpleOpt.h"

// Yoruba includes
#include "yoruba.h"
#include "yoruba_util.h"

#ifndef _YORUBA_MAIN
#define NAME "[yoruba_readgroup]"
#endif

// Functions defined in yoruba_kojopodipo.cpp
//
namespace yoruba {

int  main_kojopodipo(int argc, char* argv[]);

}  // namespace yoruba

#endif // _YORUBA_KOJOPODIPO_H_
