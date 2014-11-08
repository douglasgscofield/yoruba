// yoruba_seda.h  Copyright (c) Douglas G. Scofield, douglasgscofield@gmail.com

#ifndef _YORUBA_SEDA_H_
#define _YORUBA_SEDA_H_

// Std C/C++ includes
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <list>
#include <map>
// #ifdef C++11
// some appropriate include
// #include <unordered_map>
// #include <unordered_set>
// #else
#include <tr1/unordered_map>
#include <tr1/unordered_set>
// #endif
#include <new>

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
// #include "yoruba_lightAlignment.h"  // do I need this for 'yoruba seda'?
#include "yoruba_util.h"

#ifndef _YORUBA_MAIN
#define NAME "[yoruba_duplicate]"
#endif

// Functions defined in yoruba_seda.cpp
//
namespace yoruba {

int  main_seda(int argc, char* argv[]);

}  // namespace yoruba

#endif // _YORUBA_SEDA_H_
