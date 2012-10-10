#ifndef _YORUBA_UTIL_H_
#define _YORUBA_UTIL_H_


// Std C/C++ includes
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <string>

// #define NDEBUG  // uncomment to remove assert() code
#include <assert.h>

// BamTools includes: https://github.com/pezmaster31/bamtools
#include "api/BamReader.h"
#include "api/BamAlignment.h"
#include "api/SamHeader.h"
#include "api/SamReadGroup.h"
#include "api/SamReadGroupDictionary.h"
#include "api/SamProgram.h"
#include "api/SamProgramChain.h"
#include "ibejiAlignment.h"

#ifdef _WITH_DEBUG
#define IF_DEBUG(__lvl__) if (opt_debug >= __lvl__)
#define _DEBUG(__lvl__) if (opt_debug >= __lvl__)
#define DEBUG(__lvl__) (opt_debug >= __lvl__)
#else
#define IF_DEBUG(__lvl__) if (0)
#define _DEBUG(__lvl__) if (0)
#define DEBUG(__lvl__) (false)
#endif

namespace yoruba {

bool 
isCoordinateSorted(int32_t ref, int32_t pos, int32_t prev_ref, int32_t prev_pos);

bool 
isMateUpstream(const BamTools::BamAlignment&);

bool 
isMateDownstream(const BamTools::BamAlignment&);

void 
PrintAlignment(const BamTools::BamAlignment&);

void 
PrintAlignment(std::ostream& os, 
               const BamTools::BamAlignment& alignment);

void
printAlignmentInfo(std::ostream& os, 
               const BamTools::BamAlignment& alignment, 
               int32_t level = 0);

void
printAlignmentInfo(std::ostream& os, 
               const BamTools::BamAlignment& alignment, 
               const BamTools::RefVector& refs, 
               int32_t level = 0);

void
printAlignmentInfo_fields(std::ostream& os, 
               const BamTools::BamAlignment& alignment, 
               int32_t level = 0);

void
printAlignmentInfo_fields(std::ostream& os, 
               const BamTools::BamAlignment& alignment, 
               const BamTools::RefVector& refs, 
               int32_t level = 0);

void
printReadGroupDictionary(std::ostream& os, 
               const BamTools::SamReadGroupDictionary& rgd,
               const std::string& prefix = "",
               const std::string& rg_prefix = "@RG",
               const std::string& rg_sep = "\t",
               const std::string& rg_delim = "'",
               const std::string& rg_terminate = "\n");

void
printReadGroup(std::ostream& os, 
               const BamTools::SamReadGroup& rg,
               const std::string& rg_prefix = "@RG",
               const std::string& rg_sep = "\t",
               const std::string& rg_delim = "'",
               const std::string& rg_terminate = "\n");


} // namespace yoruba

#endif /* _YORUBA_UTIL_H_ */

