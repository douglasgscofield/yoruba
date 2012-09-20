#ifndef _YORUBA_UTIL_H_
#define _YORUBA_UTIL_H_


// Std C/C++ includes
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <string>

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
#define _DEBUG(__lvl__) if (opt_debug >= __lvl__)
#else
#define _DEBUG(__lvl__) if (0)
#endif

namespace yoruba {

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

