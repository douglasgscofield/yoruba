#ifndef _YORUBA_UTILS_H_
#define _YORUBA_UTILS_H_


// Std C/C++ includes
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <string>

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
printAlignmentInfo(std::ostream& os, const BamTools::BamAlignment& alignment, int32_t level = 0);

void
printAlignmentInfo(std::ostream& os, const BamTools::BamAlignment& alignment, const BamTools::RefVector& refs, int32_t level = 0);

void
printAlignmentInfo_fields(std::ostream& os, const BamTools::BamAlignment& alignment, int32_t level = 0);

void
printAlignmentInfo_fields(std::ostream& os, const BamTools::BamAlignment& alignment, const BamTools::RefVector& refs, int32_t level = 0);


}  // namespace yoruba


#endif /* _YORUBA_UTILS_H_ */

