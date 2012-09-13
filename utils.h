/*
 * utils.h
 *
 *  Created on: Apr 9, 2012
 *      Author: douglasgscofield
 */

#ifndef UTILS_H_
#define UTILS_H_


// Std C/C++ includes
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <string>

#include "ibejiAlignment.h"

namespace yoruba {

void 
PrintAlignment(const BamTools::BamAlignment&);

void
printAlignmentInfo(const BamTools::BamAlignment& alignment, const BamTools::RefVector& refs, int32_t level = 0);

void
printAlignmentInfo_fields(const BamTools::BamAlignment& alignment, const BamTools::RefVector& refs, int32_t level = 0);


}  // namespace yoruba


#endif /* UTILS_H_ */

