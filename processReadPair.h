/*
 * processReadPair.h
 *
 * Header file with stuff for processing a read pair.
 *
 *  Created on: Apr 9, 2012
 *      Author: douglasgscofield
 */

#ifndef PROCESSREADPAIR_H_
#define PROCESSREADPAIR_H_


// Std C/C++ includes
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <string>

#include "ibejiAlignment.h"
#include "utils.h"


namespace ibeji {

bool 
processReadPair(const BamAlignment& al1, 
        const BamAlignment& al2, 
        const BamTools::RefVector& refs, 
        const int32_t totalTail, 
        const int32_t critTail, 
        const bool diff_ref = true);

int32_t 
checkLinkPair(const BamAlignment& al1,
        const BamAlignment& al2, 
        const BamTools::RefVector& refs, 
        const int32_t totalTail, 
        const int32_t critTail, 
        const bool diff_ref = true);

int32_t
checkLinkPairCandidate(const BamAlignment&, 
        const BamTools::RefVector&, 
        const int32_t critTail);

int32_t
readTail(const BamAlignment& al, 
        const BamTools::RefVector& refs);

int32_t
readTailS(const bool mapped, const bool rev, const int32_t pos, 
        const int32_t ref_len, const int32_t aligned_len);

}  // namespace ibeji


#endif /* PROCESSREADPAIR_H_ */

