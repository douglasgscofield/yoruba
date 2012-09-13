/*
 * ibejiAlignment.h
 *
 * Map a read name to a BamAlignment.  Because of map semantics, only
 * one alignment can be kept per read name.  To hold a pair, you must
 * have two alignmentMap instantiations, one for each member of the
 * pair attached to the same read name.
 *
 *  Created on: Apr 9, 2012
 *      Author: douglasgscofield
 */

#ifndef IBEJIALIGNMENT_H_
#define IBEJIALIGNMENT_H_


#include <map>
#include <string>
#include <functional>
#include "api/BamAlignment.h"

namespace yoruba {

// first, a map using the BamAlignment class.

typedef std::map<std::string, BamTools::BamAlignment> alignmentMapFull;
typedef alignmentMapFull::iterator alignmentMapFullI;

// second, a light alignment class

class ibejiAlignment {

    public:
        ibejiAlignment(void);
        ibejiAlignment(const ibejiAlignment& other);
        ibejiAlignment(const BamTools::BamAlignment& otherBA);
        ~ibejiAlignment(void);
    public:
        std::string Name;
        int32_t     Length;
        int32_t     RefID;
        int32_t     Position;

        class { // so we have an AlignedBases.length() method
            public:
                int32_t     AlignedLength;
                int32_t length(void) { return AlignedLength; }
        } AlignedBases;
        bool        IsMapped(void) { return Mapped; };
        bool        IsReverseStrand(void) { return ReverseStrand; };
        bool        IsMateMapped(void) { return MateMapped; };
        bool        IsMateReverseStrand(void) { return MateReverseStrand; };
        bool        IsFirstMate(void) { return FirstMate; };
        bool        IsSecondMate(void) { return SecondMate; };
    private:
        bool        Mapped;
        bool        ReverseStrand;
    public:
        int32_t     MateRefID;
        int32_t     MatePosition;
    private:
        bool        MateMapped;
        bool        MateReverseStrand;
        bool        FirstMate;
        bool        SecondMate;

};  // class ibejiAlignment

typedef std::map<std::string, ibejiAlignment> alignmentMapLite;
typedef alignmentMapLite::iterator alignmentMapLiteI;

//typedef ibejiAlignment ibejiAlignment;
//typedef std::map<std::string, ibejiAlignment> alignmentMapLite;
//typedef alignmentMapLite::iterator alignmentMapLiteI;
//typedef std::map<std::string, ibejiAlignment> alignmentMap;
typedef std::map<std::string, BamTools::BamAlignment> alignmentMap;
typedef alignmentMap::iterator alignmentMapI;

};  // namespace yoruba

#endif /* IBEJIALIGNMENT_H_ */

