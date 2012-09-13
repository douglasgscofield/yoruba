/*
 * ibejiAlignment.cpp
 *
 *  Created on: Apr 10, 2012
 *      Author: douglasgscofield
 */


#include "ibejiAlignment.h"
using namespace std;
using namespace BamTools;
using namespace yoruba;

namespace yoruba {


ibejiAlignment::ibejiAlignment(void)
    : Length(0)
    , RefID(-1)
    , Position(-1)
    //, AlignedBases.AlignedLength(-1)
    , Mapped(false)
    , ReverseStrand(false)
    , MateRefID(-1)
    , MatePosition(-1)
    , MateMapped(false)
    , MateReverseStrand(false)
    , FirstMate(false)
    , SecondMate(false)
//{ }
{ AlignedBases.AlignedLength = -1; }

ibejiAlignment::ibejiAlignment(const ibejiAlignment& other)
    : Name(other.Name)
    , Length(other.Length)
    , RefID(other.RefID)
    , Position(other.Position)
    // , AlignedBases.AlignedLength(other.AlignedLength)
    , Mapped(other.Mapped)
    , ReverseStrand(other.ReverseStrand)
    , MateRefID(other.MateRefID)
    , MatePosition(other.MatePosition)
    , MateMapped(other.MateMapped)
    , MateReverseStrand(other.MateReverseStrand)
    , FirstMate(other.FirstMate)
    , SecondMate(other.SecondMate)
//{ }
{ AlignedBases.AlignedLength = other.AlignedBases.AlignedLength; }

ibejiAlignment::ibejiAlignment(const BamAlignment& otherBA)
    : Name(otherBA.Name)
    , Length(otherBA.Length)
    , RefID(otherBA.RefID)
    , Position(otherBA.Position)
    //, AlignedBases.AlignedLength(otherBA.AlignedBases.length())
    , Mapped(otherBA.IsMapped())
    , ReverseStrand(otherBA.IsReverseStrand())
    , MateRefID(otherBA.MateRefID)
    , MatePosition(otherBA.MatePosition)
    , MateMapped(otherBA.IsMateMapped())
    , MateReverseStrand(otherBA.IsMateReverseStrand())
    , FirstMate(otherBA.IsFirstMate())
    , SecondMate(otherBA.IsSecondMate())
//{ }
{ AlignedBases.AlignedLength = otherBA.AlignedBases.length(); }

ibejiAlignment::~ibejiAlignment(void)
{ }


}; // namespace yoruba


