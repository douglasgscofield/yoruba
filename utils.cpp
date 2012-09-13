/*
 * utils.cpp
 *
 * Utility functions.
 *
 *  Created on: Apr 9, 2012
 *      Author: douglasgscofield
 */

#include "utils.h"

using namespace std;
using namespace BamTools;
using namespace ibeji;

namespace ibeji {

void 
printAlignmentInfo(const BamAlignment& al, const RefVector& refs, int32_t level) {
    cout << al.Name;
    if (level > 1)
        cout << "\tprim=" << al.IsPrimaryAlignment();
    cout << "\trefid=" << al.RefID
        << "\tref=" << refs[al.RefID].RefName
        << "\tlen=" << refs[al.RefID].RefLength 
        << "\tpos=" << al.Position;
    if (level > 0) 
        cout << "\tmap=" << al.IsMapped();
    cout << "\trev=" << al.IsReverseStrand();
    if (level > 1) 
        cout << " |"
            << "\tm_refid=" << al.MateRefID
            << "\tm_pos=" << al.MatePosition
            << "\tm_map=" << al.IsMateMapped() 
            << "\tm_rev=" << al.IsMateReverseStrand();
    if (level > 0) 
        cout << "\tq_len=" << al.QueryBases.length()
            << "\ta_len=" << al.AlignedBases.length();
    cout << "\tmapq=" << al.MapQuality;
    if (level > 2)
        cout << "\tpair=" << al.IsPaired() 
            << "\tprop=" << al.IsProperPair() 
            << "\tm_1=" << al.IsFirstMate() 
            << "\tm_2=" << al.IsSecondMate() 
            << "\tisz=" << setw(7) << al.InsertSize;
    cout << endl;
}

void
printAlignmentInfo_fields(const BamAlignment& al, const RefVector& refs, int32_t level) {
    cout << setw(35) << left << al.Name << right;
    if (level > 1)
        cout << " prim " << al.IsPrimaryAlignment();
    cout << " |"
        << " refid " << setw(8) << al.RefID
        << " <" << setw(15) << refs[al.RefID].RefName
        << ",L" << setw(5) << refs[al.RefID].RefLength << ">" 
        << " pos " << setw(8) << al.Position;
    if (level > 0) 
        cout << " map " << al.IsMapped();
    cout << " rev " << al.IsReverseStrand();
    if (level > 1) 
        cout << " |"
            << " MATE refid " << setw(8) << al.MateRefID
            << " pos " << setw(8) << al.MatePosition
            << " map " << al.IsMateMapped() 
            << " rev " << al.IsMateReverseStrand();
    cout << " |";
    if (level > 0) 
        cout << " q:aln " << setw(3) << al.QueryBases.length()
            << ":" << setw(3) << al.AlignedBases.length();
    cout << " mapq " << al.MapQuality;
    if (level > 2)
        cout << " pair " << al.IsPaired() 
            << " proppair " << al.IsProperPair() 
            << " |"
            << " mat1st " << al.IsFirstMate() 
            << " mat2nd " << al.IsSecondMate() 
            << " isize " << setw(7) << al.InsertSize;
    cout << endl;
}


// Spit out basic BamAlignment data, from early BamReader code
void 
PrintAlignment(const BamAlignment& alignment) {
	cout << "---------------------------------" << endl;
	cout << "Name: "       << alignment.Name << endl;
	cout << "Aligned to: " << alignment.RefID;
	cout << ":"            << alignment.Position << endl;
}


}  // namespace ibeji
