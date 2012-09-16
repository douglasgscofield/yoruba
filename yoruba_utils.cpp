
#include "yoruba_utils.h"

using namespace std;
using namespace BamTools;
using namespace yoruba;

namespace yoruba {

void 
printAlignmentInfo(std::ostream& os, const BamAlignment& al, int32_t level) {
    os << al.Name;
    if (level > 1)
        os << "\tprim=" << al.IsPrimaryAlignment();
    os << "\trefid=" << al.RefID
        << "\tpos=" << al.Position;
    if (level > 0) 
        os << "\tmap=" << al.IsMapped();
    os << "\trev=" << al.IsReverseStrand();
    if (level > 1) 
        os << " |"
            << "\tm_refid=" << al.MateRefID
            << "\tm_pos=" << al.MatePosition
            << "\tm_map=" << al.IsMateMapped() 
            << "\tm_rev=" << al.IsMateReverseStrand();
    if (level > 0) 
        os << "\tq_len=" << al.QueryBases.length()
            << "\ta_len=" << al.AlignedBases.length();
    os << "\tmapq=" << al.MapQuality;
    if (level > 2)
        os << "\tpair=" << al.IsPaired() 
            << "\tprop=" << al.IsProperPair() 
            << "\tm_1=" << al.IsFirstMate() 
            << "\tm_2=" << al.IsSecondMate() 
            << "\tisz=" << setw(7) << al.InsertSize;
    // tags
    string tag;
    if (al.GetTag("RG", tag))
        os << "\tRG:Z:'" << tag << "'" << endl;
}

void 
printAlignmentInfo(std::ostream& os, const BamAlignment& al, const RefVector& refs, int32_t level) {
    os << al.Name;
    if (level > 1)
        os << "\tprim=" << al.IsPrimaryAlignment();
    os << "\trefid=" << al.RefID
        << "\tref=" << refs[al.RefID].RefName
        << "\tlen=" << refs[al.RefID].RefLength 
        << "\tpos=" << al.Position;
    if (level > 0) 
        os << "\tmap=" << al.IsMapped();
    os << "\trev=" << al.IsReverseStrand();
    if (level > 1) 
        os << " |"
            << "\tm_refid=" << al.MateRefID
            << "\tm_pos=" << al.MatePosition
            << "\tm_map=" << al.IsMateMapped() 
            << "\tm_rev=" << al.IsMateReverseStrand();
    if (level > 0) 
        os << "\tq_len=" << al.QueryBases.length()
            << "\ta_len=" << al.AlignedBases.length();
    os << "\tmapq=" << al.MapQuality;
    if (level > 2)
        os << "\tpair=" << al.IsPaired() 
            << "\tprop=" << al.IsProperPair() 
            << "\tm_1=" << al.IsFirstMate() 
            << "\tm_2=" << al.IsSecondMate() 
            << "\tisz=" << setw(7) << al.InsertSize;
    // tags
    string tag;
    if (al.GetTag("RG", tag))
        os << "\tRG:Z:'" << tag << "'" << endl;
}

void
printAlignmentInfo_fields(std::ostream& os, const BamAlignment& al, int32_t level) {
    os << setw(35) << left << al.Name << right;
    if (level > 1)
        os << " prim " << al.IsPrimaryAlignment();
    os << " |"
        << " refid " << setw(8) << al.RefID
        << " pos " << setw(8) << al.Position;
    if (level > 0) 
        os << " map " << al.IsMapped();
    os << " rev " << al.IsReverseStrand();
    if (level > 1) 
        os << " |"
            << " MATE refid " << setw(8) << al.MateRefID
            << " pos " << setw(8) << al.MatePosition
            << " map " << al.IsMateMapped() 
            << " rev " << al.IsMateReverseStrand();
    os << " |";
    if (level > 0) 
        os << " q:aln " << setw(3) << al.QueryBases.length()
            << ":" << setw(3) << al.AlignedBases.length();
    os << " mapq " << al.MapQuality;
    if (level > 2)
        os << " pair " << al.IsPaired() 
            << " proppair " << al.IsProperPair() 
            << " |"
            << " mat1st " << al.IsFirstMate() 
            << " mat2nd " << al.IsSecondMate() 
            << " isize " << setw(7) << al.InsertSize;
    // tags
    os << " | tags";
    string tag;
    if (al.GetTag("RG", tag))
        os << " RG:Z:'" << tag << "'" << endl;
}

void
printAlignmentInfo_fields(std::ostream& os, const BamAlignment& al, const RefVector& refs, int32_t level) {
    os << setw(35) << left << al.Name << right;
    if (level > 1)
        os << " prim " << al.IsPrimaryAlignment();
    os << " |"
        << " refid " << setw(8) << al.RefID
        << " <" << setw(15) << refs[al.RefID].RefName
        << ",L" << setw(5) << refs[al.RefID].RefLength << ">" 
        << " pos " << setw(8) << al.Position;
    if (level > 0) 
        os << " map " << al.IsMapped();
    os << " rev " << al.IsReverseStrand();
    if (level > 1) 
        os << " |"
            << " MATE refid " << setw(8) << al.MateRefID
            << " pos " << setw(8) << al.MatePosition
            << " map " << al.IsMateMapped() 
            << " rev " << al.IsMateReverseStrand();
    os << " |";
    if (level > 0) 
        os << " q:aln " << setw(3) << al.QueryBases.length()
            << ":" << setw(3) << al.AlignedBases.length();
    os << " mapq " << al.MapQuality;
    if (level > 2)
        os << " pair " << al.IsPaired() 
            << " proppair " << al.IsProperPair() 
            << " |"
            << " mat1st " << al.IsFirstMate() 
            << " mat2nd " << al.IsSecondMate() 
            << " isize " << setw(7) << al.InsertSize;
    // tags
    os << " | tags";
    string tag;
    if (al.GetTag("RG", tag))
        os << " RG:Z:'" << tag << "'" << endl;
}


// Spit out basic BamAlignment data, from early BamReader code
void 
PrintAlignment(std::ostream& os, const BamAlignment& alignment) {
	os << "---------------------------------" << endl;
	os << "Name: "       << alignment.Name << endl;
	os << "Aligned to: " << alignment.RefID;
	os << ":"            << alignment.Position << endl;
}


}  // namespace yoruba

