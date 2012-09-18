
#include "yoruba.h"
#include "yoruba_util.h"

using namespace std;
using namespace BamTools;
using namespace yoruba;


//-------------------------------------


void 
yoruba::PrintAlignment(const BamAlignment& alignment)
{
    yoruba::PrintAlignment(cerr, alignment);
}


void 
yoruba::PrintAlignment(ostream& os, 
                       const BamAlignment& alignment)
{
	os << "---------------------------------" << endl;
	os << "Name: "       << alignment.Name << endl;
	os << "Aligned to: " << alignment.RefID;
	os << ":"            << alignment.Position << endl;
}


//-------------------------------------


void
yoruba::printReadGroupDictionary(ostream& os, 
                                 const SamReadGroupDictionary& rgd,
                                 const string& prefix)
{
    if (rgd.IsEmpty()) {
        os << prefix << "read group dictionary is empty" << std::endl;
        return;
    }
    int i = 0;
    string new_prefix;
    for (SamReadGroupConstIterator rgdI = rgd.ConstBegin(); 
            rgdI != rgd.ConstEnd(); ++rgdI, ++i) {
        char new_prefix[100];
        sprintf(new_prefix, "%s@RG %d ", prefix.c_str(), i);
        yoruba::printReadGroup(os, *rgdI, new_prefix);
    }
}


//-------------------------------------


void
yoruba::printReadGroup(ostream& os, 
                       const SamReadGroup& rg,
                       const string& prefix)
{
    if (rg.HasID())
        os << prefix << "ID:'" << rg.ID << "'" << std::endl;
    if (rg.HasSequencingCenter()) 
        os << prefix << "CN:'" << rg.SequencingCenter << "'" << std::endl;
    if (rg.HasDescription())
        os << prefix << "DS:'" << rg.Description << "'" << std::endl;
    if (rg.HasProductionDate())   
        os << prefix << "DT:'" << rg.ProductionDate << "'" << std::endl;
    if (rg.HasFlowOrder()) 
        os << prefix << "FO:'" << rg.FlowOrder << "'" << std::endl;
    if (rg.HasKeySequence()) 
        os << prefix << "KS:'" << rg.KeySequence << "'" << std::endl;
    if (rg.HasLibrary()) 
        os << prefix << "LB:'" << rg.Library << "'" << std::endl;
    if (rg.HasProgram()) 
        os << prefix << "PG:'" << rg.Program << "'" << std::endl;
    if (rg.HasPredictedInsertSize()) 
        os << prefix << "PI:'" << rg.PredictedInsertSize << "'" << std::endl;
    if (rg.HasSequencingTechnology()) 
        os << prefix << "PL:'" << rg.SequencingTechnology << "'" << std::endl;
    if (rg.HasPlatformUnit()) 
        os << prefix << "PU:'" << rg.PlatformUnit << "'" << std::endl;
    if (rg.HasSample()) 
        os << prefix << "SM:'" << rg.Sample << "'" << std::endl;
}


//-------------------------------------


void 
yoruba::printAlignmentInfo(std::ostream& os, 
                           const BamAlignment& al, 
                           int32_t level)
{
    const RefVector dummy_refs;
    yoruba::printAlignmentInfo(os, al, dummy_refs, level);
}


//-------------------------------------


void 
yoruba::printAlignmentInfo(std::ostream& os, 
                           const BamAlignment& al, 
                           const RefVector& refs, 
                           int32_t level)
{
    os << al.Name;
    if (level > 1)
        os << "\tprim=" << al.IsPrimaryAlignment();
    os << "\trefid=" << al.RefID;
    if (al.IsMapped() && ! refs.empty()) {
        os << "\tref=" << refs[al.RefID].RefName
           << "\tlen=" << refs[al.RefID].RefLength;
    }
    os << "\tpos=" << al.Position;
    if (level > 0) 
        os << "\tmap=" << al.IsMapped();
    os << "\trev=" << al.IsReverseStrand();
    if (al.IsPaired() && level > 1) 
        os << " |"
            << "\tm_refid=" << al.MateRefID
            << "\tm_pos=" << al.MatePosition
            << "\tm_map=" << al.IsMateMapped() 
            << "\tm_rev=" << al.IsMateReverseStrand();
    if (level > 0) 
        os << "\tq_len=" << al.QueryBases.length()
            << "\ta_len=" << al.AlignedBases.length();
    os << "\tmapq=" << al.MapQuality;
    os << "\tpair=" << al.IsPaired();
    if (al.IsPaired() && level > 2)
        os  << "\tprop=" << al.IsProperPair() 
            << "\tm_1=" << al.IsFirstMate() 
            << "\tm_2=" << al.IsSecondMate() 
            << "\tisz=" << al.InsertSize;
    // tags
    string tag;
    if (al.GetTag("RG", tag))
        os << "\tRG:Z:'" << tag << "'";
    os << endl;
}


//-------------------------------------


void
yoruba::printAlignmentInfo_fields(std::ostream& os, 
                                  const BamAlignment& al, 
                                  int32_t level)
{
    const RefVector dummy_refs;
    yoruba::printAlignmentInfo_fields(os, al, dummy_refs, level);
}


//-------------------------------------


void
yoruba::printAlignmentInfo_fields(std::ostream& os, 
                                  const BamAlignment& al, 
                                  const RefVector& refs, 
                                  int32_t level) {
    os << setw(35) << left << al.Name << right;
    if (level > 1)
        os << " prim " << al.IsPrimaryAlignment();
    os << " |"
        << " refid " << setw(8) << al.RefID;
    if (al.IsMapped() && ! refs.empty()) 
        os << " <" << setw(15) << refs[al.RefID].RefName
           << ",L" << setw(5) << refs[al.RefID].RefLength << ">";
    os << " pos " << setw(8) << al.Position;
    if (level > 0) 
        os << " map " << al.IsMapped();
    os << " rev " << al.IsReverseStrand();
    if (al.IsPaired() && level > 1) 
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
    os << " pair " << al.IsPaired();
    if (al.IsPaired() && level > 2)
        os  << " proppair " << al.IsProperPair() 
            << " |"
            << " mat1st " << al.IsFirstMate() 
            << " mat2nd " << al.IsSecondMate() 
            << " isize " << al.InsertSize;
    // tags
    os << " | tags";
    string tag;
    if (al.GetTag("RG", tag))
        os << " RG:Z:'" << tag << "'";
    os << endl;
}


