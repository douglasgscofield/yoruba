
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
               const string& prefix,
               const string& rg_prefix,
               const string& rg_sep,
               const string& rg_delim,
               const string& rg_terminate)
{
    if (rgd.IsEmpty()) {
        os << prefix << "read group dictionary is empty" << std::endl;
        return;
    }
    SamReadGroupConstIterator rgdI = rgd.ConstBegin();
    SamReadGroupConstIterator rgd_End = rgd.ConstEnd();
    for (rgdI = rgd.ConstBegin(); rgdI != rgd_End; ++rgdI) {
        os << prefix;
        yoruba::printReadGroup(os, *rgdI, rg_prefix, rg_sep, rg_delim, rg_terminate);
    }
}


//-------------------------------------


void
yoruba::printReadGroup(ostream& os,
               const SamReadGroup& rg,
               const string& rg_prefix,
               const string& rg_sep,
               const string& rg_delim,
               const string& rg_terminate)
{
    os << rg_prefix;
    if (rg.HasID())
        os << rg_sep << "ID:" << rg_delim << rg.ID << rg_delim;
    if (rg.HasSequencingCenter())
        os << rg_sep << "CN:" << rg_delim << rg.SequencingCenter << rg_delim;
    if (rg.HasDescription())
        os << rg_sep << "DS:" << rg_delim << rg.Description << rg_delim;
    if (rg.HasProductionDate())
        os << rg_sep << "DT:" << rg_delim << rg.ProductionDate << rg_delim;
    if (rg.HasFlowOrder())
        os << rg_sep << "FO:" << rg_delim << rg.FlowOrder << rg_delim;
    if (rg.HasKeySequence())
        os << rg_sep << "KS:" << rg_delim << rg.KeySequence << rg_delim;
    if (rg.HasLibrary())
        os << rg_sep << "LB:" << rg_delim << rg.Library << rg_delim;
    if (rg.HasProgram())
        os << rg_sep << "PG:" << rg_delim << rg.Program << rg_delim;
    if (rg.HasPredictedInsertSize())
        os << rg_sep << "PI:" << rg_delim << rg.PredictedInsertSize << rg_delim;
    if (rg.HasSequencingTechnology())
        os << rg_sep << "PL:" << rg_delim << rg.SequencingTechnology << rg_delim;
    if (rg.HasPlatformUnit())
        os << rg_sep << "PU:" << rg_delim << rg.PlatformUnit << rg_delim;
    if (rg.HasSample())
        os << rg_sep << "SM:" << rg_delim << rg.Sample << rg_delim;
    os << rg_terminate;
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


