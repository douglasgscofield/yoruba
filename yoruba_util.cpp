
#include "yoruba.h"
#include "yoruba_util.h"

using namespace std;
using namespace BamTools;
using namespace yoruba;


//-------------------------------------


bool
yoruba::isCoordinateSorted(int32_t ref, int32_t pos, int32_t prev_ref, int32_t prev_pos)
{
    if (ref < prev_ref || (ref == prev_ref && pos < prev_pos))
        return false;
    return true;
}


//-------------------------------------


bool
yoruba::isMateUpstream(const BamAlignment& alignment)
{
    // assumes coordinate-sorted
    if (! alignment.IsPaired() || ! alignment.IsMateMapped()) {
        return false;
    } else if (alignment.RefID < alignment.MateRefID) {
        return false;
    } else if (alignment.RefID == alignment.MateRefID) {
        if (alignment.InsertSize > 0) {
            return false;
        } else if (alignment.InsertSize == 0) {
            cerr << "yoruba::isMateUpstream(): InsertSize == 0, strangely..." << endl;
            return false;
        } else if (alignment.InsertSize < 0) {
            return true;
        }
    } else if (alignment.RefID > alignment.MateRefID) {
        return true;
    }
    cerr << "yoruba::isMateUpstream(): unhandled case, strangely..." << endl;
    return false;
}


//-------------------------------------


bool
yoruba::isMateDownstream(const BamAlignment& alignment)
{
    // assumes coordinate-sorted
    if (! alignment.IsPaired() || ! alignment.IsMateMapped()) {
        return false;
    } else if (alignment.RefID < alignment.MateRefID) {
        return true;
    } else if (alignment.RefID == alignment.MateRefID) {
        if (alignment.InsertSize > 0) {
            return true;
        } else if (alignment.InsertSize == 0) {
            cerr << "yoruba::isMateDownstream(): InsertSize == 0, strangely..." << endl;
            return false;
        } else if (alignment.InsertSize < 0) {
            return false;
        }
    } else if (alignment.RefID > alignment.MateRefID) {
        return false;
    }
    cerr << "yoruba::isMateDownstream(): unhandled case, strangely..." << endl;
    return false;
}


//-------------------------------------


// overloaded
void
yoruba::PrintAlignment(const BamAlignment& alignment)
{
    yoruba::PrintAlignment(cerr, alignment);
}


//-------------------------------------


// overloaded
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


// overloaded
void
yoruba::printAlignmentInfo(std::ostream& os,
                           const BamAlignment& al,
                           int32_t level)
{
    const RefVector dummy_refs;
    yoruba::printAlignmentInfo(os, al, dummy_refs, level);
}


//-------------------------------------


// overloaded
void
yoruba::printAlignmentInfo(std::ostream& os,
                           const BamAlignment& al,
                           const RefVector& refs,
                           int32_t level)
{
    os << al.Name;
    if (al.IsDuplicate()) os << "\tDuplicate";
    if (level > 1)
        if (al.IsPrimaryAlignment()) os << "\tPrimary";
    os << (al.IsMapped() ? "\tMapped" : "\tUnmapped");
    os << "\tRefID=" << al.RefID;
    if (al.IsMapped() && al.RefID >= 0 && ! refs.empty()) {
        os << "[" << refs[al.RefID].RefName << ",l=" << refs[al.RefID].RefLength << "]";
    }
    os << ":Pos=" << al.Position;
    os << "\tmapQ=" << al.MapQuality;
    os << (al.IsReverseStrand() ?  "\tRev" : "\tForw");
    if (level > 0)
        os << "\tQ,A=" << al.QueryBases.length() << "," << al.AlignedBases.length();
    os << "\tPair=" << al.IsPaired();
    if (al.IsPaired() && level > 1) {
        os << " |";
        os << (al.IsMateMapped() ? "\tmMapped" : "\tmUnmapped");
        os << "\tmRefID=" << al.MateRefID;
        if (al.IsMateMapped() && al.MateRefID >= 0 && ! refs.empty()) {
            os << "[" << refs[al.MateRefID].RefName << ",l=" << refs[al.MateRefID].RefLength << "]";
        }
        os << ":mPos=" << al.MatePosition;
        os << (al.IsMateReverseStrand() ?  "\tmRev" : "\tmForw");
    }
    if (al.IsPaired() && level > 2) {
        os << "\t";
        if (al.IsProperPair()) os << "PropPair,";
        if (al.IsFirstMate()) os << "1stMate";
        if (al.IsSecondMate()) os << "2ndMate";
        os << "\tISize=" << al.InsertSize;
    }
    // tags
    string tag;
    if (al.GetTag("RG", tag))
        os << "\tRG:Z:" << tag;
    os << endl;
}


//-------------------------------------


// overloaded
void
yoruba::printAlignmentInfo_fields(std::ostream& os,
                                  const BamAlignment& al,
                                  int32_t level)
{
    const RefVector dummy_refs;
    yoruba::printAlignmentInfo_fields(os, al, dummy_refs, level);
}


//-------------------------------------


// overloaded
void
yoruba::printAlignmentInfo_fields(std::ostream& os,
                                  const BamAlignment& al,
                                  const RefVector& refs,
                                  int32_t level) {
    os << setw(35) << left << al.Name << right;
    if (al.IsDuplicate()) os << " Dup";
    if (level > 1)
        if (al.IsPrimaryAlignment()) os << " Prim";
    os << (al.IsMapped() ? " Map" : " Unmap");
    os << " |";
    os << " RefID " << setw(8) << al.RefID;
    if (al.IsMapped() && al.RefID >= 0 && ! refs.empty()) {
        os << " [" << setw(15) << refs[al.RefID].RefName << "," << 
            setw(5) << refs[al.RefID].RefLength << "]";
    }
    os << " Pos " << setw(8) << al.Position;
    os << " mapQ " << al.MapQuality;
    os << (al.IsReverseStrand() ? " Rev" : " Forw");
    if (level > 0)
        os << " Q,A " << setw(3) << al.QueryBases.length() << "," << 
            setw(3) << al.AlignedBases.length();
    os << " Pair" << al.IsPaired();
    if (al.IsPaired() && level > 1) {
        os << " |";
        os << (al.IsMateMapped() ? " Map" : " Unmap");
        os << " mRefID " << setw(8) << al.MateRefID;
        if (al.IsMateMapped() && al.MateRefID >= 0 && ! refs.empty()) {
            os << "[" << refs[al.MateRefID].RefName << ",l=" << refs[al.MateRefID].RefLength << "]";
        }
        os << " mPos " << setw(8) << al.MatePosition;
        os << (al.IsMateReverseStrand() ? " mRev" : " mForw");
    }
    if (al.IsPaired() && level > 2) {
        os << " | ";
        if (al.IsProperPair()) os << "PropPair,";
        if (al.IsFirstMate()) os << "Mate1," << al.IsFirstMate();
        if (al.IsSecondMate()) os << "Mate2" << al.IsSecondMate();
        os << " ISz " << al.InsertSize;
    }
    // tags
    os << " | tags";
    string tag;
    if (al.GetTag("RG", tag))
        os << " RG:Z:" << tag;
    os << endl;
}


