/*
 * processReadPair.cpp
 *
 *  Created on: Apr 9, 2012
 *      Author: douglasgscofield
 */

using namespace std;
#include "api/BamReader.h"
using namespace BamTools;
#include "processReadPair.h"
#include "ibejiAlignment.h"
using namespace ibeji;

namespace ibeji {

const bool debug_processReadPair = true;
const bool debug_checkLinkPairCandidate = false;
const bool debug_checkLinkPair = true;
const bool debug_readTail = false;

bool 
processReadPair(const BamAlignment& al1, 
        const BamAlignment& al2, 
        const RefVector& refs, 
        const int32_t totalTail, 
        const int32_t critTail, 
        const bool diff_ref)
{
    if ((al1.IsFirstMate() && al2.IsFirstMate())
        || (al1.IsSecondMate() && al2.IsSecondMate())) {
        cerr << "Incompatible mate orders: name1 = " << al1.Name 
             << " is1stmate " << al1.IsFirstMate() << " is2ndmate " << al1.IsSecondMate()
             << " name2 = " << al2.Name 
             << " is1stmate " << al2.IsFirstMate() << " is2ndmate " << al2.IsSecondMate()
             << endl;
        exit(1);
    }

    int32_t total_tail = -1;
    if (! (total_tail = checkLinkPair(al1, al2, refs, totalTail, critTail, diff_ref))) {
        return false;  // reject all but link pairs
        // continue;
    }
    if (critTail && ! checkLinkPairCandidate(al1, refs, critTail)
        && ! checkLinkPairCandidate(al2, refs, critTail)) {
        return false;  // neither read was a link pair candidate
    }
    if (debug_processReadPair) cout << "---------------------------------" << endl;
    int32_t lpc_tail1 = checkLinkPairCandidate(al1, refs, critTail);
    int32_t lpc_tail2 = checkLinkPairCandidate(al2, refs, critTail);
    if (debug_processReadPair) {
        printAlignmentInfo(al1, refs);
        if (lpc_tail1) {
            cout << "LINK PAIR CANDIDATE ";
            cout << ((lpc_tail1 > 0) ? "--->" : "<---") << " " << lpc_tail1 << endl;
        }
        printAlignmentInfo(al2, refs);
        if (lpc_tail2) {
            cout << "LINK PAIR CANDIDATE ";
            cout << ((lpc_tail2 > 0) ? "--->" : "<---") << " " << lpc_tail2 << endl;
        }
        cout << "TOTAL TAIL " << (abs(readTail(al1, refs)) + abs(readTail(al2, refs))) << endl;
    }

    return true;
}


int32_t
readTail(const BamAlignment& al, 
        const RefVector& refs)
{
    return readTailS(al.IsMapped(), al.IsReverseStrand(), al.Position, 
                     refs[al.RefID].RefLength, al.AlignedBases.length());
}

int32_t
readTailS(const bool mapped, const bool rev, const int32_t pos, 
        const int32_t ref_len, const int32_t aligned_len)
{
    if (! mapped) return 0;
    int32_t tail = 0;
    if (rev) {
        tail = -(pos + aligned_len);
        if (debug_readTail) cout << "readTailS: - tail = " << tail << endl;
    } else {
        tail = ref_len - pos + 1;
        if (debug_readTail) cout << "readTailS: + tail = " << tail << endl;
    }
    return tail;
}


int32_t 
checkLinkPair(const BamAlignment& al1,
        const BamAlignment& al2, 
        const RefVector& refs, 
        const int32_t totalTail, 
        const int32_t critTail, 
        const bool diff_ref)
{
    if (al1.RefID == al2.RefID && diff_ref) return 0;
    int32_t tail1 = checkLinkPairCandidate(al1, refs, critTail);
    if (! tail1) return 0;
    int32_t tail2 = checkLinkPairCandidate(al2, refs, critTail);
    if (! tail2) return 0;
    int32_t total_tail = (abs(tail1) + abs(tail2));
    if (total_tail > totalTail) return 0;
    return (total_tail);
}


int32_t
checkLinkPairCandidate(const BamAlignment& al, 
        const RefVector& refs, 
        const int32_t critTail)
{
    int32_t tail = readTail(al, refs);
    return abs(tail) <= critTail ? tail : 0;
}


// old mate-finding code, for position-sorted BAM files, didn't work well
BamAlignment 
lookForMate(BamReader& rdr, BamAlignment& al, RefVector& refs)
{
    if (! rdr.Jump(al.MateRefID, al.MatePosition)) {
        cout << "*** Could not jump to " << al.MateRefID << ":" << al.MatePosition << endl;
        exit(1);
    }
    BamAlignment al_jump;
    while (rdr.GetNextAlignment(al_jump)) {
        if (al_jump.Name == al.Name 
            && al_jump.RefID == al.MateRefID
            && al_jump.Position == al.MatePosition) {
         cout << "MATE FOUND" << endl;
            printAlignmentInfo(al_jump, refs);
            break;
        } else if (al_jump.Position > al.MatePosition) {
         cout << "NO MATE FOUND, beyond MatePosition" << endl;
            break;
        }
    }
    if (! rdr.Jump(al.RefID, al.Position)) {
        cout << "*** Could not return to " << al.RefID << ":" << al.Position << endl;
        exit(1);
    }
    // need to return to the read we were on after the jump, but how?
    rdr.GetNextAlignment(al);
    return(al_jump);
}


}  // namespace ibeji;


