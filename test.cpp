#include <stdlib.h>
#include <iostream>
#include <vector>
#include <map>
#include <tr1/unordered_map>

using namespace std;

typedef std::tr1::unordered_map<string, long> dupMap1_t;
typedef map<size_t, std::tr1::unordered_map<string, long> > dupMap2_t;

#include "MemTrack.h"
using namespace MemTrack;

int
main(int argc, char* argv[]) {
    dupMap1_t dupMap1;
    dupMap2_t dupMap2;

    cerr << "single map" << endl;
    for (long l = 0; l <= 100; ++l) {
        char buf[100];
        sprintf(buf, "chr%07ld", l);
        dupMap1[buf] = l;
        for (int i = 0; i <= 10000; ++i) {
            sprintf(buf, "%s.%d", buf, i);
            dupMap1[buf] = long(i);
        }
    }
    TrackListMemoryUsage();
    cerr << "adding map of maps" << endl;
    for (long l = 0; l <= 100; ++l) {
        char buf[100];
        sprintf(buf, "chr%07ld", l);
        dupMap2[l][buf] = l;
        for (int i = 0; i <= 10000; ++i) {
            sprintf(buf, "%s.%d", buf, i);
            dupMap2[l][buf] = long(i);
        }
    }
    TrackListMemoryUsage();
}
