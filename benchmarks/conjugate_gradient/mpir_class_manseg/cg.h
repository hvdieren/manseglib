#ifndef __cg_CG_H__
#define __cg_CG_H__

#ifndef FLOAT
#define FLOAT float
#endif

#ifndef DOUBLE
#define DOUBLE double
#endif

#ifndef FLOAT2
#define FLOAT2 DOUBLE
#endif

#ifndef MANSEG_LIB
#define MANSEG_LIB
#include "../../../mantissaSegmentation.hpp"

using HeadsArray = ManSeg::TwoSegArray<false>;
using PairsArray = ManSeg::TwoSegArray<true>;

using namespace ManSeg;

#endif

#define ALLOC(t, l) (t*)malloc((l) * sizeof(t))
#define CALLOC(t, l) (t*)calloc(sizeof(t),(l))
#define FREE(p) *p = INFINITY; free(p);

#endif
