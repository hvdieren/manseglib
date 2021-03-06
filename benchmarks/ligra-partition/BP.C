#define  MORE_ARG 1
#include "ligra-numa.h"
#include "math.h"
int maxIter=20;
using namespace std;

#define NSTATES 2
// these structs are what you'd replace with
// ManSegArrays
// the 2d array would be interesting, though
// in theory, should be fine - just an array of ManSegArrays
struct EdgeWeight
{
    double potential[NSTATES][NSTATES];
};

struct EdgeData
{
    double belief[NSTATES];
};

struct VertexInfo
{
    double potential[NSTATES];
};

struct VertexData
{
    double product[NSTATES];
};

template <class ET>
inline void writeDiv(ET *a, ET b)
{
    volatile ET newV, oldV;
    do
    {
        oldV = *a;
        newV = oldV / b;
    }
    while (!CAS(a, oldV, newV));
}

template <class ET>
inline void writeMult(ET *a, ET b)
{
    volatile ET newV, oldV;
    do
    {
        oldV = *a;
        newV = oldV * b;
    }
    while (!CAS(a, oldV, newV));
}

template <class vertex>
struct BP_F
{
    EdgeWeight *edgeW;
    EdgeData *edgeD_curr;
    EdgeData *edgeD_next;
    VertexInfo *vertI;
    VertexData *vertD_curr;
    VertexData *vertD_next;
    intT *offsets;
    static const bool use_cache = true;
    struct cache_t
    {  
        VertexData vertD_cache_next;	
    };
    BP_F(EdgeWeight *_edgeW, EdgeData *_edgeD_curr, EdgeData *_edgeD_next, VertexInfo *_vertI, VertexData *_vertD_curr, VertexData *_vertD_next, intT *_offsets) :
        edgeW(_edgeW), edgeD_curr(_edgeD_curr), edgeD_next(_edgeD_next), vertI(_vertI), vertD_curr(_vertD_curr), vertD_next(_vertD_next), offsets(_offsets) {}
    inline bool update(intT s, intT d, intT edgeIdx)
    {
        intT dstIdx = offsets[d] + edgeIdx;
        for (int i = 0; i < NSTATES; i++)
        {
            edgeD_next[dstIdx].belief[i] = 0.0;
            for (int j = 0; j < NSTATES; j++)
            {
                edgeD_next[dstIdx].belief[i] += vertI[d].potential[j] * edgeW[dstIdx].potential[i][j] * vertD_curr[d].product[j];
            }
            vertD_next[d].product[i] = vertD_next[d].product[i] * edgeD_next[dstIdx].belief[i];
        }
        return 1;
    }

    inline void create_cache(cache_t &cache, intT d)
    {
        for (int i =0; i < NSTATES; i++)
        {
            cache.vertD_cache_next.product[i]= vertD_next[d].product[i];
        }
    }
    inline bool update(cache_t &cache, intT s, intT d, intE edgeIdx)
    {
        intT dstIdx = offsets[s] + edgeIdx;
        for(int i = 0; i<NSTATES; i++)
        {
           edgeD_next[dstIdx].belief[i] = 0.0;
           for( int j=0; j< NSTATES; j++)
                edgeD_next[dstIdx].belief[i] += vertI[d].potential[j] * edgeW[dstIdx].potential[i][j] * vertD_curr[d].product[j];
           cache.vertD_cache_next.product[i]= cache.vertD_cache_next.product[i] * edgeD_next[dstIdx].belief[i]; 
        }
        return 1;
    }

    inline void commit_cache(cache_t &cache, intT d)
    {
        for (int i = 0; i < NSTATES; i++)
        {
            vertD_next[d].product[i] = cache.vertD_cache_next.product[i];
        }
    }
    // again, this will be the issue
    // since you can't do atomics..
    inline bool updateAtomic (intT s, intT d, intT edgeIdx)   //atomic Update
    {
        intT dstIdx = offsets[d] + edgeIdx;
        for (int i = 0; i < NSTATES; i++)
        {
            edgeD_next[dstIdx].belief[i] = 0.0;
            for (int j = 0; j < NSTATES; j++)
            {
                //edgeD_next[dstIdx].belief[i] += vertI[d].potential[j] * edgeW[dstIdx].potential[i][j] * vertD_curr[d].product[j];
                writeAdd(&edgeD_next[dstIdx].belief[i],vertI[d].potential[j] * edgeW[dstIdx].potential[i][j] * vertD_curr[d].product[j]);
            }
            writeMult(&(vertD_next[d].product[i]), edgeD_next[dstIdx].belief[i]);
        }
        return 1;
    }

    inline bool cond (intT d)
    {
        return 1;    //does nothing
    }
};

//resets p
struct BP_Vertex_Reset
{
    VertexData *vertD;
    BP_Vertex_Reset(VertexData *_vertD) :
        vertD(_vertD) {}
    inline bool operator () (intT i)
    {
        for (int i = 0; i < NSTATES; i++)
        {
            vertD[i].product[i] = 1.0;
        }
        return 1;
    }
};


template <class GraphType>
void Compute(GraphType &GA, long start)
{
    typedef typename GraphType::vertex_type vertex; // Is determined by GraphType
    const partitioner &part = GA.get_partitioner();
    const int perNode = part.get_num_per_node_partitions();
    //offsets
    intT n = GA.n;
    intT m = GA.m;

    mmap_ptr<intT> Degrees;
    Degrees.Interleave_allocate (n);
    mmap_ptr<intT> Offsets;
    Offsets.Interleave_allocate (n);

    parallel_for(intT j=0; j < n; ++j)
	Degrees[j] = GA.get_partition().V[j].getOutDegree();

    Offsets[0] = 0;
    
    for (intT i = 1; i < n; i++)
        Offsets[i] = Offsets[i-1] + Degrees[i-1];

    intT numEdge = Offsets[n - 1] + Degrees[n - 1];
    //create vertex data
    mmap_ptr<VertexInfo> vertI;
    vertI.part_allocate (part);
    mmap_ptr<VertexData> vertD_curr;
    vertD_curr.part_allocate (part);
    mmap_ptr<VertexData> vertD_next;
    vertD_next.part_allocate (part);

    // create edge data
    mmap_ptr<EdgeWeight> edgeW;
    edgeW.Interleave_allocate (numEdge);
    mmap_ptr<EdgeData> edgeD_curr;
    edgeD_curr.Interleave_allocate (numEdge);
    mmap_ptr<EdgeData> edgeD_next;
    edgeD_next.Interleave_allocate (numEdge);

    partitioned_vertices Frontier = partitioned_vertices::bits(part,n, m);
    int currIter=0;

    // initialise all messages to 1.0
    vertexMap(part,Frontier, BP_Vertex_Reset(vertD_curr));
    // for all iterations.. up to maxIter
    while(1&&currIter<maxIter)
    {
        currIter++;
        // initialise all messages to 1.0
        vertexMap(part,Frontier, BP_Vertex_Reset(vertD_next));
        // for all edges of the graph
        // do the main bit..

        // i have no idea, i don't understand what the belief propagation algorithm is doing
        partitioned_vertices output=edgeMap(GA, Frontier, BP_F<vertex>(edgeW, edgeD_curr, edgeD_next, vertI, vertD_curr, vertD_next, Offsets), m/20);
        output.del();
        swap(edgeD_curr, edgeD_next);
        swap(vertD_curr, vertD_next);
    }
    Frontier.del();
    vertI.del();
    vertD_curr.del();
    vertD_next.del();
    edgeD_curr.del();
    edgeD_next.del();
    edgeW.del();
    Degrees.del();
    Offsets.del();
}


