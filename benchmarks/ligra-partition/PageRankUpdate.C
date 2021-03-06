// This code is part of the project "Ligra: A Lightweight Graph Processing
// Framework for Shared Memory", presented at Principles and Practice of
// Parallel Programming, 2013.
// Copyright (c) 2013 Julian Shun and Guy Blelloch
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights (to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject to
// the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
// LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
// OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#include "ligra-numa.h"
#include "math.h"
int MaxIter=100;
template <class vertex>
struct PR_F
{
    double* p_curr, *p_next;
    double damping;
    vertex* V;
    static const bool use_cache = true;

    struct cache_t
    {
        double p_next;
    };
    PR_F(double* _p_curr, double* _p_next, double _damping, vertex* _V) :
        p_curr(_p_curr), p_next(_p_next), damping(_damping), V(_V) {}
    inline bool update(intT s, intT d)  //update function applies PageRank equation
    {
        p_next[d] += damping*(p_curr[s]/V[s].getOutDegree());
        return 1;
    }
    inline bool updateAtomic (intT s, intT d)   //atomic Update
    {
        writeAdd(&p_next[d], damping*(p_curr[s]/V[s].getOutDegree()));
        return 1;
    }

    inline void create_cache(cache_t &cache, intT d)
    {
        cache.p_next = p_next[d];
    }
    inline bool update(cache_t &cache, intT s)
    {
        cache.p_next += damping*(p_curr[s]/V[s].getOutDegree());
        return 1;
    }

    inline void commit_cache(cache_t &cache, intT d)
    {
        // Cache is used only in sequential mode
        p_next[d] = cache.p_next;
    }

    inline bool cond (intT d)
    {
        return cond_true(d);    //does nothing
    }
};

//resets p
struct PR_Vertex_Reset
{
    double* p_curr;
    PR_Vertex_Reset(double* _p_curr) :
        p_curr(_p_curr) {}
    inline bool operator () (intT i)
    {
        p_curr[i] = 0.0;
        return 1;
    }
};

double seqsum(double* a, intT n)
{
    double d = 0.;
    double err = 0.;
    for( intT i=0; i < n; ++i ) 
    {
	    //does d += a[i]; but with high accuracy
        double tmp = d;
        double y = a[i] + err;
        d = tmp + y;
        err = tmp - d;
        err += y;
    }
    return d;
}

double sumArray(const partitioner &part, double* a, intT n)
{
    double d = 0.;
	double err = 0.;
	double tmp, y;
    int p = part.get_num_partitions();
    double *psum = new double [p];
    // get sum of each partition
    map_partition( k, part, {
        intT s = part.start_of(k);
        intT e = part.start_of(k+1);
        psum[k] = seqsum( &a[s], e-s);
    });
    for( int i=0; i < p; ++i ) 
    {
        // does d += psum[i]; but with high accuracy
        tmp = d;
        y = psum[i] + err;
        d = tmp + y;
        err = tmp - d;
        err += y;
    }
    delete [] psum;
    return d;
}

double seqnormdiff(double* a, double* b, intT n)
{
    double d = 0.;
    double err = 0.;
    for( intT i=0; i < n; ++i ) 
    {
        //does d += fabs(a[i]- b[i]); but with high accuracy
        double tmp = d;
        double y = fabs(a[i] - b[i]) + err;
        d = tmp + y;
        err = tmp - d;
        err += y;
    }
    return d;
}

double normDiff(const partitioner &part, double* a, double* b, intT n)
{
    double d = 0.;
    int p = part.get_num_partitions();
    double *psum = new double [p];
    double err = 0.;
    double tmp, y;
    /* calculate sum per partition */
    map_partition( k, part, {
        intT s = part.start_of(k);
        intT e = part.start_of(k+1);
        psum[k] = seqnormdiff( &a[s], &b[s], e-s);
    } );

    for( int i=0; i < p; ++i ) 
    {
        // does d += psum[i]; but with high accuracy
        tmp = d;
        y = psum[i] + err;
        d = tmp + y;
        err = tmp - d;
        err += y;
    }
    delete [] psum;
    return d;
}

template <class GraphType>
void Compute(GraphType &GA, long start)
{
    typedef typename GraphType::vertex_type vertex; // Is determined by GraphTyp
    const partitioner &part = GA.get_partitioner();
    graph<vertex> & WG = GA.get_partition();
    const int perNode = part.get_num_per_node_partitions();
    const intT partElem = part.get_num_elements();
    intT n = GA.n;
    intT m = GA.m;
    const double damping = 0.85;
    const double epsilon = 0.0000001;
    //Data Array
    //x and y to do special allocation
    //frontier also need special node allocation
    //blocksize equal to the szie of each partitioned
    double one_over_n = 1/(double)n;

    mmap_ptr<double> p_curr;
    p_curr.part_allocate (part);
    mmap_ptr<double> p_next;
    p_next.part_allocate (part);

    double delta = 2.0;
    loop(j, part, perNode, p_curr[j] = one_over_n);
    loop(j, part, perNode, p_next[j] = 0);

    cerr << setprecision(16);
    
	int count=0;
    partitioned_vertices Frontier = partitioned_vertices::bits(part,n, m);
    while(count<MaxIter)
    {
		++count;

        // p_next[d] += damping * (p_curr[s]/V[s].getOutDegree())
        partitioned_vertices output = edgeMap(GA, Frontier, PR_F<vertex>(p_curr,p_next,damping,WG.V),m/20);

		// find value to scale PR vals by to make vector add to 1
		double scaleAdditive = (1 - sumArray(part, p_next, n))*one_over_n;
		{
			loop(j, part, perNode, p_next[j] += scaleAdditive);
		}

		// delta = abs(p_curr - p_next)
		delta = normDiff(part, p_curr, p_next, n);
        if(delta < epsilon)
        {
			cerr << count << ": delta = " << delta << "\n";
            cerr << "successfully converged in " << count << " iterations\n";
            break;
        }
        cerr << count << ": delta = " << delta << "  norm. sum = " << sumArray(part, p_next, n) << "\n";


        // reset p_curr and swap vertices
		{
			loop(j, part, perNode, p_curr[j] = 0);
		}
        swap(p_curr, p_next);
        // manage fronteir stuff
        Frontier.del();
        Frontier = output;
    }

    Frontier.del();
    p_curr.del();
    p_next.del();
}
