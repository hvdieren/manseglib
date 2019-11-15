#pragma once

/*
	Tester header class used for testing out some ideas with segmentation.
	Based on a single-element design, so is designated _old, as it will not 
	be used going forward.
*/

constexpr unsigned int TWO_SEG_FULL_TAIL_MASK = ~0;
constexpr int TWO_SEGMENT_BITS = 32;
constexpr size_t TWO_SEG_HEAD_MASK = static_cast<size_t>(TWO_SEG_FULL_TAIL_MASK) << TWO_SEGMENT_BITS;
unsigned int tailBits = 0;

// only useful for n > 0 and n < 32
constexpr int maskBits(int n) { return (static_cast<unsigned int>(TWO_SEG_FULL_TAIL_MASK) << (32 - n)); }

class ManSeg {
public:
    ManSeg()
    {
        head = 0;
        tail = 0;
    }

    ManSeg(double& d) // todo: make template ?
    {
        const size_t l = *reinterpret_cast<const size_t*>(&d);
        head = (l & TWO_SEG_HEAD_MASK) >> TWO_SEGMENT_BITS;
        tail = l & TWO_SEG_FULL_TAIL_MASK;
    }

    // todo: read head, then read n bits of tail.
    operator double() { return toDouble(*this); }

    ManSeg& operator+=(const ManSeg& rhs)
    {
        double d1 = toDouble(*this); // use toDouble, not toFullDouble
        double d2 = toDouble(rhs);
        d1 += d2;

        setSegments(d1);

        return *this; // don't return *this
    }

    friend ManSeg operator+(ManSeg lhs, const ManSeg& rhs) // this could probably be better also
    {
        lhs += rhs;
        return lhs;
    }

private:
    constexpr double toFullDouble(const ManSeg& m) { size_t l = m.head; l <<= TWO_SEGMENT_BITS; l |= m.tail; return *reinterpret_cast<double*>(&l); }
    constexpr double toHalfDouble(const ManSeg& m) { size_t l = m.head; l <<= TWO_SEGMENT_BITS; return *reinterpret_cast<double*>(&l); }
    constexpr void setSegments(const double& d) { const size_t l = *reinterpret_cast<const size_t*>(&d); head = (l & TWO_SEG_HEAD_MASK) >> TWO_SEGMENT_BITS; tail = l & TWO_SEG_FULL_TAIL_MASK; }
    double toDouble(const ManSeg& m);
   
    unsigned int head, tail;
};

double ManSeg::toDouble(const ManSeg& m)
{
    size_t l = m.head;
    l <<= TWO_SEGMENT_BITS;

    if(tailBits == TWO_SEGMENT_BITS)
        l |= m.tail;
    else if(tailBits > 0)
        l |= (m.tail & maskBits(tailBits));

    return *reinterpret_cast<double*>(&l);
}
