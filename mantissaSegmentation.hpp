#pragma once

class ManSegHead
{
public:
    ManSegHead(unsigned int& head)
        :head(head)
    {}

    ManSegHead& operator+=(const ManSegHead& rhs)
    {}

    friend ManSegHead operator+(ManSegHead lhs, const ManSegHead& rhs)
    {
        lhs += rhs;
        return lhs;
    }

    operator double()
    {
        return *reinterpret_cast<double*>((static_cast<unsigned long>(head) << segmentBits));
    }

private:
    unsigned int& head;
    static constexpr int segmentBits = 32;

};

class ManSegPair
{
public:
    ManSegPair(unsigned int head, unsigned int tail)
        :head(head), tail(tail)
    {}

    ManSegPair& operator+=(const ManSegPair& rhs)
    {

    }

    friend ManSegPair operator+(ManSegPair lhs, const ManSegPair& rhs)
    {
        lhs += rhs;
        return lhs;
    }

    /*ManSegPair operator+=(const ManSegHead& rhs)
    {}

    friend ManSegPair operator+(ManSegPair lhs, const ManSegHead& rhs)
    {
        lhs += rhs;
        return lhs;
    }*/

    // segmentation (fault)
    operator double()
    {
        unsigned long l = head;
        l <<= segmentBits;
        l |= tail;
        return *reinterpret_cast<double*>(l);
    }

private:
    unsigned int& head, tail;
    static constexpr int segmentBits = 32;

};

// for this, you need to overload the arithmetic
// operators that you would use
class ManSegArray
{
public:
    ManSegArray(size_t length)
    {
        heads = new unsigned int[length];
        tails = new unsigned int[length];
        len = length;
    }

    ~ManSegArray()
    {
        if(heads) delete[] heads;
        if(tails) delete[] tails;
    }

    // uuhh, how do we assign this [i] of things
    // the arrays are ints, so that's what you're addressing

    // void setHead(int id, unsigned int& val)
    // {
    //     heads[id] = val;
    // }

    void setHead(size_t id, unsigned int val)
    {
        heads[id] = val;
    }

    // void setPair(int id, unsigned int& head, unsigned int& tail)
    // {
    //     heads[id] = head;
    //     tails[id] = tail;
    // }

    void setPair(size_t id, unsigned int head, unsigned int tail)
    {
        heads[id] = head;
        tails[id] = tail;
    }

    // void setPair(int id, double& d)
    // {
    //     const size_t l = *reinterpret_cast<const size_t*>(&d);
    //     heads[id] = ((l & headMask) >> segmentBits);
    //     tails[id] = (l & tailMask);
    // }

    void setPair(size_t id, double d)
    {
        const size_t l = *reinterpret_cast<const size_t*>(&d);
        heads[id] = ((l & headMask) >> segmentBits);
        tails[id] = (l & tailMask);
    }

    double operator[](size_t id)
    {
        return ManSegPair(heads[id], tails[id]);
    }

    friend std::ostream& operator<<(std::ostream& os, const ManSegArray& m);

    // operator double()
    // {
    //     return ManSegPair();
    // }

// private:
    ManSegHead readHead(size_t id)
    {
        return ManSegHead(heads[id]);
    }

    ManSegPair readPair(size_t id)
    {
        return ManSegPair(heads[id], tails[id]);
    }

    static constexpr int segmentBits = 32;
    static constexpr unsigned int tailMask = ~0; // todo: should getBits()
    static constexpr unsigned long headMask = static_cast<unsigned long>(tailMask) << segmentBits;
    unsigned int* heads;
    unsigned int* tails;
    size_t len;

};

std::ostream& operator<<(std::ostream& os, const ManSegArray& m)
{
    for(size_t i = 0; i < (sizeof(m.heads)/sizeof(int)); ++i)
    {
        os << ManSegPair(m.heads[i], m.tails[i]);
    }
}
