#pragma once
#include <iostream>

#ifndef size_t
typedef unsigned long size_t;
#endif

class ManSegHead
{
public:
    ManSegHead(unsigned int& head)
        :head(head)
    {}

    ManSegHead& operator+=(const ManSegHead& rhs)
    {
        double lhs = *this;
        double rs = rhs;
        lhs += rs;

        const size_t l = *reinterpret_cast<const size_t*>(&lhs);
        head = (l >> segmentBits);

        return *this; // what do you return other than this..?
    }

    friend ManSegHead operator+(ManSegHead lhs, const ManSegHead& rhs)
    {
        lhs += rhs;
        return lhs;
    }

    // operator double()
    // {
    //     size_t l = head;
    //     l <<= segmentBits;
    //     return *reinterpret_cast<double*>(&l);
    // }

    operator double() const
    {
        size_t l = head;
        l <<= segmentBits;
        return *reinterpret_cast<double*>(&l);
    }

private:
    unsigned int& head;
    static constexpr int segmentBits = 32;

};

class ManSegPair
{
public:
    ManSegPair(unsigned int& head, unsigned int& tail)
        :head(head), tail(tail)
    {}

    // todo: other operators

    // todo: head+tail addition.
    ManSegPair& operator+=(const ManSegPair& rhs)
    {
        double lhs = *this;
        double rs = rhs;

        lhs += rs;

        const size_t l = *reinterpret_cast<const size_t*>(&lhs);
        tail = (l & tailMask);
        head = (l >> segmentBits);
    }

    friend ManSegPair operator+(ManSegPair lhs, const ManSegPair& rhs)
    {
        lhs += rhs;
        return lhs;
    }

    // do we include these?
    /*ManSegPair operator+=(const ManSegHead& rhs)
    {}

    friend ManSegPair operator+(ManSegPair lhs, const ManSegHead& rhs)
    {
        lhs += rhs;
        return lhs;
    }*/

    operator double() const
    {
        size_t l = head;
        l <<= segmentBits;
        l |= tail;
        return *reinterpret_cast<double*>(&l);
    }

private:
    unsigned int& head;
    unsigned int& tail;
    static constexpr int segmentBits = 32;
    static constexpr unsigned int tailMask = ~0;

};

// todo: for arithmetic operators on elements of this
// should use template args to work out if seghead or segpair
class ManSegArray
{
public:
    ManSegArray(size_t length)
    {
        heads = new unsigned int[length];
        tails = new unsigned int[length];
    }

    ~ManSegArray()
    {
        if(heads) delete[] heads;
        if(tails) delete[] tails;
    }

    void setHead(size_t id, unsigned int val)
    {
        heads[id] = val;
    }

    void setPair(size_t id, unsigned int head, unsigned int tail)
    {
        heads[id] = head;
        tails[id] = tail;
    }

    void setPair(size_t id, double d)
    {
        const size_t l = *reinterpret_cast<const size_t*>(&d);
        heads[id] = l >> segmentBits;
        tails[id] = l & tailMask;
    }

    double toDouble(unsigned int& head, unsigned int& tail)
    {
       size_t l = head;
        l <<= segmentBits;
        l |= tail;

        return *reinterpret_cast<double*>(&l);
    }

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
    static constexpr size_t headMask = static_cast<unsigned long>(tailMask) << segmentBits;
    unsigned int* heads;
    unsigned int* tails;

};


