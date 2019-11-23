#pragma once

#ifndef size_t
typedef unsigned long size_t;
#endif

class ManSegHead
{
public:
    ManSegHead(unsigned int& head)
        :head(head)
    {}

    ManSegHead& operator=(const ManSegHead& rhs);
    double operator=(const double& d); // todo: might be useful to do these.
    double operator+=(const ManSegHead& rhs);
    double operator-=(const ManSegHead& rhs);
    double operator*=(const ManSegHead& rhs);
    double operator/=(const ManSegHead& rhs);
    
    operator double() const
    {
        size_t l = head;
        l <<= segmentBits;
        return *reinterpret_cast<double*>(&l);
    }

    unsigned int& head;
    static constexpr int segmentBits = 32;

};

class ManSegPair
{
public:
    ManSegPair(unsigned int& head, unsigned int& tail)
        :head(head), tail(tail)
    {}

    ManSegPair& operator=(const ManSegPair& rhs);
    double operator+=(const ManSegPair& rhs);
    double operator-=(const ManSegPair& rhs);
    double operator*=(const ManSegPair& rhs);
    double operator/=(const ManSegPair& rhs);

    operator double() const
    {
        size_t l = head;
        l <<= segmentBits;
        l |= tail;
        return *reinterpret_cast<double*>(&l);
    }

    unsigned int& head;
    unsigned int& tail;
    static constexpr int segmentBits = 32;
    static constexpr unsigned int tailMask = ~0;

};

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

    void setHead(size_t id, double d);
    double readHead(size_t id);
    void setPair(size_t id, double d);
    double readPair(size_t id);
    
private:
    static constexpr int segmentBits = 32;
    static constexpr unsigned int tailMask = ~0; // todo: should getBits() ?
    static constexpr size_t headMask = static_cast<unsigned long>(tailMask) << segmentBits;
    unsigned int* heads;
    unsigned int* tails;

};

/**
 * ManSegHead functions
*/
ManSegHead& ManSegHead::operator=(const ManSegHead& rhs)
{
    unsigned int h = rhs.head;
    head = h;

    return *this;
}

double ManSegHead::operator+=(const ManSegHead& rhs)
{
    double d = *this;
    d += rhs;

    return d;
}

/*
    inline double operator+(ManSegHead lhs, const ManSegHead& rhs)
    {
        unsigned int h = lhs.head;
        ManSegHead result(h);
        result += rhs;

        return static_cast<double>(result);
    }
*/

inline double operator+(ManSegHead lhs, const ManSegHead& rhs)
{
    unsigned int h = lhs.head;
    double d = static_cast<double>(ManSegHead(h));
    d += rhs;

    return d;
}

double ManSegHead::operator-=(const ManSegHead& rhs)
{
    double d = *this;
    d -= rhs;

    return d;
}

inline double operator-(ManSegHead lhs, const ManSegHead& rhs)
{
    unsigned int h = lhs.head;
    double d = static_cast<double>(ManSegHead(h));
    d -= rhs;

    return d;
}

double ManSegHead::operator*=(const ManSegHead& rhs)
{
    double d = *this;
    d *= rhs;

    return d;
}

inline double operator*(ManSegHead lhs, const ManSegHead& rhs)
{
    unsigned int h = lhs.head;
    double d = static_cast<double>(ManSegHead(h));
    d *= rhs;

    return d;
}

double ManSegHead::operator/=(const ManSegHead& rhs)
{
    double d = *this;
    d /= rhs;

    return d;
}

inline double operator/(ManSegHead lhs, const ManSegHead& rhs)
{
    unsigned int h = lhs.head;
    double d = static_cast<double>(ManSegHead(h));
    d /= rhs;

    return d;
}

/**
 * ManSegPair functions
*/
ManSegPair& ManSegPair::operator=(const ManSegPair& rhs)
{
    unsigned int h = rhs.head;
    unsigned int t = rhs.tail;
    head = h;
    tail = t;

    return *this;
}

double ManSegPair::operator+=(const ManSegPair& rhs)
{
    double d = *this;
    d += rhs;

    return d;
}

inline double operator+(ManSegPair lhs, const ManSegPair& rhs)
{
    unsigned int h = lhs.head;
    unsigned int t = lhs.tail;
    double d = static_cast<double>(ManSegPair(h, t));
    d += rhs;

    return d;
}

double ManSegPair::operator-=(const ManSegPair& rhs)
{
    double d = *this;
    d -= rhs;

    return d;
}

inline double operator-(ManSegPair lhs, const ManSegPair& rhs)
{
    unsigned int h = lhs.head;
    unsigned int t = lhs.tail;
    double d = static_cast<double>(ManSegPair(h, t));
    d -= rhs;

    return d;
}

double ManSegPair::operator*=(const ManSegPair& rhs)
{
    double d = *this;
    d *= rhs;

    return d;
}

inline double operator*(ManSegPair lhs, const ManSegPair& rhs)
{
    unsigned int h = lhs.head;
    unsigned int t = lhs.tail;
    double d = static_cast<double>(ManSegPair(h, t));
    d *= rhs;

    return d;
}

double ManSegPair::operator/=(const ManSegPair& rhs)
{
    double d = *this;
    d /= rhs;

    return d;
}

inline double operator/(ManSegPair lhs, const ManSegPair& rhs)
{
    unsigned int h = lhs.head;
    unsigned int t = lhs.tail;
    double d = static_cast<double>(ManSegPair(h, t));
    d /= rhs;

    return d;
}

/**
 * ManSegArray functions
*/
void ManSegArray::setHead(size_t id, double d)
{
    const size_t l = *reinterpret_cast<const size_t*>(&d);
    heads[id] = l >> segmentBits;
}

double ManSegArray::readHead(size_t id)
{
    return static_cast<double>(ManSegHead(heads[id]));
}

void ManSegArray::setPair(size_t id, double d)
{
    const size_t l = *reinterpret_cast<const size_t*>(&d);
    
    heads[id] = l >> segmentBits;
    tails[id] = (l & tailMask);
}

double ManSegArray::readPair(size_t id)
{
    return static_cast<double>(ManSegPair(heads[id], tails[id]));
}
