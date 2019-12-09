#pragma once

#ifndef size_t
typedef unsigned long size_t;
#endif

class ManSegType
{};

class ManSegHead : public ManSegType
{
public:
    ManSegHead(unsigned int& head)
        :head(head)
    {}

    template<typename T>
    ManSegHead& operator=(const T& rhs);
    template<typename T>
    double operator=(const T& t);
    template<typename T>
    double operator+=(const T& rhs);
    template<typename T>
    double operator-=(const T& rhs);
    template<typename T>
    double operator*=(const T& rhs);
    template<typename T>
    double operator/=(const T& rhs);

    operator double() const
    {
        size_t l = head;
        l <<= segmentBits;
        return *reinterpret_cast<double*>(&l);
    }

    unsigned int& head;
    static constexpr int segmentBits = 32;

};

class ManSegPair : public ManSegType
{
public:
    ManSegPair(unsigned int& head, unsigned int& tail)
        :head(head), tail(tail)
    {}

    template<typename T>
    ManSegPair& operator=(const T& rhs);
    template<typename T>
    double operator+=(const T& rhs);
    template<typename T>
    double operator-=(const T& rhs);
    template<typename T>
    double operator*=(const T& rhs);
    template<typename T>
    double operator/=(const T& rhs);

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
    // how to define [] and = operators
    // still giving me trouble figuring out how to do this.

    // if we want to be able to manipulate values via references, the int& head and int& tail need to be replaced
    // by a single size_t& val
    // this *doubles* the memory requirements of heads, which probably defeats the purpose of "less bits being moved when performing arithmetic"
    // and the "less memory used" thing

    // alternatively, we try to get some kinda ManSegType inheritance working correctly, so we can access ManSeg**** functions through pointers
    // but this has the same issue of, when we have ManSegHead/ManSegPair how do we modify the underlying values??

    virtual const ManSegArray& operator=(const ManSegArray& rhs) = 0;

    virtual void set(size_t id, double d) = 0;
    virtual void set(size_t id, float f) = 0;
    virtual void set(size_t id, int i) = 0;
    virtual void set(size_t id, long l) = 0;

    virtual void setPair(size_t id, double d) = 0;
    virtual void setPair(size_t id, float f) = 0;
    virtual void setPair(size_t id, int i) = 0;
    virtual void setPair(size_t id, long li) = 0;

    virtual double read(size_t id) = 0;
    virtual ManSegArray* updoot() = 0;

    unsigned int* heads;
    unsigned int* tails;
};

template<bool useTail = false>
class ManSegBase : public ManSegArray
{
public:
    ManSegBase(size_t length)
    {
        heads = new unsigned int[length];
        tails = new unsigned int[length];
    }

    ManSegBase(unsigned int* heads, unsigned int* tails)
        :heads(heads), tails(tails)
    {}

    ~ManSegBase()
    {
        delete[] heads, tails;
    }

    // experimental: and by that i mean i can't get it working
    virtual const ManSegArray& operator=(const ManSegArray& rhs) override
    {
        unsigned int* h = rhs.heads;
        unsigned int* t = rhs.tails;
        heads = h;
        tails = t;

        return *this;
    }

    // return manseghead/tail in operator[]

    virtual void set(size_t id, double d) override;
    virtual void set(size_t id, float f) override;
    virtual void set(size_t id, int i) override;
    virtual void set(size_t id, long l) override;

    virtual void setPair(size_t id, double d) override;
    virtual void setPair(size_t id, float f) override;
    virtual void setPair(size_t id, int i) override;
    virtual void setPair(size_t id, long li) override;

    virtual double read(size_t id) override;

    virtual ManSegBase<true>* updoot() override;


    unsigned int* heads;
    unsigned int* tails;

private:
    static constexpr int segmentBits = 32;
    static constexpr unsigned int tailMask = ~0; // todo: should getBits() ?
    static constexpr size_t headMask = static_cast<unsigned long>(tailMask) << segmentBits;

};


/**
 * ManSegHead functions
*/
template<>
ManSegHead& ManSegHead::operator=(const ManSegHead& rhs)
{
    unsigned int h = rhs.head;
    head = h;

    return *this;
}

template<>
ManSegHead& ManSegHead::operator=(const ManSegPair& rhs)
{
    unsigned int h = rhs.head;
    head = h;

    return *this;
}

template<typename T>
ManSegHead& ManSegHead::operator=(const T& rhs)
{
    const size_t l = static_cast<const size_t>(rhs);
    head = static_cast<unsigned int>(l >> segmentBits);

    return *this;
}

template<typename T>
double ManSegHead::operator+=(const T& rhs)
{
    double d = *this;
    d += rhs;

    return d;
}

template<typename T>
inline double operator+(ManSegHead lhs, const T& rhs)
{
    unsigned int h = lhs.head;
    double d = static_cast<double>(ManSegHead(h));
    d += rhs;

    return d;
}

template<typename T>
double ManSegHead::operator-=(const T& rhs)
{
    double d = *this;
    d -= rhs;

    return d;
}

template<typename T>
inline double operator-(ManSegHead lhs, const T& rhs)
{
    unsigned int h = lhs.head;
    double d = static_cast<double>(ManSegHead(h));
    d -= rhs;

    return d;
}

template<typename T>
double ManSegHead::operator*=(const T& rhs)
{
    double d = *this;
    d *= rhs;

    return d;
}

template<typename T>
inline double operator*(ManSegHead lhs, const T& rhs)
{
    unsigned int h = lhs.head;
    double d = static_cast<double>(ManSegHead(h));
    d *= rhs;

    return d;
}

template<typename T>
double ManSegHead::operator/=(const T& rhs)
{
    double d = *this;
    d /= rhs;

    return d;
}

template<typename T>
inline double operator/(ManSegHead lhs, const T& rhs)
{
    unsigned int h = lhs.head;
    double d = static_cast<double>(ManSegHead(h));
    d /= rhs;

    return d;
}

/**
 * ManSegPair functions
*/
template<>
ManSegPair& ManSegPair::operator=(const ManSegPair& rhs)
{
    unsigned int h = rhs.head;
    unsigned int t = rhs.tail;
    head = h;
    tail = t;

    return *this;
}

template<>
ManSegPair& ManSegPair::operator=(const ManSegHead& rhs)
{
    unsigned int h = rhs.head;
    unsigned int t = 0U;
    head = h;
    tail = t;

    return *this;
}

template<typename T>
ManSegPair& ManSegPair::operator=(const T& rhs)
{
    const size_t l = static_cast<const size_t>(rhs);
    head = static_cast<unsigned int>(l >> segmentBits);
    tail = static_cast<unsigned int>(l & tailMask);

    return *this;
}

template <typename T>
double ManSegPair::operator+=(const T& rhs)
{
    double d = *this;
    d += rhs;

    return d;
}

template<typename T>
inline double operator+(ManSegPair lhs, const T& rhs)
{
    unsigned int h = lhs.head;
    unsigned int t = lhs.tail;
    double d = static_cast<double>(ManSegPair(h, t));
    d += rhs;

    return d;
}

template<typename T>
double ManSegPair::operator-=(const T& rhs)
{
    double d = *this;
    d -= rhs;

    return d;
}

template<typename T>
inline double operator-(ManSegPair lhs, const T& rhs)
{
    unsigned int h = lhs.head;
    unsigned int t = lhs.tail;
    double d = static_cast<double>(ManSegPair(h, t));
    d -= rhs;

    return d;
}

template<typename T>
double ManSegPair::operator*=(const T& rhs)
{
    double d = *this;
    d *= rhs;

    return d;
}

template<typename T>
inline double operator*(ManSegPair lhs, const T& rhs)
{
    unsigned int h = lhs.head;
    unsigned int t = lhs.tail;
    double d = static_cast<double>(ManSegPair(h, t));
    d *= rhs;

    return d;
}

template<typename T>
double ManSegPair::operator/=(const T& rhs)
{
    double d = *this;
    d /= rhs;

    return d;
}

template<typename T>
inline double operator/(ManSegPair lhs, const T& rhs)
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

/**
 * false specialisation
 * (heads only)
*/
template<>
void ManSegBase<false>::set(size_t id, double d)
{
    const size_t l = *reinterpret_cast<const size_t*>(&d);
    heads[id] = (l >> segmentBits);
}

template<>
void ManSegBase<false>::set(size_t id, float f)
{
    double d = f;
    const size_t l = *reinterpret_cast<const size_t*>(&d);
    heads[id] = (l >> segmentBits);
}

template<>
void ManSegBase<false>::set(size_t id, int i)
{
    double d = i;
    const size_t l = *reinterpret_cast<const size_t*>(&d);
    heads[id] = (l >> segmentBits);
}

template<>
void ManSegBase<false>::set(size_t id, long li)
{
    // double d = li;
    const size_t l = *reinterpret_cast<const size_t*>(&li);
    heads[id] = (l >> segmentBits);
}

template<>
double ManSegBase<false>::read(size_t id)
{
    return static_cast<double>(ManSegHead(heads[id]));
}


/**
 * true specialisation
 * (heads + tails)
*/
template<>
void ManSegBase<true>::set(size_t id, double d)
{
    const size_t l = *reinterpret_cast<const size_t*>(&d);
    heads[id] = l >> segmentBits;
    tails[id] = (l & tailMask);
}

template<>
void ManSegBase<true>::set(size_t id, float f)
{
    double d = f;
    const size_t l = *reinterpret_cast<const size_t*>(&d);
    heads[id] = l >> segmentBits;
    tails[id] = (l & tailMask);
}

template<>
void ManSegBase<true>::set(size_t id, int i)
{
    double d = i;
    const size_t l = *reinterpret_cast<const size_t*>(&d);
    heads[id] = l >> segmentBits;
    tails[id] = (l & tailMask);
}

template<>
void ManSegBase<true>::set(size_t id, long li)
{
    // double d = li;
    const size_t l = *reinterpret_cast<const size_t*>(&li);
    heads[id] = l >> segmentBits;
    tails[id] = (l & tailMask);
}

template<>
double ManSegBase<true>::read(size_t id)
{
    return static_cast<double>(ManSegPair(heads[id], tails[id]));
}


/**
 *  Common functions 
*/
template<bool useTail>
void ManSegBase<useTail>::setPair(size_t id, double d)
{
    const size_t l = *reinterpret_cast<const size_t*>(&d);
    heads[id] = l >> segmentBits;
    tails[id] = (l & tailMask);
}

template<bool useTail>
void ManSegBase<useTail>::setPair(size_t id, float f)
{
    double d = f;
    const size_t l = *reinterpret_cast<const size_t*>(&d);
    heads[id] = l >> segmentBits;
    tails[id] = (l & tailMask);
}

template<bool useTail>
void ManSegBase<useTail>::setPair(size_t id, int i)
{
    double d = i;
    const size_t l = *reinterpret_cast<const size_t*>(&d);
    heads[id] = l >> segmentBits;
    tails[id] = (l & tailMask);
}

template<bool useTail>
void ManSegBase<useTail>::setPair(size_t id, long li)
{
    // double d = li;
    const size_t l = *reinterpret_cast<const size_t*>(&li);
    heads[id] = l >> segmentBits;
    tails[id] = (l & tailMask);
}

template<bool useTail>
ManSegBase<true>* ManSegBase<useTail>::updoot()
{
    return new ManSegBase<true>(heads, tails);
}
