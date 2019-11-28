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

    // todo: template, but will need to be specialised.
    // manseghead/manseghead needs specialisation can stay as it is
    // if rhs = double or float or int, need to convert to manseghead then take head 
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

class ManSegPair
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
    {
    }

    ~ManSegBase()
    {
        delete[] heads, tails;
    }

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

private:
    unsigned int *heads;
    unsigned int *tails;
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

/*
    inline double operator+(ManSegHead lhs, const ManSegHead& rhs)
    {
        unsigned int h = lhs.head;
        ManSegHead result(h);
        result += rhs;

        return static_cast<double>(result);
    }
*/

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
    double d = li;
    const size_t l = *reinterpret_cast<const size_t*>(&d);
    heads[id] = (l >> segmentBits);
}

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
    double d = li;
    const size_t l = *reinterpret_cast<const size_t*>(&d);
    heads[id] = l >> segmentBits;
    tails[id] = (l & tailMask);
}

template<>
double ManSegBase<false>::read(size_t id)
{
    return static_cast<double>(ManSegHead(heads[id]));
}

template<>
ManSegBase<true>* ManSegBase<false>::updoot()
{
    return new ManSegBase<true>(this->heads, this->tails);
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
    double d = li;
    const size_t l = *reinterpret_cast<const size_t*>(&d);
    heads[id] = l >> segmentBits;
    tails[id] = (l & tailMask);
}

template<>
double ManSegBase<true>::read(size_t id)
{
    return static_cast<double>(ManSegPair(heads[id], tails[id]));
}

template<>
ManSegBase<true>* ManSegBase<true>::updoot()
{
    // this doesn't do anything but needs to be here
    // because of the template type 
    return this;
}
