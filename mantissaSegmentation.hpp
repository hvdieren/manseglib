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

    template<typename T>
    ManSegHead& operator=(const T& rhs);
    template<typename T>
    ManSegHead& operator=(const T&& other) noexcept;
    
    // don't know why this exists, might be for some chaining
    // function i don't know. function not defined anyhow
    // template<typename T>
    // double operator=(const T& t);

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
    ManSegPair& operator=(const T&& other) noexcept;

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

template<bool useTail>
class ManSegBase; // we specialise this below.

template<>
class ManSegBase<false>
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
        copiedFrom = true;
    }

    ~ManSegBase()
    {
        if(!copiedFrom)
        {
            if(heads) delete[] heads;
            if(tails) delete[] tails;
        }
    }

    const ManSegBase<false>& operator=(const ManSegBase<false>& rhs)
    {
        unsigned int* h = rhs.heads;
        unsigned int* t = rhs.tails;
        heads = h;
        tails = t;

        return *this;
    }

    // return manseghead/tail in operator[]
    ManSegHead operator[](size_t id);

    void set(size_t id, double d);
    void set(size_t id, float f);
    void set(size_t id, int i);
    void set(size_t id, long l);

    void setPair(size_t id, double d);
    void setPair(size_t id, float f);
    void setPair(size_t id, int i);
    void setPair(size_t id, long li);

    double read(size_t id);

    ManSegBase<true> updoot();

private:
    unsigned int* heads;
    unsigned int* tails;
    bool copiedFrom = false; // added this to fix dangling pointer issue. one extra bit per array is fine, right?

    static constexpr int segmentBits = 32;
    static constexpr unsigned int tailMask = ~0;
    static constexpr size_t headMask = static_cast<unsigned long>(tailMask) << segmentBits;

};

template<>
class ManSegBase<true>
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
        copiedFrom = true;
    }

    ~ManSegBase()
    {
        if(!copiedFrom)
        {
            if(heads) delete[] heads;
            if(tails) delete[] tails;
        }
    }

    const ManSegBase<true>& operator=(const ManSegBase<true>& rhs)
    {
        unsigned int* h = rhs.heads;
        unsigned int* t = rhs.tails;
        heads = h;
        tails = t;

        return *this;
    }

    // return manseghead/tail in operator[]

    ManSegPair operator[](size_t id);

    void set(size_t id, double d);
    void set(size_t id, float f);
    void set(size_t id, int i);
    void set(size_t id, long l);

    void setPair(size_t id, double d);
    void setPair(size_t id, float f);
    void setPair(size_t id, int i);
    void setPair(size_t id, long li);

    double read(size_t id);

    ManSegBase<true> updoot();

private:
    unsigned int* heads;
    unsigned int* tails;
    bool copiedFrom = false;

    static constexpr int segmentBits = 32;
    static constexpr unsigned int tailMask = ~0;
    static constexpr size_t headMask = static_cast<unsigned long>(tailMask) << segmentBits;

};


/**
 * ManSegHead functions
*/
template<>
ManSegHead& ManSegHead::operator=(const ManSegHead& other)
{
    if(this != &other)
    {
        unsigned int h = other.head;
        head = h;
    }
    return *this;
}

template<>
ManSegHead& ManSegHead::operator=(const ManSegPair& other)
{
    if(&head != &other.head)
    {
        unsigned int h = other.head;
        head = h;
    }
    return *this;
}

template<>
ManSegHead& ManSegHead::operator=(const ManSegHead&& other) noexcept
{
    if(this != &other)
    {
        unsigned int h = other.head;
        head = h;
    }
    return *this;
}

template<>
ManSegHead& ManSegHead::operator=(const ManSegPair&& other) noexcept
{
    if(&head != &other.head)
    {
        unsigned int h = other.head;
        head = h;
    }
    return *this;
}

template<typename T>
ManSegHead& ManSegHead::operator=(const T& rhs)
{
    const size_t l = *reinterpret_cast<const size_t*>(&rhs);
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
ManSegPair& ManSegPair::operator=(const ManSegPair& other)
{
    if(this != &other)
    {
        unsigned int h = other.head;
        unsigned int t = other.tail;
        head = h;
        tail = t;
    }
    return *this;
}

template<>
ManSegPair& ManSegPair::operator=(const ManSegHead& other)
{
    if(&head != &other.head)
    {
        unsigned int h = other.head;
        head = h;
        tail = 0U;
    }
    return *this;
}

template<>
ManSegPair& ManSegPair::operator=(const ManSegPair&& other) noexcept
{
    if(this != &other)
    {
        unsigned int h = other.head;
        unsigned int t = other.tail;
        head = h;
        tail = t;
    }
    return *this;
}

template<>
ManSegPair& ManSegPair::operator=(const ManSegHead&& other) noexcept
{
    if(&head != &other.head)
    {
        unsigned int h = other.head;
        head = h;
        tail = 0U;
    }
    return *this;
}

template<typename T>
ManSegPair& ManSegPair::operator=(const T& rhs)
{
    const size_t l = *reinterpret_cast<const size_t*>(&rhs);
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
void ManSegBase<false>::set(size_t id, double d)
{
    const size_t l = *reinterpret_cast<const size_t*>(&d);
    heads[id] = (l >> segmentBits);
}

void ManSegBase<false>::set(size_t id, float f)
{
    double d = f;
    const size_t l = *reinterpret_cast<const size_t*>(&d);
    heads[id] = (l >> segmentBits);
}

void ManSegBase<false>::set(size_t id, int i)
{
    double d = i;
    const size_t l = *reinterpret_cast<const size_t*>(&d);
    heads[id] = (l >> segmentBits);
}

void ManSegBase<false>::set(size_t id, long li)
{
    // double d = li;
    const size_t l = *reinterpret_cast<const size_t*>(&li);
    heads[id] = (l >> segmentBits);
}

void ManSegBase<false>::setPair(size_t id, double d)
{
    const size_t l = *reinterpret_cast<const size_t*>(&d);
    heads[id] = l >> segmentBits;
    tails[id] = (l & tailMask);
}

void ManSegBase<false>::setPair(size_t id, float f)
{
    double d = f;
    const size_t l = *reinterpret_cast<const size_t*>(&d);
    heads[id] = l >> segmentBits;
    tails[id] = (l & tailMask);
}

void ManSegBase<false>::setPair(size_t id, int i)
{
    double d = i;
    const size_t l = *reinterpret_cast<const size_t*>(&d);
    heads[id] = l >> segmentBits;
    tails[id] = (l & tailMask);
}

void ManSegBase<false>::setPair(size_t id, long li)
{
    // double d = li;
    const size_t l = *reinterpret_cast<const size_t*>(&li);
    heads[id] = l >> segmentBits;
    tails[id] = (l & tailMask);
}

double ManSegBase<false>::read(size_t id)
{
    return static_cast<double>(ManSegHead(heads[id]));
}

ManSegHead ManSegBase<false>::operator[](size_t id)
{
    return ManSegHead(heads[id]);
}

ManSegBase<true> ManSegBase<false>::updoot()
{
    return ManSegBase<true>(heads, tails);
}


/**
 * true specialisation
 * (heads + tails)
*/
void ManSegBase<true>::set(size_t id, double d)
{
    const size_t l = *reinterpret_cast<const size_t*>(&d);
    heads[id] = l >> segmentBits;
    tails[id] = (l & tailMask);
}

void ManSegBase<true>::set(size_t id, float f)
{
    double d = f;
    const size_t l = *reinterpret_cast<const size_t*>(&d);
    heads[id] = l >> segmentBits;
    tails[id] = (l & tailMask);
}

void ManSegBase<true>::set(size_t id, int i)
{
    double d = i;
    const size_t l = *reinterpret_cast<const size_t*>(&d);
    heads[id] = l >> segmentBits;
    tails[id] = (l & tailMask);
}

void ManSegBase<true>::set(size_t id, long li)
{
    // double d = li;
    const size_t l = *reinterpret_cast<const size_t*>(&li);
    heads[id] = l >> segmentBits;
    tails[id] = (l & tailMask);
}

void ManSegBase<true>::setPair(size_t id, double d)
{
    const size_t l = *reinterpret_cast<const size_t*>(&d);
    heads[id] = l >> segmentBits;
    tails[id] = (l & tailMask);
}

void ManSegBase<true>::setPair(size_t id, float f)
{
    double d = f;
    const size_t l = *reinterpret_cast<const size_t*>(&d);
    heads[id] = l >> segmentBits;
    tails[id] = (l & tailMask);
}

void ManSegBase<true>::setPair(size_t id, int i)
{
    double d = i;
    const size_t l = *reinterpret_cast<const size_t*>(&d);
    heads[id] = l >> segmentBits;
    tails[id] = (l & tailMask);
}

void ManSegBase<true>::setPair(size_t id, long li)
{
    // double d = li;
    const size_t l = *reinterpret_cast<const size_t*>(&li);
    heads[id] = l >> segmentBits;
    tails[id] = (l & tailMask);
}

double ManSegBase<true>::read(size_t id)
{
    return static_cast<double>(ManSegPair(heads[id], tails[id]));
}

ManSegPair ManSegBase<true>::operator[](size_t id)
{
    return ManSegPair(heads[id], tails[id]);
}

ManSegBase<true> ManSegBase<true>::updoot()
{
    return *this;
}
