#pragma once

#if DOS
typedef unsigned long long doublerep;
#else
typedef unsigned long doublerep;
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
        // doublerep l = head;
        // l <<= segmentBits;
        // return *reinterpret_cast<double*>(&l);
        doublerep l = 0UL;
        l |= head;
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
        doublerep l = head;
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
        len = length;
        heads = new unsigned int[length];
        tails = new unsigned int[length];
    }

    // TODO: need to update this once you figure out how to fix the ManSegBase<true> one
    ManSegBase(unsigned int* heads, unsigned int* tails)
        :heads(heads), tails(tails)
    {
    }

    ~ManSegBase()
    {
        heads = nullptr;
        tails = nullptr;
    }

    void set(size_t id, double d);
    void set(size_t id, float f);
    void set(size_t id, int i);
    void set(size_t id, long l);

    void setPair(size_t id, double d);
    void setPair(size_t id, float f);
    void setPair(size_t id, int i);
    void setPair(size_t id, long li);

    double read(size_t id);

    // return manseghead/tail in operator[]
    ManSegHead operator[](size_t id);
    const ManSegBase<false>& operator=(const ManSegBase<false>& rhs);
    const ManSegBase<false>& operator=(double* rhs);

    void toNewDoubleArray(double* d); // todo: remove

    ManSegBase<true> updoot();

    void freeMemory();

    // friend class ManSegBase<true>;

private:
    unsigned int* heads;
    unsigned int* tails;
    size_t len; // TODO: remove

    static constexpr int segmentBits = 32;
    static constexpr unsigned int tailMask = ~0;
    static constexpr doublerep headMask = static_cast<doublerep>(tailMask) << segmentBits;

};

template<>
class ManSegBase<true>
{
public:
    ManSegBase(size_t length)
    {
        len = length;
        heads = new unsigned int[length];
        tails = new unsigned int[length];
    }

    ManSegBase(unsigned int* heads, unsigned int* tails, size_t length)
        :heads(heads), tails(tails), len(length)
    {}

    ~ManSegBase()
    {
        heads = nullptr;
        tails = nullptr;
    }

    // return manseghead/tail in operator[]

    ManSegPair operator[](size_t id);
    const ManSegBase<true>& operator=(const ManSegBase<true>& rhs);
    const ManSegBase<true>& operator=(double* rhs);

    void toNewDoubleArray(double* d);

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

    void freeMemory();

private:
    unsigned int* heads;
    unsigned int* tails;
    size_t len;

    static constexpr int segmentBits = 32;
    static constexpr unsigned int tailMask = ~0;
    static constexpr doublerep headMask = static_cast<doublerep>(tailMask) << segmentBits;

};


/**
 * ManSegHead functions
*/
template<>
ManSegHead& ManSegHead::operator=(const ManSegHead& other)
{
    unsigned int h = other.head;
    head = h;
    return *this;
}

template<>
ManSegHead& ManSegHead::operator=(const ManSegHead&& other) noexcept
{
    unsigned int h = other.head;
    head = h;
    return *this;
}

template<>
ManSegHead& ManSegHead::operator=(const ManSegPair& other)
{
    unsigned int h = other.head;
    head = h;
    return *this;
}

template<>
ManSegHead& ManSegHead::operator=(const ManSegPair&& other) noexcept
{
    unsigned int h = other.head;
    head = h;
    return *this;
}

template<typename T>
ManSegHead& ManSegHead::operator=(const T& rhs)
{
    double d = rhs;
    const doublerep l = *reinterpret_cast<const doublerep*>(&d);
    head = (l >> segmentBits);
    return *this;
}

template<typename T>
ManSegHead& ManSegHead::operator=(const T&& other) noexcept
{
    double d = other;
    const doublerep l = *reinterpret_cast<const doublerep*>(&d);
    unsigned int h = (l >> segmentBits);

    head = h;
    return *this;
}

template<typename T>
double ManSegHead::operator+=(const T& rhs)
{
    double t = *this;

    t += rhs;
    *this = t;

    return t;
}

template<typename T>
inline double operator+(ManSegHead lhs, const T& rhs)
{
    // unsigned int h = lhs.head;
    // double d = static_cast<double>(ManSegHead(h));
    // d += rhs;

    double d = lhs;

    d += rhs;

    return d;
}

template<typename T>
double ManSegHead::operator-=(const T& rhs)
{
    double t = *this;

    t -= rhs;
    *this = t;

    return t;
}

template<typename T>
inline double operator-(ManSegHead lhs, const T& rhs)
{
    // unsigned int h = lhs.head;
    // double d = static_cast<double>(ManSegHead(h));
    // d -= rhs;

    double d = lhs;

    d -= rhs;

    return d;
}

template<typename T>
double ManSegHead::operator*=(const T& rhs)
{
    double t = *this;

    t *= rhs;
    *this = t;

    return t;
}

template<typename T>
inline double operator*(ManSegHead lhs, const T& rhs)
{
    // unsigned int h = lhs.head;
    // double d = ManSegHead(h);
    // d *= rhs;

    double d = lhs;

    d *= rhs;

    return d;
}

template<typename T>
double ManSegHead::operator/=(const T& rhs)
{
    double t = *this;

    t /= rhs;
    *this = t;

    return t;
}

template<typename T>
inline double operator/(ManSegHead lhs, const T& rhs)
{
    // unsigned int h = lhs.head;
    // double d = ManSegHead(h);
    // d /= rhs;

    double d = lhs;

    d /= rhs;

    return d;
}

/**
 * ManSegPair functions
*/

template<>
ManSegPair& ManSegPair::operator=(const ManSegHead& other)
{
    unsigned int h = other.head;
    head = h;
    // tail = 0U;
    return *this;
}

template<>
ManSegPair& ManSegPair::operator=(const ManSegHead&& other) noexcept
{
    unsigned int h = other.head;
    head = h;
    // tail = 0U;
    return *this;
}

ManSegPair& ManSegPair::operator=(const ManSegPair& other)
{
    unsigned int h = other.head;
    unsigned int t = other.tail;
    head = h;
    tail = t;
    return *this;
}

template<>
ManSegPair& ManSegPair::operator=(const ManSegPair&& other) noexcept
{
    unsigned int h = other.head;
    unsigned int t = other.tail;
    head = h;
    tail = t;
    return *this;
}

template<typename T>
ManSegPair& ManSegPair::operator=(const T&& other) noexcept
{
    double d = other;
    const doublerep l = *reinterpret_cast<const doublerep*>(&d);
    unsigned int h = (l >> segmentBits);
    unsigned int t = (l & tailMask);

    head = h;
    tail = t;
    return *this;
}

template<typename T>
ManSegPair& ManSegPair::operator=(const T& rhs)
{
    const doublerep l = *reinterpret_cast<const doublerep*>(&rhs);
    head = static_cast<unsigned int>(l >> segmentBits);
    tail = static_cast<unsigned int>(l & tailMask);

    return *this;
}

template <typename T>
double ManSegPair::operator+=(const T& rhs)
{
    double t = *this;

    t += rhs;
    *this = t;

    return t;
}

template<typename T>
inline double operator+(ManSegPair lhs, const T& rhs)
{
    // unsigned int h = lhs.head;
    // unsigned int t = lhs.tail;
    // double d = static_cast<double>(ManSegPair(h, t));
    double d = lhs;
    
    d += rhs;

    return d;
}

template<typename T>
double ManSegPair::operator-=(const T& rhs)
{
    double t = *this;

    t -= rhs;
    *this = t;

    return t;
}

template<typename T>
inline double operator-(ManSegPair lhs, const T& rhs)
{
    double d = lhs;
    
    d -= rhs;

    return d;
}

template<typename T>
double ManSegPair::operator*=(const T& rhs)
{
    double t = *this;

    t *= rhs;
    *this = t;

    return t;
}

template<typename T>
inline double operator*(ManSegPair lhs, const T& rhs)
{
    double d = lhs;
    
    d *= rhs;

    return d;
}

template<typename T>
double ManSegPair::operator/=(const T& rhs)
{
    double t = *this;

    t /= rhs;
    *this = t;

    return t;
}

template<typename T>
inline double operator/(ManSegPair lhs, const T& rhs)
{
    double d = lhs;

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
    const doublerep l = *reinterpret_cast<const doublerep*>(&d);
    heads[id] = (l >> segmentBits);
}

void ManSegBase<false>::set(size_t id, float f)
{
    double d = f;
    const doublerep l = *reinterpret_cast<const doublerep*>(&d);
    heads[id] = (l >> segmentBits);
}

void ManSegBase<false>::set(size_t id, int i)
{
    double d = i;
    const doublerep l = *reinterpret_cast<const doublerep*>(&d);
    heads[id] = (l >> segmentBits);
}

void ManSegBase<false>::set(size_t id, long li)
{
    // double d = li;
    const doublerep l = *reinterpret_cast<const doublerep*>(&li);
    heads[id] = (l >> segmentBits);
}

void ManSegBase<false>::setPair(size_t id, double d)
{
    const doublerep l = *reinterpret_cast<const doublerep*>(&d);
    heads[id] = l >> segmentBits;
    tails[id] = (l & tailMask);
}

void ManSegBase<false>::setPair(size_t id, float f)
{
    double d = f;
    const doublerep l = *reinterpret_cast<const doublerep*>(&d);
    heads[id] = l >> segmentBits;
    tails[id] = (l & tailMask);
}

void ManSegBase<false>::setPair(size_t id, int i)
{
    double d = i;
    const doublerep l = *reinterpret_cast<const doublerep*>(&d);
    heads[id] = l >> segmentBits;
    tails[id] = (l & tailMask);
}

void ManSegBase<false>::setPair(size_t id, long li)
{
    // double d = li;
    const doublerep l = *reinterpret_cast<const doublerep*>(&li);
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

const ManSegBase<false>& ManSegBase<false>::operator=(const ManSegBase<false>& rhs)
{
    unsigned int* h = rhs.heads;
    heads = h;

    return *this;
}

const ManSegBase<false>& ManSegBase<false>::operator=(double* rhs)
{
    for(size_t i = 0; i < len; ++i)
        operator[](i) = rhs[i];

    return *this;
}

void ManSegBase<false>::toNewDoubleArray(double* d)
{
    for(size_t i = 0; i < len; ++i)
        d[i] = ManSegHead(heads[i]);
}

ManSegBase<true> ManSegBase<false>::updoot()
{
    return ManSegBase<true>(heads, tails, len);
}

void ManSegBase<false>::freeMemory()
{
    if(heads) delete[] heads;
    if(tails) delete[] tails;
}


/**
 * true specialisation
 * (heads + tails)
*/
void ManSegBase<true>::set(size_t id, double d)
{
    const doublerep l = *reinterpret_cast<const doublerep*>(&d);
    heads[id] = l >> segmentBits;
    tails[id] = (l & tailMask);
}

void ManSegBase<true>::set(size_t id, float f)
{
    double d = f;
    const doublerep l = *reinterpret_cast<const doublerep*>(&d);
    heads[id] = l >> segmentBits;
    tails[id] = (l & tailMask);
}

void ManSegBase<true>::set(size_t id, int i)
{
    double d = i;
    const doublerep l = *reinterpret_cast<const doublerep*>(&d);
    heads[id] = l >> segmentBits;
    tails[id] = (l & tailMask);
}

void ManSegBase<true>::set(size_t id, long li)
{
    // double d = li;
    const doublerep l = *reinterpret_cast<const doublerep*>(&li);
    heads[id] = l >> segmentBits;
    tails[id] = (l & tailMask);
}

void ManSegBase<true>::setPair(size_t id, double d)
{
    const doublerep l = *reinterpret_cast<const doublerep*>(&d);
    heads[id] = l >> segmentBits;
    tails[id] = (l & tailMask);
}

void ManSegBase<true>::setPair(size_t id, float f)
{
    double d = f;
    const doublerep l = *reinterpret_cast<const doublerep*>(&d);
    heads[id] = l >> segmentBits;
    tails[id] = (l & tailMask);
}

void ManSegBase<true>::setPair(size_t id, int i)
{
    double d = i;
    const doublerep l = *reinterpret_cast<const doublerep*>(&d);
    heads[id] = l >> segmentBits;
    tails[id] = (l & tailMask);
}

void ManSegBase<true>::setPair(size_t id, long li)
{
    // double d = li;
    const doublerep l = *reinterpret_cast<const doublerep*>(&li);
    heads[id] = l >> segmentBits;
    tails[id] = (l & tailMask);
}

double ManSegBase<true>::read(size_t id)
{
    return static_cast<double>(ManSegPair(heads[id], tails[id]));
}

const ManSegBase<true>& ManSegBase<true>::operator=(const ManSegBase<true>& rhs)
{
    unsigned int* h = rhs.heads;
    unsigned int* t = rhs.tails;
    heads = h;
    tails = t;

    return *this;
}

ManSegPair ManSegBase<true>::operator[](size_t id)
{
    return ManSegPair(heads[id], tails[id]);
}

const ManSegBase<true>& ManSegBase<true>::operator=(double* rhs)
{
    for(size_t i = 0; i < len; ++i)
        operator[](i) = rhs[i];

    return *this;
}

void ManSegBase<true>::toNewDoubleArray(double* d)
{
    for(size_t i = 0; i < len; ++i)
        d[i] = ManSegPair(heads[i], tails[i]);
}

ManSegBase<true> ManSegBase<true>::updoot()
{
    return *this;
}

void ManSegBase<true>::freeMemory()
{
    if(heads) delete[] heads;
    if(tails) delete[] tails;
}


// Interim attempt
// this time, we are using doublerep as underlying type and just
// AND-ing whenever we need to access heads.

// turns out, this has a large performance impact on heads specifically
// this makes accessing and setting heads *way* worse than pairs, which are faster now
// (which should not be the case)

// since heads are also required to use doublerep, this removes the efficiency afforded by
// only reading half data from memory. which sort of defeats the purpose of the library..

// #pragma once

// #if DOS
// typedef unsigned long long doublerep;
// #else
// typedef unsigned long doublerep;
// #endif

// class ManSegHead
// {
// public:
//     ManSegHead(doublerep& head)
//         :head(head)
//     {}

//     template<typename T>
//     ManSegHead& operator=(const T& rhs);
//     template<typename T>
//     ManSegHead& operator=(const T&& other) noexcept;
    
//     // don't know why this exists, might be for some chaining
//     // function i don't know. function not defined anyhow
//     // template<typename T>
//     // double operator=(const T& t);

//     template<typename T>
//     double operator+=(const T& rhs);
//     template<typename T>
//     double operator-=(const T& rhs);
//     template<typename T>
//     double operator*=(const T& rhs);
//     template<typename T>
//     double operator/=(const T& rhs);

//     operator double() const
//     {
//         doublerep l = (head & headMask);
//         return *reinterpret_cast<double*>(&l);
//     }

//     doublerep& head;
//     static constexpr int segmentBits = 32;
//     static constexpr unsigned int tailMask = ~0;
//     static constexpr doublerep headMask = static_cast<doublerep>(tailMask) << segmentBits;

// };

// class ManSegPair
// {
// public:
//     ManSegPair(doublerep& pair)
//         :pair(pair)
//     {}


//     ManSegPair& operator=(const ManSegPair& rhs);
//     template<typename T>
//     ManSegPair& operator=(const T& rhs);
//     template<typename T>
//     ManSegPair& operator=(const T&& other) noexcept;

//     template<typename T>
//     double operator+=(const T& rhs);
//     template<typename T>
//     double operator-=(const T& rhs);
//     template<typename T>
//     double operator*=(const T& rhs);
//     template<typename T>
//     double operator/=(const T& rhs);

//     operator double() const
//     {
//         // doublerep l = head;
//         // l <<= segmentBits;
//         // l |= tail;
//         doublerep l = pair;
//         return *reinterpret_cast<double*>(&l);
//     }

//     // unsigned int& head;
//     // unsigned int& tail;
//     doublerep& pair;
//     static constexpr int segmentBits = 32;
//     static constexpr unsigned int tailMask = ~0;
//     static constexpr doublerep headMask = static_cast<doublerep>(tailMask) << segmentBits;

// };

// template<bool useTail>
// class ManSegBase; // we specialise this below.

// template<>
// class ManSegBase<false>
// {
// public:
//     ManSegBase(size_t length)
//     {
//         len = length;
//         // heads = new unsigned int[length];
//         // tails = new unsigned int[length];
//         vals = new doublerep[length];
//     }

//     // TODO: need to update this once you figure out how to fix the ManSegBase<true> one
//     // ManSegBase(unsigned int* heads, unsigned int* tails)
//     //     :heads(heads), tails(tails)
//     // {
//     // }
//     ManSegBase(doublerep* vals)
//         :vals(vals)
//     {}

//     ~ManSegBase()
//     {
//         // heads = nullptr;
//         // tails = nullptr;
//         vals = nullptr;
//     }

//     void set(size_t id, double d);
//     void set(size_t id, float f);
//     void set(size_t id, int i);
//     void set(size_t id, long l);

//     void setPair(size_t id, double d);
//     void setPair(size_t id, float f);
//     void setPair(size_t id, int i);
//     void setPair(size_t id, long li);

//     double read(size_t id);

//     // return manseghead/tail in operator[]
//     ManSegHead operator[](size_t id);
//     const ManSegBase<false>& operator=(const ManSegBase<false>& rhs);
//     const ManSegBase<false>& operator=(double* rhs);

//     void toNewDoubleArray(double* d);

//     ManSegBase<true> updoot();

//     void freeMemory();

//     // friend class ManSegBase<true>;

// private:
//     // unsigned int* heads;
//     // unsigned int* tails;

//     doublerep* vals;
//     size_t len;

//     static constexpr int segmentBits = 32;
//     static constexpr unsigned int tailMask = ~0;
//     static constexpr doublerep headMask = static_cast<doublerep>(tailMask) << segmentBits;

// };

// template<>
// class ManSegBase<true>
// {
// public:
//     ManSegBase(size_t length)
//     {
//         len = length;
//         // heads = new unsigned int[length];
//         // tails = new unsigned int[length];
//         vals = new doublerep[length];
//     }

//     // ManSegBase(unsigned int* heads, unsigned int* tails, size_t length)
//     //     :heads(heads), tails(tails), len(length)
//     // {}
//     ManSegBase(doublerep* vals, size_t len)
//         :vals(vals), len(len)
//     {}

//     ~ManSegBase()
//     {
//         // heads = nullptr;
//         // tails = nullptr;
//         vals = nullptr;
//     }

//     // return manseghead/tail in operator[]

//     ManSegPair operator[](size_t id);
//     const ManSegBase<true>& operator=(const ManSegBase<true>& rhs);
//     const ManSegBase<true>& operator=(double* rhs);

//     void toNewDoubleArray(double* d);

//     void set(size_t id, double d);
//     void set(size_t id, float f);
//     void set(size_t id, int i);
//     void set(size_t id, long l);

//     void setPair(size_t id, double d);
//     void setPair(size_t id, float f);
//     void setPair(size_t id, int i);
//     void setPair(size_t id, long li);

//     double read(size_t id);

//     ManSegBase<true> updoot();

//     void freeMemory();

// private:
//     // unsigned int* heads;
//     // unsigned int* tails;
//     doublerep* vals;
//     size_t len;

//     static constexpr int segmentBits = 32;
//     static constexpr unsigned int tailMask = ~0;
//     static constexpr doublerep headMask = static_cast<doublerep>(tailMask) << segmentBits;

// };


// /**
//  * ManSegHead functions
// */
// template<>
// ManSegHead& ManSegHead::operator=(const ManSegHead& other)
// {
//     doublerep h = other.head;
//     head = h;
//     return *this;
// }

// template<>
// ManSegHead& ManSegHead::operator=(const ManSegHead&& other) noexcept
// {
//     doublerep h = other.head;
//     head = h;
//     return *this;
// }

// template<>
// ManSegHead& ManSegHead::operator=(const ManSegPair& other)
// {
//     doublerep h = (other.pair & headMask);
//     head = h;
//     return *this;
// }

// template<>
// ManSegHead& ManSegHead::operator=(const ManSegPair&& other) noexcept
// {
//     doublerep h = (other.pair & headMask);
//     head = h;
//     return *this;
// }

// template<typename T>
// ManSegHead& ManSegHead::operator=(const T& rhs)
// {
//     double d = rhs;
//     const doublerep l = *reinterpret_cast<const doublerep*>(&d);
//     head = (l & headMask);
//     return *this;
// }

// template<typename T>
// ManSegHead& ManSegHead::operator=(const T&& other) noexcept
// {
//     double d = other;
//     const doublerep l = *reinterpret_cast<const doublerep*>(&d);
    
//     doublerep h = (l & headMask);

//     head = h;
//     return *this;
// }

// template<typename T>
// double ManSegHead::operator+=(const T& rhs)
// {
//     double t = *this;

//     t += rhs;
//     *this = t;

//     return t;
// }

// template<typename T>
// inline double operator+(ManSegHead lhs, const T& rhs)
// {
//     // unsigned int h = lhs.head;
//     // double d = static_cast<double>(ManSegHead(h));
//     // d += rhs;

//     double d = lhs;

//     d += rhs;

//     return d;
// }

// template<typename T>
// double ManSegHead::operator-=(const T& rhs)
// {
//     double t = *this;

//     t -= rhs;
//     *this = t;

//     return t;
// }

// template<typename T>
// inline double operator-(ManSegHead lhs, const T& rhs)
// {
//     // unsigned int h = lhs.head;
//     // double d = static_cast<double>(ManSegHead(h));
//     // d -= rhs;

//     double d = lhs;

//     d -= rhs;

//     return d;
// }

// template<typename T>
// double ManSegHead::operator*=(const T& rhs)
// {
//     double t = *this;

//     t *= rhs;
//     *this = t;

//     return t;
// }

// template<typename T>
// inline double operator*(ManSegHead lhs, const T& rhs)
// {
//     // unsigned int h = lhs.head;
//     // double d = ManSegHead(h);
//     // d *= rhs;

//     double d = lhs;

//     d *= rhs;

//     return d;
// }

// template<typename T>
// double ManSegHead::operator/=(const T& rhs)
// {
//     double t = *this;

//     t /= rhs;
//     *this = t;

//     return t;
// }

// template<typename T>
// inline double operator/(ManSegHead lhs, const T& rhs)
// {
//     // unsigned int h = lhs.head;
//     // double d = ManSegHead(h);
//     // d /= rhs;

//     double d = lhs;

//     d /= rhs;

//     return d;
// }

// /**
//  * ManSegPair functions
// */

// template<>
// ManSegPair& ManSegPair::operator=(const ManSegHead& other)
// {
//     doublerep h = other.head;
//     pair = h;
//     // tail = 0U;
//     return *this;
// }

// template<>
// ManSegPair& ManSegPair::operator=(const ManSegHead&& other) noexcept
// {
//     doublerep h = other.head;
//     pair = h;
//     // tail = 0U;
//     return *this;
// }

// ManSegPair& ManSegPair::operator=(const ManSegPair& other)
// {
//     doublerep p = other.pair;
//     pair = p;
//     return *this;
// }

// template<>
// ManSegPair& ManSegPair::operator=(const ManSegPair&& other) noexcept
// {
//     doublerep p = other.pair;
//     pair = p;
//     return *this;
// }

// template<typename T>
// ManSegPair& ManSegPair::operator=(const T&& other) noexcept
// {
//     double d = other;
//     const doublerep l = *reinterpret_cast<const doublerep*>(&d);
//     pair = l;

//     return *this;
// }

// template<typename T>
// ManSegPair& ManSegPair::operator=(const T& rhs)
// {
//     double d = rhs;
//     const doublerep l = *reinterpret_cast<const doublerep*>(&d);
//     pair = l;

//     return *this;
// }

// template <typename T>
// double ManSegPair::operator+=(const T& rhs)
// {
//     double t = *this;

//     t += rhs;
//     *this = t;

//     return t;
// }

// template<typename T>
// inline double operator+(ManSegPair lhs, const T& rhs)
// {
//     // unsigned int h = lhs.head;
//     // unsigned int t = lhs.tail;
//     // double d = static_cast<double>(ManSegPair(h, t));
//     double d = lhs;
    
//     d += rhs;

//     return d;
// }

// template<typename T>
// double ManSegPair::operator-=(const T& rhs)
// {
//     double t = *this;

//     t -= rhs;
//     *this = t;

//     return t;
// }

// template<typename T>
// inline double operator-(ManSegPair lhs, const T& rhs)
// {
//     double d = lhs;
    
//     d -= rhs;

//     return d;
// }

// template<typename T>
// double ManSegPair::operator*=(const T& rhs)
// {
//     double t = *this;

//     t *= rhs;
//     *this = t;

//     return t;
// }

// template<typename T>
// inline double operator*(ManSegPair lhs, const T& rhs)
// {
//     double d = lhs;
    
//     d *= rhs;

//     return d;
// }

// template<typename T>
// double ManSegPair::operator/=(const T& rhs)
// {
//     double t = *this;

//     t /= rhs;
//     *this = t;

//     return t;
// }

// template<typename T>
// inline double operator/(ManSegPair lhs, const T& rhs)
// {
//     double d = lhs;

//     d /= rhs;

//     return d;
// }


// /**
//  * ManSegArray functions
// */

// /**
//  * false specialisation
//  * (heads only)
// */
// void ManSegBase<false>::set(size_t id, double d)
// {
//     const doublerep l = *reinterpret_cast<const doublerep*>(&d);
//     // vals[id] = (l & headMask);
//     vals[id] = l;
// }

// void ManSegBase<false>::set(size_t id, float f)
// {
//     double d = f;
//     const doublerep l = *reinterpret_cast<const doublerep*>(&d);
//     vals[id] = l;
// }

// void ManSegBase<false>::set(size_t id, int i)
// {
//     double d = i;
//     const doublerep l = *reinterpret_cast<const doublerep*>(&d);
//     vals[id] = l;
// }

// void ManSegBase<false>::set(size_t id, long li)
// {
//     // double d = li;
//     const doublerep l = *reinterpret_cast<const doublerep*>(&li);
//     vals[id] = l;
// }

// void ManSegBase<false>::setPair(size_t id, double d)
// {
//     const doublerep l = *reinterpret_cast<const doublerep*>(&d);
//     vals[id] = l;
// }

// void ManSegBase<false>::setPair(size_t id, float f)
// {
//     double d = f;
//     const doublerep l = *reinterpret_cast<const doublerep*>(&d);
//     vals[id] = l;
// }

// void ManSegBase<false>::setPair(size_t id, int i)
// {
//     double d = i;
//     const doublerep l = *reinterpret_cast<const doublerep*>(&d);
//     vals[id] = l;
// }

// void ManSegBase<false>::setPair(size_t id, long li)
// {
//     // double d = li;
//     const doublerep l = *reinterpret_cast<const doublerep*>(&li);
//     vals[id] = l;
// }

// double ManSegBase<false>::read(size_t id)
// {
//     return static_cast<double>(ManSegHead(vals[id]));
// }

// ManSegHead ManSegBase<false>::operator[](size_t id)
// {
//     return ManSegHead(vals[id]);
// }

// const ManSegBase<false>& ManSegBase<false>::operator=(const ManSegBase<false>& rhs)
// {
//     for(int i = 0; i < len; ++i)
//         vals[i] = rhs.vals[i];

//     return *this;
// }

// const ManSegBase<false>& ManSegBase<false>::operator=(double* rhs)
// {
//     for(size_t i = 0; i < len; ++i)
//         operator[](i) = rhs[i];

//     return *this;
// }

// void ManSegBase<false>::toNewDoubleArray(double* d)
// {
//     for(size_t i = 0; i < len; ++i)
//         d[i] = ManSegHead(vals[i]);
// }

// ManSegBase<true> ManSegBase<false>::updoot()
// {
//     return ManSegBase<true>(vals, len);
// }

// void ManSegBase<false>::freeMemory()
// {
//     if(vals) delete[] vals;
// }


// /**
//  * true specialisation
//  * (heads + tails)
// */
// void ManSegBase<true>::set(size_t id, double d)
// {
//     const doublerep l = *reinterpret_cast<const doublerep*>(&d);
//     vals[id] = l;
// }

// void ManSegBase<true>::set(size_t id, float f)
// {
//     double d = f;
//     const doublerep l = *reinterpret_cast<const doublerep*>(&d);
//     vals[id] = l;
// }

// void ManSegBase<true>::set(size_t id, int i)
// {
//     double d = i;
//     const doublerep l = *reinterpret_cast<const doublerep*>(&d);
//     vals[id] = l;
// }

// void ManSegBase<true>::set(size_t id, long li)
// {
//     // double d = li;
//     const doublerep l = *reinterpret_cast<const doublerep*>(&li);
//     vals[id] = l;
// }

// void ManSegBase<true>::setPair(size_t id, double d)
// {
//     const doublerep l = *reinterpret_cast<const doublerep*>(&d);
//     vals[id] = l;
// }

// void ManSegBase<true>::setPair(size_t id, float f)
// {
//     double d = f;
//     const doublerep l = *reinterpret_cast<const doublerep*>(&d);
//     vals[id] = l;
// }

// void ManSegBase<true>::setPair(size_t id, int i)
// {
//     double d = i;
//     const doublerep l = *reinterpret_cast<const doublerep*>(&d);
//     vals[id] = l;
// }

// void ManSegBase<true>::setPair(size_t id, long li)
// {
//     // double d = li;
//     const doublerep l = *reinterpret_cast<const doublerep*>(&li);
//     vals[id] = l;
// }

// double ManSegBase<true>::read(size_t id)
// {
//     return static_cast<double>(ManSegPair(vals[id]));
// }

// const ManSegBase<true>& ManSegBase<true>::operator=(const ManSegBase<true>& rhs)
// {
//     for(int i = 0; i < len; ++i)
//         vals[i] = rhs.vals[i];

//     return *this;
// }

// ManSegPair ManSegBase<true>::operator[](size_t id)
// {
//     return ManSegPair(vals[id]);
// }

// const ManSegBase<true>& ManSegBase<true>::operator=(double* rhs)
// {
//     for(size_t i = 0; i < len; ++i)
//         operator[](i) = rhs[i];

//     return *this;
// }

// void ManSegBase<true>::toNewDoubleArray(double* d)
// {
//     for(size_t i = 0; i < len; ++i)
//         d[i] = ManSegPair(vals[i]);
// }

// ManSegBase<true> ManSegBase<true>::updoot()
// {
//     return *this;
// }

// void ManSegBase<true>::freeMemory()
// {
//     if(vals) delete[] vals;
// }
