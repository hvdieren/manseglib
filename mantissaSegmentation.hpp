#pragma once

namespace ManSeg
{
    // internal representation of double, so it can be easily manipulated
    typedef __UINT_FAST64_TYPE__ doublerep;

    /*
        Class representing the "head" segment of a double.
        It is defined for ease of manipulating the values in an object of type TwoSegmentArray<false>.
        It contains the sign bit, full 11 bit exponent, and 20 bits of mantissa, for a total
        of 32 bits.
        This gives less precision than IEEE-754 standard float.
        It has a maximum precision of roughly 1e-6
    */
    class Head
    {
    public:
        Head(unsigned int& head)
            :head(head)
        {}

        template<typename T>
        Head& operator=(const T& rhs);
        template<typename T>
        Head& operator=(const T&& other) noexcept;
        
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
            return *reinterpret_cast<double*>(&l);
        }

        unsigned int& head;
        static constexpr int segmentBits = 32;

    };

    /*
        Class representing the "pair" of segments of a double.
        It is defined for ease of manipulating the values in an object of type TwoSegmentArray<true>.
        Each segment (which consists of a head and tail) is 32 bits.
        The head contains the sign bit, full 11 bit exponent, and 20 bits of mantissa, 
        for a total of 32 bits.
        The tail segment contains the remaining 32 bits of mantissa.
    */
    class Pair
    {
    public:
        Pair(unsigned int& head, unsigned int& tail)
            :head(head), tail(tail)
        {}


        Pair& operator=(const Pair& rhs);
        template<typename T>
        Pair& operator=(const T& rhs);
        template<typename T>
        Pair& operator=(const T&& other) noexcept;

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

    /*
        Empty templated class for performing operations on double values, conceptually split into
        two 32 bit segments - "head" and "tail".
    */
    template<bool useTail>
    class TwoSegmentArray; // we specialise this below.

    /*
        Specialisation of TwoSegmentArray.
        Operations performed on values in the array only modify the "head" segment
        (i.e. the first 32 bits of a double value) unless specified otherwise.
    */
    template<>
    class TwoSegmentArray<false>
    {
    public:
        TwoSegmentArray(size_t length)
        {
            heads = new unsigned int[length];
            tails = new unsigned int[length];
        }

        TwoSegmentArray(unsigned int* heads, unsigned int* tails)
            :heads(heads), tails(tails)
        {}

        ~TwoSegmentArray()
        {
            heads = nullptr;
            tails = nullptr;
        }

        template<typename T>
        void set(size_t id, T t);
        template<typename T>
        void setPair(size_t id, T t);

        double read(size_t id);

        Head operator[](size_t id);

        /*
            Increases the precision of the operations performed by returning an
            object of type TwoSegmentArray<true> that has pointers to the values
            in the existing array.
        */
        TwoSegmentArray<true> increasePrecision();

        /*
            Deletes the values of the dynamic arrays used to store values in
            the array.
            NOTE: this should only be called by one object with references to
            the same set of values (such as object created using increasePrecision).
        */
        void del();

    private:
        unsigned int* heads;
        unsigned int* tails;

        static constexpr int segmentBits = 32;
        static constexpr unsigned int tailMask = ~0;
        static constexpr doublerep headMask = static_cast<doublerep>(tailMask) << segmentBits;

    };

    /*
        Specialisation of TwoSegmentArray.
        Operations performed on values in the array are exactly the same as standard IEEE-754 doubles.
    */
    template<>
    class TwoSegmentArray<true>
    {
    public:
        TwoSegmentArray(size_t length)
        {
            len = length;
            heads = new unsigned int[length];
            tails = new unsigned int[length];
        }

        TwoSegmentArray(unsigned int* heads, unsigned int* tails)
            :heads(heads), tails(tails)
        {}

        ~TwoSegmentArray()
        {
            heads = nullptr;
            tails = nullptr;
        }

        Pair operator[](size_t id);

        template<typename T>
        void set(size_t id, T t);
        template<typename T>
        void setPair(size_t id, T t);

        double read(size_t id);

        /*
            Returns *this (as we do not have any precision increase to do)
        */
        TwoSegmentArray<true> increasePrecision();

        /*
            Deletes the values of the dynamic arrays used to store values in
            the array.
            NOTE: this should only be called by one object with references to
            the same set of values (such as object created using increasePrecision).
        */
        void del();

    private:
        unsigned int* heads;
        unsigned int* tails;
        size_t len;

        static constexpr int segmentBits = 32;
        static constexpr unsigned int tailMask = ~0;
        static constexpr doublerep headMask = static_cast<doublerep>(tailMask) << segmentBits;

    };


    /**
     * Head functions
    */
    template<>
    Head& Head::operator=(const Head& other)
    {
        unsigned int h = other.head;
        head = h;
        return *this;
    }

    template<>
    Head& Head::operator=(const Head&& other) noexcept
    {
        unsigned int h = other.head;
        head = h;
        return *this;
    }

    template<>
    Head& Head::operator=(const Pair& other)
    {
        unsigned int h = other.head;
        head = h;
        return *this;
    }

    template<>
    Head& Head::operator=(const Pair&& other) noexcept
    {
        unsigned int h = other.head;
        head = h;
        return *this;
    }

    template<typename T>
    Head& Head::operator=(const T& rhs)
    {
        double d = rhs;
        const doublerep l = *reinterpret_cast<const doublerep*>(&d);
        head = (l >> segmentBits);
        return *this;
    }

    template<typename T>
    Head& Head::operator=(const T&& other) noexcept
    {
        double d = other;
        const doublerep l = *reinterpret_cast<const doublerep*>(&d);
        unsigned int h = (l >> segmentBits);

        head = h;
        return *this;
    }

    template<typename T>
    double Head::operator+=(const T& rhs)
    {
        double t = *this;

        t += rhs;
        *this = t;

        return t;
    }

    template<typename T>
    inline double operator+(Head lhs, const T& rhs)
    {
        double d = lhs;
        d += rhs;
        return d;
    }

    template<typename T>
    double Head::operator-=(const T& rhs)
    {
        double t = *this;
        t -= rhs;
        *this = t;
        return t;
    }

    template<typename T>
    inline double operator-(Head lhs, const T& rhs)
    {
        double d = lhs;
        d -= rhs;
        return d;
    }

    template<typename T>
    double Head::operator*=(const T& rhs)
    {
        double t = *this;
        t *= rhs;
        *this = t;
        return t;
    }

    template<typename T>
    inline double operator*(Head lhs, const T& rhs)
    {
        double d = lhs;
        d *= rhs;
        return d;
    }

    template<typename T>
    double Head::operator/=(const T& rhs)
    {
        double t = *this;
        t /= rhs;
        *this = t;
        return t;
    }

    template<typename T>
    inline double operator/(Head lhs, const T& rhs)
    {
        double d = lhs;
        d /= rhs;
        return d;
    }

    /**
     * ManSegPair functions
    */

    template<>
    Pair& Pair::operator=(const Head& other)
    {
        unsigned int h = other.head;
        head = h;
        return *this;
    }

    template<>
    Pair& Pair::operator=(const Head&& other) noexcept
    {
        unsigned int h = other.head;
        head = h;
        return *this;
    }

    Pair& Pair::operator=(const Pair& other)
    {
        unsigned int h = other.head;
        unsigned int t = other.tail;
        head = h;
        tail = t;
        return *this;
    }

    template<>
    Pair& Pair::operator=(const Pair&& other) noexcept
    {
        unsigned int h = other.head;
        unsigned int t = other.tail;
        head = h;
        tail = t;
        return *this;
    }

    template<typename T>
    Pair& Pair::operator=(const T&& other) noexcept
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
    Pair& Pair::operator=(const T& rhs)
    {
        const doublerep l = *reinterpret_cast<const doublerep*>(&rhs);
        head = static_cast<unsigned int>(l >> segmentBits);
        tail = static_cast<unsigned int>(l & tailMask);
        return *this;
    }

    template <typename T>
    double Pair::operator+=(const T& rhs)
    {
        double t = *this;
        t += rhs;
        *this = t;
        return t;
    }

    template<typename T>
    inline double operator+(Pair lhs, const T& rhs)
    {
        double d = lhs;
        d += rhs;
        return d;
    }

    template<typename T>
    double Pair::operator-=(const T& rhs)
    {
        double t = *this;
        t -= rhs;
        *this = t;
        return t;
    }

    template<typename T>
    inline double operator-(Pair lhs, const T& rhs)
    {
        double d = lhs;
        d -= rhs;
        return d;
    }

    template<typename T>
    double Pair::operator*=(const T& rhs)
    {
        double t = *this;
        t *= rhs;
        *this = t;
        return t;
    }

    template<typename T>
    inline double operator*(Pair lhs, const T& rhs)
    {
        double d = lhs;
        d *= rhs;
        return d;
    }

    template<typename T>
    double Pair::operator/=(const T& rhs)
    {
        double t = *this;
        t /= rhs;
        *this = t;
        return t;
    }

    template<typename T>
    inline double operator/(Pair lhs, const T& rhs)
    {
        double d = lhs;
        d /= rhs;
        return d;
    }


    /**
     * TwoSegmentArray functions
    */

    /**
     * false specialisation
     * (heads only)
    */
    template<typename T>
    void TwoSegmentArray<false>::set(size_t id, T t)
    {
        double d = t;
        const doublerep l = *reinterpret_cast<const doublerep*>(&d);
        heads[id] = (l >> segmentBits);
    }

    template<typename T>
    void TwoSegmentArray<false>::setPair(size_t id, T t)
    {
        double d = t;
        const doublerep l = *reinterpret_cast<const doublerep*>(&d);
        heads[id] = l >> segmentBits;
        tails[id] = (l & tailMask);
    }

    double TwoSegmentArray<false>::read(size_t id)
    {
        return static_cast<double>(Head(heads[id]));
    }

    Head TwoSegmentArray<false>::operator[](size_t id)
    {
        return Head(heads[id]);
    }

    TwoSegmentArray<true> TwoSegmentArray<false>::increasePrecision()
    {
        return TwoSegmentArray<true>(heads, tails);
    }

    void TwoSegmentArray<false>::del()
    {
        if(heads) delete[] heads;
        if(tails) delete[] tails;
    }


    /**
     * true specialisation
     * (heads + tails)
    */
    template<typename T>
    void TwoSegmentArray<true>::set(size_t id, T t)
    {
        double d = t;
        const doublerep l = *reinterpret_cast<const doublerep*>(&d);
        heads[id] = l >> segmentBits;
        tails[id] = (l & tailMask);
    }

    template<typename T>
    void TwoSegmentArray<true>::setPair(size_t id, T t)
    {
        double d = t;
        const doublerep l = *reinterpret_cast<const doublerep*>(&d);
        heads[id] = l >> segmentBits;
        tails[id] = (l & tailMask);
    }

    double TwoSegmentArray<true>::read(size_t id)
    {
        return static_cast<double>(Pair(heads[id], tails[id]));
    }

    Pair TwoSegmentArray<true>::operator[](size_t id)
    {
        return Pair(heads[id], tails[id]);
    }

    TwoSegmentArray<true> TwoSegmentArray<true>::increasePrecision()
    {
        return *this;
    }

    void TwoSegmentArray<true>::del()
    {
        if(heads) delete[] heads;
        if(tails) delete[] tails;
    }
}