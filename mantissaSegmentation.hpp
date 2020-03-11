#ifndef __MANSEG_LIB_H__
#define __MANSEG_LIB_H__

#include <stdint.h>

namespace ManSeg
{
    // internal representation of double, so it can be easily manipulated
    typedef uint_fast64_t doublerep;

    /* Highest achievable precision with a single segment - i.e. TwoSegArray<false> */
    constexpr double MaxSingleSegmentPrecision = 1e-5;
    /* decimal precision: num_mantissa_bits*log10(2) = 6.02... */
    constexpr double AdaptivePrecisionBound = 5e-5;

    /*
        Class representing the "head" segment of a double.
        It is defined for ease of manipulating the values in an object of type TwoSegArray<false>.
        It contains the sign bit, full 11 bit exponent, and 20 bits of mantissa, for a total
        of 32 bits.
        This gives less precision than IEEE-754 standard float.
        It has a maximum precision of roughly 1e-5, with a recommended precision bound of 5e-5.
    */
    class Head
    {
    public:
        Head(uint32_t& head)
            :head(head)
        {}

        /* 
            These equals operators must either be combined into a generic template function and
            therefore could result in extra conversions when copying data.
            If we want to keep the specialisation, then we need to either move the definition into
            a separate .cpp file and compile a .o file for the library, or we can declare them
            inline, as that allows multiple redefinitions.
            Considering it is only for a few functions, (the TwoSegArray class specialisations must be
            kept in this file) it seems like a reasonable compromise.
        */
        template<typename T>
        inline Head& operator=(const T& rhs);
        template<typename T>
        inline Head& operator=(const T&& other) noexcept;
        
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

        uint32_t& head;
        static constexpr int segmentBits = 32;
    };

    /*
        Class representing the "pair" of segments of a double.
        It is defined for ease of manipulating the values in an object of type TwoSegArray<true>.
        Each segment (which consists of a head and tail) is 32 bits.
        The head contains the sign bit, full 11 bit exponent, and 20 bits of mantissa, 
        for a total of 32 bits.
        The tail segment contains the remaining 32 bits of mantissa.
    */
    class Pair
    {
    public:
        Pair(uint32_t& head, uint32_t& tail)
            :head(head), tail(tail)
        {}


        template<typename T>
        inline Pair& operator=(const T& rhs);
        template<typename T>
        inline Pair& operator=(const T&& other) noexcept;

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

        uint32_t& head;
        uint32_t& tail;
        static constexpr int segmentBits = 32;
        static constexpr uint32_t tailMask = ~0;

    };

    /*
        Array type class for performing operations on double values, conceptually split into two 32 bit segments - "head" and "tail".
        The user is required to manually clean up memory after allocation using the del() function, as the deconstructor does not handle this. This is due to the
        way that shifting from low to high precision is handled, by just reading the extra segments in the <true> variant.
        
        <false> specialisation - Operations performed on values in the array only modify the "head" segment (i.e. the first 32 bits of a double value) unless specified otherwise.
        <true> specialisation - Operations performed on values in the array are exactly the same as standard IEEE-754 doubles, but a bit slower due to combining of head and tail segments.
    */
    template<bool useTail>
    class TwoSegArray; // explicitly specialise this below.

    /*
        Specialisation of TwoSegArray.
        User is required to manage de-allocation of memory manually, using the del()
        function, as the deconstructor does not free this space.
        Operations performed on values in the array are exactly the same as standard IEEE-754 doubles, but a bit slower due to combining of head and tail segments.
    */
    template<>
    class TwoSegArray<true>
    {
    public:
        TwoSegArray() {}

        TwoSegArray(const uint_fast64_t& length)
        {
            heads = new uint32_t[length];
            tails = new uint32_t[length];
        }

        TwoSegArray(uint32_t* heads, uint32_t* tails)
            :heads(heads), tails(tails)
        {}

        ~TwoSegArray()
        {
            heads = nullptr;
            tails = nullptr;
        }

        template<typename T>
        void set(const uint_fast64_t& id, const T& t)
        {
            double d = t;
            const doublerep l = *reinterpret_cast<const doublerep*>(&d);
            heads[id] = l >> segmentBits;
            tails[id] = (l & tailMask);
        }

        template<typename T>
        void setPair(const uint_fast64_t& id, const T& t)
        {
            double d = t;
            const doublerep l = *reinterpret_cast<const doublerep*>(&d);
            heads[id] = l >> segmentBits;
            tails[id] = (l & tailMask);
        }

        double read(const uint_fast64_t& id)
        {
            return static_cast<double>(Pair(heads[id], tails[id]));
        }

        Pair operator[](const uint_fast64_t& id)
        {
            return Pair(heads[id], tails[id]);
        }

        /*
            Returns *this (as we do not have any precision increase to do)
        */
        TwoSegArray<true> createFullPrecision()
        {
            return *this;
        }

        void alloc(const uint_fast64_t& length)
        {
            heads = new uint32_t[length];
            tails = new uint32_t[length] ();
        }

        /*
            Deletes the values of the dynamic arrays used to store values in
            the array.
            NOTE: this should only be called by one object with references to
            the same set of values (such as object created using createFullPrecision).
        */
        void del()
        {
            if(heads != nullptr) delete[] heads;
            if(tails != nullptr) delete[] tails;
        }

    private:
        uint32_t* heads;
        uint32_t* tails;

        static constexpr int segmentBits = 32;
        static constexpr uint32_t tailMask = ~0;
        static constexpr doublerep headMask = static_cast<doublerep>(tailMask) << segmentBits;

    };

    /*
        Specialisation of TwoSegArray.
        User is required to manage de-allocation of memory manually, using the del()
        function, as the deconstructor does not free this space.
        Operations performed on values in the array only modify the "head" segment
        (i.e. the first 32 bits of a double value) unless specified otherwise.
    */
    template<>
    class TwoSegArray<false>
    {
    public:
        TwoSegArray() {}

        TwoSegArray(const uint_fast64_t& length)
        {
            heads = new uint32_t[length];
            tails = new uint32_t[length] ();
        }

        TwoSegArray(uint32_t* heads, uint32_t* tails)
            :heads(heads), tails(tails)
        {}

        ~TwoSegArray()
        {
            heads = nullptr;
            tails = nullptr;
        }

        template<typename T>
        void set(const uint_fast64_t& id, const T& t)
        {
            double d = t;
            const doublerep l = *reinterpret_cast<const doublerep*>(&d);
            heads[id] = (l >> segmentBits);
        }

        template<typename T>
        void setPair(const uint_fast64_t& id, const T& t)
        {
            double d = t;
            const doublerep l = *reinterpret_cast<const doublerep*>(&d);
            heads[id] = l >> segmentBits;
            tails[id] = (l & tailMask);
        }

        double read(const uint_fast64_t& id)
        {
            return static_cast<double>(Head(heads[id]));
        }

        Head operator[](const uint_fast64_t& id)
        {
            return Head(heads[id]);
        }

        /*
            Increases the precision of the operations performed by returning an
            object of type TwoSegArray<true> that has pointers to the values
            in the existing array.
        */
        TwoSegArray<true> createFullPrecision()
        {
            return TwoSegArray<true>(heads, tails);
        }

        void alloc(const uint_fast64_t& length)
        {
            heads = new uint32_t[length];
            // we should zero tails when allocating heads array for
            // to avoid unexpected behaviour
            tails = new uint32_t[length] ();
        }

        /*
            Deletes the values of the dynamic arrays used to store values in
            the array.
            NOTE: this should only be called by one object with references to
            the same set of values (such as object created using createFullPrecision).
        */
        void del()
        {
            if(heads != nullptr) delete[] heads;
            if(tails != nullptr) delete[] tails;
        }

    private:
        uint32_t* heads;
        uint32_t* tails;

        static constexpr int segmentBits = 32;
        static constexpr uint32_t tailMask = ~0;
        static constexpr doublerep headMask = static_cast<doublerep>(tailMask) << segmentBits;

    };


    /**
     * Head functions
    */
    template<>
    inline Head& Head::operator=(const Head& other)
    {
        uint32_t h = other.head;
        head = h;
        return *this;
    }
    
    template<>
    inline Head& Head::operator=(const Head&& other) noexcept
    {
        uint32_t h = other.head;
        head = h;
        return *this;
    }

    template<>
    inline Head& Head::operator=(const Pair& other)
    {
        uint32_t h = other.head;
        head = h;
        return *this;
    }

    template<>
    inline Head& Head::operator=(const Pair&& other) noexcept
    {
        uint32_t h = other.head;
        head = h;
        return *this;
    }

    template<typename T>
    inline Head& Head::operator=(const T& rhs)
    {
        double d = rhs;
        const doublerep l = *reinterpret_cast<const doublerep*>(&d);
        head = (l >> segmentBits);
        return *this;
    }

    template<typename T>
    inline Head& Head::operator=(const T&& other) noexcept
    {
        double d = other;
        const doublerep l = *reinterpret_cast<const doublerep*>(&d);
        uint32_t h = (l >> segmentBits);

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
    inline Pair& Pair::operator=(const Head& other)
    {
        uint32_t h = other.head;
        head = h;
        return *this;
    }
    
    template<>
    inline Pair& Pair::operator=(const Head&& other) noexcept
    {
        uint32_t h = other.head;
        head = h;
        return *this;
    }

    template<>
    inline Pair& Pair::operator=(const Pair& other)
    {
        uint32_t h = other.head;
        uint32_t t = other.tail;
        head = h;
        tail = t;
        return *this;
    }
    
    template<>
    inline Pair& Pair::operator=(const Pair&& other) noexcept
    {
        uint32_t h = other.head;
        uint32_t t = other.tail;
        head = h;
        tail = t;
        return *this;
    }

    template<typename T>
    inline Pair& Pair::operator=(const T&& other) noexcept
    {
        double d = other;
        const doublerep l = *reinterpret_cast<const doublerep*>(&d);
        uint32_t h = (l >> segmentBits);
        uint32_t t = (l & tailMask);
        head = h;
        tail = t;
        return *this;
    }

    template<typename T>
    inline Pair& Pair::operator=(const T& other)
    {
        double d = other;
        const doublerep l = *reinterpret_cast<const doublerep*>(&d);
        uint32_t h = (l >> segmentBits);
        uint32_t t = (l & tailMask);
        head = h;
        tail = t;
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

    /*
        Convenient type for use of TwoSegArray<false> and TwoSegArray<true> without having to manage two separate sets of arrays.
        
        @heads - used to access only the first 32 bits of a double : [sign(1), exp(11), mantissa(20)]
        @pairs - used to access all 64 bits of a double : [sign(1), exp(11), mantissa(52)]
    */
    class ManSegArray
    {
    public:
        TwoSegArray<false> heads;
        TwoSegArray<true> pairs;

        ManSegArray() {}

        ManSegArray(const uint_fast64_t& length)
        {
            heads.alloc(length);
            pairs = heads.createFullPrecision();
        }

        /*
            Allocates length elements to the array dynamically.
            Note: this array space is used for both heads and pairs
        */
        void alloc(const uint_fast64_t& length)
        {
            heads.alloc(length);
            pairs = heads.createFullPrecision();
        }

        /* 
            Deletes space allocated to the array.
            WARNING: should only be called once, as heads and pairs share the array space.
        */
        void del() { heads.del(); }
    };
}

#endif
