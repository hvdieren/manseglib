#pragma once

#include <stdint.h>
#include <immintrin.h>

namespace ManSeg
{
    // internal representation of double, so it can be easily manipulated
    typedef std::uint_fast64_t doublerep;

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
        It has a maximum precision of roughly 1e-5
    */
    class Head
    {
    public:
        Head(float& head)
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
            __m128 seg_v = {0.0f, head};
            __m128d d_v = _mm_castps_pd(seg_v);
            return *reinterpret_cast<double*>(&d_v);
        }

        float& head;
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
        Pair(float& head, float& tail)
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
            __m128 seg_v = {tail, head};
            __m128d d_v = _mm_castps_pd(seg_v);
            return *reinterpret_cast<double*>(&d_v);
        }

        float& head;
        float& tail;
        static constexpr int segmentBits = 32;
        static constexpr std::uint32_t tailMask = ~0;

    };

    /*
        Array type class for performing operations on double values, conceptually split into two 32 bit segments - "head" and "tail".
        false specialisation - Operations performed on values in the array only modify the "head" segment (i.e. the first 32 bits of a double value) unless specified otherwise.
        true specialisation - Operations performed on values in the array are exactly the same as standard IEEE-754 doubles, but a bit slower due to combining of head and tail segments.
    */
    template<bool useTail>
    class TwoSegArray; // we specialise this below.

    /*
        Specialisation of TwoSegArray.
        Operations performed on values in the array only modify the "head" segment
        (i.e. the first 32 bits of a double value) unless specified otherwise.
    */
    template<>
    class TwoSegArray<false>
    {
    public:
        TwoSegArray() {}

        TwoSegArray(const std::uint_fast64_t& length)
        {
            heads = new float[length];
            tails = new float[length];
        }

        TwoSegArray(float* heads, float* tails)
            :heads(heads), tails(tails)
        {}

        ~TwoSegArray()
        {
            heads = nullptr;
            tails = nullptr;
        }

        template<typename T>
        void set(const std::uint_fast64_t& id, const T& t);
        template<typename T>
        void setPair(const std::uint_fast64_t& id, const T& t);

        double read(const std::uint_fast64_t& id);

        Head operator[](const std::uint_fast64_t& id);

        /*
            Increases the precision of the operations performed by returning an
            object of type TwoSegArray<true> that has pointers to the values
            in the existing array.
        */
        TwoSegArray<true> createFullPrecision();

        void alloc(const std::uint_fast64_t& length);
        /*
            Deletes the values of the dynamic arrays used to store values in
            the array.
            NOTE: this should only be called by one object with references to
            the same set of values (such as object created using createFullPrecision).
        */
        void del();

    private:
        float* heads;
        float* tails;

        static constexpr int segmentBits = 32;
        static constexpr std::uint32_t tailMask = ~0;
        static constexpr doublerep headMask = static_cast<doublerep>(tailMask) << segmentBits;

    };

    /*
        Specialisation of TwoSegArray.
        Operations performed on values in the array are exactly the same as standard IEEE-754 doubles, but a bit slower due to combining of head and tail segments.
    */
    template<>
    class TwoSegArray<true>
    {
    public:
        TwoSegArray() {}

        TwoSegArray(const std::uint_fast64_t& length)
        {
            heads = new float[length];
            tails = new float[length];
        }

        TwoSegArray(float* heads, float* tails)
            :heads(heads), tails(tails)
        {}

        ~TwoSegArray()
        {
            heads = nullptr;
            tails = nullptr;
        }

        Pair operator[](const std::uint_fast64_t& id);

        template<typename T>
        void set(const std::uint_fast64_t& id, const T& t);
        template<typename T>
        void setPair(const std::uint_fast64_t& id, const T& t);
        double read(const std::uint_fast64_t& id);

        /*
            Returns *this (as we do not have any precision increase to do)
        */
        TwoSegArray<true> createFullPrecision();

        void alloc(const std::uint_fast64_t& length);
        /*
            Deletes the values of the dynamic arrays used to store values in
            the array.
            NOTE: this should only be called by one object with references to
            the same set of values (such as object created using createFullPrecision).
        */
        void del();

    private:
        float* heads;
        float* tails;

        static constexpr int segmentBits = 32;
        static constexpr std::uint32_t tailMask = ~0;
        static constexpr doublerep headMask = static_cast<doublerep>(tailMask) << segmentBits;

    };


    /**
     * Head functions
    */
    template<>
    Head& Head::operator=(const Head& other)
    {
        float h = other.head;
        head = h;
        return *this;
    }

    template<>
    Head& Head::operator=(const Head&& other) noexcept
    {
        float h = other.head;
        head = h;
        return *this;
    }

    template<>
    Head& Head::operator=(const Pair& other)
    {
        float h = other.head;
        head = h;
        return *this;
    }

    template<>
    Head& Head::operator=(const Pair&& other) noexcept
    {
        float h = other.head;
        head = h;
        return *this;
    }

    template<typename T>
    Head& Head::operator=(const T& other)
    {
        double d = other;
        __m128d d_v = {d};
        __m128 seg_v = _mm_castpd_ps(d_v);
        head = seg_v[1];

        return *this;
    }

    template<typename T>
    Head& Head::operator=(const T&& other) noexcept
    {
        double d = other;
        __m128d d_v = {d};
        __m128 seg_v = _mm_castpd_ps(d_v);
        head = seg_v[1];

        // __m128d seg_v = {d, 0.0};
        // __m128 seg_fv = *reinterpret_cast<__m128*>(&seg_v);
        // head = seg_fv[1];

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
        float h = other.head;
        head = h;
        return *this;
    }

    template<>
    Pair& Pair::operator=(const Head&& other) noexcept
    {
        float h = other.head;
        head = h;
        return *this;
    }

    Pair& Pair::operator=(const Pair& other)
    {
        float h = other.head;
        float t = other.tail;
        head = h;
        tail = t;
        return *this;
    }

    template<>
    Pair& Pair::operator=(const Pair&& other) noexcept
    {
        float h = other.head;
        float t = other.tail;
        head = h;
        tail = t;
        return *this;
    }

    template<typename T>
    Pair& Pair::operator=(const T&& other) noexcept
    {
        double d = other;
        __m128d d_v = {d};
        __m128 seg_v = _mm_castpd_ps(d_v);
        tail = seg_v[0];
        head = seg_v[1];

        return *this;
    }

    template<typename T>
    Pair& Pair::operator=(const T& other)
    {
        double d = other;
        __m128d d_v = {d};
        __m128 seg_v = _mm_castpd_ps(d_v);
        tail = seg_v[0];
        head = seg_v[1];

        return *this;
    }

    template <typename T>
    double Pair::operator+=(const T& rhs)
    {
        double t = *this;
        double r = rhs;
        t += r;
        *this = t;
        return t;
    }

    template<typename T>
    inline double operator+(Pair lhs, const T& rhs)
    {
        double d = lhs;
        double r = rhs;
        d += r;
        return d;
    }

    template<typename T>
    double Pair::operator-=(const T& rhs)
    {
        double t = *this;
        double r = rhs;
        t -= r;
        *this = t;
        return t;
    }

    template<typename T>
    inline double operator-(Pair lhs, const T& rhs)
    {
        double d = lhs;
        double r = rhs;
        d -= r;
        return d;
    }

    template<typename T>
    double Pair::operator*=(const T& rhs)
    {
        double t = *this;
        double r = rhs;
        t *= r;
        *this = t;
        return t;
    }

    template<typename T>
    inline double operator*(Pair lhs, const T& rhs)
    {
        double d = lhs;
        double r = rhs;
        d *= r;
        return d;
    }

    template<typename T>
    double Pair::operator/=(const T& rhs)
    {
        double t = *this;
        double r = rhs;
        t /= r;
        *this = t;
        return t;
    }

    template<typename T>
    inline double operator/(Pair lhs, const T& rhs)
    {
        double d = lhs;
        double r = rhs;
        d /= r;
        return d;
    }


    /**
     * TwoSegArray functions
    */

    /**
     * false specialisation
     * (heads only)
    */
    template<typename T>
    void TwoSegArray<false>::set(const std::uint_fast64_t& id, const T& t)
    {
        double d = t;
        __m128d d_v = {d};
        __m128 seg_v = _mm_castpd_ps(d_v);
        heads[id] = seg_v[1];
    }

    template<typename T>
    void TwoSegArray<false>::setPair(const std::uint_fast64_t& id, const T& t)
    {
        double d = t;
        __m128d d_v = {d};
        __m128 seg_v = _mm_castpd_ps(d_v);
        tails[id] = seg_v[0];
        heads[id] = seg_v[1];
    }

    double TwoSegArray<false>::read(const std::uint_fast64_t& id)
    {
        return static_cast<double>(Head(heads[id]));
    }

    Head TwoSegArray<false>::operator[](const std::uint_fast64_t& id)
    {
        return Head(heads[id]);
    }

    TwoSegArray<true> TwoSegArray<false>::createFullPrecision()
    {
        return TwoSegArray<true>(heads, tails);
    }

    void TwoSegArray<false>::alloc(const std::uint_fast64_t& length)
    {
        heads = new float[length];
        tails = new float[length];
    }

    void TwoSegArray<false>::del()
    {
        if(heads) delete[] heads;
        if(tails) delete[] tails;
    }


    /**
     * true specialisation
     * (heads + tails)
    */
    template<typename T>
    void TwoSegArray<true>::set(const std::uint_fast64_t& id, const T& t)
    {
        double d = t;
        __m128d d_v = {d};
        __m128 seg_v = _mm_castpd_ps(d_v);
        tails[id] = seg_v[0];
        heads[id] = seg_v[1];
    }

    template<typename T>
    void TwoSegArray<true>::setPair(const std::uint_fast64_t& id, const T& t)
    {
        double d = t;
        __m128d d_v = {d};
        __m128 seg_v = _mm_castpd_ps(d_v);
        tails[id] = seg_v[0];
        heads[id] = seg_v[1];
    }

    double TwoSegArray<true>::read(const std::uint_fast64_t& id)
    {
        return static_cast<double>(Pair(heads[id], tails[id]));
    }

    Pair TwoSegArray<true>::operator[](const std::uint_fast64_t& id)
    {
        return Pair(heads[id], tails[id]);
    }

    TwoSegArray<true> TwoSegArray<true>::createFullPrecision()
    {
        return *this;
    }

    void TwoSegArray<true>::alloc(const std::uint_fast64_t& length)
    {
        heads = new float[length];
        tails = new float[length];
    }

    void TwoSegArray<true>::del()
    {
        if(heads) delete[] heads;
        if(tails) delete[] tails;
    }

    /*
        Convenience object - allows use of heads & pairs without having to manage
        two separate objects, or any of the eccentricities that arise from doing so.
        
        heads - single segment calculations only (i.e. the first 32 bits of a double)
        pairs - two segment calculations (i.e. the entire 64 bits of a double)
    */
   
    // todo: add this to mantissa-segmentation project
    class ManSegArray
    {
    public:
        TwoSegArray<false> heads;
        TwoSegArray<true> pairs;

        ManSegArray() {}

        ManSegArray(const std::uint_fast64_t& length)
        {
            heads.alloc(length);
            pairs = heads.createFullPrecision();
        }

        void alloc(const std::uint_fast64_t& length)
        {
            heads.alloc(length);
            pairs = heads.createFullPrecision();
        }

        /* Deletes space allocated to arrays */
        void del() { heads.del(); }
    };
}
