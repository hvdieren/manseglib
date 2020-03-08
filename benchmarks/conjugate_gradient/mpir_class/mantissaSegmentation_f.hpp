#pragma once

#include <stdint.h>
#include <immintrin.h>

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
            __m128 seg_v = _mm_set_ps(0.0f, 0.0f, head, 0.0f);
            __m128d d_v = _mm_castps_pd(seg_v);
            return d_v[0];

            // __m128 f_v = _mm_load1_ps((&head));
            // // cast to double now
            // __m128d d_v = _mm_castps_pd(f_v);
            // // cast 64-bit integer mask to double
            // __m128d d_mask = _mm_castsi128_pd(l_mask); // l_mask member of class (or library?)
            // // do AND of parts
            // __m128d res = _mm_and_pd(d_v, d_mask);

            // return res[0];
        }

        float& head;
        static constexpr int segmentBits = 32;
        static constexpr __m128i l_mask = {(long long)(0xFFFFFFFF00000000), (long long)(0x0)};

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
            __m128 seg_v = _mm_set_ps(0.0f, 0.0f, head, tail);
            __m128d d_v = _mm_castps_pd(seg_v);
            return d_v[0];

            // how do we do this properly ?

            // // load tail
            // __m128 f_v = _mm_load1_ps((&tail));
            // // how do we set the head..?
            // f_v[1] = head;   /// < ----- this is what you need to look at next in dev
            // // cast to double vector
            // __m128d d_v = _mm_castps_pd(f_v);
            // // cast 64-bit integer mask to double
            // __m128d d_mask = _mm_castsi128_pd(l_mask); // todo: make constexpr function to do this cast, so we can avoid it
            // // do AND of parts
            // __m128d res = _mm_and_pd(d_v, d_mask);

            // return res[0];
        }

        float& head;
        float& tail;
        static constexpr int segmentBits = 32;
        static constexpr __m128i l_mask = {(long long)(0xFFFFFFFFFFFFFFFF), (long long)(0x0)};

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

        TwoSegArray(const uint_fast64_t& length)
        {
            heads = new float[length];
            // we should zero tails when allocating heads array for
            // to avoid unexpected behaviour
            tails = new float[length] ();
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
        void set(const uint_fast64_t& id, const T& t);
        template<typename T>
        void setPair(const uint_fast64_t& id, const T& t);

        double read(const uint_fast64_t& id);

        Head operator[](const uint_fast64_t& id);

        /*
            Increases the precision of the operations performed by returning an
            object of type TwoSegArray<true> that has pointers to the values
            in the existing array.
        */
        TwoSegArray<true> createFullPrecision();

        void alloc(const uint_fast64_t& length);
        /*
            Deletes the values of the dynamic arrays used to store values in
            the array.
            NOTE: this should only be called by one object with references to
            the same set of values (such as object created using createFullPrecision).
        */
        void del();

    // private:
        float* heads;
        float* tails;

        static constexpr int segmentBits = 32;
        static constexpr uint32_t tailMask = ~0;
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

        TwoSegArray(const uint_fast64_t& length)
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

        Pair operator[](const uint_fast64_t& id);

        template<typename T>
        void set(const uint_fast64_t& id, const T& t);
        template<typename T>
        void setPair(const uint_fast64_t& id, const T& t);
        double read(const uint_fast64_t& id);

        /*
            Returns *this (as we do not have any precision increase to do)
        */
        TwoSegArray<true> createFullPrecision();

        void alloc(const uint_fast64_t& length);
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
        static constexpr uint32_t tailMask = ~0;
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
        __m128d d_v = _mm_set_pd(0.0, d);
        __m128 seg_v = _mm_castpd_ps(d_v);
        head = seg_v[1];

        // double d = other;
        // __m128 f_v = _mm_loadu_ps(reinterpret_cast<float*>(&d) + 1);
        // _mm_store_ss(&head, f_v);

        return *this;
    }

    template<typename T>
    Head& Head::operator=(const T&& other) noexcept
    {
        double d = other;
        __m128d d_v = _mm_set_pd(0.0, d);
        __m128 seg_v = _mm_castpd_ps(d_v);
        head = seg_v[1];

        // double d = other;
        // // i don't like this, it probably will fail at some point
        // __m128 f_v = _mm_loadu_ps(reinterpret_cast<float*>(&d) + 1);
        // _mm_store_ss(&head, f_v);

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
    Pair& Pair::operator=(const T& other)
    {
        double d = other;
        __m128d d_v = _mm_set_pd(0.0, d);
        __m128 seg_v = _mm_castpd_ps(d_v);
        tail = seg_v[0];
        head = seg_v[1];

        // double d = other;
        // // load double as float (should be free)
        // // tail of double is in f_v[0]
        // __m128 f_v = _mm_loadu_ps(reinterpret_cast<float*>(&d));
        // _mm_store_ps1(&tail, f_v);
        // // head is awkward
        // // re-load values with offset of +1 (+4 bytes)
        // // so head is in f_v[0]
        // f_v = _mm_loadu_ps(reinterpret_cast<float*>(&d) + 1); // -1 aligns tail in bottom
        // // store f_v[0] in head
        // _mm_store_ps1(&head, f_v);

        return *this;
    }

    template<typename T>
    Pair& Pair::operator=(const T&& other) noexcept
    {
        double d = other;
        __m128d d_v = _mm_set_pd(0.0, d);
        __m128 seg_v = _mm_castpd_ps(d_v);
        tail = seg_v[0];
        head = seg_v[1];

        // double d = other;
        // // load double as float (should be free)
        // // tail of double is in f_v[0]
        // __m128 f_v = _mm_loadu_ps(reinterpret_cast<float*>(&d));
        // _mm_store_ps1(&tail, f_v);
        // // head is awkward
        // // re-load values with offset of +1 (+4 bytes)
        // // so head is in f_v[0]
        // f_v = _mm_loadu_ps(reinterpret_cast<float*>(&d) + 1); // -1 aligns tail in bottom
        // // store f_v[0] in head
        // _mm_store_ps1(&head, f_v);

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
    void TwoSegArray<false>::set(const uint_fast64_t& id, const T& t)
    {
        double d = t;
        __m128d d_v = _mm_set_pd(0.0, d);
        __m128 seg_v = _mm_castpd_ps(d_v);
        heads[id] = seg_v[1];

        // double d = t;
        // // __m128d d_v = _mm_loadu_pd(&d);
        // // __m128 seg_v = _mm_castpd_ps(d_v);
        // // seg_v[0] = seg_v[1];
        // // _mm_store_ss(&heads[id], seg_v);
        // float* bits = reinterpret_cast<float*>(&d);
        // heads[id] = bits[1];
    }

    template<typename T>
    void TwoSegArray<false>::setPair(const uint_fast64_t& id, const T& t)
    {
        double d = t;
        __m128d d_v = _mm_set_pd(0.0, d);
        __m128 seg_v = _mm_castpd_ps(d_v);
        tails[id] = seg_v[0];
        heads[id] = seg_v[1];

        // double d = t;
        // __m128d d_v = _mm_loadu_pd(&d);
        // __m128 seg_v = _mm_castpd_ps(d_v);
        // _mm_store_ss(&tails[id], seg_v);
        // seg_v[0] = seg_v[1];
        // _mm_store_ss(&heads[id], seg_v);
        // float* bits = reinterpret_cast<float*>(&d);
        // tails[id] = bits[0];
        // heads[id] = bits[1];
    }

    double TwoSegArray<false>::read(const uint_fast64_t& id)
    {
        return static_cast<double>(Head(heads[id]));
    }

    Head TwoSegArray<false>::operator[](const uint_fast64_t& id)
    {
        return Head(heads[id]);
    }

    TwoSegArray<true> TwoSegArray<false>::createFullPrecision()
    {
        return TwoSegArray<true>(heads, tails);
    }

    void TwoSegArray<false>::alloc(const uint_fast64_t& length)
    {
        heads = new float[length];
        // we should zero tails when allocating heads array for
        // to avoid unexpected behaviour
        tails = new float[length] ();
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
    void TwoSegArray<true>::set(const uint_fast64_t& id, const T& t)
    {
        double d = t;
        __m128d d_v = _mm_set_pd(0.0, d);
        __m128 seg_v = _mm_castpd_ps(d_v);
        tails[id] = seg_v[0];
        heads[id] = seg_v[1];

        // double d = t;
        // // __m128d d_v = _mm_loadu_pd(&d);
        // // __m128 seg_v = _mm_castpd_ps(d_v);
        // // _mm_store_ss(&tails[id], seg_v);
        // // seg_v[0] = seg_v[1];
        // // _mm_store_ss(&heads[id], seg_v);
        // float* bits = reinterpret_cast<float*>(&d);
        // tails[id] = bits[0];
        // heads[id] = bits[1];
    }

    template<typename T>
    void TwoSegArray<true>::setPair(const uint_fast64_t& id, const T& t)
    {
        double d = t;
        __m128d d_v = _mm_set_pd(0.0, d);
        __m128 seg_v = _mm_castpd_ps(d_v);
        tails[id] = seg_v[0];
        heads[id] = seg_v[1];
        
        // double d = t;
        // __m128d d_v = _mm_loadu_pd(&d);
        // __m128 seg_v = _mm_castpd_ps(d_v);
        // _mm_store_ss(&tails[id], seg_v);
        // seg_v[0] = seg_v[1];
        // _mm_store_ss(&heads[id], seg_v);
        // float* bits = reinterpret_cast<float*>(&d);
        // tails[id] = bits[0];
        // heads[id] = bits[1];
    }

    double TwoSegArray<true>::read(const uint_fast64_t& id)
    {
        return static_cast<double>(Pair(heads[id], tails[id]));
    }

    Pair TwoSegArray<true>::operator[](const uint_fast64_t& id)
    {
        return Pair(heads[id], tails[id]);
    }

    TwoSegArray<true> TwoSegArray<true>::createFullPrecision()
    {
        return *this;
    }

    void TwoSegArray<true>::alloc(const uint_fast64_t& length)
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

        ManSegArray(const uint_fast64_t& length)
        {
            heads.alloc(length);
            pairs = heads.createFullPrecision();
        }

        void alloc(const uint_fast64_t& length)
        {
            heads.alloc(length);
            pairs = heads.createFullPrecision();
        }

        /* Deletes space allocated to arrays */
        void del() { heads.del(); }
    };
}
