// BLAS like functions
#ifndef __cg_VECTOR_H__
#define __cg_VECTOR_H__

static inline void vector_set(int n, DOUBLE a, PairsArray *x) {
    for (int i = 0; i < n; i++) (*x)[i] = a;
}

static inline void vector_rand(int n, PairsArray *x) {
    for (int i = 0; i < n; i++) (*x)[i] = rand() / (double)RAND_MAX;
}

static inline void vector_copy(int n, PairsArray *x, PairsArray *y) {
    for (int i = 0; i < n; i++) (*y)[i] = (*x)[i];
}

static inline void vector_axpy(int n, DOUBLE a, PairsArray *x, PairsArray *y) {
    for (int i = 0; i < n; i++) (*y)[i] = (*y)[i] + (*x)[i] * a;
}

static inline void vector_xpby(int n, PairsArray *x, DOUBLE b, PairsArray *y) {
    for (int i = 0; i < n; i++) (*y)[i] = (*x)[i] + b * (*y)[i];
}

static inline DOUBLE vector_dot(int n, PairsArray *x, PairsArray *y) {
    DOUBLE r = 0.0;
    for (int i = 0; i < n; i++) r += (*y)[i] * (*x)[i];
    return r;
}

static inline DOUBLE vector_norm2(int n, PairsArray *x) {
    return sqrt(vector_dot(n, x, x));
}

// limited precision version

static inline void floatm_set(int n, FLOAT a, HeadsArray *x) {
    FLOAT b = a;
    for (int i = 0; i < n; i++) (*x)[i] = b;
}

static inline void floatm_rand(int n, HeadsArray *x) {
    for (int i = 0; i < n; i++) (*x)[i] = rand() / (double)RAND_MAX;
}

static inline void floatm_copy(int n, HeadsArray *x, HeadsArray *y) {
    for (int i = 0; i < n; i++) (*y)[i] = (*x)[i];
}

static inline void floatm_axpy(int n, FLOAT a, HeadsArray *x, HeadsArray *y) {
    for (int i = 0; i < n; i++) (*y)[i] = (*y)[i] + (*x)[i] * a;
}

static inline void floatm_xpby(int n, HeadsArray *x, FLOAT b, HeadsArray *y) {
    for (int i = 0; i < n; i++) (*y)[i] = (*x)[i] + b * (*y)[i];
}

// the dot product can be computed with more precision than FLOAT

static inline FLOAT2 floatm_dot(int n, HeadsArray *x, HeadsArray *y) {
    FLOAT2 r = 0.0;
    for (int i = 0; i < n; i++) r += (*y)[i] * (*x)[i];
    return r;
}

static inline FLOAT2 floatm_norm2(int n, HeadsArray *x) {
    return sqrt(floatm_dot(n, x, x));
}

static inline FLOAT2 floatm_axpy_dot(int n, FLOAT a, HeadsArray *x, HeadsArray *y) {
    FLOAT2 r = 0.0;
    for (int i = 0; i < n; i++) {
        FLOAT2 t = (*y)[i] + (*x)[i] * a;
        (*y)[i] = t;
        r += t * t;
    }
    return r;
}

static inline FLOAT2 floatm_diff_norm2(int n, HeadsArray *x, HeadsArray *y) {
    FLOAT2 r = 0.0;
    for (int i = 0; i < n; i++) {
        FLOAT2 t = (*x)[i] - (*y)[i];
        r += t * t;
    }
    return sqrt(r);
}

// mixed precision routines

static inline void mixed_copy(int n, PairsArray *x, HeadsArray *y) {
    for (int i = 0; i < n; i++) (*y)[i] = (*x)[i];
}

static inline void mixed_axpy(int n, DOUBLE a, HeadsArray *x, PairsArray *y) {
    for (int i = 0; i < n; i++) (*y)[i] += (*x)[i] * a;
}

static inline void mixed_xpby(int n, PairsArray *x, DOUBLE b, HeadsArray *y) {
    for (int i = 0; i < n; i++) (*y)[i] = (*x)[i] + b * (*y)[i];
}

#endif
