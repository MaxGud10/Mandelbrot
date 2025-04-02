#pragma once

#include <SFML/Graphics.hpp>
#include <immintrin.h>
#include "stdlib.h"
#include <x86intrin.h>

#include <omp.h> // TODO убрать 


#define RUN_TEST(FUNC)                                                       \
    uint64_t start_ticks = __rdtsc ();                                       \
        for (size_t i = 0; i < ITERATIONS; i++)                              \
        {                                                                    \
            FUNC (set);                                                      \
        }                                                                    \
    uint64_t end_ticks   = __rdtsc ();                                       \
    printf ("CPU ticks = %ld \n", (end_ticks - start_ticks) / ITERATIONS);   \


    const int   MAX_ITERATIONS = 256;
    const float MAX_RADIUS     = 100.f;
    const int   WIDTH          = 800;
    const int   HEIGHT         = 600;
    const float D_X            = 1 / (float)WIDTH;
    const float D_Y            = 1 / (float)WIDTH;
    const float HALF_WIDTH     = (float)WIDTH  / 2;
    const float HALF_HEIGHT    = (float)HEIGHT / 2;
    const int   VECTOR_SIZE    = 8;
    const size_t ITERATIONS    = 80;

enum Mode 
{
    BY_PIXELS  =  0,
    BY_VECTOR  =  1,
    BY_SIMD    =  2
};

struct MandelBrot 
{
    float scale;
    float x_offset;
    float y_offset;

    Mode mode;

    uint32_t* pixels_array;
};

void        init_mandelbrot       (struct MandelBrot* set); 

void        mandelbrot_naive      (struct MandelBrot* set);
void        mandelbrot_vectorized (struct MandelBrot* set); 
void        mandelbrot_simd       (struct MandelBrot* set);

void        run_performance_test  (struct MandelBrot* set);

inline void get_pixels            (struct MandelBrot* set, size_t index, size_t count);
inline void get_fps               (sf::Clock &clock, sf::Text &text, struct MandelBrot* set);

void        get_mandel_brot_set   (struct MandelBrot* set);

inline  int handle_keyboard_input (struct MandelBrot* set, sf::Event &event); 

