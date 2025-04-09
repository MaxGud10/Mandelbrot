#pragma once

//#define DEBUG

#include <SFML/Graphics.hpp>
#include <immintrin.h>
#include "stdlib.h"
#include <string.h>
#include <x86intrin.h>
#include <math.h>

#include <omp.h> 

#define RUN_TEST(FUNC, mode)                                                        \
    if (mode == 1)                                                                  \
        {                                                                           \
        uint64_t start_ticks = __rdtsc();                                           \
            for (size_t i = 0; i < ITERATIONS; i++)                                 \
                {                                                                   \
                FUNC(set);                                                          \
                }                                                                   \
        uint64_t end_ticks   = __rdtsc();                                           \
        printf("CPU ticks = %ld \n", (end_ticks - start_ticks) / ITERATIONS);       \
        }                                                                           \
    else                                                                            \
        {                                                                           \
        sf::Clock clock = {};                                                       \
        clock.restart();                                                            \
        for (size_t i = 0; i < ITERATIONS; i++)                                     \
                {                                                                   \
                FUNC(set);                                                          \
                }                                                                   \
        printf("time = %f sec\n", clock.getElapsedTime().asSeconds() / ITERATIONS); \
        } 

#ifdef DEBUG 
    #define DBG(...) __VA_ARGS__
#else
    #define DBG(...)
#endif

#define FOR_VEC for (size_t i = 0; i < VECTOR_SIZE; ++i)
#define ALIGN alignas(32)  // 32 ����� ��� AVX


//=== Rendering constants ===//
const int    MAX_ITERATIONS = 256;               // ������������ ����� �������� ��� ��������
const float  MAX_RADIUS     = 100.f;             // �������� ������������ (|z|^2 >= 100)
const int    WIDTH          = 800;               // ������ ����
const int    HEIGHT         = 600;               // ������ ���� 
const float  D_X            = 1 / (float)WIDTH;  // ��� �� X � ����������� ���������
const float  D_Y            = 1 / (float)WIDTH;  // ��� �� Y
const float  HALF_WIDTH     = (float)WIDTH  / 2; // �������� ������
const float  HALF_HEIGHT    = (float)HEIGHT / 2; // �������� ������
const int    VECTOR_SIZE    = 32; // 32          // ������ ������� ��� SIMD (AVX2 = 8 float)
const size_t ITERATIONS     = 80;                // ����� ������ ��� ������ ������������������

//=== Color options ===//
const uint32_t DEFAULT_COLOR  = 0xFF000000; // ������ (��������)
const uint32_t ESCAPE_COLOR   = 0xFFFFFFFF; // ����� (����������)
const float    COLOR_SCALE    = 100.0f;     // ��������� ��� ���������

//=== Control Parameters ===//
const float ZOOM_IN_FACTOR    = 1.25f;      // ��������� ���������� (������� '+')
const float ZOOM_OUT_FACTOR   = 0.8f;       // ��������� ���������� (������� '-') 
const float MOVE_SPEED        = 0.01f;      // ��� �������� ������
const float INITIAL_SCALE     = 3.0f;       // ��������� �������
const float INITIAL_X_OFFSET  = -0.1f;      // ��������� �������� �� X
const float INITIAL_Y_OFFSET  = 0.0f;       // ��������� �������� �� Y

typedef struct MandelBrot MandelBrot_t;

typedef void (*MandelbrotFuncPtr)(MandelBrot_t*);


enum Mode 
{
    BY_PIXELS  =  1,
    BY_VECTOR  =  2,
    BY_SIMD    =  3
};

struct MandelBrot 
{
    float scale;
    float x_offset;
    float y_offset;

    Mode              mode;
    MandelbrotFuncPtr calculate;

    uint32_t* pixels_array;
    uint32_t* color_palette; 
};

void        init_mandelbrot       (MandelBrot_t* set); 
void        dtor_mandlebrot       (MandelBrot_t* set);
int         handle_keyboard_input (MandelBrot_t* set, sf::Event &event); 

void        mandelbrot_naive      (MandelBrot_t* set);
void        mandelbrot_vectorized (MandelBrot_t* set); 
void        mandelbrot_simd       (MandelBrot_t* set);
void        run_performance_test  (MandelBrot_t* set, int mode_measure);

void        get_mandel_brot_set   (MandelBrot_t* set);
void        get_pixels            (MandelBrot_t* set, size_t index, size_t count);
void        get_fps               (sf::Clock &clock, sf::Text &text, MandelBrot_t* set);
void        init_color_palette    (MandelBrot_t* set);


// TODO ��������� ������� -> gitatrid.

