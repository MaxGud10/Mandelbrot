#pragma once

//#define NDEBUG

#include <SFML/Graphics.hpp>
#include <immintrin.h>
#include "stdlib.h"
#include <x86intrin.h>

#include <omp.h> // TODO убрать 


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

#ifdef NDEBUG
    #define DBG(...) __VA_ARGS__
#else
    #define DBG(...)
#endif

//=== Rendering constants ===//
const int   MAX_ITERATIONS = 256;               // максимальное число итераций для фрактала
const float MAX_RADIUS     = 100.f;             // критерий расходимости (|z|^2 >= 100)
const int   WIDTH          = 800;               // ширина окна
const int   HEIGHT         = 600;               // высота окна 
const float D_X            = 1 / (float)WIDTH;  // шаг по X в комплексной плоскости
const float D_Y            = 1 / (float)WIDTH;  // шаг по Y
const float HALF_WIDTH     = (float)WIDTH  / 2; // половина ширины
const float HALF_HEIGHT    = (float)HEIGHT / 2; // половина высота
const int   VECTOR_SIZE    = 8;                 // размер вектора для SIMD (AVX2 = 8 float)
const size_t ITERATIONS    = 80;                // число тестов для замера производительности

//=== Color options ===//
const uint32_t DEFAULT_COLOR  = 0xFF000000; // чёрный (сходится)
const uint32_t ESCAPE_COLOR   = 0xFFFFFFFF; // белый (расходится)
const float    COLOR_SCALE    = 100.0f;     // [x]множитель для градиента

//=== Control Parameters ===//
const float ZOOM_IN_FACTOR    = 1.25f;      // множитель увеличения (клавиша '+')
const float ZOOM_OUT_FACTOR   = 0.8f;       // множитель уменьшения (клавиша '-') 
const float MOVE_SPEED        = 0.05f;      // шаг смещения камеры
const float INITIAL_SCALE     = 3.0f;       // начальный масштаб
const float INITIAL_X_OFFSET  = -0.1f;      // начальное смещение по X
const float INITIAL_Y_OFFSET  = 0.0f;       // начальное смещение по Y

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

    Mode mode;

    uint32_t* pixels_array;
};

void        init_mandelbrot       (struct MandelBrot* set); 
void        dtor_mandlebrot       (struct MandelBrot* set);

void        mandelbrot_naive      (struct MandelBrot* set);
void        mandelbrot_vectorized (struct MandelBrot* set); 
void        mandelbrot_simd       (struct MandelBrot* set);

void        run_performance_test  (struct MandelBrot* set, int mode_measure);

inline void get_pixels            (struct MandelBrot* set, size_t index, size_t count);
inline void get_fps               (sf::Clock &clock, sf::Text &text, struct MandelBrot* set);
inline  int handle_keyboard_input (struct MandelBrot* set, sf::Event &event); 

void        get_mandel_brot_set   (struct MandelBrot* set);



// TODO сделвть редми. прям как лаба (что хотел тоха конкретно уточнить в записях)
// TODO исправить русский -> gitatrid.
// TODO разделить на файлы
// TODO мб сделать typdef fot struct
///// TODO сделать больше констант 
///// TODO поменять NDEBUG и вынести в .h 

// [x]  как сделаю написать новые туду 
// TODO
// 