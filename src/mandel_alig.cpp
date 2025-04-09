#include "../include/mandlebrot.h"


void init_mandelbrot (MandelBrot_t* set) 
{
    set->scale        =  INITIAL_SCALE;
    set->x_offset     =  INITIAL_X_OFFSET;
    set->y_offset     =  INITIAL_Y_OFFSET ;

    set->pixels_array = (uint32_t*) aligned_alloc (32, HEIGHT * WIDTH * sizeof(uint32_t));

    memset(set->pixels_array, 0, HEIGHT * WIDTH * sizeof(uint32_t));

    init_color_palette(set);

    switch(set->mode) 
    {
        case BY_PIXELS:  set->calculate = mandelbrot_naive;       break;

        case BY_VECTOR:  set->calculate = mandelbrot_vectorized;  break;

        case BY_SIMD:    set->calculate = mandelbrot_simd;        break;

        default:         printf ("\nERROR: Invalid mode\n");      break;
    }
}

void get_mandel_brot_set (MandelBrot_t* set)
{
    sf::RenderWindow window (sf::VideoMode(WIDTH, HEIGHT), "Mandelbrot Set");

    sf::Image image;
    image.create (WIDTH, HEIGHT, sf::Color::Black);

    sf::Texture texture = {};
    texture.loadFromImage (image);

    sf::Sprite sprite   = {};
    sprite.setTexture     (texture);

    sf::Font font;
    font.loadFromFile ("fyodor-bold-oblique.ttf");
    sf::Text text     ("", font, 20);
    text.setFillColor (sf::Color::Blue);
    text.setPosition  (25, 5);

    sf::Clock clock = {};

    while (window.isOpen ())
    {
        sf::Event event;

        if (window.pollEvent (event))
        {
            if (event.type == sf::Event::Closed) 
            {
                DBG( printf("Window close event received\n"); )

                window.close();
            }

            handle_keyboard_input (set, event);
        }

        get_fps (clock, text, set);

        texture.update ((sf::Uint8*) set->pixels_array, WIDTH, HEIGHT, 0, 0);
        window.clear   ();
        window.draw    (sprite);
        window.draw    (text);
        window.display ();
    }
}


void mandelbrot_naive (MandelBrot_t* set)
{
    const float real_dx = set->scale * D_X;

    for (size_t y = 0; y < HEIGHT; y++)
    {
        float x0 = (-(HALF_WIDTH)            * D_X + set->x_offset) * set->scale;  // ПОПРОбОВАТЬ DOUBLE
        float y0 = (((float)y - HALF_HEIGHT) * D_Y + set->y_offset) * set->scale;

        for (size_t x = 0; x < WIDTH; x++, x0 += real_dx)
        {
            float x_n = x0;
            float y_n = y0;

            size_t count = 0;

            while (count++ < MAX_ITERATIONS)
            {
                float x2 = x_n * x_n;
                float y2 = y_n * y_n;
                float xy = x_n * y_n;

                if ((x2 + y2) >= MAX_RADIUS) break;

                x_n = x2 - y2 + x0;
                y_n = xy + xy + y0;
            }

            get_pixels (set, y * WIDTH + x, count - 1);
        }
    }
}

// версия с align 
inline void mandelbrot_vectorized(MandelBrot_t* set)
{
    const float real_dx = set->scale * D_X;
    const float real_dy = set->scale * D_Y;

    ALIGN float X0        [VECTOR_SIZE];
    ALIGN float X_N       [VECTOR_SIZE];
    ALIGN float Y_N       [VECTOR_SIZE];
    ALIGN float X2        [VECTOR_SIZE];
    ALIGN float Y2        [VECTOR_SIZE];
    ALIGN float XY        [VECTOR_SIZE];
    ALIGN int   real_count[VECTOR_SIZE];
    ALIGN int   cmp       [VECTOR_SIZE];

    for (size_t y = 0; y < HEIGHT; y++) 
    {
        const float base_x0 = (-(HALF_WIDTH) * D_X + set->x_offset) * set->scale;
        const float y0 = (((float)y - HALF_HEIGHT) * D_Y + set->y_offset) * set->scale;

        for (size_t x = 0; x < WIDTH; x += VECTOR_SIZE) 
        {

            // #pragma omp simd aligned(X0, X_N, Y_N:32)
            FOR_VEC 
            {
                X0[i] = base_x0 + (x + i) * real_dx;
                X_N[i] = X0[i];
                Y_N[i] = y0;
                real_count[i] = 0;
            }

            int iter = 0;
            while (iter++ < MAX_ITERATIONS)
            {

                // #pragma omp simd aligned(X_N, Y_N, X2, Y2, XY:32)
                FOR_VEC 
                {
                    X2[i] = X_N[i] * X_N[i];
                    Y2[i] = Y_N[i] * Y_N[i];
                    XY[i] = X_N[i] * Y_N[i];
                }

                bool all_inactive = true;
                // #pragma omp simd aligned(X2, Y2, cmp:32) reduction(||:all_inactive)
                FOR_VEC 
                {
                    cmp[i] = (X2[i] + Y2[i] <= MAX_RADIUS);
                    all_inactive = all_inactive && !cmp[i];
                }

                if (all_inactive) break;


                // #pragma omp simd aligned(X_N, Y_N, X2, Y2, XY, X0, cmp, real_count:32)
                FOR_VEC 
                {
                    real_count[i] += cmp[i];
                    X_N[i] = X2[i] - Y2[i] + X0[i];
                    Y_N[i] = XY[i] + XY[i] + y0;
                }
            }

            // #pragma omp simd aligned(real_count:32)
            for (size_t i = 0; i < VECTOR_SIZE && (x + i) < WIDTH; i++) {
                get_pixels(set, y * WIDTH + x + i, real_count[i]);
            }
        }
    }
}

// версия aligned_alloc
// inline void mandelbrot_simd(MandelBrot_t* set)
// {
//     const int SIMD_WIDTH = 8; 
//     float real_dx = D_X * set->scale;

//     __m256 MaxRadius = _mm256_set1_ps(MAX_RADIUS);
    
//     float*    x0     = (float*)    aligned_alloc (32, SIMD_WIDTH * sizeof(float)   );
//     uint32_t* counts = (uint32_t*) aligned_alloc (32, SIMD_WIDTH * sizeof(uint32_t));

//     for (size_t y = 0; y < HEIGHT; y++)
//     {
//         float base_x0 = (-(HALF_WIDTH) * D_X + set->x_offset) * set->scale;
//         float y0 = (((float)y - HALF_HEIGHT) * D_Y + set->y_offset) * set->scale;

//         for (size_t x = 0; x < WIDTH; x += VECTOR_SIZE)
//         {
//             for (size_t chunk = 0; chunk < VECTOR_SIZE; chunk += SIMD_WIDTH)
//             {
//                 size_t remaining = (chunk + SIMD_WIDTH > VECTOR_SIZE) ? 
//                                     VECTOR_SIZE - chunk : SIMD_WIDTH;
                
//                 for (size_t i = 0; i < remaining; i++) 
//                 {
//                     x0[i] = base_x0 + (x + chunk + i) * real_dx;
//                 }   

//                 for (size_t i = remaining; i < SIMD_WIDTH; i++) 
//                 {
//                     x0[i] = 0.0f;
//                 }

//                 __m256 X0 = _mm256_load_ps(x0);
//                 __m256 Y0 = _mm256_set1_ps(y0);

//                 __m256 X_N = X0;
//                 __m256 Y_N = Y0;

//                 __m256i real_count = _mm256_setzero_si256();

//                 int count = 0;
//                 while (count++ < MAX_ITERATIONS)
//                 {
//                     __m256 X2 = _mm256_mul_ps(X_N, X_N);
//                     __m256 Y2 = _mm256_mul_ps(Y_N, Y_N);
//                     __m256 XY = _mm256_mul_ps(X_N, Y_N);
//                     __m256 R2 = _mm256_add_ps(X2, Y2);

//                     __m256 mask = _mm256_cmp_ps(R2, MaxRadius, _CMP_LE_OS);
//                     if (_mm256_testz_ps(mask, mask)) break;

//                     __m256i cmp_result = _mm256_castps_si256(mask);
//                     cmp_result = _mm256_srli_epi32(cmp_result, 31);
//                     real_count = _mm256_add_epi32(real_count, cmp_result);

//                     X_N = _mm256_add_ps(_mm256_sub_ps(X2, Y2), X0);
//                     Y_N = _mm256_add_ps(_mm256_add_ps(XY, XY), Y0);
//                 }

//                 _mm256_store_si256((__m256i*)counts, real_count);

//                 for (size_t i = 0; i < remaining; i++) 
//                 {
//                     size_t px_index = y * WIDTH + x + chunk + i;
//                     if (px_index < WIDTH * HEIGHT) 
//                     { 
//                         get_pixels(set, px_index, counts[i]);
//                     }
//                 }
//             }
//         }
//     }

//     free(x0); free(counts);
// }

// версия aligned
inline void mandelbrot_simd(MandelBrot_t* set)
{
    const int SIMD_WIDTH = 8; 
    float real_dx = D_X * set->scale;

    __m256 MaxRadius = _mm256_set1_ps (MAX_RADIUS);
    
    // #pragma omp parallel for schedule(dynamic)
    // #pragma omp parallel for schedule(dynamic, 8)
    // #pragma omp parallel for schedule(guided, 8)
    for (size_t y = 0; y < HEIGHT; y++)
    {
        float base_x0 = (-(HALF_WIDTH           ) * D_X + set->x_offset) * set->scale;
        float y0      = (((float)y - HALF_HEIGHT) * D_Y + set->y_offset) * set->scale;

        for (size_t x = 0; x < WIDTH; x += VECTOR_SIZE)
        {
            for (size_t chunk = 0; chunk < VECTOR_SIZE; chunk += SIMD_WIDTH)
            {
                size_t remaining = (chunk + SIMD_WIDTH > VECTOR_SIZE) ? 
                                    VECTOR_SIZE - chunk : SIMD_WIDTH;
                
                float x0[SIMD_WIDTH] __attribute__((aligned(32))) = {};

                for (size_t i = 0; i < remaining; i++) 
                {
                    x0[i] = base_x0 + (x + chunk + i) * real_dx;
                }   

                for (size_t i = remaining; i < SIMD_WIDTH; i++) 
                {
                    x0[i] = 0.0f;
                }

                __m256 X0 = _mm256_load_ps (x0);  // загружаем выровненные данные
                __m256 Y0 = _mm256_set1_ps (y0);

                __m256 X_N = X0;
                __m256 Y_N = Y0;

                __m256i real_count = _mm256_setzero_si256 ();

                int count = 0;
                while (count++ < MAX_ITERATIONS)
                {
                    __m256 X2 = _mm256_mul_ps (X_N, X_N);
                    __m256 Y2 = _mm256_mul_ps (Y_N, Y_N);
                    __m256 XY = _mm256_mul_ps (X_N, Y_N);
                    __m256 R2 = _mm256_add_ps (X2, Y2);

                    __m256 mask = _mm256_cmp_ps(R2, MaxRadius, _CMP_LE_OS);
                    if (_mm256_testz_ps(mask, mask)) break;

                    __m256i cmp_result = _mm256_castps_si256 (mask); // float ? int
                            cmp_result = _mm256_srli_epi32   (cmp_result, 31);
                            real_count = _mm256_add_epi32    (real_count, cmp_result);

                    X_N = _mm256_add_ps (_mm256_sub_ps(X2, Y2), X0);
                    Y_N = _mm256_add_ps (_mm256_add_ps(XY, XY), Y0);
                }

                alignas(32) uint32_t counts[SIMD_WIDTH];
                _mm256_store_si256 ((__m256i*)counts, real_count);

                for (size_t i = 0; i < remaining; i++) 
                {
                    size_t px_index = y * WIDTH + x + chunk + i;
                    if (px_index < WIDTH * HEIGHT) 
                    { 
                        get_pixels(set, px_index, counts[i]);
                    }
                }
            }
        }
    } 
}

void run_performance_test (MandelBrot_t* set, int mode_measure) 
{
    if (!set->calculate) 
    {
        printf("Error: Calculation function not initialized\n");
        return;
    }

    RUN_TEST(set->calculate, mode_measure);
}


int handle_keyboard_input (MandelBrot_t* set, sf::Event &event) 
{
    if (event.type == sf::Event::KeyPressed)
    {
        switch (event.key.code)
        {            
            case sf::Keyboard::Key::Right:  set->x_offset   += MOVE_SPEED;      break;

            case sf::Keyboard::Key::Left:   set->x_offset   -= MOVE_SPEED;      break;

            case sf::Keyboard::Key::Up:     set->y_offset   -= MOVE_SPEED;      break;

            case sf::Keyboard::Key::Down:   set->y_offset   += MOVE_SPEED;      break;

            case sf::Keyboard::Key::Hyphen: set->scale      *= ZOOM_IN_FACTOR;  break;

            case sf::Keyboard::Key::Equal:  set->scale      *= ZOOM_OUT_FACTOR; break;

            default:                                                            break;
        }
    }

    return 0;
}


void get_pixels (MandelBrot_t* set,  size_t index, size_t count)
{
    if (count >= MAX_ITERATIONS) 
        set->pixels_array[index] = set->color_palette[MAX_ITERATIONS-1];
    else 
        set->pixels_array[index] = set->color_palette[count];
}

void init_color_palette(MandelBrot_t* set) 
{
    set->color_palette = (uint32_t*) calloc (MAX_ITERATIONS, sizeof(uint32_t));
    
    for (size_t i = 0; i < MAX_ITERATIONS; i++) 
    {
        // нормализация значения в диапазон [0, 1]
        float normalized = (float)i / (float)MAX_ITERATIONS;
        
        // используем тригонометрические функции для плавных переходов
        float r = sinf(2 * M_PI * normalized * 1.5f + M_PI / 3) * 0.5f + 0.5f;
        float g = sinf(2 * M_PI * normalized * 2.0f + M_PI / 2) * 0.5f + 0.5f;
        float b = sinf(2 * M_PI * normalized * 3.0f + M_PI    ) * 0.5f + 0.5f;
        
        // преобразуем в цвет формата 0xAARRGGBB
        set->color_palette[i] = 
            (0xFF << 24) |                       // Альфа-канал
            ((uint32_t) (r * 255) << 16) |       // Красный
            ((uint32_t) (g * 255) << 8)  |       // Зеленый
             (uint32_t) (b * 255);               // Синий
    }
    
    // цвет для точек внутри множества (черный)
    set->color_palette[MAX_ITERATIONS-1] = DEFAULT_COLOR;
}

void get_fps (sf::Clock &clock, sf::Text &text, MandelBrot_t* set) 
{
    if (!set->calculate) 
    {
        text.setString("Algorithm not initialized");
        return;
    }

    uint64_t start = __rdtsc();
    

    set->calculate(set);

    uint64_t       end =       __rdtsc();
    sf::Time sfml_time = clock.restart();
    
    const double cpu_ghz = 2.5; 
    double rdtsc_fps     = 1.0 / ((end - start) / (cpu_ghz * 1e9));

    double sfml_fps = 1.0 / sfml_time.asSeconds();
    
    char buffer[64];
    snprintf(buffer, sizeof(buffer), "FPS: rdtsc %.2f | SFML %.2f", rdtsc_fps, sfml_fps);
    text.setString(buffer);// -fopt-info-vec - flag
}

void dtor_mandlebrot (MandelBrot_t* set)
{
    if (!set) 
    {
        fprintf (stderr, "\nset == NULL\n"); // TODO сдеалть макросом 
        return;
    } 

    if (set->color_palette != NULL) 
    {
        free(set->color_palette);
        set->color_palette = NULL;
    }

    set->scale    = 0.0f;
    set->x_offset = 0.0f;
    set->y_offset = 0.0f;

    if (set->pixels_array != NULL) 
    {
        free(set->pixels_array);

        set->pixels_array = NULL;  

        DBG (printf("Freed pixels array\n"); ) 
    }
}