#include "../include/mandlebrot.h"


void init_mandelbrot (struct MandelBrot* set) 
{
    set->scale        =  INITIAL_SCALE;
    set->x_offset     =  INITIAL_X_OFFSET;
    set->y_offset     =  INITIAL_Y_OFFSET ;
    set->pixels_array = (uint32_t*) calloc (HEIGHT * WIDTH, sizeof(uint32_t));

    DBG( printf("Allocated pixels array: %p\n", (void*)set->pixels_array);           )
    DBG( printf("Memory allocated: %lu bytes\n", HEIGHT * WIDTH * sizeof(uint32_t)); )
}

void get_mandel_brot_set (struct MandelBrot* set)
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

        DBG( printf("Calculating FPS...\n"); )
        get_fps (clock, text, set);

        texture.update ((sf::Uint8*) set->pixels_array, WIDTH, HEIGHT, 0, 0);
        window.clear   ();
        window.draw    (sprite);
        window.draw    (text);
        window.display ();
    }
}


void mandelbrot_naive (struct MandelBrot* set)
{
    DBG( printf("Running naive Mandelbrot algorithm\n"); )

    const float real_dx = set->scale * D_X;
    DBG( printf("Real dx: %f\n", real_dx); )
    // set->scale - масштаб (увеличение/уменьшение)
    // 1 / WIDTH Ч нормировка на ширину экрана (перевод из пикселей в координаты [-1.5, 1.5]
    // real_dx - шаг по оси x в комлексной плоскости

    for (size_t y = 0; y < HEIGHT; y++)
    {
        float x0 = (-(HALF_WIDTH)            * D_X + set->x_offset) * set->scale;
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

                if (x2 + y2 >= MAX_RADIUS) 
                {
                    if (x == 0 && y == 0) 
                        printf("First point diverged at iteration %zu\n", count); 

                    break;
                }

                x_n = x2 - y2 + x0;
                y_n = xy + xy + y0;
            }

            get_pixels (set, y * WIDTH + x, count - 1);
        }
    }
    
    DBG( printf("Naive algorithm completed\n"); )
}


void mandelbrot_vectorized (struct MandelBrot* set)
{
    DBG( printf("Running vectorized Mandelbrot algorithm\n");)

    const float real_dx = set->scale * D_X;

    for (size_t y = 0; y < HEIGHT; y++) 
    {
        float x0 = (-(HALF_WIDTH)            * D_X + set->x_offset) * set->scale;
        float y0 = (((float)y - HALF_HEIGHT) * D_Y + set->y_offset) * set->scale;

        for (size_t x = 0; x < WIDTH; x += VECTOR_SIZE, x0 += VECTOR_SIZE * real_dx) 
        {
            float X0[VECTOR_SIZE] = {x0,               x0 +     real_dx, 
                                     x0 + 2 * real_dx, x0 + 3 * real_dx,
                                     x0 + 4 * real_dx, x0 + 5 * real_dx, 
                                     x0 + 6 * real_dx, x0 + 7 * real_dx};

            float X_N[VECTOR_SIZE] = {};
            float Y_N[VECTOR_SIZE] = {};

            for (size_t i = 0; i < VECTOR_SIZE; i++)
            {
                X_N[i] = X0[i]; // z0 = 0 => X_N = 0 (z1 = z0^2 + c = c => x1 = x0)
                Y_N[i] = y0;    // аналогично Y_N = 0
            }
               
            int count = 0;
            int real_count[VECTOR_SIZE] = {};

            while (count++ < MAX_ITERATIONS)
            {
                float X2[VECTOR_SIZE] = {}; 
                float Y2[VECTOR_SIZE] = {}; 
                float XY[VECTOR_SIZE] = {};

                for (size_t i = 0; i < VECTOR_SIZE; i++)
                { 
                    X2[i] = X_N[i] * X_N[i];
                    Y2[i] = Y_N[i] * Y_N[i];
                    XY[i] = X_N[i] * Y_N[i];
                }

                // провер€ем условие |z_n|^2 >= 100
                int cmp[VECTOR_SIZE] = {};
                for (size_t i = 0; i < VECTOR_SIZE; i++)
                {
                    if (X2[i] + Y2[i] <= MAX_RADIUS) 
                        cmp[i] = 1;                     // cmp[i] -> i-й бит mask
                }

                int mask = 0;
                for (size_t i = 0; i < VECTOR_SIZE; i++) 
                    mask |= (cmp[i] << i);

                if (!mask) 
                {
                    if (x == 0 && y == 0) 
                        DBG( printf("First block converged after %d iterations\n", count); )

                    break;
                }

                // обновл€ем z_n 
                for (size_t i = 0; i < VECTOR_SIZE; i++)
                { 
                    real_count[i] += cmp[i];

                    // z_{n+1} = z_n^2 + c
                    X_N[i] = X2[i] - Y2[i] + X0[i];
                    Y_N[i] = XY[i] + XY[i] + y0;
                }
            }

                for (size_t i = 0; i < VECTOR_SIZE; i++) 
                    get_pixels (set, y * WIDTH + x + i, real_count[i]);
        }
    }

    DBG (printf("Vectorized algorithm completed\n"); )
}


void mandelbrot_simd (struct MandelBrot* set) 
{
    DBG( printf("Running SIMD Mandelbrot algorithm\n"); )

    float real_dx = D_X * set->scale;

    __m256 MaxRadius = _mm256_set1_ps (MAX_RADIUS); // иницилизируем 8 элементов по 100

    __m256 DX = _mm256_set1_ps (real_dx);                       
    __m256 MUL_OFFSET = _mm256_set_ps (7, 6, 5, 4, 3, 2, 1, 0);
    DX = _mm256_mul_ps (DX, MUL_OFFSET);      // [7*dx, 6*dx, ..., 1*dx, 0*dx] - вычисл€ем начальные координаты дл€ 8 точек

    for (size_t y = 0; y < HEIGHT; y++)
    {
        float x_0 = (-(HALF_WIDTH)            * D_X + set->x_offset) * set->scale;
        float y_0 = (((float)y - HALF_HEIGHT) * D_Y + set->y_offset) * set->scale;
        
        for (size_t x = 0; x < WIDTH; x += VECTOR_SIZE, x_0 += VECTOR_SIZE * real_dx)
        {
            __m256 X0 = _mm256_set1_ps (x_0);    // [x0+7*dx, x0+6*dx, ..., x0+0*dx]
                   X0 = _mm256_add_ps  (X0, DX);

            __m256 Y0 = _mm256_set1_ps(y_0);

            __m256 X_N = X0;
            __m256 Y_N = Y0;

            int count = 0;
            __m256i real_count = _mm256_setzero_si256();

            while (count++ < MAX_ITERATIONS)
            {
                __m256 X2 = _mm256_mul_ps (X_N, X_N); // [x_n[0]?, ..., x_n[7]?]
                __m256 Y2 = _mm256_mul_ps (Y_N, Y_N);
                __m256 XY = _mm256_mul_ps (X_N, Y_N); // [x_n[0]*y_n[0], ..., x_n[7]*y_n[7]]

                __m256 R2 = _mm256_add_ps (X2, Y2);  // [x?+y?, ..., x?+y?]

                __m256 res = _mm256_cmp_ps (R2, MaxRadius, _CMP_LE_OS);

                if (!_mm256_movemask_ps (res)) break; // преобразуем векторную маску в 8-битное число (по одному биту на элемент)

                __m256i temp = _mm256_castps_si256 (res);       // интерпретируем float-маску как целые
                        temp = _mm256_srli_epi32   (temp, 31);  // cдвигаеv каждый 32-битный элемент на 31 бит, оставл€€ только старший бит (0 или 1).
                real_count   = _mm256_add_epi32    (real_count, temp);  // дл€ точек которые еще не вышли добавл€ем 1 к счетчику

                X_N = _mm256_sub_ps (X2, Y2 );
                X_N = _mm256_add_ps (X_N, X0); // (x? - y?) + x0

                Y_N = _mm256_add_ps (XY, XY );
                Y_N = _mm256_add_ps (Y_N, Y0); // 2xy + y0
            }

            uint32_t* counts = (uint32_t*) (&real_count); // преобразование AVX-регистра в массив
                                                          // [count[0], ..., count[7]]

            for (size_t i = 0; i < VECTOR_SIZE; i++)
                get_pixels (set, y * WIDTH + x + i, counts[i]); 
        } 
    }
}

void run_performance_test (struct MandelBrot* set, int mode_measure) 
{
    switch (set->mode) 
    {
        case BY_PIXELS: 
        {
            RUN_TEST(mandelbrot_naive, mode_measure);  
            break;
        }

        case BY_VECTOR:
        {
            RUN_TEST (mandelbrot_vectorized, mode_measure); 
            break;
        }

        case BY_SIMD:
        {   
            RUN_TEST (mandelbrot_simd, mode_measure);   
            break;
        }

        default:        
            break;
    }
}

inline void get_pixels (struct MandelBrot* set,  size_t index, size_t count)
{
    if (count != MAX_ITERATIONS)
        set->pixels_array[index] = 0xFFFF00FF - (count * 3) << 2; // TODO in const

    else
        set->pixels_array[index] = DEFAULT_COLOR;

    // TODO после отладки убрать эту проверку 
    if (index < 10) 
    {
        DBG( printf("Pixel %zu: count = %zu, color = 0x%08X\n", 
                     index,     count,       set->pixels_array[index]); )
    }
}


inline int handle_keyboard_input (struct MandelBrot* set, sf::Event &event) 
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


inline void get_fps (sf::Clock &clock, sf::Text &text, struct MandelBrot* set) 
{
    DBG( printf("Starting FPS calculation for mode: "); )
    switch (set->mode)
    {
        case BY_PIXELS: 
        {
            DBG( printf("Naive\n"); )

            clock.restart();
            mandelbrot_naive(set);
            break;
        }

        case BY_VECTOR:
        {
            DBG( printf("Vectorized\n"); )

            clock.restart();
            mandelbrot_vectorized(set);
            break;
        }    

        case BY_SIMD:
        { 
            DBG( printf("SIMD\n"); )

            clock.restart();
            mandelbrot_simd(set);
            break;
        }

        default:      
            printf("DEDLOH: Unknown\n");          

            break;
    }

    sf::Time elapsed_time = clock.getElapsedTime();

    char buffer[16] = {};

    sprintf (buffer, "FPS: %.2f", 1.f / elapsed_time.asSeconds ());

    text.setString (buffer);
}

void dtor_mandlebrot (struct MandelBrot* set)
{
    if (!set) 
    {
        fprintf (stderr, "\nset == NULL\n"); // TODO сдеалть макросом 
        return;
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