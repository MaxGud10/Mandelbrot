#include "../include/mandlebrot.h"


void init_mandelbrot (MandelBrot_t* set) 
{
    set->scale        =  INITIAL_SCALE;
    set->x_offset     =  INITIAL_X_OFFSET;
    set->y_offset     =  INITIAL_Y_OFFSET ;
    set->pixels_array = (uint32_t*) calloc (HEIGHT * WIDTH, sizeof(uint32_t));

    //init_color_palette(set); // [x] ?? 

    switch(set->mode) 
    {
        case BY_PIXELS:  set->calculate = mandelbrot_naive;      /*printf ( "\n[%s] set->calculate = %p\n", __func__, set->calculate);*/ break;

        case BY_VECTOR:  set->calculate = mandelbrot_vectorized; /*printf ( "\n[%s] set->calculate = %p\n", __func__, set->calculate);*/  break;

        case BY_SIMD:    set->calculate = mandelbrot_simd;       break;

        default:         printf ("\nERROR: Invalid mode\n");     break;
    }

    // DBG( printf("Allocated pixels array: %p\n", (void*)set->pixels_array);           )
    // DBG( printf("Memory allocated: %lu bytes\n", HEIGHT * WIDTH * sizeof(uint32_t)); )
}

void get_mandel_brot_set (MandelBrot_t* set)
{
    // printf ("\nSTRART IN TO get_mandel_brot_set\n");
    // printf("[%s] Selected mode: %d | calculate ptr: %p\n", 
    //     __func__, set->mode, (void*)set->calculate);

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

        //DBG( printf("Calculating FPS...\n"); )
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
    //DBG( printf("Running naive Mandelbrot algorithm\n"); )

    const float real_dx = set->scale * D_X;
    //DBG( printf("Real dx: %f\n", real_dx); )
    // set->scale - ������� (����������/����������)
    // 1 / WIDTH � ���������� �� ������ ������ (������� �� �������� � ���������� [-1.5, 1.5]
    // real_dx - ��� �� ��� x � ���������� ���������

    // �������������� �� ������� (y)
    #pragma omp parallel for schedule(dynamic)
    for (size_t y = 0; y < HEIGHT; y++)
    {
        float x0 = (-(HALF_WIDTH)            * D_X + set->x_offset) * set->scale;  // ����������� DOUBLE
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
    //DBG( printf("Naive algorithm completed\n"); )
}

void mandelbrot_vectorized (MandelBrot_t* set)
{
    //DBG( printf("Running vectorized Mandelbrot algorithm\n");)

    const float real_dx = set->scale * D_X;

    #pragma omp parallel for schedule(dynamic)
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
                Y_N[i] = y0;    // ���������� Y_N = 0
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

                // ��������� ������� |z_n|^2 >= 100
                int cmp[VECTOR_SIZE] = {};
                for (size_t i = 0; i < VECTOR_SIZE; i++)
                {
                    if (X2[i] + Y2[i] <= MAX_RADIUS) 
                        cmp[i] = 1;                     // cmp[i] -> i-� ��� mask
                }

                int mask = 0;
                for (size_t i = 0; i < VECTOR_SIZE; i++) 
                    mask |= (cmp[i] << i);

                if (!mask) 
                {
                    if (x == 0 && y == 0) 
                        //DBG( printf("First block converged after %d iterations\n", count); )

                    break;
                }

                // ��������� z_n 
                for (size_t i = 0; i < VECTOR_SIZE; i++)
                { 
                    real_count[i] += cmp[i];

                    // z_{n+1} = z_n^2 + c
                    X_N[i] = X2[i] - Y2[i] + X0[i];
                    Y_N[i] = XY[i] + XY[i] + y0;
                }
            }

                // for (size_t i = 0; i < VECTOR_SIZE; i++) 
                //     get_pixels (set, y * WIDTH + x + i, real_count[i]);
                #pragma omp critical // [x] ??? �������� ����� �� ����� 
                {
                    for (size_t i = 0; i < VECTOR_SIZE; i++) 
                        get_pixels(set, y * WIDTH + x + i, real_count[i]);
                }
        }
    }

    //DBG (printf("Vectorized algorithm completed\n"); )
}
// TODO ���������� ��� ������������ (� ������ ����). ���������� ��� ������ ��� �� ������. ����������� ���������� calloc -> aligncalloc ��� �����, ��� �� ����� ��� ������� �� ��������
// -o3 -o0, on off Asan, 3 version, Vector_size - ��� ����������� �� 2^n �������, openMd, align memory, compilers, gcc 9-13-14, clang 10-19, double-float
// ������������� fps �
// [x]-S - ���� ������� ������� asm ���� (o0 and 03)




void mandelbrot_simd (MandelBrot_t* set) 
{
    DBG( printf("Running SIMD Mandelbrot algorithm\n"); )

    float real_dx = D_X * set->scale;

    __m256 MaxRadius = _mm256_set1_ps (MAX_RADIUS); // ������������� 8 ��������� �� 100

    __m256 DX = _mm256_set1_ps (real_dx);                       
    __m256 MUL_OFFSET = _mm256_set_ps (7, 6, 5, 4, 3, 2, 1, 0);
    DX = _mm256_mul_ps (DX, MUL_OFFSET);      // [7*dx, 6*dx, ..., 1*dx, 0*dx] - ��������� ��������� ���������� ��� 8 �����

    #pragma omp parallel for schedule(dynamic)
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

                if (!_mm256_movemask_ps (res)) break; // ����������� ��������� ����� � 8-������ ����� (�� ������ ���� �� �������)

                __m256i temp = _mm256_castps_si256 (res);       // �������������� float-����� ��� �����
                        temp = _mm256_srli_epi32   (temp, 31);  // c������v ������ 32-������ ������� �� 31 ���, �������� ������ ������� ��� (0 ��� 1).
                real_count   = _mm256_add_epi32    (real_count, temp);  // ��� ����� ������� ��� �� ����� ��������� 1 � ��������

                X_N = _mm256_sub_ps (X2, Y2 );
                X_N = _mm256_add_ps (X_N, X0); // (x? - y?) + x0

                Y_N = _mm256_add_ps (XY, XY );
                Y_N = _mm256_add_ps (Y_N, Y0); // 2xy + y0
            }

            uint32_t* counts = (uint32_t*) (&real_count); // �������������� AVX-�������� � ������
                                                          // [count[0], ..., count[7]]

            for (size_t i = 0; i < VECTOR_SIZE; i++)
                get_pixels (set, y * WIDTH + x + i, counts[i]); 
        } 
    }
}

void run_performance_test (MandelBrot_t* set, int mode_measure) 
{
    if (!set->calculate) {
        printf("Error: Calculation function not initialized\n");
        return;
    }

    RUN_TEST(set->calculate, mode_measure);
    // switch (set->mode) 
    // {
    //     case BY_PIXELS: 
    //     {
    //         RUN_TEST(mandelbrot_naive, mode_measure);  
    //         break;
    //     }

    //     case BY_VECTOR:
    //     {
    //         RUN_TEST (mandelbrot_vectorized, mode_measure); 
    //         break;
    //     }

    //     case BY_SIMD:
    //     {   
    //         RUN_TEST (mandelbrot_simd, mode_measure);   
    //         break;
    //     }

    //     default:        
    //         break;
    // }
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

void dtor_mandlebrot (MandelBrot_t* set)
{
    if (!set) 
    {
        fprintf (stderr, "\nset == NULL\n"); // TODO ������� �������� 
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


// inline void get_pixels (MandelBrot_t* set,  size_t index, size_t count)
void get_pixels (MandelBrot_t* set,  size_t index, size_t count)
{
    // if (count >= MAX_ITERATIONS) {
    //     set->pixels_array[index] = set->color_palette[MAX_ITERATIONS-1];
    // } else {
    //     set->pixels_array[index] = set->color_palette[count];
    // }
    if (count != MAX_ITERATIONS)
    {
        set->pixels_array[index] = (0xFFFF00FF - (count * 3)) << 2; // TODO 0xFFFF00FF in const
    }
// [x] ������� ��� ��������� ����� ����� sin 
// ������� ��������� ����� � ���������� ��� ����� (����� � �������). ������������ � ��������  
    else
        set->pixels_array[index] = DEFAULT_COLOR;
}

void init_color_palette(MandelBrot_t* set) 
{
    set->color_palette = (uint32_t*)malloc(MAX_ITERATIONS * sizeof(uint32_t));
    
    for (size_t i = 0; i < MAX_ITERATIONS; i++) {
        // ������������ �������� � �������� [0, 1]
        float normalized = (float)i / (float)MAX_ITERATIONS;
        
        // ���������� ������������������ ������� ��� ������� ���������
        float r = sinf (2 * M_PI * normalized * 1.5f + M_PI/3) * 0.5f + 0.5f;
        float g = sinf (2 * M_PI * normalized * 2.0f + M_PI/2) * 0.5f + 0.5f;
        float b = sinf (2 * M_PI * normalized * 3.0f + M_PI)   * 0.5f + 0.5f;
        
        // ����������� � ���� ������� 0xAARRGGBB
        set->color_palette[i] = 
            (0xFF << 24) |                      // �����-�����
            ((uint32_t)(r * 255) << 16) |       // �������
            ((uint32_t)(g * 255) << 8)  |       // �������
             (uint32_t)(b * 255);               // �����
    }
    
    // ���� ��� ����� ������ ��������� (������)
    set->color_palette[MAX_ITERATIONS-1] = DEFAULT_COLOR;
}

// inline void get_fps (sf::Clock &clock, sf::Text &text, MandelBrot_t* set) 
void get_fps (sf::Clock &clock, sf::Text &text, MandelBrot_t* set) 
{
    if (!set->calculate) 
    {
        text.setString("Algorithm not initialized");
        return;
    }

    uint64_t start = __rdtsc();
    

    //printf ("\nSET->CALCULATE = %p\n", set->calculate);
    set->calculate(set);
 
    // =================
    // static bool first_run = true;  
    
    // if (first_run) 
    // {
    //     printf("[FPS] First run | calculate function: %p\n", (void*)set->calculate);
    //     first_run = false;
    // }

    // if (!set->calculate) 
    // {
    //     text.setString("Algorithm not initialized");
    //     return;
    // }
    
    //====================

    uint64_t       end =       __rdtsc();
    sf::Time sfml_time = clock.restart();
    
    const double cpu_ghz = 2.5; 
    double rdtsc_fps     = 1.0 / ((end - start) / (cpu_ghz * 1e9));

    double sfml_fps = 1.0 / sfml_time.asSeconds();
    
    char buffer[64];
    snprintf(buffer, sizeof(buffer), "FPS: rdtsc %.2f | SFML %.2f", rdtsc_fps, sfml_fps);
    text.setString(buffer);
}


