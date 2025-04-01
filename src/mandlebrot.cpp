#include "../include/mandlebrot.h"

void init_mandelbrot (struct MandelBrot* set) 
{
    set->scale        =  3.0f;
    set->x_offset     = -0.1f;
    set->y_offset     =  0.f ;
    set->pixels_array = (uint32_t*) calloc (HEIGHT * WIDTH, sizeof(uint32_t));
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
            if (event.type == sf::Event::Closed) window.close ();

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


void mandelbrot_naive (struct MandelBrot* set)
{
    const float real_dx = set->scale * (1 / (float)WIDTH);
    // set->scale - масштаб (увеличение/уменьшение)
    // 1 / WIDTH — нормировка на ширину экрана (перевод из пикселей в координаты [-1.5, 1.5]
    // real_dx - шаг по оси x в комлексной плоскости

    for (size_t y = 0; y < HEIGHT; y++)
    {
        float x0 = (-((float) WIDTH  / 2)             * (1 / (float) WIDTH) + set->x_offset) * set->scale;
        float y0 = (((float)y - ((float) HEIGHT / 2)) * (1 / (float) WIDTH) + set->y_offset) * set->scale;

        for (size_t x = 0; x < WIDTH; x++, x0 += real_dx)
        {
            float x_n = x0;
            float y_n = y0;

            size_t count = 0;

            while (count++ < 256)
            {
                float x2 = x_n * x_n;
                float y2 = y_n * y_n;
                float xy = x_n * y_n;

                if (x2 + y2 >= 100.f) break;

                x_n = x2 - y2 + x0;
                y_n = xy + xy + y0;
            }

            get_pixels (set, y * WIDTH + x, count - 1);
        }
    }
}


void mandelbrot_vectorized (struct MandelBrot* set)
{
    const float real_dx = set->scale * (1 / (float)WIDTH);

    for (size_t y = 0; y < HEIGHT; y++) 
    {
            float x0 = (-((float) WIDTH  / 2)             * (1 / (float) WIDTH) + set->x_offset) * set->scale;
            float y0 = (((float)y - ((float) HEIGHT / 2)) * (1 / (float) WIDTH) + set->y_offset) * set->scale;

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

            while (count++ < 256)
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

                // проверяем условие |z_n|^2 >= 100
                int cmp[VECTOR_SIZE] = {};
                for (size_t i = 0; i < VECTOR_SIZE; i++)
                {
                    if (X2[i] + Y2[i] <= 100.f) 
                        cmp[i] = 1;                     // cmp[i] -> i-й бит mask
                }

                int mask = 0;
                for (size_t i = 0; i < VECTOR_SIZE; i++) 
                    mask |= (cmp[i] << i);

                if (!mask) break;

                // обновляем z_n 
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
}


void run_performance_test (struct MandelBrot* set) 
{
    switch (set->mode) 
    {
        case BY_PIXELS: 
        {
            RUN_TEST (mandelbrot_naive);  
            break;
        }

        case BY_VECTOR:
        {
            RUN_TEST (mandelbrot_vectorized); 
            break;
        }

        default:        
            break;
    }
}

inline void get_pixels (struct MandelBrot* set,  size_t index, size_t count)
{
    if (count != 256)
    {
        count *= 100;
        set->pixels_array[index] = 0xFF | count << 24 | count << 16 | count << 8;
    }

    else
        set->pixels_array[index] = 0xFF000000;
}


inline int handle_keyboard_input (struct MandelBrot* set, sf::Event &event) 
{
    if (event.type == sf::Event::KeyPressed)
    {
        switch (event.key.code)
        {            
            case sf::Keyboard::Key::Right:  set->x_offset   += 0.05f; break;

            case sf::Keyboard::Key::Left:   set->x_offset   -= 0.05f; break;

            case sf::Keyboard::Key::Up:     set->y_offset   -= 0.05f; break;

            case sf::Keyboard::Key::Down:   set->y_offset   += 0.05f; break;

            case sf::Keyboard::Key::Hyphen: set->scale      *= 1.25f; break;

            case sf::Keyboard::Key::Equal:  set->scale      *= 0.8f;  break;

            default:                                                  break;
        }
    }

    return 0;
}


inline void get_fps (sf::Clock &clock, sf::Text &text, struct MandelBrot* set) 
{
    switch (set->mode)
    {
        case BY_PIXELS: 
        {
            clock.restart();
            mandelbrot_naive (set);  
            break;
        }

        case BY_VECTOR:
        {
            clock.restart();
            mandelbrot_vectorized(set); 
            break;
        }    

        default:                                     
            break;
    }

    sf::Time elapsed_time = clock.getElapsedTime();

    char buffer[16] = {};

    sprintf (buffer, "FPS: %.2f", 1.f / elapsed_time.asSeconds ());

    text.setString (buffer);
}
