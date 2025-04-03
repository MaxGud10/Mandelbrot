#include "../include/mandlebrot.h"

// Реализации функций:
// - get_mandel_brot_set
// - get_pixels
// - get_fps


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

        DBG( printf("Calculating FPS...\n"); )
        get_fps (clock, text, set);

        texture.update ((sf::Uint8*) set->pixels_array, WIDTH, HEIGHT, 0, 0);
        window.clear   ();
        window.draw    (sprite);
        window.draw    (text);
        window.display ();
    }
}


inline void get_pixels (MandelBrot_t* set,  size_t index, size_t count)
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

inline void get_fps (sf::Clock &clock, sf::Text &text, MandelBrot_t* set) 
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