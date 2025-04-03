#include "../include/mandlebrot.h"

void init_mandelbrot (MandelBrot_t* set) 
{
    set->scale        =  INITIAL_SCALE;
    set->x_offset     =  INITIAL_X_OFFSET;
    set->y_offset     =  INITIAL_Y_OFFSET ;
    set->pixels_array = (uint32_t*) calloc (HEIGHT * WIDTH, sizeof(uint32_t));

    DBG( printf("Allocated pixels array: %p\n", (void*)set->pixels_array);           )
    DBG( printf("Memory allocated: %lu bytes\n", HEIGHT * WIDTH * sizeof(uint32_t)); )
}

inline int handle_keyboard_input (MandelBrot_t* set, sf::Event &event) 
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