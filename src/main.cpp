#include "../include/mandlebrot.h"


int main (int argc, char* argv[])
{
    struct MandelBrot set = {};

    init_mandelbrot (&set);

    if (argc == 1)
    {
        printf ("You didn't specify the mode\n");
        
        return -1;
    }

    else if (argc > 1) 
    {
        sscanf (argv[1], "%d", (int*) &set.mode);

        if (set.mode > 2 || set.mode < 0)
        {
            printf ("You specified the wrong mode\n");
            return -2;
        }

        if (argc == 2) 
            get_mandel_brot_set (&set);

        else
            run_performance_test (&set);
    }

    free(set.pixels_array);

    return 0;
}