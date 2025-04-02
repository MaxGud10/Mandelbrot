#include "../include/mandlebrot.h"


int main(int argc, char* argv[])
{
    MandelBrot set = {};

    init_mandelbrot (&set);

    if (argc == 1)
    {
        printf ("You didn't specify the mode\n");
        return -1;
    }

    else if (argc > 1) 
    {
        sscanf (argv[1], "%d", (int*)&set.mode);

        if (set.mode > 3 || set.mode < 1)
        {
            printf ("You specified the wrong mode of rendering\n");
            return -2;
        }

        if (argc == 2) get_mandel_brot_set (&set);

        else 
        {
            int mode_measure = 0;
            sscanf (argv[2], "%d", &mode_measure);

            if (mode_measure < 0 || mode_measure > 1) 
            {
                printf ("You specified the wrong mode of measure\n");
                return -3;
            }

            run_performance_test (&set, mode_measure);
        }
    }

    free(set.pixels_array);

    return 0;
}