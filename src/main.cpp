#include "../include/mandlebrot.h"

int main(int argc, char* argv[])
{
    MandelBrot set = {};

    // //!!!!!
    // printf("\n[%s] Received %d argument(s):\n", __func__, argc);
    // for(int i = 0; i < argc; i++) 
    // {
    //     printf("  argv[%d] = '%s'\n", i, argv[i]);
    // }
    // //!!!!
    // printf("%s:%d:%s() - message\n", __FILE__, __LINE__, __func__);

    
    if (argc == 1)
    {
        printf("You didn't specify the mode\n"
               "Available modes:\n"
               "1 - Naive calculation\n"
               "2 - Vectorized calculation\n"
               "3 - SIMD calculation\n");

        return -1;
    }

    if (argc > 1) 
    {
        int parse_result = sscanf(argv[1], "%d", (int*)&set.mode);
        // printf("\n[%s] Parsed mode: %d (parse result: %d)\n", __func__, set.mode, parse_result);

        if (set.mode > 3 || set.mode < 1)
        {
            printf("Invalid rendering mode specified\n"
                   "Available modes: 1 (naive), 2 (vectorized), 3 (SIMD)\n");
            return -2;
        }
    }

    init_mandelbrot (&set);

    if (argc == 2) 
    {
        printf ("\n[%s] Start get_mandel_brot_set()...\n", __func__);
        get_mandel_brot_set(&set);

        printf("[%s] Selected mode: %d | calculate ptr: %p\n", 
            __func__ ,    set.mode,     (void*)set.calculate);
    }

    else if (argc == 3) 
    {
        int mode_measure = 0;
        sscanf(argv[2], "%d", &mode_measure);

        if (mode_measure < 0 || mode_measure > 1) 
        {
            printf("Invalid measurement mode specified\n");
            printf("Use 0 for seconds, 1 for CPU ticks\n");
            return -3;
        }

        run_performance_test(&set, mode_measure);
    }

    printf("[MAIN] Selected mode: %d | calculate ptr: %p\n", 
        set.mode, (void*)set.calculate);

    dtor_mandlebrot(&set);
    return 0;
}