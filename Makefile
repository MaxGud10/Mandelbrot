LIBS = -lsfml-graphics -lsfml-window -lsfml-system
SANIT_FLAGS = -fsanitize=address,alignment,undefined
OPTIMISIE_FLAG1 = -O3
OPTIMISIE_FLAG2 = -O0

all: mandle

mandle: obj/main.o obj/mandelbrot_core.o obj/mandelbrot_utils.o obj/mandelbrot_render.o
	@g++ -o mandle.exe $^ $(LIBS)

obj/mandelbrot_utils.o: src/mandelbrot_utils.cpp
	@g++ -c -mavx2 $(OPTIMISIE_FLAG1) $< -o $@

obj/mandelbrot_core.o: src/mandelbrot_core.cpp
	@g++ -c -mavx2 $(OPTIMISIE_FLAG1) $< -o $@

obj/mandelbrot_render.o: src/mandelbrot_render.cpp
	@g++ -c -mavx2 $(OPTIMISIE_FLAG1) $< -o $@

obj/main.o: src/main.cpp
	@g++ -c -mavx2 $(OPTIMISIE_FLAG1) $< -o $@

clean:
	rm obj/* mandle.exe

.PHONY: all clean