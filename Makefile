LIBS = -lsfml-graphics -lsfml-window -lsfml-system
SANIT_FLAGS = alignment,undefined
OPTIMISIE_FLAG1 = -O3
OPTIMISIE_FLAG2 = -O0 

CC = clang++
#clang++
# -fopt-info-vec
# -Rpass=vectorize

all: mandle

mandle: main.o mandlebrot.o
	@$(CC) -o mandle.exe obj/main.o obj/mandlebrot.o $(LIBS)

mandlebrot.o: src/mandlebrot.cpp     
	@$(CC) -c -mavx2 $(OPTIMISIE_FLAG1) src/mandlebrot.cpp -o obj/mandlebrot.o

main.o: src/main.cpp
	@$(CC) -c -mavx2 $(OPTIMISIE_FLAG1) src/main.cpp -o obj/main.o

clean:
	rm obj/* mandle.exe