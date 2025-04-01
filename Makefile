LIBS = -lsfml-graphics -lsfml-window -lsfml-system
SANIT_FLAGS = -fsanitize=address,alignment,undefined
OPTIMISIE_FLAG1 = -O3
OPTIMISIE_FLAG2 = -O0

all: mandle

mandle: main.o mandlebrot.o
	@g++ -o mandle.exe obj/main.o obj/mandlebrot.o $(LIBS)

mandlebrot.o: src/mandlebrot.cpp
	@g++ -c -mavx2 $(OPTIMISIE_FLAG1) src/mandlebrot.cpp -o obj/mandlebrot.o

main.o: src/main.cpp
	@g++ -c -mavx2 $(OPTIMISIE_FLAG1) src/main.cpp -o obj/main.o

clean:
	rm obj/* mandle.exe
