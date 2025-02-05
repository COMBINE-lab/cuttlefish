all: example

CC = g++
CFLAGS	= -Wall -O3 -m64 -static -Wl,--whole-archive -lpthread -Wl,--no-whole-archive -std=c++17
CLINK	= -lm -static -O3 -Wl,--whole-archive -lpthread -Wl,--no-whole-archive -std=c++17

example: example.o libs/refresh/compression/ext_lib/zstd/linux/libzstd.a
	-mkdir -p bin
	$(CC) $(CLINK) -o bin/$@ $^

example.o: example.cpp
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -rf bin
	rm example.o
