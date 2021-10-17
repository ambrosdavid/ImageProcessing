CC=gcc
CFLAGS=-Wall -lm
CFLAGS_MAIN=-Wall -ansi --pedantic -lm -g3 -O3 -fsanitize=address -fsanitize=undefined -std=gnu89 -Wextra

main_iplib : main_iplib.c bmp.o ip_lib.o 
	$(CC) $^ -o $@ $(CFLAGS_MAIN)
ip_lib.o : ip_lib.c
	$(CC) $^ -c -o $@ $(CFLAGS_MAIN)
bmp.o : bmp.c
	$(CC) $^ -c -o $@ $(CFLAGS)
test_mat : bmp.o ip_lib.o test_mat.c
	$(CC) $^ -o $@ $(CFLAGS)
mandelbrot : bmp.o ip_lib.o test_bmp.c
	$(CC) $^ -o $@ $(CFLAGS)
clean: 
	rm * main_iplib
	rm test_bmp test_mat