default: test1

test1: main.c integrate.c integrate.h interval.c interval.h
	gcc -Wall -o test1 main.c integrate.c interval.c -lm

clean:
	rm -f test1
