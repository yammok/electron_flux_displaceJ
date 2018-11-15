flux-exe: main.o load.o alloc.o init.o print.o open.o trans.o flux.o int.o util.o density.o
	gcc -Wall main.o  load.o alloc.o init.o print.o open.o trans.o flux.o int.o util.o density.o -o flux-exe -lm

clean:
	rm main.o load.o alloc.o init.o print.o open.o trans.o flux.o int.o util.o density.o
	rm flux-exe
