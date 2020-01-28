laguerre: main.o splines.o points.o laguerre.o horner.o gaus/libge.a
	$(CC) -o laguerre main.o splines.o points.o horner.o laguerre.o -L gaus -l ge 

aproks: main.o splines.o points.o aproksymator_na_bazie.o gaus/libge.a
	$(CC) -o aproks main.o splines.o points.o aproksymator_na_bazie.o -L gaus -l ge

interpol: main.o splines.o points.o interpolator.o gaus/libge.a
	$(CC) -o interpol main.o splines.o points.o interpolator.o -L gaus -l ge

prosta: main.0 splines.o points.o prosta.o
	$(CC) -o prosta main.o splines.o points.o prosta.o

aproksymator_na_bazie.o: makespl.h points.h gaus/piv_ge_solver.h
	$(CC) -I gaus -c aproksymator_na_bazie.c

laguerre.o: makespl.h points.h gaus/piv_ge_solver.h horner.h
	$(CC) -I gaus -c laguerre.c

horner.o: horner.h horner.c
	$(CC) -c horner.c horner.h

interpolator.o: makespl.h points.h gaus/piv_ge_solver.h
	$(CC) -I gaus -c interpolator.c

.PHONY: clean

clean:
	-rm *.o aproks interpol prosta horner








