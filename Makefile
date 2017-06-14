CFLAGS = -O3 -I/usr/local/include/
LFLAGS = -L/usr/local/lib/
CC = g++
OBJ = lnklst.o cmpile.o cicada.o intrpt.o bytecd.o ciclib.o userfn.o ccmain.o align3d.o

align3d: $(OBJ)
	$(CC) $(LFLAGS) -o align3d $(OBJ) /usr/local/lib/libgsl.a /usr/local/lib/libgslcblas.a -lm -lgsl -lgslcblas
	rm *.o

lnklst.o: lnklst.h
cmpile.o: lnklst.h cmpile.h
cicada.o: lnklst.h cmpile.h cicada.h
intrpt.o: lnklst.h cmpile.h cicada.h intrpt.h bytecd.h userfn.h
bytecd.o: lnklst.h cmpile.h cicada.h intrpt.h bytecd.h ciclib.h userfn.h
ciclib.o: lnklst.h cmpile.h cicada.h intrpt.h bytecd.h ciclib.h userfn.h ccmain.h
userfn.o: lnklst.h intrpt.h userfn.h ccmain.h align3d.h
	$(CC) $(CFLAGS) -c -o userfn.o userfn.cpp
ccmain.o: lnklst.h cmpile.h cicada.h intrpt.h bytecd.h ciclib.h userfn.h ccmain.h
align3d.o: lnklst.h intrpt.h userfn.h align3d.h
	$(CC) $(CFLAGS) -c -o align3d.o align3d.cpp
