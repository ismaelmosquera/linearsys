#
# Makefile to build linearsys.dll
#
# Macro to perform reordering and cleanup
#
define cleanup
	mv linearsys.dll lib
	rm *.o
endef
#
vpath %.h include
VPATH := src
CC := gcc
FLAGS := -I include
SLIB := linearsys.dll
OBJ := linearsys.o lu.o matrix.o vector.o numio.o qr.o eigen.o svd.o
$(SLIB): $(OBJ)
	$(CC) $^ -shared -lm -O2 -s -DNDEBUG -o $@ && $(cleanup)	
linearsys.o: linearsys.c linearsys.h lu.h matrix.h vector.h
	$(CC) $(FLAGS) -c $<
lu.o: lu.c lu.h matrix.h
	$(CC) $(FLAGS) -c $<
matrix.o: matrix.c matrix.h linearsys.h numio.h
	$(CC) $(FLAGS) -c $<
vector.o: vector.c vector.h numio.h
	$(CC) $(FLAGS) -c $<
numio.o: numio.c numio.h
	$(CC) $(FLAGS) -c $<
qr.o: qr.c qr.h
	$(CC) $(FLAGS) -c $<
eigen.o: eigen.c eigen.h
	$(CC) $(FLAGS) -c $<
svd.o: svd.c svd.h
	$(CC) $(FLAGS) -c $<

