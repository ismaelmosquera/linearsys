#
# Makefile to build test linearsys
#
# Macro to perform cleanup
#
define cleanup
	mv prog bin
	rm *.o
endef
#
vpath %.h ../include
CC := gcc
FLAGS := -I ../include
CCF := -L../lib -llinearsys
PROG := prog
OBJ := test_linearsys.o 
$(PROG): $(OBJ)
	$(CC) $(CCF) $^ -O2 -s -DNDEBUG -o $@ && $(cleanup)
test_linearsys.o: test_linearsys.c linearsys.h
	$(CC) $(FLAGS) -c $<

