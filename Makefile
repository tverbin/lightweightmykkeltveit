
CFLAGS+=-Ofast 
LFLAGS+=-lm
CC=clang

OBJ=decycle.o

all: decycle

decycle: $(OBJ)
	$(CC) $(CFLAGS) -o $@ $(OBJ) $(LFLAGS)

.PHONY: clean
clean:
	rm $(OBJ) decycle

