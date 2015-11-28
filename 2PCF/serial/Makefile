CC=g++
CFLAGS=-I Wall
DEPS = utils.hpp
OBJ = serial.o utils.o

%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

serial: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)

.PHONY: clean

clean:
	rm -f serial *.o *~ core
