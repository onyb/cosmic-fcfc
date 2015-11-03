CC=nvcc
CFLAGS=-I.
DEPS = utils.hpp
OBJ = serial.o utils.o

%.o: %.cu $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

grid_stride: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)

.PHONY: clean

clean:
	rm -f serial *.o *~ core
