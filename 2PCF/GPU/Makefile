CC=nvcc
CFLAGS=-I.
DEPS = utils.hpp
OBJ = grid_stride.o utils.o

%.o: %.cu $(DEPS)
        $(CC) -c -o $@ $< $(CFLAGS)

grid_stride: $(OBJ)
        $(CC) -o $@ $^ $(CFLAGS)

.PHONY: clean

clean:
      	rm -f grid_stride *.o *~ core