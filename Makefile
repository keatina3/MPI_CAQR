CC = mpicc

CFLAGS = -g -Wall

LDFLAGS = -lm

OBJECTS = prog.o utils.o dense_orthogon.o

TARGET = prog

all: $(OBJECTS)
	$(CC) $(CFLAGS) -o $(TARGET) $^ $(LDFLAGS)

.PHONY: test clean

test: all
	mpiexec --oversubscribe -n 4 ./prog -b 50 -n 10000 -m 5000

clean:
	$(RM) $(OBJECTS) $(TARGET)
