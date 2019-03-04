CC = mpicc

CFLAGS = -g -Wall

LDFLAGS = -lm

OBJECTS = prog.o utils.o dense_orthogon.o

TARGET = prog

all: $(OBJECTS)
	$(CC) $(CFLAGS) -o $(TARGET) $^ $(LDFLAGS)

.PHONY: test clean

test: all
	./test.sh

clean:
	$(RM) $(OBJECTS) $(TARGET)
