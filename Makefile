CC = mpicc

CFLAGS = -g -Wall

LDFLAGS = -lm

OBJECTS = prog.o utils.o dense_orthogon.o

TARGET = prog

all: $(OBJECTS)
	$(CC) $(CFLAGS) -o $(TARGET) $^ $(LDFLAGS)

.PHONY: clean

clean:
	$(RM) $(OBJECTS) $(TARGET)
