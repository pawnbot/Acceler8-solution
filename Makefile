CC=g++
CFLAGS=-c -O3
LDFLAGS=-O2 -pthread
SOURCES=main.c timer.c array_op.c kadane.c processor_op.c
EXECUTABLE=run
OBJS_NAMES=$(SOURCES:.c=.o)

all: $(SOURCES) $(EXECUTABLE)


$(EXECUTABLE): $(OBJS_NAMES)
	$(CC) $(LDFLAGS) $(OBJS_NAMES) -o $@

.c.o: *.h
	$(CC) $(CFLAGS) $< -o $@

clean:
	@rm -f $(OBJS_NAMES) $(EXECUTABLE)

