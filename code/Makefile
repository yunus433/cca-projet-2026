CC      := gcc
CFLAGS  := -O3 -march=native -Wall -Wextra -std=c11
INCLUDES:= -Ieuclidean-division -Igcd -Ipseudorem -Iresultant
LDFLAGS :=
LDLIBS  := -lflint -lgmp
TARGET  := a.out

SRC := \
	main.c \
	euclidean-division/euclidean.c \
	gcd/euclid.c \
	pseudorem/pseudorem.c \
	resultant/subresultant2.c \
	resultant/subresultant.c \
	resultant/resultant.c \
	resultant/sylvester.c

OBJ := $(SRC:.c=.o)

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CC) $(CFLAGS) $(OBJ) -o $@ $(LDFLAGS) $(LDLIBS)

%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

run: $(TARGET)
	./$(TARGET)

clean:
	rm -f $(OBJ) $(TARGET)

.PHONY: all run clean
