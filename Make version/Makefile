CC=g++
CC_FLAGS=-g -Wall

TARGET=NN

SRC_DIR=src
BIN_DIR=bin

SOURCES=$(SRC_DIR)/*.cpp

DEL_FILE = rm -f
LIBS_DIRS = -I./include/
LIBS = $(LIBS_DIRS)

.PHONY: clean build
all: clean build

build: $(SOURCES)
	$(CC) $(CC_FLAGS) $(LIBS) $(SOURCES) -o $(BIN_DIR)/$(TARGET)

run:
	./$(BIN_DIR)/$(TARGET)

clean:
	$(DEL_FILE) bim/*
