TARGETS = $(basename $(wildcard *.c))
INCLUDE_PATH = $(YLAB)/ycode/common/matlab

all : $(TARGETS)

% : %.c $(INCLUDE_PATH)/util.h $(INCLUDE_PATH)/numeric.h
	mex -I$(INCLUDE_PATH) -largeArrayDims $<
