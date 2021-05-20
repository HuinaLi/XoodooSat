EXTENSION = cpp
CC = g++

EXE = xoodoo
OBJ = obj
INCLUDE = -I .

LDFLAGS = -lcryptominisat5 -lm4ri
CFLAGS = -O3 -DLARGEMEM=ON -std=c++0x
# CFLAGS += -g

CXX_SOURCES = $(wildcard *.$(EXTENSION))
CXX_OBJECTS = $(patsubst  %.$(EXTENSION), $(OBJ)/%.o, $(CXX_SOURCES))
CXX_INCLUDES = $(CXX_OBJECTS:.o=.d)
CNF = $(wildcard CNF*AS*.txt)

.PHONY: all clean
all: $(OBJ) $(EXE)

$(OBJ):
	mkdir -p $(OBJ)

$(EXE): $(OBJ) $(CXX_OBJECTS)
	$(CC) $(CXX_OBJECTS) -o $(EXE) $(LDFLAGS) $(CFLAGS)

-include $(CXX_INCLUDES)

$(OBJ)/%.o: %.cpp
	$(CC) $(INCLUDE) -c $(CFLAGS) $< -o $(OBJ)/$*.o
	$(CC) $(INCLUDE) -MM $(CFLAGS) $< > $(OBJ)/$*.d
	@mv -f $(OBJ)/$*.d $(OBJ)/$*.d.tmp
	@sed -e 's|.*:|$(OBJ)/$*.o:|' < $(OBJ)/$*.d.tmp > $(OBJ)/$*.d
	@sed -e 's/.*://' -e 's/\\$$//' < $(OBJ)/$*.d.tmp | fmt -1 | \
	sed -e 's/^ *//' -e 's/$$/:/' >> $(OBJ)/$*.d
	@rm -f $(OBJ)/$*.d.tmp

clean:
	rm -rf $(CXX_OBJECTS) $(CXX_INCLUDES) $(EXE)
    ifneq ($(CNF),)
		rm CNF*AS*.txt
    else
		@echo no CNF file is generated
    endif