EXTENSION = cpp
HEADER = h

CC = g++

EXE = xoodoo
OBJ = obj
INCLUDE = -I .

LDFLAGS = -lcryptominisat5 -lm4ri
CFLAGS = -O3 -DLARGEMEM=ON -std=c++0x
CFLAGS += -g

CXX_SOURCES = $(wildcard *.$(EXTENSION))
CXX_HEADERS = $((wildcard *.$(HEADER))
CXX_OBJECTS = $(patsubst  %.$(EXTENSION), $(OBJ)/%.o, $(CXX_SOURCES))
CNF = $(wildcard CNF*.txt)

.PHONY: XoodooSAT
all: $(OBJ) $(EXE)

$(OBJ):
	mkdir -p $(OBJ)

$(EXE): $(OBJ) $(CXX_OBJECTS)
	$(CC) $(CXX_OBJECTS) -o $(EXE) $(LDFLAGS) $(CFLAGS)

-include $(OBJ)/$(CXX_OBJECTS:.o=.d)

$(OBJ)/%.o: %.cpp
	$(CC) $(INCLUDE) -c $(CFLAGS) $< -o $(OBJ)/$*.o
	$(CC) $(INCLUDE) -MM $(CFLAGS) $< > $(OBJ)/$*.d
	@mv -f $(OBJ)/$*.d $(OBJ)/$*.d.tmp
	@sed -e 's|.*:|$(OBJ)/$*.o:|' < $(OBJ)/$*.d.tmp > $(OBJ)/$*.d
	@sed -e 's/.*://' -e 's/\\$$//' < $(OBJ)/$*.d.tmp | fmt -1 | \
	sed -e 's/^ *//' -e 's/$$/:/' >> $(OBJ)/$*.d
	@rm -f $(OBJ)/$*.d.tmp

clean:
	rm -rf $(CXX_OBJECTS) $(EXE)
    ifneq ($(CNF),)
		rm $(CNF)
    else
		@echo no CNF file is generated
    endif