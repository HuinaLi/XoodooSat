EXTENSION = cpp

CC = g++

EXE = xoodoo

OBJ = obj

LDFLAGS = -lcryptominisat5 -lm4ri

CFLAGS = -DLARGEMEM=ON
CFLAGS += -g

CXX_SOURCES = main_xoo.$(EXTENSION) xoodooRound.$(EXTENSION) gen_xoodoo_Chi_cnf.$(EXTENSION) gen_xoodoo_AS_cnf.$(EXTENSION)
CXX_HEADERS = xoodooRound.h gen_xoodoo_Chi_cnf.h gen_xoodoo_AS_cnf.h
CXX_OBJECTS = $(patsubst  %.$(EXTENSION), $(OBJ)/%.o, $(CXX_SOURCES))

.PHONY: XoodooSAT
all: $(OBJ) $(EXE)

$(OBJ):
	mkdir -p $(OBJ)

$(EXE): $(CXX_OBJECTS)
	$(CC) $(CXX_OBJECTS) -o $(EXE) $(LDFLAGS) $(CFLAGS)

$(OBJ)/main_xoo.o: $(CXX_SOURCES) $(CXX_HEADERS)
	$(CC) -c $< -o $@ $(LDFLAGS) $(CFLAGS)

$(OBJ)/xoodooRound.o: xoodooRound.$(EXTENSION) gen_xoodoo_Chi_cnf.$(EXTENSION) gen_xoodoo_AS_cnf.$(EXTENSION) $(CXX_HEADERS)
	$(CC) -c $< -o $@ $(LDFLAGS) $(CFLAGS)

$(OBJ)/gen_xoodoo_Chi_cnf.o: gen_xoodoo_Chi_cnf.$(EXTENSION) gen_xoodoo_Chi_cnf.h
	$(CC) -c $< -o $@ $(LDFLAGS) $(CFLAGS)

$(OBJ)/gen_xoodoo_AS_cnf.o: gen_xoodoo_AS_cnf.$(EXTENSION) gen_xoodoo_AS_cnf.h
	$(CC) -c $< -o $@ $(LDFLAGS) $(CFLAGS)

clean:
	rm -rf $(CXX_OBJECTS) $(EXE)
	rm CNF*.txt