#ifndef GEN_XOODOO_AS_CNF_H
#define GEN_XOODOO_AS_CNF_H

#include <string.h>
#include <map>
#include <assert.h>
#include <vector>
#include <stdio.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cryptominisat5/cryptominisat.h>
#include <m4ri/m4ri.h>
using namespace CMSat;
using namespace std;

void gen_xoodoo_AS_cnf(vector<vector<string>>& cnf_clauses);
void check_xoodoo_AS_cnf();

#endif
