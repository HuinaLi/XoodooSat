//  Generate the CNF clauses of Original_DDT(17 cnf cluases)
//  gen_xoodoo_Original_Chi_cnf.cpp
//  XOODOOSAT
//
//  Created by LiHuina on 2021/11/15.
//
#include "gen_Original_xoodoo_Chi_cnf.h"

void gen_Original_xoodoo_Chi_cnf(vector<vector<string>>& cnf_clauses){
    vector<string> aClause;
    //17cnf to describe Original_Chi
    // x2 + x1' +x0' + y1' + y0'
	//7bit in all
    aClause.push_back("F");
    aClause.push_back("T");
    aClause.push_back("T");
    aClause.push_back("0");
    aClause.push_back("T");
    aClause.push_back("T");
    aClause.push_back("0");
    cnf_clauses.push_back(aClause);
    aClause.clear();


 //x2' + x1 + x0' + y2' + y0'
    aClause.push_back("T");
    aClause.push_back("F");
    aClause.push_back("T");
    aClause.push_back("T");
    aClause.push_back("0");
    aClause.push_back("T");
    aClause.push_back("0");
    cnf_clauses.push_back(aClause);
    aClause.clear();


// x2' + x1' +x0 +y2' +y1'
    aClause.push_back("T");
    aClause.push_back("T");
    aClause.push_back("F");
    aClause.push_back("T");
    aClause.push_back("T");
    aClause.push_back("0");
    aClause.push_back("0");
    cnf_clauses.push_back(aClause);
    aClause.clear();


//x1' + x0' + y2 +y1' + y0'
    aClause.push_back("0");
    aClause.push_back("T");
    aClause.push_back("T");
    aClause.push_back("F");
    aClause.push_back("T");
    aClause.push_back("T");
    aClause.push_back("0");
    cnf_clauses.push_back(aClause);
    aClause.clear();


// x2' + x0' +y2'+y1+y0'
    aClause.push_back("T");
    aClause.push_back("0");
    aClause.push_back("T");
    aClause.push_back("T");
    aClause.push_back("F");
    aClause.push_back("T");
    aClause.push_back("0");
    cnf_clauses.push_back(aClause);
    aClause.clear();


// x2' + x1' +y2'+y1'+y0
    aClause.push_back("T");
    aClause.push_back("T");
    aClause.push_back("0");
    aClause.push_back("T");
    aClause.push_back("T");
    aClause.push_back("F");
     aClause.push_back("0");
    cnf_clauses.push_back(aClause);
    aClause.clear();


// x0'+p
    aClause.push_back("0");
    aClause.push_back("0");
    aClause.push_back("T");
    aClause.push_back("0");
    aClause.push_back("0");
    aClause.push_back("0");
    aClause.push_back("F");
    cnf_clauses.push_back(aClause);
    aClause.clear();


// x1'+p
    aClause.push_back("0");
    aClause.push_back("T");
    aClause.push_back("0");
    aClause.push_back("0");
    aClause.push_back("0");
    aClause.push_back("0");
    aClause.push_back("F");
    cnf_clauses.push_back(aClause);
    aClause.clear();


// x0+y2+y1+y0'
    aClause.push_back("0");
    aClause.push_back("0");
    aClause.push_back("F");
    aClause.push_back("F");
    aClause.push_back("F");
    aClause.push_back("T");
    aClause.push_back("0");
    cnf_clauses.push_back(aClause);
    aClause.clear();


// x1+y2+y1'+y0
    aClause.push_back("0");
    aClause.push_back("F");
    aClause.push_back("0");
    aClause.push_back("F");
    aClause.push_back("T");
    aClause.push_back("F");
    aClause.push_back("0");
    cnf_clauses.push_back(aClause);
    aClause.clear();


// x2+y2'+y1+y0
    aClause.push_back("F");
    aClause.push_back("0");
    aClause.push_back("0");
    aClause.push_back("T");
    aClause.push_back("F");
    aClause.push_back("F");
    aClause.push_back("0");
    cnf_clauses.push_back(aClause);
    aClause.clear();


// x2'+p
    aClause.push_back("T");
    aClause.push_back("0");
    aClause.push_back("0");
    aClause.push_back("0");
    aClause.push_back("0");
    aClause.push_back("0");
    aClause.push_back("F");
    cnf_clauses.push_back(aClause);
    aClause.clear();


// y2+y1+y0+p'
    aClause.push_back("0");
    aClause.push_back("0");
    aClause.push_back("0");
    aClause.push_back("F");
    aClause.push_back("F");
    aClause.push_back("F");
    aClause.push_back("T");
    cnf_clauses.push_back(aClause);
    aClause.clear();


// x1+x0+y2+y0'
    aClause.push_back("0");
    aClause.push_back("F");
    aClause.push_back("F");
    aClause.push_back("F");
    aClause.push_back("0");
    aClause.push_back("T");
    aClause.push_back("0");
    cnf_clauses.push_back(aClause);
    aClause.clear();


// x2+x0+y1+y0'
    aClause.push_back("F");
    aClause.push_back("0");
    aClause.push_back("F");
    aClause.push_back("0");
    aClause.push_back("F");
    aClause.push_back("T");
    aClause.push_back("0");
    cnf_clauses.push_back(aClause);
    aClause.clear();


// x2+x1+y1'+y0
    aClause.push_back("F");
    aClause.push_back("F");
    aClause.push_back("0");
    aClause.push_back("0");
    aClause.push_back("T");
    aClause.push_back("F");
    aClause.push_back("0");
    cnf_clauses.push_back(aClause);
    aClause.clear();


 // x2+x1+x0+y0'
    aClause.push_back("F");
    aClause.push_back("F");
    aClause.push_back("F");
    aClause.push_back("0");
    aClause.push_back("0");
    aClause.push_back("T");
    aClause.push_back("0");
    cnf_clauses.push_back(aClause);
    aClause.clear();   
    
}



void check_Original_xoodoo_Chi_cnf(){
    SATSolver solver;
    solver.set_num_threads(4);//don't know why
    solver.new_vars(7);
    
    vector<vector<string>> Chi_Relation;
    gen_xoodoo_Chi_cnf(Chi_Relation);
    
    
    for (int i=0; i<1; i++) {
        
        for (int ithClause=0; ithClause<Chi_Relation.size(); ithClause++) {
            vector<Lit> clause;
            for (int ithVar=0; ithVar<7; ithVar++) {
                if (Chi_Relation[ithClause][ithVar] != "0") {
                    if (Chi_Relation[ithClause][ithVar] == "T") {
                        clause.push_back(Lit((ithVar + 3*i), true));
                    }
                    else{
                        clause.push_back(Lit((ithVar + 3*i), false));
                    }
                    
                }
            }
            solver.add_clause(clause);
 
        }
        
    }

    
    
    int solution_count = 1;
    
    while(true) {
        lbool ret = solver.solve();
        if (ret != l_True) {
            assert(ret == l_False);
            //All solutions found.
            exit(0);
        }
        
        //Use solution here. print it, for example.
        if (ret == l_True) {
            std::cout
            << solution_count <<" solution: "
            << (solver.get_model()[0] == l_False? 0:1)
            << ", " << (solver.get_model()[1] == l_False? 0:1)
            << ", " << (solver.get_model()[2] == l_False? 0:1)
             << ", " << (solver.get_model()[3] == l_False? 0:1)
            << ", " << (solver.get_model()[4] == l_False? 0:1)
            << ", " << (solver.get_model()[5] == l_False? 0:1)
             << ", " << (solver.get_model()[6] == l_False? 0:1)
            << std::endl;
            solution_count++;
        }
        else {
            //        return false;
            cout<<"no"<<endl;
        }
        //Banning found solution
        vector<Lit> ban_solution;
        
        for (uint32_t var = 0; var < solver.nVars(); var++) {
            if (solver.get_model()[var] != l_Undef) {
                ban_solution.push_back(
                                       Lit(var, (solver.get_model()[var] == l_True)? true : false));
            }
        }
        solver.add_clause(ban_solution);
    }
}

