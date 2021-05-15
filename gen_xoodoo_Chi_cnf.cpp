#include "gen_xoodoo_Chi_cnf.h"

void gen_xoodoo_Chi_cnf(vector<vector<string>>& cnf_clauses) {
    vector<string> aClause;
    //14cnf to describe Chi
    // A+B'+C'+E'+F'
	//input 3bit, output 3bit, 6bit in all
    aClause.push_back("F");
    aClause.push_back("T");
    aClause.push_back("T");
    aClause.push_back("0");
    aClause.push_back("T");
    aClause.push_back("T");
    cnf_clauses.push_back(aClause);
    aClause.clear();


 //A'+B+C'+D'+F'
    aClause.push_back("T");
    aClause.push_back("F");
    aClause.push_back("T");
    aClause.push_back("T");
    aClause.push_back("0");
    aClause.push_back("T");
    cnf_clauses.push_back(aClause);
    aClause.clear();


// A'+B'+C+D'+E'
    aClause.push_back("T");
    aClause.push_back("T");
    aClause.push_back("F");
    aClause.push_back("T");
    aClause.push_back("T");
    aClause.push_back("0");
    cnf_clauses.push_back(aClause);
    aClause.clear();


// B'+C'+D+E'+F'
    aClause.push_back("0");
    aClause.push_back("T");
    aClause.push_back("T");
    aClause.push_back("F");
    aClause.push_back("T");
    aClause.push_back("T");
    cnf_clauses.push_back(aClause);
    aClause.clear();


// A'+C'+D'+E+F'
    aClause.push_back("T");
    aClause.push_back("0");
    aClause.push_back("T");
    aClause.push_back("T");
    aClause.push_back("F");
    aClause.push_back("T");
    cnf_clauses.push_back(aClause);
    aClause.clear();


// A'+B'+D'+E'+F
    aClause.push_back("T");
    aClause.push_back("T");
    aClause.push_back("0");
    aClause.push_back("T");
    aClause.push_back("T");
    aClause.push_back("F");
    cnf_clauses.push_back(aClause);
    aClause.clear();


// A+C+E+F'
    aClause.push_back("F");
    aClause.push_back("0");
    aClause.push_back("F");
    aClause.push_back("0");
    aClause.push_back("F");
    aClause.push_back("T");
    cnf_clauses.push_back(aClause);
    aClause.clear();


// B+D+E'+F
    aClause.push_back("0");
    aClause.push_back("F");
    aClause.push_back("0");
    aClause.push_back("F");
    aClause.push_back("T");
    aClause.push_back("F");
    cnf_clauses.push_back(aClause);
    aClause.clear();


// A+B+D'+F
    aClause.push_back("F");
    aClause.push_back("F");
    aClause.push_back("0");
    aClause.push_back("T");
    aClause.push_back("0");
    aClause.push_back("F");
    cnf_clauses.push_back(aClause);
    aClause.clear();


// A'+C+D+E
    aClause.push_back("T");
    aClause.push_back("0");
    aClause.push_back("F");
    aClause.push_back("F");
    aClause.push_back("F");
    aClause.push_back("0");
    cnf_clauses.push_back(aClause);
    aClause.clear();


// A+B'+E+F
    aClause.push_back("F");
    aClause.push_back("T");
    aClause.push_back("0");
    aClause.push_back("0");
    aClause.push_back("F");
    aClause.push_back("F");
    cnf_clauses.push_back(aClause);
    aClause.clear();


// C'+D+E+F
    aClause.push_back("0");
    aClause.push_back("0");
    aClause.push_back("T");
    aClause.push_back("F");
    aClause.push_back("F");
    aClause.push_back("F");
    cnf_clauses.push_back(aClause);
    aClause.clear();


// B+C+D+F'
    aClause.push_back("0");
    aClause.push_back("F");
    aClause.push_back("F");
    aClause.push_back("F");
    aClause.push_back("0");
    aClause.push_back("T");
    cnf_clauses.push_back(aClause);
    aClause.clear();


// A+B+C+F'
    aClause.push_back("F");
    aClause.push_back("F");
    aClause.push_back("F");
    aClause.push_back("0");
    aClause.push_back("0");
    aClause.push_back("T");
    cnf_clauses.push_back(aClause);
    aClause.clear();
    
}



void check_xoodoo_Chi_cnf() {
    SATSolver solver;
    solver.set_num_threads(4);
    solver.new_vars(6);
    
    vector<vector<string>> Chi_Relation;
    gen_xoodoo_Chi_cnf(Chi_Relation);
    
    
    for (int i=0; i<1; i++) {
        
        for (int ithClause=0; ithClause<Chi_Relation.size(); ithClause++) {
            vector<Lit> clause;
            for (int ithVar=0; ithVar<6; ithVar++) {
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
            <<" -> "
            << (solver.get_model()[3] == l_False? 0:1)
            << ", " << (solver.get_model()[4] == l_False? 0:1)
            << ", " << (solver.get_model()[5] == l_False? 0:1)
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

