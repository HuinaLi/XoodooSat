//
//  gen_xoodoo_AS_cnf.cpp
//  XOODOOSAT
//
//  Created by LiHuina on 2021/4/27.
//
#include "gen_xoodoo_AS_cnf.h"


void gen_xoodoo_AS_cnf(vector<vector<string>>& cnf_clauses) {
    // ABCP四个变量
    // F1 = (C'+P)(B'+P)(A'+P)(A+B+C+P');
    vector<string> aClause;
  
    // (C'+P)
    aClause.push_back("0");
    aClause.push_back("0");
    aClause.push_back("T");
    aClause.push_back("F");
    cnf_clauses.push_back(aClause);
    aClause.clear();


    // (B'+P)
    aClause.push_back("0");
    aClause.push_back("T");
    aClause.push_back("0");
    aClause.push_back("F");
    cnf_clauses.push_back(aClause);
    aClause.clear();


    // (A'+P)
    aClause.push_back("T");
    aClause.push_back("0");
    aClause.push_back("0");
    aClause.push_back("F");
    cnf_clauses.push_back(aClause);
    aClause.clear();


    // (A+B+C+P')
    aClause.push_back("F");
    aClause.push_back("F");
    aClause.push_back("F");
    aClause.push_back("T");
    cnf_clauses.push_back(aClause);
    aClause.clear();
}

void check_xoodoo_AS_cnf() {
    int X=4, Y=3, Z=32;
    int one_round_num=X*Y*Z*2, var_num=X*Y*Z, AS_layer_num=X*Z, layer_num=4;

    SATSolver solver;
    solver.set_num_threads(4);//don't know why
    solver.new_vars(3*Y+3);
    
    vector<vector<string>> AS_Relation;
    gen_xoodoo_AS_cnf(AS_Relation);

    /*for (int ithClause=0; ithClause<AS_Relation.size(); ithClause++) {
        vector<Lit> clause;
        for (int ithVar=0; ithVar<4; ithVar++) {
            if (AS_Relation[ithClause][ithVar] != "0") {
                if (AS_Relation[ithClause][ithVar] == "T") {
                    clause.push_back(Lit((ithVar), true));
                }
                else{
                    clause.push_back(Lit((ithVar), false));
                } 
            } 
        }
        solver.add_clause(clause);
    }*/
    

    for (int x=0; x<1; x++) {
        for (int z=0; z<1; z++) {
            // 对每一个S盒check
            for (int ithClause=0; ithClause<AS_Relation.size(); ithClause++) {
                vector<Lit> clause[3];
                for (int ithVar=0; ithVar<Y; ithVar++) {
                    if (AS_Relation[ithClause][ithVar] != "0") {
                        if (AS_Relation[ithClause][ithVar] == "T") {
                            // cout<<"-"<<(ithVar*X*Z + z + Z*x) + 1<<" ";
                            clause[0].push_back(Lit((ithVar), true));
                            clause[1].push_back(Lit((Y+1 + ithVar), true));
                            clause[2].push_back(Lit(((Y+1)*2 + ithVar), true));
                        }
                        else{
                            // cout<<(ithVar*Z + z + Z*x) + 1<<" ";
                            clause[0].push_back(Lit((ithVar), false));
                            clause[1].push_back(Lit((Y+1 + ithVar), false));
                            clause[2].push_back(Lit(((Y+1)*2 + ithVar), false));
                        }
                    }
                }
                
                if (AS_Relation[ithClause][Y] != "0") {
                    if (AS_Relation[ithClause][Y] == "T") {
                        // cout<<"-"<<(var_num*layer_num + z + Z*x) + 1<<" ";
                        clause[0].push_back(Lit((Y), true));
                        clause[1].push_back(Lit((Y*2+1), true));
                        clause[2].push_back(Lit((Y*3+2), true));
                    }
                    else{
                        // cout<<(var_num*layer_num + z + Z*x) + 1<<" ";
                        clause[0].push_back(Lit((Y), false));
                        clause[1].push_back(Lit((2*Y+1), false));
                        clause[2].push_back(Lit((3*Y+2), false));
                    }
                }
                // cout<<"0"<<endl;
                solver.add_clause(clause[0]);
                solver.add_clause(clause[1]);
                solver.add_clause(clause[2]);
                clause[0].clear();
                clause[1].clear();
                clause[2].clear();
            }
        }
    }
    
    int solution_count = 0;
    lbool ret = solver.solve();
    while(ret == l_True) {
        // Use solution here. print it, for example.
        if (ret == l_True) {
            /*std::cout
            << solution_count <<" solution: "
            << (solver.get_model()[0] == l_False? 0:1)
            << ", " << (solver.get_model()[1] == l_False? 0:1)
            << ", " << (solver.get_model()[2] == l_False? 0:1)
            << ", " << (solver.get_model()[3] == l_False? 0:1)
            << std::endl;*/
            solution_count++;
            std::cout<<"solution AS: "<<solution_count<<std::endl;
        }
        else {
            // return false;
            cout<<"no"<<endl;
        }
            // Banning found solution
            vector<Lit> ban_solution;
            for (uint32_t var = 0; var < solver.nVars(); var++) {
                if (solver.get_model()[var] != l_Undef) {
                    ban_solution.push_back(
                    Lit(var, (solver.get_model()[var] == l_True)? true : false));
            }
        }
        solver.add_clause(ban_solution);
        ret = solver.solve();
    }
}
