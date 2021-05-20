//
//  KeccakRoundSAT.cpp
//  SATCollidingTrailSearch
//
//  Created by Yi Tu on 11/26/19.
//  Copyright © 2019 Liu Guozhen. All rights reserved.
//
#include <cstdlib>
#include <chrono>
#include <ctime>
#include "xoodooRound.h"
using namespace std::chrono;
using namespace XOODOOSAT;

XoodooRound::XoodooRound(int round, int weight, int thread, int mode) :solver() {
    AS_weight_num = weight;
    AS_state_num = round;
    AS_mode = mode;
    round_num = round;
    core_state_num = 2*(round_num-1);//trail a1,b1,a2,b2,...
    AS_var_num = X*Z*AS_state_num;
    AS_node_var_num = X*Z;
    round_var_num = var_num*2;
    objFilePath = default_obj_file_path;
    thread_num = thread;

    solver.set_num_threads(thread_num);//threads to use
    /*
    X=4,Y=3,Z=32
    for 3 round trail cores:
    (0-383) a1
    (384,767) b1
    (768,1151) a2
    (1152,1535) b2
    each AS has X*Z=128 varibles
    AS_a1 (1536, 1663)
    AS_a2 (1664, 1791)
    AS_b2 (1792, 1919)
     */
    solver.new_vars(var_num*core_state_num + AS_var_num);//add vars

	rho_plane[0] = 1;//which plane to rho
	rho_plane[1] = 2;
	theta_L1[0] = 1;//1st lshift in theta
	theta_L1[1] = 5;
	theta_L2[0] = 1;//2nd lshift in theta
	theta_L2[1] = 14;
	rhoW_L1[0] = 1;//1st lshift in rhoW
	rhoW_L1[1] = 0;
	rhoW_L2[0] = 0;//2nd lshift in rhoW
	rhoW_L2[1] = 11;
	rhoE_L1[0] = 0;//1st lshift in rhoE
	rhoE_L1[1] = 1;
	rhoE_L2[0] = 2;//2nd lshift in rhoE
	rhoE_L2[1] = 8;
}

void XoodooRound::theta(tXoodooState &A) {
    unsigned int x, y;
    vector<tXoodooLane> P(X, 0), E(X, 0);
    
    for(x=0; x<X; x++) {
        for(y=0; y<Y; y++)
            P[x] ^= A[indexXY(x, y)];
    }
    for(x=0; x<X; x++)
        E[x] = ROLxoo(P[(x+X-theta_L1[0])%X], theta_L1[1]) ^ ROLxoo(P[(x+X-theta_L2[0])%X], theta_L2[1]);
    for(x=0; x<X; x++)
        for(y=0; y<Y; y++)
            A[indexXY(x, y)] ^= E[x];
}

void XoodooRound::rhoW(tXoodooState &A) {
    unsigned int x, y;
    tXoodooState tempA(A);
    
    y = rho_plane[0];//A1
    for(x=0; x<X; x++)
        A[indexXY(x, y)] = ROLxoo(tempA[indexXY(x+X-rhoW_L1[0], y)], rhoW_L1[1]);//(1,0)
    y = rho_plane[1];//A2
    for(x=0; x<X; x++)
        A[indexXY(x, y)] = ROLxoo(tempA[indexXY(x+X-rhoW_L2[0], y)], rhoW_L2[1]);//(0,11)
}

void XoodooRound::rhoE(tXoodooState &A) {
    unsigned int x, y;
    tXoodooState tempA(A);
    
    y = rho_plane[0];//A1
    for(x=0; x<X; x++)
        A[indexXY(x, y)] = ROLxoo(tempA[indexXY(x+X-rhoE_L1[0], y)], rhoE_L1[1]);//(0,1)
    y = rho_plane[1];//A2
    for(x=0; x<X; x++)
        A[indexXY(x, y)] = ROLxoo(tempA[indexXY(x+X-rhoE_L2[0], y)], rhoE_L2[1]);//(2,8)
}

void XoodooRound::chi(tXoodooState &A) {
    unsigned int x, y;
    vector<tXoodooLane> B(Y,0);
    
    for(x=0; x<X; x++) {
        for(y=0; y<Y; y++)
            B[y] = A[indexXY(x, y)] ^ ((~A[indexXY(x, y+1)]) & A[indexXY(x, y+2)]);
        for(y=0; y<Y; y++)
            A[indexXY(x, y)] = B[y];
    }
}

tXoodooState XoodooRound::lambda(tXoodooState A) {
    rhoE(A);
    theta(A);
    rhoW(A);
    return A;
}

void XoodooRound::Bit2XooState(const Bitmap &A, tXoodooState &B) {//32*(4y+x) + z
    B.assign(X*Y, 0);
    unsigned int x, y, z;
    for(y = 0; y < Y; y++){
        for(x = 0; x < X; x++) {
            B[X*y+x] = 0;
        }
        for(z = 0; z < Z; z++){
            for(x = 0; x < X; x++) {
                B[X*y+x] |= ((tXoodooLane)(A[Z*X*y+Z*x+z]) << (z));
            }
        }
    }
}

void XoodooRound::XooState2Bit(const tXoodooState &A, Bitmap &B) {//32*(4y+x) + z
    B.assign(var_num,0);
    unsigned int x, y, z;
    for(y = 0; y < Y; y++){
        for(z = 0; z < Z; z++){
            for(x = 0; x < X; x++) {
                B[Z*X*y+Z*x+z] = ((A[X*y+x] >> (z)) & 0x1);
            }
        }
    }
}

void XoodooRound::Bit2State(const Bitmap &A, State &B) {
    B.clear();
    for(unsigned int i = 0; i < var_num; i++) {
        if(A[i] == 1) B.push_back(i);
    }
}

void XoodooRound::State2Bit(const State &A, Bitmap &B) {
    B.assign(var_num,0);
    for(int i = 0; i < A.size(); i++) {
        B[A[i]] = 1;
    }
}

void XoodooRound::Bit2Plane(const Bitmap &A, vector<tXoodooLane> &B) {//32*x + z
    B.assign(X,0);
    unsigned int x, z;
    
    for(x = 0; x < X; x++) {
        B[x] = 0;
    }
    for(z = 0; z < Z; z++){
        for(x = 0; x < X; x++) {
            B[x] |= ((tXoodooLane)(A[Z*x+z]) << (z));
        }
    }
}

void XoodooRound::Plane2Bit(const vector<tXoodooLane> &A, Bitmap &B) {//32*x + z
    unsigned int x, z;

    for(z = 0; z < Z; z++){
        for(x = 0; x < X; x++) {
            B[Z*x+z] = ((A[x] >> (z)) & 0x1);
        }
    }
}

void XoodooRound::State2XooState(const State &A, tXoodooState &B) {
    Bitmap BitA;//Bit
    State2Bit(A,BitA);
    Bit2XooState(BitA,B);
}

void XoodooRound::XooState2State(const tXoodooState &A, State &B) {
    Bitmap BitA(var_num, 0);
    XooState2Bit(A,BitA);
    Bit2State(BitA,B);
}


int XoodooRound::caculateXooStateWeight(const tXoodooState &A) {
    int weight = 0;
    for(int x=0; x<X; x++) {
        for(int z=0; z<Z; z++) {
            for(int y=0; y<Y; y++) {
                if(A[indexXY(x, y)]) {
                    weight+=1;
                    break;
                }
            }
        }
    }
    return weight;
}


State XoodooRound::ShiftXZ(const vector<int> &dx_dz, const State &A) {
    tXoodooState stateA(X*Y, 0);
    State2XooState(A,stateA);
    tXoodooState shiftA(stateA);

    for(int x=0; x<X; x++) {
        for(int y=0; y<Y; y++) {
            shiftA[indexXY(x, y)] = ROLxoo(stateA[indexXY(x+X-dx_dz[0], y)], dx_dz[1]);//shift (dx,dz)
        }
    }
    State res;
    XooState2State(shiftA,res);
    return res;
}

bool XoodooRound::StateEqualAfterShift(const State &A, const State &B) {//A, B is sorted
    if(A.size() != B.size()) return false;
    if(A == B) return true;
    
    for(int dx=0; dx<X; dx++) {
        for(int dz=0; dz<Z; dz++) {
            State shiftA = ShiftXZ({dx,dz},A);
            if(shiftA == B) return true;
        }
    }
    return false;
}

bool XoodooRound::isSmaller(const State &A, const State &B) {//state A<B; A, B is sorted
    int i = A.size()-1, j = B.size()-1;
    while(i>=0 && j>=0) {
        if(A[i] < B[j]) return true;
        else if(A[i] > B[j]) return false;
        else {
            i--;
            j--;
        }
    }
    return i<j;
}

vector<int> XoodooRound::genSmallestState(const State &A) {//A is sorted  
    vector<int> dx_dz = {0, 0};
    State smallest(A);
    for(int dx=0; dx<X; dx++) {
        for(int dz=0; dz<Z; dz++) {
            State shiftA = ShiftXZ({dx,dz},A);
            //after shift (dx,dz), get the smaller state
            if(isSmaller(shiftA,smallest)) {
                dx_dz = {dx,dz};
                smallest = shiftA;
            }
        }
    }
    
    return dx_dz;
}


void XoodooRound::display(ostream& fout, const State &A) {
    vector<unsigned int> tempA(X*Z,0);//State column map, 0<=tempA[i]<=3
    for(int i=0; i<A.size(); i++) {
        int y;
        y = A[i]/(X*Z);
        assert(y < Y);
        tempA[A[i]-y*X*Z] |= ((0x1) << y);
    }
    for(int x=0; x<X; x++) {
        for(int z=0; z<Z; z++) {
            if(tempA[Z*x+z] == 0) fout<<".";
            else fout<<tempA[Z*x+z];
        }
        fout<<endl;
    }
}

void XoodooRound::displayXooState(ostream& fout, const tXoodooState &A) {
    State tmp;
    XooState2State(A,tmp);
    display(fout,tmp);
}


void XoodooRound::gen_RhoW_T(map<unsigned int, unsigned int>& RhoW_index) {
    Bitmap BitA(var_num, 0);
    tXoodooState temp(X*Y, 0);
    unsigned int i, j;
    
    for(i = 0; i < BitA.size(); i++){
        BitA[i] = 1;
        
        Bit2XooState(BitA,temp);
        rhoW(temp);
        XooState2Bit(temp,BitA);
        
        for(j = 0; j < BitA.size(); j++){
            if(BitA[j] == 1){
                RhoW_index.insert(pair<unsigned int, unsigned int>(i, j));
            }
            BitA[j] = 0;
        }
    }
    return;
}

void XoodooRound::gen_RhoE_T(map<unsigned int, unsigned int>& RhoE_index, map<unsigned int, unsigned int>& inverse_RhoE_index) {
    Bitmap BitA(var_num, 0);
    tXoodooState temp(X*Y, 0);
    unsigned int i, j;
    
    for(i = 0; i < BitA.size(); i++){
        BitA[i] = 1;

        Bit2XooState(BitA,temp);
        rhoE(temp);
        XooState2Bit(temp,BitA);
        
        for(j = 0; j < BitA.size(); j++){
            if(BitA[j] == 1){
                RhoE_index.insert(pair<unsigned int, unsigned int>(i, j));
                inverse_RhoE_index.insert(pair<unsigned int, unsigned int>(j, i));
            }
            BitA[j] = 0;
        }
    }
    return;
}

void XoodooRound::gen_Theta_T(vector<vector<unsigned int>>& relation) {
    vector<vector<unsigned int>> column;
    //prepare X*Z columns
    for(unsigned int x = 0; x<X; x++) {
        for(unsigned int z = 0; z<Z; z++) {
            vector<unsigned int> aColumn = {};
            for (int y = 0; y<Y; y++) {
                aColumn.push_back( (x+y*X)*Z + z);
            }
            column.push_back(aColumn);
        }
    }
    
    for(int x = 0; x<X; x++) {
        for(int y = 0; y<Y; y++) {
            for(int z = 0; z<Z; z++) {
                vector<unsigned int> temp;
                temp.push_back(Z*(x+X*y)+z);
                temp.push_back(Z*(x+X*y)+z);
                for (int i = 0; i<Y; i++) {
                    temp.push_back( column[Z*( (x+X-theta_L1[0])%X )+ ((z+Z-theta_L1[1])%Z) ][i] );
                    temp.push_back( column[Z*( (x+X-theta_L2[0])%X )+ ((z+Z-theta_L2[1])%Z) ][i] ); 
                }
                relation.push_back(temp);
            }
        }
    }
}

void XoodooRound::gen_extend_AS_cnf_num(string &cnf_num) {
    cnf_num += to_string(AS_var_num) + " ";
    if(round_num == 2) {
        cnf_num += to_string(AS_node_var_num) + " ";
    }
}

void XoodooRound::ban_solution(const map<State,int> &A) {
    bool zero = true;
    for(auto iter=A.begin(); iter!=A.end();iter++) {
        zero = zero && (iter->first.size()==0);
    }

    for(int dx=0; dx<X; dx++) {
        for(int dz=0; dz<Z; dz++) {
            vector<Lit> ban_solutions;
            for(auto iter=A.begin(); iter!=A.end();iter++) {
                State shiftA = ShiftXZ({dx,dz},iter->first);
                Bitmap BitA(var_num,0);
                State2Bit(shiftA,BitA);
                for (uint32_t var = 0; var < BitA.size(); var++) {
                    ban_solutions.push_back(Lit(var + iter->second, (BitA[var] == 1)? true : false));
                }
            }
            solver.add_clause(ban_solutions);
            ban_solutions.clear();
            if(zero) return;
        }
    }
}

void XoodooRound::write_result(const string pathname, const int weight, const int solution_count, const States &solution) {
    ofstream out(pathname, ios::out|ios::app);
    out<<"weight "<<weight<<" solution "<< solution_count <<endl;

    for (int i=0; i<solution.size()-1; i++) {
        out<<"a"<< i+1 <<":"<<endl;
        display(out,solution[i]);
        out<<endl;
    }
    out<<"b"<< solution.size()-1 <<":"<<endl;
    display(out,solution[solution.size()-1]);
    out<<"\n"<<endl;
    out.close();
}

int XoodooRound::get_weight() {
    //get each State weight(a1,a2,...)
    vector<int> as_weight(AS_state_num,0);
    for (int i=0; i<AS_node_var_num; i++) {
        for(int j=0; j<AS_state_num; j++) {
            if (solver.get_model()[i + var_num*core_state_num + AS_node_var_num*j] == l_True) {
                as_weight[j]++;
            }
        }
    }
    //get total weight
    int weight = 0;
    if(round_num == 2) weight = as_weight[0]*2 + as_weight[1];
    else {
        for(int j=0; j<AS_state_num; j++) weight += as_weight[j];
    }
    return weight;
}


void XoodooRound::solve_and_output(const string pathname) {
    int counting = 0;
    map<int,int> solution_counts;//map<weight,solution_counts>

    while(true) {
        lbool ret = solver.solve();
        if (ret != l_True) {
            assert(ret == l_False);
            cout<<"reach end"<<endl;
            cout<<"counting: "<<counting<<endl;
            //All solutions found.
            exit(0);
        }
        
        //start processing solution
        //get weight
        int weight = get_weight();

        //get result in States format
        States solution(AS_state_num);//states to show(a1,a2,...,bn-1)
        for (unsigned int i=0; i<var_num; i++) {
            for(int j=0; j<AS_state_num; j++) {
                int offset = j*round_var_num;
                if(j == AS_state_num-1) offset -= var_num;
                if (solver.get_model()[i + offset] == l_True) {
                    solution[j].push_back(i);//State
                }
            }
        }
        //solution + 1
        auto iter = solution_counts.find(weight);
        if(iter != solution_counts.end()) {
            iter->second += 1;
        }
        else solution_counts[weight] = 1;

        //a new state, shift to smallest
        vector<int> dx_dz = genSmallestState(solution[0]);
        for (int i=0; i<AS_state_num; i++) {
            solution[i] = ShiftXZ(dx_dz,solution[i]);
        }
        
        //now extend result
        //TODO

        //writing final result
        if(round_num > 2) write_result(pathname,weight,solution_counts[weight],solution);
        if(round_num > 2 || counting%1000 == 0){cout<<counting<<": ";cout<<dx_dz[0]<<","<<dx_dz[1]<<" ";for(int i=0; i<solution[0].size(); i++)cout<<solution[0][i]<<" ";cout<<endl;}counting++;
        
        //ban found solution and shifted solutions
        map<State,int> banned;//<state,offset>
        for(int i=0; i<solution.size()-1; i++) {
            banned[solution[i]] = i*round_var_num;
        }
        ban_solution(banned);
    }
}


void XoodooRound::read2Vector(vector<vector<int>>& res, const string pathname) {
    ifstream infile;
    infile.open(pathname.data());
    assert(infile.is_open());//fail to open
    vector<int> suanz;
    string s;

    res.clear();
    while (getline(infile, s)) {
        istringstream is(s);
        double d;
        while (!is.eof()) {
            is >> d;
            suanz.push_back(d);
        }
        res.push_back(suanz);
        suanz.clear();
        s.clear();
    }
    infile.close();
}

bool XoodooRound::read2States(States& res, ifstream &infile, int mode) {
    string s;
    string start;
    if(mode == 0) {
        start = "ab";
    }
    else {
        start = "NE";
    }

    res.clear();
    while (getline(infile, s)) {
        if(s[0] == start[0] || s[0] == start[1]) {
            State temp;
            for(int x=0; x<X; x++) {
                s.clear();
                assert(getline(infile, s));
                for(int z=0; z<Z; z++) {
                    if(s[z] != '.') {
                        int column = s[z]-'0';
                        for(int y=0; y<Y; y++) {
                            if((column >> y) & 0x1) temp.push_back(x*Z + y*X*Z + z);
                        }
                    }
                }
            }
            sort(temp.begin(), temp.end());
            res.push_back(temp);
        }
        if(res.size() == AS_state_num && mode == 0) {
            return true;
        }
        if(res.size() == AS_state_num-1 && mode != 0) {
            return true;
        }
        s.clear();
    }
    return false;
}


void XoodooRound::check_trails(const string pathname, int i, int mode) {
    ifstream infile(pathname.data(), ios::in);
    assert(infile.is_open());
    States res;
    int res_size = (mode == 0)? AS_state_num:AS_state_num-1;
    int print_size = AS_state_num;
    
    for(int k=0; read2States(res,infile,mode) && k<=i; k++) {
        if(res.size() != res_size) {
            cout<<"wrong size: "<<res.size()<<endl;
            return;
        }
        if(k<i) continue;
        cout<<"\ncheck "<<i<<endl;
        vector<tXoodooState> res_xoo;
        for(int m=0; m<res.size(); m++) {
            tXoodooState tmp;
            State2XooState(res[m],tmp);
            res_xoo.push_back(tmp);
        }
        for(int m=0; m<print_size-1; m++) {
            tXoodooState tmp = lambda(res_xoo[m]);
            cout<<"a"<<m+1<<":"<<endl;
            displayXooState(cout,res_xoo[m]);
            cout<<endl;
            cout<<"b"<<m+1<<":"<<endl;
            displayXooState(cout,tmp);
            cout<<endl;
        }
        if(mode == 0) {
            cout<<"check b"<<print_size-1<<":"<<endl;
            displayXooState(cout,res_xoo[res_xoo.size()-1]);
            cout<<endl;
        }
        res.clear();
    }
    infile.close();
}

vector<vector<int>> XoodooRound::compare_trails(const string path_res1, int mode1, const string path_res2, int mode2) {
    ifstream in_res1(path_res1.data(), ios::in), in_res2(path_res2.data(), ios::in);
    States res1, res2;
    int res2_solution_num = 0;
    while(read2States(res2,in_res2,mode2)) res2_solution_num++;
    in_res2.close();

    map<int,int> same_solution, inverse_same_solution;
    for(int k=0; k<res2_solution_num; k++) inverse_same_solution[k] = -1;
    vector<vector<int>> output(2);

    for(int i=0; read2States(res1,in_res1,mode1); i++) {
        in_res2.open(path_res2.data(), ios::in);
        for(int j=0; read2States(res2,in_res2,mode2); j++) {
            bool equal = true;
            auto dx_dz1 = genSmallestState(res1[0]);
            auto dx_dz2 = genSmallestState(res2[0]);
            for(int k=0; k < min(res1.size(),res2.size()); k++) equal &= (ShiftXZ(dx_dz1,res1[k]) == ShiftXZ(dx_dz2,res2[k]));
            same_solution[i] = -1;
            if(equal) {
                same_solution[i] = j;
                inverse_same_solution[j] = i;
                cout<<i<<" in "+path_res1<<" trail equal to the "<<j<<" trail in "+path_res2<<endl;
                break;
            }
        }
        in_res2.close();
    }
    in_res1.close();

    for(auto it=same_solution.begin(); it!=same_solution.end(); it++) {
        if(it->second == -1) {
            cout<<it->first<<" trail in "+path_res1+" not found equal"<<endl;
            output[0].push_back(it->first);
        }
    }
    for(auto it=inverse_same_solution.begin(); it!=inverse_same_solution.end(); it++) {
        if(it->second == -1) {
            cout<<it->first<<" trail in "+path_res2+" not found equal"<<endl;
            output[1].push_back(it->first);
        }
    }

    return output;
}

void XoodooRound::gen_obj_T(vector<vector<int>> &obj) {
    string AS_cnf_num;//var list for pysat.card function, see pysat_card_AS.py for more information
    gen_extend_AS_cnf_num(AS_cnf_num);
    string py_cmd = "python3 pysat_card_AS.py --bound_num "+to_string(AS_weight_num)+" --var_num "+AS_cnf_num+" --mode "+to_string(AS_mode);

    FILE *fp = NULL;
    char buf[128] = {0}, result[256] = {0};
    if( (fp = popen(py_cmd.data(), "r")) != NULL ) {//run the py command and get print() output of the command
        while(fgets(buf, sizeof(buf), fp) != NULL) {
            strcat(result, buf);
        }
        pclose(fp);
        fp = NULL;
    }
    replace(AS_cnf_num.begin(),AS_cnf_num.end(),' ','_');

    int cnf_nvars = atoi(result);//the output is the number of variables of cnf
    solver.new_vars(cnf_nvars - AS_var_num);//add vars
    vector<string> mode = {"<=", ">=", "="};//atmost, atleast, euqals
    read2Vector(obj, objFilePath + "CNF_"+AS_cnf_num+"AS"+mode[AS_mode]+to_string(AS_weight_num)+".txt");
}


/*
 a0 -> b0 -> a1 -> b1 -> a2 -> b2 -> a3
 */
void XoodooRound::XoodooRound_AS() {
    map<unsigned int, unsigned int> RhoE_Relation, inverse_RhoE_Relation;
    gen_RhoE_T(RhoE_Relation,inverse_RhoE_Relation); // {i,j} {false}
    map<unsigned int, unsigned int> RhoW_Relation;
    gen_RhoW_T(RhoW_Relation); // {i,j} {false}
    vector<vector<unsigned int>> Theta_Relation = {};
    gen_Theta_T(Theta_Relation);
    cout<<"rhoE,inverse_rhoE,rhoW size: "<<RhoE_Relation.size()<<" "<<inverse_RhoE_Relation.size()<<" "<<RhoW_Relation.size()<<" "<<Theta_Relation.size()<<endl;

    //linear layer
    //a1->b1, a2->b2, ...
    cout<< "start add linear" <<endl;
    for (int round=0; round<round_num-1; round++) {
        for (int i=0; i<var_num; i++) {
            vector<unsigned int> temp={};
            temp.push_back(RhoW_Relation[Theta_Relation[RhoE_Relation[i]][0]] + var_num + round_var_num * round);//output bit position

            for (int j=1; j<Y*2+2; j++) {
                temp.push_back(inverse_RhoE_Relation[Theta_Relation[RhoE_Relation[i]][j]] + round_var_num * round);//input bit position
            }
            solver.add_xor_clause(temp, false);
        }
    }
    Theta_Relation.clear();
    RhoE_Relation.clear();
    inverse_RhoE_Relation.clear();
    RhoW_Relation.clear();
    
    /*
    chi
    generate bi --> a_i+1
    */
    vector<vector<string>> Chi_Relation;
    gen_xoodoo_Chi_cnf(Chi_Relation);
    //b1->a2
    cout<< "start add chi" <<endl;
    for (int round=0; round<round_num-2; round++) {
        for (int x=0; x<X; x++) {
           for (int z=0; z<Z; z++) {
               for (int ithClause=0; ithClause<Chi_Relation.size(); ithClause++) {
                   vector<Lit> clause;//one of the 6 clauses of a colliding Sbox
                   for (int ithVar=0; ithVar<Y; ithVar++) {
                       if (Chi_Relation[ithClause][ithVar] != "0") {//input
                           if (Chi_Relation[ithClause][ithVar] == "T") {
                               clause.push_back(Lit((ithVar*X*Z + z + Z*x + var_num + round_var_num * round), true));
                           }
                           else{
                               clause.push_back(Lit((ithVar*X*Z + z + Z*x + var_num + round_var_num * round), false)); 
                           }
                       } 
                   }
                   
                   for (int ithVar=Y; ithVar<Y*2; ithVar++) {
                       if (Chi_Relation[ithClause][ithVar] != "0") {//output
                           if (Chi_Relation[ithClause][ithVar] == "T") {
                               clause.push_back(Lit(((ithVar-Y)*X*Z + z + Z*x + var_num*2 + round_var_num * round), true));
                           }
                           else{
                               clause.push_back(Lit(((ithVar-Y)*X*Z + z + Z*x + var_num*2 + round_var_num * round), false));
                           }
                       }
                   }
                   solver.add_clause(clause);
                   clause.clear();
               }   
           }
        }
    }
    Chi_Relation.clear();
    
    vector<vector<string>> AS_Relation;
    gen_xoodoo_AS_cnf(AS_Relation);
    cout<< "start add AS" <<endl;
    //a1, a2, b2, ...
    for (int x=0; x<X; x++) {
        for (int z=0; z<Z; z++) {
            //对每一个S盒check
            for (int ithClause=0; ithClause<AS_Relation.size(); ithClause++) {
                vector<vector<Lit>> clause(AS_state_num);
                for (int ithVar=0; ithVar<Y; ithVar++) {
                    if (AS_Relation[ithClause][ithVar] != "0") {
                        for(int i=0; i<AS_state_num; i++) {
                            int offset = i*round_var_num;//a1, a2, b2, ...
                            if(i == AS_state_num-1) offset -= var_num;
                            if (AS_Relation[ithClause][ithVar] == "T") {
                                clause[i].push_back(Lit((ithVar*X*Z + z + Z*x + offset), true));
                            }
                            else {
                                clause[i].push_back(Lit((ithVar*X*Z + z + Z*x + offset), false));
                            }
                        }
                    }
                }
                
                if (AS_Relation[ithClause][Y] != "0") {
                    for(int i=0; i<AS_state_num; i++) {
                        int offset = i*AS_node_var_num;//AS_a1, AS_a2, AS_b2, ...
                        if (AS_Relation[ithClause][Y] == "T") {
                            clause[i].push_back(Lit((var_num*core_state_num + z + Z*x + offset), true));
                        }
                        else {
                            clause[i].push_back(Lit((var_num*core_state_num + z + Z*x + offset), false));
                        }
                    }
                }
                for(int i=0; i<AS_state_num; i++) {
                    solver.add_clause(clause[i]);
                    clause[i].clear();
                }
            }
        }
    }
    AS_Relation.clear();

    vector<vector<int>> obj;
    gen_obj_T(obj);
    //AS_a1a2b2
    cout<< "start add AS obj" <<endl;
    for (int ithClause=0; ithClause< obj.size(); ithClause++) {
        vector<Lit> clause; //one of the 10 clauses of a colliding Sbox

        for (int ithVar=0; ithVar<obj[ithClause].size(); ithVar++) {
            int temp = obj[ithClause][ithVar];
            if (temp > 0) {
                    clause.push_back(Lit(abs(temp) - 1 + var_num*core_state_num, false));
                }
            else{
                clause.push_back(Lit(abs(temp) - 1 + var_num*core_state_num, true));
            }
        }
        solver.add_clause(clause);
        clause.clear();
    }
    obj.clear();
}

void XoodooRound::main() {
    //print date and time
    time_point<system_clock> start = system_clock::now();
    auto st = system_clock::to_time_t(start);
    struct tm* stm = localtime(&st);
    cout<<stm->tm_mon + 1<<" "<<stm->tm_mday<<" "<<stm->tm_hour<<":"<<stm->tm_min<<":"<<stm->tm_sec<<endl;
    
    //round function, generate clauses for solver
    XoodooRound_AS();

    //mkdir result and the empty result file
    ofstream out;
    string res_prefix = objFilePath+"result/";
    if(access(res_prefix.data(), F_OK) == -1) {//check if result dir is existed
        mkdir(res_prefix.data(), S_IRWXO|S_IRWXG|S_IRWXU);
    }

    vector<string> mode = {"<=", ">=", "="};//atmost, atleast, euqals
    string res_filepath = res_prefix+"xoodoo_result_"+to_string(round_num)+"R_#AS"+mode[AS_mode]+to_string(AS_weight_num)+".txt";
    out.open(res_filepath, ios::ate);//clear result txt first
    out.close();

    //ban previous solution and solve phase
    cout<<"start banning found trails with weight"+mode[AS_mode]<<AS_weight_num<<endl;
    string ban_pre_path = objFilePath+"found"+mode[AS_mode]+to_string(AS_weight_num)+".txt";
    ifstream in_ban_pre(ban_pre_path.data(), ios::in);
    States ban_pre;
    while(read2States(ban_pre,in_ban_pre,0)) {
        map<State,int> banned;
        for(int i=0; i<ban_pre.size()-1; i++) {
           banned[ban_pre[i]] = i*round_var_num;
        }
        ban_solution(banned);
    }
    in_ban_pre.close();
    
    cout<<"start solving"<<endl;
    solve_and_output(res_filepath);
}