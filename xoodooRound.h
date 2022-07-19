#ifndef XOODOOROUND_H
#define XOODOOROUND_H

#include <stdio.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <algorithm>
#include <unordered_map>
#include <set>
#include "gen_xoodoo_Chi_cnf.h"
#include "gen_xoodoo_AS_cnf.h"

namespace XOODOOSAT {

typedef unsigned int tXoodooLane;         // a lane is 32 bits
typedef vector<tXoodooLane> tXoodooState; // a state = 3*4 lanes
typedef vector<unsigned int> Bitmap;      // size: X*Y*Z, Bitmap A represents a state in binary array, A[i] = 0 or A[i] = 1
typedef vector<unsigned int> StateColumn; // size: X*Z, StateColumn[i] is a Y-bit column
typedef vector<unsigned int> State;       // represents a state; only contains the bits that are active in a state, e.g. {122,306}
typedef vector<State> States;             // multi states
typedef vector<StateColumn> StateColumns; // multi states stored in columns

#define X 4
#define Y 3
#define Z 32          // sizeof(tXoodooLane)*8
#define var_num X*Y*Z // the number of variables(bits) of a state

#define default_obj_file_path "./" // path to AS cnf

#define indexXZ(x, z) ( Z*((x)%X) + ((z)%Z) ) /* 32*x + z */
#define indexXY(x, y) ( ((x)%X) + X*((y)%Y) ) /* x + 4*y */
#define ROLxoo(a, offset) ((offset != 0) ? ((((tXoodooLane)a) << offset) ^ (((tXoodooLane)a) >> (sizeof(tXoodooLane)*8-offset))) : a)
#define RORxoo(a, offset) ((offset != 0) ? ((((tXoodooLane)a) >> offset) ^ (((tXoodooLane)a) << (sizeof(tXoodooLane)*8-offset))) : a)

/*
 3 round: a0 -> b0 -> a1 -> b1 -> a2 -> b2 -> a3, bi=lambda(ai), ai+1=chi(bi)
 3 round trail core: a1 -> b1 -> a2 -> b2
 state to sum weight: a1, a2, b2 (a1, a2,..., an-1, bn-1)
*/

class XoodooRound
{
/*
core_state_num: how many states in the trail core
AS_weight_num: trail weight to bound
AS_state_num: how many states in the trail core to sum weight; in the example above, w(a1)+w(a2)+w(b2), 3 states to sum weight
round_num: how many rounds to analysis
AS_var_num: the number of variables(bits) of all AS nodes
AS_node_var_num: the number of variables(bits) of a AS node
AS_mode: weight sum mode, 0 for atmost, 1 for atleast, 2 for equals.
round_var_num: the number of variables(bits) of a round
*/
public:
    // 0 for differential, 1 for linear
    int analysis_mode;
    int core_state_num, AS_weight_num, AS_state_num, round_num;
    int AS_var_num, AS_node_var_num, AS_mode, round_var_num;
    int rho_plane[2], theta_L1[2], theta_L2[2], rhoW_L1[2], rhoW_L2[2], rhoE_L1[2], rhoE_L2[2];//parameters in the theta and rho process
    unsigned int theta_order;
	// path to save the AS cnf file
    string objFilePath;
    int thread_num;
    // the solver
    SATSolver solver;

    // a RhoE/RhoW relation is the input-output bit map, since RhoE/RhoW is a bijection
    // for example, RhoE_Relation[10] = 100 means: after RhoE, the 10th bit input maps to the 100th bit output
    map<unsigned int, unsigned int> RhoE_Relation, inverse_RhoE_Relation;
    map<unsigned int, unsigned int> RhoW_Relation, inverse_RhoW_Relation;
    // theta: a[x,y,z] = a[x,y,z] ^ \sum_{y=0}^{2} a[x-1,y,z-5] ^ \sum_{y=0}^{2} a[x-1,y,z-14]
    // 1 output bit is related to 7 input bits
    // Theta_Relation[i] = {i1,...,i7}
    vector<vector<unsigned int>> Theta_Relation, transpose_Theta_Relation;
    // the cnf string of chi
    vector<vector<string>> Chi_Relation;
    // the cnf string of AS
    vector<vector<string>> AS_Relation;
    // the cnf of weight
    vector<vector<int>> obj;

    /**
     * @brief Construct a new Xoodoo Round object
     * 
     * @param analy_mode int, 0 for differential, 1 for linear
     * @param rounds int, number of rounds
     * @param weight int, the weight bound
     * @param thread int, number of threads
     * @param mode int, 0 for <=weight bound, 1 for =weight bound, 2 for >=weight bound
     */
    XoodooRound(int analysis_mode, int rounds, int weight, int thread, int mode);

    // compute the order of theta transformation
    unsigned int compute_Theta_Order();

    // functions in a Xoodoo round
    void theta(tXoodooState &A);
    void transposetheta(tXoodooState &A);
    void rhoW(tXoodooState &A);
    void rhoE(tXoodooState &A);
    void inversetheta(tXoodooState &A);
    void inverserhoW(tXoodooState &A);
    void inverserhoE(tXoodooState &A);
    void chi(tXoodooState &A);
    tXoodooState lambda(tXoodooState A);
    tXoodooState inverselambda(tXoodooState A);
    tXoodooState transposelambda(tXoodooState A);
	
    // data struct converters
    void Bit2XooState(const Bitmap &A, tXoodooState &B);
    void XooState2Bit(const tXoodooState &A, Bitmap &B);
    void Bit2State(const Bitmap &A, State &B);
    void State2Bit(const State &A, Bitmap &B);
    void Bit2StateColumn(const Bitmap &A, StateColumn &B);
    void StateColumn2Bit(const StateColumn &A, Bitmap &B);
    void Bit2Plane(const Bitmap &A, vector<tXoodooLane> &B);//A,B are planes in different format
    void Plane2Bit(const vector<tXoodooLane> &A, Bitmap &B);
    void State2XooState(const State &A, tXoodooState &B);
    void XooState2State(const tXoodooState &A, State &B);
    void State2StateColumn(const State &A, StateColumn &B);
    void StateColumn2State(const StateColumn &A, State &B);
    void XooState2StateColumn(const tXoodooState &A, StateColumn &B);
    void StateColumn2XooState(const StateColumn &A, tXoodooState &B);


    // calculate the weight of a tXoodooState
    int caculateXooStateWeight(const tXoodooState &A);
    // calculate the weight of a StateColumn
    int caculateStateColumnWeight(const StateColumn &A);


    // State << (dx,dz)
    State ShiftXZ(const vector<int> &dx_dz, const State &A);
    // whether A could shiftXZ to B
    bool StateEqualAfterShift(const State &A, const State &B);
    // compare 2 states
    bool isSmaller(const State &A, const State &B);
    // get the smallest state of the ShiftXZ equivalent class of State A
    vector<int> genSmallestState(const State &A);


    // display State in a readable form, see xoodoo_result_3R_#AS<=25.txt
    void display(ostream& fout, const State &A);
    // display tXoodooState in a readable form, see xoodoo_result_3R_#AS<=25.txt
    void displayXooState(ostream& fout, const tXoodooState &A);


    // generate maps for rhow(state) and rhow^-1(state)
    void gen_RhoW_T(map<unsigned int, unsigned int>& RhoW_index, map<unsigned int, unsigned int>& inverse_RhoW_index);
    // generate maps for rhoe(state) and rhoe^-1(state)
    void gen_RhoE_T(map<unsigned int, unsigned int>& RhoE_index, map<unsigned int, unsigned int>& inverse_RhoE_index);
    // generate theta relation
    void gen_Theta_T(vector<vector<unsigned int>>& relation, vector<vector<unsigned int>>& transpose_relation);
    // generate cnf for AS weight bound
    void gen_obj_T(vector<vector<int>> &Obj, SATSolver &Solver, int weight, int as_var_num, int as_mode);
    // get var list for pysat.card function
    void gen_extend_AS_cnf_num(string &cnf_num);
    

    // add clauses to solver
    void add_lambda2solver(SATSolver &Solver, int rounds, int base_offset = 0);
    void add_chi2solver(SATSolver &Solver, int rounds, int base_offset = 0);
    void add_AS2solver(SATSolver &Solver, int rounds, int base_offset = 0, int core_var_num = 0);
    void add_ASobj2solver(SATSolver &Solver, vector<vector<int>> &Obj, int core_var_num = 0);


    // ban found solutions and shifted solutions
    void ban_solution(SATSolver &Solver, const map<State,int> &A);
    // get the trail weight by caculating AS
    int get_weight(SATSolver &Solver, int rounds, int core_var_num = 0);
    

    // writing solution to a file in a readable format
    void write_result(const string pathname, const int weight, const int solution_count, const States &solution);
    
    // solve and output the result
    void solve_and_output(SATSolver &Solver, const string pathname, int rounds, int base_offset = 0, int core_var_num = 0, const vector<Lit> &assumption = {});
	

    // read cnf generated by python-sat
    void read2Vector(vector<vector<int>>& res, const string pathname);
    // read result into States in a readable file
    // mode 0 for file generated by this programme, 1 for DC file generated by previous work(https://github.com/KeccakTeam/Xoodoo/tree/master/XooTools/Trails)
    bool read2States(States& res, ifstream &infile, int mode, const vector<string> &st = {"ab","NE","bb"});
    // read result into StateColumns in a readable file
    int read2StateColumn(StateColumns& res, ifstream &infile, int mode, const vector<string> &st = {"ab","NE","bb"});


    // check and print the i-th trail results in a readable file
    void check_trails(const string pathname, int i, int mode);
    // compare 2 trail results file, see if there's any difference
    // mode 0 for the file generated by this programme, 1 for file generated by previous work
    vector<vector<int>> compare_trails(const string path_res1, int mode1, const string path_res2, int mode2);
    

    // read assumption according to infile
    bool read2assumptions(vector<vector<Lit>> &assumptions, ifstream &infile, int mode = 2, int offset = 0, const vector<string> &st = {"ab","NE","bb"});


    // differential/linear analysis round function (TODO: trail extension), generate clauses for solver
    void XoodooRound_AS();
    // trail extension
    void extendRound_AS();
    // run the whole Xoodoo analysis process
    void main();
    // generate DDT table for chi
    void gen_chi_DDT(unordered_map<int,set<int>> &ddt);
    // extend a round by dfs
    void dfs(StateColumn &assumption, vector<int> &nonzero_col_index, unordered_map<int, set<int>> &ddt, StateColumn temp, int i, int n, int &total_weight, int &pre_weight, bool forward = true);
    // run the whole Xoodoo analysis process with trail extension
    void extend_main(bool forward = true);
};


}// end namespace

#endif
