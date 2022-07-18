//
//  xoodooRound.cpp
//
//  Created by Haochen on 05/15/21.
//
#include <cstdlib>
#include <chrono>
#include <ctime>
#include "xoodooRound.h"
using namespace std::chrono;
using namespace XOODOOSAT;

/**
 * @brief Construct a new XoodooRound object
 * 
 * @param analy_mode int, 0 for differential, 1 for linear
 * @param rounds int, number of rounds
 * @param weight int, the weight bound
 * @param thread int, number of threads
 * @param mode int, 0 for <=weight bound, 1 for =weight bound, 2 for >=weight bound
 */
XoodooRound::XoodooRound(int analy_mode, int rounds, int weight, int thread, int mode) {
    // set parameter
    analysis_mode = analy_mode;
    AS_weight_num = weight;
    AS_state_num = rounds;
    AS_mode = mode;
    round_num = rounds;

    // trail core a1,b1,a2,b2,...
    // number of states in the trail core
    core_state_num = 2 * (round_num - 1);
    // total number of AS in all states
    AS_var_num = X * Z * AS_state_num;
    // number of AS in one state(node)
    AS_node_var_num = X * Z;
    // total number of bits in all states
    round_var_num = var_num * 2;
    // 
    objFilePath = default_obj_file_path;
    thread_num = thread;

    theta_order = compute_Theta_Order();

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
    solver.new_vars(var_num*core_state_num);//add vars

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

    gen_RhoE_T(RhoE_Relation, inverse_RhoE_Relation);
    gen_RhoW_T(RhoW_Relation, inverse_RhoW_Relation);
    gen_Theta_T(Theta_Relation, transpose_Theta_Relation);
    gen_xoodoo_Chi_cnf(Chi_Relation);
    gen_xoodoo_AS_cnf(AS_Relation);
    gen_obj_T(obj, solver, AS_weight_num, AS_var_num, AS_mode);
}

/**
 * @brief compute the order of theta transformation
 * 
 * @return unsigned int, the order
 */
unsigned int XoodooRound::compute_Theta_Order() {
    unsigned int oddPartSizeX = X;
    unsigned int powerTwoPartSizeX = 1;
    while ((oddPartSizeX & 1) == 0) {
        oddPartSizeX >>= 1;
        powerTwoPartSizeX <<= 1;
    }
    unsigned int order = Z;
    if (powerTwoPartSizeX > order) order = powerTwoPartSizeX;
    switch (oddPartSizeX) {
        case 1:
            break;
        case 3:
            order *= 3;
            break;
        case 5:
            order *= 15;
            break;
        case 7:
            order *= 7;
            break;
            throw ((string)"X is not a power of two times, 3, 5 or 7");
    }
    return order;
}

/**
 * @brief theta transformation
 * 
 * @param A tXoodooState&
 */
void XoodooRound::theta(tXoodooState &A) {
    unsigned int x, y;
    vector<tXoodooLane> P(X, 0), E(X, 0);

    for (x = 0; x < X; x++) {
        for (y = 0; y < Y; y++)
            P[x] ^= A[indexXY(x, y)];
    }
    for (x = 0; x < X; x++)
        E[x] = ROLxoo(P[(x + X - theta_L1[0]) % X], theta_L1[1]) ^ ROLxoo(P[(x + X - theta_L2[0]) % X], theta_L2[1]);
    for (x = 0; x < X; x++)
        for (y = 0; y < Y; y++)
            A[indexXY(x, y)] ^= E[x];
}

/**
 * @brief transpose theta transformation
 * 
 * @param A tXoodooState&
 */
void XoodooRound::transposetheta(tXoodooState &A) {
    unsigned int x, y;
    vector<tXoodooLane> P(X, 0), E(X, 0);

    for (x = 0; x < X; x++) {
        for (y = 0; y < Y; y++)
            P[x] ^= A[indexXY(x, y)];
    }
    for (x = 0; x < X; x++)
        E[x] = RORxoo(P[(x + X + theta_L1[0]) % X], theta_L1[1]) ^ RORxoo(P[(x + X + theta_L2[0]) % X], theta_L2[1]);
    for (x = 0; x < X; x++)
        for (y = 0; y < Y; y++)
            A[indexXY(x, y)] ^= E[x];
}

/**
 * @brief the inverse of theta transformation
 * 
 * @param A tXoodooState&
 */
void XoodooRound::inversetheta(tXoodooState &A) {
    unsigned int x, y;
    vector<tXoodooLane> P(X, 0);
    unsigned int exponent = theta_order - 1;
    unsigned int powerTwo = 1;

    for (x = 0; x < X; x++) {
        for (y = 0; y < Y; y++)
            P[x] ^= A[indexXY(x, y)];
    }
    vector<tXoodooLane> E(P);
    do {
        if ((exponent & powerTwo) != 0) {
            vector<tXoodooLane> tmp(E);
            for (x = 0; x < X; x++)
                E[x] = tmp[x] ^ ROLxoo(tmp[(x + X - theta_L1[0] * powerTwo) % X], (theta_L1[1] * powerTwo) % Z) ^ ROLxoo(tmp[(x + X - theta_L2[0] * powerTwo) % X], (theta_L2[1] * powerTwo) % Z);
        }
        powerTwo <<= 1;
    } while (powerTwo <= exponent);
    for (x = 0; x < X; x++)
        E[x] ^= P[x];
    for (x = 0; x < X; x++)
        for (y = 0; y < Y; y++)
            A[indexXY(x, y)] ^= E[x];
}

/**
 * @brief rho west transformation
 * 
 * @param A tXoodooState&
 */
void XoodooRound::rhoW(tXoodooState &A) {
    unsigned int x, y;
    tXoodooState tempA(A);

    y = rho_plane[0];//A1
    for (x = 0; x < X; x++)
        A[indexXY(x, y)] = ROLxoo(tempA[indexXY(x + X - rhoW_L1[0], y)], rhoW_L1[1]);//(1,0)
    y = rho_plane[1];//A2
    for (x = 0; x < X; x++)
        A[indexXY(x, y)] = ROLxoo(tempA[indexXY(x + X - rhoW_L2[0], y)], rhoW_L2[1]);//(0,11)
}

/**
 * @brief the inverse of rho west transformation
 * 
 * @param A tXoodooState&
 */
void XoodooRound::inverserhoW(tXoodooState &A) {
    unsigned int x, y;
    tXoodooState tempA(A);

    y = rho_plane[0];//A1
    for (x = 0; x < X; x++)
        A[indexXY(x, y)] = RORxoo(tempA[indexXY(x + X + rhoW_L1[0], y)], rhoW_L1[1]);//(1,0)
    y = rho_plane[1];//A2
    for (x = 0; x < X; x++)
        A[indexXY(x, y)] = RORxoo(tempA[indexXY(x + X + rhoW_L2[0], y)], rhoW_L2[1]);//(0,11)
}

/**
 * @brief rho east transformation
 * 
 * @param A tXoodooState&
 */
void XoodooRound::rhoE(tXoodooState &A) {
    unsigned int x, y;
    tXoodooState tempA(A);

    y = rho_plane[0];//A1
    for (x = 0; x < X; x++)
        A[indexXY(x, y)] = ROLxoo(tempA[indexXY(x + X - rhoE_L1[0], y)], rhoE_L1[1]);//(0,1)
    y = rho_plane[1];//A2
    for (x = 0; x < X; x++)
        A[indexXY(x, y)] = ROLxoo(tempA[indexXY(x + X - rhoE_L2[0], y)], rhoE_L2[1]);//(2,8)
}

/**
 * @brief the inverse of rho east transformation
 * 
 * @param A tXoodooState&
 */
void XoodooRound::inverserhoE(tXoodooState &A) {
    unsigned int x, y;
    tXoodooState tempA(A);

    y = rho_plane[0];//A1
    for (x = 0; x < X; x++)
        A[indexXY(x, y)] = RORxoo(tempA[indexXY(x + X + rhoE_L1[0], y)], rhoE_L1[1]);//(0,1)
    y = rho_plane[1];//A2
    for (x = 0; x < X; x++)
        A[indexXY(x, y)] = RORxoo(tempA[indexXY(x + X + rhoE_L2[0], y)], rhoE_L2[1]);//(2,8)
}

/**
 * @brief chi transformation
 * 
 * @param A tXoodooState&
 */
void XoodooRound::chi(tXoodooState &A) {
    unsigned int x, y;
    vector<tXoodooLane> B(Y, 0);

    for (x = 0; x < X; x++) {
        for (y = 0; y < Y; y++)
            B[y] = A[indexXY(x, y)] ^ ((~A[indexXY(x, y + 1)]) & A[indexXY(x, y + 2)]);
        for (y = 0; y < Y; y++)
            A[indexXY(x, y)] = B[y];
    }
}

/**
 * @brief lambda transformation, including all linear transformations
 * 
 * @param A tXoodooState
 * @return tXoodooState 
 */
tXoodooState XoodooRound::lambda(tXoodooState A) {
    rhoE(A);
    theta(A);
    rhoW(A);
    return A;
}

/**
 * @brief transpose lambda transformation, including all linear transformations
 * 
 * @param A tXoodooState
 * @return tXoodooState 
 */
tXoodooState XoodooRound::transposelambda(tXoodooState A) {
    inverserhoW(A);
    transposetheta(A);
    inverserhoE(A);
    return A;
}

/**
 * @brief teh inverse of lambda transformation, including all linear transformations
 * 
 * @param A tXoodooState
 * @return tXoodooState 
 */
tXoodooState XoodooRound::inverselambda(tXoodooState A) {
    inverserhoW(A);
    inversetheta(A);
    inverserhoE(A);
    return A;
}

/**
 * @brief Bitmap to tXoodooState
 * 
 * @param A input const Bitmap&
 * @param B output tXoodooState&
 */
void XoodooRound::Bit2XooState(const Bitmap &A, tXoodooState &B) {
    B.assign(X*Y, 0);
    unsigned int x, y, z;
    for (y = 0; y < Y; y++) {
        for (x = 0; x < X; x++) {
            B[X*y + x] = 0;
        }
        for (z = 0; z < Z; z++) {
            for (x = 0; x < X; x++) {
                B[X*y + x] |= ((tXoodooLane)(A[Z*X*y + Z * x + z]) << (z)); // 32*(4y+x) + z
            }
        }
    }
}

/**
 * @brief tXoodooState to Bitmap
 * 
 * @param A input const tXoodooState&
 * @param B output Bitmap&
 */
void XoodooRound::XooState2Bit(const tXoodooState &A, Bitmap &B) {
    B.assign(var_num, 0);
    unsigned int x, y, z;
    for (y = 0; y < Y; y++) {
        for (z = 0; z < Z; z++) {
            for (x = 0; x < X; x++) {
                B[Z*X*y + Z * x + z] = ((A[X*y + x] >> (z)) & 0x1); // 32*(4y+x) + z
            }
        }
    }
}

/**
 * @brief Bitmap to State
 * 
 * @param A input const Bitmap&
 * @param B output State&
 */
void XoodooRound::Bit2State(const Bitmap &A, State &B) {
    B.clear();
    for (unsigned int i = 0; i < var_num; i++) {
        if (A[i] == 1) B.push_back(i);
    }
}

/**
 * @brief State to Bitmap
 * 
 * @param A input const State&
 * @param B output Bitmap&
 */
void XoodooRound::State2Bit(const State &A, Bitmap &B) {
    B.assign(var_num, 0);
    for (int i = 0; i < A.size(); i++) {
        B[A[i]] = 1;
    }
}

/**
 * @brief Bitmap to StateColumn
 * 
 * @param A input const Bitmap&
 * @param B output StateColumn&
 */
void XoodooRound::Bit2StateColumn(const Bitmap &A, StateColumn &B) {
    B.assign(X*Z, 0);
    for (int z = 0; z < Z; z++) {
        for (int x = 0; x < X; x++) {
            unsigned int col = 0;
            for (int y = 0; y < Y; y++) {
                col += (A[Z*X*y + Z * x + z] << y);
            }
            B[Z*x + z] = col;
        }
    }
}

/**
 * @brief StateColumn to Bitmap
 * 
 * @param A input const StateColumn&
 * @param B output Bitmap&
 */
void XoodooRound::StateColumn2Bit(const StateColumn &A, Bitmap &B) {
    B.assign(var_num, 0);
    for (int z = 0; z < Z; z++) {
        for (int x = 0; x < X; x++) {
            unsigned int col = A[Z*x + z];
            for (int y = 0; y < Y; y++) {
                B[Z*X*y + Z * x + z] = (col >> y) & 0x1;
            }
        }
    }
}

/**
 * @brief Bitmap to vector<tXoodooLane>
 * 
 * @param A input const Bitmap&
 * @param B output vector<tXoodooLane>&
 */
void XoodooRound::Bit2Plane(const Bitmap &A, vector<tXoodooLane> &B) {//32*x + z
    B.assign(X, 0);
    unsigned int x, z;

    for (x = 0; x < X; x++) {
        B[x] = 0;
    }
    for (z = 0; z < Z; z++) {
        for (x = 0; x < X; x++) {
            B[x] |= ((tXoodooLane)(A[Z*x + z]) << (z));
        }
    }
}

/**
 * @brief vector<tXoodooLane> to Bitmap
 * 
 * @param A input const vector<tXoodooLane>&
 * @param B output Bitmap&
 */
void XoodooRound::Plane2Bit(const vector<tXoodooLane> &A, Bitmap &B) {//32*x + z
    unsigned int x, z;

    for (z = 0; z < Z; z++) {
        for (x = 0; x < X; x++) {
            B[Z*x + z] = ((A[x] >> (z)) & 0x1);
        }
    }
}

/**
 * @brief State to tXoodooState
 * 
 * @param A input const State&
 * @param B output tXoodooState&
 */
void XoodooRound::State2XooState(const State &A, tXoodooState &B) {
    Bitmap BitA;//Bit
    State2Bit(A, BitA);
    Bit2XooState(BitA, B);
}

/**
 * @brief tXoodooState to State
 * 
 * @param A input const tXoodooState&
 * @param B output State&
 */
void XoodooRound::XooState2State(const tXoodooState &A, State &B) {
    Bitmap BitA(var_num, 0);
    XooState2Bit(A, BitA);
    Bit2State(BitA, B);
}

/**
 * @brief State to StateColumn
 * 
 * @param A input const State&
 * @param B output StateColumn&
 */
void XoodooRound::State2StateColumn(const State &A, StateColumn &B) {
    Bitmap BitA;
    State2Bit(A, BitA);
    Bit2StateColumn(BitA, B);
}

/**
 * @brief StateColumn to State
 * 
 * @param A input const StateColumn&
 * @param B output State&
 */
void XoodooRound::StateColumn2State(const StateColumn &A, State &B) {
    Bitmap BitA(var_num, 0);
    StateColumn2Bit(A, BitA);
    Bit2State(BitA, B);
}

/**
 * @brief tXoodooState to StateColumn
 * 
 * @param A input const tXoodooState&
 * @param B output StateColumn&
 */
void XoodooRound::XooState2StateColumn(const tXoodooState &A, StateColumn &B) {
    Bitmap BitA(var_num, 0);
    XooState2Bit(A, BitA);
    Bit2StateColumn(BitA, B);
}

/**
 * @brief StateColumn to tXoodooState
 * 
 * @param A input const StateColumn&
 * @param B output tXoodooState&
 */
void XoodooRound::StateColumn2XooState(const StateColumn &A, tXoodooState &B) {
    Bitmap BitA;
    StateColumn2Bit(A, BitA);
    Bit2XooState(BitA, B);
}

/**
 * @brief given a tXoodooState, calculate its weight
 * 
 * @param A(const tXoodooState&) input state
 * @return int 
 */
int XoodooRound::caculateXooStateWeight(const tXoodooState &A) {
    int weight = 0;
    for (int x = 0; x < X; x++) {
        for (int z = 0; z < Z; z++) {
            for (int y = 0; y < Y; y++) {
                // check each column for active bits
                if ((A[indexXY(x, y)] >> (z)) & 0x1) {
                    // found an active bit, weight++
                    ++weight;
                    // after weight++, exclude this column(using "break")
                    break;
                }
            }
        }
    }
    return weight;
}

/**
 * @brief given a StateColumn, calculate its weight
 * 
 * @param A (const StateColumn &) input state in column form, A[i] is a column
 * @return int 
 */
int XoodooRound::caculateStateColumnWeight(const StateColumn &A) {
    int weight = 0;
    for (int x = 0; x < X; x++) {
        for (int z = 0; z < Z; z++) {
            if (A[Z*x + z]) ++weight;
        }
    }
    return weight;
}

/**
 * @brief State << (dx, dz)
 * 
 * @param dx_dz const vector<int> &, (dx, dz), shift bits
 * @param A const State &, state to shift
 * @return State 
 */
State XoodooRound::ShiftXZ(const vector<int> &dx_dz, const State &A) {
    // first convert State to tXoodooState
    tXoodooState stateA(X*Y, 0);
    State2XooState(A, stateA);
    tXoodooState shiftA(stateA);

    for (int x = 0; x < X; x++) {
        for (int y = 0; y < Y; y++) {
            // shift (dx,dz)
            shiftA[indexXY(x, y)] = ROLxoo(stateA[indexXY(x + X - dx_dz[0], y)], dx_dz[1]);
        }
    }
    State res;
    // convert tXoodooState to State
    XooState2State(shiftA, res);
    return res;
}

/**
 * @brief check if State A and B are symmetry (shift (x,z))
 * 
 * @param A const State &
 * @param B const State &
 * @return bool
 */
bool XoodooRound::StateEqualAfterShift(const State &A, const State &B) {//A, B is sorted
    if (A.size() != B.size()) return false;
    if (A == B) return true;

    // shift (x,z)
    for (int dx = 0; dx < X; dx++) {
        for (int dz = 0; dz < Z; dz++) {
            State shiftA = ShiftXZ({ dx,dz }, A);
            if (shiftA == B) return true;
        }
    }
    return false;
}

/**
 * @brief check if State A < B
 *          value of A = \sum_{i=0}^{383} pow(2, i)*bits_of_A[i]
 * 
 * @param A const State &
 * @param B const State &
 * @return bool 
 */
bool XoodooRound::isSmaller(const State &A, const State &B) {//state A<B; A, B is sorted
    int i = A.size() - 1, j = B.size() - 1;
    while (i >= 0 && j >= 0) {
        if (A[i] < B[j]) return true;
        else if (A[i] > B[j]) return false;
        else {
            i--;
            j--;
        }
    }
    return i < j;
}

/**
 * @brief by shifting (x,z), shift State A to its smallest symmetry State
 * 
 * @param A const State &
 * @return {int, int} dx_dz, the shift position corresponding to the smallest symmetry state
 */
vector<int> XoodooRound::genSmallestState(const State &A) {//A is sorted  
    vector<int> dx_dz = { 0, 0 };
    State smallest(A);
    for (int dx = 0; dx < X; dx++) {
        for (int dz = 0; dz < Z; dz++) {
            State shiftA = ShiftXZ({ dx,dz }, A);
            //after shift (dx,dz), get the smaller state
            if (isSmaller(shiftA, smallest)) {
                dx_dz = { dx,dz };
                smallest = shiftA;
            }
        }
    }

    return dx_dz;
}

/**
 * @brief print State A
 * 
 * @param fout ostream&
 * @param A const State&
 */
void XoodooRound::display(ostream& fout, const State &A) {
    vector<unsigned int> tempA(X*Z, 0);//State column map, 0<=tempA[i]<=3
    for (int i = 0; i < A.size(); i++) {
        int y;
        y = A[i] / (X*Z);
        assert(y < Y);
        tempA[A[i] - y * X*Z] |= ((0x1) << y);
    }
    for (int x = 0; x < X; x++) {
        for (int z = 0; z < Z; z++) {
            if (tempA[Z*x + z] == 0) fout << ".";
            else fout << tempA[Z*x + z];
        }
        fout << endl;
    }
}

/**
 * @brief print tXoodooState A
 * 
 * @param fout ostream&
 * @param A const tXoodooState &
 */
void XoodooRound::displayXooState(ostream& fout, const tXoodooState &A) {
    State tmp;
    XooState2State(A, tmp);
    display(fout, tmp);
}

/**
 * @brief generate RhoW and inverse_RhoW's input-output bit relation
 * 
 * @param RhoW_index map<unsigned int, unsigned int>&, 384 size map, {input bit, output bit}
 * @param inverse_RhoW_index map<unsigned int, unsigned int>&, 384 size map, {input bit, output bit}
 */
void XoodooRound::gen_RhoW_T(map<unsigned int, unsigned int>& RhoW_index, map<unsigned int, unsigned int>& inverse_RhoW_index) {
    Bitmap BitA(var_num, 0);
    tXoodooState temp(X*Y, 0);
    unsigned int i, j;

    for (i = 0; i < BitA.size(); i++) {
        // only set the input bit to 1
        // the other input bits are 0
        // then only the corresponding output bit is 1
        // the other output bits are 0
        BitA[i] = 1;

        // convert to tXoodooState
        Bit2XooState(BitA, temp);
        // RhoW
        rhoW(temp);
        // convert back to Bitmap
        XooState2Bit(temp, BitA);
        // find the corresponding output bit, which is 1
        for (j = 0; j < BitA.size(); j++) {
            if (BitA[j] == 1) {
                // log the input-output bit relation
                // i is the input
                // j is the output
                RhoW_index.insert(pair<unsigned int, unsigned int>(i, j));
                // inverse is (j, i)
                inverse_RhoW_index.insert(pair<unsigned int, unsigned int>(j, i));
            }
            // set all bits back to 0
            BitA[j] = 0;
        }
    }
    return;
}

/**
 * @brief generate RhoE and inverse_RhoE's input-output bit relation
 * 
 * @param RhoE_index map<unsigned int, unsigned int>&
 * @param inverse_RhoE_index map<unsigned int, unsigned int>&
 */
void XoodooRound::gen_RhoE_T(map<unsigned int, unsigned int>& RhoE_index, map<unsigned int, unsigned int>& inverse_RhoE_index) {
    Bitmap BitA(var_num, 0);
    tXoodooState temp(X*Y, 0);
    unsigned int i, j;

    for (i = 0; i < BitA.size(); i++) {
        // only set the input bit to 1
        // the other input bits are 0
        // then only the corresponding output bit is 1
        // the other output bits are 0
        BitA[i] = 1;

        // convert to tXoodooState
        Bit2XooState(BitA, temp);
        // RhoE
        rhoE(temp);
        // convert back to Bitmap
        XooState2Bit(temp, BitA);
        // find the corresponding output bit, which is 1
        for (j = 0; j < BitA.size(); j++) {
            if (BitA[j] == 1) {
                // log the input-output bit relation
                // i is the input
                // j is the output
                RhoE_index.insert(pair<unsigned int, unsigned int>(i, j));
                // inverse is (j, i)
                inverse_RhoE_index.insert(pair<unsigned int, unsigned int>(j, i));
            }
            // set all bits back to 0
            BitA[j] = 0;
        }
    }
    return;
}

/**
 * @brief generate Theta and transpose Thetas input-output bit relation
 *          an output bit of Theta is related to 7 input bits of Theta
 *          a[x,y,z] = a[x,y,z] ^ \sum_{y=0}^{2} a[x-1,y,z-5] ^ \sum_{y=0}^{2} a[x-1,y,z-14]
 * 
 * @param relation vector<vector<unsigned int>>&
 * @param transpose_relation vector<vector<unsigned int>>&
 */
void XoodooRound::gen_Theta_T(vector<vector<unsigned int>>& relation, vector<vector<unsigned int>>& transpose_relation) {
    // 4*32=128 columns
    vector<vector<unsigned int>> column;
    // each column[(x,z)] contains 3 bit indexes of y, i.e. (x + y * 4)*32 + z
    for (unsigned int x = 0; x < X; x++) {
        for (unsigned int z = 0; z < Z; z++) {
            vector<unsigned int> aColumn = {};
            for (int y = 0; y < Y; y++) {
                aColumn.push_back((x + y * X)*Z + z);
            }
            column.push_back(aColumn);
        }
    }

    for (int x = 0; x < X; x++) {
        for (int y = 0; y < Y; y++) {
            for (int z = 0; z < Z; z++) {
                vector<unsigned int> temp, temp_transpose;
                // output bit a[x,y,z]
                temp.push_back(Z*(x + X * y) + z);
                temp_transpose.push_back(Z*(x + X * y) + z);
                // input bit a[x,y,z]
                temp.push_back(Z*(x + X * y) + z);
                temp_transpose.push_back(Z*(x + X * y) + z);
                // input bit \sum_{y=0}^{2} a[x-1,y,z-5] ^ \sum_{y=0}^{2} a[x-1,y,z-14]
                for (int i = 0; i < Y; i++) {
                    temp.push_back(column[Z*((x + X - theta_L1[0]) % X) + ((z + Z - theta_L1[1]) % Z)][i]);
                    temp.push_back(column[Z*((x + X - theta_L2[0]) % X) + ((z + Z - theta_L2[1]) % Z)][i]);
                    temp_transpose.push_back(column[Z*((x + X + theta_L1[0]) % X) + ((z + Z + theta_L1[1]) % Z)][i]);
                    temp_transpose.push_back(column[Z*((x + X + theta_L2[0]) % X) + ((z + Z + theta_L2[1]) % Z)][i]);
                }
                relation.push_back(temp);
                transpose_relation.push_back(temp_transpose);
            }
        }
    }
}

/**
 * @brief combine str cnf with to_string(AS_var_num)
 * 
 * @param cnf_num string &
 */
void XoodooRound::gen_extend_AS_cnf_num(string &cnf_num) {
    cnf_num += to_string(AS_var_num) + " ";
    /*if(round_num == 2) {
        cnf_num += to_string(AS_node_var_num) + " ";
    }*/
}

/**
 * @brief ban solution States and its symmetry States
 * 
 * @param Solver SATSolver &
 * @param A const map<State, int>&, A.first is the State to ban, A.second is the var offset
 */
void XoodooRound::ban_solution(SATSolver &Solver, const map<State, int> &A) {
    // check if the State is 0
    bool zero = true;
    for (auto iter = A.begin(); iter != A.end();iter++) {
        zero = zero && (iter->first.size() == 0);
    }

    for (int dx = 0; dx < X; dx++) {
        for (int dz = 0; dz < Z; dz++) {
            vector<Lit> ban_solutions;
            for (auto iter = A.begin(); iter != A.end();iter++) {
                // ban all shifted symmetric States
                State shiftA = ShiftXZ({ dx,dz }, iter->first);
                Bitmap BitA(var_num, 0);
                State2Bit(shiftA, BitA);
                for (uint32_t var = 0; var < BitA.size(); var++) {
                    ban_solutions.push_back(Lit(var + iter->second, (BitA[var] == 1) ? true : false));
                }
            }
            Solver.add_clause(ban_solutions);
            ban_solutions.clear();
            // if the State is 0, no need to shift since 0 << (x,z) = 0 for all (x,z)
            if (zero) return;
        }
    }
}

/**
 * @brief write a solution to a file in a readable form
 * 
 * @param pathname const string
 * @param weight const int
 * @param solution_count const int
 * @param solution const States&, a solution, containing r States, r is the round number
 */
void XoodooRound::write_result(const string pathname, const int weight, const int solution_count, const States &solution) {
    ofstream out(pathname, ios::out | ios::app);
    out << "weight " << weight << " solution " << solution_count << endl;

    for (int i = 0; i < solution.size() - 1; i++) {
        out << "a" << i + 1 << ":" << endl;
        display(out, solution[i]);
        out << endl;
    }
    out << "b" << solution.size() - 1 << ":" << endl;
    display(out, solution[solution.size() - 1]);
    out << "\n" << endl;
    out.close();
}

/**
 * @brief calculate the weight of the Solver's solution
 * 
 * @param Solver SATSolver &
 * @param rounds int
 * @param core_var_num int
 * @return int 
 */
int XoodooRound::get_weight(SATSolver &Solver, int rounds, int core_var_num) {
    // get each State weight(a1,a2,...)
    vector<int> as_weight(rounds, 0);
    for (int i = 0; i < AS_node_var_num; i++) {
        for (int j = 0; j < rounds; j++) {
            if (Solver.get_model()[i + core_var_num + AS_node_var_num * j] == l_True) {
                as_weight[j]++;
            }
        }
    }
    // get total weight
    int weight = 0;
    /*if(round_num == 2) weight = as_weight[0]*2 + as_weight[1];
    else {
        for(int j=0; j<rounds; j++) weight += as_weight[j];
    }*/
    for (int j = 0; j < rounds; j++) weight += as_weight[j];
    return weight;
}

/**
 * @brief solve all solutions by iteratively banning solved solutions and write solutions to a file
 * 
 * @param Solver SATSolver &
 * @param pathname const string
 * @param rounds int
 * @param base_offset int
 * @param core_var_num int
 * @param assumption const vector<Lit>&, you can solve with an assumption
 */
void XoodooRound::solve_and_output(SATSolver &Solver, const string pathname, int rounds, int base_offset, int core_var_num, const vector<Lit> &assumption) {
    int counting = 0;
    map<int, int> solution_counts;// map<weight, solution_counts>
    
    while (true) {
        lbool ret;
        if (assumption.size() == 0) {
            ret = Solver.solve();
        }
        else {
            cout << "got a assumption" << endl;
            ret = Solver.solve(&assumption);
        }
        // print date and time
        time_point<system_clock> start = system_clock::now();
        auto st = system_clock::to_time_t(start);
        struct tm* stm = localtime(&st);
        cout << '\n' << stm->tm_mon + 1 << " " << stm->tm_mday << " " << stm->tm_hour << ":" << stm->tm_min << ":" << stm->tm_sec << endl;

        if (ret != l_True) {
            assert(ret == l_False);
            cout << "reach end" << endl;
            cout << "counting: " << counting << endl;
            // All solutions found.
            exit(0);
        }

        // start processing solution
        // get weight
        int weight = get_weight(Solver, rounds, core_var_num);

        // get result in States format
        States solution(rounds);//states to show(a1,a2,...,bn-1)
        for (unsigned int i = 0; i < var_num; i++) {
            for (int j = 0; j < rounds; j++) {
                int offset = j * round_var_num;
                if (j == rounds - 1 && offset) offset -= var_num;
                if (Solver.get_model()[i + offset + base_offset] == l_True) {
                    solution[j].push_back(i);//State
                }
            }
        }
        // solution + 1
        auto iter = solution_counts.find(weight);
        if (iter != solution_counts.end()) {
            iter->second += 1;
        }
        else solution_counts[weight] = 1;

        // a new state, shift to smallest
        vector<int> dx_dz = genSmallestState(solution[0]);
        for (int i = 0; i < rounds; i++) {
            solution[i] = ShiftXZ(dx_dz, solution[i]);
        }

        // now extend result
        // TODO

        // writing final result
        if (rounds != 1) write_result(pathname, weight, solution_counts[weight], solution);
        if (rounds != 1 || counting % 1000 == 0) { cout << counting << ": ";cout << dx_dz[0] << "," << dx_dz[1] << " ";for (int i = 0; i < solution[0].size(); i++)cout << solution[0][i] << " ";cout << endl; }counting++;

        // ban found solution and shifted solutions
        map<State, int> banned;// <state,offset>
        if (rounds > 1) {
            for (int i = 0; i < rounds - 1; i++) {
                banned[solution[i]] = i * round_var_num + base_offset;
            }
        }
        else {
            banned[solution[0]] = base_offset;
        }
        ban_solution(Solver, banned);
    }
}

/**
 * @brief read cnf file into a vector
 * 
 * @param res vector<vector<int>>&
 * @param pathname const string
 */
void XoodooRound::read2Vector(vector<vector<int>>& res, const string pathname) {
    ifstream infile;
    infile.open(pathname.data());
    assert(infile.is_open());// fail to open
    vector<int> suanz;
    string s;

    res.clear();
    while (getline(infile, s)) {
        // each line is a vector<int> suanz
        istringstream is(s);
        int d;
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

/**
 * @brief read a solution from a readable solution file, read into States
 * 
 * @param res States&, output solution
 * @param infile ifstream&
 * @param mode int, 0 for XoodooSat, 1 for XooTools, 2 for single round
 * @param st const vector<string>&, default is {"ab","NE","bb"}, the start signs of a State
 * @return bool, whether found a solution in the input file
 */
bool XoodooRound::read2States(States& res, ifstream &infile, int mode, const vector<string> &st) {
    string s;
    string start = st[mode];

    res.clear();
    while (getline(infile, s)) {
        // read a line
        // if meet the start sign of a state
        if (s[0] == start[0] || s[0] == start[1]) {
            State temp;
            // read the state
            // a state is X lines
            for (int x = 0; x < X; x++) {
                s.clear();
                assert(getline(infile, s));
                // a line has Z chars
                // a char is a column
                for (int z = 0; z < Z; z++) {
                    if (s[z] != '.') {
                        int column = s[z] - '0';
                        // if the column(char) == 0, then print as '.'
                        for (int y = 0; y < Y; y++) {
                            if ((column >> y) & 0x1) temp.push_back(x*Z + y * X*Z + z);
                        }
                    }
                }
            }
            sort(temp.begin(), temp.end());
            res.push_back(temp);
        }
        // if len(res) == len(a solution), then return
        if (res.size() == AS_state_num && mode == 0) {
            return true;
        }
        if (res.size() == AS_state_num - 1 && mode == 1) {
            return true;
        }
        if (res.size() == 1 && mode == 2) {
            return true;
        }
        s.clear();
    }
    return false;
}

/**
 * @brief read a solution from a readable solution file, read into StateColumns
 * 
 * @param res StateColumns&, output solution
 * @param infile ifstream&
 * @param mode int, 0 for XoodooSat, 1 for XooTools, 2 for single round
 * @param st const vector<string>&, default is {"ab","NE","bb"}, the start signs of a State
 * @return int, the weight of the solution
 */
int XoodooRound::read2StateColumn(StateColumns& res, ifstream &infile, int mode, const vector<string> &st) {
    string s;
    string start = st[mode];
    int weight = 0;
    res.clear();
    while (getline(infile, s)) {
        // read a line
        // if meet the start sign of a state
        if (s[0] == start[0] || s[0] == start[1]) {
            StateColumn temp(X*Z, 0);
            // read the state
            // a state is X lines
            for (int x = 0; x < X; x++) {
                s.clear();
                assert(getline(infile, s));
                // a line has Z chars
                // a char is a column
                for (int z = 0; z < Z; z++) {
                    unsigned int column = 0;
                    // if the column(char) != 0, then weight++
                    if (s[z] != '.') {
                        column = s[z] - '0';
                        weight++;
                    }
                    temp[Z*x + z] = column;
                }
            }
            res.push_back(temp);
        }
        // if len(res) == len(a solution), then return
        if (res.size() == AS_state_num && mode == 0) {
            return weight;
        }
        if (res.size() == AS_state_num - 1 && mode == 1) {
            return weight;
        }
        if (res.size() == 1 && mode == 2) {
            return weight;
        }
        s.clear();
    }
    return 0;
}

/**
 * @brief read a solution from a readable solution file, read into vector<vector<Lit>>
 * 
 * @param assumptions vector<vector<Lit>>&
 * @param infile ifstream&
 * @param mode int, default is 2 for single round
 * @param offset int, var offset
 * @param st const vector<string>&, default is {"ab","NE","bb"}, the start signs of a State
 * @return bool, whether found a solution in the input file
 */
bool XoodooRound::read2assumptions(vector<vector<Lit>> &assumptions, ifstream &infile, int mode, int offset, const vector<string> &st) {
    string s;
    string start = st[mode];

    while (getline(infile, s)) {
        // read a line
        // if meet the start sign of a state
        if (s[0] == start[0] || s[0] == start[1]) {
            vector<Lit> temp;
            // read the state
            // a state is X lines
            for (int x = 0; x < X; x++) {
                s.clear();
                assert(getline(infile, s));
                // a line has Z chars
                // a char is a column
                for (int z = 0; z < Z; z++) {
                    // set column value
                    if (s[z] != '.') {
                        int column = s[z] - '0';
                        for (int y = 0; y < Y; y++) {
                            if ((column >> y) & 0x1) temp.push_back(Lit(x*Z + y * X*Z + z + offset, true));
                            else temp.push_back(Lit(x*Z + y * X*Z + z + offset, false));
                        }
                    }
                    else {
                        for (int y = 0; y < Y; y++) {
                            temp.push_back(Lit(x*Z + y * X*Z + z + offset, false));
                        }
                    }
                }
            }
            assumptions.push_back(temp);
        }
        s.clear();
    }
    return assumptions.size() > 0;
}

/**
 * @brief check if the solution is right
 * 
 * @param pathname const string
 * @param i int, check the ith trail
 * @param mode int, 0 for XoodooSat, 1 for XooTools, 2 for single round
 */
void XoodooRound::check_trails(const string pathname, int i, int mode) {
    ifstream infile(pathname.data(), ios::in);
    assert(infile.is_open());
    States res;
    int res_size; // the number of states in a solution

    if (mode == 0) res_size = AS_state_num;
    else if (mode == 1) res_size = AS_state_num - 1;
    else if (mode == 2) res_size = 1;

    int print_size; // the number of states to print
    if (mode == 0 || mode == 1) print_size = AS_state_num;
    else if (mode == 2) print_size = 1;

    // read the i-th solution
    for (int k = 0; read2States(res, infile, mode) && k <= i; k++) {
        if (res.size() != res_size) {
            cout << "wrong size: " << res.size() << endl;
            return;
        }
        // skip the 1st to (i-1)th solution
        if (k < i) continue;

        cout << "\ncheck the " << i << "th trail" << endl;
        // the solution is a1, b1, a2,...
        vector<tXoodooState> res_xoo;
        // first transform to tXoodooState
        for (int m = 0; m < res.size(); m++) {
            tXoodooState tmp;
            State2XooState(res[m], tmp);
            res_xoo.push_back(tmp);
        }
        // then run the lambda function to check if bi=lambda(ai)
        for (int m = 0; m < print_size - 1; m++) {
            tXoodooState tmp = lambda(res_xoo[m]);// DC
            /*tXoodooState tmp = transposelambda(res_xoo[m]);*/ // LC
            cout << "a" << m + 1 << ":" << endl;
            displayXooState(cout, res_xoo[m]);
            cout << endl;
            cout << "b" << m + 1 << ":" << endl;
            displayXooState(cout, tmp);
            cout << endl;
        }
        if (mode == 0) {
            cout << "check b" << print_size - 1 << ":" << endl;
            displayXooState(cout, res_xoo[res_xoo.size() - 1]);
            cout << endl;
        }
        if (mode == 2) {
            cout << "b" << AS_state_num - 1 << ":" << endl;
            displayXooState(cout, res_xoo[res_xoo.size() - 1]);
            cout << endl;
        }
        res.clear();
    }
    infile.close();
}

/**
 * @brief compare 2 solution file, check if their solutions are the same
 * 
 * @param path_res1 const string, file path 1
 * @param mode1 int, mode of file 1, 0 for XoodooSat, 1 for XooTools, 2 for single round
 * @param path_res2 const string, file path 2
 * @param mode2 int, mode of file 2, 0 for XoodooSat, 1 for XooTools, 2 for single round
 * @return vector<vector<int>>, output[0] logs the solutions in file1 not found in file2
                                output[1] logs the solutions in file2 not found in file1
 */
vector<vector<int>> XoodooRound::compare_trails(const string path_res1, int mode1, const string path_res2, int mode2) {
    ifstream in_res1(path_res1.data(), ios::in), in_res2(path_res2.data(), ios::in);
    States res1, res2;

    // get the number of solutions if file2
    int res2_solution_num = 0;
    while (read2States(res2, in_res2, mode2)) res2_solution_num++;
    in_res2.close();

    // if the i-th solution in file1 cannot be found in file2, then same_solution[i]=-1
    // same for inverse_same_solution
    map<int, int> same_solution, inverse_same_solution;
    for (int k = 0; k < res2_solution_num; k++) inverse_same_solution[k] = -1;
    vector<vector<int>> output(2);

    for (int i = 0; read2States(res1, in_res1, mode1); i++) {
        in_res2.open(path_res2.data(), ios::in);
        for (int j = 0; read2States(res2, in_res2, mode2); j++) {
            bool equal = true;
            auto dx_dz1 = genSmallestState(res1[0]);
            auto dx_dz2 = genSmallestState(res2[0]);
            // see if the two solutions are shift equal
            for (int k = 0; k < min(res1.size(), res2.size()); k++) equal &= (ShiftXZ(dx_dz1, res1[k]) == ShiftXZ(dx_dz2, res2[k]));
            same_solution[i] = -1;
            if (equal) {
                same_solution[i] = j;
                inverse_same_solution[j] = i;
                cout << i << " in " + path_res1 << " trail equal to the " << j << " trail in " + path_res2 << endl;
                break;
            }
        }
        in_res2.close();
    }
    in_res1.close();

    // logs the solutions in file1 not found in file2
    for (auto it = same_solution.begin(); it != same_solution.end(); it++) {
        if (it->second == -1) {
            cout << it->first << " trail in " + path_res1 + " not found equal" << endl;
            output[0].push_back(it->first);
        }
    }
    // logs the solutions in file2 not found in file1
    for (auto it = inverse_same_solution.begin(); it != inverse_same_solution.end(); it++) {
        if (it->second == -1) {
            cout << it->first << " trail in " + path_res2 + " not found equal" << endl;
            output[1].push_back(it->first);
        }
    }

    return output;
}

/**
 * @brief get AS and AS bound CNF
 * 
 * @param Obj vector<vector<int>>&
 * @param Solver SATSolver&
 * @param weight int
 * @param as_var_num int
 * @param as_mode int
 */
void XoodooRound::gen_obj_T(vector<vector<int>> &Obj, SATSolver &Solver, int weight, int as_var_num, int as_mode) {
    // var list for pysat.card function, see pysat_card_AS.py for more information
    string AS_cnf_num;
    // gen_extend_AS_cnf_num(AS_cnf_num);
    AS_cnf_num = to_string(as_var_num) + " ";
    string py_cmd = "python3 pysat_card_AS.py --bound_num " + to_string(weight) + " --var_num " + AS_cnf_num + " --mode " + to_string(as_mode);

    FILE *fp = NULL;
    char buf[128] = { 0 }, result[256] = { 0 };
    if ((fp = popen(py_cmd.data(), "r")) != NULL) {// run the py command and get print() output of the command
        while (fgets(buf, sizeof(buf), fp) != NULL) {
            strcat(result, buf);
        }
        pclose(fp);
        fp = NULL;
    }
    replace(AS_cnf_num.begin(), AS_cnf_num.end(), ' ', '_');

    int cnf_nvars = atoi(result);// the output is the number of variables of cnf
    Solver.new_vars(cnf_nvars);// add vars
    vector<string> mode = { "<=", ">=", "=" };// atmost, atleast, euqals
    read2Vector(Obj, objFilePath + "CNF_" + AS_cnf_num + "AS" + mode[as_mode] + to_string(weight) + ".txt");
}

/**
 * @brief add lambda relation to solver
 * 
 * @param Solver SATSolver&
 * @param rounds int
 * @param base_offset int
 */
void XoodooRound::add_lambda2solver(SATSolver &Solver, int rounds, int base_offset) {
    // lambda = RhoW(Theta(RhoE(·)))
    cout << "rhoE, inverse_rhoE, rhoW size: " << RhoE_Relation.size() << " " << inverse_RhoE_Relation.size() << " " << RhoW_Relation.size() << " " << Theta_Relation.size() << endl;
    // linear layer
    // a1->b1, a2->b2, ...
    cout << "start add linear" << endl;
    for (int round = 0; round < rounds - 1; round++) {
        for (int i = 0; i < var_num; i++) {
            vector<unsigned int> temp = {};

            // output bit position
            if(analysis_mode == 0) {
                // find the output index of RhoW
                temp.push_back(RhoW_Relation[Theta_Relation[RhoE_Relation[i]][0]] + var_num + round_var_num * round + base_offset);
            }
            else {
                temp.push_back(inverse_RhoE_Relation[transpose_Theta_Relation[inverse_RhoW_Relation[i]][0]] + var_num + round_var_num * round + base_offset);
            }

            // input bit position
            for (int j = 1; j < Y * 2 + 2; j++) {
                if(analysis_mode == 0) {
                    // find the input index of Theta, then transform back to the input of RhoE
                    temp.push_back(inverse_RhoE_Relation[Theta_Relation[RhoE_Relation[i]][j]] + round_var_num * round + base_offset);
                }
                else {
                    temp.push_back(RhoW_Relation[transpose_Theta_Relation[inverse_RhoW_Relation[i]][j]] + round_var_num * round + base_offset);
                }
            }
            Solver.add_xor_clause(temp, false);
        }
    }
}

/**
 * @brief add chi relation to solver
 * 
 * @param Solver SATSolver&
 * @param rounds int
 * @param base_offset int
 */
void XoodooRound::add_chi2solver(SATSolver &Solver, int rounds, int base_offset) {
    /*
    chi
    generate bi --> a_i+1
    */
    cout << "start add chi" << endl;
    for (int round = 0; round < rounds - 2; round++) {
        for (int x = 0; x < X; x++) {
            for (int z = 0; z < Z; z++) {
                for (int ithClause = 0; ithClause < Chi_Relation.size(); ithClause++) {
                    vector<Lit> clause;// one of the 6 clauses of a colliding Sbox
                    for (int ithVar = 0; ithVar < Y; ithVar++) {// input, 0 - Y-1
                        if (Chi_Relation[ithClause][ithVar] != "0") {
                            if (Chi_Relation[ithClause][ithVar] == "T") {
                                clause.push_back(Lit((ithVar*X*Z + z + Z * x + base_offset + round_var_num * round), true));
                            }
                            else {
                                clause.push_back(Lit((ithVar*X*Z + z + Z * x + base_offset + round_var_num * round), false));
                            }
                        }
                    }

                    for (int ithVar = Y; ithVar < Y * 2; ithVar++) {// output, Y - 2Y-1
                        if (Chi_Relation[ithClause][ithVar] != "0") {
                            if (Chi_Relation[ithClause][ithVar] == "T") {
                                clause.push_back(Lit(((ithVar - Y)*X*Z + z + Z * x + var_num + base_offset + round_var_num * round), true));
                            }
                            else {
                                clause.push_back(Lit(((ithVar - Y)*X*Z + z + Z * x + var_num + base_offset + round_var_num * round), false));
                            }
                        }
                    }
                    Solver.add_clause(clause);
                    clause.clear();
                }
            }
        }
    }
}

/**
 * @brief add AS relation to solver
 * 
 * @param Solver SATSolver&
 * @param rounds int
 * @param base_offset int
 * @param core_var_num int
 */
void XoodooRound::add_AS2solver(SATSolver &Solver, int rounds, int base_offset, int core_var_num) {
    cout << "start add AS" << endl;
    // a1, a2, b2, ...
    for (int x = 0; x < X; x++) {
        for (int z = 0; z < Z; z++) {
            // for each Sbox
            for (int ithClause = 0; ithClause < AS_Relation.size(); ithClause++) {
                vector<vector<Lit>> clause(rounds);
                for (int ithVar = 0; ithVar < Y; ithVar++) {
                    if (AS_Relation[ithClause][ithVar] != "0") {
                        for (int i = 0; i < rounds; i++) {
                            int offset = i * round_var_num;// a1, a2, b2, ...
                            if (i == rounds - 1 && offset) offset -= var_num;
                            if (AS_Relation[ithClause][ithVar] == "T") {
                                clause[i].push_back(Lit((ithVar*X*Z + z + Z * x + offset + base_offset), true));
                            }
                            else {
                                clause[i].push_back(Lit((ithVar*X*Z + z + Z * x + offset + base_offset), false));
                            }
                        }
                    }
                }

                if (AS_Relation[ithClause][Y] != "0") {
                    for (int i = 0; i < rounds; i++) {
                        // AS_a1, AS_a2, AS_b2, ...
                        int offset = i * AS_node_var_num;
                        if (AS_Relation[ithClause][Y] == "T") {
                            clause[i].push_back(Lit((core_var_num + z + Z * x + offset), true));
                        }
                        else {
                            clause[i].push_back(Lit((core_var_num + z + Z * x + offset), false));
                        }
                    }
                }
                for (int i = 0; i < rounds; i++) {
                    Solver.add_clause(clause[i]);
                    clause[i].clear();
                }
            }
        }
    }
}

/**
 * @brief add AS bound relation to solver
 * 
 * @param Solver SATSolver&
 * @param Obj vector<vector<int>>&
 * @param core_var_num int
 */
void XoodooRound::add_ASobj2solver(SATSolver &Solver, vector<vector<int>> &Obj, int core_var_num) {
    // AS_a1a2b2
    cout << "start add AS obj" << endl;
    for (int ithClause = 0; ithClause < Obj.size(); ithClause++) {
        // one of the 10 clauses of a colliding Sbox
        vector<Lit> clause;

        for (int ithVar = 0; ithVar < Obj[ithClause].size(); ithVar++) {
            int temp = Obj[ithClause][ithVar];
            if (temp > 0) {
                clause.push_back(Lit(abs(temp) - 1 + core_var_num, false));
            }
            else {
                clause.push_back(Lit(abs(temp) - 1 + core_var_num, true));
            }
        }
        Solver.add_clause(clause);
        clause.clear();
    }
}


/**
 * @brief differential/linear analysis round function
        trail: a0 -> b0 -> a1 -> b1 -> a2 -> b2 -> a3
        trail core: a1 -> b1 -> a2 -> b2 -> a3
 * 
 */
void XoodooRound::XoodooRound_AS() {
    add_lambda2solver(solver, round_num, 0);
    add_chi2solver(solver, round_num, var_num);
    add_AS2solver(solver, AS_state_num, 0, var_num*core_state_num);
    add_ASobj2solver(solver, obj, var_num*core_state_num);
}

/**
 * @brief differential/linear analysis round function, with trail extension
            ONLY extend 1 round(forward/backward)
 * 
 */
void XoodooRound::extendRound_AS() {
    // print date and time
    time_point<system_clock> start = system_clock::now();
    auto st = system_clock::to_time_t(start);
    struct tm* stm = localtime(&st);
    cout << stm->tm_mon + 1 << " " << stm->tm_mday << " " << stm->tm_hour << ":" << stm->tm_min << ":" << stm->tm_sec << endl;

    ofstream out;
    string res_prefix = objFilePath + "result/";
    if (access(res_prefix.data(), F_OK) == -1) {//check if result dir is existed
        mkdir(res_prefix.data(), S_IRWXO | S_IRWXG | S_IRWXU);
    }

    vector<vector<Lit>> assumptions;
    string pre_path = objFilePath + "extra_trails.txt";
    ifstream in_pre(pre_path.data(), ios::in);
    if (!read2assumptions(assumptions, in_pre, 2)) {
        cout << "no assumptions" << endl;
        return;
    }

    SATSolver extend_solver;
    int extend_AS_state_num = 1;
    int extend_core_var_num = var_num + round_var_num;
    int extend_AS_var_num = X * Z*extend_AS_state_num;

    int previous_weight = 24;
    int total_weight = 36;
    extend_solver.set_num_threads(thread_num);
    extend_solver.new_vars(extend_core_var_num);//add vars
    vector<vector<int>> extend_obj;
    gen_obj_T(extend_obj, extend_solver, total_weight - previous_weight, extend_AS_var_num, AS_mode);

    add_chi2solver(extend_solver, extend_AS_state_num + 2, 0);
    add_lambda2solver(extend_solver, extend_AS_state_num + 1, var_num);
    add_AS2solver(extend_solver, extend_AS_state_num, round_var_num, extend_core_var_num);
    add_ASobj2solver(extend_solver, extend_obj, extend_core_var_num);

    vector<string> mode = { "<=", ">=", "=" };//atmost, atleast, euqals
    string res_filepath = res_prefix + "xoodoo_result_4R_#AS" + mode[AS_mode] + to_string(total_weight) + ".txt";
    out.open(res_filepath, ios::ate);//clear result txt first
    out.close();
    cout << "start solving" << endl;
    if (assumptions.size())
        solve_and_output(extend_solver, res_filepath, extend_AS_state_num, round_var_num, extend_core_var_num, assumptions[0]);
}

/**
 * @brief main function, including time log and main attack fucntion
 * 
 */
void XoodooRound::main() {
    // print date and time
    time_point<system_clock> start = system_clock::now();
    auto st = system_clock::to_time_t(start);
    struct tm* stm = localtime(&st);
    cout << stm->tm_mon + 1 << " " << stm->tm_mday << " " << stm->tm_hour << ":" << stm->tm_min << ":" << stm->tm_sec << endl;

    // round function, generate clauses for solver
    XoodooRound_AS();

    // mkdir result and the empty result file
    ofstream out;
    string res_prefix = objFilePath + "result/";
    if (access(res_prefix.data(), F_OK) == -1) {//check if result dir is existed
        mkdir(res_prefix.data(), S_IRWXO | S_IRWXG | S_IRWXU);
    }

    vector<string> mode = { "<=", ">=", "=" };//atmost, atleast, euqals
    string res_filepath = res_prefix + "xoodoo_result_" + to_string(round_num) + "R_#AS" + mode[AS_mode] + to_string(AS_weight_num) + ".txt";
    out.open(res_filepath, ios::ate);//clear result txt first
    out.close();

    // ban previous solution and solve phase
    cout << "start banning found trails with weight" + mode[AS_mode] << AS_weight_num << endl;
    string ban_pre_path = objFilePath + "found" + mode[AS_mode] + to_string(AS_weight_num) + ".txt";
    ifstream in_ban_pre(ban_pre_path.data(), ios::in);
    States ban_pre;
    while (read2States(ban_pre, in_ban_pre, 0)) {
        map<State, int> banned;
        for (int i = 0; i < ban_pre.size() - 1; i++) {
            banned[ban_pre[i]] = i * round_var_num;
        }
        ban_solution(solver, banned);
    }
    in_ban_pre.close();

    cout << "start solving" << endl;
    solve_and_output(solver, res_filepath, round_num, 0, var_num*core_state_num);

}

/**
 * @brief generate DDT table for chi
 * 
 * @param ddt unordered_map<int, set<int>>&, output ddt table
 */
void XoodooRound::gen_chi_DDT(unordered_map<int, set<int>> &ddt) {
    // input of chi from 0 to 2^Y - 1
    int cnt = pow(2, Y);
    // chi map
    unordered_map<int, int> chi;

    // first calculate the chi map
    for (int i = 0; i < cnt; ++i) {
        if (i == 0) chi[i] = 0;
        else {
            // the input is i
            vector<int> A(Y), tmpa(Y);
            int y, out = 0;
            // convert the input i to Bitmap A
            for (y = 0; y < Y; ++y) {
                A[y] = (i >> y) & 0x1;
            }
            // part 1 of the chi function
            for (y = 0; y < Y; y++) {
                tmpa[y] = A[y] ^ ((~A[(y + 1) % Y]) & A[(y + 2) % Y]);
            }
            for (y = 0; y < Y; y++) {
                A[y] = tmpa[y];
            }
            // part 2 of the chi function
            for (y = 0; y < Y; ++y) {
                out += (A[y] << y);
            }
            chi[i] = out;
        }
    }

    // the input of ddt is input xor(delta)
    for (int delta = 0; delta < cnt; ++delta) {
        if (delta == 0) ddt[delta].insert(0);
        else {
            for (int a = 0; a < cnt; ++a) {
                int b = a ^ delta;
                // get all possible output xor
                int out = chi[a] ^ chi[b];
                ddt[delta].insert(out);
            }
        }
    }
}

/**
 * @brief extend a round by dfs
 * 
 * @param assumption StateColumn&, solver assumption
 * @param nonzero_col_index vector<int>&, non-zero columns
 * @param ddt unordered_map<int, set<int>>&, ddt for chi
 * @param temp StateColumn, dfs parameter
 * @param i int, dfs parameter
 * @param n int, number of non-zero columns
 * @param total_weight int&, the weight bound
 * @param pre_weight int&, previous 2-round trail core weight
 * @param forward bool, forward or backward extension
 */
void XoodooRound::dfs(StateColumn &assumption, vector<int> &nonzero_col_index, unordered_map<int, set<int>> &ddt, StateColumn temp, int i, int n, int &total_weight, int &pre_weight, bool forward) {
    if (i == n) {
        tXoodooState xooA;
        StateColumn2XooState(temp, xooA);
        if (forward) xooA = lambda(xooA);
        else xooA = inverselambda(xooA);

        int wt = caculateXooStateWeight(xooA);
        if (wt <= total_weight - pre_weight) {
            cout << "pre_weight: " << pre_weight << endl;
            cout << "sum_weight: " << pre_weight + wt << endl;
            displayXooState(cout, xooA);
            cout << endl;
        }
    }
    else {
        int index = nonzero_col_index[i];
        int col = assumption[index];
        for (int num : ddt[col]) {
            temp[index] = num;
            dfs(assumption, nonzero_col_index, ddt, temp, i + 1, n, total_weight, pre_weight, forward);
        }
    }
}

// trail analysis with trail extension
void XoodooRound::extend_main(bool forward) {
    // extendRound_AS();

    // print date and time
    time_point<system_clock> start = system_clock::now();
    auto st = system_clock::to_time_t(start);
    struct tm* stm = localtime(&st);
    cout << stm->tm_mon + 1 << " " << stm->tm_mday << " " << stm->tm_hour << ":" << stm->tm_min << ":" << stm->tm_sec << endl;

    string res_prefix = objFilePath + "result/";
    if (access(res_prefix.data(), F_OK) == -1) {// check if result dir is existed
        mkdir(res_prefix.data(), S_IRWXO | S_IRWXG | S_IRWXU);
    }

    StateColumns assumptions;
    int pre_mode = 0;
    string pre_path = objFilePath + "extra_trails.txt";
    ifstream in_pre(pre_path.data(), ios::in);

    int total_weight = AS_weight_num;
    vector<string> mode = { "back","forw" };
    string res_filepath = res_prefix + "xoodoo_result_4R_ext_" + mode[forward] + "#AS<=" + to_string(total_weight) + ".txt";
    // ofstream out;
    // out.open(res_filepath, ios::ate);// clear result txt first
    // out.close();
    cout << "start extending " + mode[forward] << endl;
    unordered_map<int, set<int>> ddt;
    gen_chi_DDT(ddt);

    int pre_weight, thres_weight = 15;
    if (pre_mode == 0 && pre_path != objFilePath + "extra_trails.txt") pre_weight = read2StateColumn(assumptions, in_pre, pre_mode);//first solution is 0
    for (int cnt = 1;pre_weight = read2StateColumn(assumptions, in_pre, pre_mode); ++cnt) {
        StateColumn assumption;
        if(forward) assumption = assumptions.back();
        else assumption = assumptions.front();
        if (pre_mode == 1) {
            tXoodooState tmp;
            StateColumn2XooState(assumptions.back(), tmp);
            tmp = lambda(tmp);
            if(forward) XooState2StateColumn(tmp, assumption);
            pre_weight += caculateXooStateWeight(tmp);
        }
        
        StateColumn temp = assumption;
        if (caculateStateColumnWeight(temp) >= thres_weight) {
            cout << cnt << " weight>=" << to_string(thres_weight) << endl;
            continue;
        }

        vector<int> nonzero_col_index;// non-zero column index
        for (int x = 0; x < X; ++x) {
            for (int z = 0; z < Z; ++z) {
                int index = Z * x + z;
                int col = assumption[index];
                if (col) {
                    nonzero_col_index.push_back(index);
                }
            }
        }
        cout << "start dfs " << cnt << endl;
        dfs(assumption, nonzero_col_index, ddt, temp, 0, nonzero_col_index.size(), total_weight, pre_weight, forward);
    }
}