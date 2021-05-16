#include <iostream>
#include <stdlib.h>
#include <chrono>
#include <assert.h>
#include <getopt.h>
#include "xoodooRound.h"

using std::cout;
using std::endl;
using XOODOOSAT::XoodooRound;

void show_usage(char* cmd) {
    cout << cmd << " [Options]" <<endl;
    cout << "Options:" << endl;
    cout << " -r, --round round_num      How many rounds to trail (default 4)." << endl;
    cout << " -w, --weight weight        The weight to bound for the trail (default 27)." << endl;
    cout << " -t, --thread thread        The number of threads for the process (default 16)." << endl;
    cout << " -m, --mode mode            The mode for weight sum (default 0), choices={0,1,2}, 0 for atmost, 1 for atleast, 2 for equals." << endl;
    cout << "                            See pysat_card_AS.py for more information." << endl;
    cout << " -h, --help                 Help information." << endl;
}

int main(int argc, char* argv[]) {
    int round_num = 4, weight = 27, thread_num = 16, mode = 0;//default: 4 rounds, max weight 27
    int c;
    bool success = true;

    while(1) {
        int opt_index = 0;
        static struct option long_options[] =
        {
            {"round", required_argument, NULL, 'r'},
            {"weight", required_argument, NULL, 'w'},
            {"thread", required_argument, NULL, 't'},
            {"mode", required_argument, NULL, 'm'},
            {"help", no_argument, NULL, 'h'}
        };
        c = getopt_long(argc, argv, "r:w:t:m:h", long_options, &opt_index);
        if(c == -1)
            break;

        switch(c) {
            case 'r':
                round_num = atoi(optarg);
                break;
            case 'w':
                weight = atoi(optarg);
                break;
            case 't':
                thread_num = atoi(optarg);
                break;
            case 'm':
                mode = atoi(optarg);
                break;
            case 'h':
                show_usage(argv[0]);
                return 0;
            default:
                success = false;
                break;
        }
    }
    if(!success) {
        show_usage(argv[0]);
        return 1;
    }
    if (optind < argc) {
        cout<<"non-option argv-elements: ";
        while (optind < argc)
            cout<<argv[optind++];
        cout<<"\n"<<endl;
        show_usage(argv[0]);
        return 1;
    }
    
    if(argc == 1) {
        show_usage(argv[0]);
        cout<<"\nnow using default args"<<endl;
    }
    cout<<round_num<<" rounds, "<<"weight "<<weight<<", "<<thread_num<<" threads"<<", in mode "<<mode<<endl;
    XoodooRound xoo(round_num,weight,thread_num,mode);

    /*auto diff = xoo.compare_trails("./test.txt","./DC.txt");
    for(int i=0; i<diff[1].size(); i++) {
        cout<<diff[1][i]<<" trail in DC"<<endl;
        xoo.check_trails("./DC.txt",diff[1][i],1);
    }*/
    //xoo.check_trails("./extra_trails.txt",1,0);
    //xoo.check_trails("./test.txt",70,0);
    xoo.XoodooRound_AS();

    return 0;
}