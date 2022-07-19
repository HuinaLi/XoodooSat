#include <iostream>
#include <stdlib.h>
#include <chrono>
#include <assert.h>
#include <getopt.h>
#include <vector>
#include "xoodooRound.h"

using std::cout;
using std::endl;
using XOODOOSAT::XoodooRound;

void show_usage(char* cmd) {
    cout << cmd << " [Options]" << endl;
    cout << "Options:" << endl;
    cout << " -a, --analysis analysis_mode       0 for differential, 1 for linear analysis (default 0)." << endl;
    cout << " -r, --round round_num              How many rounds to trail (default 3)." << endl;
    cout << " -w, --weight weight                The weight to bound for the trail (default 25)." << endl;
    cout << " -t, --thread thread                The number of threads for the process (default 16)." << endl;
    cout << " -m, --mode mode                    The mode for weight sum (default 0), choices={0,1,2}, 0 for atmost, 1 for atleast, 2 for equals." << endl;
    cout << "                                    See pysat_card_AS.py for more information." << endl;
    cout << " -h, --help                         Help information." << endl;
}

int main(int argc, char* argv[]) {
    int analysis_mode = 0, round_num = 3, weight = 25, thread_num = 16, mode = 0;// default: 3 rounds, max weight 25, differential analysis
    int c;
    bool success = true;

    // process command line parameters
    while(1) {
        int opt_index = 0;
        static struct option long_options[] =
        {
            {"analysis", required_argument, NULL, 'a'},
            {"round", required_argument, NULL, 'r'},
            {"weight", required_argument, NULL, 'w'},
            {"thread", required_argument, NULL, 't'},
            {"mode", required_argument, NULL, 'm'},
            {"help", no_argument, NULL, 'h'}
        };
        c = getopt_long(argc, argv, "a:r:w:t:m:h", long_options, &opt_index);
        if(c == -1)
            break;

        switch(c) {
            case 'a':
                analysis_mode = atoi(optarg);
                break;
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

    vector<string> analysis_print = {"differential trail analysis ", "linear trail analysis"};
    cout << analysis_print[analysis_mode] << endl;
    cout << round_num << " rounds, " << "weight " << weight << ", " << thread_num << " threads" << ", in mode " << mode << endl;
    XoodooRound xoo(analysis_mode, round_num, weight, thread_num, mode);

    /*
    // compare solution results of XoodooSat and XooTool
    auto diff = xoo.compare_trails("./dc_3_23.txt",0,"./DC23.txt",1);
    for(int i=0; i<diff[0].size(); i++) {
        cout<<diff[0][i]<<" trail in dc_3_23"<<endl;
        xoo.check_trails("./dc_3_23.txt",diff[0][i],0);
    }
    */
    // check if the trails in solution file is correct
    // xoo.check_trails("./extra_trails.txt", 1, 0);
    // xoo.check_trails("./test.txt", 70, 0);

    // run the attack
    xoo.main();

    // xoo.extend_main(true); // forward
    // xoo.extend_main(false);// backward

    return 0;
}
