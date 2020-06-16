#include <algorithm>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include "rwr.h"
#include <unordered_set>    
#include <cstdlib>
#include <cstring>


void usage() {
    cerr << "ResAcc [-d <dataset>] [-n node_count (default 20)] [-alpha alpha  (default 0.2)] [-eps epsilon (default 0.5)] [-h h_HopFWD (default 2)] [-r rmax_hop (default 0.1)] [-rf rmax_f (default 0.1)] " << endl;
}

int check_inc(int i, int max) {
    if (i == max) {
        usage();
        exit(1);
    }
    return i + 1;
}

bool maxCmp(const pair<int, double>& a, const pair<int, double>& b){
    return a.second > b.second;
}

int main(int argc, char *argv[]){
    int i = 1;
    char *endptr;
    string filename; //the name of dataset
    double alpha = 0.2;            //restart probability
    int node_count = 20;           //query node size
    double epsilon = 0.5; //relative error

    int h = 2;//parameter h: h-HopFWD; the number of hops in the h-HopFWD phase
    double r_max_hop = 0.1; // rmax_hop * 10^{-8}
    double r_max_f = 0.1; // rmax_f * 10^{-8}

    string algo = "RESACC";
    if(argc < 6){
        usage();
        exit(1);
    }
    while (i < argc) {
        if (!strcmp(argv[i], "-d")) {
            i = check_inc(i, argc);
            filename = argv[i];
        } 
        else if (!strcmp(argv[i], "-algo")) {
            i = check_inc(i, argc);
            algo = argv[i];
        }
        else if (!strcmp(argv[i], "-h")) {
            i = check_inc(i, argc);
            h = strtod(argv[i], &endptr);
            if ((h < 0) && endptr) {
                cerr << "Invalid h argument" << endl;
                exit(1);
            }
        }
        else if (!strcmp(argv[i], "-n")) {
            i = check_inc(i, argc);
            node_count = strtod(argv[i], &endptr);
            if ((node_count < 0) && endptr) {
                cerr << "Invalid node_count argument" << endl;
                exit(1);
            }
        }
        else if (!strcmp(argv[i], "-r")) {
            i = check_inc(i, argc);
            r_max_hop = strtod(argv[i], &endptr);
        }else if (!strcmp(argv[i], "-rf")) {
            i = check_inc(i, argc);
            r_max_f = strtod(argv[i], &endptr);
        }
        else if (!strcmp(argv[i], "-eps")) {
            i = check_inc(i, argc);
            epsilon = strtod(argv[i], &endptr);
            if (((epsilon < 0) || (epsilon > 1)) && endptr) {
                cerr << "Invalid rmax argument" << endl;
                exit(1);
            }
        }
        else if (!strcmp(argv[i], "-alpha")) {
            i = check_inc(i, argc);
            alpha = strtod(argv[i], &endptr);
            if (((alpha < 0) || (alpha > 1)) && endptr) {
                cerr << "Invalid rmax argument" << endl;
                exit(1);
            }
        }
        else {
            usage();
            exit(1);
        }
        i++;
    }
    
    RWR rwr = RWR(filename, alpha);
    if(algo == "GEN_QUERY"){
        ofstream outFile("dataset/" + filename + "/" +filename + ".query");
        rwr.generateQueryNode(node_count, outFile);
        outFile.close(); 
    }
    else if(algo == "TO_UNDIRECTED"){
        rwr.to_undirected_graph();
    }
    else if(algo == "GEN_GROUND_TRUTH") {
        string queryname = "dataset/" + filename + "/" +filename + ".query";
        if (!rwr.is_file_exist(queryname)) {
            cout << "please generate query file first" << endl;
        } else
            rwr.PowerMethodMulti(100, node_count, 10);/*  多线程PowerMethparameter: iteration loops, node size, thread num */
    }
    else if(algo == "RESACC"){
        string queryname = "dataset/" + filename + "/" +filename + ".query";
        if(!rwr.is_file_exist(queryname)){
            cout << "please generate query file first" << endl;
            return 0;
        }
        ifstream nodes_file("dataset/" + filename + "/" +filename + ".query");
        vector<int> test_nodes;
        while(!nodes_file.eof()){
            int temp_node;
            nodes_file >> temp_node;
            test_nodes.push_back(temp_node);
        }
        cout << "read done!" << endl;

        int h_hops = h;

        for(int t = 0; t < node_count; t++) {
            //cout << "current ite: " << t << endl;
            int test_node = test_nodes[t];
            cout << "node: " << test_node << " " << epsilon << " " << h_hops << endl;

            cout<< "---------------" <<endl;
            rwr.resacc(test_node, h_hops, r_max_hop, r_max_f, epsilon);

        }
        cout << "avg time: " << (rwr.avg_khopfwd_time  + rwr.avg_omfwd_time + rwr.avg_walk_time) / (double) node_count << endl;
        cout << "avg h-HopFWD time: " << rwr.avg_khopfwd_time / (double) node_count << endl;
        cout << "avg OMFWD time: " << rwr.avg_omfwd_time / (double) node_count << endl;
        cout << "avg remedy time: " << rwr.avg_walk_time / (double) node_count << endl;

        stringstream ss;
        ss << "estimated_rwr/" << filename << "/resacc.time";
        string outfile = ss.str();
        ofstream timefile(outfile);
        timefile <<"total time: " << (rwr.avg_khopfwd_time  + rwr.avg_omfwd_time + rwr.avg_walk_time) / (double) node_count <<endl;
        timefile <<"avg h-HopFWD time: " << rwr.avg_khopfwd_time/ (double) node_count <<endl;
        timefile <<"avg OMFWD time: " << rwr.avg_omfwd_time/ (double) node_count <<endl;
        timefile <<"avg remedy time: " << rwr.avg_walk_time/ (double) node_count <<endl;

        timefile.close();

        rwr.avg_khopfwd_time = 0;
        rwr.avg_omfwd_time = 0;
        rwr.avg_walk_time = 0;
        rwr.avg_time = 0;
    }
    return 0;
};
