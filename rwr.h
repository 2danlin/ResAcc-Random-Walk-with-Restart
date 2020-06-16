#ifndef RWR_H
#define RWR_H

#include <vector>
#include <algorithm>
#include <queue>
#include <functional>
#include <iostream>
#include <fstream>
#include <future>
#include <string>
#include <sstream>
#include "Graph.h"
#include "Random.h"
#include "alias.h"
#include <unordered_map>
#include <unordered_set>
#include <thread>
#include <sys/time.h> 
#include <time.h>
#include <cmath>
#include <limits>



bool maxScoreCmp(const pair<int, double>& a, const pair<int, double>& b){
    return a.second > b.second;
}

class pqcompare
{
  bool reverse;
public:
  pqcompare(const bool& revparam=false)
    {reverse=revparam;}
  bool operator() (const pair<int, double>& lhs, const pair<int, double>&rhs) const
  {
    if (reverse) return (lhs.second > rhs.second);
    else return (lhs.second < rhs.second);
  }
};


void RandomWalk(int walk_num, Alias &alias, Random &R, Graph& g, int* vert_count){
    for(int i = 0; i < walk_num; i++){
        int tempNode = alias.generateRandom_t(R);
        vert_count[tempNode]++;
        while(R.drand_t() > 0.2){
            int length = g.getOutSize(tempNode);
            if(length > 0){   
                int r = R.generateRandom_t() % length;
                tempNode = g.getOutVert(tempNode, r);
            }
            vert_count[tempNode]++;
        }
    }
}

class RWR
{
friend void rwr_t_PowerMethod(RWR* ppr, vector<int> nodeList, int iterations);
public:
    double avg_pre;
    double avg_recall;
    double avg_time;
    double avg_khopfwd_time;
    double avg_omfwd_time;
    double avg_walk_time;
    double avg_NDCG;

    Graph g;
    Random R;
    int vert;
    double alpha;
    string target_filename;
    double* vert_count;
    int* value_verts;
    double* r_t;
    double* r;
    unsigned NUM_CORES;
    int** multiVertCount;
    double* resultList;

    void PowerMethodMulti(int iterations, int node_count, int num_thread);
    void PowerMethodBackwardMulti(int iterations, int node_count, int num_thread);
    //const static int NUMTHREAD = 20;
    Random* Rs;

    RWR(string dataset_name, double input_alpha) {
        alpha = input_alpha;
        avg_time = 0;
        avg_khopfwd_time =0;
        avg_omfwd_time = 0;
        avg_walk_time = 0;
        avg_pre = 0 ;
        avg_recall = 0;
        avg_NDCG = 0;

        target_filename = dataset_name;
        string filename = "dataset/" + dataset_name + "/" + dataset_name + ".txt";
        g.inputGraph(filename);
        cout << "edge num: " << g.m << endl;
        vert = g.n;
        srand(unsigned(time(0)));
        R = Random(unsigned(rand()));
        vert_count = new double[vert];
        resultList = new double[vert];
        r = new double[vert];
        value_verts = new int[vert];

        for(int i =0 ; i < vert; i++){
            resultList[i] = 0;
            vert_count[i] = 0;
            value_verts[i] = -1;
            r[i] = 0;
        }
        NUM_CORES = std::thread::hardware_concurrency();
        assert(NUM_CORES >= 2);
        cout << "thread core: " << NUM_CORES << endl;
        multiVertCount = new int*[NUM_CORES];
        Rs = new Random[NUM_CORES];
        for(int i = 0; i < NUM_CORES; i++){
            Rs[i] = Random(unsigned(rand()));
            multiVertCount[i] = new int[vert];
            for(int j = 0; j < vert; j++){
                multiVertCount[i][j] = 0;
            }
        }
        cout << "init done! " << endl;
    }
    ~RWR() {
        for(int i = 0; i < NUM_CORES; i++){
            delete[] multiVertCount[i];
        }

        delete[] multiVertCount;
        delete[] vert_count;
        delete[] value_verts;
        delete[] r;
        delete[] Rs;
    }

    bool is_file_exist(string fileName)
    {
        ifstream infile(fileName);
        return infile.good();
    }

    //取s点的groundtruth
    vector<int> getRealTopK(int s, int k){
        stringstream ss;
        ss << "real_rwr/" << target_filename << "/" << s << ".txt";
        string infile = ss.str();
        ifstream real(infile);
        vector<int> realList;
        vector<double> simList;
        for(int i = 0; i < vert; i++){
            int tempId;
            double tempSim;
            real >> tempId >> tempSim;
            if(i >= k && tempSim < simList[k-1]){
                break;
             }
            if( i == 0)
                continue;
             realList.push_back(tempId);
             simList.push_back(tempSim);
         }
         real.close();
         return realList;
     }

    unordered_map<int, double> getRealTopKMap(int s, int k){
        unordered_map<int, double> answer_map;
        stringstream ss;
        ss << "real_rwr/" << target_filename << "/" << s << ".txt";
        string infile = ss.str();
        ifstream real(infile);
        double k_Sim = 0;
        for(int i = 0; i < vert; i++){
            int tempId;
            double tempSim;
            real >> tempId >> tempSim;
            if(i == k - 1){
                k_Sim = tempSim;
            }
            if(i >= k && tempSim < k_Sim){
                break;
            }
            answer_map[tempId] = tempSim;
        }
        real.close();
        return answer_map;
    }

    unordered_map<int,double> getRealRWRs(int sourceNode){
        unordered_map<int, double> real_rwr;

        stringstream ss;
        ss << "real_rwr/" << target_filename << "/" << sourceNode <<".txt";
        string infile = ss.str();
        ifstream real(infile);
        for(int i =0; i < vert; i++){
            int tempId;
            double tempSim;
            real >> tempId >> tempSim;
            real_rwr[tempId] =tempSim;
        }
        real.close();

        return real_rwr;
    }
    double getRealRWR(int sourceNode, int targetNode){
        double real_rwr;

        stringstream ss;
        ss << "real_rwr/" << target_filename << "/" << sourceNode << ".txt";
        string infile = ss.str();
        ifstream real(infile);

        for(int i =0; i < vert; i++){
            int tempId;
            double tempSim;
            real >> tempId >> tempSim;
            if(tempId == targetNode){
                real_rwr = tempSim;
                break;
            }
        }
        real.close();

        return real_rwr;
    }

    //compute the number of dangling nodes
    void getNumDeadends(){
        int num_dde = 0;
        for(int i=0; i < vert; i++){
            if(g.getOutSize(i) == 0)
                num_dde += 1;
        }
        cout<<"num of dangling nodes: " << num_dde <<endl;
        cout<<"percentage of dangling nodes: " << num_dde /(double) vert << endl;
    }

    double* MonteCarlo(int u, int k, double walk_num){
        vector<int> realList = getRealTopK(u, k);
        unordered_map<int, double> realMap = getRealTopKMap(u, k);
        double* resultList = new double[vert];
        for(int i = 0; i < vert; i++)
            resultList[i] = 0;
        if(g.getOutSize(u) == 0){
            resultList[u] = alpha;
            return resultList;
        }

        clock_t t0 = clock();
        for(double i = 0; i < walk_num; i++){
            int tempNode = u;
            resultList[tempNode] += alpha;
            while(R.drand() > alpha){
                int length = g.getOutSize(tempNode);
                if(length == 0)
                    tempNode = u;
                else{
                    int r = R.generateRandom() % length;
                    tempNode = g.getOutVert(tempNode, r);
                }
                resultList[tempNode] += alpha;
            }
        }
        clock_t t1 = clock();
        cout << "MonteCarlo time: " << (t1 - t0) / (double) CLOCKS_PER_SEC << endl;
        for(int i = 0; i < vert; i++){
            resultList[i] /= (double) walk_num;
        }
        for(int i = 0; i<k; i++){
            int tempNode = realList[i];
            //avg_top_err[i] += abs(realMap[tempNode] - resultList[tempNode]);
        }

        vector<int> newLeftVec;
        typedef priority_queue<pair<int, double>, vector<pair<int,double> >, pqcompare> pq;
        pq upper_pq(pqcompare(true));

        double UpperBound = 0;
        for(int i = 0; i < vert; i++){
            if(i < k){
                upper_pq.push(pair<int, double>(i, resultList[i]));
                if(i == k - 1)
                    UpperBound = upper_pq.top().second;
            }
            else{
                if(resultList[i] > UpperBound){
                    upper_pq.pop();
                    upper_pq.push(pair<int, double>(i, resultList[i]));
                    UpperBound = upper_pq.top().second;
                }
            }
        }
        for(int i = 0; i < k; i++){
            newLeftVec.push_back(upper_pq.top().first);
            upper_pq.pop();
        }

        cout << "precision: " << calPrecision(newLeftVec, realList, k) << endl;
        avg_pre += calPrecision(newLeftVec, realList, k);
        avg_time +=  (t1 - t0) / (double) CLOCKS_PER_SEC ;
        return resultList;
    }

    vector<int> getTopKList(double* resultList, int u, int k){
        vector<int> topKList;
        vector<pair<int,double> > sim;
        for(int i = 0; i < vert; i++){
            sim.push_back(pair<int, double>(i, resultList[i]));
        }
        sort(sim.begin(), sim.end(), maxScoreCmp);
        for(int i = 0; i < k; i++){
          topKList.push_back(sim[i].first);
      }
      return topKList;
    }

    // resacc
    void resacc(int sourceNode, int k_hops, double r_max_hop, double r_max_f, double epsilon){
        double* residual = new double[vert];
        double* map_ppr = new double[vert];
        int* hops_from_source = new int[vert];
        for(int i = 0; i < vert; i++){
            map_ppr[i] = 0.0;
            residual[i] = 0.0;
            hops_from_source[i] = numeric_limits<int>::max();
        }
        residual[sourceNode] = 1.0;
        hops_from_source[sourceNode] = 0;
        double r_sum = 0.0;
        unordered_set<int> kHopSet;
        unordered_set<int> k1HopLayer;
        //unordered_set<int> danglingNodes;
        //unordered_set<int> hubs;

        //kHopFWD from the source
        clock_t startKHopFWDTime = clock();
        cout<<"Start the h-HopFWD phase ..." <<endl;
        cout<<"----------------------------" <<endl;
        r_max_hop = epsilon * sqrt(r_max_hop) / sqrt((g.m) * 3.0 * log(2.0 * g.n) * g.n);
        //cout<<"k hops: " << k_hops <<endl;
        k1HopLayer = kHopFWD(sourceNode, k_hops, hops_from_source, kHopSet, k1HopLayer, r_max_hop, map_ppr, residual);
        //cout<< "size of k_hop_set: " << kHopSet.size() <<endl;
        //cout<<"size of k1_hop_layer: " << k1HopLayer.size() <<endl;

        clock_t endKHopFWDTime = clock();

        double kHopFWDTime = (endKHopFWDTime - startKHopFWDTime) / (double) CLOCKS_PER_SEC;
        //cout << "kHopFWD Time: " << kHopFWDTime << endl;
        avg_khopfwd_time += kHopFWDTime;

        //FWD from the k1 nodes
        clock_t startOMFWDTime = clock();
        cout<<"----------------------------" <<endl;
        cout<<"Start the OMFWD phase ..." <<endl;
        cout<<"----------------------------" <<endl;
        r_max_f = epsilon * r_max_f / sqrt( (g.m) * 3.0 * log(2.0 * g.n) * g.n);
        //cout << "r_max_f: " << r_max_ratio << endl;
        OMFWD(k1HopLayer, sourceNode, r_max_f, map_ppr, residual);

        clock_t endOMFWDTime = clock();
        double OMFWDTime = (endOMFWDTime - startOMFWDTime)/ (double) CLOCKS_PER_SEC;
        //double totalFWDTime = fwdSourceTime + remedyTime;
        //cout << "OMFWD Time: " << OMFWDTime << endl;
        avg_omfwd_time += OMFWDTime;
        //cout<< "OMFWD Process DONE " << endl;

        for(int i=0; i < vert; i++){
            if(g.getOutSize(i) == 0 && residual[i] != 0.0){
                map_ppr[i] += alpha * residual[i];
                residual[sourceNode] += (1-alpha) * residual[i];
                residual[i] = 0;
            }
        }
        for(int i=0; i < vert; i++){
            if(residual[i] > 0.0){
                r_sum += residual[i];
            }
        }
        //cout<<"residue at source: " << residual[sourceNode] <<endl;
        cout<< "After that, r_sum: " << r_sum <<endl;

        //start the random walk process
        cout<<"----------------------------" <<endl;
        cout<<"Start the Remedy phase ..." <<endl;
        cout<<"----------------------------" <<endl;
        clock_t startWalkTime = clock();
        double omega = (2 + epsilon) * log(2.0 * g.n) * g.n / epsilon / epsilon;
        unsigned long num_walks = r_sum * omega; //the parameters decides
        //cout << "# total walks: " << num_walks << endl;
        unsigned long real_total_num_walks = 0;
        for(int i = 0; i < vert; i++){
            if(residual[i] > 0.0){
                int tempNode = i;
                double residue = residual[tempNode];
                //cout<< "residue: " << i << " " << residue << endl;
                unsigned long num_tempNode_rw = ceil(residue / r_sum * num_walks);
                real_total_num_walks += num_tempNode_rw;
                //cout << "# rw for " << i << ": " << num_tempNode_rw << endl;
                double a_tempNode = residue / r_sum * num_walks /num_tempNode_rw;

                double ppr_incre = a_tempNode * r_sum / num_walks;

                for(unsigned long j=0; j<num_tempNode_rw; j++){
                    int des = tempNode;
                    while(R.drand() > alpha) {
                        int length = g.getOutSize(des);
                        if (length == 0){
                            des = sourceNode;
                        }
                        else {
                            int r = R.generateRandom() % length;
                            des = g.getOutVert(des, r);
                        }
                    }
                    //cout << "Final Destination: " << des<< endl;
                    map_ppr[des] += ppr_incre;
                }
            }
        }

        clock_t endWalkTime = clock();
        double walkTime = (endWalkTime - startWalkTime)/ (double) CLOCKS_PER_SEC;
        //cout <<"walk time: " << walkTime << endl;
        avg_walk_time += walkTime;
        double totalTime = kHopFWDTime + OMFWDTime + walkTime;

        //cout<< "walk time: " << walkTime <<endl;
        avg_time += totalTime;
        cout<<"----------------------------" <<endl;
        cout << "One Query DONE ... " << endl;
        cout<<"----------------------------" <<endl;

        stringstream ss;
        ss << "estimated_rwr/" << target_filename << "/"<< sourceNode <<".txt";
        string outfile = ss.str();
        ofstream est_ppr_file(outfile);
        for(int i =0; i < g.n; i++){
            if(map_ppr[i] > 0.0){
                est_ppr_file << i << " " << map_ppr[i] << endl;
            }
        }
        est_ppr_file.close();

        delete [] residual;
        delete [] map_ppr;
        delete [] hops_from_source;

    }

    //simple_resacc: without OMFWD
    void simple_resacc(int sourceNode, int k_hops, double r_max_hop, double epsilon){

        double* residual = new double[vert];
        double* map_ppr = new double[vert];
        int* hops_from_source = new int[vert];
        for(int i = 0; i < vert; i++){
            map_ppr[i] = 0.0;
            residual[i] = 0.0;
            hops_from_source[i] = numeric_limits<int>::max();
        }
        residual[sourceNode] = 1.0;
        hops_from_source[sourceNode] = 0;
        double r_sum = 0.0;
        unordered_set<int> kHopSet;
        unordered_set<int> k1HopLayer;
        unordered_set<int> danglingNodes;
        unordered_set<int> hubs;

        //kHopFWD from the source
        clock_t startKHopFWDTime = clock();
        r_max_hop = epsilon * sqrt(r_max_hop) / sqrt((g.m) * 3.0 * log(2.0 * g.n) * g.n);
        cout<<"k hops: " << k_hops <<endl;
        k1HopLayer = kHopFWD(sourceNode, k_hops, hops_from_source, kHopSet, k1HopLayer, r_max_hop, map_ppr, residual);
        cout<< "size of k_hop_set: " << kHopSet.size() <<endl;
        cout<<"size of k1_hop_layer: " << k1HopLayer.size() <<endl;


        clock_t endKHopFWDTime = clock();

        double kHopFWDTime = (endKHopFWDTime - startKHopFWDTime) / (double) CLOCKS_PER_SEC;
        cout << "kHopFWD Time: " << kHopFWDTime << endl;
        avg_khopfwd_time += kHopFWDTime;

        for(int i=0; i < vert; i++){
            if(g.getOutSize(i) == 0 && residual[i] != 0.0){
                map_ppr[i] += alpha * residual[i];
                residual[sourceNode] += (1-alpha) * residual[i];
                residual[i] = 0;
            }
        }
        for(int i=0; i < vert; i++){
            if(residual[i] > 0.0){
                r_sum += residual[i];
            }
        }
        cout<<"residue at source: " << residual[sourceNode] <<endl;
        cout<< "After that, r_sum: " << r_sum <<endl;

        //start the random walk process
        cout << "Start Random Walk Process" << endl;
        clock_t startWalkTime = clock();
        double omega = (2 + epsilon) * log(2.0 * g.n) * g.n / epsilon / epsilon;
        unsigned long num_walks = r_sum * omega; //the parameters decides
        cout << "# total walks: " << num_walks << endl;
        unsigned long real_total_num_walks = 0;
        for(int i = 0; i < vert; i++) {
            if (residual[i] > 0.0) {
                int tempNode = i;
                double residue = residual[tempNode];
                //cout<< "residue: " << i << " " << residue << endl;
                unsigned long num_tempNode_rw = ceil(residue / r_sum * num_walks);
                real_total_num_walks += num_tempNode_rw;
                //cout << "# rw for " << i << ": " << num_tempNode_rw << endl;
                double a_tempNode = residue / r_sum * num_walks / num_tempNode_rw;

                double ppr_incre = a_tempNode * r_sum / num_walks;

                for (unsigned long j = 0; j < num_tempNode_rw; j++) {
                    int des = tempNode;
                    while (R.drand() > alpha) {
                        int length = g.getOutSize(des);
                        if (length == 0) {
                            des = sourceNode;
                        } else {
                            int r = R.generateRandom() % length;
                            des = g.getOutVert(des, r);
                        }
                    }
                    //cout << "Final Destination: " << des<< endl;
                    map_ppr[des] += ppr_incre;
                }
            }
        }

        clock_t endWalkTime = clock();
        double walkTime = (endWalkTime - startWalkTime)/ (double) CLOCKS_PER_SEC;
        cout <<"walk time: " << walkTime << endl;
        avg_walk_time += walkTime;
        double totalTime = kHopFWDTime + walkTime;

        //cout<< "walk time: " << walkTime <<endl;
        avg_time += totalTime;
        cout << "One Query DONE ... " << endl;

        delete [] residual;
        delete [] map_ppr;
        delete [] hops_from_source;

    }

    //no_self_loop_resacc: remove the self-loop on the hopFWD
    void no_self_loop_resacc(int sourceNode, int k_hops, double r_max_f, double r_max_hop, double epsilon){
        double* residual = new double[vert];
        double* map_ppr = new double[vert];
        int* hops_from_source = new int[vert];

        for(int i = 0; i < vert; i++){
            map_ppr[i] = 0.0;
            residual[i] = 0.0;
            hops_from_source[i] = numeric_limits<int>::max();
        }
        residual[sourceNode] = 1.0;
        hops_from_source[sourceNode] = 0;
        double r_sum = 0.0;
        unordered_set<int> kHopSet;
        unordered_set<int> k1HopLayer;

        //kHopFWD from the source
        clock_t startKHopFWDTime = clock();
        r_max_hop = epsilon * sqrt(r_max_hop) / sqrt((g.m) * 3.0 * log(2.0 * g.n) * g.n);
        cout<<"k hops: " << k_hops <<endl;
        k1HopLayer = no_self_loop_kHopFWD(sourceNode, k_hops, hops_from_source, kHopSet, k1HopLayer, r_max_hop, map_ppr, residual);
        cout<< "size of k_hop_set: " << kHopSet.size() <<endl;
        cout<<"size of k1_hop_layer: " << k1HopLayer.size() <<endl;


        clock_t endKHopFWDTime = clock();

        double kHopFWDTime = (endKHopFWDTime - startKHopFWDTime) / (double) CLOCKS_PER_SEC;
        cout << "kHopFWD Time: " << kHopFWDTime << endl;
        avg_khopfwd_time += kHopFWDTime;

        //FWD from the k1 nodes
        clock_t startOMFWDTime = clock();
        r_max_f = epsilon * r_max_f / sqrt( (g.m) * 3.0 * log(2.0 * g.n) * g.n);
        OMFWD(k1HopLayer, sourceNode, r_max_f, map_ppr, residual);

        clock_t endOMFWDTime = clock();
        double OMFWDTime = (endOMFWDTime - startOMFWDTime)/ (double) CLOCKS_PER_SEC;
        //double totalFWDTime = fwdSourceTime + remedyTime;
        cout << "OMFWD Time: " << OMFWDTime << endl;
        avg_omfwd_time += OMFWDTime;
        cout<< "OMFWD Process DONE " << endl;

        for(int i=0; i < vert; i++){
            if(residual[i] > 0.0){
                r_sum += residual[i];
            }
        }
        //cout<<"residue at source: " << residual[sourceNode] <<endl;
        //cout<< "After that, r_sum: " << r_sum <<endl;

        //start the random walk process
        //cout << "Start Random Walk Process" << endl;
        clock_t startWalkTime = clock();
        double omega = (2 + epsilon) * log(2.0 * g.n) * g.n / epsilon / epsilon;
        unsigned long num_walks = r_sum * omega; //the parameters decides
        //cout << "# total walks: " << num_walks << endl;
        unsigned long real_total_num_walks = 0;
        for(int i = 0; i < vert; i++){
            if(residual[i] > 0.0){
                int tempNode = i;
                double residue = residual[tempNode];
                //cout<< "residue: " << i << " " << residue << endl;
                unsigned long num_tempNode_rw = ceil(residue / r_sum * num_walks);
                real_total_num_walks += num_tempNode_rw;
                //cout << "# rw for " << i << ": " << num_tempNode_rw << endl;
                double a_tempNode = residue / r_sum * num_walks /num_tempNode_rw;

                double ppr_incre = a_tempNode * r_sum / num_walks;

                for(unsigned long j=0; j<num_tempNode_rw; j++){
                    int des = tempNode;
                    while(R.drand() > alpha) {
                        int length = g.getOutSize(des);
                        if (length == 0){
                            des = sourceNode;
                        }
                        else {
                            int r = R.generateRandom() % length;
                            des = g.getOutVert(des, r);
                        }
                    }
                    //cout << "Final Destination: " << des<< endl;
                    map_ppr[des] += ppr_incre;
                }
            }
        }


        clock_t endWalkTime = clock();
        double walkTime = (endWalkTime - startWalkTime)/ (double) CLOCKS_PER_SEC;
        cout <<"walk time: " << walkTime << endl;
        avg_walk_time += walkTime;
        double totalTime = kHopFWDTime + OMFWDTime + walkTime;

        //cout<< "walk time: " << walkTime <<endl;
        avg_time += totalTime;
        cout << "One Query DONE ... " << endl;

        delete [] residual;
        delete [] map_ppr;
        delete [] hops_from_source;

    }

    //a new solution: without the subgraph approach
    void no_subgraph_resacc(int sourceNode, int k_hops, double r_max_f, double r_max_hop, double epsilon){
        double* residual = new double[vert];
        double* map_ppr = new double[vert];
        int* hops_from_source = new int[vert];
        for(int i = 0; i < vert; i++){
            map_ppr[i] = 0.0;
            residual[i] = 0.0;
            hops_from_source[i] = numeric_limits<int>::max();
        }
        residual[sourceNode] = 1.0;
        hops_from_source[sourceNode] = 0;
        double r_sum = 0.0;
        unordered_set<int> kHopSet;
        unordered_set<int> k1HopLayer;
        //unordered_set<int> danglingNodes;
        //unordered_set<int> hubs;


        //kHopFWD from the source
        clock_t startKHopFWDTime = clock();
        r_max_hop = epsilon * sqrt(r_max_hop) / sqrt((g.m) * 3.0 * log(2.0 * g.n) * g.n);
        cout<<"k hops: " << k_hops <<endl;
        k1HopLayer = no_subgraph_kHopFWD(sourceNode, k_hops, hops_from_source, kHopSet, k1HopLayer, r_max_hop, map_ppr, residual);
        cout<< "size of k_hop_set: " << kHopSet.size() <<endl;
        cout<<"size of k1_hop_layer: " << k1HopLayer.size() <<endl;


        clock_t endKHopFWDTime = clock();

        double kHopFWDTime = (endKHopFWDTime - startKHopFWDTime) / (double) CLOCKS_PER_SEC;
        cout << "kHopFWD Time: " << kHopFWDTime << endl;
        avg_khopfwd_time += kHopFWDTime;

        //FWD from the k1 nodes
        clock_t startOMFWDTime = clock();
        r_max_f = epsilon * r_max_f / sqrt( (g.m) * 3.0 * log(2.0 * g.n) * g.n);
        cout << "r_max_f: " << r_max_f << endl;
        OMFWD(k1HopLayer, sourceNode, r_max_f, map_ppr, residual);

        clock_t endOMFWDTime = clock();
        double OMFWDTime = (endOMFWDTime - startOMFWDTime)/ (double) CLOCKS_PER_SEC;
        //double totalFWDTime = fwdSourceTime + remedyTime;
        cout << "OMFWD Time: " << OMFWDTime << endl;
        avg_omfwd_time += OMFWDTime;
        cout<< "OMFWD Process DONE " << endl;

        for(int i=0; i < vert; i++){
            if(g.getOutSize(i) == 0 && residual[i] != 0.0){
                map_ppr[i] += alpha * residual[i];
                residual[sourceNode] += (1-alpha) * residual[i];
                residual[i] = 0;
            }
        }
        for(int i=0; i < vert; i++){
            if(residual[i] > 0.0){
                r_sum += residual[i];
            }
        }
        cout<<"residue at source: " << residual[sourceNode] <<endl;
        cout<< "After that, r_sum: " << r_sum <<endl;

        //start the random walk process
        cout << "Start Random Walk Process" << endl;
        clock_t startWalkTime = clock();
        double omega = (2 + epsilon) * log(2.0 * g.n) * g.n / epsilon / epsilon;
        unsigned long num_walks = r_sum * omega; //the parameters decides
        cout << "# total walks: " << num_walks << endl;
        unsigned long real_total_num_walks = 0;
        for(int i = 0; i < vert; i++){
            if(residual[i] > 0.0){
                int tempNode = i;
                double residue = residual[tempNode];
                //cout<< "residue: " << i << " " << residue << endl;
                unsigned long num_tempNode_rw = ceil(residue / r_sum * num_walks);
                real_total_num_walks += num_tempNode_rw;
                //cout << "# rw for " << i << ": " << num_tempNode_rw << endl;
                double a_tempNode = residue / r_sum * num_walks /num_tempNode_rw;

                double ppr_incre = a_tempNode * r_sum / num_walks;

                for(unsigned long j=0; j<num_tempNode_rw; j++){
                    int des = tempNode;
                    while(R.drand() > alpha) {
                        int length = g.getOutSize(des);
                        if (length == 0){
                            des = sourceNode;
                        }
                        else {
                            int r = R.generateRandom() % length;
                            des = g.getOutVert(des, r);
                        }
                    }
                    //cout << "Final Destination: " << des<< endl;
                    map_ppr[des] += ppr_incre;
                }
            }
        }
        cout << "real total number of walks: " << real_total_num_walks << endl;
        //updating process by using the remaining residue at source



        clock_t endWalkTime = clock();
        double walkTime = (endWalkTime - startWalkTime)/ (double) CLOCKS_PER_SEC;
        cout <<"walk time: " << walkTime << endl;
        avg_walk_time += walkTime;
        double totalTime = kHopFWDTime + OMFWDTime + walkTime;

        //cout<< "walk time: " << walkTime <<endl;
        avg_time += totalTime;
        cout << "One Query DONE ... " << endl;

        delete [] residual;
        delete [] map_ppr;
        delete [] hops_from_source;

    }


    //kHopFWD with CL: automatically find the nodes in the khopset and k1layer
    unordered_set<int> kHopFWD(int sourceNode, int k_hops, int* hops_from_source, unordered_set<int> kHopSet, unordered_set<int> k1HopLayer, double r_max_hop, double* reserve, double* residual){
        kHopSet.insert(sourceNode);

        queue<int> r_queue;
        //queue<int> r_hub_queue;

        static vector<bool> idx(g.n);
        //static vector<bool> hub_idx(g.n);
        std::fill(idx.begin(), idx.end(), false);
        //std::fill(hub_idx.begin(), hub_idx.end(),false);
        //the very first push at source
        reserve[sourceNode] += alpha * residual[sourceNode];
        int out_deg = g.getOutSize(sourceNode);

        double remain_residual = (1 - alpha) * residual[sourceNode];
        residual[sourceNode] = 0.0;

        if (out_deg == 0) {
            residual[sourceNode] += remain_residual;
        } else {
            double avg_push_residual = remain_residual / out_deg;
            for (int i = 0; i < g.getOutSize(sourceNode); i++) {
                int next = g.getOutVert(sourceNode, i);
                hops_from_source[next] = 1;
                kHopSet.insert(next);
                residual[next] += avg_push_residual;
                if(k_hops > 0 && next != sourceNode && (g.getOutSize(next) != 0) && (residual[next]/g.getOutSize(next) >= r_max_hop) && idx[next] != true) {
                    r_queue.push(next);
                    idx[next] = true;
                }else if(k_hops == 0){
                    k1HopLayer.insert(next);
                }

            }
        }

        // run the forward push for other nodes with non-zero residual

        while(r_queue.size() > 0) {
                int tempNode = r_queue.front();
                r_queue.pop();
                idx[tempNode] =false;
                // cout<<"current node: " << tempNode << endl;
                // cout<<"its residual: " << residual[tempNode] << endl;
                double current_residual = residual[tempNode];
                if(current_residual > 0.0){

                    residual[tempNode] = 0.0;
                    reserve[tempNode] += alpha * current_residual; //
                    double remain_residual = (1 - alpha) * current_residual;

                    int out_deg = g.getOutSize(tempNode);
                    if (out_deg == 0) {
                        residual[sourceNode] += remain_residual;
                        continue;
                    }

                    double avg_push_residual = remain_residual / out_deg;
                    int hops_of_tempNode = hops_from_source[tempNode];

                    for (int i = 0; i < g.getOutSize(tempNode); i++) {
                        int next = g.getOutVert(tempNode, i);
                        residual[next] += avg_push_residual;
                        //the can-pushed nodes to the queue
                        if(next != sourceNode){
                            //check the hops between next and source
                            if(hops_from_source[next] > hops_of_tempNode + 1){
                                hops_from_source[next] = hops_of_tempNode +1;
                            }
                            if(hops_from_source[next] <= k_hops){
                                // next is in the k-hop set of source, which is able to push its residue (in queue)
                                if((g.getOutSize(next) != 0) && residual[next]/ g.getOutSize(next) >= r_max_hop){
                                    if(idx[next] != true){
                                        r_queue.push(next);
                                        idx[next] =true;
                                    }
                                }

                                if(kHopSet.count(next) == 0){
                                    kHopSet.insert(next);
                                    if(k1HopLayer.count(next) > 0){
                                        k1HopLayer.erase(next);
                                    }
                                }
                            }
                            else if(hops_from_source[next] == k_hops+1){
                                if(k1HopLayer.count(next) == 0){
                                    k1HopLayer.insert(next);
                                }
                            }

                        }
                        else{
                            continue;
                        }
                    }
                }

            }
        for(int i=0; i < vert; i++){
            if(g.getOutSize(i) == 0 && residual[i] != 0.0){
                reserve[i] += alpha * residual[i];
                residual[sourceNode] += (1-alpha) * residual[i];
                residual[i] = 0;
            }
        }

        if(k_hops > 0 && residual[sourceNode] != 0) {
            double min_r_max = 0.000000000000001;
            unsigned long num_loop_from_source = (int) ceil(log(min_r_max * g.getOutSize(sourceNode)) / log(residual[sourceNode]));

            cout << "# of iteration from source: " << num_loop_from_source << endl;

            double upper_scaler_reserve =
                    (1 - pow(residual[sourceNode], num_loop_from_source - 1)) / (1 - residual[sourceNode]);
            cout << "Scaler: " << upper_scaler_reserve << endl;

            for (auto it_kHopSet = kHopSet.begin(); it_kHopSet != kHopSet.end(); ++it_kHopSet) {
                int tempNode = *it_kHopSet;
                reserve[tempNode] = reserve[tempNode] * upper_scaler_reserve;
                //check += reserve[tempNode];
                if (tempNode != sourceNode) {
                    residual[tempNode] = residual[tempNode] * upper_scaler_reserve;
                } else {
                    residual[tempNode] = pow(residual[tempNode], num_loop_from_source);
                }

            }
            for (auto it_k1Layer = k1HopLayer.begin(); it_k1Layer != k1HopLayer.end(); ++it_k1Layer) {
                int tempNode = *it_k1Layer;
                residual[tempNode] = residual[tempNode] * upper_scaler_reserve;
            }
        }

        return k1HopLayer;
    }

    //khopFWD without self-loop
    unordered_set<int> no_self_loop_kHopFWD(int sourceNode, int k_hops, int* hops_from_source, unordered_set<int> kHopSet, unordered_set<int> k1HopLayer, double r_max, double* reserve, double* residual){
        kHopSet.insert(sourceNode);

        queue<int> r_queue;

        static vector<bool> idx(g.n);
        std::fill(idx.begin(), idx.end(), false);
        //the very first push at source
        reserve[sourceNode] += alpha * residual[sourceNode];
        int out_deg = g.getOutSize(sourceNode);

        double remain_residual = (1 - alpha) * residual[sourceNode];
        residual[sourceNode] = 0.0;

        if (out_deg == 0) {
            residual[sourceNode] += remain_residual;
        } else {
            double avg_push_residual = remain_residual / out_deg;
            for (int i = 0; i < g.getOutSize(sourceNode); i++) {
                int next = g.getOutVert(sourceNode, i);
                hops_from_source[next] = 1;
                kHopSet.insert(next);
                residual[next] += avg_push_residual;
                if(k_hops > 0 && (g.getOutSize(next) != 0) && (residual[next]/g.getOutSize(next) >= r_max) && idx[next] != true) {
                    r_queue.push(next);
                    idx[next] = true;
                }else if(k_hops == 0){
                    k1HopLayer.insert(next);
                }

            }
        }

        // run the forward push for other nodes with non-zero residual

        while(r_queue.size() > 0) {
            int tempNode = r_queue.front();
            r_queue.pop();
            idx[tempNode] =false;
            double current_residual = residual[tempNode];
            if(current_residual > 0.0){

                residual[tempNode] = 0.0;
                reserve[tempNode] += alpha * current_residual; //
                double remain_residual = (1 - alpha) * current_residual;

                int out_deg = g.getOutSize(tempNode);
                if (out_deg == 0) {
                    residual[sourceNode] += remain_residual;
                    continue;
                }

                double avg_push_residual = remain_residual / out_deg;
                int hops_of_tempNode = hops_from_source[tempNode];

                for (int i = 0; i < g.getOutSize(tempNode); i++) {
                    int next = g.getOutVert(tempNode, i);
                    residual[next] += avg_push_residual;
                    //the can-pushed nodes to the queue
                        //check the hops between next and source
                        if(hops_from_source[next] > hops_of_tempNode + 1){
                            hops_from_source[next] = hops_of_tempNode +1;
                        }
                        if(hops_from_source[next] <= k_hops){
                            // next is in the k-hop set of source, which is able to push its residue (in queue)
                            if((g.getOutSize(next) != 0) && residual[next]/ g.getOutSize(next) >= r_max){
                                if(idx[next] != true){
                                    r_queue.push(next);
                                    idx[next] =true;
                                }
                            }

                            if(kHopSet.count(next) == 0){
                                kHopSet.insert(next);
                                if(k1HopLayer.count(next) > 0){
                                    k1HopLayer.erase(next);
                                }
                            }
                        }
                        else if(hops_from_source[next] == k_hops+1){
                            if(k1HopLayer.count(next) == 0){
                                k1HopLayer.insert(next);
                            }
                        }
                }
            }

        }

        return k1HopLayer;
    }

    //khopFWD without subgraph: i.e., in the whole graph
    unordered_set<int> no_subgraph_kHopFWD(int sourceNode, int k_hops, int* hops_from_source, unordered_set<int> kHopSet, unordered_set<int> k1HopLayer, double r_max, double* reserve, double* residual){
        kHopSet.insert(sourceNode);

        queue<int> r_queue;

        static vector<bool> idx(g.n);
        std::fill(idx.begin(), idx.end(), false);
        //std::fill(hub_idx.begin(), hub_idx.end(),false);
        //the very first push at source
        reserve[sourceNode] += alpha * residual[sourceNode];
        int out_deg = g.getOutSize(sourceNode);

        double remain_residual = (1 - alpha) * residual[sourceNode];
        residual[sourceNode] = 0.0;

        if (out_deg == 0) {
            residual[sourceNode] += remain_residual;
        } else {
            double avg_push_residual = remain_residual / out_deg;
            for (int i = 0; i < g.getOutSize(sourceNode); i++) {
                int next = g.getOutVert(sourceNode, i);
                hops_from_source[next] = 1;
                kHopSet.insert(next);
                residual[next] += avg_push_residual;
                if( next != sourceNode && (g.getOutSize(next) != 0) && (residual[next]/g.getOutSize(next) >= r_max) && idx[next] != true) {
                    r_queue.push(next);
                    idx[next] = true;
                }

            }
        }

        // run the forward push for other nodes with non-zero residual

        while(r_queue.size() > 0) {
            int tempNode = r_queue.front();
            r_queue.pop();
            idx[tempNode] =false;
            double current_residual = residual[tempNode];
            if(current_residual > 0.0){

                residual[tempNode] = 0.0;
                reserve[tempNode] += alpha * current_residual; //
                double remain_residual = (1 - alpha) * current_residual;

                int out_deg = g.getOutSize(tempNode);
                if (out_deg == 0) {
                    residual[sourceNode] += remain_residual;
                    continue;
                }

                double avg_push_residual = remain_residual / out_deg;
                int hops_of_tempNode = hops_from_source[tempNode];

                for (int i = 0; i < g.getOutSize(tempNode); i++) {
                    int next = g.getOutVert(tempNode, i);
                    residual[next] += avg_push_residual;
                    //the can-pushed nodes to the queue
                    if(next != sourceNode){
                        //check the hops between next and source
                        if(hops_from_source[next] > hops_of_tempNode + 1){
                            hops_from_source[next] = hops_of_tempNode +1;
                        }

                            if((g.getOutSize(next) != 0) && residual[next]/ g.getOutSize(next) >= r_max){
                                if(idx[next] != true){
                                    r_queue.push(next);
                                    idx[next] =true;
                                }
                            }

                    }
                    else{
                        continue;
                    }
                }
            }

        }
        if(k_hops > 0 && residual[sourceNode] != 0) {
            double min_r_max = 0.000000000000001; //(10^{-15})
            unsigned long num_loop_from_source = (int) ceil(log(min_r_max * g.getOutSize(sourceNode)) / log(residual[sourceNode]));

            cout << "# of iteration from source: " << num_loop_from_source << endl;

            double upper_scaler_reserve =
                    (1 - pow(residual[sourceNode], num_loop_from_source - 1)) / (1 - residual[sourceNode]);
            cout << "Scaler: " << upper_scaler_reserve << endl;


            for (int i =0; i < g.n; i++) {
                int tempNode = i;
                reserve[tempNode] = reserve[tempNode] * upper_scaler_reserve;
                //check += reserve[tempNode];
                if (tempNode != sourceNode) {
                    residual[tempNode] = residual[tempNode] * upper_scaler_reserve;
                } else {
                    residual[tempNode] = pow(residual[tempNode], num_loop_from_source);
                }

            }
        }

        for(int i=0; i<g.n; i++){
            if(residual[i]> 0.0){
                k1HopLayer.insert(i);
            }
        }

        return k1HopLayer;
    }

    //FWD by pushing all the k1 nodes into the queue
    void OMFWD(unordered_set<int> k1HopLayer, int sourceNode, double r_max, double* reserve, double* residual){
        static vector<bool> idx(g.n);
        std::fill(idx.begin(), idx.end(), false);
        queue<int> r_queue;
        //queue<int> r_hub_queue;
        //not matter how large the residue of source node would be
        //put it into the queue so that its residue could be reduced with one more forward push operations
        for(auto it_k1layer = k1HopLayer.begin(); it_k1layer != k1HopLayer.end(); ++it_k1layer){
            int k1_node = *it_k1layer;
            r_queue.push(k1_node);
            idx[k1_node] = true;
        }

        // run the forward push for all nodes with non-zero residual

        while(r_queue.size() > 0 ) {
                int tempNode = r_queue.front();
                r_queue.pop();
                idx[tempNode] = false;
                // cout<<"current node: " << tempNode << endl;
                // cout<<"its residual: " << residual[tempNode] << endl;
                double current_residual = residual[tempNode];
                if(current_residual > 0.0){

                    residual[tempNode] = 0.0;
                    reserve[tempNode] += alpha * current_residual; //
                    double remain_residual = (1 - alpha) * current_residual;

                    int out_deg = g.getOutSize(tempNode);
                    if (out_deg == 0) {
                        residual[sourceNode] += remain_residual;
                        if (residual[sourceNode] / g.getOutSize(sourceNode) >= r_max && idx[sourceNode]!= true) {
                            idx[sourceNode] = true;
                            r_queue.push(sourceNode);

                        }
                        continue;
                    }

                    double avg_push_residual = remain_residual / out_deg;

                    for (int i = 0; i < g.getOutSize(tempNode); i++) {
                        int next = g.getOutVert(tempNode, i);
                        residual[next] += avg_push_residual;

                        if((residual[next] / g.getOutSize(next) >= r_max ) && idx[next] != true){
                            idx[next] = true;
                            r_queue.push(next);
                        }
                    }
                }

            }
        }

    // the original fwd + source node must be pushed no matter whether its value satisfies the push condition
    void FWD(int sourceNode, double r_max, double* reserve, double* residual){
        static vector<bool> idx(g.n);
        std::fill(idx.begin(), idx.end(), false);
        queue<int> r_queue;
        //not matter how large the residue of source node would be
        //put it into the queue so that its residue could be reduced with one more forward push operations
        r_queue.push(sourceNode);
        idx[sourceNode] = true;

        // run the forward push for all nodes with non-zero residual

        while(r_queue.size() > 0) {
            int tempNode = r_queue.front();
            r_queue.pop();
            idx[tempNode] = false;
            // cout<<"current node: " << tempNode << endl;
            // cout<<"its residual: " << residual[tempNode] << endl;
            double current_residual = residual[tempNode];
            if(current_residual > 0.0){

                residual[tempNode] = 0.0;
                reserve[tempNode] += alpha * current_residual; //
                double remain_residual = (1 - alpha) * current_residual;

                int out_deg = g.getOutSize(tempNode);
                if (out_deg == 0) {
                    residual[sourceNode] += remain_residual;
                    if (residual[sourceNode] / g.getOutSize(sourceNode) >= r_max && idx[sourceNode]!= true) {
                        idx[sourceNode] = true;
                        r_queue.push(sourceNode);

                    }
                    continue;
                }

                double avg_push_residual = remain_residual / out_deg;

                for (int i = 0; i < g.getOutSize(tempNode); i++) {
                    int next = g.getOutVert(tempNode, i);
                    residual[next] += avg_push_residual;
                    // next is in the k-hop set of source, which is able to push its residue (in queue)
                    if (residual[next] / g.getOutSize(next) >= r_max && idx[next] != true) {
                        idx[next] = true;
                        r_queue.push(next);

                    }

                }
                }
            }

        }


    // powerMethod From the source node (forward method)
    void t_PowerMethod(vector<int> nodeList, int iterations){
        for(int i = 0; i < nodeList.size(); i++){
            int tempNode = nodeList[i];
            stringstream ss;
            ss << "real_rwr/" << target_filename << "/" << tempNode << ".txt";
            string outputFile = ss.str();
            cout << "file: " << outputFile << endl;
            PowerMethodK(iterations, outputFile, tempNode, 500);
        cout << outputFile << "done!"  << endl;
        }
    }

    void PowerMethodK(int iterations, string outputFile, int u, int k){
        unordered_map<int, double> map_residual;
        map_residual.clear();
        map_residual[u] = 1.0;

        int num_iter=0;
        double* map_ppr = new double[vert];
        for(int i = 0; i < vert; i++){
            map_ppr[i] = 0;
        }
        while( num_iter < iterations ){
            cout << u << ": iter " << num_iter << endl;
            num_iter++;

            vector< pair<int,double> > pairs(map_residual.begin(), map_residual.end());
            map_residual.clear();
            for(auto &p: pairs){
                if(p.second > 0){
                    map_ppr[p.first] += alpha*p.second;
                    int out_deg = g.getOutSize(p.first);

                    double remain_residual = (1-alpha)*p.second;
                    if(out_deg==0){
                        map_residual[u] += remain_residual;
                    }
                    else{
                        double avg_push_residual = remain_residual / out_deg;
                        for(int i = 0; i < g.getOutSize(p.first); i++){
                            int next = g.getOutVert(p.first, i);
                            map_residual[next] += avg_push_residual;
                        }
                    }
                }
            }
        }
        ofstream fout(outputFile);
        vector<pair<int, double> > pprs;
        for(int j = 0; j < vert; j++){
            pprs.push_back(pair<int, double>(j, map_ppr[j]));
        }
        sort(pprs.begin(), pprs.end(), maxScoreCmp);
        for(int j = 0; j < vert; j++){
            if(pprs[j].second >= 0){
                fout << pprs[j].first << " " << pprs[j].second << "\n";
            }

        }
        fout.close();
        delete[] map_ppr;
    }

    //calculate Precision
    double calPrecision(vector<int> topK1, vector<int> realList, int k, bool isShowMissing = false){
        int size = realList.size();
        int size2 = topK1.size();
        int hitCount = 0;
        for(int i = 0; i < size2; i++) {
            bool isFind = false;
            for (int j = 0; j < size; j++) {
                if (topK1[i] == realList[j]) {
                    hitCount++;
                    isFind = true;
                    break;
                }
            }
            /*if(!isFind){
               cout << "useless node: " << topK1[i] << endl;
            }*/
        }
        cout << "hit Count: " << hitCount << endl;
        double result = hitCount / (double) k;
        return result < 1 ? result : 1;

    }

    double calNDCG(vector<int> candidates, int k, int s){
        vector<int> topK = getRealTopK(s, k);
        unordered_map<int, double> realMap = getRealTopKMap(s, k);

        double correct = 0;
        for(int i = 0; i < k; i++){
            if(realMap[candidates[i]] == realMap[topK[i]])
                correct++;
            else{
                cout << "misMatch : " << candidates[i] << ", " << topK[i] << endl;
            }
        }
        return correct / (double)k;

        double Zp = 0;
        for(int i = 1; i <= k; i++){
            Zp += (pow(2, realMap[topK[i-1]]) - 1) / (log(i+1) / log(2));
        }
        double NDCG = 0;
        for(int i = 1; i <= k; i++){
            NDCG += (pow(2, realMap[candidates[i-1]]) - 1) / (log(i+1) / log(2));
        }
        return NDCG / Zp;
    }

    //generate random query node
    void generateQueryNode(int nodeNum, ofstream& fout){
        for(int i = 0; i < nodeNum; i++){
            int tempNode = R.generateRandom() % vert;
            if(g.getOutSize(tempNode) == 0){
                i--;
                continue;
            }
            fout << tempNode << endl;
        }
    }


    //the later algorithms for community detection
    //find the seeds by spread hubs
    void find_seeds_by_spread_hubs(int num_seeds, vector<int>& set_seeds){
        cout<<"# of seeds required: " << num_seeds <<endl;
        vector<bool> marked;
       //queue<int> unmarked;
        vector<pair<int, double>> sorted_deg;
        for(int i=0; i<vert; i++){
            marked.push_back(false);
            //unmarked.push(i);
            sorted_deg.push_back(make_pair(i, g.getOutSize(i)));
        }

        //sort all the degrees
        cout << "sorting degrees..." << endl;
        sort(sorted_deg.begin(), sorted_deg.end(),
             [](pair<int, double> const &l, pair<int, double> const &r) { return l.second > r.second; });

        int flag_max_deg = 0;
        while(set_seeds.size() < num_seeds){
            int node_max_deg;
            while(marked[sorted_deg[flag_max_deg].first] != false && flag_max_deg < vert){
                flag_max_deg++;
            }
            if(flag_max_deg < vert){
                node_max_deg = sorted_deg[flag_max_deg].first;
                marked[node_max_deg] = true;
                set_seeds.push_back(node_max_deg);
                for(int i=0; i< g.getOutSize(node_max_deg); i++){
                    int out_neighbour = g.outAdjList[node_max_deg][i];
                    marked[out_neighbour] = true;
                }
            }
            else{
                cout<<"cannot find enough seeds "<<endl;
                break;
            }

        }

        cout<<"finding done.." <<endl;
        cout<<"# of seeds: " << set_seeds.size()<<endl;
        //for(int i=0; i<set_seeds.size(); i++){
        //    cout<<set_seeds[i]<<endl;
        //}
    }

    //seed_expansion_by_ResAcc
    void seed_expansion_by_resacc(int sourceNode, int k_hops, double r_max_hop, double r_max_f, double epsilon, vector<pair<int, double>>& pprs){
        double* residual = new double[vert];
        double* map_ppr = new double[vert];
        //double* temp_ppr_v = new double[vert];
        //double* temp_residual_v = new double[vert];
        int* hops_from_source = new int[vert];
        //vector<pair<int, int>> node_outdeg(vert);
        for(int i = 0; i < vert; i++){
            map_ppr[i] = 0.0;
            residual[i] = 0.0;
            //temp_ppr_v[i] = 0.0;
            //temp_residual_v[i] = 0.0;
            hops_from_source[i] = numeric_limits<int>::max();
            //node_outdeg[i] = make_pair(i, g.getOutSize(i));
        }
        residual[sourceNode] = 1.0;
        hops_from_source[sourceNode] = 0;
        double r_sum = 0.0;
        unordered_set<int> kHopSet;
        unordered_set<int> k1HopLayer;

        //kHopFWD from the source
      //  clock_t startKHopFWDTime = clock();
        r_max_hop = epsilon * sqrt(r_max_hop) / sqrt((g.m) * 3.0 * log(2.0 * g.n) * g.n);
        cout<<"k hops: " << k_hops <<endl;
        k1HopLayer = kHopFWD(sourceNode, k_hops, hops_from_source, kHopSet, k1HopLayer, r_max_hop, map_ppr, residual);
        cout<< "size of k_hop_set: " << kHopSet.size() <<endl;
        cout<<"size of k1_hop_layer: " << k1HopLayer.size() <<endl;


       // clock_t endKHopFWDTime = clock();

       // double kHopFWDTime = (endKHopFWDTime - startKHopFWDTime) / (double) CLOCKS_PER_SEC;
       // cout << "kHopFWD Time: " << kHopFWDTime << endl;
       // avg_khopfwd_time += kHopFWDTime;

        //FWD from the k1 nodes
       /// clock_t startOMFWDTime = clock();
        r_max_f = epsilon * r_max_f / sqrt( (g.m) * 3.0 * log(2.0 * g.n) * g.n);
        cout << "r_max_f: " << r_max_f << endl;
        OMFWD(k1HopLayer, sourceNode, r_max_f, map_ppr, residual);

       // clock_t endOMFWDTime = clock();
       // double OMFWDTime = (endOMFWDTime - startOMFWDTime)/ (double) CLOCKS_PER_SEC;
        //double totalFWDTime = fwdSourceTime + remedyTime;
        //cout << "OMFWD Time: " << OMFWDTime << endl;
       // avg_omfwd_time += OMFWDTime;
        cout<< "OMFWD Process DONE " << endl;

        for(int i=0; i < vert; i++){
            if(g.getOutSize(i) == 0 && residual[i] != 0.0){
                map_ppr[i] += alpha * residual[i];
                residual[sourceNode] += (1-alpha) * residual[i];
                residual[i] = 0;
            }
        }
        for(int i=0; i < vert; i++){
            if(residual[i] > 0.0){
                r_sum += residual[i];
            }
        }
        cout<<"residue at source: " << residual[sourceNode] <<endl;
        cout<< "After that, r_sum: " << r_sum <<endl;

        //start the random walk process
        cout << "Start Random Walk Process" << endl;

        //clock_t startWalkTime = clock();
        double omega = (2 + epsilon) * log(2.0 * g.n) * g.n / epsilon / epsilon;
        unsigned long num_walks = r_sum * omega; //the parameters decides
        cout << "# total walks: " << num_walks << endl;
        unsigned long real_total_num_walks = 0;
        for(int i = 0; i < vert; i++){
            if(residual[i] > 0.0){
                int tempNode = i;
                double residue = residual[tempNode];
                //cout<< "residue: " << i << " " << residue << endl;
                unsigned long num_tempNode_rw = ceil(residue / r_sum * num_walks);
                real_total_num_walks += num_tempNode_rw;
                //cout << "# rw for " << i << ": " << num_tempNode_rw << endl;
                double a_tempNode = residue / r_sum * num_walks /num_tempNode_rw;

                double ppr_incre = a_tempNode * r_sum / num_walks;

                for(unsigned long j=0; j<num_tempNode_rw; j++){
                    int des = tempNode;
                    while(R.drand() > alpha) {
                        int length = g.getOutSize(des);
                        if (length == 0){
                            des = sourceNode;
                        }
                        else {
                            int r = R.generateRandom() % length;
                            des = g.getOutVert(des, r);
                        }
                    }
                    //cout << "Final Destination: " << des<< endl;
                    map_ppr[des] += ppr_incre;
                }
            }
        }
        cout << "real total number of walks: " << real_total_num_walks << endl;
        //updating process by using the remaining residue at source



        //clock_t endWalkTime = clock();
       // double walkTime = (endWalkTime - startWalkTime)/ (double) CLOCKS_PER_SEC;
       // cout <<"walk time: " << walkTime << endl;
        //avg_walk_time += walkTime;
      //  double totalTime = kHopFWDTime + OMFWDTime + walkTime;

        //cout<< "walk time: " << walkTime <<endl;
        //avg_time += totalTime;
        cout << "One Query DONE ... " << endl;

        for(int i=0; i<vert; i++){
            double ppr_divide_deg = map_ppr[i] / (double) g.getOutSize(i);
            pprs.push_back(make_pair(i, ppr_divide_deg));
        }

        delete [] residual;
        delete [] map_ppr;
        delete [] hops_from_source;
    }

    //seed_expansion_by_power
    void seed_expansion_by_power(int sourceNode, int iterations, vector<pair<int,double>>& pprs){
        unordered_map<int, double> map_residual;
        map_residual.clear();
        map_residual[sourceNode] = 1.0;

        int num_iter=0;
        double* map_ppr = new double[vert];
        for(int i = 0; i < vert; i++){
            map_ppr[i] = 0;
        }
        while( num_iter < iterations ){
            cout << sourceNode << ": iter " << num_iter << endl;
            num_iter++;

            vector< pair<int,double> > pairs(map_residual.begin(), map_residual.end());
            map_residual.clear();
            for(auto &p: pairs){
                if(p.second > 0){
                    map_ppr[p.first] += alpha*p.second;
                    int out_deg = g.getOutSize(p.first);

                    double remain_residual = (1-alpha)*p.second;
                    if(out_deg==0){
                        map_residual[sourceNode] += remain_residual;
                    }
                    else{
                        double avg_push_residual = remain_residual / out_deg;
                        for(int i = 0; i < g.getOutSize(p.first); i++){
                            int next = g.getOutVert(p.first, i);
                            map_residual[next] += avg_push_residual;
                        }
                    }
                }
            }
        }
        cout << "One Query DONE ... " << endl;

        for(int i=0; i<vert; i++){
            double ppr_divide_deg = map_ppr[i] / (double) g.getOutSize(i);
            pprs.push_back(make_pair(i, ppr_divide_deg));
        }
        delete[] map_ppr;
    }

    void find_communities_by_RWR(string method, vector<int> set_seeds, int k_hop, double r_max_ratio, double r_max, double epsilon){
        double avg_time_expansion = 0;
        double avg_time_community = 0;
        double avg_normalized_cut = 0;
        double avg_conductance = 0;
        for(int i=0; i<set_seeds.size(); i++){
            clock_t startExpansionTime = clock();
            int temp_seed = set_seeds[i];
            //start the seed expansion by ppr for each seed
            //vector<double> pprs_seed;
            vector<pair<int, double>> pprs_seed;
            if(method == "resacc"){
                seed_expansion_by_resacc(temp_seed, k_hop, r_max_ratio, r_max, epsilon, pprs_seed);
            }
            else if(method == "power"){
                seed_expansion_by_power(temp_seed, 1000, pprs_seed);
            }


            clock_t endExpansionTime = clock();
            avg_time_expansion += (endExpansionTime - startExpansionTime)/ (double) CLOCKS_PER_SEC;

            cout<<"starting computing conductance..." <<endl;
            clock_t startCommunityTime = clock();
            //sort their ppr values (divided by deg(x))
            sort(pprs_seed.begin(), pprs_seed.end(),
                 [](pair<int, double> const &l, pair<int, double> const &r) { return l.second > r.second; });

            //compute the communities
            //vector<int> commun_set;
            unordered_set<int> commun_set;
            long volum_commun_set = 0;
            long cut_commun_set = 0;
            double best_conductance = 1.0;
            unordered_set<int> best_communset;
            best_communset.insert(pprs_seed[0].first); //the first node

            for(int i = 0; i<vert * 0.05; i++){
                int temp_node = pprs_seed[i].first;
                volum_commun_set += g.getOutSize(temp_node);
                for(int j=0; j<g.getOutSize(temp_node); j++){
                    //check where an out-neighbour of this node is in the set
                    int temp_out_neighbour = g.outAdjList[temp_node][j];
                    if(commun_set.count(temp_out_neighbour) == 0){
                        cut_commun_set +=1;
                    }
                    else{
                        cut_commun_set -= 1;
                    }
                }
                //cout<<"v: " <<temp_node <<" cut: " << cut_commun_set <<" vol: " << volum_commun_set <<endl;
                commun_set.insert(temp_node);
                double current_conductance = cut_commun_set / (double) min(volum_commun_set, g.m- volum_commun_set);

                if(current_conductance < best_conductance){
                    best_conductance = current_conductance;
                    best_communset = commun_set;
                }
            }

            clock_t endCommunityTime = clock();
            avg_time_community += (endCommunityTime - startCommunityTime)/ (double) CLOCKS_PER_SEC;
            cout<<"best conductance: " << best_conductance <<endl;
            avg_conductance += best_conductance;
            avg_normalized_cut += cut_commun_set / (double) volum_commun_set;
        }

        stringstream ss;
        if(method =="resacc"){
            ss << "community_detection_results/" << target_filename << "/resacc.txt";
        }
        else if(method =="power"){
            ss << "community_detection_results/" << target_filename << "/power.txt";
        }

        string outfile = ss.str();
        ofstream results_file(outfile);

        results_file<<"avg time expansion: " << avg_time_expansion / (double) set_seeds.size()<<endl;
        results_file<<"avg time community: " << avg_time_community / (double) set_seeds.size()<<endl;
        results_file<<"total time: " <<(avg_time_community+avg_time_expansion) <<endl;
        results_file<<"avg conductance: " << avg_conductance / (double) set_seeds.size()<<endl;
        results_file<<"avg normalized cut: " << avg_normalized_cut / (double) set_seeds.size()<<endl;

        results_file.close();

    }

    void find_communities_without_RWR(string method, vector<int> set_seeds, int k_hop, double r_max_ratio, double r_max, double epsilon){
        double avg_time_expansion = 0;
        double avg_time_community = 0;
        double avg_normalized_cut = 0;
        double avg_conductance = 0;
        for(int i=0; i<set_seeds.size(); i++){
            clock_t startExpansionTime = clock();
            int temp_seed = set_seeds[i];
            //start the seed expansion by ppr for each seed
            //vector<double> pprs_seed;
            queue<int> queue_for_traversal;
            vector<int> pprs_seeds;//with distance
            vector<int> visited;
            for(int i = 0; i<vert; i++){
                visited.push_back(false);
            }
            queue_for_traversal.push(temp_seed);
            pprs_seeds.push_back(temp_seed);
            while(queue_for_traversal.size() > 0){
                int temp_node = queue_for_traversal.front();
                queue_for_traversal.pop();
                for(int i =0; i<g.getOutSize(temp_node); i++){
                    int out_neighbor = g.outAdjList[temp_node][i];
                    if(visited[out_neighbor] == false){
                        queue_for_traversal.push(out_neighbor);
                        pprs_seeds.push_back(out_neighbor);
                        visited[out_neighbor] = true;
                    }
                }
            }


            cout<<"starting computing conductance..." <<endl;
            clock_t startCommunityTime = clock();
            //sort their ppr values (divided by deg(x);
            //compute the communities
            //vector<int> commun_set;
            unordered_set<int> commun_set;
            long volum_commun_set = 0;
            long cut_commun_set = 0;
            double best_conductance = 1.0;
            unordered_set<int> best_communset;
            best_communset.insert(temp_seed); //the first node



            for(int i = 0; i<vert * 0.05; i++){
                int temp_node = pprs_seeds[i];
                volum_commun_set += g.getOutSize(temp_node);
                for(int j=0; j<g.getOutSize(temp_node); j++){
                    //check where an out-neighbour of this node is in the set
                    int temp_out_neighbour = g.outAdjList[temp_node][j];
                    if(commun_set.count(temp_out_neighbour) == 0){
                        cut_commun_set +=1;
                    }
                    else{
                        cut_commun_set -= 1;
                    }
                }
                //cout<<"v: " <<temp_node <<" cut: " << cut_commun_set <<" vol: " << volum_commun_set <<endl;
                commun_set.insert(temp_node);
                double current_conductance = cut_commun_set / (double) min(volum_commun_set, g.m- volum_commun_set);

                if(current_conductance < best_conductance){
                    best_conductance = current_conductance;
                    best_communset = commun_set;
                }
            }

            clock_t endCommunityTime = clock();
            avg_time_community += (endCommunityTime - startCommunityTime)/ (double) CLOCKS_PER_SEC;
            cout<<"best conductance: " << best_conductance <<endl;
            avg_conductance += best_conductance;
            avg_normalized_cut += cut_commun_set / (double) volum_commun_set;
        }

        stringstream ss;
        if(method =="resacc"){
            ss << "community_detection_results/" << target_filename << "/resacc.txt";
        }
        else if(method =="power"){
            ss << "community_detection_results/" << target_filename << "/power.txt";
        }

        string outfile = ss.str();
        ofstream results_file(outfile);

        results_file<<"avg time expansion: " << avg_time_expansion / (double) set_seeds.size()<<endl;
        results_file<<"avg time community: " << avg_time_community / (double) set_seeds.size()<<endl;
        results_file<<"total time: " <<(avg_time_community+avg_time_expansion) <<endl;
        results_file<<"avg conductance: " << avg_conductance / (double) set_seeds.size()<<endl;
        results_file<<"avg normalized cut: " << avg_normalized_cut / (double) set_seeds.size()<<endl;

        results_file.close();

    }

    void to_undirected_graph(){
        stringstream ss;
        ss << "dataset/" << target_filename << "/undirected.txt";
        string outfile = ss.str();
        ofstream graph_file(outfile);

        for(int i=0; i< vert; i++){
            int temp_node = i;
            for(int j=0; j<g.getOutSize(temp_node); j++){
                int temp_out_neighbour = g.outAdjList[temp_node][j];


                if(temp_out_neighbour<temp_node){
                    unordered_set<int> out_of_temp_out_neighbour;
                    for(int l=0; l<g.getOutSize(temp_out_neighbour); l++){
                        out_of_temp_out_neighbour.insert(g.outAdjList[temp_out_neighbour][l]);
                    }
                    if(out_of_temp_out_neighbour.count(temp_node) == 0){
                        graph_file<<temp_node<<" "<< temp_out_neighbour<<endl;
                        graph_file<<temp_out_neighbour<<" "<<temp_node<<endl;
                    }
                }
                else{
                    graph_file<<temp_node<<" "<< temp_out_neighbour<<endl;
                    graph_file<<temp_out_neighbour<<" "<<temp_node<<endl;
                }

            }
        }
        graph_file.close();
    }
};

void rwr_t_PowerMethod(RWR* rwr, vector<int> nodeList, int iterations){
    return rwr->t_PowerMethod(nodeList, iterations);
}


void RWR::PowerMethodMulti(int iterations, int node_count, int num_thread){
    struct timeval t_start,t_end; 
    gettimeofday(&t_start, NULL); 
    long start = ((long)t_start.tv_sec)*1000+(long)t_start.tv_usec/1000; 
    string inputFile = "dataset/" + target_filename + "/" + target_filename + ".query";
    ifstream node_file(inputFile);
    vector<int> nodes;
    for(int i = 0; i < node_count; i++){
        int temp_node;
        node_file >> temp_node;
        if(g.getOutSize(temp_node) == 0){
            i--;
            cout << "illegal : " << temp_node << endl;
            continue;
        }
        nodes.push_back(temp_node);
    }
    node_file.close();
    if(node_count < num_thread){
        num_thread = node_count;
    }
    vector<thread> threads;
    for(int i = 0; i < num_thread-1; i++){
        vector<int> t_nodes;
        for(int j = 0; j < node_count / num_thread; j++){
            t_nodes.push_back(nodes[i * node_count / num_thread + j]);
        }
        threads.push_back(thread(rwr_t_PowerMethod, this, t_nodes, iterations));
    }
    vector<int> t_nodes;
    for(int j = 0; j < node_count / num_thread; j++){
        t_nodes.push_back(nodes[(num_thread-1) * node_count / num_thread + j]);
    }
    t_PowerMethod(t_nodes, iterations);
    for (int i = 0; i < num_thread - 1; i++){
        threads[i].join();
    }
    gettimeofday(&t_end, NULL); 
    long end = ((long)t_end.tv_sec)*1000+(long)t_end.tv_usec/1000; 
    int cost_time = end - start;

    cout << "cost: " << cost_time / (double) 1000 << endl;
    stringstream ss;
    ss << "estimated_rwr/" << target_filename << "/power.time";
    string outfile = ss.str();
    ofstream timefile(outfile);
    timefile << cost_time / (double) 1000 << endl;
    timefile.close();

}


#endif
