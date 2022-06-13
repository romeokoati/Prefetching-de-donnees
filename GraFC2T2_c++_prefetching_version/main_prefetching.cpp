#include "Expgen.hpp"
#include "utils.hpp"



using namespace std;

bool sortcol( const vector<double>& v1,const vector<double>& v2 )
{
    return v1[0] < v2[0];
}

int main(int argc, char *argv[]) {
    string input_dataset_file = argv[1];
    string input_trust_network = argv[2];
    string graph_t = argv[3];

    //# ---->  evaluation metrics

    //# you can add any metric in the form prec@N (Precision), recall@N, map@N, hr@N (Hit Ratio)
    //# with N in {1, 2, 3, 5, 10, 15, 20, 30, 40, 50, 100}
    vector<string> eval_metric = {"prec@10", "recall@10", "map@10", "hr@10"};

    //# ---->  parameters of experiment

    //# the time is divided in time slices of equal length
    int number_of_time_slices = 8;

    List_param* rs1_param = new List_param(graph_t, "CIU", "LDF", "IT", 0, 0.0, 180, 0, 0.9, 0.9);
//    List_param* rs2_param = new List_param("STG", "CIU", "LDF", "IT", 540, 0.1, 365, 10, 0.7, 0.9);

    vector<List_param*> rs_list; //# list of configurated recsys
    rs_list.push_back(rs1_param);
//    rs_list.push_back(rs2_param);
    clock_t start;
    start = clock();

    //# ---->  Read data
    vector<vector<double>> linkstream(read_data(input_dataset_file, 5));
    vector<vector<double>> trust_network(read_data(input_trust_network, 2));

    // make the time field the first column of linkstream
    for (int i=0; i<linkstream.size(); i++){
        double t = linkstream[i][4], u = linkstream[i][0], i_i = linkstream[i][1], c = linkstream[i][3], r = linkstream[i][2];
        linkstream[i][0] = t;
        linkstream[i][1] = u;
        linkstream[i][2] = i_i;
        linkstream[i][3] = c;
        linkstream[i][4] = r;
    }
    sort(linkstream.begin(), linkstream.end(), sortcol);
    //std::cout << linkstream.size() << std::endl;
    main_run(rs_list, linkstream, trust_network, eval_metric, number_of_time_slices);

    start = clock() - start;

    std::cout << "Spending time: "<< (double)start/(double)CLOCKS_PER_SEC << std::endl;
    return 0;
}



