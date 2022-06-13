#ifndef DEF_EXPGEN
#define DEF_EXPGEN

#include "utils.hpp"
#include "Evaluation.hpp"


/**
#
###############################################################################
###############################################################################
#
# EXPGEN
#
###############################################################################
########################################################
*/
class Expgen {

    public:
    Expgen(std::vector<std::vector<double>> linkstream, int slice_length,
            std::vector<GraphRecsys*> recsys_list, std::vector<std::vector<double>> trust_network, double rating_max);
    void save(std::string filename, std::string resdir="");
    std::map<std::string, std::map<std::string, double>> get_exp_results();
    void run();

    //private:
    std::vector<std::vector<double>> m_linkstream;
    int m_slice_length;
    std::vector<GraphRecsys*> m_recsys_list;

    int m_tmin, m_tmax;

    std::vector<std::vector<double>> m_trust_network;
    double m_rating_max, m_rating_median;

    std::map<std::string, std::map<std::string, double>> m_exp_results;
};


void main_run(std::vector<List_param*> rs_list, std::vector<std::vector<double>> linkstream, std::vector<std::vector<double>> trust_network,
              std::vector<std::string> eval_metric, int number_of_time_slices);


#endif // DEF_EXPGEN


