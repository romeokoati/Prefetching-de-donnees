#ifndef DEF_EVALUATION
#define DEF_EVALUATION

#include "RecSys.hpp"



/**###############################################################################

###############################################################################
#
# RECSYSBASIC
#
###############################################################################
#
#  EVALUATION OF RECOMMENDER SYSTEMS
#
*/

/*
typedef struct Eval_result {

} Eval_result;
*/

class Evaluation {

    public:
    Evaluation(std::map<int, std::vector<int>> links_to_rec, std::map<int, std::vector<int>> rec_links);
    std::vector<std::string> get_evaluation_metrics();
    std::map<std::string, double> get_result_values();
    std::map<std::string, double> get_result_weights();
    void compute_evaluation_results();
    ~Evaluation();

    /**#
    #
    ###########################################################################
    # private methods
    #*/
    //private:

    std::map<int, std::vector<int>> m_links_to_rec;
    std::map<int, std::vector<int>> m_rec_links;
    std::map<int, std::vector<int>> m_rec_links_binary;
    std::map<std::string, double> m_result_values;
    std::map<std::string, double> m_result_weights;
    std::vector<std::string> m_evaluation_metrics;

    /**###########################################################################
    # return [hr@5, hr@10, hr@15, hr@20, hr@30, hr@40, hr@50]
    */
    std::map<std::string, std::map<int, double>> _get_hit_ratio();
    /**#
    ###########################################################################
    # return [precision@5, precision@10, precision@15, precision@20, precision@30, precision@40, precision@50]
    */
    std::map<std::string, std::map<int, double>> _get_precision();
    /**#
    ###########################################################################
    # return [recall@5, recall@10, recall@15, recall@20, recall@30, recall@40, recall@50]
    */
    std::map<std::string, std::map<int, double>> _get_recall();
    /**#
    ###########################################################################
    # return [map@5, map@10, map@15, map@20, map@30, map@40, map@50]
    */
    std::map<std::string, std::map<int, double>> _get_map();

    double _get_average_precision(std::vector<int> user_rec_links_binary);
    /**#
    ###########################################################################
    # Mean Reciprocal Rank
    # return [mrr@5, mrr@10, mrr@15, mrr@20, mrr@30, mrr@40, mrr@50]
    */
    std::map<std::string, std::map<int, double>> _get_mrr();

    double _get_reciprocal_rank(std::vector<int> user_rec_links_binary);
    /**#
    ###########################################################################
    # F-Measure
    # return [fm@5, fm@10, fm@15, fm@20, fm@30, fm@40, fm@50]
    */
    std::map<std::string, std::map<int, double>> _get_fmeasure(double b);
};

std::vector<std::string> evaluation_metric_list();


#endif // DEF_EVALUATION
