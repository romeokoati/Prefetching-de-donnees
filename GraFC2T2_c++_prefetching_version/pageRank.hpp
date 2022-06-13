#ifndef DEF_PAGERANK
#define DEF_PAGERANK

#include "GraphMatrix.hpp"
//#include "utils.h"


/**#
###############################################################################
###############################################################################
#
# PAGERANK
#
###############################################################################
###############################################################################
#
**/
std::map<std::string, double> pagerank_scipy (Graph* G, std::map<std::string, double> personalization, std::map<std::string, double> dangling,
                                   int max_iter=100, double alpha=0.85, double tol=1.0e-6, std::string weight="weight");


#endif // DEF_PAGERANK

