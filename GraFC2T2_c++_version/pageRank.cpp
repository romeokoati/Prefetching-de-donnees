#include "pageRank.hpp"


std::map<std::string, double> pagerank_scipy (Graph* G, std::map<std::string, double> personalization, std::map<std::string, double> dangling,
                                   int max_iter, double alpha, double tol, std::string weight) {
    std::map<std::string, double> result;
    int N = G->nodeNum, M_e = G->edgeNum;
    if(N == 0)
        return result;

    std::vector<std::string> nodes_name;

    std::vector<Node*> nodes = G->nodes;

    for (int i=0; i<N; i++){
        nodes_name.push_back(nodes[i]->name);
    }

    std::vector<double> S(N, 0.0);

    SparseMatrix M = to_sparse_matrix(G);

    SparseMatrix product = M.normalize();

    int size_product = product.size_();

    for (int i = 0; i < size_product; i++) 
    {
        S[product.data[i]->row] += 1;
    }

    //# initial vector
    std::vector<double> x(N, 1.0/N), p(N, 1.0/N), dangling_weights;
    //# Personalization vector
    if(personalization.empty())
    {

    }
    else
    {
        double sum = 0.0;
        for(int i=0; i<N; i++) {
            std::string node = nodes_name[i];
            double node_perso = personalization[node];
            if(personalization.find(node) == personalization.end()) {
                p[i] = 0.0;
            }
            else {
                p[i] = node_perso;
                sum += node_perso;
            }
        }
        for (int i = 0; i < N; i++){
            p[i] = p[i] / sum;
        }
    }

    //# Dangling nodes
    if(dangling.empty())
    {
        dangling_weights = p;
    }
    else
    {
        //# Convert the dangling dictionary into an array in nodelist order
        double sum = 0.0;
        for(std::map<std::string, double>::iterator it = dangling.begin(); it!=dangling.end(); ++it){
            if(std::find(nodes_name.begin(), nodes_name.end(), it->first) != nodes_name.end()){
                dangling_weights.push_back(it->second);
                sum += it->second;
            } else {
                dangling_weights.push_back(0.0);
            }
        }
        for (int i = 0; i < dangling_weights.size(); i++){
            dangling_weights[i] = dangling_weights[i] / sum;
        }
    }

    std::vector<int> is_dangling;
    for (int i=0; i<N; i++){
        if(S[i] == 0.0){
            is_dangling.push_back(i);
        }
    }
    
    //# power iteration: make up to max_iter iterations
    for (int i = 0; i < max_iter; i++)
    {
        std::vector<double> xlast(N, 0.0), tmp(N, 0.0);

        for(int l=0; l<N; l++)
        {
            xlast[l] = x[l];
        }
        
        double sum = 0.0, temp;
        for (int j = 0; j < is_dangling.size(); j++)
        {
            sum += x[is_dangling[j]];
        }
        
        for (int k = 0; k < size_product; k++)
        {
            int ind_ligne = product.data[k]->row;
            int ind_colone = product.data[k]->col;
            tmp[ind_colone] += x[ind_ligne] * product.data[k]->value;
        }

        for (int k = 0; k < tmp.size(); k++)
        {
            x[k] = alpha * (tmp[k] + (sum * dangling_weights[k])) + ((1.0 - alpha) * p[k]);
        }

        //# check convergence, l1 norm
        double err = 0.0;
        double tmp_er = 0.0;
        for (int k = 0; k < N; k++)
        {
            tmp_er = x[k] - xlast[k];
            err += std::abs(tmp_er);
        }

        if(err < N * tol)
        {
            for (int k = 0; k < N; k++)
            {
                result[nodes[k]->name] = x[k];
            }
            return result;
        }
    }

    for (int k = 0; k < N; k++){
        result[nodes[k]->name] = x[k];
    }
    return result;
}
