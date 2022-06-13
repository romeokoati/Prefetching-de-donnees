#include "Evaluation.hpp"

using namespace std;


std::vector<std::string> metrics = {"hr","prec","recall","map","mrr","f0.5","f1","f2"};
std::vector<int> tops_n = {1,2,3,5,10,15,20,30,40,50,100};
std::vector<std::string> eval_metrics = {
    "hr@1","hr@2","hr@3","hr@5","hr@10","hr@15","hr@20","hr@30","hr@40","hr@50","hr@100",
    "prec@1","prec@2","prec@3","prec@5","prec@10","prec@15","prec@20","prec@30","prec@40","prec@50","prec@100",
    "recall@1","recall@2","recall@3","recall@5","recall@10","recall@15","recall@20","recall@30","recall@40","recall@50","recall@100",
    "map@1","map@2","map@3","map@5","map@10","map@15","map@20","map@30","map@40","map@50","map@100",
    "mrr@1","mrr@2","mrr@3","mrr@5","mrr@10","mrr@15","mrr@20","mrr@30","mrr@40","mrr@50","mrr@100",
    "f0.5@1","f0.5@2","f0.5@3","f0.5@5","f0.5@10","f0.5@15","f0.5@20","f0.5@30","f0.5@40","f0.5@50","f0.5@100",
    "f1@1","f1@2","f1@3","f1@5","f1@10","f1@15","f1@20","f1@30","f1@40","f1@50","f1@100",
    "f2@1","f2@2","f2@3","f2@5","f2@10","f2@15","f2@20","f2@30","f2@40","f2@50","f2@100"
};

vector<string> evaluation_metric_list () {
    return eval_metrics;
}


Evaluation::Evaluation(map<int, vector<int>> links_to_rec, map<int, vector<int>> rec_links){
    if (links_to_rec.size() < 0)
        cerr << "Evaluation : links_to_rec bad value!";
    if (rec_links.size() < 0)
        cerr << "Evaluation : rec_links bad value!";
    m_links_to_rec = links_to_rec;
    m_rec_links = rec_links;
    //m_rec_links_binary = NULL;
    //m_result_values = map<string, double>;
    //m_result_weights = map<string, double>;
    m_evaluation_metrics = eval_metrics;
}
Evaluation::~Evaluation(){}


vector<string> Evaluation::get_evaluation_metrics(){
    return m_evaluation_metrics;
}
map<string, double> Evaluation::get_result_values(){
    return m_result_values;
}
map<string, double> Evaluation::get_result_weights(){
    return m_result_weights;
}

void Evaluation::compute_evaluation_results(){
    //# transform all recommendation list to binary list
    vector<string> keyDist;
    map<int, vector<int>>::iterator im;

    for(im=m_links_to_rec.begin() ; im!=m_links_to_rec.end() ; im++){
        int u = (*im).first;
        std::vector<int>().swap(m_rec_links_binary[u]);
        int taille = (m_rec_links[u]).size();
        for(int i=0; i<taille; i++){
            vector<int> links;
            links = m_links_to_rec[u];
            int link = m_rec_links[u][i];
            if(find(links.begin(), links.end(), link) != links.end()){
                m_rec_links_binary[u].push_back(1);
            } else {
                m_rec_links_binary[u].push_back(0);
            }
        }
    }
//    for(im=m_links_to_rec.begin() ; im!=m_links_to_rec.end() ; im++){
//        int u = (*im).first;
//        int taille = (m_rec_links[u]).size();
//        for(int i=0; i<taille; i++)
//            std::cout << m_rec_links_binary[u][i] << "; ";
//    }
//    std::cout << "\n ";
    //# return all evaluation
    map<string, map<int, double>> eval_results, eval_result_weights;
    map<string, map<int, double>> result1 = _get_hit_ratio();
    eval_results["hr"] = result1["hr"]; eval_result_weights["hr"] = result1["w"];
    map<string, map<int, double>> result2 = _get_precision();
    eval_results["prec"]= result2["prec"];  eval_result_weights["prec"] = result2["w"];
    map<string, map<int, double>> result3 = _get_recall();
    eval_results["recall"] = result3["recall"];  eval_result_weights["recall"] = result3["w"];
    map<string, map<int, double>> result4 = _get_map();
    eval_results["map"] = result4["map"];  eval_result_weights["map"] = result4["w"];
    map<string, map<int, double>> result5 = _get_mrr();
    eval_results["mrr"] = result5["mrr"];  eval_result_weights["mrr"] = result5["w"];
    map<string, map<int, double>> result6 = _get_fmeasure(0.5);
    eval_results["f0.5"] = result6["fm"];  eval_result_weights["f0.5"] = result6["w"];
    map<string, map<int, double>> result7 = _get_fmeasure(1);
    eval_results["f1"] = result7["fm"];  eval_result_weights["f1"] = result7["w"];
    map<string, map<int, double>> result8 = _get_fmeasure(2);
    eval_results["f2"] = result8["fm"];  eval_result_weights["f2"] = result8["w"];

//     for(map<int, double>::iterator it=result2["prec"].begin() ; it!=result2["prec"].end() ; it++){
//        double u = (*it).second;
//        std::cout << "prec : " << u << std::endl;
//     }
//     for(map<int, double>::iterator it=result2["w"].begin() ; it!=result2["w"].end() ; it++){
//        double u = (*it).second;
//        std::cout << "w : " << u << std::endl;
//     }

    for(int i=0; i<metrics.size(); i++){
        string metric = metrics[i];
        for(int j=0; j<tops_n.size(); j++){
            int n = tops_n[j];
            m_result_values[metric + "@" + to_string(n)] = eval_results[metric][n];
            m_result_weights[metric + "@" + to_string(n)] = eval_result_weights[metric][n];
//            std::cout << "m_result_values : " << m_result_values[metric + "@" + to_string(n)] << std::endl;
//            std::cout << "m_result_weights : " << m_result_weights[metric + "@" + to_string(n)] << std::endl;
        }
    }
}

map<string, map<int, double>> Evaluation::_get_hit_ratio(){
    map<int, double> hri, wi;
    double nb_u = 1.0 * m_rec_links_binary.size();
    int top_n_size = tops_n.size();
    for(int i=0; i<top_n_size; i++){
        int n = tops_n[i];
        hri[n] = 0.0;
        wi[n] = nb_u;
    }
    if(nb_u > 0.0){
        for(map<int, vector<int>>::iterator im=m_rec_links_binary.begin(); im!=m_rec_links_binary.end(); im++){
            int idx = (*im).first;
            for(int j=0; j<top_n_size; j++){
                double sum = 0.0;
                int n = tops_n[j], rec_size = m_rec_links_binary[idx].size();
                int end_rec = (n < rec_size) ? n : rec_size;
                for(int top=0; top<end_rec; top++){
                    sum += m_rec_links_binary[idx][top];
                }
//                std::cout << "sum : " << sum << std::endl;
                if(sum >= 1.0)
                    hri[tops_n[j]] += 1.0;
            }
        }
        for(int i=0; i<top_n_size; i++){
            int n = tops_n[i];
            hri[n] = (hri[n] * 1.0) / nb_u;
//            std::cout << "hit ratio : " << hri[n] << std::endl;
        }
    }
    map<string, map<int, double>> result;
    result["w"] = wi;
    result["hr"] = hri;
    return result;
}
map<string, map<int, double>> Evaluation::_get_precision(){
    map<int, double> preci, wi, deno_preci, nume_preci;
    int top_n_size = tops_n.size();
    for(int i=0; i<top_n_size; i++){
        int n = tops_n[i];
        preci[n] = 0.0;
        deno_preci[n] = 0.0;
        nume_preci[n] = 0.0;
        wi[n] = 0.0;
    }
    for(map<int, vector<int>>::iterator im=m_rec_links_binary.begin(); im!=m_rec_links_binary.end(); im++){
        int idx = (*im).first;
        for(int j=0; j<top_n_size; j++){
            int n = tops_n[j];
            //# denominator [number of recommendations]
            vector<int> tmp;
            int rec_links_size = m_rec_links[idx].size();
            int end_rec_links = (n < rec_links_size) ? n : rec_links_size;
            for(int i=0; i<end_rec_links; i++)
                tmp.push_back(m_rec_links[idx][i]);
            deno_preci[n] += tmp.size();
//            std::cout << "deno_preci : " << deno_preci[n] << std::endl;

            //# numerator [number of good recommendations]
            double sum = 0.0;
            int rec_size = m_rec_links_binary[idx].size();
            int end_rec = (n < rec_size) ? n : rec_size;
            for(int top=0; top<end_rec; top++){
                sum += m_rec_links_binary[idx][top];
            }
//            std::cout << "sum : " << sum << std::endl;
            nume_preci[n] += sum;
        }
    }
    for(int i=0; i<top_n_size; i++){
        int n = tops_n[i];
        preci[n] = (deno_preci[n] > 0.0) ? ((1.0 * nume_preci[n]) / (1.0 * deno_preci[n])) : 0.0;
        wi[n] = deno_preci[n];
//        std::cout << "preci : " << preci[n] << std::endl;
//        std::cout << "wi : " << wi[n] << std::endl;
    }
    map<string, map<int, double>> result;
    result["w"] = wi;
    result["prec"] = preci;
    return result;
}

map<string, map<int, double>> Evaluation::_get_recall(){
    map<int, double> recalli, wi, deno_recalli, nume_recalli;
    int top_n_size = tops_n.size();
    for(int i=0; i<top_n_size; i++){
        int n = tops_n[i];
        recalli[n] = 0.0;
        deno_recalli[n] = 0.0;
        nume_recalli[n] = 0.0;
        wi[n] = 0.0;
    }
    for(map<int, vector<int>>::iterator im=m_rec_links_binary.begin(); im!=m_rec_links_binary.end(); im++){
        int idx = (*im).first;
        for(int j=0; j<top_n_size; j++){
            //# denominator [number of links observed]
            int n = tops_n[j];
            int u_nb_links_to_rec = m_links_to_rec[idx].size();
            deno_recalli[n] = deno_recalli[n] + u_nb_links_to_rec;

            //# numerator [number of good recommendations]
            double sum = 0.0;
            int rec_size = m_rec_links_binary[idx].size();
            int end_rec = (n < rec_size) ? n : rec_size;
            for(int top=0; top<end_rec; top++){
                sum += m_rec_links_binary[idx][top];
            }
            nume_recalli[n] += sum;
        }
    }
    for(int i=0; i<top_n_size; i++){
        //# compute the recall metric
        int n = tops_n[i];
        recalli[n] = (deno_recalli[n] > 0.0) ? ((1.0 * nume_recalli[n]) / (1.0 * deno_recalli[n])) : 0.0;
        wi[n] = deno_recalli[n];
    }
    map<string, map<int, double>> result;
    result["w"] = wi;
    result["recall"] = recalli;
    return result;
}

map<string, map<int, double>> Evaluation::_get_map(){
    double nb_u = 1.0 * (m_rec_links_binary.size());
    int top_n_size = tops_n.size();
    map<int, double> mapi, wi, nume_mapi;
    for(int i=0; i<top_n_size; i++){
        int n = tops_n[i];
        mapi[n] = 0.0;
        nume_mapi[n] = 0.0;
        wi[n] = nb_u;
    }
    if(nb_u > 0.0){
        for(map<int, vector<int>>::iterator im=m_rec_links_binary.begin(); im!=m_rec_links_binary.end(); im++){
            int idx = (*im).first;
            for(int j=0; j<top_n_size; j++){
                int n = tops_n[j];
                vector<int> tmp;
                int rec_size = m_rec_links_binary[idx].size();
                int end_rec = (n < rec_size) ? n : rec_size;
                for(int top=0; top<end_rec; top++){
                    tmp.push_back(m_rec_links_binary[idx][top]);
                }
//                for(int i=0; i<j; i++)
//                    tmp.push_back(m_rec_links_binary[idx][i]);
                nume_mapi[n] += _get_average_precision(tmp);
            }
        }
        for(int i=0; i<top_n_size; i++){
            //# compute the map metric
            int n = tops_n[i];
            mapi[n] = (1.0 * nume_mapi[n]) / (1.0 * nb_u);
//            std::cout << "mapi[n] : " << mapi[n] << std::endl;
//            std::cout << "wi : " << wi[n] << std::endl;
        }
    }
    map<string, map<int, double>> result;
    result["w"] = wi;
    result["map"] = mapi;
    return result;
}

double Evaluation::_get_average_precision(vector<int> user_rec_links_binary){
    double average_precision = 0.0;
    vector<int> indexes_good_rec;
    for(int i=0; i<user_rec_links_binary.size(); i++){
        if(user_rec_links_binary[i] == 1){
            indexes_good_rec.push_back(i);
//            std::cout << "indexes_good_rec : " << i << std::endl;
        }
    }
    for(int i=0; i<indexes_good_rec.size(); i++){
        average_precision += (double)(i + 1) / (1.0 * (double)(indexes_good_rec[i] + 1));
    }

    average_precision = (indexes_good_rec.size() > 0) ? ((1.0 * average_precision) / (1.0 * (double)indexes_good_rec.size())) : 0.0;
    return average_precision;
}

map<string, map<int, double>> Evaluation::_get_mrr(){
    double nb_u = 1.0 * (m_rec_links_binary.size());
    int top_n_size = tops_n.size();
    map<int, double> mrri, wi, nume_mrri;
    for(int i=0; i<tops_n.size(); i++){
        int n = tops_n[i];
        mrri[n] = 0.0;
        nume_mrri[n] = 0.0;
        wi[n] = nb_u;
    }
    if(nb_u > 0){
        for(map<int, vector<int>>::iterator im=m_rec_links_binary.begin(); im!=m_rec_links_binary.end(); im++){
            int idx = (*im).first;
            for(int j=0; j<top_n_size; j++){
                int n = tops_n[j];
                vector<int> tmp;
                int rec_size = m_rec_links_binary[idx].size();
                int end_rec = (n < rec_size) ? n : rec_size;
                for(int top=0; top<end_rec; top++){
                    tmp.push_back(m_rec_links_binary[idx][top]);
                }
                nume_mrri[n] += _get_reciprocal_rank(tmp);
            }
        }
        for(int i=0; i<top_n_size; i++){
            //# compute the mrr metric
            int n = tops_n[i];
            mrri[n] = nume_mrri[n]/(1.0 * nb_u);
//            std::cout << "mrri[n] : " << mrri[n] << std::endl;
//            std::cout << "wi : " << wi[n] << std::endl;
        }
    }
    map<string, map<int, double>> result;
    result["w"] = wi;
    result["mrr"] = mrri;
    return result;
}

double Evaluation::_get_reciprocal_rank(vector<int> user_rec_links_binary){
    double reciprocal_rank = 0.0;
    auto it = find(user_rec_links_binary.begin(), user_rec_links_binary.end(), 1);
    if(it != user_rec_links_binary.end()){
        int index = it - user_rec_links_binary.begin();
        reciprocal_rank = 1.0 / (1.0 * (double)(index + 1));
    }
    return reciprocal_rank;
}

map<string, map<int, double>> Evaluation::_get_fmeasure(double b){
    map<int, double> fmi, wi, true_positivei, false_negativei, false_positivei;
    int top_n_size = tops_n.size();
    for(int i=0; i<top_n_size; i++){
        int n = tops_n[i];
        fmi[n] = 0.0;
        true_positivei[n] = 0.0;
        false_negativei[n] = 0.0;
        false_positivei[n] = 0.0;
        wi[n] = 0.0;
    }
    for(map<int, vector<int>>::iterator im=m_rec_links_binary.begin(); im!=m_rec_links_binary.end(); im++){
        int idx = (*im).first;
        for(int j=0; j<top_n_size; j++){
            int n = tops_n[j];
            //# True positive [number of good recommendations]
            double sum = 0.0;
            int rec_size = m_rec_links_binary[idx].size();
            int end_rec = (n < rec_size) ? n : rec_size;
            for(int top=0; top<end_rec; top++){
                sum += m_rec_links_binary[idx][top];
            }
//            std::cout << "sum : " << sum << std::endl;
            true_positivei[n] += sum;

            //# False negative [number of links observed but which are not predicted]
            false_negativei[n] += (m_links_to_rec[idx].size() - true_positivei[n]);

            //# False positive [predicted True but which are False]
            false_positivei[n] += (n - true_positivei[n]);
        }
    }
    for(int i=0; i<top_n_size; i++){
        //# numerator of F-measure
        int n = tops_n[i];
        fmi[n] = (1 + b * b) * true_positivei[n];
        //# denominator of F-measure
        wi[n] = ((1 + b * b) * true_positivei[n]) + (b * b * false_negativei[n]) + false_positivei[n];
        fmi[n] = (wi[n] > 0) ? (1.0 * fmi[n])/(1.0 * wi[n]) : 0.0;
//        std::cout << "fmi : " << fmi[n] << std::endl;
//        std::cout << "wi : " << wi[n] << std::endl;
    }
    map<string, map<int, double>> result;
    result["w"] = wi;
    result["fm"] = fmi;
    return result;
}

