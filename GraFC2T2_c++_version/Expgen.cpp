#include "Expgen.hpp"


using namespace std;

Expgen::Expgen (vector<vector<double>> linkstream, int slice_length,
                vector<GraphRecsys*> recsys_list, vector<vector<double>> trust_network, double rating_max) {
   
    m_linkstream = linkstream;
    m_slice_length = slice_length;
    m_recsys_list = recsys_list;
    m_tmin = linkstream[0][0];
    m_tmax = linkstream[linkstream.size()-1][0];
    m_trust_network = trust_network;
    m_rating_max = rating_max;
    m_rating_median = rating_max/2.0;

}

void Expgen::save(string filename, string resdir) {
    Out save_file = Out(filename, resdir);
    string exp_results = "";

    save_file.write(to_string(m_linkstream.size()));
    save_file.write(to_string(m_recsys_list.size()));
    save_file.write(to_string(m_tmin));
    save_file.write(to_string(m_tmax));
}

map<string, map<string, double>> Expgen::get_exp_results () {
    return m_exp_results;
}

void Expgen::run(){
    int exp_tmin = m_tmin, slice_iter = 1, link_iter = 0, nb_all_links = m_linkstream.size();
    std::vector<int> user_list, item_list, content_list;
    std::vector<double> last_user_rating_mean;
    std::map<int, std::vector<int>> last_user_item;
    std::vector<int> last_item_list;
    std::map<int, int> last_user_list_id;

    std::map<std::string, std::map<int, std::map<std::string, std::map<std::string, double>>>> results_iter;
    std::map<int, std::vector<int>> user_trust_global;
    std::map<int, std::vector<int>>::iterator im, im2;
    std::vector<std::vector<double>> cumulate_substream;

    //# The initial time corresponds to the time of appearance of the first link of linkstream
    for(int i=0; i<m_recsys_list.size(); i++){
        m_recsys_list[i]->m_tbegin = exp_tmin;
        //results_iter[m_recsys_list[i].m_name] = {}
    }
    for(int i=0 ; i<m_trust_network.size() ; i++){
        user_trust_global[m_trust_network[i][0]].push_back(m_trust_network[i][1]);
    }

    while(link_iter<nb_all_links){
        //# Extract substream and links to recommend
        std::vector<std::vector<double>> substream;
        int t, u, i, c;
        double r;
        int slice_tmax = exp_tmin + (slice_iter * m_slice_length), number_of_links_rec = 0;
        std::vector<int> slice_new_user_list, slice_new_item_list, slice_new_content_list;
        std::map<int, std::vector<int>> links_to_rec, slice_new_user_item;

        while(link_iter<nb_all_links){
            std::vector<double> link = m_linkstream[link_iter];
            bool b_old_u = false, b_old_i = false;
            t = (int)link[0];
            u = (int)link[1];
            i = (int)link[2];
            c = (int)link[3];
            r = link[4];
            if(t > slice_tmax)
                break;
            substream.push_back(link);
            cumulate_substream.push_back(link);
            if(std::find(user_list.begin(), user_list.end(), u) != user_list.end())
                b_old_u = true;
            else
                slice_new_user_list.push_back(u);

            if(std::find(item_list.begin(), item_list.end(), i) != item_list.end())
                b_old_i = true;
            else
                slice_new_item_list.push_back(i);

            if(std::find(content_list.begin(), content_list.end(), c) == content_list.end())
                slice_new_content_list.push_back(c);

            if((b_old_u==true) && (b_old_i==true)){
                if (last_user_item.find(u) != last_user_item.end()){
                    if(std::find(last_user_item[u].begin(), last_user_item[u].end(), i) == last_user_item[u].end()){
                        if((r >= m_rating_median) && (r >= last_user_rating_mean[last_user_list_id[u]])){
                            links_to_rec[u].push_back(i);
                            number_of_links_rec += 1;
                        }
                    }
                }
            }
            slice_new_user_item[u].push_back(i);
            link_iter++;
        }
        /* ###################################################################
            # we have the substream of the the current time_slice
        */
        std::cout << "len(substream) " << substream.size() <<" | len(cumulate_substream) " << cumulate_substream.size() << std::endl;

        user_list.reserve(user_list.size() + distance(slice_new_user_list.begin(), slice_new_user_list.end()));
        user_list.insert(user_list.end(), slice_new_user_list.begin(), slice_new_user_list.end());
        item_list.reserve(item_list.size() + distance(slice_new_item_list.begin(), slice_new_item_list.end()));
        item_list.insert(item_list.end(), slice_new_item_list.begin(), slice_new_item_list.end());
        content_list.reserve(content_list.size() + distance(slice_new_content_list.begin(), slice_new_content_list.end()));
        content_list.insert(content_list.end(), slice_new_content_list.begin(), slice_new_content_list.end());

        std::set<int> s(user_list.begin(), user_list.end());
        user_list.assign(s.begin(), s.end());
//        reverse(user_list.begin(), user_list.end());
        std::set<int> s2(item_list.begin(), item_list.end());
        item_list.assign(s2.begin(), s2.end());
//        reverse(item_list.begin(), item_list.end());
        std::set<int> s3(content_list.begin(), content_list.end());
        content_list.assign(s3.begin(), s3.end());

        //# creation of integer user_id and integer item_id for similarity_matrix
        std::map<int, int> user_list_id, id_user_list, item_list_id, id_item_list, content_list_id, id_content_list;

        int u_id = 0; int user_list_i = 0;
        for(int i=0; i<user_list.size(); i++) {
            user_list_i = user_list[i];
            user_list_id[user_list_i] = u_id;
            id_user_list[u_id] = user_list_i;
            u_id++;
        }
        int i_id = 0; int item_list_i = 0;
        for(int i=0; i<item_list.size(); i++) {
            item_list_i = item_list[i];
            item_list_id[item_list_i] = i_id;
            id_item_list[i_id] = item_list_i;
            i_id++;
        }
        int c_id = 0; int content_list_i = 0;
        for(int i=0; i<content_list.size(); i++) {
            content_list_i = content_list[i];
            content_list_id[content_list_i] = c_id;
            id_content_list[c_id] = content_list_i;
            c_id++;
        }
        // # creation of information matrix
        int nb_user = user_list.size(), nb_item = item_list.size();
        std::vector<std::vector<double>> rating_matrix(nb_user, std::vector<double>(nb_item, 0.0));
        std::map<int, std::vector<int>> item_content;
        for(int i=0; i<cumulate_substream.size(); i++) {
            int u = (int)cumulate_substream[i][1], i_i = (int)cumulate_substream[i][2];
            rating_matrix[user_list_id[u]][item_list_id[i_i]] = cumulate_substream[i][4];

            // # for item contents
            item_content[i_i].push_back((int)cumulate_substream[i][3]);
        }
        // # for user trust
        std::map<int, std::vector<int>> user_trust;
        for(std::map<int,std::vector<int>>::iterator it = user_trust_global.begin(); it != user_trust_global.end(); ++it) {
            std::set_intersection((it->second).begin(),(it->second).end(),
                                    user_list.begin(),user_list.end(),
                                    back_inserter(user_trust[it->first]));
        }
        //  # user information
        std::vector<double> user_rating_mean(nb_user, 0.0), divide(nb_user, 0.0);
        for (int i=0; i<nb_user; i++){
            for (int j=0; j<nb_item; j++){
                user_rating_mean[i] += rating_matrix[i][j];
                divide[i] = (rating_matrix[i][j] == 0.0) ? divide[i] : (divide[i] + 1.0);
            }
        }
        for (int j=0; j<nb_user; j++)
            user_rating_mean[j] = (user_rating_mean[j] == 0.0) ? user_rating_mean[j] : (user_rating_mean[j] / divide[j]);

        std::vector<std::vector<double>> implicit_rating_matrix(nb_user, std::vector<double>(nb_item, 0.0));
        std::vector<std::vector<double>> user_jaccard_similarity(nb_user, std::vector<double>(nb_user, 0.0));
        int u_i, i_i;
        // compute implicit Rating matrix to be used for Jaccard similarity computation.
        for(u_i=0; u_i<nb_user; u_i++){
            for(i_i=0; i_i<nb_item; i_i++){
                double rrr = rating_matrix[u_i][i_i];
                if((rrr >= m_rating_median) && (rrr >= user_rating_mean[u_i]))
                    implicit_rating_matrix[u_i][i_i] = 1.0;
                else
                    implicit_rating_matrix[u_i][i_i] = 0.0;
            }
        }
        // Compute user jaccard similarity matrix.
        for(int u_i_1=0; u_i_1<nb_user; u_i_1++)
        {
            for(int u_i_2=u_i_1; u_i_2<nb_user; u_i_2++)
            {
                double sim_u1_u2 = jaccard_similarity_score(implicit_rating_matrix[u_i_1],implicit_rating_matrix[u_i_2]);
                
                user_jaccard_similarity[u_i_1][u_i_2] = sim_u1_u2;
                user_jaccard_similarity[u_i_2][u_i_1] = sim_u1_u2;
            }
        }

        //# global information on users and items
        Global_info global_info;
        //# users global information
        global_info.user_trust_map = user_trust;
        global_info.user_rating_mean = user_rating_mean;
        global_info.user_jaccard_similarity = user_jaccard_similarity;

        //# rating info [max, median, min]
        global_info.rating_info[0] = m_rating_max;
        global_info.rating_info[1] = m_rating_max/2.0;
        global_info.rating_info[2] =  0.0;

        //# user_list_id , item_list_id and content_list_id
        global_info.user_list_id = user_list_id;
        global_info.item_list_id = item_list_id;
        global_info.id_user_list = id_user_list;
        global_info.id_item_list = id_item_list;

        //# number of distinct ratings
        int count_ = 0;
        for (int i=0; i<nb_user; i++)
            for (int j=0; j<nb_item; j++)
                if(rating_matrix[i][j] != 0.0)
                    count_++;
        global_info.nb_ratings = count_;
        
        global_info.rating_matrix = rating_matrix;
        /**#
        ###################################################################
        #
        # applied recommender systems for the current time slice
        #**/
        for(int i=0; i<m_recsys_list.size(); i++) {
            std::vector<int> keys;
            for(std::map<int,std::vector<int>>::iterator it = links_to_rec.begin(); it != links_to_rec.end(); ++it) {
              keys.push_back(it->first);
            }

            std::map<int, std::vector<int>> rec_links = m_recsys_list[i]->get_recommended_list(keys, last_user_item, last_item_list);

            Evaluation recsys_eval(links_to_rec, rec_links);
            recsys_eval.compute_evaluation_results();
            results_iter[m_recsys_list[i]->m_name][slice_iter]["value"] = recsys_eval.get_result_values();
            results_iter[m_recsys_list[i]->m_name][slice_iter]["weight"] = recsys_eval.get_result_weights();

            // # update recsys for next step

            m_recsys_list[i]->update_recsys(substream, cumulate_substream, global_info);

        }
        /*#
        ###################################################################
        #
        # preparation of the next time slice
        #*/
        last_user_rating_mean = user_rating_mean;
        last_user_list_id = user_list_id;
        last_item_list = item_list;
        std::vector<int> keys_item;
        
        for(std::map<int,std::vector<int>>::iterator it = slice_new_user_item.begin(); it != slice_new_user_item.end(); ++it) {
            int u = it->first;
            if(last_user_item.find(u) != last_user_item.end()){

                last_user_item[u].reserve(last_user_item[u].size() + distance(slice_new_user_item[u].begin(), slice_new_user_item[u].end()));
                last_user_item[u].insert(last_user_item[u].end(), slice_new_user_item[u].begin(), slice_new_user_item[u].end());
                std::set<int> s3(last_user_item[u].begin(), last_user_item[u].end());
                last_user_item[u].assign(s3.begin(), s3.end());
            } else{
                last_user_item[u] = slice_new_user_item[u];
            }
        }
        slice_iter += 1;
        //# supression of objects that were create during this iteration
        std::vector<std::vector<double>>().swap(implicit_rating_matrix);
        std::vector<std::vector<double>>().swap(user_jaccard_similarity);

        std::vector<std::vector<double>>().swap(substream);
        std::map<int,std::vector<int>>().swap(links_to_rec);
        number_of_links_rec = 0;
        std::vector<int>().swap(slice_new_user_list);
        std::vector<int>().swap(slice_new_item_list);
        std::map<int,std::vector<int>>().swap(slice_new_user_item);
    }
    /**#
    #######################################################################
    #
    #  compute time averaged values of any evaluation metric
    #*/
    std::vector<std::string> eval_metrics(evaluation_metric_list());
    for(int i=0; i<m_recsys_list.size(); i++) {
        std::string r_name = m_recsys_list[i]->m_name;
        for(int j=0; j<eval_metrics.size(); j++) {
            double numerator = 0.0, denominator = 0.0;
            string metric_j = eval_metrics[j];
            for(std::map<int, std::map<std::string, std::map<std::string, double>>>::iterator it = results_iter[r_name].begin();
                    it != results_iter[r_name].end(); ++it) {
                numerator += (results_iter[r_name][it->first]["value"][metric_j] *
                              results_iter[r_name][it->first]["weight"][metric_j]);
                denominator += results_iter[r_name][it->first]["weight"][metric_j];

            }
            m_exp_results[r_name][metric_j] = (denominator > 0.0) ? (numerator/denominator) : 0.0;
        }
    }

    //# supression of created dobjects
    std::vector<int>().swap(user_list);
    std::vector<int>().swap(item_list);
    std::map<int,std::vector<int>>().swap(last_user_item);
    std::map<std::string, std::map<int, std::map<std::string, std::map<std::string, double>>>>().swap(results_iter);

}

void main_run(std::vector<List_param*> rs_list, std::vector<std::vector<double>> linkstream, std::vector<std::vector<double>> trust_network,
              std::vector<std::string> eval_metric, int number_of_time_slices){
    double tmin = linkstream[0][0];
    double tmax = linkstream[linkstream.size()-1][0];
    if(number_of_time_slices <= 2)
        number_of_time_slices = 2;
    int slice_length = (int)((tmax - tmin + (double)number_of_time_slices) / (double)number_of_time_slices);

    std::map<std::string, int> dict_graph, dict_content, dict_time, dict_trust;
    dict_graph["BIP"] = 0; dict_graph["STG"] = 1; dict_graph["LSG"] = 2;
    dict_content["NA"] = 0; dict_content["CI"] = 1; dict_content["CIU"] = 3;
    dict_time["NA"] = 0; dict_time["EDF"] = 1; dict_time["LDF"] = 2;
    dict_trust["NA"] = 0; dict_trust["ET"] = 1; dict_trust["IT"] = 2;

    std::vector<GraphRecsys*> recsys_list;
    for(int i=0; i<rs_list.size(); i++){
        std::string rs_name = rs_list[i]->graph + "-" + rs_list[i]->content + "-" + rs_list[i]->time + "-" + rs_list[i]->trust;
        std::vector<int> rs_id(4, 0);
        rs_id[0] = dict_graph[rs_list[i]->graph];
        rs_id[1] = dict_content[rs_list[i]->content];
        rs_id[2] = dict_time[rs_list[i]->time];
        rs_id[3] = dict_trust[rs_list[i]->trust];

        RecSys_result* result = get_recsys(rs_id);
        GraphRecsys* recsys = result->recsys;
        std::vector<std::string> param_list = result->param_list;
        std::string suffix = "";

        for(int j=0; j<param_list.size(); j++){
            if(param_list[j]=="delta"){
                recsys->m_delta = rs_list[i]->delta * 86400;
                suffix += "-Delta-" + std::to_string(rs_list[i]->delta);
            } else if(param_list[j]=="beta"){
                recsys->m_beta = rs_list[i]->beta;
                suffix += "-Beta-" + std::to_string(rs_list[i]->beta);
            } else if(param_list[j]=="nt"){
                recsys->m_nt = rs_list[i]->t0 * 86400;
                suffix += "-T0-" + std::to_string(rs_list[i]->t0);
            } else if(param_list[j]=="ta"){
                recsys->m_ta = rs_list[i]->k;
                suffix += "-K-" + std::to_string(rs_list[i]->k);
            } else if(param_list[j]=="k"){
                recsys->m_k = rs_list[i]->gamma;
                suffix += "-Gamma-" +  std::to_string(rs_list[i]->gamma);
            } else if(param_list[j]=="alpha"){
                recsys->m_alpha = rs_list[i]->alpha;
                suffix += "-Alpha-" +  std::to_string(rs_list[i]->alpha);
            }
        }
        recsys->m_name = rs_name + suffix;
        recsys_list.push_back(recsys);
    }
    int timestamp = (int)(time(0));
    std::string resdir = "GraFC2T2";
    std::string savefile = resdir + "-" + std::to_string(timestamp);

    std::cout << "BEGIN : " << savefile << std::endl;

    Expgen param_exp(linkstream, slice_length, recsys_list, trust_network, 5.0);

    param_exp.run();

    auto recsys_results = param_exp.get_exp_results();

    std::string workdir = get_current_dir_name();
    std::string recsys_results_file_save = workdir + "/" + resdir + "-RESULTS-" + std::to_string(timestamp) + ".txt";

    std::ofstream monFlux(recsys_results_file_save.c_str(), std::ios::app);

    for(auto it = recsys_results.cbegin(); it != recsys_results.cend(); ++it) {
        std::string line = "\n----> " + (it->first) + "\n";
         monFlux << line.c_str() << std::endl;
        for (int i = 0; i < eval_metric.size(); i++){
            line = eval_metric[i] + " = " + to_string(recsys_results[it->first][eval_metric[i]]) + "\n";
             monFlux << line.c_str() << std::endl;
        }
    }
    std::cout << "END : " << savefile << std::endl;
}




