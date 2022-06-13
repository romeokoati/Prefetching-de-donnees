#ifndef DEF_RECSYS
#define DEF_RECSYS


#include "pageRank.hpp"


#define LINK_WEIGHT int t = (int)substream[j][0]; \
                    string u = to_string(int(substream[j][1])); \
                    string i = to_string(int(substream[j][2])); \
                    string c = to_string(int(substream[j][3])); \
                    double r = substream[j][4]; \
                    string ustr = "u"+u, istr = "i"+i, cstr = "c"+c; \
                    double ui_w = rating_to_link_weight(u, r, m_user_rating_mean, m_user_list_id, m_rating_max);

#define ADD_TO_GRAPH_NO_CONTENT m_graph->add_node(ustr, "u"); \
								m_graph->add_node(istr, "i"); \
								m_graph->add_edge_name(ustr, istr, 1, 1, t); \
								m_graph->add_edge_name(istr, ustr, 1, 1, t);

#define ADD_TO_GRAPH_NCI m_graph->add_node(cstr, "c"); \
						 m_graph->add_edge_name(cstr, istr, 1, 1, t); \
						 m_graph->add_edge_name(istr, cstr, 1, 1, t);

#define ADD_TO_GRAPH_NCU m_graph->add_node(cstr, "c"); \
						 m_graph->add_edge_name(cstr, ustr, 1, 1, t); \
						 m_graph->add_edge_name(ustr, cstr, 1, 1, t);

#define ADD_TO_GRAPH_SSTR_ISTR m_graph->add_node(sstr, "s"); \
								m_graph->add_edge_name(sstr, istr, 1, 1, t); \
								m_graph->add_edge_name(istr, sstr, 1, 1, t);

#define ADD_TO_GRAPH_CSTR_SSTR m_graph->add_node(sstr, "s"); \
								m_graph->add_edge_name(sstr, cstr, 1, 1, t); \
								m_graph->add_edge_name(cstr, sstr, 1, 1, t);

#define ADD_TO_GRAPH_LSG_NO_CONTENT m_graph->add_node(ut, "u"); \
						            m_graph->add_node(it, "i"); \
						            m_graph->add_edge_name(ut, it, 1, 1, t); \
						            m_graph->add_edge_name(it, ut, 1, 1, t); \
						            \
                                    if(m_user_last_sessions.find(ustr) != m_user_last_sessions.end()){ \
                                        string last_ut = m_user_last_sessions[ustr]; \
                                        m_graph->add_edge_name(ut, last_ut, 1, 1, t); \
                                        m_graph->add_edge_name(last_ut, ut, 1, 1, t); \
                                    } \
                                    if(m_item_last_sessions.find(istr) != m_item_last_sessions.end()){ \
                                        string last_it = m_item_last_sessions[istr]; \
                                        m_graph->add_edge_name(it, last_it, 1, 1, t); \
                                        m_graph->add_edge_name(last_it, it, 1, 1, t); \
                                    }\
                                    m_user_last_sessions[ustr] = ut; \
                                    m_item_last_sessions[istr] = it; \

#define ADD_TO_GRAPH_LSG_NCI m_graph->add_node(ct, "c"); \
								 m_graph->add_edge_name(ct, it, 1, 1, t); \
								 m_graph->add_edge_name(it, ct, 1, 1, t); \
                                if(m_content_last_sessions.find(cstr) != m_content_last_sessions.end()){ \
                                    string last_ct = m_content_last_sessions[cstr]; \
                                    m_graph->add_edge_name(ct, last_ct, 1, 1, t); \
                                    m_graph->add_edge_name(last_ct, ct, 1, 1, t); \
                                }\
                                m_content_last_sessions[cstr] = ct; \

#define ADD_TO_GRAPH_LSG_NCU m_graph->add_node(ct, "c"); \
                             m_graph->add_edge_name(ct, ut, 1, 1, t); \
                             m_graph->add_edge_name(ut, ct, 1, 1, t); \
                            if(m_content_last_sessions.find(cstr) != m_content_last_sessions.end()){ \
                                string last_ct = m_content_last_sessions[cstr]; \
                                m_graph->add_edge_name(ct, last_ct, 1, 1, t); \
                                m_graph->add_edge_name(last_ct, ct, 1, 1, t); \
                            }\
                            m_content_last_sessions[cstr] = ct; \

#define REC_LIST_BEGIN for(int i=0; i<all_items.size(); i++) {\
                            all_items_istr.push_back("i" + to_string(all_items[i])); \
                            } \
                        for (int i=0; i<nodes.size(); i++) {\
                            d[nodes[i]->name] = 0.0; \
                        } \
                        for (int j=0; j<users_to_rec.size(); j++){ \
                            vector<string> new_items, new_items_final; \
                            map<int, double> new_items_rank; \
                            vector<int> u_trusted; \
                            int u = users_to_rec[j]; \
                            string ustr = "u" + to_string(u); \
                

#define REC_LIST_TRUST if(m_user_trust.find(u) != m_user_trust.end())\
                            u_trusted = m_user_trust[u];

#define REC_LIST_RESET_PERS_VECTOR d[ustr] = 0.0; \
                    d[sstr] = 0.0; \
                    for (int i=0; i<u_trusted.size(); i++){ \
                        string uistr = "u"+to_string(u_trusted[i]); \
                        d[uistr] = 0.0; \
                        if(m_user_last_sessions.find(uistr) != m_user_last_sessions.end()){ \
                            string suistr = m_user_last_sessions[ustr]; \
                            d[suistr] = 0; \
                        } \
                    }

#define REC_LIST_SUM_D double sum=0.0; \
                        for(map<string, double>::iterator it = d.begin(); it != d.end(); ++it) \
                            sum += it->second; \
                        double sum_d = sum * 1.0; \
                        for(map<string, double>::iterator it = d.begin(); it!=d.end(); ++it){ \
                            d[it->first] = (d[it->first] * 1.0) / sum_d;\
                        }

#define REC_LIST_END    for (int i=0; i<m_rating_matrix[0].size(); i++){ \
                            if(m_rating_matrix[m_user_list_id[u]][i] == 0.0){ \
                                new_items.push_back("i" + to_string(m_id_item_list[i])); \
                            } \
                        } \
                        for (int k=0; k<new_items.size(); k++){ \
                            std::string item = new_items[k];\
                            if(rank_.find(item) != rank_.end()){ \
                                new_items_final.push_back(item); \
                            }\
                        }\
                        for (int k=0; k<new_items_final.size(); k++){ \
                            std::string item = new_items_final[k];\
                            if(rank_[item] > 1e-15){ \
                                std::string str = item.substr(1, item.length() - 1 );\
                                int id = atoi(str.c_str());\
                                new_items_rank[id] = rank_[item]; \
                            }\
                        }\
                        int top_size = (new_items_rank.size() < 100) ? new_items_rank.size() : 100; \
                        vector<pair<int, double>> top_100(top_size); \
                        partial_sort_copy(new_items_rank.begin(), new_items_rank.end(), \
                                            top_100.begin(), top_100.end(), \
                                            [](pair<const int, double> const& l, pair<const int, double> const& r) \
                                            { return l.second > r.second; }); \
                        for (int i=0; i<top_size; i++) {\
                            rec[u].push_back(top_100[i].first); \
                        } \
                        vector<string>().swap(new_items); \
                        map<int,double>().swap(new_items_rank); \
                        vector<string>().swap(new_items_final); \
                    } \
                    return rec;

#define REC_LIST_END_LSG vector<int> new_items; \
                for (int i=0; i<all_items.size(); i++) {\
                    if(find(user_item_list[u].begin(), user_item_list[u].end(), all_items[i]) != \
                       user_item_list[u].end()) {\
                       }\
                    else \
                         new_items.push_back(all_items[i]); \
                }\
                map<int, double> new_items_rank; \
                for(map<string, double>::iterator it = rank_.begin(); it!=rank_.end(); ++it){ \
                    string str = it->first;\
                    char str_i = str[0];\
                    if(str_i == 'i'){\
                        int item = stoi((tokenizer(str, ",",1))[0]); \
                        if(find(new_items.begin(), new_items.end(), item)!=new_items.end()){ \
                            if(new_items_rank.find(item) != new_items_rank.end()) {\
                                new_items_rank[item] += rank_[str]; \
                            }\
                            else {\
                                new_items_rank[item] = rank_[str]; \
                            }\
                        } \
                    }\
                }\
                for(map<int, double>::iterator it = new_items_rank.begin(); it!=new_items_rank.end(); ++it){ \
                        std::cout<<it->first<<","<<it->second<<std::endl;\
                }\
                map<int, double> new_items_rank_final; \
                for(map<int, double>::iterator it2 = new_items_rank.begin(); it2!=new_items_rank.end(); ++it2){ \
                    int item = it2->first; \
                    if(new_items_rank[item] > 1e-16) \
                        new_items_rank_final[item] = new_items_rank[item]; \
                } \
                int top_size = (new_items_rank_final.size() < 100) ? new_items_rank_final.size() : 100; \
                vector<pair<int, double>> top_100(top_size); \
                partial_sort_copy(new_items_rank_final.begin(), new_items_rank_final.end(), \
                                    top_100.begin(), top_100.end(), \
                                    [](pair<const int, double> const& l, pair<const int, double> const& r) \
                                    { return l.second > r.second; }); \
                for (int i=0; i<top_size; i++) {\
                    rec[u].push_back(top_100[i].first); \
                 } \
                vector<int>().swap(new_items); \
                map<int,double>().swap(new_items_rank_final); \
                map<int,double>().swap(new_items_rank); \
            }\
            return rec;



#define cwd getcwd
#define cd chdir



/**#
###############################################################################
###############################################################################
#
#  RECOMMENDER SYSTEMS BASIC GENERIC CLASS
#
###############################################################################
#*/

struct User_trust_map {
	int id_user;
	std::vector<int> id_trusted_users;
};


//# global information on users and items
struct Global_info {
    std::map<int, std::vector<int>> user_trust_map;
    //std::vector<User_trust_map> user_trust_map;
    std::vector<double> user_rating_mean;
    std::vector<std::vector<double>> user_jaccard_similarity;
    double rating_info[3];				// rating info [max, median, min]
    std::map<int, int> user_list_id, item_list_id;
    std::map<int, int> id_user_list, id_item_list;
    int nb_ratings;
    std::vector<std::vector<double>> rating_matrix;

};

struct List_param {
    std::string graph;
    std::string content;
    std::string time;
    std::string trust;
    int delta;
    double beta;
    int t0;
    double k;
    double gamma;
    double alpha;
    List_param(std::string graph, std::string content, std::string time, std::string trust,
               int delta, double beta, int t0, double k, double gamma, double alpha){
        this->graph = graph;
        this->content = content;
        this->time = time;
        this->trust = trust;
        this->delta = delta;
        this->beta = beta;
        this->t0 = t0;
        this->k = k;
        this->gamma = gamma;
        this->alpha = alpha;
    }
};

class Recsysgen {

    public:
    Recsysgen(int tbegin, std::vector<int> recsys_id, std::string name);
    void update_recsys(std::vector<std::vector<double>> substream, std::vector<std::vector<double>> cumulate_substream, Global_info global_info);
    std::map<int, std::vector<int>> get_recommended_list(std::vector<int> users_to_rec, std::map<int, std::vector<int>> user_item_list, std::vector<int> all_items);
    ~Recsysgen();

    //private:
    int m_tbegin;
    std::vector<int> m_recsys_id;
    std::string m_name;
};
/**#
###############################################################################
###############################################################################
#
# GRAPHRS (RECSYS GRAPH)
#
###############################################################################
###############################################################################
#
*/


class GraphRecsys : public Recsysgen {

    public:
    GraphRecsys();
    GraphRecsys(int tbegin, std::vector<int> recsys_id, std::string name, int graph_type=0, int content=0,
                 int time=0, int kp=0, int nt=0, int ta=NULL, int delta=0,
                double alpha=0, double beta=0, double k=0);
    std::string __str__();
    void update_recsys(std::vector<std::vector<double>> substream, std::vector<std::vector<double>> cumulate_substream, Global_info global_info);
    std::map<int, std::vector<int>> get_recommended_list(std::vector<int> users_to_rec, std::map<int, std::vector<int>> user_item_list, std::vector<int> all_items);
    void _time_weight(int tnow=0);
    ~GraphRecsys();

    public:

    double m_alpha=0, m_beta=0, m_k=0;
    int m_graph_type=0, m_content=0, m_time=0, m_nt=0, m_ta=NULL, m_delta=0, m_kp=0;
    Graph* m_graph;

    std::map<std::string, std::string> m_user_last_sessions, m_item_last_sessions, m_content_last_sessions;
    std::map<std::string, std::string> m_users_info, m_items_info;

    std::map<int, std::vector<int>> m_user_trust;
    //std::vector<User_trust_map> m_user_trust;
    std::vector<double> m_user_rating_mean;
    std::vector<std::vector<double>> m_user_pearson_similarity;

    double m_rating_max, m_rating_median, m_rating_min;
    std::map<int, int> m_user_list_id, m_item_list_id;
    std::map<int, int> m_id_user_list, m_id_item_list;

    std::vector<std::vector<double>> m_rating_matrix;

};

// for get_recSys result return
struct RecSys_result {
	GraphRecsys* recsys;
	std::vector<std::string> param_list;

	RecSys_result(GraphRecsys* recsys, std::vector<std::string> param_list){
        this->recsys = recsys;
		this->param_list = param_list;
	}
	//List_param param_list;
};


double rating_to_link_weight(std::string u, double r, std::vector<double> user_rating_mean, std::map<int, int> user_list_id, double max_rating);

RecSys_result* get_recsys(std::vector<int> recsys_id);






#endif
