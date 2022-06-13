#include "utils.hpp"
#include <iostream> 
#include <string> 
#include <regex> 
#include <iterator> 
#include <bits/stdc++.h>
using namespace std; 

Out::Out (std::string file, std::string newdir) {
    std::string workdir = get_current_dir_name();
    std::string outdir = workdir + "/out/";
    m_filename = "";
    if (newdir != ""){
        m_filename = outDir + newdir;
        mkdir(m_filename.c_str(), 777);
        m_filename = m_filename + file + ".txt";
    } else
        m_filename = outDir + file + ".txt";
    m_filetmp = m_filename;

}

Out::~Out (){}

void Out::write (std::string text) {
    std::ofstream monFlux(m_filetmp.c_str(), std::ios::app);
    monFlux << text.c_str() << std::endl;
}
void Out::writewt(std::string text){ // write with time
    //std::string time;
    time_t curr_time;
	tm * curr_tm;
	char date_[100];
	char time_[100];

	time(&curr_time);
	curr_tm = localtime(&curr_time);

	strftime(date_, 50, "%d/%m/%Y", curr_tm);
	strftime(time_, 50, " %T", curr_tm);

    std::ofstream monFlux(m_filetmp.c_str(), std::ios::app);
    monFlux << date_ << time_ << " ; " << text.c_str() << std::endl;
}
void Out::close(){}
void Out::copyFile(std::string srcdir, std::string dstdir, std::string filename){
    std::string outsrcdir = outDir + srcdir;
    std::string outsrcfile = outsrcdir + m_filename;
    std::string outdstdir = outDir + dstdir;
    mkdir(outdstdir.c_str(), 777);
    std::ifstream src(outsrcfile, std::ios::binary);
    std::ofstream dst(outdstdir, std::ios::binary);

    dst << src.rdbuf();
}
DataDistCdfCcdf Out::dataDistCdfCcdf(std::vector<double> tableau, std::map<double, double> dist, std::map<double, double> cdf,
                                        std::map<double, double> ccdf){
    sort(tableau.begin(), tableau.end());
    double minv = tableau[0], maxv = tableau[tableau.size()-1];
    std::vector<double> valDist;
    std::map<double, double>::iterator im ;
    int i = 0, nb = 0, taille = tableau.size();

    for(i=0; i<taille; i++){
        ++dist[tableau[i]];
    }
    for(im=dist.begin() ; im!=dist.end() ; im++)
        valDist.push_back((double)(*im).first);

    sort(valDist.begin(), valDist.end());

    /* #CDF*/
    int len_valDist = valDist.size(), key;
    for(i=0; i<len_valDist; i++){
        key = valDist[i];
        nb += dist[key];
        cdf[key] = (100.0 * nb)/taille;
    }
    for(im=dist.begin(); im!=dist.end(); im++)
        (*im).second = (100.0 * (double)((*im).first))/taille;

    /* #CCDF*/
    for(i=0; i<len_valDist; i++){
        key = valDist[i];
        ccdf[key] = (100.0 * nb)/taille;
        nb -= dist[key];
    }
    DataDistCdfCcdf dataDistCdfCcdf;
    dataDistCdfCcdf.cdf = cdf;
    dataDistCdfCcdf.ccdf = ccdf;
    dataDistCdfCcdf.dist = dist;
    dataDistCdfCcdf.maxv = maxv;
    dataDistCdfCcdf.minv = minv;
    return dataDistCdfCcdf;

}
DistAndCcdf Out::distAndCcdf(std::vector<double> tableau, std::map<double, double> dist, std::map<double, double> ccdf){
    std::sort(tableau.begin(), tableau.end());
    double minv = tableau[0], maxv = tableau[tableau.size()-1];
    std::vector<double> valDist;
    std::map<double, double>::iterator im;
    int i = 0, nb = 0, taille = tableau.size();

    for(i=0; i<taille; i++)
        ++dist[tableau[i]];
    for(im=dist.begin() ; im!=dist.end() ; im++)
        valDist.push_back((double)(*im).first);

    std::sort(valDist.begin(), valDist.end());
    int len_valDist = valDist.size(), key;
    nb = len_valDist;

    /* #CCDF*/
    for(i=0; i<len_valDist; i++){
        key = valDist[i];
        ccdf[key] = (100.0 * nb)/taille;
        nb -= dist[key];
    }
    DistAndCcdf distAndCcdf;
    distAndCcdf.ccdf = ccdf;
    distAndCcdf.dist = dist;
    distAndCcdf.maxv = maxv;
    distAndCcdf.minv = minv;
    return distAndCcdf;
}

double tfunction_identity(double weight_init, double Dt, double nt, double ta){
    return weight_init;
}

double tfunction_half_life(double weight_init, double Dt, double nt, double ta){
    if(nt > 0.0)
        return (weight_init * exp(-(log(2) * Dt * 1.0 / nt)));
    return (weight_init * exp(-1 * Dt * log(2)));
}

double tfunction_logistic(double weight_init, double Dt, double nt, double ta){
    double K = (nt > 0.0) ? (ta * 1.0) / nt : ta;
    return (1.0 - (1.0 / (1.0 + exp(-1.0 * K * (Dt - nt)))));
}

double tfunction_constant_decay(double weight_init, double Dt, double nt, double ta){
    double tw =  (nt > 0 && Dt <= nt) ? (1 - ((Dt * 1.0)/(nt * 1.0))) : 0.0;
    return tw;
}

double tfunction_short_term(double weight_init, double Dt, double nt, double ta){
    double tw = (Dt < nt) ? 1.0 : 0.0;
    return tw;
}

std::vector<std::vector<double>> read_data(std::string path, int n) {
	std::ifstream in(path);
    std::vector<std::vector<double>> fields;
    if (in) {
        std::string line;
        while (getline(in, line)) {
            std::stringstream sep(line);
            std::string field;
            std::vector<double> row;
            while (getline(sep, field, ';')) {
                row.push_back(stod(field));
            }
            if(row.size() >= n)
                fields.push_back(row);
        }
    }
	return fields;
}

// Function to return the Jaccard distance
double jaccard_similarity_score(std::vector<double> v1, std::vector<double> v2)
{
    int size_v1 = v1.size();
    int a=0,v1_i=0, v2_i=0, i = 0, b = 0;
    for(i=0; i<size_v1; i++)
    {
        v1_i = v1[i];
        v2_i = v2[i];
        if((v1_i==v2_i) && (v1_i==1.0))
        {
            a++;
        } 
        else if (v1_i != v2_i) 
        {
            b++;
        }

    }

    // Calculate the Jaccard index using the formula
    double score = (a+b > 0) ? (double)a / (double)(a + b) : 1.0;

    return score;
}

std::vector<string> tokenizer(string s,string del = " ",int useregex=1)
{
    
    if (useregex==1)
    {
        regex regexp("(i\\()|(\\))"); 
        s =regex_replace(s, regexp, "");
    }
    
    
    std::vector<string> result ;
    int start = 0;
    int end = s.find(del);

    while (end != -1) 
    {
        // cout << s.substr(start, end - start) << endl;
        result.push_back(s.substr(start,end-start));
        
        start = end + del.size();
        end = s.find(del, start);
    }

    result.push_back(s.substr(start, end - start));
    
    return result;
}

