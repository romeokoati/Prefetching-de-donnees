#include <vector>

void jaccard_similarity_score(double* v1, double* v2, double* intersect, int n){
    for(int i=0; i<n; i++)
        for(int j=0; j<n; j++)
            if(v1[i]==v2[j])
                intersect[i] = v1[i];
}
