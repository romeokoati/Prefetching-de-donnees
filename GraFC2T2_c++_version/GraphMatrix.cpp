#include "GraphMatrix.hpp"


bool operator==(const Edge& lhs, const Edge& rhs){
    return (lhs.head->name == rhs.head->name) && (lhs.tail->name == rhs.tail->name);
}

SparseMatrix::SparseMatrix(int n){
    if(n < 1)
        throw "Matrix dimentions must be positive integer.";
    this->m_rows = n;
    this->m_cols = n;
}

SparseMatrix::SparseMatrix(int rows, int cols){
    if(rows < 1 || cols < 1)
        throw "Matrix dimentions must be positive integer.";
    this->m_rows = rows;
    this->m_cols = cols;
}

double SparseMatrix::get(int row, int col) const{
    for(int i=0; i<data.size(); i++)
        if(data[i]->row == row && data[i]->col == col)
            return data[i]->value;
    return 0;
}

bool SparseMatrix::set_(int row, int col, double value){
    if(row > this->m_rows || col > this->m_cols){
        throw "Cannot insert this element. Index is too big!";
        return false;
    }
    Element* elt = new Element(row, col, value);
    data.push_back(elt);
    return true;
}

SparseMatrix SparseMatrix::operator*(const SparseMatrix& mat) const{
    if(this->cols() != mat.rows())
        throw "For multiplication, rows and cols of first matrix must match cols and rows of second matrix, respectively!";

    SparseMatrix result(this->m_rows, mat.cols());
    auto mat_t = mat.copy_();
    mat_t.transpose();

    for(int i = 0; i < this->data.size();){
        int nRow = this->data[i]->row;
        for(int j = 0; j < mat_t.data.size();){
            int nCol = mat_t.data[j]->row;
            int i_curr = i, j_curr = j;    // take a temporary pointer for each matrix and calculate this entry
            double sum = 0;
            while((i_curr < this->data.size()) && ((this->data[i_curr]->row) == nRow) &&
                  (j_curr < mat_t.data.size()) && ((mat_t.data[j_curr]->row) == nCol)) {
                if((this->data[i_curr]->col) < (mat_t.data[j_curr]->col)){
                    i_curr++;
                } else if((this->data[i_curr]->col) > (mat_t.data[j_curr]->col)){
                    j_curr++;
                } else{
                    sum += (this->data[i_curr]->value) * (mat_t.data[j_curr]->value);
//                    if(nRow==1 && nCol==2)
//                        std::cout << "mat1[i] " << this->data[i_curr]->value << ", mat2[j]" << mat_t.data[j_curr]->value << ", " << sum  << std::endl;
                    i_curr++;
                    j_curr++;
                }
            }
//            std::cout << "( " << nRow << "," << nCol << ", " << sum << ")" << std::endl;
            if (sum != 0.0)   {
                result.set_(nRow, nCol, sum);
            }
            while ((j < mat_t.data.size()) && ((mat_t.data[j]->row) == nCol))
                j++;    // the rest of matrix
        }
        while ((i < this->data.size()) && ((this->data[i]->row) == nRow))
            i++;        // the rest of other matrix
    }
    return result;
}

int SparseMatrix::rows() const {
    return this->m_rows;
}

int SparseMatrix::cols() const {
    return this->m_cols;
}

int SparseMatrix::size_() const {
    return data.size();
}

void SparseMatrix::print(bool full) const{
    if(full){
        for(int i = 0; i < m_rows; i++){
            for(int j = 0; j < m_cols; j++){
                if(this->get(i, j) != 0)
                    std::cout << this->get(i,j) << "    ";
                else
                    std::cout << "0" << "   ";
            }
            std::cout << std::endl;
        }
    }
    else{
        std::cout << "(ROW, COL, VAlUE)" << std::endl;
        for(unsigned int i=0; i< this->data.size(); ++i){
            std::cout << "(" << this->data[i]->row << ", " << this->data[i]->col << ", " << this->data[i]->value << ")" << std::endl;
        }
    }
}

void SparseMatrix::transpose()
{
    for(unsigned int i=0; i< this->data.size(); ++i)
        std::swap(this->data[i]->row, this->data[i]->col);
    std::sort(this->data.begin(), this->data.end(), [](const Element* e1, const Element* e2) { return e1->row < e2->row; } );
}

SparseMatrix SparseMatrix::copy_() const{
    SparseMatrix _copy(this->m_rows, this->m_cols);
    _copy.data = this->data;
    return _copy;
}

SparseMatrix SparseMatrix::normalize ()
{
    int M = this->size_();
    SparseMatrix result(this->rows(), this->cols());
    std::vector<double> S(this->cols(), 0.0);
    // Sum element over the column axis.
    for (int i=0; i<M; i++)
    {
        Element* elt = this->data[i];
        int n_1_id = elt->row;
        S[n_1_id] += elt->value;
    }
    for (int i=0; i<M; i++)
    {
        Element* elt = this->data[i];
        int n_1_id = elt->row;
        int n_2_id = elt->col;
        result.set_(n_2_id, n_1_id, (elt->value / S[n_2_id]));
    }
    return result;
}

SparseMatrix to_sparse_matrix(Graph* G)
{
    int N = G->nodeNum, M = G->edgeNum;
    SparseMatrix agacency_mat(N);
    for (int i=0; i<M; i++)
    {
        Edge* edge = G->edges[i];
        int n_1_id = edge->head->id;
        int n_2_id = edge->tail->id;
        agacency_mat.set_(n_1_id, n_2_id, edge->weight);
    }
    return agacency_mat;
}
