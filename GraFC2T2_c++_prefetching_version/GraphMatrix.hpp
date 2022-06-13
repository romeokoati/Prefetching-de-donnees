#ifndef DEF_GRAPHMATRIX
#define DEF_GRAPHMATRIX


#include <iostream>
#include <vector>
#include <set>
#include <algorithm>
#include <time.h>
#include <fstream>
#include <sstream>
#include <string>
#include <ctime>
#include <map>
#include "unistd.h"

//Vertex node ---------------------------------------------- -------------------------
struct Node {
	std::string name;								//Vertex name
	int id;                                     //Vertex number (starting at 0)
	std::string node_type;							//Vertex type (u=user, i=item, c=content)
	Node() {
		this->name = "";
		this->node_type = "";
		this->id = -1;
	}
	Node(std::string name, std::string node_type, int id) {
		this->name = name;
		this->node_type = node_type;
		this->id = id;
	}
	void set_(std::string name, std::string node_type, int id) {
		this->name = name;
		this->node_type = node_type;
		this->id = id;
	}
};
//Vertex node ---------------------------------------------- -------------------------

//Arc node ---------------------------------------------- -------------------------
struct Edge {
    int w_init;                                 //
    int t;                                   //
	double weight;									//The weight of the arc
	Node* tail;									//The tail node of the arc
	Node* head;									//Arc big head node

	Edge() {
		this->weight = 0.0;
		this->w_init = 0;
		this->t = 0;
		this->tail = NULL;
		this->head = NULL;
	}
	Edge(double weight, int w_init, int t, Node* head, Node* tail) {
		this->w_init = w_init;
		this->t = t;
		this->weight = weight;
		this->tail = tail;
		this->head = head;
	}
	void set_(double weight, int w_init, int t, Node* head, Node* tail) {
	    this->w_init = w_init;
		this->t = t;
		this->weight = weight;
		this->tail = tail;
		this->head = head;
	}
};
//Arc node ---------------------------------------------- -------------------------

bool operator==(const Edge& lhs, const Edge& rhs);
//Figure----------------------------------------------- ----------------------
class Graph {						//Number of sides

public:
    int nodeNum;								//Vertex number
	int edgeNum;

    std::vector<Node*> nodes;						//Vertex array
	std::vector<Edge*> edges;						//side array

public:
	//Construction method
	Graph() {
		// Initialize the property
		nodeNum = 0;
		edgeNum = 0;
	}
	// Add a node
	void add_node(std::string name, std::string node_type) {
		auto result = std::find_if(nodes.begin(), nodes.end(),
            [&name](Node* const& n) { return n->name == name; }
        );
        if (result == nodes.end()) {
            Node* node = new Node(name, node_type, nodeNum);
            nodes.push_back(node);
            this->nodeNum++;
        }
	}
	void add_edge_name(std::string head_name, std::string tail_name, double weight, int w_init, int t) {
	    bool trouve = false;
	    int k=0;
		Edge* tmp_e = new Edge();
	    while ((k<edges.size()) && (!trouve)){

			tmp_e = edges[k];
            if ((tmp_e->head->name == head_name) && (tmp_e->tail->name == tail_name)) {
                trouve = true;
            }
            k++;
        }
        if (!trouve) {
            Edge* edge = new Edge(weight, w_init, t, getNode(head_name), getNode(tail_name));
            edges.push_back(edge);
            this->edgeNum++;
        } else {
            edges[k-1]->t = t;
        }
	}
	//Add arc
	void add_edge(Node* head, Node* tail, double weight, int w_init, int t) {
        std::string head_name = head->name;
        std::string tail_name = tail->name;
        add_edge_name(head_name, tail_name, weight, w_init, t);
	}

	// Remove the node
	void removeNode(std::string name) {
		Node* node = this->getNode(name);
		if (!node) {
			return;
		}
		for (int i = 0; i < nodeNum; i++) {							// Delete all arcs
			this->removeEdge(node, nodes[i]);
		}
		for (int i = 0; i < nodeNum; i++) {							// Delete all the arcs
			this->removeEdge(nodes[i], node);
		}
		for (int i = 0; i < nodeNum; i++) {							// Move out of the array
			if (nodes[i] == node) {
				nodes.erase(nodes.begin() + i);
				break;
			}
		}
		// release the memory
		delete node;
		this->nodeNum--;
	}
	// Remove the tail to the arc of the head
	void removeEdge(Node * tail, Node * head) {
		for (int i = 0; i < edgeNum; i++) {
			//1. Find this arc
			if ((edges[i]->tail == tail) && (edges[i]->head == head)) {
				Edge* temp = edges[i];
				//2. Remove from array (use iterator to remove the elements of the specified table)
				edges.erase(edges.begin() + i);
				this->edgeNum--;
				delete temp;
			}
		}
	}
	//getter......
	Node* getNode(std::string name) {
		for (int i = 0; i < nodeNum; i++) {
			if (nodes[i]->name == name) {
				return nodes[i];
			}
		}
		return NULL;
	}
	Node* getNode(int id) {
		for (int i = 0; i < nodeNum; i++) {
			if (nodes[i]->id == id) {
				return nodes[i];
			}
		}
		return NULL;
	}

	Edge* getEdge(Node * tail, Node * head) {
		for (int i = 0; i < edgeNum; i++) {
			if (edges[i]->tail == tail && edges[i]->head == head) {
				return edges[i];
			}
		}
		return NULL;
	}
	//print
	void display() {
		//Print arc
		std::cout << "Arc tail" << " --> " << "arc head" << "   " << "weight" << "time" << "   " << std::endl;
		for (int i = 0; i < edgeNum; i++) {
			std::cout << edges[i]->tail->name << "    --> " << edges[i]->head->name << "      " << edges[i]->weight <<
			"       " << edges[i]->t << std::endl;
		}
		std::cout << std::endl;
		std::cout << std::endl;
	}
	//@Override
	void to_string() {
		std::cout << "nodeNum = " << this->nodeNum << std::endl;
		std::cout << "edgeNum = " << this->edgeNum << std::endl;
		this->display();
	}
};



struct Element {
    int row, col;
    double value;

    Element(int row, int col, double value) {
		this->row = row;
		this->col = col;
		this->value = value;
	}
};

class SparseMatrix{

    public:
    SparseMatrix(int n);                   // Square Matrix
    SparseMatrix(int rows, int cols);     // Generic Matrix

    // Get and Set
    double get(int row, int col) const;
    bool set_(int row, int col, double value);

    // Operator overloads
    SparseMatrix operator*(const SparseMatrix& mat) const;
    SparseMatrix normalize ();

    // Get information about matrix
    int rows() const;
    int cols() const;
    int size_() const;

    // Other functions
    void print(bool full = false) const;

    void transpose();                       // transpose of matrix
    SparseMatrix copy_() const;              // deep copy of matrix
    std::vector<Element*> data;          // data of matrix
private:
    int m_rows, m_cols;

};

SparseMatrix to_sparse_matrix(Graph* G);

#endif // !GraphMatrix





