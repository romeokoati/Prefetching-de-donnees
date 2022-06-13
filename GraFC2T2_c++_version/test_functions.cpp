#include <iostream>
#include <string>
#include <map>


using namespace std;

static long string_hash(pair<string, int> node_hash)
{
    register int len;
    register unsigned char *p;
    register long x;
    char ob_sval[1];
    ob_sval[0] = 3;

    if (node_hash.second != -1)
        return node_hash.second;

    len = node_hash.first.size();
    if (len == 0) {
        node_hash.second = 0;
        return 0;
    }
    p = (unsigned char *)(ob_sval);
    std::cout << "node_hash.first.size(): " << node_hash.first.size() << std::endl;
    x = *p << 7;
    std::cout << "node_hash : " << *p << std::endl;
    while (--len >= 0)
        x = (1000003*x) ^ *p++;
    x ^= node_hash.first.size();
    if (x == -1)
        x = -2;
    node_hash.second = x;
    return x;
}




