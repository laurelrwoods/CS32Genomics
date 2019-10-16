#ifndef TRIE_INCLUDED
#define TRIE_INCLUDED

#include <string>
#include <vector>
#include <algorithm>
#include <iostream>


template<typename ValueType>
class Trie
{
public:
    Trie();
    ~Trie();
    void reset();
    void insert(const std::string& key, const ValueType& value);
    std::vector<ValueType> find(const std::string& key, bool exactMatchOnly) const;

      // C++11 syntax for preventing copying and assignment
    Trie(const Trie&) = delete;
    Trie& operator=(const Trie&) = delete;

private:
    struct Child;
    
    struct Node {
        std::vector<ValueType> values;
        std::vector<Child> children;
    };
    
    struct Child {
        Child(char nm, Node* ptr){
            label = nm;
            nodeptr = ptr;
        }
        char label;
        Node * nodeptr;
    };
    
    Node * m_root;
    
    void clear(Node * root);
    void insertHelper(Node* root, const std::string& key, const ValueType& value);
    void findHelper(Node* root, std::vector<ValueType>& vals,
                                      const std::string& key, bool exactMatchOnly) const;
    void printHelper(Node* root, int depth);
};


/////PUBLIC MEMBER FUNCTIONS//////

template <typename ValueType>
Trie<ValueType>::Trie(){
    m_root = new Node;    //empty root node
}

template <typename ValueType>
Trie<ValueType>::~Trie(){
    clear(m_root);
}

template <typename ValueType>
void Trie<ValueType>::reset(){
    clear(m_root);    //reclaim memory
    m_root = new Node;    //empty root node
}

template <typename ValueType>
void Trie<ValueType>::insert(const std::string& key, const ValueType& value){
    insertHelper(m_root, key, value);
}

template <typename ValueType>
std::vector<ValueType> Trie<ValueType>::find(const std::string& key, bool exactMatchOnly) const{
    std::vector<ValueType> values;
    findHelper(m_root, values, key, exactMatchOnly);
    return values;
}


/////PRIVATE MEMBER FUNCTIONS////


template <typename ValueType>
void Trie<ValueType>::clear(Node * root){
    if (root == nullptr)    //done!
        return;
    for (int i = 0; i < root->children.size(); i++){
        clear(root->children[i].nodeptr);   //delete all the children's subtrees
    }
    delete root;
}

template <typename ValueType>
void Trie<ValueType>::insertHelper(Node* root, const std::string& key, const ValueType& value){
    if (key.empty()){
        root->values.push_back(value);
        return;
    }
    //check for child node with right label
    for(int i = 0; i < root->children.size(); i++){
        Child & ch = root->children[i];
        if (ch.label == key[0]){    //next node already exists, no need to create new one
            //recurse with child node and substring of key
            insertHelper(ch.nodeptr, key.substr(1), value);
            return;
        }
    }
    //no children have correct label so create new node
    Node * n = new Node;
    Child c(key[0], n);
    //add to vector of children labels
    root->children.push_back(c);
    insertHelper(n, key.substr(1), value);
}

template <typename ValueType>
void Trie<ValueType>::findHelper(Node* root, std::vector<ValueType>& vals,
                                 const std::string& key, bool exactMatchOnly) const{
    if (key.empty()){
        //found match or snip, so insert values associated with it
        vals.insert(vals.end(), root->values.begin(), root->values.end());
        return;
    }
    
    bool firstChar = false;
    if (root == m_root)
        firstChar = true;   //cant have a mismatch on first char

    
    //iterate through and find matches
    for(int i = 0; i < root->children.size(); i++){
        Child & ch = root->children[i];
        if (ch.label == key[0]){
            findHelper(ch.nodeptr, vals, key.substr(1), exactMatchOnly);
        }else if (!firstChar && !exactMatchOnly){
            //look for snips-- if this char is wrong, all the others must be right so pass in true
            findHelper(ch.nodeptr, vals, key.substr(1), true);
        }
    }
}

#endif // TRIE_INCLUDED
