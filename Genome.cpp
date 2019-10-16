#include "provided.h"
#include <string>
#include <vector>
#include <iostream>
#include <istream>
using namespace std;

class GenomeImpl
{
public:
    GenomeImpl(const string& nm, const string& sequence);
    static bool load(istream& genomeSource, vector<Genome>& genomes);
    int length() const;
    string name() const;
    bool extract(int position, int length, string& fragment) const;
private:
    string m_name;
    string m_sequence;
};

GenomeImpl::GenomeImpl(const string& nm, const string& sequence)
{
    m_name = nm;
    m_sequence = sequence;
}

bool GenomeImpl::load(istream& genomeSource, vector<Genome>& genomes) 
{
    if (!genomeSource)
        return false;
    
    char first;
    string name;
    string sequence;
    
    genomeSource.get(first);    //check that file starts with a name
    if (first != '>')
        return false;

    while (genomeSource){
        name.clear();
        sequence.clear();
        
        getline(genomeSource, name);    //fill out name field
        if (name.empty())
            return false;   //no name- invalid
        
        char c;
        bool newline = false;
        while (genomeSource.get(c)){
            if (newline && c == '>')        //starting another genome
                break;
            if (newline && c == '\n')       //empty line is invalid
                return false;
            newline = false;    //reset
            
            switch(c){
                case 'A':
                case 'T':
                case 'C':
                case 'G':
                case 'N':
                    sequence += c;  //valid characters in sequence
                    break;
                case '\n':
                    newline = true;
                    break;
                default:
                    return false;   //invalid format
            }
        }
        
        if (sequence.empty())
            return false;
        
        Genome g(name, sequence);   //create new genome and add to vector
        genomes.push_back(g);
    }
    
    return true;  // valid file
}

int GenomeImpl::length() const
{
    return m_sequence.size();
}

string GenomeImpl::name() const
{
    return m_name;
}

bool GenomeImpl::extract(int position, int length, string& fragment) const
{
    //invalid input
    if (position < 0 || length < 0 || position + length > this->length())
        return false;
    fragment = m_sequence.substr(position, length);
    return true;    //succesfully extracted fraggment
}

//******************** Genome functions ************************************

// These functions simply delegate to GenomeImpl's functions.
// You probably don't want to change any of this code.

Genome::Genome(const string& nm, const string& sequence)
{
    m_impl = new GenomeImpl(nm, sequence);
}

Genome::~Genome()
{
    delete m_impl;
}

Genome::Genome(const Genome& other)
{
    m_impl = new GenomeImpl(*other.m_impl);
}

Genome& Genome::operator=(const Genome& rhs)
{
    GenomeImpl* newImpl = new GenomeImpl(*rhs.m_impl);
    delete m_impl;
    m_impl = newImpl;
    return *this;
}

bool Genome::load(istream& genomeSource, vector<Genome>& genomes) 
{
    return GenomeImpl::load(genomeSource, genomes);
}

int Genome::length() const
{
    return m_impl->length();
}

string Genome::name() const
{
    return m_impl->name();
}

bool Genome::extract(int position, int length, string& fragment) const
{
    return m_impl->extract(position, length, fragment);
}
