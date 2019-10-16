#include "provided.h"
#include "Trie.h"
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <unordered_map>
using namespace std;

class GenomeMatcherImpl
{
public:
    GenomeMatcherImpl(int minSearchLength);
    void addGenome(const Genome& genome);
    int minimumSearchLength() const;
    bool findGenomesWithThisDNA(const string& fragment, int minimumLength, bool exactMatchOnly, vector<DNAMatch>& matches) const;
    bool findRelatedGenomes(const Genome& query, int fragmentMatchLength, bool exactMatchOnly, double matchPercentThreshold, vector<GenomeMatch>& results) const;

private:
    int m_minSearchLength;
    
    struct GenomeFragment{
        GenomeFragment(int gn, int p){
            genomeNum = gn;
            pos = p;
        }
        int genomeNum;
        int pos;
    };
    
    Trie<GenomeFragment> m_trie;
    vector<Genome> m_genomes;
    
    
};

bool compare( GenomeMatch g1, GenomeMatch g2) {
    if (g1.percentMatch > g2.percentMatch)
        return true;
    else if (g1.percentMatch < g2.percentMatch)
        return false;
    
    return g1.genomeName < g2.genomeName;
}

GenomeMatcherImpl::GenomeMatcherImpl(int minSearchLength)
{
    m_minSearchLength = minSearchLength;
    
}

void GenomeMatcherImpl::addGenome(const Genome& genome)
{
    m_genomes.push_back(genome);
    for (int i = 0; i <= genome.length() - m_minSearchLength; i++){
        string frag;
        genome.extract(i, m_minSearchLength, frag);
        GenomeFragment gf(m_genomes.size() - 1, i);
        m_trie.insert(frag, gf);
    }

}

int GenomeMatcherImpl::minimumSearchLength() const
{
    return m_minSearchLength;
}

bool GenomeMatcherImpl::findGenomesWithThisDNA(const string& fragment, int minimumLength, bool exactMatchOnly, vector<DNAMatch>& matches) const
{
    if (fragment.size() < minimumLength || minimumLength < m_minSearchLength)
        return false;   //invalid input
    
    unordered_map<int, int> matchMap;   //maps genome index to dnaMatch index
    
    //find all fragments matching first few characters
    vector<GenomeFragment> potentials = m_trie.find(fragment.substr(0, m_minSearchLength), exactMatchOnly);
    for (int i = 0; i < potentials.size(); i++){
        
        bool snipsAllowed = !exactMatchOnly;        //track whether all characters must match exactly
        int gIndex = potentials[i].genomeNum;
        int gPos = potentials[i].pos;
        int start = gPos + m_minSearchLength;   //we already looked at first m_minSearchLength chars
        int matchLength = m_minSearchLength;    //track how many chars match
        
        string prefix;
        m_genomes[gIndex].extract(gPos, m_minSearchLength, prefix);
        //check if the potential match is already a snip
        if (prefix != fragment.substr(0, m_minSearchLength))
            snipsAllowed = false;   //already had our mismatch character- no more allowed
        
        //check the rest of the chars in the fragment string
        for (int j = start; j < gPos + fragment.size(); j++){
            if (j >= m_genomes[gIndex].length())
                break;  //went past end of genome sequence
            string f;
            m_genomes[gIndex].extract(j, 1, f);
            if (f == fragment.substr(j - gPos, 1))  //check one char at a time
                matchLength++;
            else if (snipsAllowed){     //one mismatch is ok
                matchLength++;
                snipsAllowed = false;   //all the rest must be the same
            }else
                break;  //more 
        }

        if (matchLength >= minimumLength){  //long enough to add to matches vector
            //check if this genome already has a match
            auto it = matchMap.find(gIndex);
            if (it == matchMap.end()){
                //genome is not currently in matches, we must add
                DNAMatch d;
                d.genomeName = m_genomes[gIndex].name();
                d.length = matchLength;
                d.position = gPos;
                matches.push_back(d);
                matchMap[gIndex] = matches.size() - 1;   //add to map for future reference
            }else if (matches[it->second].length < matchLength){   //replace with longer match
                    matches[it->second].length = matchLength;
                    matches[it->second].position = gPos;
            }
        }
    }
    
    return (matches.size() > 0);    //found at least one match
}

bool GenomeMatcherImpl::findRelatedGenomes(const Genome& query, int fragmentMatchLength, bool exactMatchOnly, double matchPercentThreshold, vector<GenomeMatch>& results) const
{
    if (fragmentMatchLength < m_minSearchLength)
        return false;
    
    unordered_map<string, int> matchCountMap;  //from genome name to num matches
    
    //number of sequences we will look at
    double s = query.length() / fragmentMatchLength;
    
    for (int i = 0; i < s; i++){
        int pos = i * fragmentMatchLength;
        //get sequence from query
        string seq;
        query.extract(pos, fragmentMatchLength, seq);
        vector<DNAMatch> matches;
        //get matches for this sequence
        findGenomesWithThisDNA(seq, fragmentMatchLength, exactMatchOnly, matches);
        //for each match, increment count
        for (int j = 0; j < matches.size(); j++){
            auto it = matchCountMap.find(matches[j].genomeName);
            if (it == matchCountMap.end()){
                //not in map yet
                string nm = matches[j].genomeName;
                matchCountMap[nm] = 1;  //found one match
            }else{
                it->second++;
            }
        }
    }
    
    //iterate through map and compute percentages
    for (auto it = matchCountMap.begin(); it != matchCountMap.end(); it++){
        double percentage = it->second / s * 100.0;
        if (percentage >= matchPercentThreshold){
            GenomeMatch gmatch;
            gmatch.genomeName = it->first;
            gmatch.percentMatch = percentage;
            results.push_back(gmatch);
        }
    }
    
    
    sort(results.begin(), results.end(), compare);
    return results.size() > 0;
}



//******************** GenomeMatcher functions ********************************

// These functions simply delegate to GenomeMatcherImpl's functions.
// You probably don't want to change any of this code.

GenomeMatcher::GenomeMatcher(int minSearchLength)
{
    m_impl = new GenomeMatcherImpl(minSearchLength);
}

GenomeMatcher::~GenomeMatcher()
{
    delete m_impl;
}

void GenomeMatcher::addGenome(const Genome& genome)
{
    m_impl->addGenome(genome);
}

int GenomeMatcher::minimumSearchLength() const
{
    return m_impl->minimumSearchLength();
}

bool GenomeMatcher::findGenomesWithThisDNA(const string& fragment, int minimumLength, bool exactMatchOnly, vector<DNAMatch>& matches) const
{
    return m_impl->findGenomesWithThisDNA(fragment, minimumLength, exactMatchOnly, matches);
}

bool GenomeMatcher::findRelatedGenomes(const Genome& query, int fragmentMatchLength, bool exactMatchOnly, double matchPercentThreshold, vector<GenomeMatch>& results) const
{
    return m_impl->findRelatedGenomes(query, fragmentMatchLength, exactMatchOnly, matchPercentThreshold, results);
}
