//
//  main.cpp
//  genomics
//
//  Created by Laurel Woods on 3/7/19.
//  Copyright © 2019 Laurel Woods. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <string>
#include "Trie.h"
#include "provided.h"

using namespace std;

void loadFiles(){
    int numFiles = 8;
    string files[] = {
        "Desulfurococcus_mucosus.txt", "Ferroglobus_placidus.txt", "Ferroplasma_acidarmanus.txt", "Halobacterium_jilantaiense.txt", "Halorientalis_persicus.txt", "Halorientalis_regularis.txt",  "Halorubrum_californiense.txt", "Halorubrum_chaoviator.txt"
    };
    vector<Genome> vg;
    for (int i = 0; i < numFiles; i++){
        string filename = "/Users/laurelwoods/Desktop/data/" + files[i];
        ifstream strm(filename);
        if (!strm){
            cout << "cannot open " << filename << endl;
            exit(1);
        }
        
        bool success = Genome::load(strm, vg);
        
        if (success){
            cout << "Loaded " << vg.size() << " genomes successfully from " << files[i] << endl;
        }else
            cout << "error loading genome data from " << files[i] << endl;
    }
    string f;
    //vg[0].extract(0, 200, f);
    //cout << f;
}

void testTrie(){
    Trie<int> t;
    t.insert("hi", 1);
    t.insert("him", 2);
    t.insert("hit", 3);
    t.insert("hit", 7);
    t.insert("hip", 2);
    t.insert("hip", 10);
    t.insert("lit", 0);
    t.insert("himm", 4);
    t.insert("lol", 2);
    t.insert("lo", 5);
    t.insert("lox", 0);
    t.insert("la", 9);
    
    vector<int> v = t.find("hit", false);
    for (int i = 0; i < v.size(); i++)
        cout << v[i] << endl;
    
    t.reset();
    vector<int> vec = t.find("hit", false);
    for (int i = 0; i < vec.size(); i++)
        cout << vec[i] << endl;
}

void testAddGenome(){
    Genome g1("a", "ACTG");
    Genome g2("b", "TCGACT");
    Genome g3("c", "TCTCG");
    GenomeMatcher gm(3);
    gm.addGenome(g1);
    gm.addGenome(g2);
    gm.addGenome(g3);
}

void nathan(){
    GenomeMatcher genes(10);
    vector<Genome> vec;
    ifstream input("/Users/laurelwoods/Desktop/data/Ferroplasma_acidarmanus.txt");
    Genome::load(input, vec);
    genes.addGenome(vec[0]);
}

void testFindGenome(){
    Genome g1("1", "ACTG");
    Genome g2("2", "TCGACT");
    Genome g3("3", "TCTCG");
    GenomeMatcher gm(3);
    gm.addGenome(g1);
    gm.addGenome(g2);
    gm.addGenome(g3);
    vector<DNAMatch> matches;
    gm.findGenomesWithThisDNA("CGACG", 4, false, matches);
    for (int i = 0; i < matches.size(); i++){
        cout << "Match " << i << ": " << endl;
        cout << matches[i].genomeName << " at " << matches[i].position << " of length "
        << matches[i].length << endl;
        
    }
}

void testFindGenome2(){
    Genome g1("1", "CGGTGTACNACGACTGGGGATAGAATATCTTGACGTCGTACCGGTTGTAGTCGTTCGACCGAAGGGTTCCGCGCCAGTAC");
    Genome g2("2", "TAACAGAGCGGTNATATTGTTACGAATCACGTGCGAGACTTAGAGCCAGAATATGAAGTAGTGATTCAGCAACCAAGCGG");
    Genome g3("3", "TTTTGAGCCAGCGACGCGGCTTGCTTAACGAAGCGGAAGAGTAGGTTGGACACATTNGGCGGCACAGCGCTTTTGAGCCA");
    GenomeMatcher gm(4);
    gm.addGenome(g1);
    gm.addGenome(g2);
    gm.addGenome(g3);
    
    vector<DNAMatch> matches;
    cout << gm.findGenomesWithThisDNA("GTATAT", 6, false, matches) << endl;
    for (int i = 0; i < matches.size(); i++){
        cout << "Match " << i << ": " << endl;
        cout << matches[i].genomeName << " at " << matches[i].position << " of length "
        << matches[i].length << endl;
        
    }
}

void testRelatedGenomes(){
    Genome g1("1", "ATGATCG");
    Genome g2("2", "TCGACTCGG");
    Genome g3("3", "TCTCGTTA");
    GenomeMatcher gm(2);
    gm.addGenome(g1);
    gm.addGenome(g2);
    gm.addGenome(g3);
    vector<GenomeMatch> matches;
    Genome q("query", "ACTTATTCG");
    cout << gm.findRelatedGenomes(q, 3, true, .1, matches) << endl;
    for (int i = 0; i < matches.size(); i++){
        cout << matches[i].genomeName << ": " << matches[i].percentMatch << endl;
    }
}

void test(){
    GenomeMatcher gm(3);
    Genome g1("fish", "CTGAGAGTTTCAGTC");
    Genome g2("salmon", "GGCTCATGCTCG");
    gm.addGenome(g1);
    gm.addGenome(g2);
    vector<DNAMatch> matches;
    cout << gm.findGenomesWithThisDNA("CTC", 3, false, matches) << endl;
    for (int i = 0; i < matches.size(); i++){
        cout << "Match " << i << ": " << endl;
        cout << matches[i].genomeName << " at " << matches[i].position << " of length "
        << matches[i].length << endl;
        
    }
}

//int main() {
//    //testAddGenome();
//    //loadFiles();
//    //nathan();
//    //testFindGenome2();
//    //testRelatedGenomes();
//    test();
//    cout << "\nProgram finished." << endl;
//}
