#ifndef TREE_COMPARER_HPP
#define TREE_COMPARER_HPP

#include "anc.hpp"
#include "tree_builder.hpp"
#include "anc_builder.hpp"

float DistanceUsingPearsonCorrelation(Tree& tr1, Tree& tr2);
float PartitionMetric(Tree& tr1, Tree& tr2, float threshold = 1.0);
float BranchScoreMetric(Tree& tr1, Tree& tr2, float threshold = 1.0);
float TimeWhileKAncestorsDistance(Tree& tr1, Tree& tr2);
float GetTotalBranchLength(Tree& tr);
float GetTMRCA(Tree& tr);

void PairwiseTMRCA(Tree& tr, std::vector<float>& pairwiseTMRCA);
float GetPairwiseTMRCA(Node& n, std::vector<float>& pairwiseTMRCA, std::vector<Leaves>& leaves);

#endif //TREE_COMPARER_HPP

