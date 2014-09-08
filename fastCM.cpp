//============================================================================
// Name        : fastCM.cpp
// Author      : Peter Windridge
// Version     :
// Copyright   : Your copyright notice
// Description : Implement the algorithm described at
//               http://mathoverflow.net/questions/180301/finding-loops-and-double-edges-asap-in-configuration-model-random-graph
//				 At present it is setup for a d regular graph.
//============================================================================

/* R code
cm<-read.csv("/tmp/CM.csv",header=F)
max(cm)
simples <- (cm$V1 == 750)
mean(cm$V1[!simples])
hist(cm$V1[!simples])
sum(simples)/length(cm$V1)  # should be around exp( (1-d^2)/4)
*/


#include <iostream>
#include <cassert>
#include <array>
#include <algorithm>
#include <numeric>
#include <random>
#include <fstream>

using namespace std;

// undefine SCORING to SELECT vertices in lexicographic order..
#define SCORING

// number of vertices
const unsigned int n = 250;
const unsigned int reps = 10000;


int main() {
	// MT random number generator
	random_device rd;
	default_random_engine mt(rd());
	uniform_real_distribution<> unif(0.0, 1.0);

	ofstream csv("/tmp/CM.csv");

	array <unsigned int,n> half_edges;

	for (unsigned int r = 0; r < reps; r++) {

		// just use a regular graph (more standard results available!)
		half_edges.fill(4);
		unsigned int total_half_edges = accumulate(half_edges.begin(), half_edges.end(), 0);
		assert (total_half_edges % 2 == 0);

		// Adjacency matrix
		unsigned int A[n][n]= {};

		// M is number of pairs completed
		unsigned int M = 0;
		unsigned int N = total_half_edges / 2;

		unsigned int SELECTd_vertex;
		unsigned int SAMPLEd_vertex;
		double U_sample;
		unsigned int cumsum = 0;

		// array to hold the badness scores (the initial bad scores are just the degrees as there are no edges formed)
		array <unsigned int,n> bad_scores(half_edges);

		unsigned int max_bad_score; //, bad_score;

		while (total_half_edges > 0)
		{
#ifdef SCORING
			// find vertex with the largest badness
			max_bad_score = 0;
			for (unsigned int i = 0; i < n; i++) {
				if (bad_scores[i]> max_bad_score) {
					max_bad_score = bad_scores[i];
					SELECTd_vertex = i;
				}
			}
#else
			// just find first vertex with at least one half-edge
			//SELECTd_vertex = find_if(half_edges.begin(), half_edges.end(), [](unsigned int h){return (h > 0);} );
			for (unsigned int i = 0; i < n; i++) {
				if (half_edges[i] > 0) {
					SELECTd_vertex = i;
					break;
				}
			}
#endif

			// if (half_edges[SELECTd_vertex] == 0) cout << SELECTd_vertex << " bad score = " << max_bad_score << endl;

			// remove a half-edge from the selected vertex
			assert(half_edges[SELECTd_vertex] > 0);
			half_edges[SELECTd_vertex]--;

			// now choose another half-edge at random
			// to do this we first sample a U(0,half-edges - 1) random variable
			U_sample = unif(mt) * (total_half_edges - 1);
			SAMPLEd_vertex = 0, cumsum = 0;
			while (cumsum < U_sample) {
				cumsum += half_edges[SAMPLEd_vertex];
				SAMPLEd_vertex ++;
				assert(SAMPLEd_vertex <= n);
			}

			SAMPLEd_vertex--;

			assert(SAMPLEd_vertex >= 0);
			assert(half_edges[SAMPLEd_vertex] > 0);
			half_edges[SAMPLEd_vertex] --;

			if (A[SELECTd_vertex][SAMPLEd_vertex] > 0) {
			//	cout << "double edge!" << endl;
				break;
			}

			if(SELECTd_vertex == SAMPLEd_vertex){
			//	cout << "loop!" << endl;
				break;
			}

#ifdef SCORING
			// update the scores
			bad_scores[SELECTd_vertex] --; // for the half-edge selected
			bad_scores[SELECTd_vertex] += half_edges[SAMPLEd_vertex]; // +half_edges at the new vertex
			bad_scores[SAMPLEd_vertex] --; // half-edge sampled
			bad_scores[SAMPLEd_vertex] += half_edges[SELECTd_vertex];

			// now for any vertex already connected to SELECTd or SAMPLEd vertex,
			// we need to decrement its score by 1.

			for(unsigned int i =0; i < n; i++) {
				if (bad_scores[i] > 0) {
					if (A[i][SAMPLEd_vertex] > 0) {
						bad_scores[i]--;
					}
					if (A[i][SELECTd_vertex] > 0) {
						bad_scores[i]--;
					}
				}
			}


			if (half_edges[SELECTd_vertex] == 0) bad_scores[SELECTd_vertex]=0;
			if (half_edges[SAMPLEd_vertex] == 0) bad_scores[SAMPLEd_vertex]=0;
#endif

			// update the adjacency matrix
			A[SELECTd_vertex][SAMPLEd_vertex] = 1;
			A[SAMPLEd_vertex][SELECTd_vertex] = 1;


			// remove the selected half-edge from the pool
			M++;
			total_half_edges -= 2;
		}

		csv << M << endl;

	}

	return 0;
}
