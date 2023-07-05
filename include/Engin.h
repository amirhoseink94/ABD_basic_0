#ifndef ENGIN_H
#define ENGIN_H

#include <set>
#include <iostream>
#include <stdio.h>
#include <armadillo>
#include <bits/stdc++.h>

#include "Particle.h"
#include "Segment.h"
#include "Face.h"
#include "Body.h"

#include "../CCD/ccd.hpp"
#include "../CCD/timer.hpp"
#include "../CCD/interval_root_finder.hpp"

using namespace std;
using namespace arma;


namespace ABD{
	class Engin
	{
		public:
		// All objects, as the body type, that would be nice :D
		vector<Body*> dynamic_objects;
		vector<Body*> static_objects;
		unsigned int N;
		vec Q;												// stacked Q =(q_1, q_2,..., q_\ell)
		vec dQ;												// derivitive of Q
		mat H;												// hessian mattix

		double Dl_t;									// time step, TODO: moving it to a better plave
		double energy_tresh_hold;
		double distance_tresh_hold;


		// initialization phase needs these two:
		void init();										// the initilization,

		void add_object(Body*);

		// First step:

		void apply_force();							// apply external forces to each object, for now, it is just the particles applyforce


		vec calculate_next_Q();

		// Thirs spte: moving
		void update_A_p();							// update the affine trnasformation A and translation p based on q, it seems we do not need this
		void apply_tranformation();

		void apply_tranformation_temp(const vec&);

		double calculate_energy(const vec&);
		double calculate_energy_collision(const vec&, vec&);


		// tools
		mat PSPD(const mat&);					// normalizing the hessian matrix to be semipositive

		double distance_VF(Particle&, const Face&, vec&, vec&);
		double distance_EE(const Segment&, const Segment&, vec&, vec&);

		double distance_PF(vec&, const Face&);

		float contact_detection_topological(vec&, vec&);
		float contact_detection_topological2(vec&, vec&);
		float contact_detection_exact(vec&, vec&);

	private:
		bool contact_edge_edge(const vec&, const vec&, const vec&, const vec&,
                           const vec&, const vec&, const vec&, const vec&, float&);
		bool contact_vertex_face(const vec&, const vec&, const vec&, const vec&,
                             const vec&, const vec&, const vec&, const vec&, float&);
	};
}
#endif
