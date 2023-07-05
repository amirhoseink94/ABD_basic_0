#ifndef BODY_H
#define BODY_H

#include <set>
#include <iostream>
#include <stdio.h>
#include <armadillo>
#include <bits/stdc++.h>

#include "Particle.h"
#include "Segment.h"
#include "Face.h"


using namespace std;
using namespace arma;


namespace ABD{
	enum SHAPE{SPHERE, PLANE, PLANE_STIP, PLANE_CURVED};
	class Body
	{
		public:
		// all points, segments, faces
		// points are the essential material
		vector<Particle*> points;
		// faces and segments are used for aimation and later on for the
		// collision detection, etc
		vector<Segment> segments;
		vector<Face> faces;

		Particle* center; // the centere of the sphere --just for the topological cases
		double radius;
		// all the informations for the ABD
		Mat<double> M; 									// the mass matrix
		Mat<double> M_inv;								// the inversion of mass matrix
		vec q;										// the q vector which is 12x1 vector,
		vec q_dot;								// the q dot vector which plays the role of velocity, 12x1 vector
		Mat<double> A;										// the affine transformation
		Mat<double> p;										// the translation

		double k;												// stifness parameter
		double v;												// volume parameter
		double energy_tresh_hold;				// the tresh hold for the minimum energy, we want to find the zero,

		vec q_mad;								// the current configuration in one parameter, will be used in energy calulation

		double Dl_t;									// time step, TODO: moving it to a better plave
		double Dl_t_r;

		bool movable;
		SHAPE sh;
		// initialization phase needs these two:
		Body(bool, ABD::SHAPE, vec);



		// First step:

		void apply_force();							// apply external forces to each particle
		void update_q_mad();						// update the current q_mad, to use it in calulation

		// Secons step
		/*
		Note: we have to update q_mad before using E_der and E_der_der, because the
		whole calulation is based on the q_mad.
		*/
		vec E_der(const Mat<double>&);			// calculating the derivite of the current energy
		Mat<double> E_der_der(const Mat<double>&);	//calculating the hessian matrix of current energy


		vec calculate_next_q();

		// Thirs spte: moving
		void update_A_p();							// update the affine trnasformation A and translation p based on q
		void apply_tranformation();


		//============================================================================
		bool operator==(const Face& rhs) const;
		Body& operator= (const Body& obj);

		// drawing materials for now

		void draw_body();
		void print();
		//  Geometry
		void reshape();


		double calculate_energy(const Mat<double>&);
		Mat<double> PSPD(const Mat<double>&);

		// file handling
		private:
		void load_from_file_sphere(vec);
		void load_from_file_plane(vec);
		void load_from_file_plane_stip(vec, double);
		void load_from_file_plane_curved(vec);

		void init(bool);								// the initilization, for now is just setting up dynamics
		void set_up_dynamic();					// calculate the mass matrix, initiate the q vector based on identity for A, and not moving traslation p


	};

}

#endif
