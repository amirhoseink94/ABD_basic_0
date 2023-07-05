// openGL libraries
#include <GL/gl.h>
#include <GLFW/glfw3.h>
#include <GL/glut.h>
// armadillo library
#include <armadillo>
// system and c++ libraries
#include <stdio.h>
#include <iostream>
#include <cmath>
#include <time.h>
#include <vector>

// my libraries
#include "Pair.h"
#include "Vector.h"
#include "Shape.h"

using namespace std;
typedef Pair<Vec3f> Segment;
typedef Vector<Vec3f> Face;




// global variables

int width = 500, height = 500;   // Size of the drawing area, to be set in reshape().
int frameNumber = 0;     // For use in animation.
int rotate = 0;

int zoom = 45.;
int scr_width = 800;
int scr_height = 800;

int x_camera = 0.;
int y_camera = 0.;
int z_camera = 0.;

int x_rotate_angle = 0;
int y_rotate_angle = 0;

int azimuth = -300.;
int elevation = 0.;

bool right_button_flag = false;
bool left_button_flag = false;
int oldX = 0.;
int oldY = 0.;

// Animations settings
int frame_no = 0;
double dl_t = 1/2.;
double G = 0.0001;
double tresh_hold = 5;
// key handlers
bool UP_key_flag = false;
bool DOWN_key_flag = false;
bool RIGHT_key_flag = false;
bool LEFT_key_flag = false;

bool n_key_flag = false;
bool m_key_flag = false;
bool w_key_flag = false;
bool s_key_flag = false;
bool a_key_flag = false;
bool d_key_flag = false;
bool z_key_flag = false;

void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods);
void key_handler();



// rednders and drawing
void render();
//void render_particles(vector<Particle>&);
void draw_cube();
void draw_space();

// body maker
Shape make_triangle(Vec3f, Vec3f, Vec3f);
// mesh maker
Shape make_mesh(Shape start, int n);
Shape make_mesh_handler(Shape start, int itr, int n);
Shape merge_shapes(Shape s, Shape v);



int main()
{
	if ( !glfwInit() )
		return -1;
	cout << "(the glfw just initilized!)" << endl;

	Vec3f A(-20, 0 , -20);
	Vec3f B(20, 0 , 20); //
	Vec3f C(-20, 0 , 20); //
	Vec3f D(20, 0, -20);
	Vec3f E(20, 0, 20); //
	Vec3f F(-20, 0, 20); //


	Shape ABC = make_triangle(A, B, C);
	Shape ABD = make_triangle(A, B, D);
	Shape ABC_m = make_mesh(ABC, 4);
	Shape ABD_m = make_mesh(ABD, 4);
	Shape s = merge_shapes(ABC_m, ABD_m);

	/*Vec3f A(0, 10 , 0);
	Vec3f B(-10, 0 , -10); //
	Vec3f C(10, 0 , -10); //
	Vec3f D(0, -10, 0);
	Vec3f E(10, 0, 10); //
	Vec3f F(-10, 0, 10); //
	hape ABC = make_triangle(A, B, C);
	Shape BCD = make_triangle(B, C, D);
	Shape ABF = make_triangle(A, B, F);
	Shape BDF = make_triangle(B, D, F);

	Shape AEC = make_triangle(A, E, C);
	Shape ECD = make_triangle(E, C, D);
	Shape AEF = make_triangle(A, E, F);
	Shape EDF = make_triangle(E, D, F);




	Shape ABC_m = make_mesh(ABC, 2);
	Shape BCD_m = make_mesh(BCD, 2);
	Shape ABF_m = make_mesh(ABF, 2);
	Shape BDF_m = make_mesh(BDF, 2);

	Shape AEC_m = make_mesh(AEC, 2);
	Shape ECD_m = make_mesh(ECD, 2);
	Shape AEF_m = make_mesh(AEF, 2);
	Shape EDF_m = make_mesh(EDF, 2);

	Shape up = merge_shapes(merge_shapes(ABC_m, BCD_m), merge_shapes(ABF_m, BDF_m));
	Shape down = merge_shapes(merge_shapes(AEC_m, ECD_m), merge_shapes(AEF_m, EDF_m));

	Shape s = merge_shapes(up, down);*/
	//s.write_to_file();

	cout<<"we are good to go"<<endl;


	GLFWwindow* window = glfwCreateWindow(width, height, "test", NULL, NULL);
	if( !window )
	{
		glfwTerminate();
		return -1;
	}

	glfwMakeContextCurrent(window);
	glfwSwapInterval(1);

	// keyboard handles
	glfwSetKeyCallback(window, key_callback);

	//the main loop
	render();
	glfwSwapBuffers(window);

	while( !glfwWindowShouldClose( window) )
	{
		glfwPollEvents();
		key_handler();
		render();
		draw_space();
		s.draw_shape();

		if(z_key_flag)
		{

			cout<<"Z is pressed"<<endl;
			s.write_to_file();
			cout<<"The file is ready!"<<endl;
		}

		glfwSwapBuffers(window);
	}

	glfwTerminate();

	return 0;
}






// render camera and scene
void render()
{
    glEnable(GL_DEPTH_TEST);
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
    glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );

    // reshape display window function, locating camera
    glMatrixMode(GL_PROJECTION); //these apply directy to my view, what I see in the clip!
    glLoadIdentity();
    gluPerspective(zoom, scr_width / scr_height, 1.0, 1000.0);


    // Use global parameters to rotate camera;
    glMatrixMode(GL_MODELVIEW); // this is all about the camera
    glLoadIdentity();
    gluLookAt(0, 0, 195.0, 0.0, 50.0, 0.0, 0.0, 1.0, 0.0);

    glTranslatef(-x_camera, -y_camera, -z_camera);
    glRotatef(x_rotate_angle, 1., 0., 0.);
    glRotatef(y_rotate_angle, 0., 1., 0.);

    // drawing the floor
    glColor3ub(255, 255, 255);

    // draw the bones
    glPushMatrix();
    glColor3ub(0, 255, 255);
    glPopMatrix();
}

// rendering objects at each time
/*void render_particles(vector<Particle>& particles)
{
	long unsigned int N = particles.size();
	for(long unsigned int i=0; i<N; i++)
	{
		vec F_tot(3, 1, fill::zeros);

		for(long unsigned int j=0; j<N; j++)
		{
			if( i!=j )
			{
				vec F = particles[j].pos - particles[i].pos;
				//cout<<"The force direction is: "<<i<<" "<<F<<endl;
				double len = vec_length(F);
				if(len < tresh_hold)
				{
					//cout<<"warning! we are in the dark age!"<<endl;
					continue;
				}
				F = F*(1./len);
				//cout<<"##"<<i<<len<<" "<<F<<endl;
				F = F*(G * particles[i].m * particles[j].m / pow(len, 2));
				F_tot = F_tot + F;
			}
		}
		//cout<<particles[i].F<<"|"<<F_tot<<endl;
		particles[i].F = F_tot;

	}
	/*particles[1].pos.x = 0;
	particles[1].pos.y = 50;
	particles[1].pos.y = 50;*/
	/*for(long unsigned int i=0; i<N; i++)
	{
		//Vec3f a =
	//	particles[i].pos = particles[i].pos + dl_t * particles[i].v;
		particles[i].v = particles[i].v + particles[i].F * ((1./particles[i].m) * dl_t);
		particles[i].pos = particles[i].pos + particles[i].v * (dl_t);
		//cout<<"pos: "<<i<<";"<<particles[i].pos<<endl;
		//cout<<"vel: "<<i<<";"<<particles[i].v<<endl;
	}

}*/
// --------------------------------- Key handling --------------------------------
void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
	if ( key == GLFW_KEY_Z )
    {
        if ( action == GLFW_PRESS )
            z_key_flag = true;
        else if ( action == GLFW_RELEASE )
            z_key_flag = false;
        else if ( action == GLFW_REPEAT)
            z_key_flag = true;
	}

    if ( key == GLFW_KEY_N )
    {
        if ( action == GLFW_PRESS )
            n_key_flag = true;
        else if ( action == GLFW_RELEASE )
            n_key_flag = false;
        else if ( action == GLFW_REPEAT)
            n_key_flag = true;
	}

	if ( key == GLFW_KEY_M )
    {
        if ( action == GLFW_PRESS )
            m_key_flag = true;
        else if ( action == GLFW_RELEASE )
            m_key_flag = false;
        else if ( action == GLFW_REPEAT)
            m_key_flag = true;
	}

	if ( key == GLFW_KEY_W )
    {
        if ( action == GLFW_PRESS )
            w_key_flag = true;
        else if ( action == GLFW_RELEASE )
            w_key_flag = false;
        else if ( action == GLFW_REPEAT)
            w_key_flag = true;
	}

	if ( key == GLFW_KEY_S )
    {
        if ( action == GLFW_PRESS )
            s_key_flag = true;
        else if ( action == GLFW_RELEASE )
            s_key_flag = false;
        else if ( action == GLFW_REPEAT)
            s_key_flag = true;
	}

	if ( key == GLFW_KEY_A )
    {
        if ( action == GLFW_PRESS )
            a_key_flag = true;
        else if ( action == GLFW_RELEASE )
            a_key_flag = false;
        else if ( action == GLFW_REPEAT)
            a_key_flag = true;
	}

	if ( key == GLFW_KEY_D )
    {
        if ( action == GLFW_PRESS )
            d_key_flag = true;
        else if ( action == GLFW_RELEASE )
            d_key_flag = false;
        else if ( action == GLFW_REPEAT)
            d_key_flag = true;
	}

	if ( key == GLFW_KEY_UP )
    {
        if ( action == GLFW_PRESS )
            UP_key_flag = true;
        else if ( action == GLFW_RELEASE )
            UP_key_flag = false;
        else if ( action == GLFW_REPEAT)
            UP_key_flag = true;
	}

	if ( key == GLFW_KEY_DOWN )
    {
        if ( action == GLFW_PRESS )
            DOWN_key_flag = true;
        else if ( action == GLFW_RELEASE )
            DOWN_key_flag = false;
        else if ( action == GLFW_REPEAT)
            DOWN_key_flag = true;
	}

	if ( key == GLFW_KEY_RIGHT )
    {
        if ( action == GLFW_PRESS )
            RIGHT_key_flag = true;
        else if ( action == GLFW_RELEASE )
            RIGHT_key_flag = false;
        else if ( action == GLFW_REPEAT)
            RIGHT_key_flag = true;
	}

	if ( key == GLFW_KEY_LEFT )
    {
        if ( action == GLFW_PRESS )
            LEFT_key_flag = true;
        else if ( action == GLFW_RELEASE )
            LEFT_key_flag = false;
        else if ( action == GLFW_REPEAT)
            LEFT_key_flag = true;
	}
}

void key_handler()
{
	if(n_key_flag)
		z_camera++;

	if(m_key_flag)
		z_camera--;


	if(w_key_flag)
		y_camera++;

	if(s_key_flag)
		y_camera--;

	if(a_key_flag)
		x_camera--;

	if(d_key_flag)
		x_camera++;

	if(UP_key_flag)
		x_rotate_angle++;

	if(DOWN_key_flag)
		x_rotate_angle--;

	if(LEFT_key_flag)
		y_rotate_angle++;

	if(RIGHT_key_flag)
		y_rotate_angle--;
}

// --------------------------------- drawing objects --------------------------------

void draw_space()
{
    glBegin(GL_LINES);
    // x
    glColor3ub(255, 0, 0);
    glVertex3f(0, 0, 0);
    glVertex3f(500, 0, 0);
    // y
    glColor3ub(0, 255, 0);
    glVertex3f(0, 0, 0);
    glVertex3f(0, 500, 0);
    // z
    glColor3ub(0, 0, 255);
    glVertex3f(0, 0, 0);
    glVertex3f(0, 0, 500);
    glEnd();
}
// ----------------------------------------


Shape make_triangle(Vec3f A, Vec3f B, Vec3f C)
{
	Segment AB(A, B);
	Segment BC(B, C);
	Segment CA(C, A);


	Face ABC(A, B, C);



	Shape s;
	s.points.insert(A);
	s.points.insert(B);
	s.points.insert(C);


	s.segments.insert(AB);
	s.segments.insert(BC);
	s.segments.insert(CA);

	s.faces.insert(ABC);

	return s;
}
// ----------------------------------------
Shape make_mesh(Shape start, int n)
{
    Shape res = make_mesh_handler(start, 0, n);
    return res;
}

Shape make_mesh_handler(Shape start, int itr, int n)
{
  if(itr == n)
  {
    return start;
  }
  //cout<<itr<<endl;
  set<Vec3f>::iterator set_itr = start.points.begin();
  Vec3f A((*set_itr).x, (*set_itr).y, (*set_itr).z);
  set_itr++;
  Vec3f B((*set_itr).x, (*set_itr).y, (*set_itr).z);
  set_itr++;
  Vec3f C((*set_itr).x, (*set_itr).y, (*set_itr).z);

  Vec3f m_AB = (A + B).num_multi(0.5);
  Vec3f m_BC = (B + C).num_multi(0.5);
  Vec3f m_CA = (C + A).num_multi(0.5);
  //cout<<C<<A<<C + A<<endl;
  Segment m_AB_A(m_AB, A);
  Segment m_AB_B(m_AB, B);
  Segment m_BC_B(m_BC, B);
  Segment m_BC_C(m_BC, C);
  Segment m_CA_C(m_CA, C);
  Segment m_CA_A(m_CA, A);
  Segment m_AB_m_CA(m_AB, m_CA);
  Segment m_AB_m_BC(m_AB, m_BC);
  Segment m_BC_m_CA(m_BC, m_CA);

  Face m_AB_m_CA_A(m_AB, m_CA, A);
  Face m_AB_m_BC_m_CA(m_AB, m_CA, m_BC);
  Face m_AB_m_BC_B(m_AB, m_BC, B);
  Face m_BC_m_CA_C(m_BC, m_CA, C);

  Shape up;
  up.points.insert(A);
  up.points.insert(m_AB);
  up.points.insert(m_CA);

  up.segments.insert(m_AB_A);
  up.segments.insert(m_CA_A);
  up.segments.insert(m_AB_m_CA);

  up.faces.insert(m_AB_m_CA_A);
  // ----------------------------

  Shape left;
  left.points.insert(B);
  left.points.insert(m_AB);
  left.points.insert(m_BC);

  left.segments.insert(m_AB_B);
  left.segments.insert(m_BC_B);
  left.segments.insert(m_AB_m_CA);

  left.faces.insert(m_AB_m_BC_B);
  // ----------------------------

  Shape right;
  right.points.insert(C);
  right.points.insert(m_BC);
  right.points.insert(m_CA);

  right.segments.insert(m_BC_C);
  right.segments.insert(m_CA_C);
  right.segments.insert(m_BC_m_CA);

  right.faces.insert(m_BC_m_CA_C);
  // ----------------------------

  Shape middle;
  middle.points.insert(m_AB);
  middle.points.insert(m_BC);
  middle.points.insert(m_CA);

  middle.segments.insert(m_AB_m_CA);
  middle.segments.insert(m_BC_m_CA);
  middle.segments.insert(m_AB_m_BC);

  middle.faces.insert(m_AB_m_BC_m_CA);
  //----
  Shape res_up = make_mesh_handler(up, itr+1, n);
  Shape res_middle = make_mesh_handler(middle, itr+1, n);
  Shape res_right = make_mesh_handler(right, itr+1, n);
  Shape res_left = make_mesh_handler(left, itr+1, n);
  Shape res = merge_shapes(merge_shapes(res_left, res_right), merge_shapes(res_up, res_middle));


  return res;
}

Shape merge_shapes(Shape s, Shape v)
{
  Shape m;
  m.points = s.points;
  m.segments = s.segments;
  m.faces = s.faces;
  for(set<Vec3f>::iterator itr=v.points.begin(); itr != v.points.end(); itr++)
    m.points.insert((*itr));

  for(set<Segment>::iterator itr=v.segments.begin(); itr != v.segments.end(); itr++)
    m.segments.insert((*itr));

  for(set<Face>::iterator itr=v.faces.begin(); itr !=v.faces.end(); itr++)
    m.faces.insert((*itr));

  return m;
}
