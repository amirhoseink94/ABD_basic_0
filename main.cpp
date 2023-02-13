// openGL libraries
#include <GL/gl.h>
#include <GLFW/glfw3.h>
#include <GL/glut.h>
// system and c++ libraries
#include <stdio.h>
#include <iostream>
#include <cmath>
#include <time.h>
#include <vector>
// my libraries
#include "include/Particle.h"
#include "include/Vector.h"
#include "include/Pair.h"
#include "include/Shape.h"
#include "include/Utilities.h"
//#include "utilities.cpp"

// -----------------------------------------------------------------------------------
typedef Pair<Vec3D> Segment;
// --------------------------------- global variables --------------------------------

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
float dl_t = 10.;
float G = 0.0001;
float tresh_hold = 5;
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


using namespace std;
// rednders and drawing
void render();
void render_particles(vector<Particle>&);
void draw_cube();
void draw_space();
void draw_particles(vector<Particle>);

// tool functions
float distance(Vec3D, Vec3D);

//Shape make_mesh(Shape start, int itr, int n);


int main( int argc, char** argv)
{
	srand(time(NULL));
	srand(rand());
	vector<Particle> particles;

	Particle p1(Vec3D(50, 50, 50), Vec3D(0, 0, 0), 10, Vec3D(0, 0, 0) );
	Particle p2(Vec3D(0, 50, 50), Vec3D(0, 0, 0), 100, Vec3D(0, 0, 0) );

	particles.push_back(p1);
	particles.push_back(p2);
	cout<<p1.pos<<endl;
	cout<<p2.pos<<endl;
	Shape a;

	Vec3D A(10,10,10);
	Vec3D B(10,20,10);
	Vec3D C(10,10,20);
	Segment AB(A, B);
	Segment BC(B, C);
	Segment CA(C, A);
	Face ABC(A, B, C);
	Segment s(Vec3D(50, 50, 50),Vec3D(0, 0, 0));

	a.points.insert(A);
	a.points.insert(B);
	a.points.insert(C);

	a.segments.insert(AB);
	a.segments.insert(BC);
	a.segments.insert(CA);

	a.faces.insert(ABC);

	set<Vec3D>::iterator set_itr = a.points.begin();
  Vec3D AA((*set_itr).x, (*set_itr).y, (*set_itr).z);
  set_itr++;
  Vec3D BB((*set_itr).x, (*set_itr).y, (*set_itr).z);
  set_itr++;
  Vec3D CC((*set_itr).x, (*set_itr).y, (*set_itr).z);
  cout<<"A"<<AA<<endl;
  cout<<"B"<<BB<<endl;
  cout<<"C"<<CC<<endl;
	cout<<a.points.size()<<(A==B)<<endl;

	Util::make_mesh(a, 1);




	if ( !glfwInit() )
		return -1;
	cout << "(the glfw just initilized!)" << endl;
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

	//render();
	while( !glfwWindowShouldClose( window) )
	{
		glfwPollEvents();
		key_handler();
		render();
		draw_space();
		draw_cube();
		a.draw_shape();
		render_particles(particles);
		//cout<<"z is pressed"<<endl;
		/*if(z_key_flag)
		{
		}*/
		//render_particles(particles);
		draw_particles(particles);
		glfwSwapBuffers(window);
	}

	glfwTerminate();

	return 0;
}
// ---
Shape make_mesh(Shape start, int itr, int n)
{
  if(itr == n)
  {
    return start;
  }
  cout<<"we are here!"<<endl;
	return start;
  /*set<Vec3D>::iterator set_itr = start.points.begin();
  Vec3D A((*set_itr).x, (*set_itr).y, (*set_itr).z);
  itr++;
  Vec3D B((*set_itr).x, (*set_itr).y, (*set_itr).z);
  itr++;
  Vec3D C((*set_itr).x, (*set_itr).y, (*set_itr).z);
  cout<<"A"<<A<<endl;
  cout<<"B"<<B<<endl;
  cout<<"C"<<C<<endl;
*/
}


// --------------------------------- rendering --------------------------------

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
void render_particles(vector<Particle>& particles)
{
	long unsigned int N = particles.size();
	for(long unsigned int i=0; i<N; i++)
	{
		Vec3D F_tot(0, 0, 0);
		for(long unsigned int j=0; j<N; j++)
		{
			if( i!=j )
			{
				Vec3D F = particles[j].pos - particles[i].pos;
				cout<<"#"<<i<<F<<endl;
				float len = F.length();
				if(len < tresh_hold)
				{
					cout<<"warning! we are in the dark age!"<<endl;
					continue;
				}
				F = F.num_multi(1./len);
				cout<<"##"<<i<<len<<" "<<F<<endl;
				F = F.num_multi(G * particles[j].m / pow(len, 2));
				F_tot = F_tot + F;
			}
		}
		cout<<particles[i].F<<"|"<<F_tot<<endl;
		int x;
		//cin>>x;
		particles[i].F = F_tot;

	}
	/*particles[1].pos.x = 0;
	particles[1].pos.y = 50;
	particles[1].pos.y = 50;*/
	for(long unsigned int i=0; i<particles.size(); i++)
	{
		//Vec3D a =
	//	particles[i].pos = particles[i].pos + dl_t * particles[i].v;
		particles[i].v = particles[i].v + particles[i].F.num_multi(dl_t);
		particles[i].pos = particles[i].pos + particles[i].v.num_multi(dl_t);
		cout<<"pos: "<<i<<";"<<particles[i].pos<<endl;
		cout<<"vel: "<<i<<";"<<particles[i].v<<endl;
	}

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

void draw_particles(vector<Particle> p)
{
	for(long unsigned int i=0; i<p.size(); i++)
	{
		glPointSize(2);
		glBegin(GL_POINTS);
			if(i%2==0)
				glColor3ub(255, 255, 0);
			else
				glColor3ub(0, 255, 255);
			glVertex3f(p[i].pos.x, p[i].pos.y, p[i].pos.z);
		glEnd();
	}
}


void draw_cube()
{
	float r = 3;
	for(int i=0; i<=180; i++)
	{
		for(int j=0; j<360; j++)
		{
			glBegin(GL_POINTS);
				glColor3ub(255, i, j);
				glVertex3f(r * sin(i) * cos(j), r * sin(i) * sin(j), r * cos(i));
			glEnd();
		}
	}
}

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
