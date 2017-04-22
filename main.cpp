// Victor Priano
// Computer Graphics (Loop Subdivision)

////////////////////////////////////////////////////////////////////////
// A simple wrapper for to store 3D vectors
// This file is not my original work but was supplied to us.
// Functionality: This struct allows me to store the coordinates
//                of any vertex specified on any subdivision level
////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <cmath>
using namespace std;

struct Vector3
{
        float x,y,z;

        Vector3() : x(0.0), y(0.0), z(0.0)
        {}

        Vector3(float x, float y, float z)
                : x(x), y(y), z(z)
        {}

        Vector3(const Vector3 & v)
                : x(v.x), y(v.y), z(v.z)
        {}
        
        Vector3 operator+(const Vector3 & rhs) const
        { return Vector3(x + rhs.x, y + rhs.y, z + rhs.z); }
        Vector3 operator-(const Vector3 & rhs) const
        { return Vector3(x - rhs.x, y - rhs.y, z - rhs.z); }
        Vector3 operator*(float rhs) const
        { return Vector3(x * rhs, y * rhs, z * rhs); }
        Vector3 operator/(float rhs) const
        { return Vector3(x / rhs, y / rhs, z / rhs); }
        Vector3 operator+=(const Vector3 & rhs)
        { x += rhs.x; y += rhs.y; z += rhs.z; return *this; }
        Vector3 operator-=(const Vector3 & rhs)
        { x -= rhs.x; y -= rhs.y; z -= rhs.z; return *this; }
        Vector3 operator*=(float rhs)
        { x *= rhs; y *= rhs; z *= rhs; return *this; }
        Vector3 operator/=(float rhs)
        { x /= rhs; y /= rhs; z /= rhs; return *this; }

        float magnitude() const
        { return sqrt(x * x + y * y + z * z); }
        void normalize()
        { *this /= magnitude(); }
        float dot(const Vector3 & rhs) const
        {
                return x * rhs.x + y * rhs.y + z * rhs.z;
        }
        Vector3 cross(const Vector3 & rhs) const
        {
                return Vector3(y * rhs.z - z * rhs.y,
                                        z * rhs.x - x * rhs.z,
                                        x * rhs.y - y * rhs.x);
        }
        
        void print() 
        { cout << x << " " << y << " " << z << endl, cout.flush(); }
};
////////////////////////////////////////////////////////////////////////
//                        subdivision main                            //
////////////////////////////////////////////////////////////////////////
#include <GL/glut.h>
#include <cmath>
#include <fstream>
#include <vector>
#include <limits>

// define constants to help with calculations
const int WINDOW_WIDTH = 800;
const int WINDOW_HEIGHT = 800;
#define WEIGHT_1 3.0/8.0
#define WEIGHT_2 1.0/8.0
#define PI 3.14159265

// variables dealing directly with the triangle mesh data
// vector structs to hold the coordinates and connectivity amongst them
vector< Vector3 > coordinates;
vector< Vector3 > vertices;
vector< vector<int> > connect;
vector< Vector3 > new_coordinates;
vector< Vector3 > new_vertices;
vector< Vector3 > tmp_coordinates;
vector< Vector3 > tmp_vertices;
vector< Vector3 > update_coordinates;
vector< Vector3 > update_vertices;
int num_coordinates = 0;
int faces = 0;

// since OpenGl runs an infinite loop, I need flags to activate 
// certain functions
int view=0;
int divide=4;
int test_print=0;
int once=0;
int old_vertices=0;
int subdivision_level=1;

////////////////////////////////////////////////////////////////////////
// A storage container for an edge that specifies what vertex is
// adjacent and opposite of the edge
// Functionality: easier access to vertices that give weight to a
//                certain edge
////////////////////////////////////////////////////////////////////////
struct edges 
{
        float opp1;
        float opp2;
        float adj1;
        float adj2;
};
vector <edges> edgez;
////////////////////////////////////////////////////////////////////////
// Renders a quad at cell (x, y) with dimensions CELL_LENGTH          //
////////////////////////////////////////////////////////////////////////
void renderPixel(float x, float y, float z, float r = 0.0, float g = 0.0, float b = 0.0)
{
        // ...
        // Complete this function
        // ...
        glBegin(GL_POINTS);
        glColor3f(r,g,b);
        glVertex3f(x,y,0.0);        
        glEnd();
}
////////////////////////////////////////////////////////////////////////
//                      debugging                                     //
////////////////////////////////////////////////////////////////////////
// Simple debugging helper functions
// Functionality: Renders a point at (x,y) in yellow to
//                distinguish the newly inserted coordinate
void test_point(float x, float y)
{
        glPointSize(4.0);
        glColor3f(1.0, 1.0, 0.0);
        glEnable(GL_POINT_SMOOTH);
        renderPixel(x, y, 0.0);
        glPointSize(1.0);
        glColor3f(1.0, 1.0, 1.0);
        glDisable(GL_POINT_SMOOTH);
}
// Functionality: Prints to the terminal the values of the vectors
//                containing the coordinates, face list, and adjcency
//                                  list, mostly used for debugging purposes
void print()
{
        cout << "VERTICES\n", cout.flush();
        for (int i=0; i<new_vertices.size();++i)
                cout << new_vertices[i].x << " " 
                     << new_vertices[i].y << " " 
                     << new_vertices[i].z << endl, cout.flush();
        cout << "COORDINATES\n", cout.flush();
        for (int j=0; j<new_coordinates.size();++j)
                cout << "( " << new_coordinates[j].x << " " 
                     << new_coordinates[j].y << " " << new_coordinates[j].z 
                     << " )" << endl, cout.flush();
        cout << "ADJACENCY LIST\n", cout.flush();
        for (int k=0; k<connect.size();++k)
        {        
                cout << connect[k][0] << ": ", cout.flush();
                for (int l=1; l<connect[k].size();++l)
                        cout << connect[k][l] << " ", cout.flush();
                cout << endl, cout.flush();
        }
}
////////////////////////////////////////////////////////////////////////
//                    Graphics rendering                              //
////////////////////////////////////////////////////////////////////////
// Rendering algorithms (completed in lab)
// Functionality: DDA renders a line from (x1, y1) to (x2, y2)
void DDA (float x1, float y1, float x2, float y2)
{
        float dx = x2-x1;
        float dy = y2-y1;
        float m = dy/dx;
        float steps = (abs(dx) > abs(dy)) ? abs(dx) : abs(dy);                        
        float x=x1, y=y1;

        float xinc = dx/steps;
        float yinc = dy/steps;
        
        for (int i = 0; i < (int)steps; ++i)
        {
                x+=xinc;
                y+=yinc;
                renderPixel(x, y, 0);
        }
}
// Functionality: Renders part of a circle using an initial start (x,y)
//                                   with radius r and moving along the arc calculated
//                              through a formula
//                NOT USED
void MidpointCircle(int x, int y, int r)
{
        int i=0, j=r;
        float midpoint = 0;
        while (i <= j)
        {
                renderPixel(x+i, y+j, 0);
                renderPixel(x-i, y+j, 0);
                renderPixel(x+i, y-j, 0);
                renderPixel(x-i, y-j, 0);
                renderPixel(x+j, y+i, 0);
                renderPixel(x-j, y+i, 0);
                renderPixel(x+j, y-i, 0);
                renderPixel(x-j, y-i, 0);
                ++i;

                midpoint = (r*r) - (i*i) - ((j+0.5)*(j+0.5));
                j = (midpoint < 0) ? j-1 : j;
        }
}
////////////////////////////////////////////////////////////////////////
//                          read data                                 //
////////////////////////////////////////////////////////////////////////
// Reads a text file through command line
// Functionality: Parses the specified file and stores the files 
//                contents in the appropriate container
void read_file(char* filename)
{
        float data = 0.0;
        ifstream file;
        file.open(filename);
        if (file == NULL) 
        { 
                cout << "read_file() failed: FILE DOES NOT EXIST.\n";
                exit(1); 
        }
        file >> data;
        num_coordinates = (int)data;
        file >> data;
        faces = (int)data;
        Vector3 point;
        for (int i=0; i < num_coordinates; ++i)
        {
                file >> data;
                point.x = data;
                file >> data;
                point.y = data;
                file >> data;
                point.z = data;
                coordinates.push_back(point);
        }
        Vector3 edge;
        for (int i=0; i < faces; ++i)
        {
                file >> data;
                edge.x = data;
                file >> data;
                edge.y = data;
                file >> data;
                edge.z = data;
                vertices.push_back(edge);
        }
}
////////////////////////////////////////////////////////////////////////
//                       helper functions                             //
////////////////////////////////////////////////////////////////////////
// Simple helper functions relating to the subdividing functions
// The applied weight to an existing vertex
//Function
float alpha(int n)
{
        float a = 1.0/4.0;
        float b = cos(2*PI/n)*a;
        b += WEIGHT_1;
        b *= b;
        float c = 5.0/8.0;
        c -= b;
        return (1.0/float(n))*c;
}

float beta(int n)
{
        return 1.0 - (float)n*alpha(n);
}
////////////////////////////////////////////////////////////////////////
//                     helper functions                               //
////////////////////////////////////////////////////////////////////////
bool is_contained_int(vector<int> v, float n)
{
        for (int i=0; i< v.size(); ++i) if ((float)v[i] == n) return true;
        return false;
}

bool is_contained_Vector3(vector<Vector3> v, Vector3 t)
{
        for (int i=0; i< v.size(); ++i) if (v[i].x == t.x && v[i].y == t.y && v[i].z == t.z) return true;
        return false;
}

int find(vector<Vector3>v, Vector3 t)
{
        for (int i=0; i< v.size(); ++i) if (v[i].x == t.x && v[i].y == t.y && v[i].z == t.z) return i;
        return -1;
}
////////////////////////////////////////////////////////////////////////
//                      main subdivide functions                      //
////////////////////////////////////////////////////////////////////////
// Creates the adjacency list
// Functionality: Updates the connect vector which stores a vertex 
//                  index and the indices of all vertices that connect to
//                that specific vertex
void connectivity(vector<Vector3> v, vector<Vector3> c)
{
        vector<int> tmp;
        connect.clear();
        for(int i=0; i < c.size(); ++i) tmp.push_back(i), connect.push_back(tmp), tmp.clear();
        for (int j=0; j < c.size(); ++j)
        {
                for (int k=0; k < v.size(); ++k)
                {
                        if (v[k].x == connect[j][0])
                        {
                                if (!is_contained_int(connect[j], v[k].y)) connect[j].push_back((int)v[k].y);
                                if (!is_contained_int(connect[j], v[k].z)) connect[j].push_back((int)v[k].z);
                        }
                        if (v[k].y == connect[j][0])
                        {
                                if (!is_contained_int(connect[j], v[k].x)) connect[j].push_back((int)v[k].x);
                                if (!is_contained_int(connect[j], v[k].z)) connect[j].push_back((int)v[k].z);
                        }
                        if (v[k].z == connect[j][0])
                        {
                                if (!is_contained_int(connect[j], v[k].x)) connect[j].push_back((int)v[k].x);
                                if (!is_contained_int(connect[j], v[k].y)) connect[j].push_back((int)v[k].y);
                        }
                }
        }
}
// Functionality: applies weight to existing vertices using neighbors
void update_even(vector<Vector3> c)
{
        float alpha_constant=0.0;
        float beta_constant=0.0;
        for (int i=0; i < old_vertices; ++i)
        {
                Vector3        update = c[i];
                Vector3 tmp;
                alpha_constant = alpha(connect[i].size()-1);
                beta_constant = beta(connect[i].size()-1);
                update*= beta_constant;
                for (int j=1; j < connect[i].size(); ++j) { tmp+=c[connect[i][j]]; }
                tmp*=alpha_constant;
                update+=tmp;
                new_coordinates[i] = update;
        }
}
// Functionality: returns the midpoint that has been already weighed
Vector3 update_odd(vector<Vector3>c, float n, float x)
{
        Vector3 odd;
        for (int i=0; i<edgez.size();++i)
        {
                if ( (edgez[i].adj1==n && edgez[i].adj2==x) || (edgez[i].adj1==x && edgez[i].adj2==n) )
                {
                        odd.x=WEIGHT_1*c[(int)n].x+WEIGHT_1*c[(int)x].x
                              +WEIGHT_2*c[(int)edgez[i].opp1].x+WEIGHT_2*c[(int)edgez[i].opp2].x;
                        odd.y=WEIGHT_1*c[(int)n].y+WEIGHT_1*c[(int)x].y
                              +WEIGHT_2*c[(int)edgez[i].opp1].y+WEIGHT_2*c[(int)edgez[i].opp2].y;
                        odd.z=WEIGHT_1*c[(int)n].z+WEIGHT_1*c[(int)x].z
                              +WEIGHT_2*c[(int)edgez[i].opp1].z+WEIGHT_2*c[(int)edgez[i].opp2].z;
                        return odd;
                }
        }
        return odd;
}
// Functionality: groups edges into a premad edge struct that stores
//                each edges adjacent and opposite points
void group_edges(vector<Vector3> v)
{
        float side_1=0.0, side_2=0.0,side_3=0.0;
        float next1=0.0, next2=0.0, next3=0.0;
        edgez.clear();
        for (int i=0; i < v.size()-1; ++i)
        {
                side_1 = v[i].x;
                side_2 = v[i].y;
                side_3 = v[i].z;
                for (int j=i+1; j < v.size(); ++j)
                {
                        edges line;
                        next1 = v[j].x;
                        next2 = v[j].y;
                        next3 = v[j].z;
                        if ( (side_1 == next1 || side_1 == next2 || side_1 == next3) &&
                             (side_2 == next1 || side_2 == next2 || side_2 == next3) )
                        {
                                line.adj1 = side_1, line.adj2 = side_2, line.opp1 = side_3;
                                if (next1 != side_1 && next1 != side_2) line.opp2 = next1;
                                else if (next2 != side_1 && next2 != side_2) line.opp2 = next2;
                                else line.opp2 = next3;
                                edgez.push_back(line);
                      }
                        else if ( (side_1 == next1 || side_1 == next2 || side_1 == next3) && 
                                       (side_3 == next1 || side_3 == next2 || side_3 == next3))
                        {
                                line.adj1 = side_1, line.adj2 = side_3, line.opp1 = side_2;
                                if (next1 != side_1 && next1 != side_3) line.opp2 = next1;
                                else if (next2 != side_1 && next2 != side_3) line.opp2 = next2;
                                else line.opp2 = next3;
                                edgez.push_back(line);
                        }
                        else if ( (side_2 == next1 || side_2 == next2 || side_2 == next3) && 
                                      (side_3 == next1 || side_3 == next2 || side_3 == next3))
                        {
                                line.adj1 = side_2, line.adj2 = side_3, line.opp1 = side_1;
                                if (next1 != side_2 && next1 != side_3) line.opp2 = next1;
                                else if (next2 != side_2 && next2 != side_3) line.opp2 = next2;
                                else line.opp2 = next3;
                                edgez.push_back(line);
                        }
                        else {}
                }
        }
}
// Functionality: The meat of the loop subdivision code. It generates
//                the midpoint location and ask for its new position.
//                It updates the coordinates and vertices vectors
void subdivide(vector<Vector3> v, vector<Vector3> c)
{
        new_vertices.clear();
        new_coordinates.clear();
        new_coordinates=c;
        old_vertices=c.size();
        Vector3 update, mid1, mid2, mid3;
        for (int i=0; i<v.size();++i)
        {

                mid1 = update_odd(c, v[i].x, v[i].y);
                if (!is_contained_Vector3(new_coordinates, mid1)) new_coordinates.push_back(mid1);
                mid2 = update_odd(c, v[i].x, v[i].z);
                if (!is_contained_Vector3(new_coordinates, mid2)) new_coordinates.push_back(mid2);
                mid3 = update_odd(c, v[i].y, v[i].z);
                if (!is_contained_Vector3(new_coordinates, mid3)) new_coordinates.push_back(mid3);
                update.x=v[i].x;
                update.y=find(new_coordinates, mid1);
                update.z=find(new_coordinates, mid2);
                new_vertices.push_back(update);
                update.x=v[i].y;
                update.y=find(new_coordinates, mid1);
                update.z=find(new_coordinates, mid3);
                new_vertices.push_back(update);
                update.x=v[i].z;
                update.y=find(new_coordinates, mid2);
                update.z=find(new_coordinates, mid3);
                new_vertices.push_back(update);
                update.x=find(new_coordinates, mid1);
                update.y=find(new_coordinates, mid2);
                update.z=find(new_coordinates, mid3);
                new_vertices.push_back(update);
        }
}
////////////////////////////////////////////////////////////////////////
//                      print                                         //
////////////////////////////////////////////////////////////////////////
// Functionality:  prints whatever is in the coordinates vector based on
//                 the adjacency list
void print_mesh(vector<Vector3> c, int n)
{
        for (int i = 0; i < connect.size();++i)
                for (int j=1; j < connect[i].size(); ++j)
                {
                        if (!n) DDA(c[i].x*8, c[i].y*8, c[connect[i][j]].x*8, c[connect[i][j]].y*8);
                        else if (n==1) DDA(c[i].x*8, c[i].z*8,c[connect[i][j]].x*8, c[connect[i][j]].z*8);
                        else DDA(c[i].y*8, c[i].z*8, c[connect[i][j]].y*8, c[connect[i][j]].z*8);
                }
}
// basic GL_render() call
void GL_render()
{
        glClear(GL_COLOR_BUFFER_BIT);

        
        print_mesh(new_coordinates, view);
        glutPostRedisplay();
        
        glutSwapBuffers();
}
// Functionality: updates program based on keyboard input
void KeyboardInput(unsigned char key, int x, int y)
{
        switch (key)
        {
                //    T    //
                case 84:
                case 116:
                        subdivision_level=1;
                        ++view;
                        view = (view > 2) ? 0 : view;
                        new_coordinates.clear(), tmp_coordinates.clear();
                        new_vertices.clear(), tmp_vertices.clear();
                        new_coordinates=coordinates, tmp_coordinates=coordinates;
                        new_vertices=vertices, tmp_vertices=vertices;
                        connectivity(new_vertices, new_coordinates);
                        group_edges(new_vertices);
                        break;
                //    L    //
                case 76:
                case 108:
                        ++subdivision_level;
                        if (subdivision_level<5)
                        {
                                tmp_coordinates=new_coordinates;
                                tmp_vertices=new_vertices;
                                subdivide(tmp_vertices, tmp_coordinates);
                                update_even(tmp_coordinates);
                                connectivity(new_vertices, new_coordinates);                        
                                group_edges(new_vertices);
                        }
                        else
                        { 
                                subdivision_level=1;
                                new_coordinates=update_coordinates;
                                new_vertices=update_vertices; 
                                connectivity(new_vertices, new_coordinates);                        
                                group_edges(new_vertices);
                        }
                        break;
                default:
                        break;
        }
}
////////////////////////////////////////////////////////////////////////
//                      initialize original                           // 
////////////////////////////////////////////////////////////////////////
void subdivide_init()
{
        new_vertices=vertices, tmp_vertices=vertices, update_vertices=vertices;
        new_coordinates=coordinates, tmp_coordinates=coordinates, update_coordinates=coordinates;
        connectivity(vertices, coordinates);
        group_edges(vertices);
}
////////////////////////////////////////////////////////////////////////
//Initializes OpenGL attributes
void GLInit(int* argc, char** argv)
{
        glutInit(argc, argv);
        glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);
        glutInitWindowSize(WINDOW_WIDTH, WINDOW_HEIGHT);

        // ...
        // Complete this function
        // ...
        glutCreateWindow(“Victor Priano“);

        // The default view coordinates is (-1.0, -1.0) bottom left & (1.0, 1.0) top right.
        // For the purposes of this lab, this is set to the number of pixels
        // in each dimension.
        glMatrixMode(GL_PROJECTION_MATRIX);
        glOrtho(0, WINDOW_WIDTH, 0, WINDOW_HEIGHT, -1, 1);
        glClearColor(1.0f, 1.0f, 1.0f, 0.0f);

        glutDisplayFunc(GL_render);
        glutKeyboardFunc(KeyboardInput);
}
////////////////////////////////////////////////////////////////////////
//                       main code                                    //
////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv)
{        
        read_file(argv[1]);
        subdivide_init();
        GLInit(&argc, argv);
        glutMainLoop();

        return 0;
}
