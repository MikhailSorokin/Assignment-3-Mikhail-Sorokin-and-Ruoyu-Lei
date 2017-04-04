#ifndef __MESH_HPP__

#include <QtGui>
#include <QtOpenGL>

#include <iostream>
#include <map>
#include <unordered_map>
using namespace std;

struct Mesh_Face {
    Mesh_Face() {
        vert[0] = vert[1] = vert[2] = -1;
    }

    Mesh_Face(long v0, long v1, long v2) {
        vert[0] = v0; vert[1] = v1; vert[2] = v2;
    }
    long vert[3]; // indices (in the vertex array) of all vertices (mesh_vertex)
    QVector3D faceNormal;

    bool operator==(const Mesh_Face &other) const
        {
            return vert[0] == other.vert[0] &&
            vert[1] == other.vert[1] &&
            vert[2] == other.vert[2] ;
        }
};


struct Mesh_Edge {

    // adding this default constructor to make the nested vector initialization work
    Mesh_Edge() {
        startVertexID = -1;
        endVertexID = -1;
    }

    Mesh_Edge(int start, int end) {
        startVertexID = start;
        endVertexID = end;
    }

    bool operator<(const Mesh_Edge &other) const
        {
        if (startVertexID < other.startVertexID)
        {
            return true;
        } else if (startVertexID == other.startVertexID)
        {
            if (endVertexID < other.endVertexID)
            {
                return true;
            }
            else if (endVertexID == other.endVertexID)
            {
                return false;
            }
        }
            return false;
    }

    bool operator==(const Mesh_Edge &other) const
      { if (startVertexID == other.startVertexID && endVertexID == other.endVertexID) {
        return true;
      } else if (startVertexID == other.endVertexID && endVertexID == other.startVertexID) {
        return true;
      } else {
        return false;
      }
      }

    int startVertexID;
    int endVertexID;
};


struct Mesh_Vertex {

    Mesh_Vertex(float x, float y, float z) {
        position = QVector3D(x,y,z);
    }

    QVector3D position;
    float avgEdgeLength;
    QVector3D normal;

    vector<Mesh_Edge> edges;
};

struct Mesh {
    vector<Mesh_Face> faces; // Mesh faces.
    vector<Mesh_Vertex> vertices; //Mesh Vertices
    vector<vector<Mesh_Face>> facesAdjVertex; //Faces connected to a vertex

    QOpenGLBuffer vertexBuffer, baryBuffer;



    bool load_obj(QString filename);
    void storeVBO();
    void recenter();
    void add_face(const vector<int> &cur_vert);
    void process_example();

    //Helper functions
    void debug_print(string str);
    float length(QVector3D edgeVector);
    float gaussian(float x, float y, float z, float u, float v, float w, float sigma);
    Mesh_Vertex* get_midpoint_vertex(QVector3D a, QVector3D b);
    int check_and_add(map<Mesh_Edge,int>& edgeToMidpointMap, vector<Mesh_Edge>& newEdges, Mesh_Edge midpointEdge, Mesh_Vertex* midpoint);
    vector<int> intersected_face(vector<Mesh_Face> f1,vector<Mesh_Face> f2, int v1 ,int v2);

    //MUST IMPLEMENT
    void compute_average_edge_lengths();
    void compute_vertex_normals();

    void inflate(float factor);
    void random_noise(float factor);
    void smooth();
    void sharpen(float factor);
    void split_faces();
    void split_long_edges();
    void loop_subdivision();
};

#endif // __MESH_HPP__
