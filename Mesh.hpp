#ifndef __MESH_HPP__

#include <QtGui>
#include <QtOpenGL>

#include <iostream>
#include <map>
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
};


struct Mesh_Edge {

    Mesh_Edge(int start, int end) {
        startVertexID = start;
        endVertexID = end;
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
    float length(QVector3D edgeVector);
    float gaussian(float x, float y, float z, float u, float v, float w, float sigma);

    //MUST IMPLEMENT
    void compute_average_edge_lengths();
    void compute_vertex_normals();

    void inflate(float factor);
    void random_noise(float factor);
    void smooth();
    void sharpen();
};

#endif // __MESH_HPP__
