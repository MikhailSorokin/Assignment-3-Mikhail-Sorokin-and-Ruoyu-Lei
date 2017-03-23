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

struct Mesh_Vertex {

    Mesh_Vertex(int x, int y, int z) {
        position = QVector3D(x,y,z);
    }
    QVector3D position;
    QVector3D normal;
};

/*struct Mesh_Edge {

    Mesh_Edge(int start, int end) {
        startVertexID = start;
        endVertexID = end;
    }
    int startVertexID;
    int endVertexID;
};*/

struct Mesh {
    vector<Mesh_Face> faces; // Mesh faces.
    vector<Mesh_Vertex> vertices; //Mesh Vertices
    //vector<Mesh_Edge> edges; //Mesh edges
    vector<vector<Mesh_Face>> facesAdjVertex; //Faces connected to a vertex

    QOpenGLBuffer vertexBuffer, baryBuffer;



    bool load_obj(QString filename);
    void storeVBO();
    void recenter();
    void add_face(const vector<int> &cur_vert);
    void process_example();
};

#endif // __MESH_HPP__
