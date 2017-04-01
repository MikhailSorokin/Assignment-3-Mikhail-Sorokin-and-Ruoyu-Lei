#include <iostream>
#include <fstream>
#include <iomanip>
#include <unordered_map>
#include <ctime>

#include "Mesh.hpp"

using namespace std;

void Mesh::debug_print(string str) {
    cout<<"debug info: "<<str<<"\n";
}

void Mesh::add_face(const vector<int> &cur_vert) {
    if(cur_vert.size() > 3) {
        // If number of edges in face is greater than 3,
        // decompose into triangles as a triangle fan.
        int v0 = cur_vert[0], v1 = cur_vert[1], v2 = cur_vert[2];

        faces.push_back(Mesh_Face(v0, v1, v2));     // First face

        // all subsequent faces
        for( size_t i = 3; i < cur_vert.size(); i++ ) {
            v1 = v2; v2 = cur_vert[i];
            faces.push_back(Mesh_Face(v0, v1, v2));
        }
    }
    else if(cur_vert.size() == 3) {
        faces.push_back(Mesh_Face(cur_vert[0], cur_vert[1], cur_vert[2]));
    }
}

bool Mesh::load_obj(QString filename) {
    QFile objfile(filename);
    if (!objfile.open(QIODevice::ReadOnly | QIODevice::Text)) {
        return false; //error
    }
    QTextStream in(&objfile);
    long face_cnt = 0;

    while (!in.atEnd()) {
        QString line = in.readLine();
        line = line.trimmed();
        line = line.replace("\t", " ");

        QStringList tokens = line.trimmed().split(' ', QString::SkipEmptyParts);
        if(tokens.size() == 0) continue;

        if(tokens[0] == "v") {
            if(tokens.size() < 4) return false; // eror
            float x = tokens[1].toFloat();
            float y = tokens[2].toFloat();
            float z = tokens[3].toFloat();
            vertices.push_back(Mesh_Vertex(x,y,z));
        }
        if(tokens[0] == "f") {
            vector<int> cur_vert;
            for(int i = 1; i < tokens.size(); i++) {
                QStringList indexes = tokens[i].split("/");
                if(indexes.size() >= 1) {
                    if(indexes[0].toLong() < 0) {  cur_vert.push_back(vertices.size() + indexes[0].toLong()); }
                    else { cur_vert.push_back(indexes[0].toLong() - 1); }
                }
            }
            face_cnt++;
            add_face(cur_vert);
        }
    }

    cout << "face_cnt=" << face_cnt << endl;
    cout << "faces.size()=" << faces.size() << endl;
    cout << "vertices.size()=" << vertices.size() << endl;

    recenter();
    return true;
}

void Mesh::compute_average_edge_lengths() {
    //Construct the edges vector
    for (Mesh_Face& face : faces) {
        for (int vi = 0; vi < 3; vi++) {
            //cout << "(" << vertices[face.vert[vi]].position[0] << "," << vertices[face.vert[vi]].position[1] << "," << vertices[face.vert[vi]].position[2] << ")";
            vertices[face.vert[vi]].edges.push_back(Mesh_Edge(face.vert[vi], face.vert[(vi + 1) % 3]));
            vertices[face.vert[vi]].edges.push_back(Mesh_Edge(face.vert[vi], face.vert[(vi + 2) % 3]));
            //cout << "Size: " << facesAdjVertex[face.vert[vi]].size() << endl;
            //cout << "(" << facesAdjVertex[face.vert[vi]][0].faceNormal[0] << "," << facesAdjVertex[face.vert[vi]][0].faceNormal[1] << "," << facesAdjVertex[face.vert[vi]][0].faceNormal[2] << ")";
        }
    }

    //Compute normals for every vertex
    for (int vi = 0; vi < vertices.size(); vi++) {
        float vertexAdjEdgeSum = 0.f;

        for (Mesh_Edge& edge : vertices[vi].edges) {
            vertexAdjEdgeSum += length(vertices[edge.endVertexID].position - vertices[edge.startVertexID].position);
        }

        vertices[vi].avgEdgeLength = vertexAdjEdgeSum / vertices[vi].edges.size();
        //cout << vertices[vi].avgEdgeLength << endl;
    }
}

float Mesh::length(QVector3D e) {
    return sqrt(e[0] * e[0] + e[1] * e[1] + e[2] * e[2]);
}

void Mesh::compute_vertex_normals() {

    //Compute normals for every face
    for (Mesh_Face& face : faces) {
        QVector3D P0P1 = vertices[face.vert[1]].position - vertices[face.vert[0]].position;
        QVector3D P0P2 = vertices[face.vert[2]].position - vertices[face.vert[0]].position;
        face.faceNormal = QVector3D::crossProduct(P0P1, P0P2);
    }

    /* Adding adjacent faces to every vertex */
    for (int i = 0; i < vertices.size(); i++) {
        facesAdjVertex.push_back(vector<Mesh_Face>());
    }

    for (Mesh_Face& face : faces) {
        for (int vi = 0; vi < 3; vi++) {
            //cout << "(" << vertices[face.vert[vi]].position[0] << "," << vertices[face.vert[vi]].position[1] << "," << vertices[face.vert[vi]].position[2] << ")";
            facesAdjVertex[face.vert[vi]].push_back(face);
            //cout << "Size: " << facesAdjVertex[face.vert[vi]].size() << endl;
            //cout << "(" << facesAdjVertex[face.vert[vi]][0].faceNormal[0] << "," << facesAdjVertex[face.vert[vi]][0].faceNormal[1] << "," << facesAdjVertex[face.vert[vi]][0].faceNormal[2] << ")";
        }
    }

    //Compute normals for every vertex
    for (int vi = 0; vi < vertices.size(); vi++) {
        QVector3D normalSum(0,0,0);

        vector<Mesh_Face> adjFaces = facesAdjVertex[vi];
        for (Mesh_Face& adjFace : adjFaces) {
            normalSum += adjFace.faceNormal;
        }

        vertices[vi].normal = normalSum.normalized();
        //cout << "(" << vertices[vi].normal[0] << "," << vertices[vi].normal[1] << "," << vertices[vi].normal[2] << ")" << endl;
    }
}


void Mesh::recenter() {
    if( vertices.size() < 1) return;
    QVector3D maxPoint = vertices[0].position;
    QVector3D minPoint = vertices[0].position;

    // Find the AABB
    for( uint i = 0; i < vertices.size(); ++i ) {
        QVector3D & point = vertices[i].position;
        if( point[0] > maxPoint[0] ) maxPoint[0] = point[0];
        if( point[1] > maxPoint[1] ) maxPoint[1] = point[1];
        if( point[2] > maxPoint[2] ) maxPoint[2] = point[2];
        if( point[0] < minPoint[0] ) minPoint[0] = point[0];
        if( point[1] < minPoint[1] ) minPoint[1] = point[1];
        if( point[2] < minPoint[2] ) minPoint[2] = point[2];
    }

    // Center of the AABB
    QVector3D center = QVector3D( (maxPoint[0] + minPoint[0]) / 2.0f,
    (maxPoint[1] + minPoint[1]) / 2.0f,
    (maxPoint[2] + minPoint[2]) / 2.0f );

    // Translate center of the AABB to the origin
    for( uint i = 0; i < vertices.size(); ++i ) {
        QVector3D & point = vertices[i].position;
        point = point - center;
    }
}


void Mesh::process_example() {
    for(size_t v = 0; v < vertices.size(); v++) {
        if(vertices[v].position[0] > 0) {
            vertices[v].position[0] += 3.5;
        }
    }
}

void Mesh::inflate(float factor) {
    //facesAdjVertex.clear();

    compute_average_edge_lengths();
    compute_vertex_normals();

    for(size_t v = 0; v < vertices.size(); v++) {
        float displacement = vertices[v].avgEdgeLength * factor;
        vertices[v].position = vertices[v].position + vertices[v].normal * displacement;
    }
}

void Mesh::random_noise(float factor) {
    compute_average_edge_lengths();

    for(size_t v = 0; v < vertices.size(); v++) {
        float displacement = vertices[v].avgEdgeLength * factor;
        float clampedDisplacement = (displacement >= 0) ? displacement : 0;
        srand(rand()%500);
        vertices[v].position = vertices[v].position + QVector3D(rand()%2,rand()%2,rand()%2) * clampedDisplacement;
    }
}


// 3d gaussian
float Mesh::gaussian(float x, float y, float z, float u, float v, float w, float sigma){
    float std = -0.5 * (x - u)*(x - u) / (sigma * sigma) +
            -0.5 * (y - v)*(y - v) / (sigma * sigma) +
            -0.5 * (z - w)*(z - w) / (sigma * sigma);
    float wght = exp(std);

    return wght;
}

void Mesh::smooth() {
    compute_average_edge_lengths();

    // the new vertices positions
    vector<QVector3D> smoothedVertices;

    for (size_t iv = 0; iv < vertices.size(); iv++){
        QVector3D currPst = vertices[iv].position;
        float sigma = vertices[iv].avgEdgeLength;
        float norm_const = 0;
        QVector3D point(0,0,0);

        for (size_t iu = 0; iu < vertices[iv].edges.size(); iu++){
            int neighborId = vertices[iv].edges[iu].startVertexID == iv ? vertices[iv].edges[iu].endVertexID : vertices[iv].edges[iu].startVertexID;

            QVector3D neibPst = vertices[neighborId].position;

            float gauss = gaussian(currPst.x(), currPst.y(), currPst.z(), neibPst.x(), neibPst.y(),neibPst.z(),sigma);
            norm_const += gauss;
            point += gauss * neibPst;
        }

        QVector3D pos = (point + currPst) / (norm_const + 1);
        // QVector3D pos = vertices[iv].position;
        smoothedVertices.push_back(pos);
    }

    // replace with new positions
    for (size_t i = 0; i < vertices.size();i++){
        vertices[i].position = smoothedVertices[i];
    }
}

void Mesh::sharpen(){
    // ???
}

int Mesh::check_and_add(map<Mesh_Edge*,Mesh_Vertex*>& edgeToMidpointMap, vector<Mesh_Edge>& newEdges, Mesh_Edge* midpointEdge, Mesh_Vertex* midpoint) {
    if (edgeToMidpointMap.find(midpointEdge) == edgeToMidpointMap.end()) {
        edgeToMidpointMap[midpointEdge] = midpoint;
        vertices.push_back(*midpoint);
        vertices[midpointEdge->startVertexID].edges.push_back(*midpointEdge);
        vertices[midpointEdge->endVertexID].edges.push_back(*midpointEdge);
        newEdges.push_back(*midpointEdge);

        cout << "AGH" << endl;

        return vertices.size() - 1;
    }

    return midpointEdge->endVertexID; //endVertexID will always be our new ID
}

Mesh_Vertex* Mesh::get_midpoint_vertex(QVector3D a, QVector3D b) {
    Mesh_Vertex* mid = new Mesh_Vertex((a.x() + b.x())/2, (a.y() + b.y())/2, (a.z() + b.z())/2);
    return mid;
}

void Mesh::add_edges(int start, int end, int index, string key, unordered_map<string,int>& h_edges, vector<vector<Mesh_Edge>>& n_edges) {
    // string key = start < end ? to_string(start) + "," + to_string(end) : to_string(end) + "," + to_string(start);
    auto f = h_edges.find(key);

    if (f == h_edges.end()) {
        // it's a new edge
        Mesh_Edge e1 = Mesh_Edge(start,index);
        Mesh_Edge e2 = Mesh_Edge(index,end);

        n_edges[start].push_back(e1);
        n_edges[index].push_back(e1);

        n_edges[end].push_back(e2);
        n_edges[index].push_back(e2);

        // put new vertex into h_edges
        h_edges[key] = index;
    }
}

void Mesh::split_faces(){
    clock_t begin = clock();

    //1. A hash table of to keep track of new edges added to the original vertices.
    // This does not include new edges between new vertices because they must be new.
    map<Mesh_Edge*,Mesh_Vertex*> midEdgeToPoint;

    //2. Add new faces to this vector array.
    vector<Mesh_Face> newFaces;

    //3. Clear all edges within the vertices
    for (int vi = 0; vi < vertices.size(); vi++) {
        vertices[vi].edges.clear();
    }

    //4. Iterate through all faces and add all new faces with vedges to the newFaces vector.
    for (Mesh_Face& face: faces) {
        //Original vertex indices
        int og1 = face.vert[0];
        int og2 = face.vert[1];
        int og3 = face.vert[2];

        //This will be the indices of all new vertices
        int new1 = vertices.size();
        int new2 = vertices.size() + 1;
        int new3 = vertices.size() + 2;

        //5. Get the midpoints between every vertex in current face.
        vector<Mesh_Edge> newEdges;

        //Midpoints of OG 1 and 2
        auto midPoint1 = get_midpoint_vertex(vertices[og1].position, vertices[og2].position);
        Mesh_Edge* m11 = new Mesh_Edge(og1, new1);
        Mesh_Edge* m12 = new Mesh_Edge(og2, new1);
        new1 = check_and_add(midEdgeToPoint, newEdges, m11, midPoint1);
        check_and_add(midEdgeToPoint, newEdges, m12, midPoint1);

        //Midpoints of OG 1 and 3
        auto midPoint2 = get_midpoint_vertex(vertices[og1].position, vertices[og3].position);
        Mesh_Edge* m21 = new Mesh_Edge(og1, new2);
        Mesh_Edge* m22 = new Mesh_Edge(og3, new2);
        new2 = check_and_add(midEdgeToPoint, newEdges, m21, midPoint2);
        check_and_add(midEdgeToPoint, newEdges, m22, midPoint2);

        //Midpoints of OG 2 and 3
        auto midPoint3 = get_midpoint_vertex(vertices[og2].position, vertices[og3].position);
        Mesh_Edge* m31 = new Mesh_Edge(og2, new3);
        Mesh_Edge* m32 = new Mesh_Edge(og3, new3);
        new3 = check_and_add(midEdgeToPoint, newEdges, m31, midPoint3);
        check_and_add(midEdgeToPoint, newEdges, m32, midPoint3);

        //6. Connect all midpoints

        // create 3 new edges: xy, xz, yz
        // put them into n_edges
        Mesh_Edge new12 = Mesh_Edge(new1,new2);
        Mesh_Edge new13 = Mesh_Edge(new1,new3);
        Mesh_Edge new23 = Mesh_Edge(new2,new3);

        vertices[new1].edges.push_back(new12);
        vertices[new1].edges.push_back(new13);

        vertices[new2].edges.push_back(new12);
        vertices[new2].edges.push_back(new23);

        vertices[new3].edges.push_back(new13);
        vertices[new3].edges.push_back(new23);


        //create new faces: one, two, three, four
        // a x y => one
        // x b z => two
        // x y z => three
        // y z c => four
        Mesh_Face one = Mesh_Face(og1,new1,new2);
        Mesh_Face two = Mesh_Face(new1,og2,new3);
        Mesh_Face three = Mesh_Face(new1,new2,new3);
        Mesh_Face four = Mesh_Face(new2,new3,og3);

        newFaces.push_back(one);
        newFaces.push_back(two);
        newFaces.push_back(three);
        newFaces.push_back(four);
    }

    faces = newFaces;

    //TODO: Delete dynamically allocated stuff from memory
    //for (auto& p : midEdgeToPoint) { delete p.first; delete p.second; }

    clock_t endtime = clock();

    double elapsed_secs = double(endtime - begin) / CLOCKS_PER_SEC;
    cout << "split faces run for " + to_string(elapsed_secs) << " seconds" << endl;


}

void Mesh::storeVBO() {
    vector<QVector3D> tri_vert, tri_bary;

    for(long f = 0; f < (long)faces.size(); f++) {
        tri_vert.push_back(vertices.at(faces[f].vert[0]).position);
        tri_vert.push_back(vertices.at(faces[f].vert[1]).position);
        tri_vert.push_back(vertices.at(faces[f].vert[2]).position);

        tri_bary.push_back(QVector3D(1,0,0));;
        tri_bary.push_back(QVector3D(0,1,0));;
        tri_bary.push_back(QVector3D(0,0,1));;
    }

    if(vertexBuffer.isCreated()) vertexBuffer.destroy();
    vertexBuffer.create();
    vertexBuffer.setUsagePattern( QOpenGLBuffer::StaticDraw );
    vertexBuffer.bind();
    vertexBuffer.allocate(&tri_vert[0] , sizeof( QVector3D ) * tri_vert.size());

    if(baryBuffer.isCreated()) baryBuffer.destroy();
    baryBuffer.create();
    baryBuffer.setUsagePattern( QOpenGLBuffer::StaticDraw );
    baryBuffer.bind();
    baryBuffer.allocate(&tri_bary[0] , sizeof( QVector3D ) * tri_bary.size());
}
