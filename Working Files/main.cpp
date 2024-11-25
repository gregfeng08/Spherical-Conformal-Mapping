#include <fstream>
#include <vector>
#include <set>

#include "OBJFileReader.h"
#include "Solid.h"
#include "iterators.h"
#include "SolidDelegate.h"
#include "Trait.h"

using namespace MeshLib;

double totalEnergy(Solid &mesh, bool isHarmonic) {
	SolidEdgeIterator eiter(&mesh);
	double ret = 0;
	for(; !eiter.end();++eiter) {
		Edge* e = *eiter;
		HalfEdge* he = e->halfedge(0);
		Point v1 = he->source()->point();
		Point v2 = he->target()->point();

		double energy = (v1-v2).norm2();
		if(isHarmonic) {
			energy*=e->kuv;
		}
		ret+=energy;
	}
	return ret;
}

void gauss_map(Solid &mesh) {
	SolidVertexIterator viter(&mesh);
    for(; !viter.end(); ++viter) {
        Vertex* v = *viter;

        Point averageNormal(0, 0, 0);
        int faceCount = 0;

        HalfEdge* start = v->halfedge();

        HalfEdge* he = start;
        do {
            Face* face = he->face();
            if (face != nullptr) {
                averageNormal += face->norm();
                ++faceCount;
            }
            he = he->he_sym()->he_next();
        } while (he != start && he != nullptr);

        if (faceCount > 0) {
            averageNormal /= faceCount;
            averageNormal /= averageNormal.norm(); 
        }
        v->point() = averageNormal;
    }
}

void computeAbsDeriv(Solid &mesh, bool isHarmonic) {
	SolidVertexIterator viter(&mesh);
	for(; !viter.end();++viter) {
		Vertex* v = *viter;
		Point laplacian = Point(0,0,0);
		VertexVertexIterator vviter(v);
		for(; !vviter.end();++vviter) {
			Vertex* adjVert = *vviter;
			Point vertLap = Point(0,0,0);
			vertLap = v->point()-adjVert->point();
			if(isHarmonic) {
				vertLap=vertLap*(&mesh)->vertexEdge(v,adjVert)->kuv;
			}
			laplacian+=(v->point()-adjVert->point());
		}
		v->absoluteDeriv = (laplacian-(v->point()*(laplacian*v->point())));
	}
}

void computeKUV(Solid &mesh) {
	SolidEdgeIterator eiter(&mesh);
	for(; !eiter.end();++eiter) {
		Edge* e = *eiter;
		HalfEdge* he1 = e->halfedge(0);
		HalfEdge* he2 = e->halfedge(1);

		Point p1 = he1->source()->point();
		Point p2 = he1->he_next()->target()->point();
		Point p3 = he1->target()->point();
		Point p4 = he1->he_sym()->he_next()->target()->point();

		double alphaUV = (p1-p2)*(p3-p2)/(2*((p1-p2)^(p3-p2)).norm());
		double betaUV = (p1-p4)*(p3-p4)/(2*((p1-p4)^(p3-p4)).norm());

		e->kuv = alphaUV+betaUV;
	}
}

void computeFaceAreas(Solid &mesh) {
	SolidFaceIterator fiter(&mesh);
	for(; !fiter.end();++fiter) {
		Face* face = *fiter;
		HalfEdge* he1 = face->halfedge();
		Point crossArg1 = he1->target()->point()-he1->source()->point();
		Point crossArg2 = he1->he_prev()->source()->point()-he1->source()->point();

		double crossNorm = (crossArg1^crossArg2).norm()/2;
		face->area = crossNorm;
	}
}

void recenterMesh(Solid &mesh) {
	computeFaceAreas(mesh);

	SolidVertexIterator viter(&mesh);
	for(; !viter.end();++viter) {
		Vertex* v = *viter;
		double sumFaceAreas = 0;
        int faceCount = 0;

        HalfEdge* start = v->halfedge();
        HalfEdge* he = start;

        do {
            Face* face = he->face();
			sumFaceAreas+=face->area;
            if (face != nullptr) {
                ++faceCount;
            }
            he = he->he_sym()->he_next();
        } while (he != start && he != nullptr);

		v->area = sumFaceAreas/3;
	}

	Point massCenter = Point(0,0,0);
	double mass = 0;
	SolidVertexIterator vviter(&mesh);
	for(; !vviter.end();++vviter) {
		Vertex* v = *vviter;
		massCenter+=v->point()*v->area;
		mass+=v->area;
	}
	massCenter/=mass;

	SolidVertexIterator vvviter(&mesh);
	for(; !vvviter.end();++vvviter) {
		Vertex* v = *vvviter;
		v->point()-=massCenter;
		v->point()/=v->point().norm();
	}

}

int main(int argc, char *argv[]) {
    // Read in the obj file
	Solid mesh;
	Solid newMesh;
	OBJFileReader of;
	std::ifstream in(argv[1]);
	of.readToSolid(&mesh, in);
	mesh.copyto(newMesh);

    //MAIN PROGRAM

	//Tuette Map Calculation
	gauss_map(newMesh);
	computeKUV(newMesh);

	double prevEnergy = totalEnergy(newMesh, false);
	std::cout << "Init Tuette:" << prevEnergy << std::endl;
	double currEnergy = 0;
	double threshold = 1e-5;
	double stepSize = 1e-4;
	int ctr = 0;
	int maxIter = 2e9;
	while(true) {
		std::cout << ctr++ << std::endl;
		computeAbsDeriv(newMesh, false);
		SolidVertexIterator viter(&newMesh);
		for(; !viter.end();++viter) {
			Vertex* v = *viter;
			v->point() = v->point()-v->absoluteDeriv*stepSize;
			v->point() /= v->point().norm();
		}
		currEnergy = totalEnergy(newMesh, false);
		if(std::abs(currEnergy-prevEnergy)<threshold||ctr>maxIter) {
			break;
		} else {
			prevEnergy = currEnergy;
		}
		std::cout << currEnergy << "," << prevEnergy << "," << stepSize << std::endl;
	}

	//Spherical Conformal Mapping Calculation
	double prevEnergyHarmonic = totalEnergy(newMesh, true);
	std::cout << "Init Harmonic:" << prevEnergy << std::endl;
	double currEnergyHarmonic = 0;
	double thresholdHarmonic = 1e-5;
	double stepSizeHarmonic = 1e-4;

	int ctrHarmonic = 0;
	int maxIterHarmonic = 2e9;

	while(true) {
		std::cout << ctrHarmonic++ << std::endl;
		computeAbsDeriv(newMesh, true);
		SolidVertexIterator viter(&mesh);
		for(; !viter.end();++viter) {
			Vertex* v = *viter;
			v->point() = v->point()-v->absoluteDeriv*stepSizeHarmonic;
			v->point() /= v->point().norm();
		}
		currEnergyHarmonic = totalEnergy(newMesh, true);
		recenterMesh(newMesh);
		if(std::abs(currEnergyHarmonic-prevEnergyHarmonic)<thresholdHarmonic||ctrHarmonic>maxIterHarmonic) {
			break;
		} else {
			prevEnergyHarmonic = currEnergyHarmonic;
		}
		std::cout << "Harmonic:" << currEnergyHarmonic << "," << prevEnergyHarmonic << "," << stepSizeHarmonic << std::endl;
	}

    // Write out the resultant obj file
	int vObjID = 1;
	std::map<int, int> vidToObjID;

	std::ofstream os(argv[2]);

	SolidVertexIterator iter(&newMesh);

	for(; !iter.end(); ++iter)
	{
		Vertex *v = *iter;
		Point p = v->point();
		os << "v " << p[0] << " " << p[1] << " " << p[2] << std::endl;
		vidToObjID[v->id()] = vObjID++;
	}
	os << "# " << (unsigned int)newMesh.numVertices() << " vertices" << std::endl;

	float u = 0.0, v = 0.0;
	for(iter.reset(); !iter.end(); ++iter)
	{
		Vertex *vv = *iter;
		std::string key( "uv" );
		std::string s = Trait::getTraitValue (vv->string(), key );
		if( s.length() > 0 )
		{
			sscanf( s.c_str (), "%f %f", &u, &v );
		}
		os << "vt " << u << " " << v << std::endl;
	}
	os << "# " << (unsigned int)newMesh.numVertices() << " texture coordinates" << std::endl;

	SolidFaceIterator fiter(&newMesh);
	for(; !fiter.end(); ++fiter)
	{
		Face *f = *fiter;
		FaceVertexIterator viter(f);
		os << "f " ;
		for(; !viter.end(); ++viter)
		{
			Vertex *v = *viter;
			os << vidToObjID[v->id()] << "/" << vidToObjID[v->id()] << " ";
		}
		os << std::endl;
	}
	os.close();
}