#define _CRT_SECURE_NO_WARNINGS 1

#include <vector>
#include <iostream>
#include <algorithm>
using namespace std;
#include <random>

#include <string>
#include <stdio.h>
#include <algorithm>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

thread_local default_random_engine engine;
thread_local uniform_real_distribution<double> uniform(0, 1);

#define M_PI 3.141592653589793238

#define LUM_SPHERIQUE false

class Vector {
public:
	Vector(double x = 0, double y = 0, double z = 0) :x(x), y(y), z(z) {}
	double norm2() { return x * x + y * y + z * z; }
	void normalize() {
		double n = sqrt(norm2());
		x /= n;
		y /= n;
		z /= n;
	}
	double x, y, z;
};

Vector operator+(const Vector& a, const Vector& b) {
	return Vector(a.x + b.x, a.y + b.y, a.z + b.z);
}
Vector operator-(const Vector& a, const Vector& b) {
	return Vector(a.x - b.x, a.y - b.y, a.z - b.z);
}
Vector operator/(const Vector& a, const double b) {
	return Vector(a.x / b, a.y / b, a.z / b);
}
Vector operator*(const Vector& a, const double b) {
	return Vector(a.x * b, a.y * b, a.z * b);
}
Vector operator*(const double b, const Vector& a) {
	return Vector(a.x * b, a.y * b, a.z * b);
}
Vector operator*(const Vector& a, const Vector& b) {
	return Vector(a.x * b.x, a.y * b.y, a.z * b.z);
}
bool operator!=(const Vector& a, const Vector& b) {
	if (a.x == b.x && a.y == b.y && a.z == b.z) {
		return false;
	}
	else {
		return true;
	}
}
double dot(const Vector& a, const Vector& b) {
	return a.x*b.x + a.y*b.y + a.z*b.z;
}
Vector cross(const Vector& a, const Vector& b) {
	return Vector(a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x);
}
double absolute(double a) {
	return max(a, -a);
}

double calcul_intensite(double I, Vector L, Vector P, Vector n) {
	Vector l = L - P; // Le vecteur de P vers la source de lumière L
	double d2 = l.norm2(); // La distance au carré entre P et la source de lumière L
	l.normalize();
	double intensite = I * max(0., dot(n, l)) / (4 * M_PI * M_PI * d2);
	return intensite;
}

Vector random_cos(const Vector& N) {
	double r1 = uniform(engine);
	double r2 = uniform(engine);
	double s = sqrt(1 - r2);
	Vector v(cos(2 * M_PI*r1)*s, sin(2 * M_PI*r1)*s, sqrt(r2));
	Vector T1;

	if (N.x <= N.y && N.x <= N.z) {
		T1 = Vector(0, -N.z, N.y);
	}
	else if (N.y <= N.x && N.y <= N.z) {
		T1 = Vector(-N.z, 0, N.x);
	}
	else {
		T1 = Vector(-N.y, N.x, 0);
	}
	T1.normalize();
	Vector T2 = cross(N, T1);

	return T1 * v.x + T2 * v.y + N * v.z;
}

class Ray {
public:
	Ray(const Vector& C = Vector(), const Vector& u = Vector(1, 1, 1)) : C(C), u(u) {}

	Vector C, u;
};

class Object {
public:
	Object(const Vector& albedo = Vector(1., 1., 1.), bool miroir = false, bool transparent = false) : albedo(albedo), miroir(miroir), transparent(transparent) {}
	Vector albedo; // couleur du triangle
	bool miroir; // le triangle est-il miroir ?
	bool transparent; // le triangle est-il transparent ?

	virtual bool intersect(Ray ray, Vector& P, Vector& N, double& t, int& ids) = 0;
};

class TriangleIndices {
public:
	TriangleIndices(int vtxi = -1, int vtxj = -1, int vtxk = -1, int ni = -1, int nj = -1, int nk = -1, int uvi = -1, int uvj = -1, int uvk = -1, int group = -1, bool added = false) : vtxi(vtxi), vtxj(vtxj), vtxk(vtxk), uvi(uvi), uvj(uvj), uvk(uvk), ni(ni), nj(nj), nk(nk), group(group) {
	};
	int vtxi, vtxj, vtxk; // indices within the vertex coordinates array
	int uvi, uvj, uvk;  // indices within the uv coordinates array
	int ni, nj, nk;  // indices within the normals array
	int group;       // face group
};

class Triangle : public Object {
public:
	Triangle(const Vector& A, const Vector& B, const Vector& C, const Vector& albedo = Vector(1., 1., 1.), bool miroir = false, bool transparent = false) : A(A), B(B), C(C) {
		this->albedo = albedo;
		this->miroir = miroir;
		this->transparent = transparent;
	}
	Vector A;
	Vector B;
	Vector C;

	bool intersect(Ray ray, Vector& P, Vector& N, double& t, int& ids) {

		N = cross(B - A, C - A);
		N.normalize();

		t = dot(A - ray.C, N) / dot(ray.u, N);

		if (t >= 0) {
			P = ray.C + ray.u * t;
			Vector barycentric = get_barycentric_coordinates(P, A, B, C);
			if (barycentric.x >= 0 && barycentric.x <= 1 && barycentric.y >= 0 && barycentric.y <= 1 && barycentric.z >= 0 && barycentric.z <= 1) {
				return true;
			}
			else {
				return false;
			}
		}
		else {
			return false;
		}
	}

	Vector get_center() {
		return (A + B + C) / 3;
	}

	Vector get_barycentric_coordinates(Vector& P, Vector& A, Vector& B, Vector& C) {
		double beta = (dot(P - A, B - A) * (C - A).norm2() - dot(C - A, B - A) * dot(P - A, C - A)) / ((B - A).norm2() * (C - A).norm2() - dot(B - A, C - A) * dot(B - A, C - A));
		double gamma = (dot(P - A, C - A) * (B - A).norm2() - dot(B - A, C - A) * dot(P - A, B - A)) / ((C - A).norm2() * (B - A).norm2() - dot(C - A, B - A) * dot(C - A, B - A));
		double alpha = 1 - beta - gamma;

		return Vector(alpha, beta, gamma);
	}
};

class Bbox {
public:
	Bbox() : P1(Vector()), P2(Vector()) {}
	Bbox(const Vector& A, const Vector& B) : P1(A), P2(B) {}
	Vector P1;
	Vector P2;

	void setP1(const Vector& A) {
		P1 = A;
	}

	void setP2(const Vector& B) {
		P2 = B;
	}

	bool isInside(Triangle& t) {
		Vector center = t.get_center();
		if (center.x > P1.x && center.x < P2.x && center.y > P1.y && center.y < P2.y && center.z > P1.z && center.z < P2.z) {
			return true;
		}
		else {
			return false;
		}
	}

	vector<Bbox> divide() {
		Vector diag = P2 - P1;
		double xSize = absolute(P2.x - P1.x);
		double ySize = absolute(P2.y - P1.y);
		double zSize = absolute(P2.z - P1.z);

		Bbox b1, b2;

		if (xSize >= ySize && xSize >= zSize) {
			b1 = Bbox(Vector(P1.x, P1.y, P1.z), Vector(diag.x / 2, P2.y, P2.z));
			b2 = Bbox(Vector(diag.x / 2, P1.y, P1.z), Vector(P2.x, P2.y, P2.z));
		}
		else if (ySize >= xSize && ySize >= zSize) {
			b1 = Bbox(Vector(P1.x, P1.y, P1.z), Vector(P2.x, diag.y / 2, P2.z));
			b2 = Bbox(Vector(P1.x, diag.y / 2, P1.z), Vector(P2.x, P2.y, P2.z));
		}
		else {
			b1 = Bbox(Vector(P1.x, P1.y, P1.z), Vector(P2.x, P2.y, diag.z / 2));
			b2 = Bbox(Vector(P1.x, P1.y, diag.z / 2), Vector(P2.x, P2.y, P2.z));
		}

		vector<Bbox> Bboxes;
		Bboxes.push_back(b1);
		Bboxes.push_back(b2);

		return Bboxes;
	}
};

class BVHNode {
public:
	BVHNode() : debut(0), fin(0), b(Bbox()), isElementary(false) {}
	BVHNode(int debut, int fin, Bbox b) : debut(debut), fin(fin), b(b), isElementary(false) {}
	BVHNode(int debut, int fin, Bbox b, bool isElementary) : debut(debut), fin(fin), b(b), isElementary(isElementary) {}
	BVHNode *fg, *fd;
	Bbox b;
	int debut, fin;
	bool isElementary;

	void quick_sort(vector<TriangleIndices> indices, vector<Vector> vertices, int max_it) {
		if (fin - debut <= 4 || max_it == 0) { // condition initiale : le BVHNode de gauche contient tout et le BVHNode de droite est vide
			fg = new BVHNode(debut, fin, b, true);
			fd = new BVHNode(fin, fin, Bbox(), true);
		}
		else {
			vector<Bbox> dividedBboxes = b.divide();
			int pivot = debut;
			vector<TriangleIndices> newIndices = indices;
			for (int i = debut; i <= fin; i++) {
				TriangleIndices ti = indices[i];
				Vector A = vertices[ti.vtxi];
				Vector B = vertices[ti.vtxj];
				Vector C = vertices[ti.vtxk];

				Triangle face(A, B, C);

				if (dividedBboxes[0].isInside(face)) {
					newIndices[i] = indices[pivot];
					newIndices[pivot] = indices[i];
					pivot++;
				}
			}

			fg = new BVHNode(debut, pivot - 1, dividedBboxes[0]);
			fd = new BVHNode(pivot, fin, dividedBboxes[1]);

			fg->quick_sort(indices, vertices, max_it - 1);
			fd->quick_sort(indices, vertices, max_it - 1);
		}
	}
};

class Geometry : public Object {
public:
	~Geometry() {}
	Geometry() {};
	Geometry(double scale = 1, Vector trans = Vector(), const Vector& albedo = Vector(1., 1., 1.), bool miroir = false, bool transparent = false) : scale(scale), trans(trans), bvhnode(BVHNode()) {
		this->albedo = albedo;
		this->miroir = miroir;
		this->transparent = transparent;
	}
	double scale;
	Vector trans;

	BVHNode bvhnode;

	void initialize_BVHNode() {
		Bbox b = getBbox(indices);
		bvhnode = BVHNode(0, indices.size() - 1, b);
	}

	void quick_sort_BVHNode(int max_it) {
		initialize_BVHNode();
		bvhnode.quick_sort(indices, vertices, max_it);
	}

	Bbox getBbox(vector<TriangleIndices> triangleIndices) {

		double xmin = 0, ymin = 0, zmin = 0, xmax = 0, ymax = 0, zmax = 0;

		for (int i = 0; i < vertices.size(); i++) {
			Vector vertex = vertices[i];

			xmin = min(vertex.x, xmin);
			xmax = max(vertex.x, xmax);
			ymin = min(vertex.y, ymin);
			ymax = max(vertex.y, ymax);
			zmin = min(vertex.z, zmin);
			zmax = max(vertex.z, zmax);

		}

		Bbox b(Vector(xmin, ymin, zmin), Vector(xmax, ymax, zmax));

		return b;
	}

	bool intersect(Ray ray, Vector& P, Vector& N, double& t, int& ido) {
		Vector nearestP, nearestN;
		bool firstRound = true;
		for (int i = 0; i < indices.size(); i++) {

			TriangleIndices ti = indices[i];
			Vector A = vertices[ti.vtxi];
			Vector B = vertices[ti.vtxj];
			Vector C = vertices[ti.vtxk];

			Triangle face(A, B, C);

			if (face.intersect(ray, P, N, t, ido)) {
				if (firstRound) {
					nearestP = P;
					nearestN = N;
					ido = i;
					firstRound = false;
				}
				else if ((P - C).norm2() < (nearestP - C).norm2()) {
					nearestP = P;
					nearestN = N;
					ido = i;
				}
			}
		}
		if (!firstRound) {
			P = nearestP;
			N = nearestN;
			return true;
		}
		else {
			return false;
		}
	}

	void readOBJ(const char* obj, bool load_textures) {
		// fonction originalement écrite par Nicolas BONNEEL pour lire un fichier .obj
		char matfile[255];
		char grp[255];

		FILE* f;
		f = fopen(obj, "r");
		int curGroup = -1;
		while (!feof(f)) {
			char line[255];
			if (!fgets(line, 255, f)) break;

			std::string linetrim(line);
			linetrim.erase(linetrim.find_last_not_of(" \r\t") + 1);
			strcpy(line, linetrim.c_str());

			if (line[0] == 'u' && line[1] == 's') {
				sscanf(line, "usemtl %[^\n]\n", grp);
				curGroup++;
			}

			if (line[0] == 'v' && line[1] == ' ') {
				Vector vec;

				Vector col;
				if (sscanf(line, "v %lf %lf %lf %lf %lf %lf\n", &vec.x, &vec.y, &vec.z, &col.x, &col.y, &col.z) == 6) {
					col.x = std::min(1., std::max(0., col.x));
					col.y = std::min(1., std::max(0., col.y));
					col.z = std::min(1., std::max(0., col.z));

					vertices.push_back(vec);
					vertexcolors.push_back(col);

				}
				else {
					sscanf(line, "v %lf %lf %lf\n", &vec.x, &vec.y, &vec.z);
					vertices.push_back(vec);
				}
			}
			if (line[0] == 'v' && line[1] == 'n') {
				Vector vec;
				sscanf(line, "vn %lf %lf %lf\n", &vec.x, &vec.y, &vec.z);
				normals.push_back(vec);
			}
			if (line[0] == 'v' && line[1] == 't') {
				Vector vec;
				sscanf(line, "vt %lf %lf\n", &vec.x, &vec.y);
				uvs.push_back(vec);
			}
			if (line[0] == 'f') {
				TriangleIndices t;
				int i0, i1, i2, i3;
				int j0, j1, j2, j3;
				int k0, k1, k2, k3;
				int nn;
				t.group = curGroup;

				char* consumedline = line + 1;
				int offset;

				nn = sscanf(consumedline, "%u/%u/%u %u/%u/%u %u/%u/%u%n", &i0, &j0, &k0, &i1, &j1, &k1, &i2, &j2, &k2, &offset);
				if (nn == 9) {
					if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
					if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
					if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
					if (j0 < 0) t.uvi = uvs.size() + j0; else	t.uvi = j0 - 1;
					if (j1 < 0) t.uvj = uvs.size() + j1; else	t.uvj = j1 - 1;
					if (j2 < 0) t.uvk = uvs.size() + j2; else	t.uvk = j2 - 1;
					if (k0 < 0) t.ni = normals.size() + k0; else	t.ni = k0 - 1;
					if (k1 < 0) t.nj = normals.size() + k1; else	t.nj = k1 - 1;
					if (k2 < 0) t.nk = normals.size() + k2; else	t.nk = k2 - 1;
					indices.push_back(t);
				}
				else {
					nn = sscanf(consumedline, "%u/%u %u/%u %u/%u%n", &i0, &j0, &i1, &j1, &i2, &j2, &offset);
					if (nn == 6) {
						if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
						if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
						if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
						if (j0 < 0) t.uvi = uvs.size() + j0; else	t.uvi = j0 - 1;
						if (j1 < 0) t.uvj = uvs.size() + j1; else	t.uvj = j1 - 1;
						if (j2 < 0) t.uvk = uvs.size() + j2; else	t.uvk = j2 - 1;
						indices.push_back(t);
					}
					else {
						nn = sscanf(consumedline, "%u %u %u%n", &i0, &i1, &i2, &offset);
						if (nn == 3) {
							if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
							if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
							if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
							indices.push_back(t);
						}
						else {
							nn = sscanf(consumedline, "%u//%u %u//%u %u//%u%n", &i0, &k0, &i1, &k1, &i2, &k2, &offset);
							if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
							if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
							if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
							if (k0 < 0) t.ni = normals.size() + k0; else	t.ni = k0 - 1;
							if (k1 < 0) t.nj = normals.size() + k1; else	t.nj = k1 - 1;
							if (k2 < 0) t.nk = normals.size() + k2; else	t.nk = k2 - 1;
							indices.push_back(t);
						}
					}
				}

				consumedline = consumedline + offset;

				while (true) {
					if (consumedline[0] == '\n') break;
					if (consumedline[0] == '\0') break;
					nn = sscanf(consumedline, "%u/%u/%u%n", &i3, &j3, &k3, &offset);
					TriangleIndices t2;
					t2.group = curGroup;
					if (nn == 3) {
						if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
						if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
						if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
						if (j0 < 0) t2.uvi = uvs.size() + j0; else	t2.uvi = j0 - 1;
						if (j2 < 0) t2.uvj = uvs.size() + j2; else	t2.uvj = j2 - 1;
						if (j3 < 0) t2.uvk = uvs.size() + j3; else	t2.uvk = j3 - 1;
						if (k0 < 0) t2.ni = normals.size() + k0; else	t2.ni = k0 - 1;
						if (k2 < 0) t2.nj = normals.size() + k2; else	t2.nj = k2 - 1;
						if (k3 < 0) t2.nk = normals.size() + k3; else	t2.nk = k3 - 1;
						indices.push_back(t2);
						consumedline = consumedline + offset;
						i2 = i3;
						j2 = j3;
						k2 = k3;
					}
					else {
						nn = sscanf(consumedline, "%u/%u%n", &i3, &j3, &offset);
						if (nn == 2) {
							if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
							if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
							if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
							if (j0 < 0) t2.uvi = uvs.size() + j0; else	t2.uvi = j0 - 1;
							if (j2 < 0) t2.uvj = uvs.size() + j2; else	t2.uvj = j2 - 1;
							if (j3 < 0) t2.uvk = uvs.size() + j3; else	t2.uvk = j3 - 1;
							consumedline = consumedline + offset;
							i2 = i3;
							j2 = j3;
							indices.push_back(t2);
						}
						else {
							nn = sscanf(consumedline, "%u//%u%n", &i3, &k3, &offset);
							if (nn == 2) {
								if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
								if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
								if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
								if (k0 < 0) t2.ni = normals.size() + k0; else	t2.ni = k0 - 1;
								if (k2 < 0) t2.nj = normals.size() + k2; else	t2.nj = k2 - 1;
								if (k3 < 0) t2.nk = normals.size() + k3; else	t2.nk = k3 - 1;
								consumedline = consumedline + offset;
								i2 = i3;
								k2 = k3;
								indices.push_back(t2);
							}
							else {
								nn = sscanf(consumedline, "%u%n", &i3, &offset);
								if (nn == 1) {
									if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
									if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
									if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
									consumedline = consumedline + offset;
									i2 = i3;
									indices.push_back(t2);
								}
								else {
									consumedline = consumedline + 1;
								}
							}
						}
					}
				}

			}

		}
		fclose(f);

	}

	vector<TriangleIndices> indices;
	vector<Vector> vertices;
	vector<Vector> normals;
	vector<Vector> uvs;
	vector<Vector> vertexcolors;

};

class Sphere : public Object {
public:
	Sphere(const Vector& O = Vector(), double R = 10, Vector albedo = Vector(1, 1, 1), const bool miroir = false,
		const bool transparent = false, const double n = 1, const bool isLum = false) :
		O(O), R(R), n(n), isLum(isLum) {
		this->albedo = albedo;
		this->miroir = miroir;
		this->transparent = transparent;
	}
	Vector O; // centre de la sphère
	double R; // rayon de la sphère
	double n; // indice de réfraction de la sphère
	bool isLum; // la sphère est-elle une sphère lumineuse ?

	bool intersect(Ray ray, Vector& P, Vector& N, double& t, int& ids) {
		double a = 1;
		double b = 2 * dot(ray.u, ray.C - O);
		double c = (ray.C - O).norm2() - R * R;
		double delta = b * b - 4 * a*c;
		if (delta >= 0) {
			if (-b - sqrt(delta) > 0) {
				t = (-b - sqrt(delta)) / (2 * a);
			}
			else if (-b + sqrt(delta) > 0) {
				t = (-b + sqrt(delta)) / (2 * a);
			}
			else {
				return false;
			}
			P = ray.C + ray.u * t;
			N = P - O;
			N.normalize();
			return true;
		}
		else {
			return false;
		}
	}
};

class Scene {
public:
	Scene(const vector<Object*> objects = vector<Object*>(), const Vector& C = Vector(),
		const Vector& L = Vector(), const double I = 0, const double n = 1) :
		objects(objects), C(C), L(L), I(I), n(n) {}
	vector<Object*> objects;
	Vector C; // position de la caméra / de l'oeil du spectateur
	Vector L; // source de lumière
	double I; // luminosité de la source de lumière
	double n; // indice de réfraction de la scène

	void addObject(Object* o) {
		objects.push_back(o);
	}

	bool intersect(Ray ray, Vector& P, Vector& N, double& t, int& ids) {
		Vector nearestP, nearestN;
		bool firstRound = true;
		vector<double> intersections;
		for (int unsigned i = 0; i < objects.size(); i++) {
			Object* object = objects[i];
			if (object->intersect(ray, P, N, t, ids)) {
				if (firstRound) {
					nearestP = P;
					nearestN = N;
					ids = i;
					firstRound = false;
				}
				else if ((P - C).norm2() < (nearestP - C).norm2()) {
					nearestP = P;
					nearestN = N;
					ids = i;
				}
			}
		}
		if (!firstRound) {
			P = nearestP;
			N = nearestN;
			return true;
		}
		else {
			return false;
		}
	}

	Vector get_color(Vector& P, Vector& N, Ray ray, int nbRebonds) {
		Vector rayColor;
		double t;
		int ids;
		double n1, n2;
		n1 = n;
		if (nbRebonds == 0) {
			return Vector(0., 0., 0.);
		}

		/* ---------- Si la lumière est sphérique ------------- */
		if (LUM_SPHERIQUE) {}

		/* ---------- Si la lumière est ponctuelle ------------- */
		else { 
			if (intersect(ray, P, N, t, ids)) {
				Object* o = objects[ids];
				/* ---------- 1. Cas de la sphère MIROIR ------------- */
				if (o->miroir) { // La sphère est une sphère miroir, on enclanche le processus récursif
					if (nbRebonds == 0) {
						return Vector(0., 0., 0.);
					}
					else {
						Vector uMiroir = ray.u - N * 2 * dot(ray.u, N);
						Vector PMiroir = P + N * 0.0001;
						Ray rayMiroir(PMiroir, uMiroir);
						return get_color(PMiroir, uMiroir, rayMiroir, nbRebonds - 1);
					}
				}

				/* ---------- 2. Cas de la sphère TRANSPARENTE ------------- */
				/*else if (s.transparent) {
					if (nbRebonds == 0) {
						return Vector();
					}
					else {
						int sensN;
						if (dot(N, ray.u) < 0) { // On est en train de rentrer dans la sphère, donc N et u ont un produit scalaire négatif
							n2 = s.n;
							sensN = -1;
						}
						else { // On est à l'intérieur de la sphère, en train d'en sortir
							n1 = s.n;
							n2 = n;
							sensN = 1;
						}
						// if ce qu'il y a dans sqrt est positif (traiter l'autre cas)
						Vector Tn = N * sensN * sqrt(1 - (n1 / n2) * (n1 / n2) * (1 - dot(ray.u, N) * dot(ray.u, N)));
						Vector Tt = (ray.u - N * dot(ray.u, N)) * (n1 / n2);
						Vector uRefraction = Tn + Tt;
						Vector PRefraction = P + N * sensN * 0.0001;
						Ray rayRefraction(PRefraction, uRefraction);
						nbRebonds--;
						return get_color(PRefraction, uRefraction, rayRefraction, nbRebonds);
					}
				}*/

				else {
					/* ---------- Contribution directe ------------- */
					double intensite = calcul_intensite(I, L, P, N);
					Vector PL = L - P;
					double distance2_PL = PL.norm2();
					PL.normalize();
					Ray rayPrime(P + PL * 0.0001, PL);
					Vector PPrime, NPrime;
					double tPrime;
					int idsPrime;
					bool ombre = intersect(rayPrime, PPrime, NPrime, tPrime, idsPrime); // On doit rajouter t la distance entre P (origine de rayPrime) et PPrime : c'est directement la solution t résolue par intersect
					double distance2_PPPrime = (PPrime - P).norm2();
					if (ombre && distance2_PPPrime < distance2_PL) {
						/* ---------- On est dans une zone d'ombre ------------- */
						rayColor = Vector(0, 0, 0);
					}
					else {
						rayColor = o->albedo * intensite;
					}

					/* ---------- Contribution indirecte ------------- */
					nbRebonds--;
					Vector reflechi = random_cos(N);
					Ray ray_reflechi(P + N * 0.001, reflechi);
					rayColor = rayColor + o->albedo * get_color(P, N, ray_reflechi, nbRebonds); // aussi / M_PI ? 
				}
			}
			else { // Sinon, mettre les pixels en noir
				rayColor = Vector(0, 0, 0);
			}
		}

		return rayColor;
	}
};

int main() {
	int resolution = 3;
	int W = 100 * resolution;
	int H = 100 * resolution;
	double fov = 60 * M_PI / 180.;

	Vector rouge(1, 0, 0), vert(0, 1, 0), bleu(0, 0, 1), cyan(0, 1, 1), violet(1, 0, 1), jaune(1, 1, 0), blanc(1, 1, 1), noir(0, 0, 0);

	Vector C(0, 0, 55);
	Vector L(-10, -20, 40);
	double I = 50000000000;

	vector<unsigned char> image(W * H * 3, 0);
	vector<Object*> objects;

	Scene scene = Scene(objects, C, L, I);

	if (LUM_SPHERIQUE) {
		scene.addObject(new Sphere(Vector(-20, -40, 50), 15., blanc, false, false, 1, true));
	}

	/* ----- Ajouts de sphères à la scène ------------ */
	//spheres.push_back(Sphere(Vector(0, -3, 10), 12., blanc, false, false));
	//scene.addObject(new Sphere(Vector(-13, 0, 3), 10., blanc, true, false));
	// scene.addObject(new Sphere(Vector(0, -15, 5), 7., rouge, true, false));
	//spheres.push_back(Sphere(Vector(13, 0, 0), 10., blanc, false, true, 1.4));
	//scene.addObject(new Sphere(Vector(13, 0, 10), 5., blanc, false, false));
	//spheres.push_back(Sphere(Vector(-20, 3, 25), 2., jaune, false, false));
	//spheres.push_back(Sphere(Vector(10, -10, -25), 6., cyan, false, false));
	//scene.addObject(new Sphere(Vector(0, -1000, 0), 950., vert, false, false));
	//scene.addObject(new Sphere(Vector(0, 1000, 0), 990., bleu, false, false));
	/*scene.addObject(new Sphere(Vector(0, 0, -5000), 4990., rouge, false, false));
	//spheres.push_back(Sphere(Vector(0, 0, 1000), 960., jaune, false, false));
	scene.addObject(new Sphere(Vector(1000, 0, 0), 965., jaune, false, false));
	scene.addObject(new Sphere(Vector(-1000, 0, 0), 965., cyan, false, false));

	/* ----- Ajouts d'un triangle miroir à la scène ------------ */
	/*Triangle* t0 = new Triangle(Vector(0, -20, 3), Vector(13, -20, 20), Vector(0, 20, 3), blanc, true);
	scene.addObject(t0);*/

	/* ----- Ajouts du maillage à la scène ------------ */
	Geometry* g0 = new Geometry(30, Vector(0, 25, 0));
	g0->readOBJ("Beautiful Girl.obj", false);
	cout << g0->vertices.size() << endl;
	for (int i = 0; i < g0->vertices.size(); i++) {
		g0->vertices[i] = Vector(g0->vertices[i].x, -g0->vertices[i].z, g0->vertices[i].y) * g0->scale + g0->trans;
	}

	/* ----- Création et Quicksort des BVHNodes et des Bbox ------------ */
	/*g0->quick_sort_BVHNode(100);*/
	scene.addObject(g0);

	int k = 0;
#pragma omp parallel for
	for (int i = 0; i < H; i++) {
		cout << k << endl;
		k++;
		for (int j = 0; j < W; j++) {
			Vector P; // Point d'intersection (s'il existe) entre le rayon et la sphère
			Vector N; // La normale au point P sur la sphère
			Vector V; // Vecteur directeur du rayon à construire
			Ray ray; // Rayon envoyé vers le pixel (i,j)
			double x;
			double y;
			double R;
			double u;
			double v; // u et v, construits à partir de x, y et R, permettent d'appliquer un antialiasing

			int nbRays = 1; // nombre de rayon de réflexion pour la lumière diffuse
			Vector meanRayColor(0, 0, 0);
			for (int k = 0; k < nbRays; k++) {
				x = uniform(engine);
				y = uniform(engine);
				R = sqrt(-2 * log(x));
				u = R * cos(2 * M_PI * y) * 0.5;
				v = R * sin(2 * M_PI * y) * 0.5;
				V = Vector(j - W / 2 - 0.5 + u, i - H / 2 - 0.5 + v, -H / (2 * tan(fov / 2.)));
				V.normalize();
				ray = Ray(C, V);
				Vector rayColor = scene.get_color(P, N, ray, 5);
				meanRayColor = meanRayColor + rayColor;
			}
			meanRayColor = meanRayColor / nbRays;

			image[(i*W + j) * 3 + 0] = min(255., pow(meanRayColor.x, 0.45));
			image[(i*W + j) * 3 + 1] = min(255., pow(meanRayColor.y, 0.45));
			image[(i*W + j) * 3 + 2] = min(255., pow(meanRayColor.z, 0.45));
		}
	}
	stbi_write_png("maillage_beautiful_girl_boxes_2.png", W, H, 3, &image[0], 0);

	return 0;
}
