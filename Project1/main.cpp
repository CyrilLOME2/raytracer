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

class Sphere {
public:
	Sphere(const Vector& O = Vector(), double R = 10, Vector albedo = Vector(1, 1, 1), const bool miroir = false,
		const bool transparent = false, const double n = 1, const bool isLum = false) :
		O(O), R(R), albedo(albedo), miroir(miroir), transparent(transparent), n(n), isLum(isLum) {}
	Vector O; // centre de la sphère
	double R; // rayon de la sphère
	bool miroir; // la sphère est-elle miroir ?
	bool transparent; // la sphère est-elle transparente ?
	double n; // indice de réfraction de la sphère
	bool isLum; // la sphère est-elle une sphère lumineuse ?
	Vector albedo; // couleur de la sphere

	bool intersect(Ray ray, Vector& P, Vector& N, double& t) {
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
	Scene(const std::vector<Sphere> spheres = std::vector<Sphere>(), const Vector& C = Vector(),
		const Vector& L = Vector(), const double I = 0, const double n = 1) :
		spheres(spheres), C(C), L(L), I(I), n(n) {}
	std::vector<Sphere> spheres;
	Vector C; // position de la caméra / de l'oeil du spectateur
	Vector L; // source de lumière
	double I; // luminosité de la source de lumière
	double n; // indice de réfraction de la scène

	void addSphere(Sphere &s) {}

	bool intersect(Ray ray, Vector& P, Vector& N, double& t, int& ids) {

		Vector nearestP, nearestN;
		bool firstRound = true;
		std::vector<double> intersections;
		for (int unsigned i = 0; i < spheres.size(); i++) {
			Sphere sphere = spheres[i];
			if (sphere.intersect(ray, P, N, t)) {
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
		if (LUM_SPHERIQUE) { 
			if (intersect(ray, P, N, t, ids)) {
				Sphere s = spheres[ids];
				/* ---------- 1. Cas de la sphère MIROIR ------------- */
				if (s.miroir) { // La sphère est une sphère miroir, on enclanche le processus récursif
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
					if (s.isLum) {
						/* ---------- Si c'est la sphère lumineus, on effectue seulement une atténuation de distance ------------- */
						rayColor = s.albedo * I / (4 * M_PI * s.R * s.R);
					}
					else {
						/* ---------- Eclairage direct ------------- */
						Sphere sLum = spheres[0];
						Vector OP = P - sLum.O;
						OP.normalize();

						Vector NPrime = random_cos(OP);
						Vector PPrime = NPrime * (sLum.R + 0.001) + sLum.O;
						Vector wi = PPrime - P;
						double d2 = wi.norm2();
						wi.normalize();

						Ray rayIntersect(P + wi * 0.001, wi);
						Vector PIntersect, NIntersect;
						double tIntersect;
						int idsIntersect;
						bool intersection = intersect(rayIntersect, PIntersect, NIntersect, tIntersect, idsIntersect);

						/* ---------- On ne met pas d'intensité s'il y a un obstacle ------------- */
						if (intersection && (P - PIntersect).norm2() < d2 && idsIntersect != 0) {
							rayColor = Vector(0, 0, 0);
						}
						/* ---------- Sinon, on calcule l'intensité ------------- */
						else {
							double intensite = (I / (4 * M_PI * M_PI * d2) * dot(OP * (-1), wi) * dot(NPrime, wi * (-1)) / dot(OP, NPrime));
							rayColor = s.albedo * intensite;
						}

						/* ---------- Eclairage indirecte à l'aide d'une BRDF diffuse ------------- */
						Vector reflechi = random_cos(N);
						Ray ray_reflechi(P + N * 0.001, reflechi);
						rayColor = rayColor + s.albedo * get_color(P, N, ray_reflechi, nbRebonds - 1); // aussi / M_PI ?
					}
				}
			}
			else { // Sinon, mettre les pixels en noir
				rayColor = Vector(0, 0, 0);
			}
		}

		/* ---------- Si la lumière est ponctuelle ------------- */
		else {
			if (intersect(ray, P, N, t, ids)) {
				Sphere s = spheres[ids];
				/* ---------- 1. Cas de la sphère MIROIR ------------- */
				if (s.miroir) { // La sphère est une sphère miroir, on enclanche le processus récursif
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
						rayColor = s.albedo * intensite;
					}

					/* ---------- Contribution indirecte ------------- */
					nbRebonds--;
					Vector reflechi = random_cos(N);
					Ray ray_reflechi(P + N * 0.001, reflechi);
					rayColor = rayColor + s.albedo * get_color(P, N, ray_reflechi, nbRebonds); // aussi / M_PI ? 
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
	int precision = 2;
	int W = 256 * precision;
	int H = 256 * precision;
	double fov = 60 * M_PI / 180.;

	Vector rouge(1, 0, 0), vert(0, 1, 0), bleu(0, 0, 1), cyan(0, 1, 1), violet(1, 0, 1), jaune(1, 1, 0), blanc(1, 1, 1);

	Vector C(0, 0, 55);
	Vector L(-10, -20, 40);
	double I = 10000000000;

	std::vector<unsigned char> image(W * H * 3, 0);
	vector<Sphere> spheres;

	if (LUM_SPHERIQUE) {
		spheres.push_back(Sphere(Vector(-20, -40, 50), 10., blanc, false, false, 1, true));
	}

	/* ----- Ajouts de sphères à la scène ------------ */
	//spheres.push_back(Sphere(Vector(0, -3, 10), 12., blanc, false, false));
	spheres.push_back(Sphere(Vector(-13, 0, 3), 10., blanc, false, false));
	spheres.push_back(Sphere(Vector(0, -15, 5), 7., rouge, true, false));
	//spheres.push_back(Sphere(Vector(13, 0, 0), 10., blanc, false, true, 1.4));
	spheres.push_back(Sphere(Vector(13, 0, 10), 5., blanc, false, false));
	//spheres.push_back(Sphere(Vector(-20, 3, 25), 2., jaune, false, false));
	//spheres.push_back(Sphere(Vector(10, -10, -25), 6., cyan, false, false));
	spheres.push_back(Sphere(Vector(0, -1000, 0), 950., vert, false, false));
	spheres.push_back(Sphere(Vector(0, 1000, 0), 990., bleu, false, false));
	spheres.push_back(Sphere(Vector(0, 0, -5000), 4990., rouge, false, false));
	//spheres.push_back(Sphere(Vector(0, 0, 1000), 960., jaune, false, false));
	spheres.push_back(Sphere(Vector(1000, 0, 0), 965., jaune, false, false));
	spheres.push_back(Sphere(Vector(-1000, 0, 0), 965., cyan, false, false));
	Scene scene = Scene(spheres, C, L, I);

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

			int nbRays = 1000; // nombre de rayon de réflexion pour la lumière diffuse
			Vector meanRayColor(0, 0, 0);
			for (int k = 0; k < nbRays; k++) {
				/* --------- ANTIALIASING ----------- */
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
	stbi_write_png("lum_ponct_1000rays.png", W, H, 3, &image[0], 0);

	return 0;
}
