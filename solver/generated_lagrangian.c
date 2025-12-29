
#include <math.h>
#include <stddef.h>

typedef struct {
	float x;
	float y;
	float z;
} Vector3D;
void elastic_pendulum_step(Vector3D* q, Vector3D* dq, Vector3D* _dq, Vector3D* _ddq, float t, size_t N) {
    // Auto-generated Euler-Lagrange Equations using sympy.physics.mechanics
    // Constants have been collapsed into their values.
    float g = 9.81;
    float k = 50.0;
    float l0 = 1.0;
    float m = 1.0;
    _dq[0].x = dq[0].x;
    _ddq[0].x = (-k*q[0].x*(-l0 + sqrt(pow(q[0].x, 2) + pow(q[0].y, 2) + pow(q[0].z, 2)))/sqrt(pow(q[0].x, 2) + pow(q[0].y, 2) + pow(q[0].z, 2)) - k*(-l0 + sqrt(pow(-q[0].x + q[1].x, 2) + pow(-q[0].y + q[1].y, 2) + pow(-q[0].z + q[1].z, 2)))*(q[0].x - q[1].x)/sqrt(pow(-q[0].x + q[1].x, 2) + pow(-q[0].y + q[1].y, 2) + pow(-q[0].z + q[1].z, 2)))/m;
    _dq[0].y = dq[0].y;
    _ddq[0].y = (-k*q[0].y*(-l0 + sqrt(pow(q[0].x, 2) + pow(q[0].y, 2) + pow(q[0].z, 2)))/sqrt(pow(q[0].x, 2) + pow(q[0].y, 2) + pow(q[0].z, 2)) - k*(-l0 + sqrt(pow(-q[0].x + q[1].x, 2) + pow(-q[0].y + q[1].y, 2) + pow(-q[0].z + q[1].z, 2)))*(q[0].y - q[1].y)/sqrt(pow(-q[0].x + q[1].x, 2) + pow(-q[0].y + q[1].y, 2) + pow(-q[0].z + q[1].z, 2)))/m;
    _dq[0].z = dq[0].z;
    _ddq[0].z = (-g*m - k*q[0].z*(-l0 + sqrt(pow(q[0].x, 2) + pow(q[0].y, 2) + pow(q[0].z, 2)))/sqrt(pow(q[0].x, 2) + pow(q[0].y, 2) + pow(q[0].z, 2)) - k*(-l0 + sqrt(pow(-q[0].x + q[1].x, 2) + pow(-q[0].y + q[1].y, 2) + pow(-q[0].z + q[1].z, 2)))*(q[0].z - q[1].z)/sqrt(pow(-q[0].x + q[1].x, 2) + pow(-q[0].y + q[1].y, 2) + pow(-q[0].z + q[1].z, 2)))/m;
    _dq[1].x = dq[1].x;
    _ddq[1].x = (-k*(-l0 + sqrt(pow(-q[0].x + q[1].x, 2) + pow(-q[0].y + q[1].y, 2) + pow(-q[0].z + q[1].z, 2)))*(-q[0].x + q[1].x)/sqrt(pow(-q[0].x + q[1].x, 2) + pow(-q[0].y + q[1].y, 2) + pow(-q[0].z + q[1].z, 2)) - k*(-l0 + sqrt(pow(-q[1].x + q[2].x, 2) + pow(-q[1].y + q[2].y, 2) + pow(-q[1].z + q[2].z, 2)))*(q[1].x - q[2].x)/sqrt(pow(-q[1].x + q[2].x, 2) + pow(-q[1].y + q[2].y, 2) + pow(-q[1].z + q[2].z, 2)))/m;
    _dq[1].y = dq[1].y;
    _ddq[1].y = (-k*(-l0 + sqrt(pow(-q[0].x + q[1].x, 2) + pow(-q[0].y + q[1].y, 2) + pow(-q[0].z + q[1].z, 2)))*(-q[0].y + q[1].y)/sqrt(pow(-q[0].x + q[1].x, 2) + pow(-q[0].y + q[1].y, 2) + pow(-q[0].z + q[1].z, 2)) - k*(-l0 + sqrt(pow(-q[1].x + q[2].x, 2) + pow(-q[1].y + q[2].y, 2) + pow(-q[1].z + q[2].z, 2)))*(q[1].y - q[2].y)/sqrt(pow(-q[1].x + q[2].x, 2) + pow(-q[1].y + q[2].y, 2) + pow(-q[1].z + q[2].z, 2)))/m;
    _dq[1].z = dq[1].z;
    _ddq[1].z = (-g*m - k*(-l0 + sqrt(pow(-q[0].x + q[1].x, 2) + pow(-q[0].y + q[1].y, 2) + pow(-q[0].z + q[1].z, 2)))*(-q[0].z + q[1].z)/sqrt(pow(-q[0].x + q[1].x, 2) + pow(-q[0].y + q[1].y, 2) + pow(-q[0].z + q[1].z, 2)) - k*(-l0 + sqrt(pow(-q[1].x + q[2].x, 2) + pow(-q[1].y + q[2].y, 2) + pow(-q[1].z + q[2].z, 2)))*(q[1].z - q[2].z)/sqrt(pow(-q[1].x + q[2].x, 2) + pow(-q[1].y + q[2].y, 2) + pow(-q[1].z + q[2].z, 2)))/m;
    _dq[2].x = dq[2].x;
    _ddq[2].x = -k*(-l0 + sqrt(pow(-q[1].x + q[2].x, 2) + pow(-q[1].y + q[2].y, 2) + pow(-q[1].z + q[2].z, 2)))*(-q[1].x + q[2].x)/(m*sqrt(pow(-q[1].x + q[2].x, 2) + pow(-q[1].y + q[2].y, 2) + pow(-q[1].z + q[2].z, 2)));
    _dq[2].y = dq[2].y;
    _ddq[2].y = -k*(-l0 + sqrt(pow(-q[1].x + q[2].x, 2) + pow(-q[1].y + q[2].y, 2) + pow(-q[1].z + q[2].z, 2)))*(-q[1].y + q[2].y)/(m*sqrt(pow(-q[1].x + q[2].x, 2) + pow(-q[1].y + q[2].y, 2) + pow(-q[1].z + q[2].z, 2)));
    _dq[2].z = dq[2].z;
    _ddq[2].z = (-g*m - k*(-l0 + sqrt(pow(-q[1].x + q[2].x, 2) + pow(-q[1].y + q[2].y, 2) + pow(-q[1].z + q[2].z, 2)))*(-q[1].z + q[2].z)/sqrt(pow(-q[1].x + q[2].x, 2) + pow(-q[1].y + q[2].y, 2) + pow(-q[1].z + q[2].z, 2)))/m;
return;
}