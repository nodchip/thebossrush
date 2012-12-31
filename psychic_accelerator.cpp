#include <bitset>
#include <deque>
#include <list>
#include <map>
#include <queue>
#include <set>
#include <stack>
#include <vector>
#include <algorithm>
#include <functional>
#include <iterator>
#include <locale>
#include <memory>
#include <stdexcept>
#include <utility>
#include <string>
#include <fstream>
#include <ios>
#include <iostream>
#include <iosfwd>
#include <iomanip>
#include <istream>
#include <ostream>
#include <sstream>
#include <streambuf>
#include <complex>
#include <numeric>
#include <valarray>
#include <exception>
#include <limits>
#include <new>
#include <typeinfo>
#include <cassert>
#include <cctype>
#include <cerrno>
#include <cfloat>
#include <climits>
#include <cmath>
#include <csetjmp>
#include <csignal>
#include <cstdlib>
#include <cstddef>
#include <cstdarg>
#include <ctime>
#include <cstdio>
#include <cstring>
#include <cwchar>
#include <cwctype>
using namespace std;

static const double PI = atan(1.0) * 4.0;
const double EPS = 1e-8;
const double INF = 1e12;
typedef complex<double> P;
namespace std {
  bool operator < (const P& a, const P& b) {
    return real(a) != real(b) ? real(a) < real(b) : imag(a) < imag(b);
  }
}
double cross(const P& a, const P& b) {
  return imag(conj(a)*b);
}
double dot(const P& a, const P& b) {
  return real(conj(a)*b);
}

struct L : public vector<P> {
  L(const P &a, const P &b) {
    push_back(a); push_back(b);
  }
};

int ccw(P a, P b, P c) {
  b -= a; c -= a;
  if (cross(b, c) > 0)   return +1;       // counter clockwise
  if (cross(b, c) < 0)   return -1;       // clockwise
  if (dot(b, c) < 0)     return +2;       // c--a--b on line
  if (norm(b) < norm(c)) return -2;       // a--b--c on line
  return 0;
}

bool intersectLL(const L &l, const L &m) {
  return abs(cross(l[1]-l[0], m[1]-m[0])) > EPS || // non-parallel
         abs(cross(l[1]-l[0], m[0]-l[0])) < EPS;   // same line
}
bool intersectLS(const L &l, const L &s) {
  return cross(l[1]-l[0], s[0]-l[0])*       // s[0] is left of l
         cross(l[1]-l[0], s[1]-l[0]) < EPS; // s[1] is right of l
}
bool intersectLP(const L &l, const P &p) {
  return abs(cross(l[1]-p, l[0]-p)) < EPS;
}
bool intersectSS(const L &s, const L &t) {
  return ccw(s[0],s[1],t[0])*ccw(s[0],s[1],t[1]) <= 0 &&
         ccw(t[0],t[1],s[0])*ccw(t[0],t[1],s[1]) <= 0;
}
bool intersectSP(const L &s, const P &p) {
  return abs(s[0]-p)+abs(s[1]-p)-abs(s[1]-s[0]) < EPS; // triangle inequality
}

P projection(const L &l, const P &p) {
  double t = dot(p-l[0], l[0]-l[1]) / norm(l[0]-l[1]);
  return l[0] + t*(l[0]-l[1]);
}
P reflection(const L &l, const P &p) {
  return p + 2.0 * (projection(l, p) - p);
}
double distanceLP(const L &l, const P &p) {
  return abs(p - projection(l, p));
}
double distanceLL(const L &l, const L &m) {
  return intersectLL(l, m) ? 0 : distanceLP(l, m[0]);
}
double distanceLS(const L &l, const L &s) {
  if (intersectLS(l, s)) return 0;
  return min(distanceLP(l, s[0]), distanceLP(l, s[1]));
}
double distanceSP(const L &s, const P &p) {
  const P r = projection(s, p);
  if (intersectSP(s, r)) return abs(r - p);
  return min(abs(s[0] - p), abs(s[1] - p));
}
double distanceSS(const L &s, const L &t) {
  if (intersectSS(s, t)) return 0;
  return min(min(distanceSP(s, t[0]), distanceSP(s, t[1])),
             min(distanceSP(t, s[0]), distanceSP(t, s[1])));
}
P crosspoint(const L &l, const L &m) {
  double A = cross(l[1] - l[0], m[1] - m[0]);
  double B = cross(l[1] - l[0], l[1] - m[0]);
  if (abs(A) < EPS && abs(B) < EPS) return m[0]; // same line
  if (abs(A) < EPS) assert(false); // !!!PRECONDITION NOT SATISFIED!!!
  return m[0] + B / A * (m[1] - m[0]);
}

P calculateCenter(const P& a, const P& b, const P& c, const P& d) {
	if (!intersectLL(L(a, b), L(c, d))) {
		return (b + c) * 0.5;
	}

	const P p = b;
	const P q = b + (a - b) * polar(1.0, PI * 0.5);
	const P r = c;
	const P s = c + (d - c) * polar(1.0, PI * 0.5);
	return crosspoint(L(p, q), L(r, s));
}

int main()
{
	int N;
	double a;
	while (scanf("%d%lf", &N, &a) == 2 && N){
		vector<L> segments;
		vector<double> s;
		for (int i = 0; i < N; ++i){
			double x0, y0, x1, y1;
			scanf("%lf%lf%lf%lf", &x0, &y0, &x1, &y1);
			const P a(x0, y0);
			const P b(x1, y1);
			segments.push_back(L(a, b));
			s.push_back(abs(a - b));
		}

		vector<double> curveVelocities;
		curveVelocities.push_back(0.0);

		//加速度より各カーブの限界速度Viを計算
		vector<double> curveLengths;
		vector<double> curveRadius;
		vector<double> curveTheta;
		for (int i = 0; i < N - 1; ++i){
			const P& p0 = segments[i][0];
			const P& p1 = segments[i][1];
			const P& p2 = segments[i + 1][0];
			const P& p3 = segments[i + 1][1];
			P center = calculateCenter(p0, p1, p2, p3);
			const double r = abs(p1 - center);
			curveRadius.push_back(r);
			double d = dot(p0 - p1, p2 - p3) / abs(p0 - p1) / abs(p3 - p2);
			d = max(-1.0, d);
			d = min(1.0, d);
			double theta = acos(d);
			
			// 200 * 200あれば届く
			if (intersectSS(L(p1, p1 + (p0 - p1) * 1e5), L(p2, p2 + (p3 - p2) * 1e5))) {
				theta = PI * 2.0 - theta;
			}

			curveTheta.push_back(theta);
			curveVelocities.push_back(sqrt(a * r));
			curveLengths.push_back(r * theta);
		}

		curveVelocities.push_back(0.0);

		//(最後の点から最初の点に向かって)各カーブについて次のカーブに入る速度V_{i+1}から
		//現在のカーブを出るときの限界速度を計算し、V_{i}を①②のminに更新
		for (int i = N - 1; i >= 1; --i){
			curveVelocities[i] = min(curveVelocities[i], sqrt(curveVelocities[i + 1] * curveVelocities[i + 1] + 2 * a * s[i]));
		}

		//(最初の点から最後の点に向かって)各カーブについて前のカーブから出たときの速度V_{i-1}から
		//現在のカーブに入る事ができる限界速度を計算し、V_{i}を①②③のminに更新
		for (int i = 1; i < N; ++i){
			curveVelocities[i] = min(curveVelocities[i], sqrt(curveVelocities[i - 1] * curveVelocities[i - 1] + 2 * a * s[i - 1]));
		}

		double answer = 0;
		for (int i = 0; i < N - 1; ++i){
			answer += curveLengths[i] / curveVelocities[i + 1];
		}

		for (int i = 0; i < N; ++i){
			const double v0 = curveVelocities[i];
			const double v1 = curveVelocities[i + 1];
			const double x = (v1 * v1 - v0 * v0 + 2 * a * s[i]) / (4 * a);
			const double v = sqrt(v0 * v0 + 2 * a * x);
			const double t0 = (sqrt(v0 * v0 + 2 * a * x) - v0) / a;
			const double t1 = (v - sqrt(max(0.0, v * v - 2.0 * a * (s[i] - x)))) / a;

			//fprintf(stderr, "v0:%f v1:%f x:%f v:%f t0:%f t1:%f\n", v0, v1, x, v, t0, t1);
			//if (i < N - 1) {
			//	fprintf(stderr, "\tr:%f v:%f theta:%f\n", curveRadius[i], curveVelocities[i], curveTheta[i]);
			//}
			
			answer += t0;
			answer += t1;
		}

		printf("%.10lf\n", answer);
	}
}
