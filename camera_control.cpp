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
static const double EPS = 1e-8;
static const double PI = 4.0 * atan(1.0);
static const double PI2 = 8.0 * atan(1.0);
typedef long long ll;

#define ALL(c) (c).begin(), (c).end()
#define CLEAR(v) memset(v,0,sizeof(v))
#define MP(a,b) make_pair((a),(b))
#define REP(i,n) for(int i=0;i<(int)n;++i)

typedef complex<double> P;

static const double kMaxTime = 1e4;

// 以下ユーティリティ
struct EventPoint {
	vector<int> vocalIndexes;
	P representativePoint;
	double time;
	// 他のイベントポイントにマージされた場合
	bool merged;

	EventPoint() : merged(false) { }
};

struct VocalPosition {
	P position;
	double timeBegin;
	VocalPosition(const P& position, double timeBegin) : position(position), timeBegin(timeBegin) { }
};

double dot(const P& a, const P& b) {
	return a.real() * b.real() + a.imag() * b.imag();
}

double cross(const P& a, const P& b) {
	return a.real() * b.imag() - a.imag() * b.real();
}

bool isSameAngle(const P& a, const P& b) {
	// EPSの値に注意
	return abs(a / abs(a) - b / abs(b)) < EPS;
}

// 相似形で移動しているかどうか
// 四点が同一方向に存在する場合は含めない
bool isHomothetic(const P& a, const P& b, const P& c, const P& d) {
	return !isSameAngle(a, b) && !isSameAngle(c, d) &&
		isSameAngle(a, c) && isSameAngle(b, d) &&
		abs(cross(a - b, c - d)) < EPS;
}

// 四点が同一方向に存在するか？
bool isSameAngle(const P& a, const P& b, const P& c, const P& d) {
	return isSameAngle(a, b) && isSameAngle(a, c) && isSameAngle(a, d);
}

// ボーカル
struct Vocal {
	int index;
	// 移動スケジュールリスト
	// 「この時間にはここにいる」
	// 最後の要素は番兵として無限時間における位置を置く
	vector<VocalPosition> vocalPositions;
	// パートスケジュールリスト
	vector<pair<double, double> > vocalParts;
	// イベント番号
	vector<int> eventPointIndexes;
	// TODO: ノードリスト

	void inputPositions() {
		int M;
		cin >> M;
		REP(m, M) {
			double x, y, t;
			cin >> x >> y >> t;
			vocalPositions.push_back(VocalPosition(P(x, y), t));
		}
		vocalPositions.push_back(VocalPosition(vocalPositions.back().position, kMaxTime));
	}
	P getPosition(double time) const {
		// 整数用二分探索
		// TODO: timeがkMaxTimeだった場合に処理が正しく行われるかどうか
		int left = 0;
		int right = vocalPositions.size();
		while (left + 1 < right){
			const int middle = (left + right) / 2;
			if (vocalPositions[middle].timeBegin < time) {
				left = middle;
			} else {
				right = middle;
			}
		}

		const int index = min(left, (int)vocalPositions.size() - 2);

		// 比の計算
		const double t0 = vocalPositions[index].timeBegin;
		const double t1 = vocalPositions[index + 1].timeBegin;
		const double r = (time - t0) / (t1 - t0);

		// 位置の計算
		const P& p0 = vocalPositions[index].position;
		const P& p1 = vocalPositions[index + 1].position;
		const P p = p0 * (1.0 - r) + p1 * r;
		return p;
	}
	void inputVocalParts() {
		int L;
		cin >> L;
		REP(l, L) {
			pair<double, double> part;
			cin >> part.first >> part.second;
			vocalParts.push_back(part);
		}

		// TODO: ソート＆区間併合を必要とさせるかどうか？
		//       本質的ではないので予めソート+区間併合なしで入力させたい
	}
	void translate(const P& p) {
		REP(i, vocalPositions.size()) {
			vocalPositions[i].position -= p;
		}
	}
	double calculateCross(const Vocal& other, double time) {
		const P pointMy = getPosition(time);
		const P pointOther = other.getPosition(time);
		return abs(cross(pointMy, pointOther));
	}
	void registerEventPoints(const Vocal& other, vector<EventPoint>& eventPoints, double time) {
		const P pointMy = getPosition(time);
		const P pointOther = other.getPosition(time);

		if (dot(pointMy, pointOther) < 0) {
			// 原点を挟んで反対側にある
			return;
		}

		{
			const double r = abs(pointMy / abs(pointMy) - pointOther / abs(pointOther));
			assert(abs(r) < EPS);
		}

		eventPoints.push_back(EventPoint());
		EventPoint& eventPoint = eventPoints.back();
		eventPoint.representativePoint = pointMy;
		eventPoint.time = time;
		eventPoint.vocalIndexes.push_back(index);
		eventPoint.vocalIndexes.push_back(other.index);
	}
	void extractEventPoints(Vocal& other, vector<EventPoint>& eventPoints) {
		// タイミングを列挙する
		vector<double> times;
		REP(i, vocalPositions.size()) {
			times.push_back(vocalPositions[i].timeBegin);
		}
		REP(i, other.vocalPositions.size()) {
			times.push_back(other.vocalPositions[i].timeBegin);
		}

		sort(ALL(times));

		times.erase(unique(ALL(times)), times.end());

		// スケジュールの区間内で視線が重なるタイミングをイベントポイントに追加する
		// 視線が重なるタイミングはcross(a,b)==0かつdot(a,b)>0のとき
		REP(i, times.size() - 1) {
			const double timeBegin = times[i];
			const double timeEnd = times[i + 1];
			const P aBegin = getPosition(timeBegin);
			const P aEnd = getPosition(timeEnd);
			const P bBegin = other.getPosition(timeBegin);
			const P bEnd = other.getPosition(timeEnd);

			// TODO: カメラを挟んで反対側にあり、かつ動かない場合

			// TODO: A=B=C=0のパターン
			// 平行移動している場合、又は四点とも同じ方向にある場合
			if (isHomothetic(aBegin, aEnd, bBegin, bEnd) ||
				isSameAngle(aBegin, aEnd, bBegin, bEnd))
			{
				continue;
			}

			// Maximaで計算したが大丈夫か？
			// 大丈夫だ、問題ない。
			const double a = aBegin.real();
			const double b = aBegin.imag();
			const double c = aEnd.real();
			const double d = aEnd.imag();
			const double e = bBegin.real();
			const double f = bBegin.imag();
			const double g = bEnd.real();
			const double h = bEnd.imag();

			const double A = (c-a)*h+(b-d)*g+(a-c)*f+(d-b)*e;
			const double B = a*h-b*g+(c-2*a)*f+(2*b-d)*e;
			const double C = a*f-b*e;

			// 縮退して一次方程式となる場合
			if (abs(A) < EPS) {

				// さらに縮退する場合
				if (abs(B) < EPS) {
					// カメラを挟んで反対側におり、かつ動かない場合C=0?
					continue;
				}

				const double r = -C / B;
				if (-EPS < r && r < 1.0 + EPS) {
					registerEventPoints(other, eventPoints, timeBegin * (1.0 - r) + timeEnd * r);
				}
				continue;
			}

			double D = B * B - 4 * A * C;
			if (abs(D) < EPS) {
				const double r = -B / (2.0 * A);
				if (-EPS < r && r < 1.0 + EPS) {
					registerEventPoints(other, eventPoints, timeBegin * (1.0 - r) + timeEnd * r);
				}
			} else if (D > 0.0) {
				double r0 = (-B + sqrt(D)) / (2.0 * A);
				if (-EPS < r0 && r0 < 1.0 + EPS) {
					registerEventPoints(other, eventPoints, timeBegin * (1.0 - r0) + timeEnd * r0);
				}

				double r1 = (-B - sqrt(D)) / (2.0 * A);
				if (-EPS < r1 && r1 < 1.0 + EPS) {
					registerEventPoints(other, eventPoints, timeBegin * (1.0 - r1) + timeEnd * r1);
				}
			}
		}
	}
	void extractVocalParts(double timeBegin, double timeEnd, vector<pair<double, double> >& extractedVocalParts) const {
		REP(i, vocalParts.size()) {
			if (vocalParts[i].second < timeBegin) {
				continue;
			}

			if (timeEnd < vocalParts[i].first) {
				continue;
			}

			const double begin = max(timeBegin, vocalParts[i].first);
			const double end = min(timeEnd, vocalParts[i].second);
			extractedVocalParts.push_back(MP(begin, end));
		}
	}
};

bool compareEventPointByTime(const EventPoint& lh, const EventPoint& rh) {
	return lh.time < rh.time;
}

bool isMerged(const EventPoint& ep) {
	return ep.merged;
}

void mergeEventPoints(vector<EventPoint>& eventPoints) {
	const int n = eventPoints.size();

	for (int i = 0; i < n; ++i) {
		EventPoint& ep0 = eventPoints[i];

		// 既に他のイベントポイントに併合されていた場合はスキップ
		if (ep0.merged) {
			continue;
		}

		for (int j = i + 1; j < n && eventPoints[j].time < ep0.time + EPS; ++j) {
			EventPoint& ep1 = eventPoints[j];

			// 既に他のイベントポイントに併合されていた場合はスキップ
			if (ep1.merged) {
				continue;
			}

			if (!isSameAngle(ep0.representativePoint, ep1.representativePoint)) {
				continue;
			}

			assert(abs(ep0.time - ep1.time) < EPS);
			ep1.merged = true;
			ep0.vocalIndexes.insert(ep0.vocalIndexes.end(), ALL(ep1.vocalIndexes));
		}
	}

	REP(i, n) {
		EventPoint& ep = eventPoints[i];
		sort(ALL(ep.vocalIndexes));
		ep.vocalIndexes.erase(unique(ALL(ep.vocalIndexes)), ep.vocalIndexes.end());
	}

	eventPoints.erase(remove_if(ALL(eventPoints), isMerged), eventPoints.end());
}

struct EventEdge {
	int src, dst, vocalIndex;
	double timeBegin, timeEnd;
	bool merged;
};

bool compareEventEdgeByTime(const EventEdge& lh, const EventEdge& rh) {
	return lh.timeBegin != rh.timeBegin ? lh.timeBegin < rh.timeBegin : lh.timeEnd < rh.timeEnd;
}

// Spaghetti Source - 各種アルゴリズムの C++ による実装 http://www.prefield.com/algorithm/
typedef double Weight;
struct Edge {
	int src, dst;
	Weight weight;
	Edge(int src, int dst, Weight weight) :
	src(src), dst(dst), weight(weight) { }
};
bool operator < (const Edge &e, const Edge &f) {
	return e.weight != f.weight ? e.weight > f.weight : // !!INVERSE!!
		e.src != f.src ? e.src < f.src : e.dst < f.dst;
}
typedef vector<Edge> Edges;
typedef vector<Edges> Graph;

bool visit(const Graph &g, int v, vector<int> &order, vector<int> &color) {
	color[v] = 1;
	for (Edges::const_iterator e = g[v].begin(), eEnd = g[v].end(); e != eEnd; ++e) {
		if (color[e->dst] == 2) continue;
		if (color[e->dst] == 1) return false;
		if (!visit(g, e->dst, order, color)) return false;
	}
	order.push_back(v); color[v] = 2;
	return true;
}
bool topologicalSort(const Graph &g, vector<int> &order) {
	int n = g.size();
	vector<int> color(n);
	REP(u, n) if (!color[u] && !visit(g, u, order, color))
		return false;
	reverse(ALL(order));
	return true;
}

bool isMergedVocalPart(const pair<double, double>& part) {
	return part.first == 0 && part.second == 0;
}

void mergeVocalParts(vector<pair<double, double> >& parts) {
	sort(ALL(parts));
	int lastIndex = 0;
	for (int i = 1; i < (int)parts.size(); ++i) {
		if (parts[i].first < parts[lastIndex].second) {
			parts[lastIndex].second = max(parts[lastIndex].second, parts[i].second);
			parts[i].first = 0.0;
			parts[i].second = 0.0;
		} else {
			lastIndex = i;
		}
	}

	parts.erase(remove_if(ALL(parts), isMergedVocalPart), parts.end());
}

int main() {
	for (int N; cin >> N && N; ) {
		vector<Vocal> vocals(N);

		double cameraX, cameraY;
		cin >> cameraX >> cameraY;

		// 入力
		REP(n, N) {
			vocals[n].inputPositions();
			vocals[n].inputVocalParts();
			vocals[n].index = n;
		}

		// カメラを原点に移動する
		const P cameraPosition(cameraX, cameraY);
		REP(n, N) {
			vocals[n].translate(cameraPosition);
		}

		set<double> eventTimings;
		eventTimings.insert(0);
		eventTimings.insert(kMaxTime);
		REP(n, N) {
			const vector<pair<double, double> >& vocalParts = vocals[n].vocalParts;
			REP(i, vocalParts.size()) {
				eventTimings.insert(vocalParts[i].first);
				eventTimings.insert(vocalParts[i].second);
			}
		}

		// 2人ずつ視線が交差するポイントを検出し、イベントポイントに加える
		vector<EventPoint> eventPoints;
		for (int i = 0; i < N; ++i) {
			Vocal& vocal = vocals[i];

			for (set<double>::iterator it = eventTimings.begin(), itEnd = eventTimings.end(); it != itEnd; ++it) {
				vocal.registerEventPoints(vocal, eventPoints, *it);
			}

			for (int j = i + 1; j < N; ++j) {
				vocals[j].extractEventPoints(vocal, eventPoints);
			}
		}

		// 時間順にソート
		sort(ALL(eventPoints), compareEventPointByTime);

		// 同じイベントポイントを持つボーカルの集合を一つのノードにまとめる
		mergeEventPoints(eventPoints);

		// エッジを列挙する
		REP(eventPointIndex, eventPoints.size()) {
			const EventPoint& ep = eventPoints[eventPointIndex];
			REP(j, ep.vocalIndexes.size()) {
				const int vocalIndex = ep.vocalIndexes[j];
				vocals[vocalIndex].eventPointIndexes.push_back(eventPointIndex);
			}
		}

		vector<EventEdge> eventEdges;
		REP(vocalIndex, vocals.size()) {
			const Vocal& vocal = vocals[vocalIndex];
			REP(i, vocal.eventPointIndexes.size() - 1) {
				EventEdge eventEdge;
				eventEdge.src = vocal.eventPointIndexes[i];
				eventEdge.dst = vocal.eventPointIndexes[i + 1];
				eventEdge.timeBegin = eventPoints[eventEdge.src].time;
				eventEdge.timeEnd = eventPoints[eventEdge.dst].time;
				eventEdge.vocalIndex = vocalIndex;
				eventEdge.merged = false;
				eventEdges.push_back(eventEdge);
			}
		}

		// DAGを生成する
		sort(ALL(eventEdges), compareEventEdgeByTime);

		Graph graph(eventPoints.size());

		for (int i = 0; i < (int)eventEdges.size(); ++i) {
			const EventEdge& eventEdge = eventEdges[i];

			// 既に他のイベントポイント間に併合されていた場合はスキップ
			if (eventEdge.merged) {
				continue;
			}

			const double timeBegin = eventEdge.timeBegin;
			const double timeEnd = eventEdge.timeEnd;
			vector<pair<double, double> > vocalParts;
			const Vocal& vocal = vocals[eventEdge.vocalIndex];
			vocal.extractVocalParts(timeBegin, timeEnd, vocalParts);
			const P positionBegin = vocal.getPosition(timeBegin);
			const P positionEnd = vocal.getPosition(timeEnd);

			for (int j = i + 1; j < (int)eventEdges.size() && eventEdges[j].timeBegin < timeBegin + EPS; ++j) {
				// 既に他のイベントポイント間に併合されていた場合はスキップ
				if (eventEdges[j].merged) {
					continue;
				}

				// 開始時刻と終了時刻が等しいかどうか
				if (abs(eventEdge.timeEnd - eventEdges[j].timeEnd) > EPS) {
					continue;
				}

				// 相似に移動する場合はボーカルパートを併合する
				const Vocal& vocalOther = vocals[eventEdges[j].vocalIndex];
				const P positionOtherBegin = vocalOther.getPosition(timeBegin);
				const P positionOtherEnd = vocalOther.getPosition(timeEnd);

				if (isHomothetic(positionBegin, positionEnd, positionOtherBegin, positionOtherEnd) ||
					isSameAngle(positionBegin, positionEnd, positionOtherBegin, positionOtherEnd))
				{
					vocalOther.extractVocalParts(timeBegin, timeEnd, vocalParts);
					eventEdges[j].merged = true;
				}
			}

			// イベントノード間のコスト求める
			mergeVocalParts(vocalParts);
			double timeSum = 0.0;
			REP(i, vocalParts.size()) {
				timeSum += vocalParts[i].second - vocalParts[i].first;
			}

			graph[eventEdge.src].push_back(Edge(eventEdge.src, eventEdge.dst, timeSum));
		}

		// DAGの最長パスを求める
		vector<int> order;
		topologicalSort(graph, order);

		vector<double> dp(eventPoints.size());
		REP(i, order.size()) {
			const int src = order[i];
			REP(j, graph[src].size()) {
				const int dst = graph[src][j].dst;
				const double weight = graph[src][j].weight;
				if (dp[dst] < dp[src] + weight) {
					dp[dst] = dp[src] + weight;
				}
			}
		}

		const double answer = *max_element(ALL(dp));
		printf("%.20f\n", answer);
	}
}
