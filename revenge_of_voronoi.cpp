//書いていて死にたくなった
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <set>
#include <functional>
#include <utility>
#include <map>
#include <climits>

using namespace std;

//#define DEBUG
//#define CHECK

#ifdef DEBUG
#define CHECK
#endif

#ifdef CHECK
static int DISTANCE[1024 * 1024];
static char RESULT[1024 * 1024];
#endif

static const int EMPTY = -1;
int TABLE[1024 * 1024];
int WIDTH, HEIGHT;
//static const int MAX_LETTER = CHAR_MAX;

struct DIC
{
	int to_index[CHAR_MAX + 1];
	vector<char> to_char;
	char strongest;
	DIC() : strongest(127){
		fill(to_index, to_index + 128, EMPTY);
	}
	void Insert(char x){
		to_char.push_back(x);
	}
	void MakeDic(){
		sort(to_char.begin(), to_char.end());
		to_char.erase(unique(to_char.begin(), to_char.end()), to_char.end());
		for (int i = 0; i < to_char.size(); ++i){
			to_index[to_char[i]] = i;
		}
	}
	int ToIndex(char x){
		return to_index[x];
	}
	int ToLetter(int index) const{
		return to_char[index];
	}
	int GetCount(){
		return to_char.size();
	}
};

int PositionToIndex(int x, int y)
{
	return x * HEIGHT+ y;
}

void IndexToPosition(int index, int& x, int& y)
{
	x = index / HEIGHT;
	y = index % HEIGHT;
}

int Distance(int x1, int y1, int x2, int y2)
{
	return abs(x1 - x2) + abs(y1 - y2);
}

//最終段階
//適合点の探索
map<pair<int, int>, map<int, set<int> > > EDGES;	//(文字インデックス,文字インデックス)->元場所インデックス->先場所インデックス
vector<pair<int, int> > ROUTE;
set<pair<int, int> > VISITED;

void RecRoute(int letter_index, vector<int>* graph)
{
	vector<int>& edge = graph[letter_index];
	for (vector<int>::iterator it = edge.begin(); it != edge.end(); ++it){
		const pair<int, int> p1 = make_pair(letter_index, *it);
		const pair<int, int> p2 = make_pair(*it, letter_index);
		if (VISITED.find(p1) != VISITED.end() || VISITED.find(p2) != VISITED.end()){
			continue;
		}
		VISITED.insert(p1);
		VISITED.insert(p2);

		ROUTE.push_back(p1);
		RecRoute(*it, graph);
	}
}

//再帰開始時にはMARK[0]にcharcters[0]を順に入れてから回す。
vector<int> MARK[CHAR_MAX + 1];
bool RecPair(int route_index)
{
#ifdef DEBUG
	for (int i = 0; i < route_index; ++i){
		int x1, y1, x2, y2;
		IndexToPosition(MARK[ROUTE[i].first].back(), x1, y1);
		IndexToPosition(MARK[ROUTE[i].second].back(), x2, y2);
		cout << "(" << x1 << "," << y1 << ")-(" << x2 << "," << y2 << ") ";
	}
	cout << endl;
#endif

	//適合する条件がある場合は終了
	if (route_index == ROUTE.size()){
		return true;
	}

	const int letter_index = ROUTE[route_index].first;
	const int position_index = MARK[letter_index].back();

	const int next_letter_index = ROUTE[route_index].second;
	const map<int, set<int> >& route = EDGES[ROUTE[route_index]];
	map<int, set<int> >::const_iterator it_route;
	if ((it_route = route.find(position_index)) == route.end()){
		return false;
	}
	const set<int>& edge = it_route->second;

	for (set<int>::const_iterator it = edge.begin(); it != edge.end(); ++it){
		const int next_position_index = *it;

		if (!MARK[next_letter_index].empty() && MARK[next_letter_index].back() != next_position_index){
			continue;
		}

		MARK[next_letter_index].push_back(next_position_index);

		if (RecPair(route_index + 1)){
			return true;
		}

		MARK[next_letter_index].pop_back();
	}

	return false;
}

void InsertNeighbor(map<pair<int, int>, vector<pair<int, int> > >& neighbors, int l1, int l2, int p1, int p2)
{
	if (l1 > l2){
		swap(l1, l2);
		swap(p1, p2);
	}
	neighbors[make_pair(l1, l2)].push_back(make_pair(p1, p2));
}

bool isGeneratedBoundary(int strongX, int strongY, int weakX, int weakY, int x1, int y1, int x2, int y2)
{
	return Distance(strongX, strongY, x1, y1) <= Distance(strongX, strongY, x2, y2) &&
		Distance(weakX, weakY, x1, y1) > Distance(weakX, weakY, x2, y2);
}

//領域同士の境界から点のペアの候補を探索する際、条件に合わない候補を削除する
struct RemovePosPairs : public unary_function<pair<int, int>, bool>
{
	int strongX, strongY;
	int weakX, weakY;
	RemovePosPairs(int strong, int weak){
		IndexToPosition(strong, strongX, strongY);
		IndexToPosition(weak, weakX, weakY);
	}
	bool operator()(const pair<int, int>& rh)const{
		int x1, y1, x2, y2;
		IndexToPosition(rh.first, x1, y1);
		IndexToPosition(rh.second, x2, y2);
		return !isGeneratedBoundary(strongX, strongY, weakX, weakY, x1, y1, x2, y2);
	}
};

int main()
{
	int testCase = 0;
	while (cin >> WIDTH >> HEIGHT && !(WIDTH == 0 && HEIGHT == 0)){
		//cout << "Case " << ++testCase << ":" << endl;

		//入力
		string line;
		getline(cin, line);

		DIC dic;
		for (int y = 0; y < HEIGHT; ++y){
			getline(cin, line);
			for (int x = 0; x < WIDTH; ++x){
				const char ch = line[x];
				TABLE[PositionToIndex(x, y)] = ch;
				dic.Insert(ch);
			}
		}
		dic.MakeDic();

		if (dic.GetCount() == 1){
			cout << line[0] << " 0 0" << endl << endl;
			continue;
		}

		for (int y = 0; y < HEIGHT; ++y){
			for (int x = 0; x < WIDTH; ++x){
				TABLE[PositionToIndex(x, y)] = dic.ToIndex(TABLE[PositionToIndex(x, y)]);
			}
		}

		const int number_of_letters = dic.GetCount();

		//境界情報検出
		//領域検出
		//隣接領域検出
		vector<int> letter_positions[CHAR_MAX + 1];	//文字インデックス -> 場所インデックス
		vector<int> graph[CHAR_MAX + 1];	//接続元文字インデックス -> 接続先文字インデックス
		map<pair<int, int>, vector<pair<int, int> > > neighbors;	//(文字インデックス,文字インデックス)->(場所インデックス,場所インデックス)

		for (int y = 0; y < HEIGHT; ++y){
			for (int x = 0; x < WIDTH; ++x){
				//領域
				const int position_index = PositionToIndex(x, y);
				const int letter_index = TABLE[position_index];
				letter_positions[letter_index].push_back(position_index);

				//右方向
				if (x != WIDTH - 1){
					const int position_index_right = PositionToIndex(x + 1, y);
					const int letter_index_right = TABLE[position_index_right];
					if (letter_index != letter_index_right){
						graph[letter_index].push_back(letter_index_right);
						graph[letter_index_right].push_back(letter_index);
						InsertNeighbor(neighbors, letter_index, letter_index_right, position_index, position_index_right);
					}
				}

				//下方向
				if (y != HEIGHT - 1){
					const int position_index_bottom = PositionToIndex(x, y + 1);
					const int letter_index_bottom = TABLE[position_index_bottom];
					if (letter_index != letter_index_bottom){
						graph[letter_index].push_back(letter_index_bottom);
						graph[letter_index_bottom].push_back(letter_index);
						InsertNeighbor(neighbors, letter_index, letter_index_bottom, position_index, position_index_bottom);
					}
				}
			}
		}

		//領域のペアを調べて点のペアの候補を作成する
		EDGES.clear();
		vector<int> wholeEdges[100 * 100];
		for (map<pair<int, int>, vector<pair<int, int> > >::iterator it = neighbors.begin(); it != neighbors.end(); ++it) {
			const int strongAreaIndex = it->first.first;
			const int weakAreaIndex = it->first.second;
			const vector<int>& positions1 = letter_positions[strongAreaIndex];
			const vector<int>& positions2 = letter_positions[weakAreaIndex];
			const vector<pair<int, int> >& boundaries = it->second;
			vector<pair<int, int> > pos_pairs;

			//一個目は候補の全列挙
			{
				int strong_x, strong_y, weak_x, weak_y;
				IndexToPosition(boundaries[0].first, strong_x, strong_y);
				IndexToPosition(boundaries[0].second, weak_x, weak_y);
				for (vector<int>::const_iterator it_position1 = positions1.begin(); it_position1 != positions1.end(); ++it_position1){
					int x1, y1;
					IndexToPosition(*it_position1, x1, y1);
					for (vector<int>::const_iterator it_posiiton2 = positions2.begin(); it_posiiton2 != positions2.end(); ++it_posiiton2){
						int x2, y2;
						IndexToPosition(*it_posiiton2, x2, y2);

						if (Distance(strong_x, strong_y, x1, y1) <= Distance(strong_x, strong_y, x2, y2) &&
							Distance(weak_x, weak_y, x1, y1) > Distance(weak_x, weak_y, x2, y2))
						{
							pos_pairs.push_back(make_pair(*it_position1, *it_posiiton2));
						}
					}
				}
			}

			//二個目以降は候補削り
			for (vector<pair<int, int> >::const_iterator it_boundaries = boundaries.begin() + 1; it_boundaries < boundaries.end(); ++it_boundaries){
				const int strong = it_boundaries->first;
				const int weak = it_boundaries->second;
				pos_pairs.erase(remove_if(pos_pairs.begin(), pos_pairs.end(), RemovePosPairs(strong, weak)), pos_pairs.end());
			}

			//候補の保存
			{
				map<int, set<int> >& v = EDGES[make_pair(strongAreaIndex, weakAreaIndex)];
				for (vector<pair<int, int> >::iterator it_pos_pair = pos_pairs.begin(); it_pos_pair < pos_pairs.end(); ++it_pos_pair){
					v[it_pos_pair->first].insert(it_pos_pair->second);
					wholeEdges[it_pos_pair->first].push_back(it_pos_pair->second);
				}
			}
			{
				map<int, set<int> >& v = EDGES[make_pair(weakAreaIndex, strongAreaIndex)];
				for (vector<pair<int, int> >::iterator it_pos_pair = pos_pairs.begin(); it_pos_pair < pos_pairs.end(); ++it_pos_pair){
					v[it_pos_pair->second].insert(it_pos_pair->first);
					wholeEdges[it_pos_pair->second].push_back(it_pos_pair->first);
				}
			}

			//const int strongAreaIndex = it->first.first;
			//const int weakAreaIndex = it->first.second;
			//const vector<int>& strongPositions = letter_positions[strongAreaIndex];
			//const vector<int>& weakPositions = letter_positions[weakAreaIndex];
			//const vector<pair<int, int> >& boundaries = it->second;

			//for (vector<int>::const_iterator strongPosition = strongPositions.begin(); strongPosition != strongPositions.end(); ++strongPosition){
			//	int x1, y1;
			//	IndexToPosition(*strongPosition, x1, y1);

			//	for (vector<int>::const_iterator weakPosition = weakPositions.begin(); weakPosition != weakPositions.end(); ++weakPosition){
			//		int x2, y2;
			//		IndexToPosition(*weakPosition, x2, y2);

			//		bool ok = true;
			//		for (vector<pair<int, int> >::const_iterator boundary = boundaries.begin(); ok && boundary != boundaries.end(); ++boundary){
			//			int strongX, strongY, weakX, weakY;
			//			IndexToPosition(boundary->first, strongX, strongY);
			//			IndexToPosition(boundary->second, weakX, weakY);

			//			ok = isGeneratedBoundary(strongX, strongY, weakX, weakY, x1, y1, x2, y2);
			//		}

			//		if (ok){
			//			EDGES[make_pair(strongAreaIndex, weakAreaIndex)][*strongPosition].insert(*weakPosition);
			//			EDGES[make_pair(weakAreaIndex, strongAreaIndex)][*weakPosition].insert(*strongPosition);
			//			wholeEdges[*strongPosition].push_back(*weakPosition);
			//			wholeEdges[*weakPosition].push_back(*strongPosition);
			//		}
			//	}
			//}
		}

		//隣接領域探索順番の決定
		ROUTE.clear();
		for (int i = 0; i < dic.GetCount(); ++i){
			sort(graph[i].begin(), graph[i].end());
			graph[i].erase(unique(graph[i].begin(), graph[i].end()), graph[i].end());
		}
		VISITED.clear();
		RecRoute(0, graph);

		//枝狩り(両端がつながっていないエッジを削除)
		bool cleanuped = false;
		while (!cleanuped){
			cleanuped = true;

#ifdef DEBUG
			cout << "*";
#endif

			vector<pair<int, int> > areaPairsToRemove;
			//(文字インデックス,文字インデックス)->元場所インデックス->先場所インデックス
			for (map<pair<int, int>, map<int, set<int> > >::iterator edgeMap = EDGES.begin(); edgeMap != EDGES.end(); ++edgeMap){
				vector<pair<int, int> > edgesToRemove;
				const map<int, set<int> >& edges = edgeMap->second;
				for (map<int, set<int> >::const_iterator edge = edges.begin(); edge != edges.end(); ++edge){
					const int srcAreaIndex = TABLE[edge->first];

					//枝の先に接続先が有るかどうかチェック
					const int srcPositionIndex = edge->first;
					const set<int>& destPositionIndexes = edge->second;
					for (set<int>::const_iterator destPositionIndex = destPositionIndexes.begin(); destPositionIndex != destPositionIndexes.end(); ++destPositionIndex){
						const vector<int>& wholeDests = wholeEdges[*destPositionIndex];
						set<int> destAreas;
						for (vector<int>::const_iterator wholeDest = wholeDests.begin(); wholeDest != wholeDests.end(); ++wholeDest){
							destAreas.insert(TABLE[*wholeDest]);
						}

						if (graph[TABLE[*destPositionIndex]].size() == destAreas.size()){
							continue;
						}

						edgesToRemove.push_back(make_pair(srcPositionIndex, *destPositionIndex));
					}
				}

				for (vector<pair<int, int> >::iterator edgeToRemove = edgesToRemove.begin(); edgeToRemove != edgesToRemove.end(); ++edgeToRemove){
					edgeMap->second[edgeToRemove->first].erase(edgeToRemove->second);
					if (edgeMap->second[edgeToRemove->first].empty()){
						edgeMap->second.erase(edgeToRemove->first);
					}
				}

				if (edgeMap->second.empty()){
					areaPairsToRemove.push_back(edgeMap->first);
				}
			}

			for (vector<pair<int, int> >::iterator areaPairToRemove = areaPairsToRemove.begin(); areaPairToRemove != areaPairsToRemove.end(); ++areaPairToRemove){
				EDGES.erase(*areaPairToRemove);
			}
		}

		//ペアAND
		for (int i = 0; i < number_of_letters; ++i){
			MARK[i].clear();
		}
		for (vector<int>::iterator it = letter_positions[0].begin(); it != letter_positions[0].end(); ++it){
			MARK[0].push_back(*it);
			if (RecPair(0)){
				break;
			}
			MARK[0].pop_back();
		}

		for (int i = 0; i < number_of_letters; ++i){
			int x, y;
			IndexToPosition(MARK[i].back(), x, y);
			cout << (char)dic.ToLetter(i) << " " << x << " " << y << endl;
		}

		//cout << endl;

#ifdef CHECK
		fill(DISTANCE, DISTANCE + sizeof(DISTANCE) / sizeof(DISTANCE[0]), INT_MAX);
		for (int y = 0; y < HEIGHT; ++y){
			for (int x = 0; x < WIDTH; ++x){
				char ch = 0;
				for (int i = 0; i < number_of_letters; ++i){
					char c = dic.ToLetter(i);
					int xx, yy;
					IndexToPosition(MARK[i].back(), xx, yy);
					int d = Distance(xx, yy, x, y);
					if (DISTANCE[PositionToIndex(x, y)] > d){
						DISTANCE[PositionToIndex(x, y)] = d;
						ch = c;
					}
				}
				RESULT[PositionToIndex(x, y)] = ch;
				//cout << ch;
			}
			//cout << endl;
		}

		for (int y = 0; y < HEIGHT; ++y){
			for (int x = 0; x < WIDTH; ++x){
				int index = PositionToIndex(x, y);
				if (dic.ToLetter(TABLE[index]) != RESULT[index]){
					cout << "error (" << x << "," << y << ") " << (char)dic.ToLetter(TABLE[index]) << " -> " << RESULT[index] << endl;
				}
			}
		}
#endif
	}
}
