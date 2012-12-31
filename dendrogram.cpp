#include <vector>
#include <algorithm>
#include <functional>
#include <climits>
#include <algorithm>
#include <functional>
#include <set>
#include <cstdio>

using namespace std;

static const int ROOT = -1;
static const int EMPTY = -1;

struct BRANCH
{
	int y;
	vector<int> after;
	bool operator<(const BRANCH& rh) const{
		return y < rh.y;
	}
};

struct NODE
{
	vector<int> children;
	set<int> childrenTemp;
	int parent;
	int minLeaf;	//自分が持っている子供たちの中で最も番号の小さい葉
};

struct NodeSorter : public binary_function<int, int, bool>
{
	const vector<NODE>& graph;
	NodeSorter(const vector<NODE>& graph) : graph(graph){}
	bool operator()(int lh, int rh) const{
		return graph[lh].minLeaf > graph[rh].minLeaf;
	}
};

int main()
{
	for (int n, m, q; scanf("%d%d%d", &n, &m, &q) != EOF && n;){
		vector<BRANCH> branches;
		for (int i = 0; i < m; ++i){
			BRANCH branch;
			int numberOfChildren;
			scanf("%d%d", &branch.y, &numberOfChildren);
			for (int j = 0; j < numberOfChildren; ++j){
				int child;
				scanf("%d", &child);
				--child;
				branch.after.push_back(child);
			}
			branches.push_back(branch);
		}

		sort(branches.begin(), branches.end());

		//グラフ構築
		vector<NODE> graph(n);	//0-(n-1):葉 n-(2n-1): 分岐点
		for (int i = 0; i < n; ++i){
			graph[i].minLeaf = i;
			graph[i].parent = ROOT;
		}

		//葉と分岐をつなぐ
		for (vector<BRANCH>::iterator it = branches.begin(); it != branches.end(); ++it){
			const int parent = graph.size();
			NODE node;
			node.childrenTemp.insert(it->after.begin(), it->after.end());
			node.minLeaf = EMPTY;
			node.parent = ROOT;
			for (vector<int>::iterator jt = it->after.begin(); jt != it->after.end(); ++jt){
				if (graph[*jt].parent != ROOT){
					graph[graph[*jt].parent].childrenTemp.erase(*jt);
					graph[graph[*jt].parent].childrenTemp.insert(parent);
					node.parent = graph[*jt].parent;
				}
				graph[*jt].parent = parent;
			}

			graph.push_back(node);
		}

		//minLeafを書き換える
		for (int start = 0; start < n; ++start){
			int current = start;
			do {
				graph[current].minLeaf = start;
				current = graph[current].parent;
			} while (current != ROOT && graph[current].minLeaf == EMPTY);
		}

		//枝の順番を並び変える
		for (int i = n; i < graph.size(); ++i){
			NODE& node = graph[i];
			node.children.insert(node.children.end(), node.childrenTemp.begin(), node.childrenTemp.end());
			sort(node.children.begin(), node.children.end(), NodeSorter(graph));
		}

		//DFS開始
		vector<int> answer;
		vector<bool> visited(graph.size());
		for (int start = 0; start < n; ++start){
			if (visited[start]){
				continue;
			}

			//根を探す
			int root = start;
			while (graph[root].parent != ROOT){
				root = graph[root].parent;
			}

			//トポロジカルソート
			vector<int> stk;
			stk.push_back(root);
			while (!stk.empty()){
				const int currentNode = stk.back();
				stk.pop_back();
				if (visited[currentNode]){
					continue;
				}
				visited[currentNode] = true;

				if (currentNode < n){
					answer.push_back(currentNode);
				}

				stk.insert(stk.end(), graph[currentNode].children.begin(), graph[currentNode].children.end());
			}
		}

		for (int i = 0; i < q; ++i){
			int query;
			scanf("%d", &query);
			printf("%d\n", answer[query - 1] + 1);
		}
	}
}
