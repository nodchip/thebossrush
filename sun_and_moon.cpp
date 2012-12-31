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
typedef unsigned long long ull;

#define ALL(c) (c).begin(), (c).end()
#define CLEAR(v) memset(v,0,sizeof(v))
#define MP(a,b) make_pair((a),(b))
#define REP(i,n) for(int i=0;i<(int)n;++i)
#define ABS(a) ((a)>0?(a):-(a))
template<class T> T MIN(const T& a, const T& b) { return a < b ? a : b; }
template<class T> T MAX(const T& a, const T& b) { return a > b ? a : b; }
template<class T> void MIN_UPDATE(T& a, const T& b) { if (a > b) a = b; }
template<class T> void MAX_UPDATE(T& a, const T& b) { if (a < b) a = b; }

ll normalize(ll x, ll m) {
  return ((x % m) + m) % m;
}

bool isPrime(int n) {
	for (int i = 2; i * i <= n; ++i) {
		if (n % i == 0) {
			return false;
		}
	}
	return true;
}

// éüêî, åWêî
map<ll, ll> derive(const map<ll, ll>& polynomial, ll m) {
	map<ll, ll> derived;
	for (map<ll, ll>::const_iterator it = polynomial.begin(), itEnd = polynomial.end(); it != itEnd; ++it) {
		if (it->first == 0) {
			continue;
		}
		derived[it->first - 1] = normalize(normalize(it->first, m) * normalize(it->second, m), m);
	}
	return derived;
}

ll powMod(ll x, ll k, ll m) {
  if (k == 0)     return 1;
  if (k % 2 == 0) return powMod(x*x % m, k/2, m);
  else            return x*powMod(x, k-1, m) % m;
}

// éüêî, åWêî
ll calc(const map<ll, ll>& polynomial, ll x, ll m) {
	ll res = 0;
	for (map<ll, ll>::const_iterator it = polynomial.begin(), itEnd = polynomial.end(); it != itEnd; ++it) {
		res += normalize(powMod(x, it->first, m) * normalize(it->second, m), m);
    res %= m;
	}
	return res;
}

ll trySqrt(ll i) {
  ll sq = sqrt(i);
  for (int di = -2; di <= 2; ++di) {
    if (0 < sq + di && (sq + di) * (sq + di) == i) {
      return sq + di;
    }
  }
  return 0;
}

void addFactors(ll n, set<ll>& factors) {
	if (n < 0) {
		n = -n;
	}

	for (ll i = 1; i * i * i <= n; ++i) {
		if (n % i == 0) {
			factors.insert(i);
			factors.insert(trySqrt(n / i));
		}
	}
}

ll solve(const map<ll, ll>& polynomial) {
  bool zero = true;
  for (map<ll, ll>::const_iterator it = polynomial.begin(), itEnd = polynomial.end(); zero && it != itEnd; ++it) {
    zero = (it->second == 0);
  }
  if (zero) {
    return 1;
  }

	vector<ll> primes;
	for (int i = 1000000000; primes.size() < 10; ++i) {
		if (isPrime(i)) {
			primes.push_back(i);
		}
	}

  set<ll> factors;
	addFactors(polynomial.begin()->second, factors);
  factors.erase(0);

  set<ll>::iterator it = factors.begin();
  for (; it != factors.end(); ++it) {
    bool ok = true;
    for (int i = 0; ok && i < primes.size(); ++i) {
      ll m = primes[i];
      map<ll, ll> derived = derive(polynomial, m);

      ll lh = calc(polynomial, *it, m);
			ll rh = calc(derived, *it, m);
			ok = (lh == 0 && rh == 0);
    }

    if (ok) {
      break;
    }
  }

  if (it != factors.end()) {
    return *it;
	} else {
    return 0;
  }
}

int main() {
	std::ios::sync_with_stdio(false);

	int N;
	cin >> N;

	// éüêî, åWêî
	map<ll, ll> polynomial;
	REP(n, N) {
    ll degree, coefficient;
		cin >> degree >> coefficient;
    polynomial[degree] += coefficient;
	}

  ll answer = solve(polynomial);

  if (answer) {
		cout << "Yes " << answer << endl;
	} else {
    cout << "No" << endl;
  }
}
