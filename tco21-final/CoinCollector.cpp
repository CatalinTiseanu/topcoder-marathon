/* 
  Solution starts for real at line ~550 
  First hundred lines are macros and useful functions like random sampling
  Vast majority is not used by the actual solution
*/

#pragma GCC target "avx"

#include <algorithm>
#include <bitset>
#include <cassert>
#include <cctype>
#include <climits>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <deque>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <list>
#include <map>
#include <numeric>
#include <queue>
#include <set>
#include <sstream>
#include <stack>
#include <string>
#include <utility>
#include <vector>
#include <array>

#include <sys/time.h>
#include <random>
#include <chrono>
#include <unordered_map>
#include <unordered_set>

#include <x86intrin.h>
#include <emmintrin.h>

using namespace std;
//#define DEBUG              // general information//

template<class T, class V> std::ostream&  operator <<(std::ostream& stream, const pair<T,V> & p) {
  stream << p.first << "," << p.second << " ";
  return stream;
}

template<class T> std::ostream&  operator <<(std::ostream& stream, const vector<T> & v) {
  for (auto el : v)
    stream << el << " ";
  return stream;
}

template<class T> std::ostream&  operator << (std::ostream& stream, const vector< vector<T> > & v) {
  for (auto line : v) {
    for (auto el : line)
      stream << el << " ";
    stream << "\n";
  }
  return stream;
}

class debugger {
public:
  template<class T> void output (const T& v) {
    cerr << v << " ";
  }
  
  template<class T> debugger& operator , (const T& v) {
    output(v);
    return *this;
  }
} dbg;

// Adapted from http://codeforces.com/blog/entry/15643
vector<string> debug_split(const string& s, char c) {
    vector<string> v; stringstream ss(s); string x;
    while (getline(ss, x, c)) v.emplace_back(x); return move(v);}

//void debug_err(vector<string>::iterator it) {cerr << "\n";}
template<typename T, typename... Args>
void debug_err(vector<string>::iterator it, T a, Args... args) {
    cerr << it -> substr((*it)[0] == ' ', it -> length()) << " = " << a << " | ";
    debug_err(++it, args...);
}

// End debug

////////////////
// templates  //
////////////////

// general
//template size
template<typename T> int size(const T& c) { return int(c.size()); }
//template abs
template<typename T> T abs(T x) { return x < 0 ? -x : x; }
//template sqr
template<typename T> T sqr(T x) { return x*x; }
//template remin
template<typename T> bool remin(T& x, T y) { if (x <= y) return false; x = y; return true; }
//template remax
template<typename T> bool remax(T& x, T y) { if (x >= y) return false; x = y; return true; }
//template toStr
template<typename T> string toStr(T x) { stringstream ss; ss << x; return ss.str(); }
//helper len
inline int len(const string & x){ return x.length(); }
//helper print_mask
inline void print_mask(int mask, int n){for (int i = 0; i < n; ++i) cerr << (char)('0' + (mask >> i & 1));}
//helper stoi
//int stoi(const string & s){stringstream ss(s); int res; ss >> res; return res;}
//helper itos
string itos(const int & x){stringstream ss; ss << x; return ss.str();}

// types
//template int64
typedef long long            int64 ;
//template uint64
typedef unsigned long long   uint64 ;
typedef unsigned char uchar;
typedef unsigned short ushort;
typedef int64 hash_type;

// shortcuts
#define all(_xx)             _xx.begin(), _xx.end()

#define pb                   push_back
#define vl                   vector<long long>
#define vs                   vector<string>
#define vvs                  vector<vector<string>>
#define vvi                  vector<vector<int>>
#define vvl                  vector<vector<long long>>
#define vd                   vector<double>
#define vvd                  vector<vector<double>>
#define vc                   vector<char>
#define vvc                  vector<vector<char>>
#define vpdd                 vector<pair<double,double> >

#define vi                   vector<int>
#define vpii                 vector<pair<int,int>> 
//#define vpii                 vector<pair<short,short> >


#define pii                  pair<int, int>
#define pcc                  pair<char, char>
#define pucc                 pair<uchar, uchar>
#define pll                  pair<long long, long long>
#define pdd                  pair<double, double>
#define mp(XX, YY)           make_pair(XX, YY)
#define fi                   first
#define se                   second

#define ll                   long long
#define ull                  unsigned long long
#define SS                   stringstream

#define ptype unsigned short

// for loops
#define re(II, NN)           for (int II(0), _NN(NN); (II) < (_NN); ++(II))
#define fod(II, XX, YY)      for (int II(XX), _YY(YY); (II) >= (_YY); --(II))
#define fo(II, XX, YY)       for (int II(XX), _YY(YY); (II) <= (_YY); ++(II))
#define foreach(it,c) for(__typeof((c).begin()) it=(c).begin();it!=(c).end();++it)

// Useful hardware instructions
#define bitcount             __builtin_popcount
#define gcd                  __gcd

// Useful all around
#define checkbit(n,b)        ((n >> (b)) & 1)
#define pt2(b)               (1LL<<(b))
#define DREP(a)              sort(a.begin(), a.end());a.erase(unique(a.begin(), a.end()),a.end())
#define INDEX(arr,ind)       (lower_bound(arr.begin(), arr.end(),ind)-arr.begin())
#define PAUSE cerr << "Press any key to continue ..."; char xxxx; scanf("%c", &xxxx);

#define v_fill(xx,val) fill(all(xx), val)
#define m_fill(xx,val) memset(xx, val, sizeof(xx))

#define SUM(vv) accumulate(vv.begin(), vv.end(), 0)
#define oo (10001)

#define ALWAYS_INLINE inline __attribute__((always_inline))
#define ALIGNED __attribute__ ((aligned(16)))

#define likely(x)  __builtin_expect(!!(x),1)
#define unlikely(x) __builtin_expect(!!(x),0)

#define SSE_LOAD(a)    _mm_load_si128((const __m128i*)&a)
#define SSE_STORE(a, b) _mm_store_si128((__m128i*)&a, b)


ull getTicks() {
#ifdef __i386
    ull time;
    __asm__ volatile ("rdtsc" : "=A" (time));
    return time;
#else
    ull timelo, timehi;
    __asm__ volatile ("rdtsc" : "=a" (timelo), "=d" (timehi));
    return (timehi << 32) + timelo;
#endif
}

#define SUBMIT

// changed NON SUBMIT from 3.6e9 * 1.32

const double inv_frequency = 1.0 / 2.8e9;

double convertTicks(ull time) {
#ifndef SUBMIT  
    return time / 3.3e9 * 1.23;
#else
    return time * inv_frequency;
    //return time / 3.6e9 * 1.46;
#endif
}

double getTime() {
  return convertTicks(getTicks()) * 1.0;
}

struct Timer {
  ull ticks;
  const string name;
  bool running;
  int count;
  
  Timer() {;}

  Timer(string _name) : name(_name) {
    ticks = 0;
    count = 0;
    running = false;
  }

  void reset() {
    ticks = 0;
    count = 0;
    running = false; 
  }

  double elapsed() {
    return convertTicks(running ? getTicks() - ticks : ticks);
  }

  void start() {
    if (running) return;
    ticks = getTicks() - ticks;
    running = true;
  }
  
  void stop() {
    if (!running) return;
    ticks = getTicks() - ticks;
    running = false;
    count++;
  }
  
  void check() {
    double delta = elapsed();
    cerr << name << ": " << delta << " count: " << count;
    if (count) cerr << " mean per sec: " << (1 / delta) * count;
    cerr << endl;
  }

//  double elapsed() {;}
//  void startTimer() {;}
//  void stopTimer() {;}
//  void checkTimer() {;}

};

// class ChronoTimer {
// public:
//   double start, end, total;
  
//   double getTime()
//   { 
//     timeval tv;
//     gettimeofday(&tv, 0);
//     return (tv.tv_sec + tv.tv_usec * 0.000001);
//   }

//   void startTimer() {
//     start = getTime();
//   }
  
//   double stopTimer() {
//     double elapsed_seconds = checkTimer();
//     // total += elapsed_seconds;
//     // cerr << "ChronoTimer: " << elapsed_seconds << "\n";
//     return elapsed_seconds;
//   }
  
//   double checkTimer() {
//     timeval tv;
//     gettimeofday(&tv, 0);
//     return tv.tv_sec + tv.tv_usec * 0.000001 - start;
//   }
  
//   ChronoTimer() {total = 0.0; startTimer();}
// };

namespace Color {
  enum Code {
    FG_RED      = 31,
    FG_GREEN    = 32,
    FG_YELLOW   = 33,
    FG_BLUE     = 34,
    FG_DEFAULT  = 39,
    BG_RED      = 41,
    BG_GREEN    = 42,
    BG_BLUE     = 44,
    BG_DEFAULT  = 49
  };
  
  class Modifier {
    Code code;
  public:
    Modifier(Code pCode) : code(pCode) {}
    friend std::ostream&
    operator<<(std::ostream& os, const Modifier& mod) {
      return os << "\033[" << mod.code << "m";
    }
  };
}

const uchar
v_white[3]  = { 255, 255, 255 }, v_black[3] = { 0, 0, 0 },     v_red[3] = { 120, 50, 80 },
v_yellow[3] = { 200, 155, 0 },   v_green[3] = { 30, 200, 70 }, v_purple[3] = { 175, 32, 186 },
v_blue[3]   = { 55, 140, 185 },  v_grey[3] = { 127, 127, 127 };

Color::Modifier green(Color::FG_GREEN);
Color::Modifier red(Color::FG_RED);
Color::Modifier yellow(Color::FG_YELLOW);
Color::Modifier blue(Color::FG_BLUE);
Color::Modifier def(Color::FG_DEFAULT);

///////////////////////////////////////////////////////////

template <class T>
pair<T, T> operator - (const pair<T, T> & a, const pair<T, T> & b) {
    return mp(a.fi - b.fi, a.se - b.se);
}

template <class T>
pair<T, T> operator + (const pair<T, T> & a, const pair<T, T> & b) {
    return mp(a.fi + b.fi, a.se + b.se);
}

template <class T, class U>
pair<T, T> operator / (const pair<T, T> & a, const U d) {
    return mp(a.fi / d, a.se / d);
}

template <class T, class U>
pair<T, T> operator * (const U d, const pair<T, T> & a) {
    return mp(a.fi * d, a.se * d);
}

template <typename T>
vector<T> join_vector(const vector<T> & left, const vector<T> & right) {
    if (!size(right)) return left;
    auto res = vector<T>(left);
    res.insert(res.end(), all(right));
    return res;
}

template <typename T>
void concat_vector(vector<T> & left, const vector<T> & right) {
  if (!size(right)) return;
    left.insert(left.end(), all(right));
}

template <typename T>
T sum_vector(vector<T> & v) {
    return accumulate(all(v), 0);
}

template <typename T>
T mean_vector(vector<T> & v) {
    return sum_vector(v) / size(v);
}

vi string_to_vector(const string & s) {
  vi res; for (auto c : s) res.pb(c); return res;
}

template <typename T>
vector<T> operator + (const vector<T>& left, const vector<T>& right) {
  return join_vector(left, right);
}

/////////////////////////////////////////////////////////////

mt19937 random_engine;

double nextDouble() {return 1.0 * (random_engine() % (1<<30) + 1) / ((1<<30)-1);}
int nextInt() {return random_engine();}
int nextInt(int range) {return random_engine() % range;}
int nextInt(int first, int last) {return random_engine() % (last - first + 1) + first;}

/////////////////////////////////////////////////////////////

struct RNG {
  unsigned int x = 123456789;
  unsigned int y = 362436069;
  unsigned int z = 521288629;
  
  unsigned int rand() {
    x ^= x << 16;
    x ^= x >> 5;
    x ^= x << 1;

    unsigned int t = x;
    x = y;
    y = z;
    z = t ^ x ^ y;

    return z;   
  }
    
  inline int nextInt(int x) {
    return rand() % x;
  }
    
  inline int nextInt(int a, int b) {
    return a + (rand() % (b - a + 1));
  }
    
  inline double nextDouble() {
    return (rand() + 0.5) * (1.0 / 4294967296.0);
  }

  void seed(int val) {
    random_engine.seed(val);
    x = nextInt((1<<30));
    y = nextInt((1<<30));
    z = nextInt((1<<30));
  }
};

static RNG rng;

template<typename T>
T ChooseRandom (const vector<T>& vec) {
  assert (size(vec));
  return vec[rng.nextInt(size(vec))];
}

template<typename T>
void FastShuffle(vector<T>& vec) {
    fod (i, size(vec) - 1, 0) {
        swap(vec[i], vec[rng.nextInt(0, i)]);
    }
}
////////////////////////////////////////////////////////////
double dist_l2 (const pii& A, const pii& B) {
  return hypot (1.0*A.fi - B.fi, 1.0 * A.se - B.se);
}

template<typename T>
pair<double, double> computeMeanVariance (const vector<T>& v) {
  pair<double, double> result;
  for (const auto& x : v) {
    result.fi += x;
  } 
  result.fi /= size(v);
  for (const auto& x : v) {
    result.se += sqr(x - result.fi);
  }
  result.se = sqrt(result.se / size(v));
  return result;
}

template<typename T>
double computePercentage (vector<T>& v, int index) {
  return 1.0 * v[index] / sum_vector(v);
}

template<typename T>
vector<double> computeAllPercentages (vector<T>& v) {
  vd result (size(v));
  auto S = sum_vector(v);
  re (i, size(v))
    result[i] = 1.0 * v[i] / S;
  return result;
}

///////////////////////////////////////////////////////////

#ifdef __APPLE__
#define HOME_COMPUTER
#define LOCAL
#endif

//#define LOCAL

#ifdef LOCAL
#define WARN
#define HIGHLIGHT
#define DEBUG
#define INFO
#endif

#ifdef DEBUG
#define debug(args...)            {dbg,args; cerr<<"\n"; cerr.flush();}
#define trace(args...) { vector<string> _v = debug_split(#args, ','); debug_err(_v.begin(), args); }
#define echo(arg) {cerr << #arg << ": "; dbg,arg; cerr << "\n"; cerr.flush();}
#else
#define debug(args...) {}
#define echo(arg) {}
#define trace(args...) {} 
#endif

#ifdef WARN
#define warn(args...)             {dbg,args; cerr<<"\n"; cerr.flush();}
#else
#define warn(args...)
#endif

#ifdef INFO
#define info(args...)             {cerr << blue; dbg,args; cerr<<def<<"\n"; cerr.flush();}  
#else
#define info(args...)
#endif

#ifdef HIGHLIGHT
#define highlight(args...)        {cerr << red; dbg,args; cerr<<def<<"\n"; cerr.flush();}
#else
#define highlight(args...)
#endif

#ifdef LOCAL
#define MAX_TIME 9.0
#define TIME_BUFFER 0.5
#else
// DON'T CHANGE
#define MAX_TIME 9.0
#define TIME_BUFFER 0.5
#endif

#define MAX_SAFE_TIME 4.0
#define EPS 1e-7

#define delta_type double
#define double_type double

//////////////////////////////////////////////

// ACTUAL SOLUTION STARTS HERE

const int INF = (1<<30);
const int MAX_A = 30 * 30;
const int MAX_D = 6;

const int DR[4] = {-1, 0, 1, 0};
const int DC[4] = {0, 1, 0, -1};

const int MAX_LOOKAHEAD = 4;

//////////////////////////////////////////////
//////////////////////////////////////////////
Timer solution_timer("solution");
//////////////////////////////////////////////
unordered_map<string, int64> Counter;

const vvi DICES = {
    {0}, // STAY
    {1},
    {2}, 
    {3},
    {4},
    {5},
    {6},
    {1, 3, 5},
    {2, 4, 6},
    {1, 2, 3},
    {4, 5, 6},
    {1, 2, 3, 4, 5, 6}
};

const char DIRS[] = {'U','R','D','L'};   

unordered_map<string, int> name2id;
unordered_map<int, string> id2name;

int go[MAX_A][4][7];
double DICE_WEIGHT[12] = {0.0};

// 1 = coin, -1 = spike, 0 = nothing
int HIST[MAX_A + 10];

vc is_spike;
vc is_coin;
//////////////////////////////////////////////
int A;
int N;
int D;
int max_score;

int POS_SCORE;

vi REAL_SCORE;

pii O;
vs dice;
vi real_dice;

vvc ogrid;
vc grid;

double timeElapsed;
double timeRemaining;

double ratioS; // percent spikes
double ratioC; // percent coins

vi rng2dice;

struct Move {
    int dir;
    int id;
};

const Move FINISH = {-1, -1};


int p2id(pii a) {
    return a.fi * N + a.se;
}

int p2id(int r, int c) {
    return r * N + c;
}

pii id2p(int id) {
    return pii {id / N, id % N};
}

struct State {
    int score;
    int turn;
    int pos; // player pos
    
    int rem_coins;
    int col_coins;

    int cnt[12]; // might be more efficient to store dices directly

    vi dice; // id's of dices

    vc is_coin;

    State () {
        score = 0;
        turn = 0;
        is_coin = ::is_coin;
        dice = ::real_dice;

        rem_coins = accumulate(all(is_coin), (int)0);
        col_coins = 0;
        m_fill(cnt, 0);
        for (auto d : dice)
            cnt[d]++;
    }

    void RemoveDice(int id) { 
        assert(cnt[id] >= 1);
        cnt[id]--;
    }

    void DoMove(Move & m, int val) {
        int npos = go[pos][m.dir][val];

        RemoveDice(m.id);

        if (is_spike[npos]) {
            score /= 2;
        } else if (is_coin[npos]) {
            is_coin[npos] = 0;
            score += 100;
            col_coins++;
            rem_coins--;
        }

        pos = npos;
        turn++;
    }

    void AddDice(int id) {
        cnt[id]++;
    }
};
Timer rec_timer("rec");

namespace Solver {
    double efficiency;
    double per_turn;

    int start_depth;

    // This linearizes the grid and precomputes which cell one reaches after a given 
    // direction and die value
    // Linearization = (r, c) goes to (r * N + c) which allows viewing verything as an array
    // Directions: UP=0, RIGHT=1, DOWN=2, LEFT=3
    // DICE_WEIGHT[i] = likelihood of dice with id i of occuring as replacement
    // eg. DICE_WEIGHT[0] = 0.04 since id=0 represents STAY and DICE_WEIGHT[11] = 0.24 since id=11 represents
    // RANDOM
    void Init() {
        m_fill(go, -1);

        is_spike.assign(A, 0);
        is_coin.assign(A, 0);
        grid.assign(A, ' ');

        re (r, N) {
            re (c, N) {
                int pid = p2id(r, c);

                is_spike[pid] = ogrid[r][c] == 'S';
                is_coin[pid] = ogrid[r][c] == 'C';

                grid[pid] = ogrid[r][c];

                re (val, 7) {
                    re (d, 4) {
                        int nr = (r + val * DR[d] + N) % N;
                        int nc = (c + val * DC[d] + N) % N;

                        int nid = p2id(nr, nc);

                        go[pid][d][val] = nid;
                    }
                }
            }
        }

        rng2dice.clear();
        re (i, 12) {
            for (auto v : DICES[i]) {
                rng2dice.pb(i);
                DICE_WEIGHT[i] += 1.0 / 25.0;
            }

        }

        assert(rng2dice.size() == 25);
    }

    int depth;
    int max_depth;
    State state;

    double fudgefactor;

    /*
      Recursive lookahead
    */
    pair<double, Move> Rec3() {
        depth++;

        // If reached max search depth return
        if (depth >= max_depth) {
            depth--;
            return pair<double, Move> {state.score, FINISH};
        }

        // Compute what percent of the score gathered so far we need to collect in the future
        // to justify stepping on a spike and incurrign a 50% reduction in score
        // The smaller D, the pessimistic we are about the future and therefore the more we penalize
        // stepping on a spike

        fudgefactor = 0.3 + 0.1 * D;
        if (D == 2)
            fudgefactor = 0.3;
        
        // Estimate how many more turns we have to collect score
        int rem_turns = (A - max_depth);
        // What's the maximum amount of score (coins) we can still gather?
        double eremscore = min(per_turn * rem_turns, (double)state.rem_coins * 100);
        // What's the minimum score we need in order to justify hitting a spike?
        // Obviously if the expected value in the future is higher than right now we continue
        double min_score = state.score - eremscore * fudgefactor;  


        Move bst_move = FINISH;
        double bst_score = min_score;

        // How much more depth do we have ? (eg recursive calls)
        int rem_look = max_depth - depth;

        // Add random noise in order to break ties
        auto jitter = [&] (void) {
            return rng.nextDouble();
        };

        int buf = 0;

        // Simulate all die for D <= 3 (test out for just D = 2)
        // This part is key and got a huge boost for D=2
        // For the very first replacement, instead of assuming a RANDOM die
        // try all 12 replacement possibilities (weighted by they % of occuring, eg. STAY has a 4% chance)
        // and from that continue

        if (D <= 3) {
            buf = 1;
        } else {
            buf = 0;
        }

        // This is critical for pruning
        // Since spikes are so detrimental to the score
        // We compute the most optimistic scenario after hitting a spike which is every 
        // turn after we collect a coin
        // If this is still not enough, then we don't even check this possibility
        // This helps a lot with pruning
        auto compute_max_pos_score = [&] (int dir, int i) {
            int nspikes = 0;
            int NVALS = DICES[i].size();
            for (auto val : DICES[i]) {
                nspikes += is_spike[go[state.pos][dir][val]];
            }

            double max_pos_score = 1.0 * (
                (NVALS - nspikes) * (state.score + min(rem_look, state.rem_coins) * 100) +
                nspikes * (state.score / 2 + min(rem_look-1, state.rem_coins) * 100)
            ) / NVALS;
            return max_pos_score;
        };

        pair<double, Move> ret;

        // For every die we have in hand (state.cnt[i] > 0)
        re (i, 12) if (state.cnt[i]) {
            re (dir, 4) { // For every direction
                double cur_score = 0;
                double max_pos_score = compute_max_pos_score(dir, i); 
                int NVALS = DICES[i].size();

                // If the estimate for the max possible score is an improvement over what we have
                if (max_pos_score > bst_score) {
                    for (auto val : DICES[i]) { // for every value this dice could take
                        // Compute how the next state would look like
                        int npos = go[state.pos][dir][val];

                        state.cnt[i]--;
                        state.turn++;

                        int opos = state.pos;
                        int oscore = state.score;
                        int delta_coin = state.is_coin[npos];

                        if (state.is_coin[npos]) {
                            state.score += 100;
                            state.is_coin[npos] = 0;
                        } else if (is_spike[npos]) {
                            state.score /= 2;
                        }
                        
                        state.pos = npos;

                        // Here is where we do the recursion


                        if (depth + 1 == max_depth) { // Don't replace if we're reaching max depth
                            ret = Rec3();
                            cur_score += (ret.fi + jitter());
                        } else if (depth + 1 <= start_depth + buf) { // if D <= 3 and first replacement
                            fo (die, 0, 11) {
                                state.AddDice(die);
                                ret = Rec3();
                                // Try all possible dice as replacement weighted by their % of occuring
                                cur_score += (ret.fi + jitter()) * DICE_WEIGHT[die];
                                state.RemoveDice(die);
                            }
                        } else { // Just assume we're going to replace the dice with RANDOM (id=11)
                            state.AddDice(11);
                            ret = Rec3();
                            cur_score += (ret.fi + jitter());
                            state.RemoveDice(11);
                        }
                        // Undo the recursive call
                        state.score = oscore;
                        
                        if (delta_coin) {
                            state.is_coin[npos] = delta_coin;
                            state.rem_coins += delta_coin;
                            state.col_coins -= delta_coin;
                        }
                        state.pos = opos;
                    
                        // TUNED: is we're hitting a spike decrease expected score
                        // Not clear how much it helped
                        if (D >= 4 && is_spike[npos]) {
                            cur_score -= 100;
                        }

                        state.cnt[i]++;
                        state.turn--;
                    }

                    cur_score /= NVALS;
                } else {
                    cur_score = -1e3;
                }

                // HACK: I did tune this a lot and it helped
                // Basically thus prioritizes using the RANDOM dices (id=11) first
                // Since the more random a dice, the more unpredictable, there the 'worse' 
                if (i == 11) {
                    // TUNED
                    if (D >= 4) {
                        cur_score += 25;
                    }
                }

                if (cur_score > bst_score) {
                    bst_score = cur_score;
                    bst_move = Move {dir, i};
                }

                // Always do STAY
                // If we have a STAY die (id=0) then we always choose that
                if (!i && depth < A - MAX_LOOKAHEAD) break;
            }

            // Only do dice=STAY for one dir (STAY with dir=0 is the same as for dir=1)
            if (!i) {
                break;
            }
        }
        depth--;
        return pair<double, Move> {bst_score, bst_move};
    }

    Move Solve(const State & real_state) {
        state = real_state;

        efficiency = 1.0 * state.score / (0.1 + state.col_coins * 100);
        per_turn = 1.0 * state.score / (state.turn + 1);

        // We're starting with lookahead (max_depth) = 4
        // As we run out of time we reduce the lookahead to 3
        // Finally, once we hit panic mode we only have a search depth of 2

        int lookahead = 4;

        // Adjust lookahead
        double avg_time = timeElapsed / (state.turn + 1);
        double need_time = avg_time * (A - state.turn);

        // HACK
        if (1.0 * need_time > timeRemaining) {
            lookahead = 3;
        }

        // TIME PANIC (DOUBLE CHECK?)
        if (timeRemaining <= 0.1) {
            lookahead = 2;

        }

        // FOR FINAL SUBMISSION
        // Note: I should have figured out earlier that the safety play is to just finish the game ... :)
        // Would have made my worries around hitting TLE much smaller
        if (timeRemaining < -0.1) {
            lookahead = 1;
        }

        depth = real_state.turn - 1;
        max_depth = min(A, real_state.turn + lookahead);
        start_depth = real_state.turn;
        warn(real_state.turn, "Calling Rec:", lookahead, "->", real_state.score, "eff:", efficiency, "per_turn:", per_turn);
        int delta = max_depth - real_state.turn;
        
        auto ret = Rec3();        
        auto bst_move = ret.second;

        return bst_move;
    }
}

//////////////////////////////////////////////

class CoinCollector {
public:
    int real_turn;
    Move real_last_move;
    State real_state;

    void Log() {
        ofstream fout("logXXX.aux", ofstream::out | ofstream::app);
        fout << N << " " << D << " " << ratioC << " " << ratioS << "\n";
        fout.close();
        exit(0);
    }

    void Info() {
        warn("N:", N, "D:", D, "ratioC:", ratioC, "ratioS:", ratioS);
    }

    CoinCollector() {
        name2id["STAY"] = 0;
        name2id["ONE"] = 1;
        name2id["TWO"] = 2;
        name2id["THREE"] = 3;
        name2id["FOUR"] = 4;
        name2id["FIVE"] = 5;
        name2id["SIX"] = 6;
        name2id["ODD"] = 7;
        name2id["EVEN"] = 8;
        name2id["LOW"] = 9;
        name2id["HIGH"] = 10;
        name2id["RANDOM"] = 11;

        for (auto & e : name2id) {
            id2name[e.se] = e.fi;
        }
    }

    void Init (int _N, int _D, int _r, int _c, vs & _dice, vvc & _grid) {
        N = _N;
        A = sqr(N);
        D = _D;
        O = {_r, _c};
        dice = _dice;
        real_dice.clear();
        for (auto d : dice) {
            real_dice.pb(name2id[d]);
        }
        ogrid = _grid;

        real_turn = 0;

        random_engine.seed(0);
        rng.seed(0); 

        Solver::Init();

        real_state = State();
        real_state.pos = p2id(O);

        ratioC = 0;
        ratioS = 0;

        re (p, A) {
            ratioC += grid[p] == 'C';
            ratioS += grid[p] == 'S';
        }

        POS_SCORE = ratioC * 100;

        ratioC /= A;
        ratioS /= A;
        max_score = 0;

        m_fill(HIST, 0);
        //Log(); 

        Info();
    }

    void Update(int _val, string _name) {
        int bef_score = real_state.score;
        real_state.DoMove(real_last_move, _val);
        int aft_score = real_state.score;

        if (aft_score > bef_score) {
            HIST[real_state.turn - 1] = 1;
        } else if (aft_score < bef_score) {
            HIST[real_state.turn - 1] = -1;
        }

        real_state.AddDice(name2id[_name]);

        REAL_SCORE.pb(real_state.score);
    }

    Move SampleMove() { 
        int diceID = real_turn % D;
        int did = real_state.dice[diceID];
        int dir = real_turn % 4;   

        return Move {dir, did};
    }

    Move GetMove(int _timeElapsed) {
        solution_timer.start();
    
        timeElapsed = _timeElapsed / 1000.0;
        timeRemaining = MAX_TIME - timeElapsed;

        if (timeElapsed >= 9.8) {
            solution_timer.stop();
            return FINISH;
        }

        max_score = max(real_state.score, max_score);
        //Log();

        warn(real_turn, "/", A-1, "[GetMove]", real_turn, "score:", real_state.score, "/", max_score, "elapsed:", timeElapsed);
  
        Move final_move;

        // Always use STAY unless near the end of the game
        if (real_state.rem_coins == 0) {
            final_move = FINISH;
        } else if (real_state.cnt[0] && real_state.turn < A - MAX_LOOKAHEAD) {
            final_move = Move {0, 0};
        } else {
            final_move = Solver::Solve(real_state);
        }
        warn("T:", real_turn, "Move:", final_move.dir, " / ", final_move.id);
        // rec_timer.check();

        real_turn++;
        real_last_move = final_move;

        solution_timer.stop();
        return final_move;
    } 
};

////////////////////

// Below is the boilerplate needed to run the solution
#define PERSONAL

template<class T> void getVector(vector<T>& v) {
    for (int i = 0; i < v.size(); ++i)
        cin >> v[i];
}

void RunTester() {
    std::ios::sync_with_stdio(false);
    cin.tie(0);

    CoinCollector COIN;
    

    int N;
    int D;
    int PersonR;
    int PersonC;
    int timeElapsed;
    int diceValue;

    cin >> N;
    cin >> D;
    cin >> PersonR;
    cin >> PersonC;

    vector<string> dice(D);
    for (int i=0; i< D; i++)
        cin >> dice[i];

    vector< vector<char> > grid(N, vector<char>(N, -1));  
    for (int r=0; r<N; r++)
        for (int c=0; c<N; c++)
            cin >> grid[r][c];

    COIN.Init(N, D, PersonR, PersonC, dice, grid);

    timeElapsed = 0;

    for (int turn=1; turn<=N*N; turn++) {
        auto m = COIN.GetMove(timeElapsed);

        //print the move 

        // Early exit
        if (m.dir == -1) {
            cout << "X" << "\n";
            cout.flush();
            break;
        }

        cout << DIRS[m.dir] << "\n";
        cout << id2name[m.id] << "\n";
        cout.flush();    

        //read the updated information
        string new_dice;
        cin >> timeElapsed;
        cin >> diceValue;
        cin >> new_dice;
    
        COIN.Update(diceValue, new_dice);
    }
}

int main(int argc, char ** argv) {
    RunTester();
    return 0;
}