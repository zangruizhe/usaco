//
// Created by racer on 2022/8/23.
//
#include <algorithm>
#include <cmath>
#include <deque>
#include <fstream>
#include <iostream>
#include <iterator>
#include <limits>
#include <map>
#include <numeric>
#include <queue>
#include <set>
#include <sstream>
#include <stack>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using namespace std;

#define caf const auto &
#define rep(i, a, b) for (int i = a; i < (b); ++i)
#define all(x) begin(x), end(x)
#define sz(x) (int)(x).size()
#define mp make_pair
#define f first
#define s second
#define getunique(v)                                                           \
  {                                                                            \
    sort(v.begin(), v.end());                                                  \
    v.erase(unique(v.begin(), v.end()), v.end());                              \
  }

using lld = long double;
using ull = unsigned long long;
using ll = long long;

using Index = pair<int, int>;

using pii = pair<int, int>;
using pll = pair<ll, ll>;
using vi = vector<int>;
using vvi = vector<vector<int>>;
using vvvi = vector<vector<vector<int>>>;

using vb = vector<bool>;
using vvb = vector<vector<bool>>;
using vvvb = vector<vector<vector<bool>>>;

double dmax = numeric_limits<double>::max();
int imax = numeric_limits<int>::max();

void usaco(string filename) {
  // #pragma message("be careful, freopen may be wrong")
  freopen((filename + ".in").c_str(), "r", stdin);
  freopen((filename + ".out").c_str(), "w", stdout);
}

auto isPrime = [](int v) {
  if (v < 2)
    return false;
  if (v <= 3)
    return true;
  for (int i = 2; i * i <= v; i++) {
    if (v % i == 0)
      return false;
  }
  return true;
};

vector<string> ReadInput() {
  vector<string> rst;
  std::string line;

  while (std::getline(std::cin, line)) {
    if (line.empty())
      break;
    rst.push_back(line);
  }

  return rst;
}

vector<string> ReadLines(const string &file) {
  std::ifstream input(file);
  vector<string> lines;
  for (std::string line; getline(input, line);)
    lines.push_back(line);
  return lines;
}

template <typename T> vector<T> SplitLine(const string &s) {
  vector<T> rst;

  T t;
  stringstream ss(s);
  while (ss >> t) {
    rst.push_back(t);
  }

  return rst;
}

template <typename T>
void WriteLine(const vector<T> &lines, const string &file) {
  ofstream fout(file);
  for (const auto &l : lines) {
    fout << l << endl;
  }
  fout.close();
}

template <typename T> string Join(const vector<T> &input, const string &split) {
  stringstream ss;
  for (size_t i = 0; i < input.size(); ++i) {
    if (i != 0)
      ss << split;
    ss << input[i];
  }

  return ss.str();
}

template <typename T>
string JoinV2(const vector<T> &lines, string delim = " ") {

  auto strings = vector<string>{};
  std::transform(lines.begin(), lines.end(), back_inserter(strings),
                 [](const T &v) { return to_string(v); });

  return std::accumulate(
      strings.begin(), strings.end(), string(),
      [delim](const std::string &a, const std::string &b) -> std::string {
        return a + (a.empty() ? "" : delim) + b;
      });
}

struct ride {
  vector<string> Solution(const vector<string> &ins) {
    int rst = 1;
    int rst2 = 1;
    for (auto c : ins[0])
      rst *= (c - 'A' + 1);
    for (auto c : ins[1])
      rst2 *= (c - 'A' + 1);
    if ((rst % 47) == (rst2 % 47))
      return vector<string>{"GO"};
    else
      return vector<string>{"STAY"};
  }
}; // namespace ride

struct gift1 {
  vector<string> Solution(const vector<string> &ins) {
    unordered_map<string, int> names;
    int i = 1;
    for (; i <= stoi(ins[0]); ++i) {
      names[ins[i]] = 0;
    }
    while (i < ins.size() - 1) {
      auto name = ins[i];

      int amount;
      int num;
      stringstream ss(ins[++i]);
      ss >> amount;
      ss >> num;

      if (num > 0) {
        names[name] -= amount;
        names[name] += amount % num;
        int give = amount / num;
        while (num > 0) {
          names[ins[++i]] += give;
          num--;
        }
      }
      i++;
    }
    vector<string> rst;
    for (int j = 1; j <= stoi(ins[0]); ++j) {
      rst.push_back(ins[j] + " " + to_string(names[ins[j]]));
    }
    return rst;
  }
}; // namespace gift1

struct friday {
  vector<string> Solution(const vector<string> &lines) {
    int n = stoi(lines[0]);
    int day = 0;
    vector<int> leap{31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
    vector<int> no_leap{31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
    unordered_map<int, int> rst;
    for (int y = 1900; y < 1900 + n; y++) {
      auto month = no_leap;
      if (y % 100 == 0 && y % 400 == 0)
        month = leap;
      if (y % 100 != 0 && y % 4 == 0)
        month = leap;
      for (auto &m : month) {
        rst[(day + 13) % 7]++;
        day += m;
      }
    }

    string ans;
    ans = (to_string(rst[6]));
    ans = ans + " " + (to_string(rst[0]));
    ans = ans + " " + (to_string(rst[1]));
    ans = ans + " " + (to_string(rst[2]));
    ans = ans + " " + (to_string(rst[3]));
    ans = ans + " " + (to_string(rst[4]));
    ans = ans + " " + (to_string(rst[5]));

    return vector<string>{ans};
  }
}; // namespace friday

struct beads {
  int search(const vector<pair<char, int>> &msg, int i) {
    if (i >= msg.size())
      return 0;
    auto start = msg[i];
    auto s = next(msg.begin(), i);

    if (start.first == 'b') {
      auto l = find_if(s, msg.end(),
                       [](const pair<char, int> &v) { return v.first == 'r'; });
      auto r = find_if(l, msg.end(),
                       [](const pair<char, int> &v) { return v.first == 'b'; });

      int rst = 0;
      for (auto j = s; j < r; j++) {
        rst += j->second;
      }
      return rst;
    } else if (start.first == 'r') {
      auto l = find_if(s, msg.end(),
                       [](const pair<char, int> &v) { return v.first == 'b'; });
      auto r = find_if(l, msg.end(),
                       [](const pair<char, int> &v) { return v.first == 'r'; });

      int rst = 0;
      for (auto j = s; j < r; j++) {
        rst += j->second;
      }
      return rst;
    }

    return start.second + search(msg, i + 1);
  }

  vector<string> Solution(const vector<string> &lines) {
    vector<pair<char, int>> msg;
    auto cur = lines[1][0];
    auto cur_n = 0;
    for (auto c : lines[1]) {
      if (c == cur)
        cur_n++;
      else {
        msg.emplace_back(cur, cur_n);
        cur = c;
        cur_n = 1;
      }
    }
    if (cur_n != 0)
      msg.emplace_back(cur, cur_n);

    auto len = msg.size();

    if (len <= 3) {
      int rst = 0;
      for (auto v : msg) {
        rst += v.second;
      }
      return vector<string>{to_string(rst)};
    }

    msg.insert(msg.end(), msg.begin(), msg.end());

    int rst = 0;
    for (int i = 0; i < msg.size(); ++i) {
      rst = max(rst, search(msg, i));
    }

    return vector<string>{to_string(rst)};
  }
}; // namespace beads

struct milk2 {
  vector<string> Solution(const vector<string> &lines) {
    vector<pair<int, int>> nums;
    for (int i = 1; i < lines.size(); ++i) {
      int l;
      int r;
      stringstream ss(lines[i]);
      ss >> l;
      ss >> r;
      nums.emplace_back(l, r);
    }
    sort(nums.begin(), nums.end(), [](auto l, auto r) {
      if (l.first == r.first)
        return l.second > r.second;
      else
        return l.first < r.first;
    });

    vector<pair<int, int>> tmp;
    int l = nums[0].first;
    int r = nums[0].second;
    for (const auto &n : nums) {
      if (n.first == l)
        continue;
      auto li = n.first;
      auto ri = n.second;
      if (li <= r && ri >= l) {
        r = max(r, ri);
      }
      if (li > r) {
        tmp.emplace_back(l, r);
        l = li;
        r = ri;
      }
    }
    tmp.emplace_back(l, r);

    int feed = tmp[0].second - tmp[0].first;
    int no_feed = 0;
    for (int i = 1; i < tmp.size(); i++) {
      feed = max(feed, tmp[i].second - tmp[i].first);
      no_feed = max(no_feed, tmp[i].first - tmp[i - 1].second);
    }

    return vector<string>{to_string(feed) + " " + to_string(no_feed)};
  }
};

struct transform {
  void Flop(vector<vector<char>> &nums) {
    int n = nums.size();
    for (int i = 0; i < n; i++) {
      for (int j = i + 1; j < n; j++) {
        auto tmp = nums[i][j];
        nums[i][j] = nums[j][i];
        nums[j][i] = tmp;
      }
    }
  }

  vector<vector<char>> Rotation90(const vector<vector<char>> &nums) {
    int n = nums.size();
    vector<vector<char>> rst(n, vector<char>());
    for (int i = 0; i < n; i++) {
      rst[i] = nums[n - 1 - i];
    }

    Flop(rst);
    return rst;
  }

  vector<vector<char>> Rotation180(const vector<vector<char>> &nums) {
    return Rotation90(Rotation90(nums));
  }

  vector<vector<char>> Rotation270(const vector<vector<char>> &nums) {
    return Rotation180(Rotation90(nums));
  }

  vector<vector<char>> Reflection(const vector<vector<char>> &nums) {
    int n = nums.size();
    vector<vector<char>> rst(n, vector<char>());
    for (int i = 0; i < n; i++) {
      rst[i].insert(rst[i].begin(), nums[i].crbegin(), nums[i].crend());
    }
    return rst;
  }

  bool isMatch(const vector<vector<char>> &src,
               const vector<vector<char>> &dst) {
    int n = src.size();
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        if (src[i][j] != dst[i][j])
          return false;
      }
    }
    return true;
  }

  vector<string> Solution(const vector<string> &lines) {
    int n = stoi(lines.front());
    vector<vector<char>> src(n, vector<char>());
    vector<vector<char>> dst(n, vector<char>());

    for (int i = 1; i <= n; i++) {
      for (auto &c : lines[i])
        src[i - 1].push_back(c);
    }

    for (int i = 1; i <= n; i++) {
      for (auto &c : lines[i + n])
        dst[i - 1].push_back(c);
    }

    if (isMatch(Rotation90(src), dst)) {
      return vector<string>{"1"};
    }

    if (isMatch(Rotation180(src), dst)) {
      return vector<string>{"2"};
    }

    if (isMatch(Rotation270(src), dst)) {
      return vector<string>{"3"};
    }

    if (isMatch(Reflection(src), dst)) {
      return vector<string>{"4"};
    }

    if (isMatch(src, dst)) {
      return vector<string>{"6"};
    }

    src = Reflection(src);

    if (isMatch(Rotation90(src), dst)) {
      return vector<string>{"5"};
    }

    if (isMatch(Rotation180(src), dst)) {
      return vector<string>{"5"};
    }

    if (isMatch(Rotation270(src), dst)) {
      return vector<string>{"5"};
    }

    return vector<string>{"7"};
  }
};

struct namenum {
  unordered_map<char, vector<char>> choose;
  unordered_set<string> names;
  vector<string> rst;

  void backtrack(string nums, int start, string path) {
    if (path.size() == nums.size()) {
      if (names.find(path) != names.end())
        rst.push_back(path);
      return;
    }
    for (auto &v : choose[nums[start]]) {
      backtrack(nums, start + 1, path + v);
    }
  }

  vector<string> Solution(const vector<string> &lines) {
    std::ifstream infile("dict.txt");
    std::string line;
    while (std::getline(infile, line)) {
      names.insert(line);
    }
    choose['2'] = vector<char>{'A', 'B', 'C'};
    choose['5'] = vector<char>{'J', 'K', 'L'};
    choose['8'] = vector<char>{'T', 'U', 'V'};
    choose['3'] = vector<char>{'D', 'E', 'F'};
    choose['6'] = vector<char>{'M', 'N', 'O'};
    choose['9'] = vector<char>{'W', 'X', 'Y'};
    choose['4'] = vector<char>{'G', 'H', 'I'};
    choose['7'] = vector<char>{'P', 'R', 'S'};
    backtrack(lines[0], 0, "");
    if (rst.empty()) {
      rst.push_back("NONE");
    }
    return rst;
  }
};

struct palsquare {
  bool isPalindrome(const string &num) {
    int l = 0;
    int r = num.size() - 1;
    while (l <= r) {
      if (num[l] != num[r])
        return false;
      l++;
      r--;
    }
    return true;
  }

  string Convert(int y, int base) {
    string num = "0123456789ABCDEFGHIJK";

    string tmp;
    while (y > 0) {
      tmp = num[y % base] + tmp;
      y /= base;
    }
    return tmp;
  }

  vector<string> Solution(const vector<string> &lines) {
    int base = stoi(lines[0]);
    vector<string> rst;
    for (int i = 1; i <= 300; i++) {
      auto y = Convert(i * i, base);

      if (isPalindrome(y)) {
        auto x = Convert(i, base);
        rst.push_back(x + " " + y);
      }
    }
    return rst;
  }
};

struct dualpal {
  bool isPalindrome(const string &num) {
    int l = 0;
    int r = num.size() - 1;
    while (l <= r) {
      if (num[l] != num[r])
        return false;
      l++;
      r--;
    }
    return true;
  }

  string Convert(int y, int base) {
    string num = "0123456789ABCDEFGHIJK";

    string tmp;
    while (y > 0) {
      tmp = num[y % base] + tmp;
      y /= base;
    }
    return tmp;
  }

  vector<string> Solution(const vector<string> &lines) {
    int n;
    int s;
    stringstream ss(lines[0]);
    ss >> n;
    ss >> s;
    s++;

    vector<string> rst;

    while (n > 0) {
      int good = 0;
      for (int i = 2; i <= 10; i++) {
        auto tmp = Convert(s, i);
        if (isPalindrome(tmp))
          good++;
        if (good == 2) {
          rst.push_back(to_string(s));
          n--;
          break;
        }
      }
      s++;
    }

    return rst;
  }
};

struct milk {
  vector<string> Solution(const vector<string> &lines) {
    vector<pair<int, int>> tmp;
    for (auto &v : lines) {
      int n;
      int s;
      stringstream ss(v);
      ss >> n;
      ss >> s;
      tmp.emplace_back(n, s);
    }

    int need = tmp[0].first;
    vector<pair<int, int>> farm(next(tmp.begin()), tmp.end());
    sort(farm.begin(), farm.end(),
         [](pair<int, int> a, pair<int, int> b) { return a.first < b.first; });

    int rst = 0;
    for (auto &v : farm) {
      if (need > v.second) {
        need -= v.second;
        rst += v.first * v.second;
      } else {
        rst += v.first * need;
        need = 0;
        break;
      }
    }
    return vector<string>{to_string(rst)};
  }
};

struct barn1 {
  vector<string> Solution(const vector<string> &lines) {
    int M;
    int S;
    int C;
    stringstream ss(lines[0]);
    ss >> M;
    ss >> S;
    ss >> C;

    vector<int> nums;
    for (int i = 1; i < lines.size(); i++) {
      nums.push_back(stoi(lines[i]));
    }

    sort(nums.begin(), nums.end());

    priority_queue<pair<int, int>> pq;
    for (int i = 1; i < nums.size(); i++) {
      pq.push({nums[i] - nums[i - 1], i - 1});
    }

    unordered_set<int> index;
    for (int i = 0; i < M - 1; i++) {
      if (!pq.empty()) {
        index.insert(pq.top().second);
        pq.pop();
      }
    }

    int rst = 0;
    int pre = 0;
    for (int i = 0; i < nums.size(); i++) {
      if (index.find(i) != index.end()) {
        rst += (nums[i] - nums[pre] + 1);
        pre = i + 1;
      }
    }
    rst += nums.back() - nums[pre] + 1;

    return vector<string>{to_string(rst)};
  }
};

struct crypt1 {
  vector<int> cache3;
  vector<int> cache2;

  void backtrack(const vector<int> &nums, int need_num, int start, int cur,
                 string cur_s) {
    if (need_num == cur_s.length()) {
      if (need_num == 3)
        cache3.push_back(cur);
      if (need_num == 2)
        cache2.push_back(cur);
      return;
    }

    for (int i = start; i < nums.size(); i++) {
      backtrack(nums, need_num, start, cur * 10 + nums[i],
                cur_s + to_string(nums[i]));
    }
  }

  bool IsLegal(const unordered_set<int> &nums, int num, int len) {
    int i = 0;
    while (num > 0) {
      if (nums.find(num % 10) == nums.end())
        return false;
      num /= 10;
      i++;
    }

    return i == len;
  }

  vector<string> Solution(const vector<string> &lines) {
    cache2 = vector<int>();
    cache3 = vector<int>();
    vector<int> nums;
    unordered_set<int> set_nums;
    int n;
    stringstream ss(lines[1]);
    while (ss >> n) {
      set_nums.insert(n);
      nums.push_back(n);
    }
    backtrack(nums, 3, 0, 0, "");
    backtrack(nums, 2, 0, 0, "");
    //        for (auto &v: nums) {
    //            backtrack(nums, 3, 0, v, to_string(v));
    //            backtrack(nums, 2, 0, v, to_string(v));
    //        }

    int rst = 0;

    for (auto &v2 : cache2) {
      // xy
      auto x = v2 % 10;
      auto y = v2 / 10;
      for (auto &v3 : cache3) {
        if (IsLegal(set_nums, x * v3, 3) && IsLegal(set_nums, y * v3, 3) &&
            IsLegal(set_nums, v2 * v3, 4)) {
          //                    cout << v3 << "*" << v2 << "=" << v2 * v3 <<
          //                    endl;
          rst++;
        }
      }
    }

    return vector<string>{to_string(rst)};
  }
};

struct combo {
  unordered_set<string> rst;

  vector<int> Range(int num, int N) {
    int l = num - 2;
    int r = num + 2;

    vector<int> tmp;

    for (int i = l; i <= r; i++) {
      int t = i;
      while (t <= 0)
        t += N;
      while (t > N) {
        int x = t % N;
        if (x == 0)
          t = 1;
        else
          t = x;
      }
      tmp.push_back(t);
    }

    return tmp;
  }

  void backtrack(const vector<int> &key, int N, string cur, int start) {
    if (start >= key.size()) {
      rst.insert(cur);
      return;
    }

    for (auto &v : Range(key[start], N)) {
      backtrack(key, N, cur + ":" + to_string(v), start + 1);
    }
  }

  vector<string> Solution(vector<string> lines) {
    rst = unordered_set<string>();
    vector<int> master;
    vector<int> host;

    int n;
    stringstream ss(lines[1]);
    while (ss >> n) {
      host.push_back(n);
    }

    stringstream ss2(lines[2]);
    while (ss2 >> n) {
      master.push_back(n);
    }

    int N = stoi(lines[0]);
    backtrack(master, N, "", 0);
    backtrack(host, N, "", 0);
    return vector<string>{to_string(rst.size())};
  }
};

struct wormhole {
#define MAX_N 12
  int N;
  vector<int> X;
  vector<int> Y;
  vector<int> partner;
  vector<int> next_on_right;

  bool cycle_exists(void) {
    for (int start = 1; start <= N; start++) {
      // does there exist a cylce starting from start
      int pos = start;
      for (int count = 0; count < N; count++)
        pos = next_on_right[partner[pos]];
      if (pos != 0)
        return true;
    }
    return false;
  }

  // count all solutions
  int solve() {
    // find first unpaired wormhole
    int i, total = 0;
    for (i = 1; i <= N; i++)
      if (partner[i] == 0)
        break;

    // everyone paired?
    if (i > N) {
      if (cycle_exists())
        return 1;
      else
        return 0;
    }

    // try pairing i with all possible other wormholes j
    for (int j = i + 1; j <= N; j++)
      if (partner[j] == 0) {
        // try pairing i & j, let recursion continue to
        // generate the rest of the solution
        partner[i] = j;
        partner[j] = i;
        total += solve();
        partner[i] = partner[j] = 0;
      }
    return total;
  }

  vector<string> Solution(vector<string> lines) {
    N = stoi(lines[0]);
    X = vector<int>(MAX_N + 1);
    Y = vector<int>(MAX_N + 1);
    partner = vector<int>(MAX_N + 1);
    next_on_right = vector<int>(MAX_N + 1);

    for (int i = 1; i < lines.size(); i++) {
      stringstream ss(lines[i]);
      ss >> X[i] >> Y[i];
    }

    for (int i = 1; i <= N; i++) // set next_on_right[i]...
      for (int j = 1; j <= N; j++)
        if (X[j] > X[i] && Y[i] == Y[j]) // j right of i...
          if (next_on_right[i] == 0 || X[j] - X[i] < X[next_on_right[i]] - X[i])
            next_on_right[i] = j;

    return vector<string>{to_string(solve())};
  }
};

struct skidesign {
  vector<string> Solution(vector<string> lines) {
    vector<int> nums;
    for (int i = 1; i < lines.size(); i++)
      nums.push_back(stoi(lines[i]));
    sort(nums.begin(), nums.end());
    int rst = numeric_limits<int>::max();
    for (int i = 0; i < 83; i++) {
      auto l = i;
      auto h = l + 17;
      auto tmp = 0;
      for (auto &v : nums) {
        if (v < l)
          tmp += (l - v) * (l - v);
        else if (v > h)
          tmp += (v - h) * (v - h);
      }
      rst = min(tmp, rst);
    }
    return vector<string>{to_string(rst)};
  }

  vector<string> skidesignV3(vector<string> lines) {
    vector<int> nums;
    for (int i = 1; i < lines.size(); i++)
      nums.push_back(stoi(lines[i]));
    sort(nums.begin(), nums.end());
    auto mid = nums.front() + (nums.back() - nums.front()) / 2;
    auto l = mid - 8;
    auto h = l + 17;
    int rst = 0;
    for (auto &v : nums) {
      if (v < l)
        rst += (l - v) * (l - v);
      else if (v > h)
        rst += (v - h) * (v - h);
    }
    return vector<string>{to_string(rst)};
  }

  vector<string> skidesignV2(vector<string> lines) {
    vector<int> nums;
    for (int i = 1; i < lines.size(); i++)
      nums.push_back(stoi(lines[i]));
    sort(nums.begin(), nums.end());
    int i = 0;
    int j = nums.size() - 1;
    while (i <= j) {
      auto t1 = nums[j] - nums[i];
      auto t2 = nums[j - 1] - nums[i];
      auto t3 = nums[j] - nums[i + 1];
      if (t1 <= 17) {
        break;
      } else if (t2 <= 17 && t3 <= 17) {
        if (t2 < t3) {
          i++;
          break;
        } else {
          j--;
          break;
        }
      } else if (t2 <= 17) {
        j--;
        break;
      } else if (t3 <= 17) {
        i++;
        break;
      } else {
        i++;
        j--;
      }
    }

    int low = nums[i];
    int high = nums[i] + 17;
    cout << low << ":" << high << endl;
    cout << i << ":" << j << endl;
    int rst = 0;
    for (int k = 0; k < i; k++) {
      rst += (low - nums[k]) * (low - nums[k]);
      cout << "rst:" << rst << endl;
    }

    for (int k = j + 1; k < nums.size(); k++) {
      if (nums[k] > high) {
        rst += (nums[k] - high) * (nums[k] - high);
        cout << "rst:" << rst << endl;
      }
    }
    return vector<string>{to_string(rst)};
  }
};

struct ariprog {

  vector<string> Solution(vector<string> lines) {
    int n = stoi(lines[0]);
    int m = stoi(lines[1]);
    vector<int> tmp;
    for (int i = 0; i <= m; i++)
      tmp.push_back(i * i);

    vector<bool> has_v(tmp.back() * 2 + 1, false);
    vector<int> nums;
    for (int i = 0; i < tmp.size(); i++) {
      for (int j = i; j < tmp.size(); j++) {
        int t = tmp[i] + tmp[j];
        nums.push_back(t);
        has_v[t] = true;
      }
    }
    sort(nums.begin(), nums.end());

    priority_queue<pair<int, int>, vector<pair<int, int>>, greater<>> pq;

    for (int j = 0; j < nums.size() - n; j++) {
      if (nums[j] == nums[j + 1])
        continue;

      for (int k = j + 1; k < nums.size(); k++) {
        int i = nums[k] - nums[j];

        if (i < nums.back() / (n - 2)) {
          if (has_v[nums[j] + (n - 1) * i]) {
            int next = nums[j];
            int num = 1;

            while (has_v[next] && num < n) {
              next += i;
              num++;
            }
            if (num == n) {
              pq.push({i, nums[j]});
            }
          }
        } else {
          break;
        }
      }
    }

    vector<string> rst;

    while (!pq.empty()) {
      auto tmp = pq.top();
      pq.pop();
      string s = to_string(tmp.second) + " " + to_string(tmp.first);
      if (!rst.empty()) {
        if (s != rst.back())
          rst.push_back(s);
      } else {
        rst.push_back(s);
      }
    }

    if (rst.empty())
      rst.push_back("NONE");

    return rst;
  }
};

struct milk3 {
  vector<vector<vector<int>>> mem;
  vector<int> cap;
  vector<int> contain;
  unordered_set<int> save;

  void Backtrack(int from, int to) {
    if (mem[contain[0]][contain[1]][contain[2]] == true) {
      return;
    }
    mem[contain[0]][contain[1]][contain[2]] = true;

    if (contain[0] == 0) {
      save.insert(contain[2]);
    }

    for (int i = 0; i < contain.size(); i++) {
      for (int j = 0; j < contain.size(); j++) {
        if (j == i)
          continue;
        if (i == to && j == from)
          continue;
        int tmp = min(cap[j] - contain[j], contain[i]);
        if (tmp != 0) {
          contain[i] -= tmp;
          contain[j] += tmp;

          Backtrack(i, j);

          contain[i] += tmp;
          contain[j] -= tmp;
        }
      }
    }
  }

  vector<string> Solution(vector<string> lines) {
    save = unordered_set<int>();
    int A, B, C;

    stringstream ss(lines[0]);
    ss >> A;
    ss >> B;
    ss >> C;

    mem = vector<vector<vector<int>>>(
        C + 1, vector<vector<int>>(C + 1, vector<int>(C + 1, false)));
    contain = vector<int>(3, 0);
    contain[2] = C;
    cap = vector<int>{A, B, C};

    save.insert(C);

    int tmp = min(cap[0] - contain[0], contain[2]);
    contain[2] -= tmp;
    contain[0] += tmp;
    Backtrack(2, 0);

    vector<int> rst(save.begin(), save.end());

    sort(rst.begin(), rst.end());

    string r = to_string(rst.front());
    for (int i = 1; i < rst.size(); i++) {
      r += " " + to_string(rst[i]);
    }

    return vector<string>{r};
  }
};

struct numtri {
  vector<vector<int>> mem;
  vector<vector<int>> nums;

  int dp(int i, int j) {
    if (i < 0 || j < 0 || j >= nums[i].size()) {
      return 0;
    }
    if (mem[i][j] != -1)
      return mem[i][j];
    int tmp = 0;
    tmp = max(tmp, nums[i][j] + dp(i - 1, j - 1));
    tmp = max(tmp, nums[i][j] + dp(i - 1, j));
    mem[i][j] = tmp;
    return mem[i][j];
  };

  vector<string> Solution(const vector<string> &ins) {
    for (int i = 1; i < ins.size(); i++) {
      nums.push_back(SplitLine<int>(ins[i]));
    }

    for (auto &v : nums) {
      mem.emplace_back(vector<int>(v.size(), -1));
    }

    int i = nums.size() - 1;
    for (int j = 0; j < nums.back().size(); j++) {
      dp(i, j);
    }

    int rst = 0;
    for (int j = 0; j < mem.back().size(); j++) {
      rst = max(mem[i][j], rst);
    }

    return vector<string>{to_string(rst)};
  }
};

struct pprime {
  auto Solution(const vector<string> &ins) {
    vector<int> nums = SplitLine<int>(ins.front());
    int a = nums.front();
    int b = nums.back();

    auto add = [](string &s) {
      // add mid by 1 then carry to up
      bool carry = true;
      int i = (s.size() - 1) / 2;
      while (carry && i >= 0) {
        if (s[i] == '9') {
          s[i] = '0';
          i--;
        } else {
          s[i]++;
          carry = false;
        }
      }

      if (carry) {
        s = "1" + s;
      }
      return;
    };

    auto mirror = [](string &s) {
      for (int i = 0; i < s.size() / 2; i++) {
        s[s.size() - 1 - i] = s[i];
      }
    };

    auto isPrime = [](int v) {
      if (v < 2)
        return false;
      if (v <= 3)
        return true;
      for (int i = 2; i * i <= v; i++) {
        if (v % i == 0)
          return false;
      }
      return true;
    };

    vector<int> rst;
    // gen palindrome number
    auto tmp = to_string(a);
    auto b_s = to_string(b);
    while (tmp.size() < b_s.size() || tmp <= b_s) {
      mirror(tmp);
      auto v = stoi(tmp);
      if (isPrime(v)) {
        rst.push_back(v);
      }
      add(tmp);
    }

    return rst;
  }

  vector<string> SolutionV1(const vector<string> &ins) {
    vector<int> nums = SplitLine<int>(ins.front());
    int a = nums.front();
    int b = nums.back();

    vector<bool> check(b - a + 1, true);

    int j = b / 2 + 1;
    for (int m = 2; m <= j; m++) {
      if (m * m > b)
        break;

      int k = m;
      int q = j;

      while (k <= q) {
        int mid = k + (q - k) / 2;
        if (mid * m < a) {
          k = mid + 1;
        } else {
          q = mid - 1;
        }
      }

      k -= 1;
      for (; k <= j; k++) {
        auto tmp = k * m;
        if (a <= tmp && tmp <= b) {
          check[tmp - a] = false;
        }
        if (tmp > b)
          break;
      }
    }

    auto IsGood = [](int num) {
      vector<int> n;
      while (num > 0) {
        n.push_back(num % 10);
        num /= 10;
      }

      int i = 0;
      int j = n.size() - 1;
      while (i < j) {
        if (n[i] != n[j])
          return false;
        i++;
        j--;
      }
      return true;
    };

    vector<string> rst;
    for (int i = 0; i < check.size(); i++) {
      if (check[i] && IsGood(i + a)) {
        rst.push_back(to_string(i + a));
      }
    }

    return rst;
  }
};

struct sprime {
  vector<int> dp(int n, vector<int> nums) {
    if (n == 0)
      return nums;

    vector<int> tmp;
    for (auto &v : nums) {
      v *= 10;
      for (auto i : vector<int>{1, 2, 3, 5, 7, 9}) {
        if (isPrime(v + i)) {
          tmp.push_back(v + i);
        }
      }
    }
    return dp(n - 1, tmp);
  }

  auto Solution(const vector<string> &ins) {
    // 1 digit primme
    // then add next digit to check
    // continue
    return dp(stoi(ins.front()), vector<int>{0});
  }
};

struct castle {
  const int Wwest = 1;
  const int Wnorth = 2;
  const int Weast = 4;
  const int Wsouth = 8;

  vector<vector<int>> matrix;
  vector<vector<int>> rooms;
  unordered_map<int, int> rooms_size;

  void dfs(int i, int j, int num) {
    if (i < 0 || i >= matrix.size() || j < 0 || j >= matrix.front().size()) {
      return;
    }

    if (rooms[i][j] != -1) {
      return;
    }

    rooms[i][j] = num;
    rooms_size[num]++;

    auto wall = matrix[i][j];
    if (!(wall & Wwest))
      dfs(i, j - 1, num);
    if (!(wall & Wnorth))
      dfs(i - 1, j, num);
    if (!(wall & Weast))
      dfs(i, j + 1, num);
    if (!(wall & Wsouth))
      dfs(i + 1, j, num);
  }

  auto Solution(const vector<string> &ins) {
    for (int i = 1; i < ins.size(); i++) {
      matrix.emplace_back(SplitLine<int>(ins[i]));
    }

    // go through matrix, market every room with number
    // add room size for number
    int m = matrix.size();
    int n = matrix.front().size();
    rooms = vector<vector<int>>(m, vector<int>(n, -1));
    int room_n = 0;
    for (int i = 0; i < m; i++) {
      for (int j = 0; j < n; j++) {
        if (rooms[i][j] == -1) {
          dfs(i, j, room_n);
          room_n++;
        }
      }
    }
    // go through matrix;
    // check if remove every wall add room on both side
    // do not add the same room id
    tuple<int, int, char> rst;
    int cur_max = 0;
    for (int j = 0; j < n; j++) {
      for (int i = m - 1; i >= 0; i--) {
        auto cur_room = rooms[i][j];
        auto cur_room_s = rooms_size[rooms[i][j]];

        if (j > 0 && rooms[i][j - 1] != cur_room) {
          auto next = rooms_size[rooms[i][j - 1]] + cur_room_s;
          if (next > cur_max) {
            cur_max = next;
            rst = tuple<int, int, char>{i, j - 1, 'E'};
          }
        }

        if (i + 1 < m && rooms[i + 1][j] != cur_room) {
          auto next = rooms_size[rooms[i + 1][j]] + cur_room_s;
          if (next > cur_max) {
            cur_max = next;
            rst = tuple<int, int, char>{i + 1, j, 'N'};
          }
        }
      }
    }

    vector<string> rst_v;
    rst_v.push_back(to_string(rooms_size.size()));
    int large = 0;
    for (auto &v : rooms_size) {
      large = max(v.second, large);
    }
    rst_v.push_back(to_string(large));
    rst_v.push_back(to_string(cur_max));
    rst_v.push_back(to_string(get<0>(rst) + 1) + " " +
                    to_string(get<1>(rst) + 1) + " " + string(1, get<2>(rst)));
    return rst_v;
  }
};

struct castleV2 {
  int whole_wall = 15;
  unordered_map<int, vector<pair<int, int>>> walls;

  void BuildDoor(vector<pair<int, pair<int, int>>> &choose, int start, int acc,
                 vector<pair<int, int>> &path) {
    walls[acc] = path;
    if (start >= choose.size()) {
      return;
    }
    for (int i = start; i < choose.size(); i++) {
      const auto &[num, dir] = choose[i];
      path.push_back(dir);
      BuildDoor(choose, i + 1, acc + num, path);
      path.pop_back();
    }
  }

  vector<vector<int>> matrix;
  vector<vector<int>> rooms;
  unordered_map<int, int> rooms_size;
  //    vector<vector<bool>> visited;

  void dfs(int i, int j, int num) {
    if (i < 0 || i >= matrix.size() || j < 0 || j >= matrix.front().size()) {
      return;
    }

    if (rooms[i][j] != -1) {
      return;
    }

    rooms[i][j] = num;
    rooms_size[num]++;

    for (auto &door : walls[15 - matrix[i][j]]) {
      int m = i + door.first;
      int n = j + door.second;
      dfs(m, n, num);
    }
  }

  auto Solution(const vector<string> &ins) {
    for (int i = 1; i < ins.size(); i++) {
      matrix.emplace_back(SplitLine<int>(ins[i]));
    }

    vector<pair<int, pair<int, int>>> choose;
    choose.emplace_back(pair<int, pair<int, int>>{1, pair<int, int>{0, -1}});
    choose.emplace_back(pair<int, pair<int, int>>{2, pair<int, int>{-1, 0}});
    choose.emplace_back(pair<int, pair<int, int>>{4, pair<int, int>{0, 1}});
    choose.emplace_back(pair<int, pair<int, int>>{8, pair<int, int>{1, 0}});

    // build number to the door and wall
    vector<pair<int, int>> path;
    BuildDoor(choose, 0, 0, path);

    // go through matrix, market every room with number
    // add room size for number
    int m = matrix.size();
    int n = matrix.front().size();
    rooms = vector<vector<int>>(m, vector<int>(n, -1));
    int room_n = 0;
    for (int i = 0; i < m; i++) {
      for (int j = 0; j < n; j++) {
        if (rooms[i][j] == -1) {
          dfs(i, j, room_n);
          room_n++;
        }
      }
    }

    // go through matrix;
    // check if remove every wall add room on both side
    // do not add the same room id

    auto dir = [](pair<int, int> wall) -> char {
      if (wall == pair<int, int>{0, -1})
        return 'W';
      if (wall == pair<int, int>{-1, 0})
        return 'N';
      if (wall == pair<int, int>{0, 1})
        return 'E';
      return 'S';
    };

    tuple<int, int, char> rst;
    int cur_max = 0;
    unordered_set<int> visited;
    for (int i = 0; i < m; i++) {
      for (int j = 0; j < n; j++) {
        auto cur_room = rooms[i][j];

        visited.clear();
        visited.insert(cur_room);
        for (auto &wall : walls[matrix[i][j]]) {
          int x = i + wall.first;
          int y = j + wall.second;

          if (x < 0 || x >= matrix.size() || y < 0 ||
              y >= matrix.front().size()) {
            continue;
          }

          auto next_room = rooms[x][y];

          if (visited.find(next_room) == visited.end()) {
            visited.insert(cur_room);

            auto tmp = rooms_size[next_room] + rooms_size[cur_room];
            auto cur = tuple<int, int, char>{i, j, dir(wall)};
            if (tmp > cur_max) {
              cur_max = tmp;
              rst = cur;
            } else if (tmp == cur_max) {
              auto &[old_i, old_j, old_dir] = rst;
              if (j < old_j) {
                rst = cur;
              } else if (j == old_j) {
                if (i > old_i) {
                  rst = cur;
                } else if (i == old_i) {
                  if (dir(wall) == 'N') {
                    rst = cur;
                  }
                }
              }
            }
          }
        }
      }
    }

    vector<string> rst_v;
    rst_v.push_back(to_string(rooms_size.size()));
    int large = 0;
    for (auto &v : rooms_size) {
      large = max(v.second, large);
    }
    rst_v.push_back(to_string(large));
    rst_v.push_back(to_string(cur_max));
    rst_v.push_back(to_string(get<0>(rst) + 1) + " " +
                    to_string(get<1>(rst) + 1) + " " + string(1, get<2>(rst)));
    return rst_v;
  }
};

struct frac1 {
  unordered_map<double, pair<int, int>> mem;
  vector<vector<bool>> visited;

  void BackTrack(int N, int up, int down) {
    if (up > N || down > N)
      return;
    if (up > down)
      return;

    if (visited[up][down])
      return;
    visited[up][down] = true;

    auto tmp = (double)(up) / down;
    if (mem.find(tmp) == mem.end() || mem[tmp].first > up) {
      mem[tmp] = pair<int, int>{up, down};
    }

    BackTrack(N, up, down + 1);
    BackTrack(N, up + 1, down);
  }

  auto Solution(const vector<string> &lines) {
    auto N = stoi(lines.front());
    // try all combination
    visited = decltype(visited)(N + 1, vector<bool>(N + 1, false));
    BackTrack(N, 1, 1);

    vector<string> rst;
    rst.push_back("0/1");

    // push pq to rst;
    vector<double> key;
    for (auto &v : mem) {
      key.push_back(v.first);
    }
    sort(key.begin(), key.end());
    for (auto &v : key) {
      auto ans = mem[v];
      rst.push_back(to_string(ans.first) + "/" + to_string(ans.second));
    }
    //        rst.push_back("1/1");
    return rst;
  }
};

struct sort3 {
  using pq_type = priority_queue<int, vector<int>, greater<>>;

  pq_type q1{};
  pq_type q2{};
  pq_type q3{};

  void Swap(vector<int> &nums, int i, int j) {
    auto tmp = nums[i];
    nums[i] = nums[j];
    nums[j] = tmp;
  }

  int Sort(vector<int> &nums, int last) {
    auto SwapByQueue = [&](pq_type &q1, pq_type &q2, int cur) -> int {
      auto next = q2.top();
      cout << next << "<-" << cur << endl;
      for (int i = 0; i < nums.size(); i++) {
        if (i == next || i == cur) {
          cout << "[" << nums[i] << "]"
               << " ";
        } else {
          cout << nums[i] << " ";
        }
      }
      cout << endl;

      Swap(nums, cur, next);

      q1.push(next);

      q2.pop();
      q2.push(cur);

      return 1 + Sort(nums, cur) + Sort(nums, next);
    };

    // put cur to the right pos
    auto v = nums[last];
    if (v == 1) {
      if (!q3.empty() && q3.top() < last)
        return SwapByQueue(q1, q3, last);
      else if (!q2.empty() && q2.top() < last)
        return SwapByQueue(q1, q2, last);
    } else if (v == 2) {
      if (!q3.empty() && q3.top() < last)
        return SwapByQueue(q2, q3, last);
    }
    return 0;
  }

  auto Solution(const vector<string> &line) {
    auto nums = vector<int>{};
    for (int i = 1; i < line.size(); i++) {
      nums.push_back(stoi(line[i]));
    }

    // put index in queue
    for (int i = 0; i < nums.size(); i++) {
      auto v = nums[i];
      if (v == 1)
        q1.push(i);
      else if (v == 2)
        q2.push(i);
      else
        q3.push(i);
    }

    for (auto &v : nums) {
      cout << v << " ";
    }
    cout << endl;

    auto rst = int{0};
    for (int i = nums.size() - 1; i >= 0; i--) {
      rst += Sort(nums, i);
    }
    return vector<int>{rst};
  }
};

struct sort3V2 {
  void Swap(vector<int> &nums, int i, int j) {

    cout << i << "<-" << j << endl;
    for (int k = 0; k < nums.size(); k++) {
      if (k == i || k == j) {
        cout << "[" << nums[k] << "]"
             << " ";
      } else {
        cout << nums[k] << " ";
      }
    }
    cout << endl;

    auto tmp = nums[i];
    nums[i] = nums[j];
    nums[j] = tmp;
  }

  auto Solution(const vector<string> &line) {
    vector<int> nums;
    for (int i = 1; i < line.size(); i++) {
      nums.push_back(stoi(line[i]));
    }

    // collect 2 position
    vector<int> index3;
    vector<int> index2;
    vector<int> index1;
    for (int i = 0; i < nums.size(); i++) {
      if (nums[i] == 2)
        index2.push_back(i);

      if (nums[i] == 1)
        index1.push_back(i);

      if (nums[i] == 3)
        index3.push_back(i);
    }

    int rst = 0;

    // loop nums,
    // if find 1 > first 2, swap
    // if find 3 < end 2, swap
    // after swap make sure 2 index is sorted
    for (int i = 0; i < nums.size(); i++) {
      auto v = nums[i];
      if (index2.empty()) {
        if (index1.empty())
          continue;
        else {
          if (v == 3 && i < index1.back()) {
            rst++;

            Swap(nums, i, index1.back());
            index1.erase(prev(index1.end()));
            index1.push_back(i);
            sort(index1.begin(), index1.end());
          }
        }
      } else {
        if (v == 1 && i > index2.front()) {
          rst++;

          Swap(nums, i, index2.front());
          index2.erase(index2.begin());
          index2.push_back(i);
          sort(index2.begin(), index2.end());

        } else if (v == 3 && i < index2.back()) {

          rst++;
          Swap(nums, i, index2.back());
          index2.erase(prev(index2.end()));
          index2.push_back(i);
          sort(index2.begin(), index2.end());
        }
      }
    }
    if (nums.back() == 2 && index3.empty() != true) {
      for (int i = 0; i < nums.size(); i++) {
        if (i == 3) {
          rst++;
          Swap(nums, i, nums.size() - 1);
          break;
        }
      }
    }

    return vector<int>{rst};
  }
};

struct sort3V3 {

  vector<int> nums;
  vector<int> final_nums;
  vector<bool> visited;
  vector<vector<int>> graph;

  int rst = 0;
  bool dfs(int start_i, int next_i, unordered_set<int> &path, int &rst) {
    // find the cycle and add the item num in the cycle
    if (start_i == next_i) {
      // path is good, add rst; and mark visited;
      for (auto &i : path) {
        visited[i] = true;
      }
      rst += (path.size() - 1);
      cout << start_i << " ; rst=" << rst << endl;
      return true;
    }
    // walk through graph by next_v
    for (auto &i : graph[final_nums[next_i]]) {
      if (visited[i] == false && path.find(i) == path.end() && i >= start_i) {
        path.insert(i);
        if (dfs(start_i, i, path, rst)) {
          return true;
        }
        path.erase(i);
      }
    }
    return false;
  }

  auto Solution(const vector<string> &line) {
    for (int i = 1; i < line.size(); i++) {
      nums.push_back(stoi(line[i]));
    }
    final_nums = nums;
    sort(final_nums.begin(), final_nums.end());
    visited = vector<bool>(nums.size(), false);

    graph = vector<vector<int>>(4, vector<int>{});
    for (int i = 0; i < nums.size(); i++) {
      graph[nums[i]].push_back(i);
    }

    unordered_set<int> path;
    for (int i = 0; i < nums.size(); i++) {
      if (final_nums[i] == nums[i]) {
        visited[i] = true;
      } else if (!visited[i]) {
        path.clear();
        for (auto j : graph[final_nums[i]]) {
          if (visited[j])
            continue;

          path.insert(j);
          if (dfs(i, j, path, rst)) {
            break;
          }
          path.erase(j);
        }
      }
    }

    return vector<int>{rst};
  }
};

struct sort3V4 {
  vector<int> nums;
  vector<int> final_nums;
  vector<bool> visited;
  vector<vector<int>> graph;
  int rst = 0;

  auto Solution(const vector<string> &line) {
    for (int i = 1; i < line.size(); i++) {
      nums.push_back(stoi(line[i]));
    }
    final_nums = nums;
    sort(final_nums.begin(), final_nums.end());

    graph = vector<vector<int>>(4, vector<int>{});
    for (int i = 0; i < nums.size(); i++) {
      graph[nums[i]].push_back(i);
    }

    visited = vector<bool>(nums.size(), false);

    for (int i = 0; i < nums.size(); i++) {
      if (visited[i])
        continue;
      if (nums[i] == final_nums[i])
        continue;
      else {
        // find swap dir index
        bool find = false;
        for (auto &j : graph[final_nums[i]]) {
          if (visited[j] == false && final_nums[j] == nums[i]) {
            rst++;
            visited[i] = true;
            visited[j] = true;
            find = true;
            break;
          }
        }
        if (!find) {
          for (auto &j : graph[final_nums[i]]) {
            if (visited[j] == false) {
              for (auto &k : graph[final_nums[j]]) {
                if (visited[k] == false && final_nums[k] == nums[i]) {
                  rst += 2;
                  visited[i] = true;
                  visited[j] = true;
                  visited[k] = true;
                  break;
                }
              }
              break;
            }
          }
        }
      }
    }
    cout << "rst=" << rst << endl;

    return vector<int>{rst};
  }
};

struct holstein {
  vector<int> vita;

  vector<vector<int>> foods;

  int rst_waste = numeric_limits<int>::max();
  vector<int> rst;

  void dfs(int food_i, int vita_i, vector<int> &need_vita, vector<int> &path) {

    if (vita_i >= need_vita.size() || all_of(need_vita.begin(), need_vita.end(),
                                             [](int v) { return v <= 0; })) {

      if (rst.empty() || path.size() < rst.size()) {
        rst = path;
        return;
      } else if (path.size() == rst.size()) {
        if (path.back() - path.front() < rst.back() - rst.front()) {
          rst = path;
          return;
        }
      }
      return;
    }

    if (!rst.empty() && path.size() > rst.size()) {
      return;
    }

    if (need_vita[vita_i] <= 0) {
      dfs(food_i, vita_i + 1, need_vita, path);
    }

    for (int i = food_i; i < foods.size(); i++) {
      const auto &food = foods[i];

      // choose with food need to use
      if (foods[i][vita_i] <= 0)
        continue;

      for (int j = 0; j < need_vita.size(); j++) {
        need_vita[j] -= food[j];
      }
      path.push_back(i);

      dfs(i + 1, vita_i, need_vita, path);

      path.pop_back();

      for (int j = 0; j < need_vita.size(); j++) {
        need_vita[j] += food[j];
      }
    }
  }

  auto Solution(const vector<string> &lines) {
    vita = SplitLine<int>(lines[1]);
    for (int i = 3; i < lines.size(); i++) {
      foods.push_back(SplitLine<int>(lines[i]));
    }

    vector<int> path{};
    dfs(0, 0, vita, path);
    for (auto &v : rst) {
      v++;
    }
    return vector<string>{to_string(rst.size()) + " " + Join<int>(rst, " ")};
  }
};

struct hamming {
  vector<int> rst;
  ll rst_sum = -1;

  vector<vector<bool>> isLegal;
  // Hamming distance of num i j is > D

  bool dfs(int N, int right_limit, int start, vector<int> &path, ll path_sum) {
    // pick up from num=0, try to build N codeword
    // get the smallest combine

    if (path.size() == N) {
      if (rst.empty() || rst_sum > path_sum) {
        rst = path;
        rst_sum = path_sum;
      }
      return true;
    }

    for (int i = start; i <= right_limit; i++) {
      bool is_good = all_of(path.begin(), path.end(),
                            [this, i](int v) { return this->isLegal[i][v]; });

      if (is_good) {
        path.push_back(i);

        if (dfs(N, right_limit, i + 1, path, path_sum + i)) {
          return true;
        }

        path.pop_back();
      }
    }
    return false;
  }

  int HammingDistance(int i, int j) {
    auto tmp = i ^ j;
    auto rst = int{0};
    while (tmp > 0) {
      rst += tmp % 2;
      tmp /= 2;
    }
    return rst;
  }

  auto Solution(const vector<string> &lines) {
    auto input = SplitLine<int>(lines[0]);
    auto N = input[0];
    auto B = input[1];
    auto D = input[2];
    auto right_limit = (1 << (B + 1)) - 1;

    // build isLegal
    isLegal = decltype(isLegal)(right_limit + 1,
                                vector<bool>(right_limit + 1, false));
    for (int i = 0; i <= right_limit; i++) {
      for (int j = 0; j <= right_limit; j++) {
        if (HammingDistance(i, j) >= D) {
          isLegal[i][j] = true;
        }
      }
    }

    cout << isLegal[0][7] << endl;

    // pick up from num=0, try to build N codeword
    // get the smallest combine
    auto path = vector<int>{};
    dfs(N, right_limit, 0, path, 0);

    vector<string> rst_str;
    int i = 0;
    for (; i + 10 < rst.size(); i += 10) {
      rst_str.push_back(Join(
          vector<int>(next(rst.begin(), i), next(rst.begin(), i + 10)), " "));
    }

    rst_str.push_back(Join(vector<int>(next(rst.begin(), i), rst.end()), " "));
    return rst_str;
  }
};

struct preface {
  vector<string> Roman{"I", "V", "X", "L", "C", "D", "M"};
  vector<int> Roman_i{1, 5, 10, 50, 100, 500, 1000};

  vector<pair<int, string>> Roman_choice;
  // vector<int, string> choice;

  vector<string> nums;

  string Backtrack(int target) {
    if (target < 1)
      return "";

    if (target == 1)
      return "I";

    if (!nums[target].empty())
      return nums[target];

    auto it = prev(upper_bound(
        Roman_choice.begin(), Roman_choice.end(), target,
        [](int t, const pair<int, string> &v) { return t < v.first; }));

    auto it_str = it->second;

    auto tmp = string{};
    if (it->first == target) {
      tmp = it_str;
    } else {
      tmp = it_str + Backtrack(target - it->first);
    }
    nums[target] = tmp;
    return tmp;
  }

  auto Solution(const vector<string> &lines) {
    auto N = stoi(lines.front());
    nums = vector<string>(N + 1, "");

    // build choice
    for (int i = 0; i < Roman.size(); i++) {
      auto s = Roman[i];
      auto v = Roman_i[i];

      Roman_choice.push_back({v, s});
      if (v == 1 || v == 10 || v == 100 || v == 1000) {
        Roman_choice.push_back({v + v, s + s});
        Roman_choice.push_back({v + v + v, s + s + s});
        if (v != 1000) {
          Roman_choice.push_back({Roman_i[i + 1] - v, s + Roman[i + 1]});
          Roman_choice.push_back({Roman_i[i + 2] - v, s + Roman[i + 2]});
        }
      }
    }

    sort(Roman_choice.begin(), Roman_choice.end(),
         [](const pair<int, string> &l, pair<int, string> &r) {
           return l.first < r.first;
         });

    unordered_map<string, int> tmp;
    for (auto i = N; i >= 1; i--) {
      string rst_str = Backtrack(i);
      for (auto &c : rst_str) {
        tmp[string(1, c)]++;
      }
    }

    auto rst = vector<string>{};
    for (auto &v : Roman) {
      if (tmp[v] != 0) {
        rst.push_back(v + " " + to_string(tmp[v]));
      }
    }

    return rst;
  }
};
struct subset {
  vector<vector<int>> mem;

  int dp(int limit, int target) {
    if (target == 0)
      return 1;
    if (mem[limit][target] != -1)
      return mem[limit][target];

    auto tmp = int{0};
    for (int i = limit; i >= 1; i--) {
      if (target - i >= 0) {
        tmp += dp(i - 1, target - i);
      }
    }
    mem[limit][target] = tmp;
    return tmp;
    // num to combine target under limit
  }

  auto Solution(const vector<string> &lines) {
    auto N = stoi(lines.front());
    if (N == 1 || (N * (N + 1) / 2) % 2 != 0)
      return vector<int>{0};

    auto need = N * (N + 1) / 4;
    mem = vector<vector<int>>(N + 1, vector<int>(need + 1, -1));

    return vector<int>{dp(N - 1, need - N)};
  }
};

struct runround {
  bool GoodNum(const string &M) {
    if (M.find('0') != string::npos)
      return false;
    auto tmp = M;
    getunique(tmp);
    if (tmp.size() != M.size())
      return false;

    auto visited = vector<bool>(M.size(), false);
    auto next = (M[0] - '0') % M.size();
    while (visited[next] == false) {
      visited[next] = true;
      next = (next + (M[next] - '0')) % M.size();
    }
    return all_of(visited.begin(), visited.end(), [](bool v) { return v; });
  }

  auto Solution(const vector<string> &lines) {
    auto M = stoi(lines.front());
    M++;

    while (M < numeric_limits<int>::max()) {
      if (GoodNum(to_string(M))) {
        return vector<int>{M};
      }
      M++;
    }

    auto rst = vector<int>{};
    return rst;
  }
};

struct lamps {
  unordered_set<string> rst{};
  int N;
  int C;
  vector<int> ON;
  vector<int> OFF;

  unordered_map<string, unordered_set<int>> mem;

  bool isGood(const string &state) {
    auto is_on = all_of(ON.begin(), ON.end(),
                        [&state](int i) { return state[i] == '1'; });
    if (!is_on)
      return false;

    auto is_off = all_of(OFF.begin(), OFF.end(),
                         [&state](int i) { return state[i] == '0'; });
    if (!is_off)
      return false;

    return true;
  }

  bool HaveState(int N, const string &state) {
    const auto &cur = mem[state];
    return cur.find(N) != cur.end();
  }

  string op1(string state) {
    for (int i = 1; i < state.size(); i++) {
      state[i] = state[i] == '1' ? '0' : '1';
    }
    return state;
  };
  string op2(string state) {
    for (int i = 0; i * 2 + 1 < state.size(); i++) {
      auto v = i * 2 + 1;
      state[v] = state[v] == '1' ? '0' : '1';
    }
    return state;
  };
  string op3(string state) {
    for (int i = 1; i * 2 < state.size(); i++) {
      auto v = i * 2;
      state[v] = state[v] == '1' ? '0' : '1';
    }
    return state;
  };

  string op4(string state) {
    for (int i = 0; i * 3 + 1 < state.size(); i++) {
      auto v = i * 3 + 1;
      state[v] = state[v] == '1' ? '0' : '1';
    }
    return state;
  };

  void dfs(int N, string &&state) {
    if (N == 0) {
      if (isGood(state)) {
        auto tmp = state;
        tmp.erase(tmp.begin());
        rst.insert(tmp);
      }
      return;
    }

    if (HaveState(N, state))
      return;

    mem[state].insert(N);

    dfs(N - 1, op1(state));
    dfs(N - 1, op2(state));
    dfs(N - 1, op3(state));
    dfs(N - 1, op4(state));
  }

  auto Solution(const vector<string> &lines) {
    N = stoi(lines[0]);
    C = stoi(lines[1]);
    ON = SplitLine<int>(lines[2]);
    ON.pop_back();
    OFF = SplitLine<int>(lines[3]);
    OFF.pop_back();

    dfs(C, string(N + 1, '1'));

    if (rst.empty())
      return vector<string>{"IMPOSSIBLE"};

    auto ans = vector<string>(rst.begin(), rst.end());
    sort(ans.begin(), ans.end());
    return ans;
  }
};

struct prefix {
  vector<string> choice;
  vector<int> mem;
  int rst = 0;
  int dfs(const string &target, int start) {
    if (mem[start] != -1)
      return mem[start];

    if (start >= target.size())
      return 0;

    int tmp = 0;
    for (const auto &v : choice) {
      bool can_use = true;
      for (int i = 0; i < v.size(); i++) {
        if (v[i] != target[start + i]) {
          can_use = false;
          break;
        }
      }
      if (can_use) {
        tmp = max(tmp, int(dfs(target, start + v.size()) + v.size()));
      }
    }

    mem[start] = tmp;
    return tmp;
  }

  auto Solution(const vector<string> &lines) {
    int i = 0;
    while (lines[i] != ".") {
      for (const auto &v : SplitLine<string>(lines[i])) {
        choice.push_back(v);
      }
      i++;
    }

    auto target = string{};
    i++;
    while (i < lines.size()) {
      target += lines[i];
      i++;
    }

    // cout << lines << endl;
    // cout << "choice=" << choice << endl;
    mem = decltype(mem)(target.size() + 1, -1);
    auto rst = dfs(target, 0);

    return vector<int>{rst};
  }
};

struct nocows {
  const int mod = 9901;

  vector<vector<ll>> less_high;
  vector<vector<ll>> with_high;

  vector<ll> with_num;

  ll BuildTreeLessThanHigh(int N, int high) {

    if (less_high[N][high] != -1)
      return less_high[N][high];

    if (N == 0 && high == 0) {
      less_high[N][high] = 1;
      return 1;
    }

    if ((N == 0 || high == 0) || (N > with_num[high - 1])) {
      less_high[N][high] = 0;
      return 0;
    }

    ll tmp = 0;
    tmp = BuildTreeWithHigh(N, high - 1) + BuildTreeLessThanHigh(N, high - 1);
    tmp %= mod;

    less_high[N][high] = tmp;
    // cout << "BuildTreeLessThanHigh: N=" << N << " ; high=" << high << endl;
    return tmp;
  }

  ll BuildTreeWithHigh(int N, int high) {
    // ways to build tree with high
    if ((N < (2 * high - 1)) || (N > with_num[high])) {
      with_high[N][high] = 0;
      return 0;
    }

    if (with_high[N][high] != -1)
      return with_high[N][high];

    if (high == 1 && N == 1) {
      with_high[N][high] = 1;
      return 1;
    }

    if ((high == 1 || N == 1)) {
      with_high[N][high] = 0;
      return 0;
    }

    N--;
    ll tmp = 0;
    for (int i = 0; i * 2 + 1 <= N / 2; i++) {
      auto v = i * 2 + 1;
      auto c = (2 * v == N) ? 1 : 2;

      tmp += c * (BuildTreeWithHigh(v, high - 1) *
                      BuildTreeLessThanHigh(N - v, high - 1) +
                  // l == high , r < high
                  BuildTreeLessThanHigh(v, high - 1) *
                      BuildTreeWithHigh(N - v, high - 1) +
                  // l < high, r == high
                  BuildTreeWithHigh(v, high - 1) *
                      BuildTreeWithHigh(N - v, high - 1)
                  // l == high r == high
                 );
      tmp %= mod;
    }
    tmp %= mod;
    with_high[N][high] = tmp;
    // cout << "BuildTreeWithHigh: N=" << N << " ; high=" << high << endl;
    return tmp;
  }

  auto Solution(const vector<string> &lines) {
    auto tmp = SplitLine<int>(lines.front());
    auto N = tmp.front();
    auto K = tmp.back();

    less_high = decltype(less_high)(N + 1, vector<ll>(K + 1, -1));
    with_high = decltype(less_high)(N + 1, vector<ll>(K + 1, -1));
    for (int i = 0; i <= K; i++) {
      with_num.push_back((1 << i) - 1);
    }

    auto rst = BuildTreeWithHigh(N, K);

    return vector<ll>{rst};
  }

  auto SolutionV2(const vector<string> &lines) {
    int table[101][202];
    int smalltrees[101][202];
    auto tmp = SplitLine<int>(lines.front());
    auto N = tmp.front();
    auto K = tmp.back();

    table[1][1] = 1;
    for (int i = 2; i <= K; i++) {
      for (int j = 1; j <= N; j += 2)
        for (int k = 1; k <= j - 1 - k; k += 2) {
          auto c = (k != j - 1 - k) ? 2 : 1;

          table[i][j] += c * (smalltrees[i - 2][k] * table[i - 1][j - 1 - k] +
                              table[i - 1][k] * smalltrees[i - 2][j - 1 - k] +
                              table[i - 1][k] * table[i - 1][j - 1 - k]);
          table[i][j] %= mod;
        }
      for (int k = 0; k <= N; k++) {
        smalltrees[i - 1][k] += table[i - 1][k] + smalltrees[i - 2][k];
        smalltrees[i - 1][k] %= mod;
      }
    }
    int rst = table[K][N];
    cout << rst << endl;
    return vector<int>();
  }
};

struct zerosum {
  vector<int> nums;
  vector<string> rst;

  int Cal(const string &path) {
    int rst = 0;
    int i = 0;
    int cur = 0;
    auto op = '+';
    for (const auto &c : path) {
      // cout << "rst=" << rst << "; cur=" << cur << endl;
      if (c == '+' || c == '-') {
        rst += (op == '+' ? cur : -cur);
        // reset rst
        op = c;
        cur = 0;
      } else if (c != ' ') {
        cur = cur * 10 + (c - '0');
      }
      i++;
    }
    rst += (op == '+' ? cur : -cur);
    return rst;
  }

  void dfs(int cur_i, const string &path) {
    if (cur_i == nums.size() - 1) {
      auto tmp = path;
      tmp += to_string(nums.back());
      if (Cal(tmp) == 0) {
        rst.push_back(tmp);
      }
      return;
    }
    dfs(cur_i + 1, path + to_string(nums[cur_i]) + "+");
    dfs(cur_i + 1, path + to_string(nums[cur_i]) + "-");
    dfs(cur_i + 1, path + to_string(nums[cur_i]) + " ");
  }

  auto Solution(const vector<string> &lines) {
    auto N = stoi(lines.front());

    for (int i = 1; i <= N; i++) {
      nums.push_back(i);
    }

    auto nums = stack<int>{};
    dfs(0, "");
    sort(rst.begin(), rst.end());
    return rst;
  }
};

struct money {
  vector<int> coins;
  vector<vector<ll>> mem;

  ll dp(int N, int cur) {
    if (N < 0 || cur >= coins.size())
      return 0;

    if (N == 0)
      return 1;

    if (mem[N][cur] != -1)
      return mem[N][cur];

    ll tmp = 0;
    for (int i = cur; i < coins.size(); i++) {
      auto v = N - coins[i];
      if (v >= 0)
        tmp += dp(v, i);
    }

    mem[N][cur] = tmp;
    return tmp;
  }

  ll DP_Array(int N) {
    mem = vector<vector<ll>>(coins.size() + 1, vector<ll>(N + 1, 0));
    // i = coin index
    // j = money needed
    for (int i = 0; i < mem.size(); i++) {
      mem[i][0] = 1;
    }

    for (int i = 1; i <= coins.size(); i++) {
      for (int j = 1; j <= N; j++) {
        // do not use this type of coin
        mem[i][j] = mem[i - 1][j];

        // use one of this type of coin
        auto c = coins[i - 1];
        if (j >= c) {
          mem[i][j] += mem[i][j - c];
        }
      }
    }
    return mem[coins.size()][N];
  }

  auto Solution(const vector<string> &lines) {
    auto tmp = SplitLine<int>(lines.front());
    auto V = tmp.front();
    auto N = tmp.back();
    for (int i = 1; i < lines.size(); i++) {
      for (auto &v : SplitLine<int>(lines[i])) {
        coins.push_back(v);
      };
    }

    // mem table
    // mem = vector<vector<ll>>(N + 1, vector<ll>(coins.size() + 1, -1));
    // auto rst = dp(N, 0);

    // dp table
    auto rst = DP_Array(N);

    return vector<ll>{rst};
  }
};

struct concom {
  vector<vector<pair<int, int>>> graph;

  void dfs(vector<int> &sub_company, unordered_map<int, int> &cnt_company,
           vector<int> &rst) {
    if (sub_company.empty())
      return;
    // get connect c from group[sub] add c to rst
    for (auto &s : sub_company) {
      for (auto &c : graph[s]) {
        // add to cnt
        // if cnt duplicate then add percent
        cnt_company[c.first] += c.second;
      }
      rst.push_back(s);
    }
    sub_company.clear();

    // remove percent > 50, add to sub_cpy
    for (auto &v : cnt_company) {
      if (v.second > 50 && find(rst.begin(), rst.end(), v.first) == rst.end()) {
        sub_company.push_back(v.first);
      }
    }
    for (auto &v : sub_company) {
      cnt_company.erase(v);
    }
    dfs(sub_company, cnt_company, rst);
  }

  auto Solution(const vector<string> &lines) {
    auto tmp_company = unordered_set<int>{};
    graph = decltype(graph)(101, vector<pair<int, int>>{});

    // build graph
    for (int i = 1; i < lines.size(); i++) {
      auto tmp = SplitLine<int>(lines[i]);
      if (tmp[0] != tmp[1]) {
        graph[tmp[0]].push_back({tmp[1], tmp[2]});
      }

      tmp_company.insert(tmp[0]);
      tmp_company.insert(tmp[1]);
    }

    auto companies = vector<int>(tmp_company.begin(), tmp_company.end());
    sort(companies.begin(), companies.end());

    auto rst = vector<pair<int, int>>{};
    for (int i = 0; i < companies.size(); i++) {
      auto sub_cpy = vector<int>{};
      auto cnt_cpy = unordered_map<int, int>{};
      auto start = companies[i];
      for (auto &v : graph[start]) {
        if (v.second > 50) {
          sub_cpy.push_back(v.first);
        } else {
          cnt_cpy[v.first] += v.second;
        }
      }

      auto sub_rst = vector<int>{start};

      dfs(sub_cpy, cnt_cpy, sub_rst);

      for (const auto &v : sub_rst) {
        if (start != v) {
          rst.push_back({start, v});
        }
      }
    }

    sort(rst.begin(), rst.end());
    auto tmp = vector<string>{};
    for (auto &v : rst) {
      tmp.push_back(to_string(v.first) + " " + to_string(v.second));
    }
    return tmp;
  }
};

struct ttwo {
  int rst = 0;

  vector<vector<char>> board;
  vector<pair<int, int>> dirs{
      {-1, 0},
      {0, 1},
      {1, 0},
      {0, -1},
  };

  pair<Index, int> Move(const Index &src, const int &dir) {
    auto i = src.first + dirs[dir].first;
    auto j = src.second + dirs[dir].second;
    if (i < 0 || i >= board.size() || j < 0 || j >= board.front().size() ||
        board[i][j] == '*') {
      return {src, (dir + 1) % 4};
    }
    return {make_pair(i, j), dir};
  }

  void dfs(Index F, int F_dir, Index C, int C_dir, int min) {
    // cout << "F = " << F << " ; " << F_dir << " ; C=" << C << " ; " << C_dir
    //      << " ; min=" << min << endl;
    if (min > 100 * 3) {
      rst = 0;
      return;
    }

    if (F == C) {
      rst = min;
      return;
    }

    auto [N_F, N_F_dir] = Move(F, F_dir);
    auto [N_C, N_C_dir] = Move(C, C_dir);
    dfs(N_F, N_F_dir, N_C, N_C_dir, min + 1);
  }

  auto Solution(const vector<string> &lines) {
    for (int i = 0; i < 10; i++) {
      board.push_back(SplitLine<char>(lines[i]));
    }
    // get F/c
    auto F = Index{0, 0};
    auto C = Index{0, 0};
    for (int i = 0; i < board.size(); i++) {
      for (int j = 0; j < board.front().size(); j++) {
        if (board[i][j] == 'F') {
          F = make_pair(i, j);
        } else if (board[i][j] == 'C') {
          C = make_pair(i, j);
        }
      }
    }

    dfs(F, 0, C, 0, 0);

    return vector<int>{rst};
  }
};

struct maze1 {
  vector<string> maze;
  vector<Index> GetDoor() {
    vector<Index> rst;
    for (int j = 1; j + 2 < maze.front().size(); j += 2) {
      if (maze[0][j] == ' ') {
        rst.push_back({1, j});
      }
      if (maze[maze.size() - 1][j] == ' ') {
        rst.push_back({maze.size() - 2, j});
      }
    }

    for (int i = 1; i < maze.size(); i += 2) {
      if (maze[i][0] == ' ') {
        rst.push_back({i, 1});
      }
      if (maze[i][maze.front().size() - 1] == ' ') {
        rst.push_back({i, maze.front().size() - 2});
      }
    }

    return rst;
  }

  vector<Index> dir{{-1, 0}, {1, 0}, {0, -1}, {0, 1}};
  pair<bool, Index> TryMove(Index start, Index v) {
    auto i = v.first * 2 + start.first;
    auto j = v.second * 2 + start.second;
    if (i < 0 || i >= maze.size() || j < 0 || j >= maze.front().size()) {
      return {false, {}};
    }
    if (maze[v.first + start.first][v.second + start.second] == ' ') {
      return {true, make_pair(i, j)};
    }

    return {false, {}};
  }

  vector<vector<int>> visited;

  void dfs(Index start, int dist) {
    visited[start.first][start.second] =
        min(dist, visited[start.first][start.second]);
    auto rst = 0;
    for (auto &v : dir) {
      auto [ok, x] = TryMove(start, v);
      if (ok && visited[x.first][x.second] > dist + 1) {
        dfs(x, dist + 1);
      }
    }
    return;
  }

  auto Solution(const vector<string> &lines) {
    auto tmp = SplitLine<int>(lines.front());
    maze = vector<string>(tmp.back() * 2 + 1, string(tmp.front() * 2 + 1, ' '));
    for (int i = 1; i < lines.size(); i++) {
      for (int j = 0; j < lines[i].size(); j++) {
        maze[i - 1][j] = lines[i][j];
      }
    }

    visited = vector<vector<int>>(
        maze.size(),
        vector<int>(maze.front().size(), numeric_limits<int>::max()));

    auto doors = GetDoor();
    dfs(doors.front(), 1);
    dfs(doors.back(), 1);

    auto rst = 0;
    for (auto &v : visited) {
      for (auto &k : v) {
        if (k != numeric_limits<int>::max()) {
          rst = max(rst, k);
        }
      }
    }
    return vector<int>{rst};
  }
};

struct cowtour {
  vector<Index> points;
  vector<vector<double>> graph_points_distance;
  vector<vector<int>> graph;
  vector<int> group;

  void dfs(int cur, int group_num) {
    if (group[cur] != -1)
      return;
    group[cur] = group_num;
    for (auto &v : graph[cur]) {
      dfs(v, group_num);
    }
  }

  double distance(int a, int b) {
    const auto &l = points[a];
    const auto &r = points[b];

    return sqrt(pow(l.first - r.first, 2) + pow(l.second - r.second, 2));
  }

  auto Solution(const vector<string> &lines) {
    auto N = stoi(lines.front());
    graph_points_distance = vector<vector<double>>(
        N, vector<double>(N, numeric_limits<double>::max()));

    // get points
    for (int i = 1; i <= N; i++) {
      auto tmp = SplitLine<int>(lines[i]);
      points.push_back({tmp.front(), tmp.back()});
    }

    // get graph
    graph = decltype(graph)(N, vector<int>{});
    for (int i = N + 1; i < lines.size(); i++) {
      for (int j = 0; j < N; j++) {
        auto k = i - N - 1;
        if (lines[i][j] == '1' && k != j) {
          graph_points_distance[k][j] = distance(k, j);
          graph[k].push_back(j);
        }
      }
    }

    // mark group
    group = vector<int>(N, -1);
    int group_num = 0;
    for (int i = 0; i < points.size(); i++) {
      caf p = points[i];
      if (group[i] == -1) {
        dfs(i, group_num);
        group_num++;
      }
    }

    // cal distance
    for (int k = 0; k < N; k++) {
      for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
          if (i == j) {
            graph_points_distance[i][j] = 0;
          } else if (graph_points_distance[i][j] >
                     graph_points_distance[i][k] +
                         graph_points_distance[k][j]) {
            graph_points_distance[i][j] =
                graph_points_distance[i][k] + graph_points_distance[k][j];
          }
        }
      }
    }

    auto diam = vector<double>(N, 0);
    auto group_diam = vector<double>(group_num, 0);
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
        if (graph_points_distance[i][j] != dmax &&
            diam[i] < graph_points_distance[i][j]) {
          diam[i] = graph_points_distance[i][j];
        }
      }

      if (diam[i] > group_diam[group[i]]) {
        group_diam[group[i]] = diam[i];
      }
    }

    auto rst = numeric_limits<double>::max();

    for (int i = 0; i < points.size(); i++) {

      auto ans = dmax;
      for (int j = 0; j < points.size(); j++) {
        if (group[i] == group[j])
          continue;

        auto tmp = diam[i] + distance(i, j) + diam[j];
        tmp = max(tmp, group_diam[group[i]]);
        tmp = max(tmp, group_diam[group[j]]);
        ans = min(ans, tmp);
      }

      rst = min(rst, ans);
    }

    return vector<string>{to_string(rst)};
  }
};

struct comehome {
  vector<vector<int>> group;

  void bfs(int start, vector<int> &dist) {
    auto pq =
        priority_queue<pair<int, int>, vector<pair<int, int>>, greater<>>{};
    pq.push({0, start});

    while (!pq.empty()) {
      auto cur = pq.top();
      pq.pop();

      auto c_i = cur.second;
      auto c_dst = cur.first;
      auto c_c = (char)('A' + c_i);

      if (dist[c_i] != -1 && dist[c_i] <= c_dst)
        continue;

      dist[c_i] = cur.first;
      for (int i = 0; i < group.front().size(); i++) {
        if (group[c_i][i] == -1)
          continue;

        auto tmp = c_dst + group[c_i][i];
        if (dist[i] != -1 && dist[i] < tmp)
          continue;
        pq.push({tmp, i});
      }
    }

    return;
  }

  void dfs(int start, int dst, vector<int> &dist) {
    if (dist[start] != -1 && dist[start] <= dst) {
      return;
    }

    dist[start] = dst;

    // if ('A' + start == 'Z') {
    //   return;
    // }

    for (int i = 0; i < group.size(); i++) {
      if (group[start][i] == -1)
        continue;
      dfs(i, dst + group[start][i], dist);
    }
  }

  pair<int, char> bruteforce() {
    auto rst = imax;
    auto rst_c = char{};
    auto dst = group;

    for (int k = 0; k < group.size(); k++) {
      for (int i = 0; i < group.size(); i++) {
        for (int j = 0; j < group.size(); j++) {
          if (i == j) {
            group[i][j] = 0;
            dst[i][j] = 0;
            continue;
          }
          if (dst[k][i] == -1 || dst[k][j] == -1)
            continue;

          auto tmp = dst[k][i] + dst[k][j];
          if (dst[i][j] == -1 || dst[i][j] > tmp) {
            dst[i][j] = tmp;
          }
        }
      }
    }

    for (auto i = 0; i < group.size(); i++) {
      auto c = char('A' + i);
      if ('A' <= c && c < 'Z') {
        auto t = dst[i]['Z' - 'A'];
        if (t >= 0 && rst > t) {
          rst = t;
          rst_c = c;
        }
      }
    }

    return {rst, rst_c};
  }

  auto Solution(const vector<string> &lines) {
    group = decltype(group)('z' - 'A' + 1, vector<int>('z' - 'A' + 1, -1));
    // build graph
    for (int i = 1; i < lines.size(); i++) {
      auto s = stringstream(lines[i]);
      char f, t;
      int dist;
      s >> f >> t >> dist;
      if (f == t) {
        group[f - 'A'][t - 'A'] = 0;
        group[t - 'A'][f - 'A'] = 0;

      } else if (group[f - 'A'][t - 'A'] == -1 ||
                 group[f - 'A'][t - 'A'] > dist) {
        group[f - 'A'][t - 'A'] = dist;
        group[t - 'A'][f - 'A'] = dist;
      }
    }

    // auto [rst, rst_c] = bruteforce();

    // get whole cow
    auto rst = imax;
    auto rst_c = char{};
    auto tmp = vector<int>(group.size(), -1);
    dfs('Z' - 'A', 0, tmp);
    // bfs('Z' - 'A', tmp);
    for (int i = 0; i < 'Z' - 'A'; i++) {

      auto t = tmp[i];
      if (t >= 0 && rst > t) {
        rst = t;
        rst_c = 'A' + i;
      }
    }

    return vector<string>{string(1, rst_c) + " " + to_string(rst)};
  }
};

struct fracdec {
  void dfs(int N, int D, string &post, vector<int> &remainder) {
    auto iter = find(remainder.begin(), remainder.end(), N);
    if (iter != remainder.end()) {
      auto i = distance(remainder.begin(), iter);
      post.insert(i, "(");
      post += ")";
      return;
    }

    remainder.push_back(N);

    auto tmp = (N * 10) / D;
    post += to_string(tmp);

    auto r = (N * 10) % D;
    if (r == 0) {
      return;
    }

    dfs(r, D, post, remainder);
  }

  auto Solution(const vector<string> &lines) {
    auto tmp = SplitLine<int>(lines.front());
    auto N = tmp.front();
    auto D = tmp.back();

    auto pre = to_string(N / D);
    N %= D;

    auto post = string{""};
    auto remainder = vector<int>{};
    dfs(N, D, post, remainder);
    // cout << remainder << endl;

    auto rst = vector<string>{};
    int i = 0;
    auto rst_str = pre + "." + post;
    for (auto &c : rst_str) {
      if (rst.empty() || rst.back().size() >= 76) {
        rst.push_back(string(1, c));
      } else {
        rst.back() += c;
      }
    }

    return rst;
  }
};

struct agrinet {
  vvi matrix;
  vb InTree;
  priority_queue<pii, vector<pii>, greater<>> pq;

  auto Solution(const vector<string> &lines) {
    auto N = stoi(lines.front());
    auto tmp = vector<int>();
    for (int i = 1; i < lines.size(); i++) {
      if (tmp.size() == N) {
        matrix.push_back(tmp);
        tmp.clear();
      }

      auto t = SplitLine<int>(lines[i]);
      for (auto &v : t) {
        tmp.push_back(v);
      }
    }

    matrix.push_back(tmp);

    pq = priority_queue<pii, vector<pii>, greater<>>{};
    InTree = vector<bool>(N, false);
    auto cost = 0;

    auto Cut = [this](int start) {
      for (int i = 0; i < matrix.size(); i++) {
        if (matrix[start][i] != 0 && InTree[i] == false) {
          pq.push({matrix[start][i], i});
        }
      }
    };

    pq.push({0, 0});

    while (!pq.empty()) {
      auto cur = pq.top();
      pq.pop();
      if (InTree[cur.second])
        continue;
      else {
        cost += cur.first;
        InTree[cur.second] = true;
        Cut(cur.second);
      }
    }

    auto rst = vector<int>{cost};
    return rst;
  }
};

struct inflate {
  vector<pii> question;
  vvi mem;
  vi dp_array;
  int dp(int time, int end) {
    if (time <= 0)
      return 0;
    if (mem[time][end] != -1)
      return mem[time][end];

    auto tmp = 0;
    for (int i = end; i > 0; i--) {
      const auto &v = question[i - 1];
      if (time >= v.second) {
        tmp = max(tmp, v.first + dp(time - v.second, i));
      }
    }

    mem[time][end] = tmp;
    return tmp;
  }

  int dp_loop2(int M) {
    for (int i = 0; i < question.size(); i++) {
      for (int j = 1; j <= M; j++) {
        caf q = question[i];

        if (j >= q.second) {
          dp_array[j] = max(dp_array[j], dp_array[j - q.second] + q.first);
        }
      }
    }
    return dp_array.back();
  }

  int dp_loop(int M) {
    for (int i = 1; i <= question.size(); i++) {
      for (int j = 1; j <= M; j++) {
        caf q = question[i - 1];

        mem[i][j] = max(mem[i][j], mem[i - 1][j]);
        if (j >= q.second) {
          mem[i][j] = max(mem[i][j], mem[i][j - q.second] + q.first);
        }
      }
    }

    return mem.back().back();
  }

  int brute_force2(int M) {
    auto rst = 0;
    while (M > 0) {
      auto tmp = 0;
      auto tmp_M = pii{};
      for (caf q : question) {
        if (M < q.second)
          continue;

        auto v = (M / q.second) * q.first;
        if (v > tmp) {
          tmp = v;
          tmp_M = q;
        }
      }

      if (tmp == 0)
        break;

      M -= tmp_M.second;
      rst += tmp_M.first;
    }
    return rst;
  }

  int brute_force(int M) {
    auto cmp = [](pii a, pii b) {
      auto A = double(a.first) / a.second;
      auto B = double(b.first) / b.second;
      if (abs(A - B) <= 0.00001) {
        return a.second > b.second;
      } else {
        return A < B;
      }
    };
    auto pq = priority_queue<pii, vector<pii>, decltype(cmp)>(cmp);
    for (auto &q : question) {
      pq.push(q);
    }

    auto ans = 0;
    while (M > 0 && !pq.empty()) {
      auto q = pq.top();
      pq.pop();
      ans += (M / q.second) * q.first;
      M %= q.second;
    }
    return ans;
  }

  auto Solution(const vector<string> &lines) {
    auto tmp = SplitLine<int>(lines.front());
    auto M = tmp.front();
    auto N = tmp.back();
    question = vector<pii>(N);

    for (int i = 1; i < lines.size(); i++) {
      const auto &tmp = SplitLine<int>(lines[i]);
      question[i - 1] = make_pair(tmp[0], tmp[1]);
    }
    sort(question.begin(), question.end(),
         [](pii a, pii b) { return a.second < b.second; });

    // solve 1
    // mem = vvi(M + 1, vi(N + 1, -1));
    // auto ans = dp(M, N);

    // solve 2
    // mem = vvi(M + 1, vi(N + 1, -1));
    // auto ans = brute_force(M);
    // solve 3
    // auto ans = brute_force2(M);

    // solve 4
    // mem = vvi(N + 1, vi(M + 1, 0));
    // auto ans = dp_loop(M);

    // solve 5 best
    dp_array = vi(M + 1, 0);
    auto ans = dp_loop2(M);

    auto rst = vector<int>{ans};
    return rst;
  }
};

struct humble {
  vector<int> choice;

  auto Solution(const vector<string> &lines) {
    auto tmp = SplitLine<int>(lines.front());
    auto k = tmp.front();
    auto N = tmp.back();
    choice = SplitLine<int>(lines.back());
    sort(choice.begin(), choice.end());
    auto i_mark = vi(choice.size(), 0);
    auto rst = vector<ll>{1};

    while (rst.size() <= N) {
      auto up = rst.back();
      auto next = numeric_limits<ll>::max();
      // check every c
      // cal c with rst[j], j is from last j c use
      // find the first j in rst, that let c *j > up
      // get min up, that is the next
      // next

      rep(i, 0, choice.size()) {
        caf c = choice[i];

        rep(j, i_mark[i], rst.size()) {
          if (c * rst[j] > up) {
            i_mark[i] = j;
            next = min(next, c * rst[j]);
            break;
          }
        }
      }

      rst.push_back(next);
    }

    return vector<ll>{rst[N]};
  }
};

struct contact {
  auto Solution(const vector<string> &lines) {
    auto tmp = SplitLine<int>(lines.front());
    auto A = tmp[0];
    auto B = tmp[1];
    auto N = tmp[2];
    auto nums = string();
    rep(i, 1, lines.size()) { nums += lines[i]; }

    auto mem = unordered_map<string_view, int>{};

    rep(i, A, (B + 1)) {
      if (nums.size() < i)
        break;
      rep(j, 0, i) {
        for (int k = j; k < nums.size(); k += i) {
          if (k + i - 1 >= nums.size())
            continue;
          mem[string_view(&nums[k], i)]++;
        }
      }
    }

    // loop through map get value set
    using Item = pair<int, string_view>;
    auto cmp = [](Item a, Item b) -> bool {
      if (a.first == b.first) {
        if (a.second.size() == b.second.size()) {
          return a.second > b.second;
        } else {
          return a.second.size() > b.second.size();
        }
      } else {
        return a.first < b.first;
      }
    };

    auto items = priority_queue<Item, vector<Item>, decltype(cmp)>(cmp);
    for (auto &v : mem) {
      items.push({v.second, v.first});
    }

    // loop val set
    auto cache = map<int, vector<string>>{};
    while (!items.empty()) {
      auto v = items.top();
      items.pop();
      cache[v.first].push_back(string(v.second));
    }

    auto rst = vector<string>{};

    int n = 0;
    for (auto i = cache.crbegin(); i != cache.crend() && n < N; i++, n++) {
      rst.push_back(to_string(i->first));
      caf v = i->second;
      for (int j = 0; j < v.size(); j++) {
        if (j % 6 == 0)
          rst.push_back(v[j]);
        else
          rst.back() += " " + v[j];
      }
    }

    return rst;
  }
};

struct stamps {
  vi nums;

  int dp_table(int k) {
    auto mem = vi(k * nums.back() + 1, numeric_limits<int>::max());
    mem[0] = 0;
    for (int i = 0; i < nums.size(); i++) {
      for (int j = 0; j < mem.size(); j++) {
        if (mem[j] < k) {
          mem[j + nums[i]] = min(mem[j + nums[i]], 1 + mem[j]);
        }
      }
    }

    int i = 0;
    while (i < mem.size() && mem[i] <= k) {
      i++;
    }

    return i - nums.front();
  }

  auto Solution(const vector<string> &lines) {
    auto tmp = SplitLine<int>(lines.front());
    auto k = tmp.front();

    for (int i = 1; i < lines.size(); i++) {
      auto tmp = SplitLine<int>(lines[i]);
      for (auto &v : tmp) {
        nums.push_back(v);
      }
    }
    sort(all(nums));
    auto rst = dp_table(k);

    return vector<int>{rst};
  }
};

struct fact4 {
  auto Solution(const vector<string> &lines) {
    auto N = stoi(lines.front());

    auto GetNum5 = [](int N) {
      int rst = 0;
      while (N >= 5) {
        rst += (N / 5);
        N /= 5;
      }
      return rst;
    };

    auto num_2 = GetNum5(N);

    auto rst = ll{1};
    for (int i = 1; i <= N; i++) {
      auto tmp = i;
      while (tmp % 5 == 0)
        tmp /= 5;

      while (tmp % 2 == 0 && num_2 > 0) {
        tmp /= 2;
        num_2--;
      }
      rst = (rst * tmp) % 10;
    }
    return vector<ll>{rst};
  }
};
struct kimbits {
  vector<int> nums;

  bool dp(ll &need, ll num_1, int first_1, string &cur) {
    if (need == 0)
      return true;

    if (num_1 == 0)
      return false;

    for (int i = cur.size() - 1; i > first_1; i--) {
      if (cur[i] != '1' && num_1 > 0) {
        cur[i] = '1';

        if (dp(--need, num_1 - 1, first_1, cur)) {
          return true;
        }

        cur[i] = '0';
      } else {
        break;
      }
    }
    return false;
  }

  string dp_table(ll N, ll L, ll I) {
    auto CreatNum = [&]() {
      auto i = int{int(N) - 1};
      auto t = ll{1};
      while (i >= 0) {
        nums[i] = t;
        t *= 2;
        i--;
      }
    };

    CreatNum();

    auto str = string(N, '0');
    for (int i = str.size() - 1; i >= 0; i--) {
      str[i] = '1';

      if (dp(--I, L - 1, i, str)) {
        break;
      }
      str[i] = '0';
    }
    return str;
  }

  string dp_table2(int N, int L, int I) {
    auto cache = priority_queue<int, vector<int>, less<>>{};

    auto GetNum = [&](const string &s) {
      auto tmp = 0;
      for (int i = 0; i < s.size(); i++) {
        if (s[i] == '1') {
          tmp += nums[i];
        }
      }
      return tmp;
    };

    auto toString = [N](int i) -> string {
      auto tmp = string{};
      while (i > 0) {
        int m = i % 2;
        tmp = to_string(m) + tmp;
        i /= 2;
      }
      while (tmp.size() < N) {
        tmp = "0" + tmp;
      }
      return tmp;
    };

    cache.push(0);
    auto q = queue<int>();
    q.push(0);

    auto Loop = [&]() {
      auto len = q.size();
      for (int j = 0; j < len; j++) {
        auto s = q.front();
        q.pop();
        // cout << toString(s) << endl;

        for (int k = N - 1; k >= 0; k--) {
          caf v = nums[k];
          if ((s & v) == v)
            break;

          auto tmp = s | v;
          q.push(tmp);
          cache.push(tmp);
          if (cache.size() > I) {
            cache.pop();
          }
        }
      }
    };

    for (int i = 1; i <= L; i++) {
      Loop();
    }

    return toString(cache.top());
  }

  string toString(int N, int i) {
    auto tmp = string{""};
    while (i > 0) {
      int m = i % 2;
      tmp = to_string(m) + tmp;
      i /= 2;
    }
    while (tmp.size() < N) {
      tmp = "0" + tmp;
    }
    return tmp;
  };

  string TrySolve4(int N, int L, int I) {
    auto rst = string(N, '0');
    auto dp = vvi(N, vi(L + 1, 0));

    for (int i = N - 1; i >= 0; i--) {
      for (int j = 1; j <= L; j++) {
        if (i == N - 1) {
          dp[i][j] = 1;
        } else {
          dp[i][j] = dp[i + 1][j - 1] + dp[i + 1][j];
        }
      }
    }

    I--;
    for (int l = L; l > 0; l--) {
      for (int i = N - 1; i >= 0; i--) {
        auto find = false;
        if (I == 0) {
          return rst;
        }

        if (I > dp[i][l]) {
          I -= dp[i][l];
        } else {
          I--;
          rst[i] = '1';
          break;
        }
      }
    }

    return rst;
  }

  auto Solution(const vector<string> &lines) {
    auto tmp = SplitLine<ll>(lines.front());
    auto N = tmp[0];
    auto L = tmp[1];
    auto I = tmp[2];
    nums = vector<int>(N, 0);

    if (I == 1) {
      return vector<string>{string(N, '0')};
    }
    if (N == L) {
      return vector<string>{toString(N, I - 1)};
    }

    // auto str = dp_table(N, L, I - 1);
    // auto str = dp_table2(N, L, I);
    auto str = TrySolve4(N, L, I);

    auto rst = vector<string>{str};
    return rst;
  }
};

struct spin {
  vector<pair<int, vector<pii>>> wheel{};

  vector<pii> Join(const vector<pii> &l, const vector<pii> &r) {
    auto tmp = vector<pii>{};
    auto parse = [](const vector<pii> &l) {
      auto tmp_l = vector<pii>{};
      for (auto l_v : l) {
        if (l_v.first > l_v.second) {
          tmp_l.push_back(make_pair(l_v.first, 0));
          tmp_l.push_back(make_pair(0, l_v.second));
        } else {
          tmp_l.push_back(l_v);
        }
      }
      return tmp_l;
    };

    auto tmp_l = parse(l);
    auto tmp_r = parse(r);

    for (auto l_t : tmp_l) {
      for (auto r_t : tmp_r) {
        if (l_t.first <= r_t.second && l_t.second >= r_t.first) {
          tmp.push_back(make_pair(max(l_t.first, r_t.first),
                                  min(l_t.second, r_t.second)));
        }
      }
    }
    return tmp;
  }

  bool CanThrough() {
    auto pre = wheel.front().second;

    for (int i = 1; i < wheel.size(); i++) {
      caf r = wheel[i].second;
      pre = Join(pre, r);
      if (pre.empty())
        return false;
    }

    return (!pre.empty());
  }

  auto Solution(const vector<string> &lines) {
    // get wheel
    rep(i, 0, lines.size()) {
      auto tmp = SplitLine<int>(lines[i]);
      wheel.push_back(make_pair(tmp[0], vector<pii>{}));
      for (int i = 2; i < tmp.size(); i += 2) {
        wheel.back().second.push_back(
            make_pair(tmp[i], (tmp[i] + tmp[i + 1]) % 360));
      }
    }

    // calculate for every second
    for (int i = 0; i <= 360; i++) {
      if (CanThrough()) {
        return vector<string>{to_string(i)};
      }

      for (auto &v : wheel) {
        auto sp = v.first;
        for (auto &w : v.second) {
          w.first = (sp + w.first) % 360;
          w.second = (sp + w.second) % 360;
        }
      }
    }

    auto rst = vector<string>{"none"};
    return rst;
  }
};

struct ratios {
  vi goal;
  vvi type;

  vector<int> rst;
  vi Cal(int x, int y, int z) {
    vi tmp;
    for (int i = 0; i < type.front().size(); i++) {
      tmp.push_back(x * type[0][i] + y * type[1][i] + z * type[2][i]);
    }
    return tmp;
  }

  bool Check(vi cal, int &multi) {
    for (int i = 0; i < goal.size(); i++) {
      if (goal[i] == 0 && cal[i] == 0)
        continue;

      if (goal[i] == 0 || cal[i] == 0)
        return false;

      if (goal[i] != 0 && cal[i] != 0 && (cal[i] % goal[i]) != 0) {
        return false;
      }
    }

    if (cal[0] == 0) {
      multi = 0;
    } else {
      multi = cal[0] / goal[0];
    }
    for (int i = 1; i < goal.size(); i++) {
      auto tmp = 0;

      if (cal[i] != 0)
        tmp = cal[i] / goal[i];

      if (multi == 0 || tmp == 0) {
        multi = max(tmp, multi);
      } else {
        if (tmp != multi) {
          return false;
        }
      }
    }
    return true;
  }

  vvvb visited;

  void dfs(int x, int y, int z) {
    if (x >= 100 || y >= 100 || z >= 100) {
      return;
    }
    if (visited[x][y][z])
      return;

    visited[x][y][z] = true;

    auto cal = Cal(x, y, z);
    auto multi = 1;
    if (Check(cal, multi)) {
      if (rst.back() > multi) {
        rst.back() = multi;
        rst[0] = x;
        rst[1] = y;
        rst[2] = z;
      }
    }

    dfs(x + 1, y, z);
    dfs(x, y + 1, z);
    dfs(x, y, z + 1);
  }

  void dp_table() {
    for (int x = 1; x < 100; x++) {
      for (int y = 1; y < 100; y++) {
        for (int z = 1; z < 100; z++) {
          auto cal = Cal(x, y, z);
          auto multi = 1;
          if (Check(cal, multi)) {
            if (rst.back() > multi) {
              rst.back() = multi;
              rst[0] = x;
              rst[1] = y;
              rst[2] = z;
            }
          }
        }
      }
    }
  }

  auto Solution(const vector<string> &lines) {
    rst = vector<int>(4, numeric_limits<int>::max());
    visited = vvvb(100, vvb(100, vb(100, false)));

    goal = SplitLine<int>(lines.front());
    for (int i = 1; i < lines.size(); i++) {
      type.push_back(SplitLine<int>(lines[i]));
    }
    // cout << goal << endl;
    // for (auto &t : type)
    //   cout << t << endl;

    dfs(1, 0, 0);
    dfs(0, 1, 0);
    dfs(0, 0, 1);

    // dp_table();
    if (rst.back() == numeric_limits<int>::max()) {
      return vector<string>{"NONE"};
    } else {
      auto tmp = vector<string>{};
      for (auto &v : rst) {
        tmp.push_back(to_string(v));
      }
      return vector<string>{Join(tmp, " ")};
    }
  }
};
struct msquare {
  string target;
  string start = string{"12345678"};
  unordered_map<char, int> cur_map;
  unordered_map<char, int> t_map;

  int rst_n = numeric_limits<int>::max();
  string rst_s = "";
  string OP_A(string cur) {
    int i = 0;
    int j = cur.size() - 1;
    while (i < j) {
      swap(cur[i], cur[j]);
      i++;
      j--;
    }
    return cur;
  }
  string OP_B(string cur) {
    auto l = string(1, cur[3]);
    auto r = string(1, cur[4]);
    cur.erase(next(cur.begin(), 3));
    cur.erase(next(cur.begin(), 3));
    return l + cur + r;
  }

  string OP_C(string cur) {
    swap(cur[1], cur[2]);
    swap(cur[1], cur[5]);
    swap(cur[1], cur[6]);
    return cur;
  }

  unordered_set<string> visited;

  void Print(const string &cur, const string &path) {
    cout << endl;
    for (int i = 0; i < 4; i++) {
      cout << cur[i];
    }
    cout << endl;

    for (int i = 7; i >= 4; i--) {
      cout << cur[i];
    }
    cout << endl;

    cout << path << endl;
  }

  void dfs(int n_A, int n_B, string cur, string path) {
    Print(cur, path);
    UpdateMap(cur_map, cur);

    if (path.size() > rst_n)
      return;

    if (visited.find(cur) != visited.end()) {
      return;
    }

    if (cur == target) {
      if (path.size() < rst_n) {
        rst_s = path;
        rst_n = path.size();
      } else if (path.size() == rst_n && path < rst_s) {
        rst_s = path;
      }
      return;
    }

    visited.insert(cur);

    auto good_up_down = all_of(cur.begin(), cur.end(), [this, &cur](char c) {
      return CheckUpDown(c, cur);
    });

    auto good_in_row = all_of(cur.begin(), cur.end(),
                              [this](char c) { return CheckInRightRow(c); });

    if (good_in_row && good_up_down) {
      dfs(n_A, n_B + 1, OP_B(cur), path + "B");
    } else {
      if (!good_up_down) {
        if ((!CheckUpDown(cur[1], cur)) && (!CheckUpDown(cur[2], cur))) {
          dfs(n_A, n_B, OP_C(cur), path + "C");
        } else {
          dfs(n_A, n_B + 1, OP_B(cur), path + "B");
        }
      }

      if (!good_in_row) {
        if (!CheckInRightRow('1')) {
          dfs(n_A + 1, n_B, OP_A(cur), path + "A");
        } else {
          if (CheckInRightRow(cur[1]) || CheckInRightRow(cur[2])) {
            dfs(n_A, n_B + 1, OP_B(cur), path + "B");
          } else {
            dfs(n_A, n_B, OP_C(OP_C(cur)), path + "CC");
          }
        }
      }
    }
  }

  void UpdateMap(unordered_map<char, int> &m, const string &s) {
    for (int i = 0; i < s.size(); i++) {
      m[s[i]] = i;
    }
  }

  bool CheckUpDown(char c, const string &cur) {
    return cur[7 - cur_map[c]] == target[7 - t_map[c]];
  }

  bool CheckInRightRow(char c) {
    if (cur_map[c] < 4) {
      return t_map[c] < 4;
    } else {
      return t_map[c] >= 4;
    }
  }

  string bfs() {
    auto mem = unordered_map<string, string>{};

    auto cmp = [](const pair<string, string> &a,
                  const pair<string, string> &b) {
      if (a.second.size() == b.second.size()) {
        return a > b;
      }
      return a.second.size() > b.second.size();
    };

    auto IsBetterPath = [&](const pair<string, string> &v) {
      if (mem.find(v.first) == mem.end())
        return true;

      if (mem[v.first].size() > v.second.size())
        return true;

      if (mem[v.first].size() == v.second.size() && mem[v.first] > v.second)
        return true;

      return false;
    };

    auto pq = priority_queue<pair<string, string>, vector<pair<string, string>>,
                             decltype(cmp)>(cmp);
    pq.push(make_pair(start, ""));

    while (!pq.empty()) {
      auto v = pq.top();
      pq.pop();
      if (IsBetterPath(v)) {
        mem[v.first] = v.second;

        auto a = make_pair(OP_A(v.first), v.second + "A");
        auto b = make_pair(OP_B(v.first), v.second + "B");
        auto c = make_pair(OP_C(v.first), v.second + "C");
        for (auto &t : {a, b, c}) {
          if (IsBetterPath(t)) {
            pq.push(t);
          }
        }
      }
    }

    return mem[target];
  }

  auto Solution(const vector<string> &lines) {
    auto tmp = SplitLine<int>(lines.front());
    for (auto &v : tmp) {
      target += to_string(v);
    }

    // UpdateMap(t_map, target);
    // visited.clear();
    // dfs(1, 0, OP_A(start), "A");
    // visited.clear();
    // dfs(0, 0, start, "");

    rst_s = bfs();
    auto rst = vector<string>{to_string(rst_s.size()), rst_s};
    return rst;
  }
};

struct butter {
  int N, P, C;
  vector<int> cows;
  vector<vector<pll>> path;

  vector<int> GetShortPath(int cow) {
    auto cow_path = vector<int>(P, numeric_limits<int>::max());
    // pair<dist, pasture>
    auto pq = priority_queue<pll, vector<pll>, greater<>>{};
    pq.push(make_pair(0, cows[cow]));

    while (!pq.empty()) {
      auto t = pq.top();
      pq.pop();

      if (cow_path[t.second] <= t.first)
        continue;

      cow_path[t.second] = t.first;

      caf pasture_path = path[t.second];
      for (auto &v : pasture_path) {
        auto tmp = v.second + t.first;
        if (tmp < cow_path[v.first]) {
          pq.push(make_pair(tmp, v.first));
        }
      }
    }
    return cow_path;
  }

  void dfs(int cur_p, int dst, vector<int> &cow_path) {
    if (cow_path[cur_p] <= dst)
      return;

    cow_path[cur_p] = dst;
    for (auto &v : path[cur_p]) {
      auto tmp = dst + v.second;
      if (tmp < cow_path[v.first]) {
        dfs(v.first, tmp, cow_path);
      }
    }
  }

  auto Solution(const vector<string> &lines) {
    auto tmp = SplitLine<int>(lines.front());
    N = tmp[0];
    P = tmp[1] + 1;
    C = tmp[2];

    for (int i = 1; i <= N; i++) {
      cows.push_back(stoi(lines[i]));
    }

    path = vector<vector<pll>>(P);
    for (int i = N + 1; i < lines.size(); i++) {
      auto tmp = SplitLine<int>(lines[i]);
      path[tmp[0]].push_back(make_pair(tmp[1], tmp[2]));
      path[tmp[1]].push_back(make_pair(tmp[0], tmp[2]));
    }

    auto mini_paths = vvi{};
    for (int i = 0; i < cows.size(); i++) {
      mini_paths.push_back(GetShortPath(i));

      // dfs
      // auto tmp = vector<int>(P, numeric_limits<int>::max());
      // dfs(cows[i], 0, tmp);
      // mini_paths.push_back(tmp);
    }

    auto rst = numeric_limits<ll>::max();
    for (int i = 1; i < mini_paths.front().size(); i++) {
      ll tmp = 0;
      for (int j = 0; j < mini_paths.size(); j++) {
        tmp += mini_paths[j][i];
      }
      rst = min(rst, tmp);
    }

    return vector<ll>{rst};
  }
};

struct fence {
  map<int, map<int, int>> graph;
  map<int, int> degree;
  auto Solution(const vector<string> &lines) {
    auto rst = vector<int>{};

    for (int i = 1; i < lines.size(); i++) {
      auto tmp = SplitLine<int>(lines[i]);
      graph[tmp[0]][tmp[1]]++;
      degree[tmp[0]]++;
      graph[tmp[1]][tmp[0]]++;
      degree[tmp[1]]++;
    }

    auto st = stack<int>{};
    auto start = degree.begin()->first;
    for (caf v : degree) {
      if (v.second % 2 != 0) {
        start = v.first;
        break;
      }
    }

    st.push(start);

    while (!st.empty()) {
      auto n = st.top();

      auto &cnt = graph[n];
      if (cnt.empty()) {
        rst.push_back(n);
        st.pop();
      } else {
        auto next = cnt.begin()->first;
        st.push(next);

        cnt[next]--;
        if (cnt[next] == 0)
          cnt.erase(next);

        graph[next][n]--;

        if (graph[next][n] == 0)
          graph[next].erase(n);
      }
    }

    reverse(all(rst));
    return rst;
  }
};

struct shopping {
  vi goods;
  vi buys;
  vi buys_id;
  vector<pair<int, unordered_map<int, int>>> offers;
  unordered_map<string, int> key_map;

  ll rst_v = numeric_limits<ll>::max();

  ll CalLeft(const vi cur) {
    ll tmp = 0;
    for (int i = 0; i < cur.size(); i++) {
      tmp += cur[i] * goods[buys_id[i]];
    }
    return tmp;

    return 0;
  }

  void CalLeft(ll cur_v) {
    auto tmp = cur_v;
    for (int i = 0; i < buys.size(); i++) {
      tmp += buys[i] * goods[i];
    }

    rst_v = min(rst_v, tmp);
  };

  void dfs(int start, ll cur_v) {
    if (all_of(all(buys), [](auto &v) { return v == 0; })) {
      rst_v = min(rst_v, cur_v);
    }

    if (start >= offers.size()) {

      CalLeft(cur_v);
      return;
    }

    for (int i = start; i < offers.size(); i++) {
      if (all_of(all(offers[i].second),
                 [this](auto &v) { return buys[v.first] >= v.second; })) {
        for (caf v : offers[i].second)
          buys[v.first] -= v.second;

        dfs(i, cur_v + offers[i].first);

        for (caf v : offers[i].second)
          buys[v.first] += v.second;
      }
    }

    CalLeft(cur_v);
  }

  ll bfs() {
    struct Item {
      ll pay = 0;
      int state = 0;
      vi cur_item;
    };

    auto cmp = [](const Item &a, const Item &b) { return a.pay < b.pay; };

    auto pq = priority_queue<Item, vector<Item>, decltype(cmp)>(cmp);
    auto first = Item{};
    first.pay = 0;
    first.state = 0;
    first.cur_item = buys;

    pq.push(first);

    auto CalLeft = [this](const Item &cur) {
      ll tmp = cur.pay;
      for (int i = 0; i < cur.cur_item.size(); i++) {
        tmp += cur.cur_item[i] * goods[i];
      }
      return tmp;
    };

    auto GetKey = [this](const Item &cur) {
      auto tmp = int{0};
      for (int i = 0; i < cur.cur_item.size(); i++) {
        if (cur.cur_item[i] > 0) {
          tmp = tmp * 10 + cur.cur_item[i];
        }
      }
      return tmp;
    };

    auto mem = unordered_map<int, ll>{};

    while (!pq.empty()) {
      auto cur = pq.top();
      pq.pop();
      auto val = CalLeft(cur);
      auto key = cur.state;
      if (mem.find(key) != mem.end() && mem[key] < val)
        continue;

      mem[key] = val;

      for (int i = 0; i < offers.size(); i++) {
        if (all_of(all(offers[i].second), [this, &cur](auto &v) {
              return cur.cur_item[v.first] >= v.second;
            })) {

          auto tmp = cur;

          for (caf v : offers[i].second)
            tmp.cur_item[v.first] -= v.second;

          tmp.pay += offers[i].first;

          tmp.state = GetKey(tmp);
          auto tmp_v = CalLeft(tmp);
          if (mem.find(tmp.state) != mem.end() && mem[tmp.state] < tmp_v)
            continue;

          pq.push(tmp);
        }
      }
    }

    for (caf v : mem) {
      rst_v = min(rst_v, v.second);
    }
    return rst_v;
  }

  ll dp_table() {
    ll dp[6][6][6][6][6][100];

    auto mem = unordered_map<int, int>{};
    auto default_v = numeric_limits<int>::max();

    for (int a1 = 0; a1 <= buys[buys_id[0]]; a1++) {
      for (int a2 = 0; a2 <= buys[buys_id[1]]; a2++) {
        for (int a3 = 0; a3 <= buys[buys_id[2]]; a3++) {
          for (int a4 = 0; a4 <= buys[buys_id[3]]; a4++) {
            for (int a5 = 0; a5 <= buys[buys_id[4]]; a5++) {
              mem[buys_id[0]] = a1;
              mem[buys_id[1]] = a2;
              mem[buys_id[2]] = a3;
              mem[buys_id[3]] = a4;
              mem[buys_id[4]] = a5;

              for (int r = 0; r <= offers.size(); r++) {
                // cout << a1 << a2 << a3 << a4 << a5 << endl;

                if (r == 0) {
                  auto tmp = CalLeft(vi{a1, a2, a3, a4, a5});
                  dp[a1][a2][a3][a4][a5][r] = tmp > 0 ? tmp : default_v;
                  // cout << "default=" << dp[a1][a2][a3][a4][a5][r] << endl;
                  continue;
                }

                auto tmp = dp[a1][a2][a3][a4][a5][r - 1];
                // cout << "r-1 tmp=" << tmp << endl;
                caf offer = offers[r - 1];
                if (all_of(all(offer.second), [this, &mem](auto &v) {
                      return mem[v.first] >= v.second;
                    })) {

                  // cout << "loop:" << mem << endl;
                  for (caf v : offer.second)
                    mem[v.first] -= v.second;

                  auto pre =
                      dp[mem[buys_id[0]]][mem[buys_id[1]]][mem[buys_id[2]]]
                        [mem[buys_id[3]]][mem[buys_id[4]]][r];

                  if (pre != default_v) {
                    tmp = min(tmp, pre + offer.first);
                  } else {
                    tmp = min(tmp, ll(offer.first));
                  }

                  // cout << "use offer =" << r << " , tmp=" << tmp << endl;
                  for (caf v : offer.second)
                    mem[v.first] += v.second;
                }

                dp[a1][a2][a3][a4][a5][r] = tmp;
              }
            }
          }
        }
      }
    }

    auto rst = dp[buys[buys_id[0]]][buys[buys_id[1]]][buys[buys_id[2]]]
                 [buys[buys_id[3]]][buys[buys_id[4]]][offers.size()];
    return rst == numeric_limits<int>::max() ? 0 : rst;
  }

  auto Solution(const vector<string> &lines) {
    int N = stoi(lines[0]);
    int i = 1;
    for (; i <= N; i++) {
      offers.emplace_back();
      caf l = SplitLine<int>(lines[i]);

      offers.back().first = l.back();
      for (int j = 1; j <= 2 * l.front(); j += 2) {
        offers.back().second[l[j]] = l[j + 1];
      }
    }
    i++;
    buys = vi(1000, 0);
    goods = vi(1000, 0);
    for (; i < lines.size(); i++) {
      caf l = SplitLine<int>(lines[i]);
      buys[l.front()] = l[1];
      goods[l.front()] = l[2];
      buys_id.push_back(l.front());
    }

    while (buys_id.size() < 5) {
      buys_id.push_back(0);
    }

    // dfs(0, 0);
    // bfs();
    rst_v = dp_table();

    return vector<ll>{rst_v};
  }
};

struct camelot {
  int R, C;
  vector<pii> knights;
  pii king;
  vvi king_path;

  vector<pii> Next(pii start) {
    auto i = start.first;
    auto j = start.second;
    return {{i - 1, j - 2}, {i - 2, j - 1}, {i - 1, j + 2}, {i - 2, j + 1},
            {i + 1, j - 2}, {i + 2, j - 1}, {i + 1, j + 2}, {i + 2, j + 1}};
  }

  vector<pii> KingNext(pii start) {
    auto i = start.first;
    auto j = start.second;
    return {{i, j - 1},     {i, j + 1},     {i - 1, j - 1}, {i - 1, j},
            {i - 1, j + 1}, {i + 1, j - 1}, {i + 1, j},     {i + 1, j + 1}};
  }

  vvi king_dfs(pii start) {
    auto map = vvi(R, vi(C, imax));

    struct Item {
      int step;
      pii point;
    };

    auto cmp = [](const Item &a, const Item &b) { return a.step > b.step; };
    auto pq = priority_queue<Item, vector<Item>, decltype(cmp)>(cmp);
    auto tmp = Item{};
    tmp.step = 0;
    tmp.point = start;
    pq.push(tmp);

    while (!pq.empty()) {
      auto t = pq.top();
      pq.pop();

      if (map[t.point.first][t.point.second] <= t.step)
        continue;

      map[t.point.first][t.point.second] = t.step;

      auto next_p = KingNext(t.point);
      for (caf n : next_p) {
        if (n.first < 0 || n.first >= R || n.second < 0 || n.second >= C) {
          continue;
        }
        auto v = t.step + 1;
        if (v < map[n.first][n.second]) {
          auto tmp = Item{};
          tmp.step = t.step + 1;
          tmp.point = n;
          pq.push(tmp);
        }
      }
    }

    return map;
  }

  using PATH = pair<vvi, vvi>;
  vector<PATH> paths;

  PATH dfs(pii start) {
    auto map = vvi(R, vi(C, imax));
    auto map_path = vvi(R, vi(C, imax));

    struct Item {
      int step;
      int path;
      pii point;
    };

    auto cmp = [](const Item &a, const Item &b) { return a.step > b.step; };
    auto pq = priority_queue<Item, vector<Item>, decltype(cmp)>(cmp);
    auto tmp = Item{};
    tmp.step = 0;
    tmp.path = king_path[start.first][start.second];
    tmp.point = start;
    pq.push(tmp);

    while (!pq.empty()) {
      auto t = pq.top();
      pq.pop();

      if (map[t.point.first][t.point.second] < t.step)
        continue;

      if (map[t.point.first][t.point.second] == t.step) {

        map_path[t.point.first][t.point.second] =
            min(map_path[t.point.first][t.point.second], t.path);

        continue;
      }

      map[t.point.first][t.point.second] = t.step;
      map_path[t.point.first][t.point.second] = t.path;

      auto next_p = Next(t.point);
      for (caf n : next_p) {
        if (n.first < 0 || n.first >= R || n.second < 0 || n.second >= C) {
          continue;
        }
        auto v = t.step + 1;
        if (v <= map[n.first][n.second]) {
          auto tmp = Item{};
          tmp.step = v;
          tmp.point = n;
          tmp.path = min(t.path, king_path[n.first][n.second]);
          pq.push(tmp);
        }
      }
    }

    return {map, map_path};
  }

  auto Solution(const vector<string> &lines) {
    // get positions
    auto tmp = SplitLine<int>(lines[0]);
    R = tmp.front();
    C = tmp.back();
    auto k = SplitLine<string>(lines[1]);
    king.first = stoi(k.back()) - 1;
    king.second = k.front()[0] - 'A';

    for (int i = 2; i < lines.size(); i++) {
      auto tmp = SplitLine<string>(lines[i]);
      for (int j = 0; j < tmp.size(); j += 2) {
        auto v = make_pair(stoi(tmp[j + 1]) - 1, tmp[j][0] - 'A');
        knights.push_back(v);
      }
    }

    if (knights.empty())
      return vector<int>{0};

    king_path = king_dfs(king);

    for (caf k : knights) {
      paths.push_back(dfs(k));
    }

    auto map = vvi(R, vi(C, -1));
    auto min_map = vvi(R, vi(C, imax));
    rep(r, 0, R) {
      rep(c, 0, C) {
        auto tmp = 0;
        auto good_to_go = true;
        rep(k, 0, knights.size()) {
          if (paths[k].first[r][c] == imax) {
            good_to_go = false;
            break;
          }
          tmp += paths[k].first[r][c];
        }
        if (good_to_go) {
          map[r][c] = tmp;
        }
      }
    }

    auto rst = numeric_limits<int>::max();
    rep(r, 0, R) {
      rep(c, 0, C) {
        if (map[r][c] == -1)
          continue;

        if (rst <= map[r][c])
          continue;

        auto tmp = king_path[r][c];

        rep(k, 0, knights.size()) { tmp = min(tmp, paths[k].second[r][c]); }

        rst = min(rst, tmp + map[r][c]);
      }
    }

    return vector<int>{rst};
  }
};

struct range {
  vector<string> graph;
  using Item = pair<pii, int>;
  vector<pii> GetNext(const Item &v) {
    auto rst = vector<pii>{};
    auto n_i = v.first.first + v.second;
    auto n_j = v.first.second + v.second;

    rep(i, v.first.first, n_i) { rst.push_back({i, n_j}); }

    rep(j, v.first.second, n_j) { rst.push_back({n_i, j}); }

    rst.push_back({n_i, n_j});
    return rst;
  }

  int GetMaxRange(Item t) {
    auto next = GetNext(t);
    while (!next.empty() && all_of(all(next), [this](const pii &v) {
      if (v.first < 0 || v.first >= graph.size() || v.second < 0 ||
          v.second >= graph.front().size())
        return false;
      return graph[v.first][v.second] == '1';
    })) {
      t.second++;
      next = GetNext(t);
    }
    return t.second - 1;
  };

  map<int, int> dp_table() {
    auto dp = vvi(graph.size(), vi(graph.front().size(), 0));
    auto rst = map<int, int>{};
    rep(i, 0, graph.size()) {
      rep(j, 0, graph.front().size()) {
        if (graph[i][j] == '0')
          continue;

        auto t = 1;
        if (i != 0) {
          t = max(t, dp[i - 1][j]);
        }
        if (j != 0) {
          t = max(t, dp[i][j - 1]);
        }
        dp[i][j] = GetMaxRange({{i, j}, t});
      }
    }

    rep(i, 0, graph.size()) {
      rep(j, 0, graph.front().size()) {
        for (int k = 1; k <= dp[i][j]; k++) {
          rst[k + 1]++;
        }
      }
    }
    return rst;
  }

  auto Solution(const vector<string> &lines) {
    for (int i = 1; i < lines.size(); i++) {
      graph.push_back(lines[i]);
    }
    auto tmp_mem = dp_table();
    auto rst = vector<string>{};
    for (caf v : tmp_mem) {
      rst.push_back(to_string(v.first) + " " + to_string(v.second));
    }
    return rst;
  }
};

struct game1 {
  vi nums;
  vi acc_nums;
  vvvi mem;
  int N;

  int GetSum(int i, int j) {
    if (i == 0)
      return acc_nums[j];
    return acc_nums[j] - acc_nums[i - 1];
  }

  int dp(int x, int i, int j) {
    if (x == 0)
      return 0;

    if (i == j)
      return nums[i];

    if (mem[x][i][j] != -1)
      return mem[x][i][j];

    auto v = max(nums[i] + GetSum(i + 1, j) - dp(x - 1, i + 1, j),
                 nums[j] + GetSum(i, j - 1) - dp(x - 1, i, j - 1));

    mem[x][i][j] = v;
    return v;
  }

  int dp_table() {
    auto dp_t = vector<vector<int>>(N + 1, vector<int>(N, 0));

    for (int x = 1; x <= N; x++) {
      for (int i = 0; i + x - 1 < N; i++) {
        auto r = i + x - 1;
        if (x == 1) {
          dp_t[x][i] = nums[i];
        } else {
          dp_t[x][i] = max(nums[i] + GetSum(i + 1, r) - dp_t[x - 1][i + 1],
                           nums[r] + GetSum(i, r - 1) - dp_t[x - 1][i]);
        }
      }
    }

    return dp_t[N][0];
  }

  auto Solution(const vector<string> &lines) {
    for (int i = 1; i < lines.size(); i++) {
      for (caf v : SplitLine<int>(lines[i])) {
        nums.push_back(v);
      }
    }
    N = nums.size();

    acc_nums.push_back(nums.front());
    for (int i = 1; i < nums.size(); i++) {
      acc_nums.push_back(acc_nums.back() + nums[i]);
    }

    // dp function
    //  mem = vvvi(N + 1, vvi(N + 1, vi(N + 1, -1)));
    //  auto A = dp(nums.size(), 0, nums.size() - 1);

    // dp_table
    auto A = dp_table();

    auto rst = vector<string>{to_string(A) + " " +
                              to_string(GetSum(0, nums.size() - 1) - A)};
    return rst;
  }
};

struct heritage {
  template <typename T> struct TreeNode {
    T val;
    TreeNode *left;
    TreeNode *right;

    TreeNode() : val(T{}), left(nullptr), right(nullptr) {}

    TreeNode(T x) : val(x), left(nullptr), right(nullptr) {}

    TreeNode(T x, TreeNode *left, TreeNode *right)
        : val(x), left(left), right(right) {}
  };

  string IN;
  string PRE;
  TreeNode<char> root;

  auto build(int in_l, int in_r, int pre_l, int pre_r) -> TreeNode<char> * {
    if (pre_l == pre_r)
      return nullptr;

    auto r = PRE[pre_l];
    auto c_r = new TreeNode(r);
    auto dst = distance(IN.begin() + in_l,
                        find(IN.begin() + in_l, IN.begin() + in_r, r));
    c_r->left = build(in_l, in_l + dst, pre_l + 1, pre_l + dst + 1);
    c_r->right = build(in_l + dst + 1, in_r, pre_l + dst + 1, pre_r);
    c_r->val = r;
    return c_r;
  }

  string PostOrder(TreeNode<char> *root) {
    if (root == nullptr)
      return "";
    auto rst = string{};
    rst += PostOrder(root->left);
    rst += PostOrder(root->right);
    rst += string(1, root->val);
    return rst;
  }

  auto Solution(const vector<string> &lines) {
    IN = lines.front();
    PRE = lines.back();
    root = *build(0, IN.size(), 0, PRE.size());

    auto rst = vector<string>{PostOrder(&root)};
    return rst;
  }
};

struct fence9 {
  auto InTriangleV1(pii v, pii v0, pii v1, pii v2) -> bool {
    auto det = [](pii v1, pii v2) -> int {
      return v1.first * v2.second - v1.second * v2.first;
    };

    auto a = (det(v, v2) - det(v0, v2)) / double(det(v1, v2));
    auto b = -(det(v, v1) - det(v0, v1)) / double(det(v1, v2));
    if (a > 0 && b > 0 && a + b < 1)
      return true;

    return false;
  }

  auto inTriangle(pii v, pii v0, pii v1, pii v2) -> bool {
    // do not finish
    auto CrossProduct = [](pii a, pii b) {
      return a.first * b.second - a.second * b.first;
    };
    auto SameSide = [](pii p1, pii p2, pii a, pii b) { return false; };

    return false;
  }

  auto Solution(const vector<string> &lines) {
    auto tmp = SplitLine<int>(lines.front());
    auto v0 = make_pair(0, 0);
    auto v1 = make_pair(tmp[0], tmp[1]);
    auto v2 = make_pair(tmp[2], 0);
    double a1 = double(v1.second) / v1.first;
    double a2 = double(v1.second) / (v1.first - v2.first);
    double b2 = double(v1.second * v2.first) / (v2.first - v1.first);

    auto acc = int{0};

    if (v1.first <= v2.first) {
      for (int i = 1; i < v2.first; i++) {
        auto v = int{0};
        if ((i <= v1.first)) {
          v = int(a1 * i);
        } else {
          v = int(a2 * i + b2);
        }

        if (v <= 0)
          continue;
        if (InTriangleV1(make_pair(i, v), v0, v1, v2)) {
          acc += v;
        } else {
          acc += v - 1;
        }
      }
    } else {
      for (int i = 1; i < v1.first; i++) {
        if ((i < v2.first)) {
          // acc += int(ceil(a1 * i - 1));
          auto v = int(a1 * i);
          if (v <= 0)
            continue;
          if (InTriangleV1(make_pair(i, v), v0, v1, v2)) {
            acc += v;
          } else {
            acc += v - 1;
          }
        } else
          for (int y = int(a2 * i + b2); y <= int(a1 * i); y++) {
            if (y <= 0)
              continue;
            if (InTriangleV1(make_pair(i, y), v0, v1, v2)) {
              acc++;
            }
          }
      }
    }

    return vector<int>{acc};
  }
};

struct rockers {
  vvvi dp;

  auto Solution(const vector<string> &lines) {
    auto tmp = SplitLine<int>(lines.front());
    auto N = tmp[0];
    auto T = tmp[1];
    auto M = tmp[2];
    auto Len = SplitLine<int>(lines.back());
    dp = vvvi(N + 1, vvi(M + 1, vi(T + 1, 0)));
    for (int i = 0; i <= N; i++) {
      for (int j = 0; j <= M; j++) {
        for (int k = 0; k <= T; k++) {
          if (i == 0 || j == 0) {
            dp[i][j][k] = 0;
            continue;
          }
          auto l = Len[i - 1];

          dp[i][j][k] = dp[i - 1][j][k];
          if (l <= k) {
            dp[i][j][k] = max(dp[i][j][k], 1 + dp[i - 1][j][k - l]);
          } else if (j > 1 && l <= T) {
            dp[i][j][k] = max(dp[i][j][k], 1 + dp[i - 1][j - 1][T - l]);
          }
        }
      }
    }

    auto rst = vector<int>{dp[N][M][T]};
    return rst;
  }
};

struct nuggets {
  vector<int> nums;
  vector<vector<bool>> nums_dp;
  const int LIMIT = 65024;

  auto Solution(const vector<string> &lines) {
    rep(i, 1, lines.size()) { nums.push_back(stoi(lines[i])); }
    // cout << nums << endl;
    nums_dp =
        vector<vector<bool>>(nums.size(), vector<bool>(LIMIT + 1000, false));
    for (int i = 0; i < nums_dp.size(); i++) {
      caf v = nums[i];
      for (int j = 0; j < nums_dp.front().size(); j++) {
        if (j == v) {
          nums_dp[i][j] = true;
          continue;
        } else {
          if (j >= v) {
            nums_dp[i][j] = nums_dp[i][j] || nums_dp[i][j - v];
          }
          if (i > 0) {
            nums_dp[i][j] = nums_dp[i][j] || nums_dp[i - 1][j];
          }
        }
      }
    }
    auto last = int{};
    for (int j = 0; j < nums_dp.front().size(); j++) {
      if (nums_dp.back()[j] == false) {
        last = j;
      }
    }
    last = last > LIMIT ? 0 : last;
    auto rst = vector<int>{last};
    return rst;
  }
};

struct fence6 {
  vector<unordered_map<int, int>> graph{};
  vector<int> lens{};

  auto Solution(const vector<string> &lines) {
    auto N = stoi(lines.front());

    graph = vector<unordered_map<int, int>>(N);
    lens = vi(N);

    for (int i = 1; i < lines.size(); i += 3) {
      auto tmp = SplitLine<int>(lines[i]);
      auto index = tmp[0] - 1;
      lens[index] = tmp[1];
      for (caf v : SplitLine<int>(lines[i + 1])) {
        graph[index][v - 1] = 0;
      };

      for (caf v : SplitLine<int>(lines[i + 2])) {
        graph[index][v - 1] = 1;
      };
    }

    auto rst = imax;

    rep(i, 0, graph.size()) {
      auto visited = vi(graph.size(), imax);

      auto pq =
          priority_queue<pair<int, pii>, vector<pair<int, pii>>, greater<>>{};

      for (caf l : graph[i]) {
        pq.push(mp(lens[l.f], mp(i, l.f)));

        while (!pq.empty()) {
          auto t = pq.top();
          pq.pop();

          if (t.s.s == i)
            rst = min(rst, t.first);
          else {
            if (visited[t.s.s] <= t.f)
              continue;
            visited[t.s.s] = t.f;

            auto next_i = graph[t.s.s][t.s.f] == 0 ? 1 : 0;
            for (caf v : graph[t.s.s]) {
              if (v.s != next_i)
                continue;

              auto tmp = t.first + lens[v.f];

              if (tmp >= visited[v.f] || tmp >= rst)
                continue;

              pq.push(mp(tmp, mp(t.s.s, v.f)));
            }
          }
        }
      }
    }
    return vector<int>{rst};
  }
};

/*
ID: racer
TASK: hamming
LANG: C++
*/

#define TASK "hamming"
#define INFile TASK ".in"
#define OUTFile TASK ".out"

int main() {
  auto lines = vector<string>{};

#ifdef __clang__
  lines = ReadInput();
#else
  lines = ReadLines(INFile);
#endif

  auto rst = hamming().Solution(lines);

#ifdef __clang__
  cout << "result=" << endl;
  for (auto &v : rst) {
    cout << v << endl;
  }
#else
  WriteLine(rst, OUTFile);
#endif
  return 0;
}
