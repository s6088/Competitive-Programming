//dsu
int par(int x)
{
  if (p[x] != x)
    return p[x] = par(p[x]);
  return p[x];
}
void dsu(int x, int y)
{
  int a = par(x), b = par(y);
  if (a == b)
    return;
  if (z[a] < z[b])
    swap(a, b);
  z[a] += z[b];
  p[b] = a;
}

//trie
int newWord(string s)
{
  int index = 0;
  for (int i = 0; i < s.size(); i++)
  {
    int x = s[i] - 'a';
    if (!trie[index][x])
      trie[index][x] = ++cntNode;
    index = trie[index][x];
  }
  isEnd[index] = true;
  return index;
}

bool isPresent(string s)
{
  int index = 0;
  for (int i = 0; i < s.size(); i++)
  {
    int x = s[i] - 'a';
    if (!trie[index][x])
      return false;
    index = trie[index][x];
  }
  return isEnd[index];
}

void newWord(int n)
{
  int index = 0;
  for (int i = 29; i >= 0; i--)
  {
    int x = (n & (1 << i));
    if (x)
      x = 1;
    if (!trie[index][x])
      trie[index][x] = ++cntNode;
    index = trie[index][x];
    cnt[index]++;
  }
}

int mn(int n)
{
  int index = 0, ans = 0;
  for (int i = 29; i >= 0; i--)
  {
    int x = (n & (1 << i));
    if (x)
      x = 1;
    if (cnt[trie[index][x]])
    {
      cnt[trie[index][x]]--;
      index = trie[index][x];
    }
    else
    {
      ans |= (1 << i);
      cnt[trie[index][!x]]--;
      index = trie[index][!x];
    }
  }
  return ans;
}

//Segment Tree

size of st 4n before built st[] << -inf int h = ceil(log2((double)n));

//O(n)
void built(int l, int h, int p)
{
  if (l > h)
    return;
  if (l == h)
  {
    st[p] = a[l];
    return;
  }
  int m = (l + h) / 2;
  built(l, m, 2 * p);
  built(m + 1, h, 2 * p + 1);
  st[p] = min(st[2 * p], st[2 * p + 1]);
}

//O(logn)
int rQ(int ql, int qh, int l, int h, int p)
{
  if (l > h || ql > h || qh < l)
    return inf; //0 in sum
  if (lazy[p])
  {
    st[p] += lazy[p]; ///sum->(h-l+1)*lazy[p]
    if (l != h)
      lazy[2 * p] += lazy[p], lazy[2 * p + 1] += lazy[p];
    lazy[p] = 0;
  }
  if (ql <= l && qh >= h)
    return st[p];
  int m = (l + h) / 2;
  return min(rQ(ql, qh, l, m, 2 * p), rQ(ql, qh, m + 1, h, 2 * p + 1));
}

//init lazy <- 0

void update(int ql, int qh, int v, int l, int h, int p)
{

  if (lazy[p])
  {
    st[p] += lazy[p]; ///sum->(h-l+1)*lazy[p]
    if (l != h)
      lazy[2 * p] += lazy[p], lazy[2 * p + 1] += lazy[p];
    lazy[p] = 0;
  }

  if (l > h || l > qh || h < ql)
    return;

  if (ql <= l && qh >= h)
  {
    st[p] += v; ///(h-l+1)*v
    if (l != h)
      lazy[2 * p] += v, lazy[2 * p + 1] += v;
    return;
  }

  int m = (l + h) / 2;
  update(ql, qh, v, l, m, 2 * p);
  update(ql, qh, v, m + 1, h, 2 * p + 1);
  st[p] = min(st[2 * p], st[2 * p + 1]);
}

///lca root change each query
int lca(int r, int u, int v)
{
  if (r == 0)
    return 0;
  if (r == u || r == v)
    return r;
  int z = 0, count = 0;
  for (int i = 0; i < graph[r].size(); i++)
  {
    int l = lca(graph[r][i], u, v);
    if (l)
      count++, z = l;
  }
  if (count == 2)
    return r;
  return z;
}
//binary raised lca
void dfs(int u, int a)
{
  l[u] = l[a] + 1;
  P[u][0] = a;
  for (int i = 0; i < graph[u].size(); i++)
  {
    int v = graph[u][i];
    if (l[v] == -1)
      dfs(v, u);
  }
}

void built(int s, int n)
{
  dfs(s, 0);
  for (int j = 1; j <= 30; j++)
    for (int i = 1; i <= n; i++)
      if (P[i][j - 1] != -1)
        P[i][j] = P[P[i][j - 1]][j - 1];
}

int lca(int u, int v)
{

  if (l[u] < l[v])
    swap(u, v);

  for (int i = 30; i >= 0; i--)
    if (l[u] - l[v] >= (1 << i))
      u = P[u][i];

  if (u == v)
    return u;

  for (int i = 30; i >= 0; i--)
    if (P[u][i] != -1 && P[u][i] != P[v][i])
      u = P[u][i], v = P[v][i];

  return P[u][0];
}

//bit 1D
int prfSum(int i)
{
  return i == 0 ? 0 : bit[i] + prfSum(i - (-i & i));
}
void update(int n, int x, ll v)
{
  for (int i = x; i <= n; i += i & -i)
    bit[i] += v;
}

//bit 2D
void update(int n, int m, int x, int y, int v)
{
  for (int i = x; i <= n; i += i & -i)
    for (int j = y; j <= m; j += j & -j)
      BIT[i][j] += v;
}

int prfSum(int x, int y)
{
  int sum = 0;
  for (int i = x; i > 0; i -= i & -i)
    for (int j = y; j > 0; j -= j & -j)
      sum += BIT[i][j];
  return sum;
}

int prfSum(int r1, int c1, int r2, int c2)
{
  return prfSum(r2, c2) - prfSum(r2, c1 - 1) - prfSum(r1 - 1, c2) + prfSum(r1 - 1, c1 - 1);
}

void built(int n, int m)
{
  memset(BIT, 0, sizeof BIT);
  for (int i = 1; i <= n; i++)
    for (int j = 1; j <= m; j++)
      BIT[i][j] = AUX[i][j] - prfSum(1, 1, i, j);
}

//sparse table 1D 1 indexing no cleaning req
void built(int n)
{
  for (int i = 1; i <= n; i++)
    S[i][0] = i;
  for (int j = 1; (1 << j) <= n; j++)
  {
    for (int i = 1; i + (1 << j) - 1 <= n; i++)
    {
      if (a[S[i][j - 1]] < a[S[i + (1 << (j - 1))][j - 1]])
        S[i][j] = S[i][j - 1];
      else
        S[i][j] = S[i + (1 << (j - 1))][j - 1];
    }
  }
}

int rQ(int l, int r)
{
  int len = r - l + 1;
  int k = log2(len);
  return min(a[S[l][k]], a[S[l + len - (1 << k)][k]]);
}

///2D rmq
///0 index + S[LG][N][LG][N]
void built(int n, int m)
{
  for (int ir = 0; ir < n; ir++)
  {
    for (int ic = 0; ic < m; ic++)
      S[0][ir][0][ic] = a[ir][ic];
    for (int jc = 1; jc <= log2(m) + 1; jc++)
      for (int ic = 0; ic + (1 << (jc - 1)) < m; ic++)
        S[0][ir][jc][ic] = max(S[0][ir][jc - 1][ic], S[0][ir][jc - 1][ic + (1 << (jc - 1))]);
  }
  for (int jr = 1; jr <= log2(n); jr++)
    for (int ir = 0; ir < n; ir++)
      for (int jc = 0; jc <= log2(m); jc++)
        for (int ic = 0; ic < m; ic++)
          S[jr][ir][jc][ic] = max(S[jr - 1][ir][jc][ic], S[jr - 1][ir + (1 << (jr - 1))][jc][ic]);
}

int rQ(int x1, int y1, int x2, int y2)
{

  int kx = log2(x2 - x1 + 1), ky = log2(y2 - y1 + 1), pkx = (1 << kx), pky = (1 << ky);
  int r1 = max(S[kx][x1][ky][y1], S[kx][x1][ky][y2 + 1 - pky]);
  int r2 = max(S[kx][x2 + 1 - pkx][ky][y1], S[kx][x2 + 1 - pkx][ky][y2 + 1 - pky]);
  return max(r1, r2);
}

/*
reset the jth bit of S = S &= ~(1 << j)
toggle the jth bit S = S ^= (1 << j)
the value of the least significant bit that is on = S & (-S)
To turn on all bits in a set of size n,
use S = (1 << n) - 1 S is a power of 2 then S & (S - 1) will be 0. 
Turn off the last bit in S = S & (S − 1)
Turn on the last zero in S = S | (S + 1)
Turn off the last consecutive run of ones in S = S & (S + 1)                                                                                              
Turn on the last consecutive run of zeroes in S = S | (S − 1) //bit 1 based indexing
*/