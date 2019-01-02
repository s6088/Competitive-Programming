/*
d(0) = 1, d(1) = 0
d(n) = (n-1) * ( d(n-1) + d(n-2) )
shakib -> tamim,  shakib <-> tamim
*/
long long dea(ll n)
{
    if (n == 0)
        return 1;
    if (n == 1)
        return 0;
    if (dp[n] != -1)
        return dp[n];
    return dp[n] = ((n - 1) * (dea(n - 1) + dea(n - 2)));
}

long long comb(int n, int k)
{
    if (k == 0 || n == k)
        return 1;
    if (dp[n][k] == -1)
        dp[n][k] = (comb(n - 1, k) + comb(n - 1, k - 1)) % mod;
    return dp[n][k];
}
// finds how many permutations of a size-n set does not have any of the firts k elements
// in their correspondient position (kth position).
//Number of permutations of n elements such that exactly k of them are in their correct position.
//f(n, k) = comb (n, k) * D(n - k)
long long dea(int n, int k)
{
    long long ans = fact[n];
    for (int i = 1; i <= k; ++i)
    {
        ans -= ((i & 1) ? 1 : -1) * (comb(k, i) * fact[n - i]) % mod;
        ans = (ans + mod) % mod;
    }
    return ans;
}

/*the number of lattice points between ab = 1 + __gcd(abs(ax - bx), abs(ay - by))

If d is a divisor of n, then there are Phi(n/d) numbers i <= n for which gcd(i,n)=d

a number n have log10(n) + 1 digit
number of digit in fact(n) in base b = floor(cumlog[1 to n]/log(b)) + 1

sum i to n (sum of real divisor i) there exist a pattern

a ^ n + b ^ n = (a + b)(a ^ n - 1 + b ^ n - 1) - ab(a ^ n - 2 + b ^ n - 2)
mod of negative number(mod + (a % mod)) % mod
                
Derangement : 1, 0, 1, 2, 9, 44, 265, 1854, 14833, 133496, ... 
Catalan numbers : 1, 1, 2, 5, 14, 42, 132, 429, 1430, 4862, ... 
pow(2, p - 1) * (pow(2, p) - 1) is a perfect number whenever(pow(2, p) - 1) is prime there can be only a prime factor of x greater than sqrt(x)

cat(n) = 2n ! / (n !*(n + 1) !) cat(n + 1) = (2n + 2) * (2n + 1) * cat(n) / (n + 2) * (n + 1)

                */

//string mod
int strmod(string num, int a)
{
    int res = 0;
    for (int i = 0; i < num.length(); i++)
        res = (res * 10 + (int)num[i] - '0') % a;
    return res;
}

//predicate binary search
int bs(int l, int r, int f)
{
    if (l > r)
        return -1;
    int m = (l + r) / 2;
    int x = prfSum(m), y = prfSum(m - 1);
    if (x == f && y < f)
        return m;
    if (x >= f)
        r = m - 1;
    if (x < f)
        l = m + 1;
    return bs(l, r, f);
}

void solve()
{
    for (i = 2; i * i <= n; ++i)
    {
        j = n / i;
        ans += (i + j) * (j - i + 1) / 2;
        ans += i * (j - i);
    }
}

inline int trialingZerosInFactN(int n)
{
    if (n)
        return n / 5 + zero(n / 5);
    return 0;
}
inline int firstThreeDigitN_K(int n, int k)
{
    double x = k * log10(n);
    return int(pow(10, x - LL(x)) * 100);
}
template <class T>
inline T mpow(T p, T e, T M)
{
    if (e == 0)
        return 1;
    if (e % 2 == 0)
    {
        ll t = mpow(p, e / 2, M);
        return (T)((t * t) % M);
    }
    return (T)((ll)mpow(p, e - 1, M) * (ll)p) % M;
}

template <class T>
inline T modinv(T a, T M) { return mpow(a, M - 2, M); }

///Derangement

!n = (n - 1)(!(n - 1) + !(n - 2)) !n = n ! - nC1(n - 1) ! + nC2(n - 2) ! - ....nCn * 0 !

                                                                               struct matrix
{
    ll m[10][10];
};

matrix Mul(matrix aa, matrix bb, int n)
{
    matrix ret;
    memset(ret.m, 0, sizeof ret.m);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
        {
            for (int k = 0; k < n; k++)
            {
                ret.m[i][j] += (aa.m[i][k] * bb.m[k][j]);
                ret.m[i][j] %= mod;
            }
        }
    return ret;
}

matrix exp(matrix a, ll b, int n)
{
    if (b == 1)
        return a;
    if (b & 1)
    {
        return Mul(a, exp(a, b - 1, n), n);
    }
    else
    {
        matrix tmp = exp(a, b / 2, n);
        return Mul(tmp, tmp, n);
    }
}

Sn = n / 2 * (2a + (n - 1) * d) Sn = a * (r ^ n - 1) / (r - 1) void extEuclid(int a, int b, int &x, int &y, int &gcd)
{
    x = 0;
    y = 1;
    gcd = b;
    int m, n, q, r;
    for (int u = 1, v = 0; a != 0; gcd = a, a = r)
    {
        q = gcd / a;
        r = gcd % a;
        m = x - u * q;
        n = y - v * q;
        x = u;
        y = v;
        u = m;
        v = n;
    }
}

void gcd(int a, int b, int &x, int &y, int &g)
{
    if (b == 0)
    {
        x = 1, y = 0, g = a;
        return;
    }
    gcd(b, a % b, x, y, g);
    int x1 = y, y1 = x - (a / b) * y;
    x = x1, y = y1;
}

void gcd(int a, int b)
{
    if (b == 0)
    {
        x = 1, y = 0, g = a;
        return;
    }
    gcd(b, a % b);
    int x1 = y, y1 = x - (a / b) * y;
    x = x1, y = y1;
}

long long modInv(long long n, long long m)
{
    int x, y, gcd;
    extEuclid(n, m, x, y, gcd);
    if (gcd == 1)
        return x % m;
    return 0;
}
void gcd(int a, int b)
{
    if (b == 0)
    {
        x = 1, y = 0, g = a;
        return;
    }
    gcd(b, a % b);
    int x1 = y, y1 = x - (a / b) * y;
    x = x1, y = y1;
}

long long nCr(long long n, long long r)
{
    if (dp[n][r] != -1)
        return dp[n][r];
    if (n == r)
        return 1;
    if (r == 1)
        return n;
    return dp[n][r] = nCr(n - 1, r - 1) + nCr(n - 1, r);
}

bool isprime(long long x)
{
    for (long long i = 2; i * i <= x; i++)
        if (x % i == 0)
            return false;
    return true;
}
//segmented sieve
void simpleSieve(ll limit)
{
    bool mark[limit + 1];
    memset(mark, true, sizeof(mark));

    for (ll p = 2; p * p <= limit; p++)
        if (mark[p] == true)
            for (ll i = p * 2; i <= limit; i += p)
                mark[i] = false;
    for (ll p = 2; p <= limit; p++)
        if (mark[p] == true)
            prime.push_back(p);
}

ll segmentedSieve(ll low, ll high)
{
    low = max(2ll, low);
    ll limit = high - low + 1, cnt = 0;
    bool mark[limit + 1];
    memset(mark, true, sizeof(mark));

    for (ll i = 0; i < prime.size() && prime[i] * prime[i] <= high; i++)
    {
        ll loLim = max((low / prime[i]), 2ll) * prime[i];
        if (loLim < low)
            loLim += prime[i];
        for (ll j = loLim; j <= high; j += prime[i])
            mark[j - low] = false;
    }

    for (ll i = low; i <= high; i++)
        if (mark[i - low] == true)
            cnt++;

    return cnt;
}

int minSwaps(int n)
{

    int ans = 0;

    for (int i = 0; i < n; i++)
    {
        if (is[i] || a[i].S == i)
            continue;
        int cS = 0;
        for (int j = i; !is[j]; j = a[j].S)
        {
            is[j] = 1, cS++;
        }
        ans += (cS - 1);
    }

    return ans;
}

int kadane(int l, int r)
{
    int msf = INT_MIN, meh = 0;
    for (int i = l; i <= r; i++)
    {
        meh = meh + a[i];
        if (msf < meh)
            msf = meh;
        if (meh < 0)
            meh = 0;
    }
    return msf;
}

//merging
void merge(int l, int r)
{
    int m = (l + r) / 2;
    for (int i = l, j = m + 1, k = 0; i <= m || j <= r; k++)
    {
        int x = (i <= m) ? a[i] : INT_MAX, y = (j <= r) ? a[j] : INT_MAX;
        if (x < y)
            ta[k] = a[i++];
        else
            ta[k] = a[j++];
    }
    for (int i = l, k = 0; i <= r; i++, k++)
    {
        a[i] = ta[k];
    }
}

void msort(int l, int r)
{
    if (l >= r)
        return;
    int m = (l + r) / 2;
    msort(l, m);
    msort(m + 1, r);
    merge(l, r);
}

//general sieve

void sieve()
{
    primes.clear();
    for (int i = 4; i < N; i += 2)
        isPrime[i] = 1;
    isPrime[0] = isPrime[1] = 1;
    int sq = sqrt(N);
    for (int i = 3; i <= sq; i += 2)
    {
        if (!isPrime[i])
        {
            for (int j = i + i; j < N; j += i)
                isPrime[j] = 1;
        }
    }
    primes.push_back(2);
    for (int i = 3; i < N; i += 2)
        if (!isPrime[i])
            primes.push_back(i);
}

///
int cntDiv(int x)
{
    if (nd[x] || x == 0)
        return nd[x];
    int i = x;
    map<int, int> cnt;
    while (i > 1)
        cnt[sd[i]]++, i /= sd[i];
    int ans = 1;
    for (map<int, int>::iterator it = cnt.begin(); it != cnt.end(); it++)
        ans *= (it->second) + 1;
    return nd[x] = ans;
}
///main
for (int i = 0; i < N; i++)
    sd[i] = i;
for (int i = 2; i * i < N; i++)
{
    if (sd[i] != i)
        continue;
    for (int j = i * i; j < N; j += i)
        if (sd[j] == j)
            sd[j] = i;
}
///

/*Euler's totient function  applied to a positive integer n is defined to be
the number of positive integers less than or equal to n that are relatively prime to n*/
note 1 : dont cont gcd(1, 1) = 1 if your graph has no self loop or similar, nirob vi

    void phi()
{
    for (int i = 0; i < M; i++)
        relPrimeLE[i] = i;
    for (int i = 2; i < M; i++)
    {
        if (relPrimeLE[i] != i)
            continue;
        for (int j = i; j < M; j += i)
            relPrimeLE[j] -= (relPrimeLE[j] / i);
    }
}
