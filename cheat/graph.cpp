/* ' ' has mejor effect in ordering all i(1 to n|not 0) -> graph[0].pb(i) and s[i]+=" "*/

vector<int> topologicalSort(int V)
{
    vector<int> inDegree(V, 0);

    for (int i = 0; i < V; i++)
    {
        for (int j = 0; j < graph[i].size(); j++)
            inDegree[graph[i][j]]++;
    }

    queue<int> q;
    for (int i = 0; i < V; i++)
        if (inDegree[i] == 0)
            q.push(i);

    int cnt = 0;
    vector<int> top_order;
    while (!q.empty())
    {
        int u = q.front();
        q.pop();
        top_order.push_back(u);
        for (int i = 0; i < graph[u].size(); i++)
            if (--inDegree[graph[u][i]] == 0)
                q.push(graph[u][i]);
        cnt++;
    }
    if (cnt != V)
        top_order.clear();
    return top_order;
}

int rec(int cur, int par)
{
    //ring
    //foom main
    //min(rec(u1, 1) + w1, rec(u2, 1) + w2)
    if (cur == 1)
        return 0;
    if (graph[cur].size() > 2)
        cout << "fuck" << endl;
    for (int i = 0; i < graph[cur].size(); i++)
    {
        int v = graph[cur][i].first;
        int w = graph[cur][i].second;
        if (v != par)
        {
            return w + rec(v, cur);
        }
    }
}

struct Edge
{
    int v, flow, C, rev;
};

int level[N], start[N];
vector<Edge> adj[N];

void addEdge(int u, int v, int C)
{
    Edge a{v, 0, C, adj[v].size()};
    Edge b{u, 0, C, adj[u].size()};
    adj[u].push_back(a);
    adj[v].push_back(b);
}

bool bfs(int s, int t)
{
    memset(level, -1, sizeof level);
    level[s] = 0;
    queue<int> q;
    q.push(s);
    while (!q.empty())
    {
        int u = q.front();
        q.pop();
        for (int i = 0; i < adj[u].size(); i++)
        {
            Edge &e = adj[u][i];
            if (level[e.v] < 0 && e.flow < e.C)
            {
                level[e.v] = level[u] + 1;
                q.push(e.v);
            }
        }
    }
    return level[t] < 0 ? false : true;
}

int sendFlow(int u, int flow, int t)
{
    if (u == t)
        return flow;

    for (; start[u] < adj[u].size(); start[u]++)
    {
        Edge &e = adj[u][start[u]];

        if (level[e.v] == level[u] + 1 && e.flow < e.C)
        {
            int curr_flow = min(flow, e.C - e.flow);
            int temp_flow = sendFlow(e.v, curr_flow, t);
            if (temp_flow > 0)
            {
                e.flow += temp_flow;
                adj[e.v][e.rev].flow -= temp_flow;
                return temp_flow;
            }
        }
    }

    return 0;
}

int DinicMaxflow(int s, int t)
{
    if (s == t)
        return -1;
    int total = 0;
    while (bfs(s, t))
    {
        memset(start, 0, sizeof start);
        while (int flow = sendFlow(s, INT_MAX, t))
            total += flow;
    }
    return total;
}

//O(ve2)
int isAugPath(int s, int t, int n)
{

    memset(c, -1, sizeof c);
    memset(p, -1, sizeof p);
    queue<int> q;
    q.push(s);

    while (!q.empty())
    {
        int u = q.front();
        if (u == t)
            goto path;
        q.pop();
        c[u] = 0;
        for (int i = 1; i <= n; i++)
        {
            int v = i;

            if (c[v] == -1 && r[u][v] > 0)
            {
                p[v] = u;
                c[v] = 0;
                q.push(v);
            }
        }
    }

    return 0;

path:
    int u = t, mn = INT_MAX;
    while (p[u] != -1)
    {
        mn = min(mn, r[p[u]][u]);
        u = p[u];
    }

    u = t;
    while (p[u] != -1)
    {
        r[p[u]][u] -= mn;
        r[u][p[u]] += mn;
        u = p[u];
    }
    return mn;
}

32
    /*K(m,n) is hamilton iff m=n>=2
*/

    void
    dfs(int u)
{
    d[u] = 1;
    for (int i = 0; i < graph[u].size(); i++)
    {
        int v = graph[u][i];
        if (!d[v])
            dfs(v);
    }
}

int dx[] = {-1, 0, 0, 1};
int dy[] = {0, 1, -1, 0};

int dfs(int x, int y, int n, int m)
{
    c[x][y] = 1;
    for (int i = 0; i < 4; i++)
    {
        int X = x + dx[i], Y = y + dy[i];
        if (X < n && Y < m && X >= 0 && Y >= 0 && c[X][Y] == 0 && s[X][Y] != 'D')
            dfs(X, Y, n, m);
    }
}

void bfs(int u)
{
    c[u] = 0;
    queue<int> q;
    q.push(u);
    while (!q.empty())
    {
        u = q.front();
        q.pop();
        for (int i = 0; i < graph[u].size(); i++)
        {
            int v = graph[u][i];
            if (c[v] == -1)
            {
                c[v] = c[u] + 1;
                q.push(v);
            }
        }
    }
}

void tarjan(int u)
{
    r[u] = l[u] = ticks++;
    s.push(u);
    c[u] = 1;
    for (int i = 0; i < graph[u].size(); i++)
    {
        int v = graph[u][i];
        if (r[v] == -1)
        {
            tarjan(v);
            l[u] = min(l[u], l[v]);
        }
        else if (c[v])
        {
            l[u] = min(l[u], l[v]);
        }
    }
    if (r[u] == l[u])
    {
        int v;
        cur++;
        do
        {
            v = s.top();
            s.pop();
            c[v] = 0;
            scc[v] = cur;
        } while (u != v);
    }
}

int l[N], r[N], c[N], scc, col[N], cnt;

vector<int> graph[N], ans;

void tarjan(int u)
{
    l[u] = r[u] = ++cnt;
    ans.push_back(u);
    c[u] = 1;

    for (int i = 0; i < graph[u].size(); i++)
    {
        int v = graph[u][i];
        if (r[v] == 0)
            tarjan(v);
        if (c[v])
            l[u] = min(l[u], l[v]);
    }

    if (l[u] == r[u])
    {
        ++scc;
        for (int i = 0; i < ans.size(); i++)
            col[ans[i]] = scc, c[ans[i]] = 0;
        ans.clear();
    }
}

void findCut(int u, int p, int cnt)
{
    low[u] = num[u] = cnt;
    for (int i = 0; i < graph[u].size(); i++)
    {
        int v = graph[u][i];
        if (num[v] == -1)
        {
            if (p == -1)
                rootChildren++, isCut[u] = (rootChildren > 1);
            findCut(v, u, cnt + 1);
            if (low[v] >= num[u])
                isCut[u] = true;
            if (low[v] > num[u])
                ans.push_back(ii(min(u, v), max(u, v)));
            low[u] = min(low[u], low[v]);
        }
        else if (p != v)
            low[u] = min(low[u], num[v]);
    }
}

void findcut(int u, int par)
{
    int child = 0;
    vis[u] = low[u] = ++dfsTime;
    for (int i = 0; i < graph[u].size(); i++)
    {
        int v = graph[u][i];
        if (v == par)
            continue;
        if (vis[v] != -1)
            low[u] = min(low[u], vis[v]);
        else
        {
            child++;
            findcut(v, u);
            low[u] = min(low[u], low[v]);
            if (low[v] >= vis[u])
                isCutPoint[u] = 1;
        }
    }
    if (par == -1)
        isCutPoint[u] = (child > 1);
}

void dijkstra(int s)
{
    for (int i = 0; i < M; i++)
        dist[i] = inf;
    dist[s] = 0;
    priority_queue<ii, vector<ii>, greater<ii>> pq;
    pq.push(ii(0, s));
    while (!pq.empty())
    {
        int D = pq.top().first, u = pq.top().second;
        pq.pop();
        if (dist[u] < D)
            continue;
        for (int i = 0; i < graph[u].size(); i++)
        {
            int w = graph[u][i].first, v = graph[u][i].second;
            if (dist[v] > dist[u] + w)
            {
                dist[v] = dist[u] + w;
                pq.push(ii(dist[v], v));
            }
        }
    }
}

void dfs(int x, int m)
{
    color[x] = true;
    // cout << x << endl;
    for (int i = 0; i < m; i++)
    {
        int u = edges[i].second.first, v = edges[i].second.second;
        if (x == u && color[v] == false)
        {
            dfs(v, m);
        }
    }
}

void bf(int n, int m)
{
    for (int i = 1; i < n; i++)
    {
        for (int j = 0; j < m; j++)
        {
            int u = edges[j].second.first, v = edges[j].second.second;
            int w = edges[j].first;
            d[v] = min(d[v], w + d[u]);
        }
    }
}

void rec(int n, int m)
{
    for (int j = 0; j < m; j++)
    {
        int u = edges[j].second.first, v = edges[j].second.second;
        int w = edges[j].first;
        if (d[v] > d[u] + w && color[u] == false)
            dfs(u, m);
    }
}

///is bipartite;
int bfs(int s)
{
    color[s] = 0;
    queue<int> q;
    q.push(s);
    while (!q.empty())
    {
        int u = q.front();
        q.pop();
        for (int i = 0; i < graph[u].size(); i++)
        {
            int v = graph[u][i];
            if (color[v] == -1)
            {
                color[v] = 1 + color[u];
                color[v] %= 2;
                q.push(v);
            }
            else if (color[u] == color[v])
                return 0;
        }
    }
    return 1;
}

///isDag
bool dfs(int u)
{
    color[u] = gray;
    bool ans = true;
    for (int i = 0; i < graph[u].size(); i++)
    {
        int v = graph[u][i];
        if (color[v] == white)
            ans &= dfs(v);
        else if (color[v] == gray)
            ans &= 0;
    }
    color[u] = black;
    return ans;
}

floydWarshells
        detect neg cycle /
    min cycle sum start at i end at i : set the main diagonal = inf,
                                                     graph[i][i]

    void
    init(int x, int y)
{
    for (int i = x; i <= y; i++)
    {
        for (int j = x; j <= y; j++)
        {
            if (i != j)
                d[i][j] = inf; //INT_MAX OVERFLOW
            else
                d[i][j] = 0;
            p[i][j] = i;
        }
    }
}

void assp(int x, int y)
{

    for (int k = x; k <= y; k++)
    {
        for (int i = x; i <= y; i++)
        {
            for (int j = x; j <= y; j++)
            {
                if (d[i][j] > d[i][k] + d[k][j])
                {
                    d[i][j] = d[i][k] + d[k][j];
                    p[i][j] = p[k][j];
                }
            }
        }
    }
}

void path(int x, int y)
{
    if (x != y)
        path(x, p[x][y]);
    printf(" %d", y);
}
