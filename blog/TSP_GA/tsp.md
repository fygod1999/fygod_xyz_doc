## 遗传算法
### 遗传算法简介
遗传算法(`Genetic Algorithms`，简称 `GA`)是一种基于自然选择原理和自然遗传机制的搜索(寻优)算法，它是模拟自然界中的生命进化机制，在人工系统中实现特定目标的优化。遗传算法的实质是通过群体搜索技术，根据适者生存的原则逐代进化，最终得到最优解或准最优解。它必须做以下操作:初始群体的产生、求每一个体的适应度、根据适者生存的原则选择优良个体、被选出的优良个体两两配对，通过随机交叉其染色体的基因并随机变异某些染色体的基因后生成下一代群体，按此方法使群体逐代进化，直到满足进化终止条件。其实现方法如下：
- 根据具体问题确定可行解域，确定一种编码方法，能用数值串或字符串表示 可行解域的每一解。
- 对每一解应有一个度量好坏的依据，它用一函数表示，叫做适应度函数。
- 确定进化参数群体规模M，交叉概率pc、变异概率pm、进化终止条件。

### 遗传算法求解TSP问题
问题：给定平面上20个点的名称与坐标，两个点之间的距离为它们的欧几里得距离。求一条路径，刚好经过每个点1次，使其路径长度最短。
参数设定如下：
```
种群大小：M=50
最大代数：G=1000
交叉率：pc=1,交叉率为1能保证种群的充分进化
变异率：pm=0.1,一般而言，变异发生的可能性较小
```

在该问题中，每一条路径就是所谓的染色体（解的编码），每条路径的长度就是该个体的适应性（路径长度越短，适应性越强）。交叉操作就是选择两条路径，取一个分界点`k`，将两条路径分别以分界点`k`分成前后两段，并且将两条路径重新组合得到新的两条路径。这里的交叉操作蕴含了变异操作，但是能够让子代继承父代的优良特性。变异操作也是实现群体多样性的一种手段，也是全局寻优的保证，具体实现为，按照给定的变异率，对选定的变异的个体，随机的选取三个整数`u`。

`in.txt`

```txt
34
北京 116.46 39.92
天津 117.2 39.13
上海 121.48 31.22
重庆 106.54 29.59
拉萨 91.11 29.97
乌鲁木齐 87.68 43.77
银川 106.27 38.47
呼和浩特 111.65 40.82
南宁 108.33 22.84
哈尔滨 126.63 45.75
长春 125.35 43.88
沈阳 123.38 41.8
石家庄 114.48 38.03
太原 112.53 37.87
西宁 101.74 36.56
济南 117 36.65
郑州 113.6 34.76
南京 118.78 32.04
合肥 117.27 31.86
杭州 120.19 30.26
福州 119.3 26.08
南昌 115.89 28.68
长沙 113 28.21
武汉 114.31 30.52
广州 113.23 23.16
台北 121.5 25.05
海口 110.35 20.02
兰州 103.73 36.03
西安 108.95 34.27
成都 104.06 30.67
贵阳 106.71 26.57
昆明 102.73 25.04
香港 114.1 22.2
澳门 113.33 22.13
```

```cpp
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <deque>
#include <iostream>
#include <map>
#include <queue>
#include <set>
#include <stack>
#include <string>
#include <vector>
#include <ctime>
using namespace std;
typedef long long LL;
const int maxn = 1e2 + 7;
const int INF = 0x7fffffff;
const double PI = acos(-1);
struct Point
{ //点类
    string name;
    double x, y;
    int i; //编号
};
vector<Point> p;
double d[maxn][maxn]; //距离矩阵
int n;
double sum = 0; //当前最短路径长度

double dist(Point a, Point b)
{ //计算两点距离
    return sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y));
}

double get_sum(vector<Point> a)
{ //返回路径长度
    double sum = 0;
    for (int i = 1; i < a.size(); i++)
    {
        sum += d[a[i].i][a[i - 1].i];
    }
    sum += d[a[0].i][a[a.size() - 1].i];
    return sum;
}

void init()
{                                //初始化
    srand((unsigned)time(NULL)); //设置随机数种子
    cin >> n;
    p.clear();
    for (int i = 0; i < n; i++)
    {
        Point t;
        cin >> t.name >> t.x >> t.y;
        t.i = i;
        p.push_back(t);
    }
    for (int i = 0; i < n; i++)
    {
        for (int j = i + 1; j < n; j++)
        {
            d[i][j] = d[j][i] = dist(p[i], p[j]);
        }
    }
    sum = get_sum(p);
}

void show()
{ //显示当前结果
    cout << "路径长度: " << sum << endl;
    cout << "路径:";
    for (int i = 0; i < n; i++)
        cout << ' ' << p[i].name;
    puts("");
}

int w = 100;                 //限定种群只能活100个个体
vector<vector<Point>> group; //种群，也就是染色体列表

void Improve_Circle()
{ //改良圈法得到初始序列
    vector<Point> cur = p;
    for (int t = 0; t < w; t++)
    { //重复50次
        for (int i = 0; i < n; i++)
        { //构造随机顺序
            int j = rand() % n;
            swap(cur[i], cur[j]);
        }
        int flag = 1;
        while (flag)
        {
            flag = 0;
            //不断选取uv子串，尝试倒置uv子串的顺序后解是否更优，如果更优则变更
            for (int u = 1; u < n - 2; u++)
            {
                for (int v = u + 1; v < n - 1; v++)
                {
                    if (d[cur[u].i][cur[v + 1].i] + d[cur[u - 1].i][cur[v].i] <
                        d[cur[u].i][cur[u - 1].i] + d[cur[v].i][cur[v + 1].i])
                    {
                        for (int k = u; k <= (u + v) / 2; k++)
                        {
                            swap(cur[k], cur[v - (k - u)]);
                            flag = 1;
                        }
                    }
                }
            }
        }
        group.push_back(cur);
        double cur_sum = get_sum(cur);
        if (cur_sum < sum)
        {
            sum = cur_sum;
            p = cur;
        }
    }
}

vector<int> get_randPerm(int n)
{ //返回一个随机序列
    vector<int> c;
    for (int i = 0; i < n; i++)
    {
        c.push_back(i);
    }
    for (int i = 0; i < n; i++)
    {
        swap(c[i], c[rand() % n]);
    }
    return c;
}

//排序时用到的比较函数
bool cmp(vector<Point> a, vector<Point> b) { return get_sum(a) < get_sum(b); }

int dai = 200; //一共进行200代的进化选择
int c[maxn];
double bylv = 0.1; //变异率

void genetic_algorithm()
{ //遗传算法
    vector<vector<Point>> A = group, B, C;
    // A：当前代的种群  B:交配产生的子代  C：变异产生的子代
    for (int t = 0; t < dai; t++)
    {
        B = A;
        vector<int> c = get_randPerm(A.size());
        for (int i = 0; i + 1 < c.size(); i += 2)
        {
            int F = rand() % n; //基因划分分界点
            int u = c[i], v = c[i + 1];
            for (int j = F; j < n;
                 j++)
            { //交换随机选的2个个体的基因后半段，也就是交配
                swap(B[u][j], B[v][j]);
            }
            //交换后可能发生冲突，需要解除冲突
            //保留F前面的部分不变，F后面的部分有冲突则交换
            int num1[1000] = {0}, num2[1000] = {0};
            for (int j = 0; j < n; j++)
            {
                num1[B[u][j].i]++;
                num2[B[v][j].i]++;
            }
            vector<Point> v1;
            vector<Point> v2;
            for (int j = 0; j < n; j++)
            {
                if (num1[B[u][j].i] == 2)
                {
                    v1.push_back(B[u][j]);
                }
            }
            for (int j = 0; j < n; j++)
            {
                if (num2[B[v][j].i] == 2)
                {
                    v2.push_back(B[v][j]);
                }
            }
            int p1 = 0, p2 = 0;
            for (int j = F; j < n; j++)
            {
                if (num1[B[u][j].i] == 2)
                {
                    B[u][j] = v2[p2++];
                }
                if (num2[B[v][j].i] == 2)
                {
                    B[v][j] = v1[p1++];
                }
            }
        }
        C.clear();
        int flag = 1;
        for (int i = 0; i < A.size(); i++)
        {
            if (rand() % 100 >= bylv * 100)
                continue;
            //对于变异的个体,取3个点u<v<w,把子串[u,v]插到w后面
            int u, v, w;
            u = rand() % n;
            do
            {
                v = rand() % n;
            } while (u == v);
            do
            {
                w = rand() % n;
            } while (w == u || w == v);
            if (u > v)
                swap(u, v);
            if (v > w)
                swap(v, w);
            if (u > v)
                swap(u, v);

            vector<Point> vec;
            for (int j = 0; j < u; j++)
                vec.push_back(A[i][j]);
            for (int j = v; j < w; j++)
                vec.push_back(A[i][j]);
            for (int j = u; j < v; j++)
                vec.push_back(A[i][j]);
            for (int j = w; j < n; j++)
                vec.push_back(A[i][j]);
            C.push_back(vec);
        }
        //合并A，B，C
        for (int i = 0; i < B.size(); i++)
        {
            A.push_back(B[i]);
        }
        for (int i = 0; i < C.size(); i++)
        {
            A.push_back(C[i]);
        }
        sort(A.begin(), A.end(), cmp); //从小到大排序
        vector<vector<Point>> new_A;
        for (int i = 0; i < w; i++)
        {
            new_A.push_back(A[i]);
        }
        A = new_A;
    }
    group = A;
    sum = get_sum(group[0]);
    p = group[0];
}

int main()
{
#ifndef ONLINE_JUDGE
    freopen("in.txt", "r", stdin);
#endif
    init();
    cout << "初始";
    show();
    cout << "改良圈法";
    Improve_Circle();
    show();
    cout << "遗传算法";
    genetic_algorithm();
    show();
    return 0;
}
```