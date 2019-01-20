#include <bits/stdc++.h>

#define SIMPLE_TEST 0
#define RANDOM_TEST 1
#define TWO_OPT_TEST 2

#define MAX 105
#define MAX_CITY_NO 500

using namespace std;

struct Point {
    double x;
    double y;

    Point() {

    }

    Point(double a, double b) {
        x = a;
        y = b;
    }

    double deg2rad(double deg) {
        return deg * (3.14159 / 180);
    }


    double getDistanceFromLatLonInKm(double lat1, double lon1, double lat2, double lon2) {
        double R = 6371; // Radius of the earth in km
        double dLat = deg2rad(lat2 - lat1);  // deg2rad below
        double dLon = deg2rad(lon2 - lon1);
        double a = sin(dLat / 2.0) * sin(dLat / 2.0) +
                   cos(deg2rad(lat1)) * cos(deg2rad(lat2)) * sin(dLon / 2) * sin(dLon / 2);
        double c = 2 * atan2(sqrt(a), sqrt(1 - a));
        double d = R * c; // Distance in km
        return d;
    }


    double dist(Point other) {
//        return getDistanceFromLatLonInKm(x, y, other.x, other.y);
        return hypot(x - other.x, y - other.y);
    }
};

int n;
Point pt[MAX];
bool Vis[MAX];
vector<int> Tour;
int bestStartNNSimple;


class distAndPoint {
public:
    Point p;
    int d;
    int realIdx;

    distAndPoint() {

    }

    distAndPoint(Point pp, int x, int ridx) {
        p = pp;
        d = x;
        realIdx = ridx;
    }

    bool operator<(const distAndPoint &dp) {
        return (d < dp.d);
    }

};


double Cost(vector<int> Tour) {
    double Ans = 0;
    int Size = (int) Tour.size();
    for (int i = 0; i < Size; i++) {
        int xx = Tour[i];
        int yy = Tour[(i + 1) % Size];
        Ans += pt[xx].dist(pt[yy]);
    }
    return Ans;
}

void PrintSolution(vector<int> Tour) {
    printf("Solution Cost : %lf\n", Cost(Tour));
    printf("Tour Description : ");
    for (int i = 0; i < Tour.size(); i++)
        printf("%d ", Tour[i]);
    printf("\n\n");
}


class nnTourObject {
public:
    vector<int> tourList;
    int cost;

    nnTourObject() {}

    double Cost() {
        double Ans = 0;
        int Size = static_cast<int>(tourList.size());
        for (int i = 0; i < Size; i++) {
            int xx = tourList[i];
            int yy = tourList[(i + 1) % Size];
            Ans += pt[xx].dist(pt[yy]);
        }
        return Ans;
    }

    nnTourObject(vector<int> v) {
        tourList = v;
        cost = static_cast<int>(Cost());
    }

    bool operator<(const nnTourObject &nt) {
        return (cost < nt.cost);
    }
};

vector<nnTourObject> nnGreedySimple_5_Tours;
vector<nnTourObject> nnGreedyRandom_5_Tours;

vector<nnTourObject> shGreedySimple_5_Tours;
vector<nnTourObject> shGreedyRandom_5_Tours;


vector<nnTourObject> twoOptFirst_4_Tours_NN;
vector<nnTourObject> twoOptFirst_4_Tours_SH;

vector<nnTourObject> twoOptBest_4_Tours_NN;
vector<nnTourObject> twoOptBest_4_Tours_SH;


void printVector(vector<int> vec) {
    for (int i = 0; i < vec.size(); i++) {
        cout << vec[i] << " ";
    }
    cout << endl;
}

/*****************NN Simple Start**********************/

int FindNearestUnvisited(int x) {
    double Min = LLONG_MAX;
    int Idx = -1;
    for (int i = 0; i < n; i++) {
        if (Vis[i])
            continue;
        if (pt[x].dist(pt[i]) < Min) {
            Min = pt[x].dist(pt[i]);
            Idx = i;
        }
    }
    return Idx;
}


void nnGreedySimple(int startIdx) {
    Tour.clear();
    memset(Vis, 0, sizeof(Vis));

    int Start = startIdx;
    int Last = Start;
    Vis[Last] = 1;
    Tour.push_back(Start);

    while (Tour.size() < n) {
        int Idx = FindNearestUnvisited(Last);
        Vis[Idx] = 1;
        Tour.push_back(Idx);
    }
}

/*****************NN Simple End**********************/


/*****************Savings Simple Start ****************/

int adj[MAX_CITY_NO][MAX_CITY_NO] = {};
int parent[MAX_CITY_NO] = {};

void SavingAddEdge(int u, int v) {
    ++adj[u][v];
    ++adj[v][u];
}

void SavingRemoveEdge(int u, int v) {
    --adj[u][v];
    --adj[v][u];
}

bool isSavingSameRoute(int u, int v) {
    int par_u = parent[u], par_v = parent[v];
    return par_u == par_v;
}

void SavingMergeRoute(int u, int v) {
    int par_u = parent[u], par_v = parent[v];
    int p = min(par_u, par_v);

    for (int i = 0; i < n; ++i) {
        if (parent[i] == par_u || parent[i] == par_v)
            parent[i] = p;
    }
}

class savingObject {
public:
    int u, v;
    double savingsVal;

    bool operator<(const savingObject &obj) const {
        return savingsVal < obj.savingsVal;
    }


    savingObject(int x, int y, int val) {
        u = x;
        v = y;
        savingsVal = val;
    }
};

void SHGreedySimple(int randStart) {
    Tour.clear();
    memset(Vis, 0, sizeof(Vis));

    double dist[n][n];

    for (int u = 0; u < n; u++) {
        for (int v = 0; v < n; v++) {
            dist[u][v] = pt[u].dist(pt[v]);
        }
    }

    int d = randStart;


    priority_queue<savingObject> pq;

    for (int u = 0, i = 0; u < n; ++u) {
        for (int v = u + 1; v < n; ++v, ++i) {
            if (u != d && v != d) {
                savingObject obj(u, v, static_cast<int>(dist[u][d] + dist[d][v] - dist[u][v]));
                pq.push(obj);
            }
        }
    }

    memset(adj, 0, sizeof(adj)); // init to 0
    memset(parent, -1, sizeof(parent));

    for (int u = 0; u < n; ++u) {
        parent[u] = u;
    }

    for (int u = 0; u < n; ++u) {
        if (u != d) {
            SavingAddEdge(u, d);
            SavingAddEdge(u, d);
        }
    }


    while (!pq.empty()) {
        int u = pq.top().u;
        int v = pq.top().v;
        pq.pop();

        if (!isSavingSameRoute(u, v) && adj[u][d] != 0 && adj[d][v] != 0) {
            SavingAddEdge(u, v);
            SavingMergeRoute(u, v);
            SavingRemoveEdge(u, d);
            SavingRemoveEdge(v, d);
        }
    }


    int s = d;
    do {
        for (int i = 0; i < n; ++i) {
            if (adj[s][i] == 1) {
                Tour.push_back(s);
                SavingRemoveEdge(s, i); // removes the edge so that same edge doesn't come twice
                s = i;
            }
        }
    } while (s != d);

}


// 1 = NN , 2 = SH
vector<nnTourObject> SimpleUtil(int NNorSH) {
    string name = "";
    vector<nnTourObject> nn_5_Tours;

    int t = 0;
    while (t < 5) {
        t++;

        int randomStart = rand() % n;

        if (NNorSH == 1) {
            nnGreedySimple(randomStart);
            name = "Nearest Neighbor Greedy Simple Tour: ";
        } else if (NNorSH == 2) {
            SHGreedySimple(randomStart);
            name = "Savings Heuristic Greedy Simple Tour: ";
        }

        nnTourObject nnTour(Tour);
        nn_5_Tours.push_back(nnTour);
    }

    sort(nn_5_Tours.begin(), nn_5_Tours.end());

#ifdef SIMPLE_TEST
    cout << name << endl;
    cout << "Best Case: " << endl;
    cout << "   Start Node: " << nn_5_Tours[0].tourList[0] << endl;
    cout << "   Tour : ";
    printVector(nn_5_Tours[0].tourList);
    cout << "   Cost: " << nn_5_Tours[0].cost << endl;
#endif

    bestStartNNSimple = nn_5_Tours[0].tourList[0];

#ifdef SIMPLE_TEST
    cout << "Worst Case: " << endl;
    cout << "   Start Node: " << nn_5_Tours[4].tourList[0] << endl;
    cout << "   Tour : ";
    printVector(nn_5_Tours[4].tourList);
    cout << "   Cost: " << nn_5_Tours[4].cost << endl;
#endif

    double sum = 0;
    for (int i = 0; i < nn_5_Tours.size(); i++) {
        sum += nn_5_Tours[i].cost;
    }
    sum = sum / nn_5_Tours.size();

#ifdef SIMPLE_TEST
    cout << "Average Cost: " << sum << endl;
#endif
    nnGreedySimple_5_Tours = nn_5_Tours;

    return nn_5_Tours;
}

/******************* Savings Simple End *******************/



/***********NN Random Start***************/

int FindNearestRandomUnvisited(int x) {
    vector<distAndPoint> nearest_nodes;

    for (int i = 0; i < n; i++) {
        if (!Vis[i]) {
            distAndPoint tem(pt[i], pt[x].dist(pt[i]), i);
            nearest_nodes.push_back(tem);
        }
    }

    sort(nearest_nodes.begin(), nearest_nodes.end());

    int divisor = 5;
    if (nearest_nodes.size() < 5) divisor = static_cast<int>(nearest_nodes.size());
    int randomIdx = rand() % divisor;

    distAndPoint tem = nearest_nodes[randomIdx];
    int Idx = tem.realIdx;

    return Idx;
}


void nnGreedyRandomized(int startIdx) {
    Tour.clear();
    memset(Vis, 0, sizeof(Vis));

    int Start = startIdx;
    int Last = Start;
    Vis[Last] = true;
    Tour.push_back(Start);

    while (Tour.size() < n) {

        int Idx = FindNearestRandomUnvisited(Last);
        Vis[Idx] = true;
        Tour.push_back(Idx);
    }
}



/***********NN Random End***************/


/*********Savings Random Start *****************/

void shGreedyRandomized(int Start) {
    Tour.clear();
    memset(Vis, 0, sizeof(Vis));

    double dist[n][n];

    for (int u = 0; u < n; u++) {
        for (int v = 0; v < n; v++) {
            dist[u][v] = pt[u].dist(pt[v]);
        }
    }

    int d = Start;


    priority_queue<savingObject> pq;

    for (int u = 0, i = 0; u < n; ++u) {
        for (int v = u + 1; v < n; ++v, ++i) {
            if (u != d && v != d) {
                savingObject obj(u, v, static_cast<int>(dist[u][d] + dist[d][v] - dist[u][v]));
                pq.push(obj);
            }
        }
    }

    memset(adj, 0, sizeof(adj)); // init to 0
    memset(parent, -1, sizeof(parent));

    for (int u = 0; u < n; ++u) {
        parent[u] = u;
    }

    for (int u = 0; u < n; ++u) {
        if (u != d) {
            SavingAddEdge(u, d);
            SavingAddEdge(u, d);
        }
    }


    while (!pq.empty()) {


        vector<savingObject> vec;
        int highest = 5;
        if (pq.size() < 5) {
            highest = static_cast<int>(pq.size());
        }
        for (int i = 0; i < highest; i++) {
            vec.push_back(pq.top());
            pq.pop();
        }

        int randomIdx = rand() % highest;

        int u = vec[randomIdx].u;
        int v = vec[randomIdx].v;

        for (int i = 0; i < highest; i++) {
            if (i == randomIdx) continue;
            pq.push(vec[i]);
        }
        vec.clear();

        if (!isSavingSameRoute(u, v) && adj[u][d] != 0 && adj[d][v] != 0) {
            SavingAddEdge(u, v);
            SavingMergeRoute(u, v);
            SavingRemoveEdge(u, d);
            SavingRemoveEdge(v, d);
        }
    }


    int s = d;
    do {
        for (int i = 0; i < n; ++i) {
            if (adj[s][i] == 1) {
                Tour.push_back(s);
                SavingRemoveEdge(s, i); // removes the edge so that same edge doesn't come twice
                s = i;
            }
        }
    } while (s != d);

}


// NN = 1 , SH = 2
vector<nnTourObject> RandomizedUtil(int NNorSH) {
    string name = "";
    vector<nnTourObject> nn_5_Tours;
    vector<nnTourObject> simpleNNTour;

    int t = 0;

    while (t < 10) {
        t++;
        if (NNorSH == 1) {
            nnGreedyRandomized(bestStartNNSimple);
            name = "Nearest Neighbor Greedy Randomized Tour: ";
        } else if (NNorSH == 2) {
            shGreedyRandomized(bestStartNNSimple);
            name = "Savings Heuristic Greedy Randomized Tour: ";
        }
        nnTourObject nnTour(Tour);
        nn_5_Tours.push_back(nnTour);
    }

    sort(nn_5_Tours.begin(), nn_5_Tours.end());

#ifdef RANDOM_TEST
    cout << name << endl;
    cout << "Best Case: " << endl;
    cout << "   Start Node: " << nn_5_Tours[0].tourList[0] << endl;
    cout << "   Tour : ";
    printVector(nn_5_Tours[0].tourList);
    cout << "   Cost: " << nn_5_Tours[0].cost << endl;

    cout << "Worst Case: " << endl;
    cout << "   Start Node: " << nn_5_Tours[4].tourList[0] << endl;
    cout << "   Tour : ";
    printVector(nn_5_Tours[4].tourList);
    cout << "   Cost: " << nn_5_Tours[4].cost << endl;
#endif

    double sum = 0;
    for (int i = 0; i < nn_5_Tours.size(); i++) {
        sum += nn_5_Tours[i].cost;
    }
    sum = sum / nn_5_Tours.size();

#ifdef RANDOM_TEST
    cout << "Average Cost: " << sum << endl;
#endif

    return nn_5_Tours;
}


/********** Savings Random End ****************/






/************2-Opt First improve Start***************/
void TwoOptHeuristicFirstImprovement(vector<int> vec) {
    Tour.clear();
    memset(Vis, 0, sizeof(Vis));

    Tour = vec;

    while (true) {
        double Curr = Cost(Tour);
        bool Changed = false;

        for (int i = 0; i < Tour.size(); i++) {
            for (int j = i + 2; j < Tour.size(); j++) {
                reverse(Tour.begin() + i + 1, Tour.begin() + j + 1);
                double NewCost = Cost(Tour);
                if (NewCost < Curr) {
                    Changed = true;
                    break;
                }
                reverse(Tour.begin() + i + 1, Tour.begin() + j + 1);
            }
            if (Changed)
                break;
        }
        if (!Changed)
            break;
    }
}

/***************2-Opt best improve start***********/
void TwoOptHeuristicBestImprovement(vector<int> vec) {
    int ii = 0, jj = 0;
    Tour.clear();
    memset(Vis, 0, sizeof(Vis));

    Tour = vec;

    while (true) {
        double Curr = Cost(Tour);
        bool Changed = false;

        for (int i = 0; i < Tour.size(); i++) {
            for (int j = i + 2; j < Tour.size(); j++) {
                reverse(Tour.begin() + i + 1, Tour.begin() + j + 1);
                double NewCost = Cost(Tour);
                if (NewCost < Curr) {
                    Changed = true;
                    ii = i;
                    jj = j;
                    Curr = NewCost;
                }
                reverse(Tour.begin() + i + 1, Tour.begin() + j + 1);
            }
            if (Changed) {
                reverse(Tour.begin() + ii + 1, Tour.begin() + jj + 1);
                break;
            }
        }

        if (!Changed)
            break;
    }
}


// twoOptType == 1 -> First , 2 -> Best
// NNorSH == 1 -> NN , 2 -> SH

vector<nnTourObject> twoOptUtil(int NNorSH, int twoOptType) {
    string name = "";
    vector<nnTourObject> myTourResults;
    vector<nnTourObject> bestToursOfNNHandSH;

    if (NNorSH == 1) {
        bestToursOfNNHandSH.push_back(nnGreedyRandom_5_Tours[0]);
        bestToursOfNNHandSH.push_back(nnGreedyRandom_5_Tours[1]);
        bestToursOfNNHandSH.push_back(nnGreedyRandom_5_Tours[2]);
        bestToursOfNNHandSH.push_back(nnGreedySimple_5_Tours[0]);

        if (twoOptType == 1) {
            for (int i = 0; i < bestToursOfNNHandSH.size(); i++) {
                TwoOptHeuristicFirstImprovement(bestToursOfNNHandSH[i].tourList);
                nnTourObject tem(Tour);
                myTourResults.push_back(tem);
            }
            name = "2-Opt First for NNH Greedy Random.";
        } else if (twoOptType == 2) {
            for (int i = 0; i < bestToursOfNNHandSH.size(); i++) {
                TwoOptHeuristicBestImprovement(bestToursOfNNHandSH[i].tourList);
                nnTourObject tem(Tour);
                myTourResults.push_back(tem);
            }
            name = "2-Opt Best for NNH Greedy Random.";
        }

    } else if (NNorSH == 2) {
        //SH vectors go here.
        bestToursOfNNHandSH.push_back(shGreedyRandom_5_Tours[0]);
        bestToursOfNNHandSH.push_back(shGreedyRandom_5_Tours[1]);
        bestToursOfNNHandSH.push_back(shGreedyRandom_5_Tours[2]);
        bestToursOfNNHandSH.push_back(shGreedySimple_5_Tours[0]);

        if (twoOptType == 1) {
            for (int i = 0; i < bestToursOfNNHandSH.size(); i++) {
                TwoOptHeuristicFirstImprovement(bestToursOfNNHandSH[i].tourList);
                nnTourObject tem(Tour);
                myTourResults.push_back(tem);
            }
            name = "2-Opt First for SH Greedy Random.";
        } else if (twoOptType == 2) {
            for (int i = 0; i < bestToursOfNNHandSH.size(); i++) {
                TwoOptHeuristicBestImprovement(bestToursOfNNHandSH[i].tourList);
                nnTourObject tem(Tour);
                myTourResults.push_back(tem);
            }
            name = "2-Opt First for SH Greedy Random.";
        }
    }


    sort(myTourResults.begin(), myTourResults.end());

#ifdef TWO_OPT_TEST
    cout << name << " : " << endl;
    cout << "Best Case: " << endl;
    cout << "   Start Node: " << myTourResults[0].tourList[0] << endl;
    cout << "   Tour : ";
    printVector(myTourResults[0].tourList);
    cout << "   Cost: " << myTourResults[0].cost << endl;

    cout << "Worst Case: " << endl;
    cout << "   Start Node: " << myTourResults[myTourResults.size() - 1].tourList[0] << endl;
    cout << "   Tour : ";
    printVector(myTourResults[myTourResults.size() - 1].tourList);
    cout << "   Cost: " << myTourResults[myTourResults.size() - 1].cost << endl;
#endif

    double sum = 0;
    for (int i = 0; i < myTourResults.size(); i++) {
        sum += myTourResults[i].cost;
    }
    sum = sum / myTourResults.size();

#ifdef TWO_OPT_TEST
    cout << "Average Cost: " << sum << endl;
#endif

    return myTourResults;
}


int main() {
//    freopen("pr76.tsp", "r", stdin);
//    freopen("pr76_out.txt", "w", stdout);
//    cout << "Output For pr76.tsp Dataset." << endl << endl;


//    freopen("berlin52.tsp", "r", stdin);
//    freopen("berlin52_out.txt", "w", stdout);
//    cout << "Output For berlin52.tsp Dataset." << endl << endl ;


    freopen("st70.tsp", "r", stdin);
    freopen("st70_out.txt", "w", stdout);
    cout << "Output For st70.tsp Dataset." << endl << endl ;

    srand(time(NULL));

    scanf("%d", &n);
    for (int i = 0; i < n; i++) { scanf("%lf %lf", &pt[i].x, &pt[i].y); }

    /***************Run NN_Greedy_Simple()*****************/

    cout << endl << endl << "##############################################################################" << endl
         << endl;

    nnGreedySimple_5_Tours = SimpleUtil(1);
    cout << endl << endl << "##############################################################################" << endl
         << endl;

    shGreedySimple_5_Tours = SimpleUtil(2);
    cout << endl << endl << "##############################################################################" << endl
         << endl;

    nnGreedyRandom_5_Tours = RandomizedUtil(1);
    cout << endl << endl << "##############################################################################" << endl
         << endl;

    shGreedyRandom_5_Tours = RandomizedUtil(2);
    cout << endl << endl << "##############################################################################" << endl
         << endl;

    /* (nnsh, twoOptType)..
     * NNorSH == 1 -> NN , 2 -> SH
     * twoOptType == 1 -> First , 2 -> Best
     */
    twoOptFirst_4_Tours_NN = twoOptUtil(1, 1);
    cout << endl << endl << "##############################################################################" << endl
         << endl;

    twoOptFirst_4_Tours_SH = twoOptUtil(2, 1);
    cout << endl << endl << "##############################################################################" << endl
         << endl;

    twoOptBest_4_Tours_NN = twoOptUtil(1, 2);
    cout << endl << endl << "##############################################################################" << endl
         << endl;

    twoOptBest_4_Tours_SH = twoOptUtil(2, 2);
    cout << endl << endl << "##############################################################################" << endl
         << endl;

    return 0;
}
