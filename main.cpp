#include <iostream>
#include<vector>
#include<math.h>
#include<cstring>
#include<algorithm>
#include <numeric>
#include <ctime>
#include <direct.h>
#include <set>
#include<fstream>

using namespace std;
#pragma GCC optimize(3)

const int n = 29;
double bestRoute = 27603; 
const char *dir = "Sahara.txt";

//const int n = 38;
//double bestRoute = 6659.432;
//const char* dir = "Djibouti.txt";

//const int n = 194;
//double bestRoute = 9353.55;
//const char *dir = "Qatar.txt";

bool PRE_EA = false; // Whether initialize population with pure Genetic Algorithn
double TIME_LIMIT = 150; 
int hybrid_generation = 100;
int popuSize = 100;
double m_rate = 0.05;
int experiment_times = 1;
int robin_constant = popuSize / 10;

double x[n], y[n];
double G[n][n];

double eval(vector<int> &a) {
    double ans = 0;
    for (int i = 0; i < n; i++) {
        if (i != n - 1) {
            ans += G[a[i]][a[i + 1]];
        } else ans += G[a[i]][a[0]];
    }
    return 1 / ans;
}

vector<double> eval_vector(vector<vector<int>> &p) {
    vector<double> temp;
    for (int i = 0; i < p.size(); i++) {
        double fit = eval(p[i]);
        temp.push_back(fit);
    }
    return temp;
}

double dist(int i, int j) {

    return sqrt(1.0 * (x[i] - x[j]) * (x[i] - x[j]) + (y[i] - y[j]) * (y[i] - y[j]));
}

vector<int> ans;
double bestv = 1.0 / 2047483648.0;

double getBest(vector<vector<int>> p) {
    vector<double> score = eval_vector(p);
    vector<int> t;
    double temp = 0;
    for (int i = 0; i < score.size(); ++i) {
        if (score[i] > temp) {
            temp = score[i];
            t = p[i];
        }
    }
    if (bestv < temp) {
        bestv = temp;
        ans = t;
    }
    return temp;
}

vector<double> Stat(vector<double> score) {
    double mean = 0;
    for (auto num:score)mean += 1.0 / num;
    mean /= (score.size() * 1.0);
    double var = 0;
    for (auto num:score)var += (1.0 / num - mean) * (1.0 / num - mean);
    var /= (score.size() * 1.0);
    double b = 0;
    for (auto num:score)b = max(num, b);
    printf("In population：best = %f mean = %f var = %f\n", 1.0 / b, mean, var);
    return {mean, var};
}


bool vis[n];

vector<int> order_crossover(vector<int> &a, vector<int> b) {
//    b = calibra(a,b);
    memset(vis, false, sizeof(vis));
    int pos1 = rand() % n;
    int pos2 = pos1 + rand() % (n - pos1);
    vector<int> temp;
    for (int i = pos1; i <= pos2; i++) {
        vis[a[i]] = true;
        temp.push_back(a[i]);
    }
    for (int i = 0; i < n; i++) {
        if (!vis[b[i]])temp.push_back(b[i]);
    }
    return temp;
}

vector<int> mutation(vector<int> a) {
    int pos1 = rand() % n;
    int pos2 = rand() % n;
    while (pos2 == pos1)pos2 = rand() % n;
    swap(a[pos1], a[pos2]);
    return a;
}

vector<int> shuffle_mutation(vector<int> a, int k) {
    int start = rand() % k;
    random_shuffle(a.begin() + start, a.begin() + start + k);
    return a;
}

vector<int> edge_mutation(vector<int> a) {
    int pos1 = rand() % n;
    int pos2 = rand() % n;
    if (pos1 > pos2)swap(pos1, pos2);
    reverse(a.begin() + pos1, a.begin() + pos2);
    return a;
}


int wheel(vector<double> prob) {
    double r = rand() % 10000 * 1.0 / 10000.0;
    int i = 0;
    while (r > 0 and i < prob.size() - 1) {
        r -= prob[i];
        i += 1;
    }
    return i;
}

vector<double> robin(vector<double> score) {
    vector<double> cnt(score.size(), 0);
    vector<int> pos;
    for (int i = 0; i < score.size(); ++i)pos.push_back(i);
    for (int i = 0; i < score.size(); i++) {
        random_shuffle(pos.begin(), pos.end());
        for (int j = 0; j < robin_constant; ++j) {
            if (score[i] >= score[pos[j]])cnt[i]++;
        }
    }
    return cnt;

}


vector<double> rank_score(vector<int> pos) {
    vector<double> ans(pos.size(), 0.0);
    double tot = pos.size() * (pos.size() - 1) / 2;
//    printf("%f %d\n",tot,pos.size());
    for (int i = 0; i < pos.size(); i++) {
//        printf("%d %d %d\n",pos[i],n,i);
        ans[pos[i]] = (pos.size() - i) * 1.0 / tot;
    }
    return ans;
}

bool inplace_MN_Search(vector<int> &a) {
    double ori = 1.0 / eval(a);
    int swapi = 0, swapj = 0;
    double swap_minv = ori;
    for (int i = 1; i < n; i++) {
        for (int j = i + 3; j < n - 1; j++) {
            double temp = ori - G[a[i - 1]][a[i]] - G[a[i]][a[i + 1]] - G[a[j - 1]][a[j]] - G[a[j]][a[j + 1]] +
                          G[a[i - 1]][a[j]] + G[a[j]][a[i + 1]] + G[a[j - 1]][a[i]] + G[a[i]][a[j + 1]];
            if (temp < swap_minv) {
                swap_minv = temp;
                swapi = i;
                swapj = j;
            }
        }
    }
    double RL_minv = ori;
    int RL_i = 0, RL_j = 0;
    for (int i = 1; i < n; i++) {
        for (int j = i + 3; j < n - 1; j++) {
//            printf("%d %d\n",i,j);
            int mid = (i + j) / 2;
            double temp =
                    ori - G[a[mid]][a[mid + 1]] - G[a[i - 1]][a[i]] - G[a[j]][a[j + 1]] + G[a[i - 1]][a[mid + 1]] +
                    G[a[j]][a[i]] + G[a[mid]][a[j + 1]];
            if (temp < RL_minv) {
                RL_minv = temp;
                RL_i = i, RL_j = j;
            }
        }
    }
    double inverse_minv = ori;
    int inverse_i = 0, inverse_j = 0;
    for (int i = 1; i < n; i++) {
        for (int j = i + 3; j < n - 1; j++) {
            double temp = ori - G[a[i - 1]][a[i]] - G[a[j]][a[j + 1]] + G[a[i - 1]][a[j]] + G[a[i]][a[j + 1]];

            if (temp < inverse_minv) {
                inverse_minv = temp;
                inverse_i = i;
                inverse_j = j;
            }
        }
    }

    double swap2_min = ori;
    int swap2_i = 0;
    for (int i = 0; i < n - 1; i++) {
        double temp = ori - G[a[(i - 1 + n) % n]][a[i]] - G[a[(i + 1) % n]][a[(i + 2) % n]] +
                      G[a[(i - 1 + n) % n]][a[(i + 1) % n]] + G[a[i]][a[(i + 2) % n]];
        if (temp < swap2_min) {
            swap2_min = temp;
            swap2_i = i;
        }
    }

//    if(rand()%n<10)
//    printf("%.1f %.1f %.1f %.1f %.1f\n", swap2_min, inverse_minv, RL_minv, swap_minv, ori);

    double minv = min(swap2_min, min(inverse_minv, min(RL_minv, swap_minv)));
    if (minv == ori) {
//        shuffle_mutation(a,rand()%10 + 10);
        return false;
    }
    if (inverse_minv == minv) {
//        printf("inverse: %d %d\n",inverse_i,inverse_j);
        reverse(a.begin() + inverse_i, a.begin() + inverse_j + 1);
        inplace_MN_Search(a);
    }
    if (swap_minv == minv) {
//        printf("swap: %d %d %.1f\n",swapi,swapj,1.0/eval(a));
        swap(a[swapi], a[swapj]);
//        printf("after swap: %.1f\n",1.0/eval(a));
        inplace_MN_Search(a);
    }
    if (RL_minv == minv) {
        //        printf("ROTATE: %d %d\n",RL_i,RL_j);
        vector<int> temp;
        for (int i = 0; i < RL_i; i++)temp.push_back(a[i]);
        for (int i = (RL_i + RL_j) / 2 + 1; i <= RL_j; i++) {
            temp.push_back(a[i]);
        }
        for (int i = RL_i; i <= (RL_i + RL_j) / 2; i++)
            temp.push_back(a[i]);
        for (int i = RL_j + 1; i < n; i++)
            temp.push_back(a[i]);
        a = temp;
        inplace_MN_Search(a);
    }
    if (swap2_min == minv) {
//        printf("swap2: %d\n",swap2_i);
        swap(a[swap2_i], a[(swap2_i + 1) % n]);
        inplace_MN_Search(a);
    }
    return true;
}


vector<int> diverse(vector<int> a, vector<vector<int>> P) {
    vector<int> cnt(P.size(), 0);
    for (int i = 0; i < P.size(); i++) {
        int s = 0;
        while(P[i][s]!=a[0])s++;
        for (int j = 0; j < n; j++) {
            if (a[j] != P[i][(s+j)%n])cnt[i]++;
        }
    }
    vector<int> pos;
    for (int i = 0; i < P.size(); i++)pos.push_back(i);
    sort(pos.begin(), pos.end(), [&](int i, int j) {
        return cnt[i] > cnt[j];
    });
    return pos;
}


vector<vector<int>> hybrid_nextGeneration(vector<vector<int>> P, double mutate_rate) {

    vector<vector<int>> ans;
    for (int i = 0; i < P.size(); i++)
        if (!inplace_MN_Search(P[i]) && ans.size() < popuSize) {
            //            P[i] = shuffle_mutation(P[i],n/10);
            ans.push_back(shuffle_mutation(P[i], n / 10));

        }
    getBest(P);

    vector<double> score = eval_vector(P);
    double best = 0;
    for (auto v:score)best = max(best, v);
    if (best >= 1.0 / bestRoute)return P;
    int n_p = popuSize;

    vector<int> pos;
    for (int i = 0; i < n_p; i++)pos.push_back(i);

    sort(pos.begin(), pos.end(), [&](int i, int j) {
        return score[i] > score[j];
    });
//    classic elitism
    for (int i = 0; i < n_p / 10; i++)ans.push_back(P[pos[i]]);
//      diverse elitism
//    auto d_pos = diverse(P[pos[0]], P);
//    for (int i = 0; i < n_p / 10; i++)ans.push_back(P[d_pos[i]]);

    score = rank_score(pos);

    int cnt = ans.size();

    while (cnt < 2 * n_p) {

        int i = wheel(score);
        int j = wheel(score);
        while (j == i)j = wheel(score);

        ans.push_back(order_crossover(P[i], P[j]));
        ans.push_back(order_crossover(P[j], P[i]));
        cnt += 2;

        if (rand() % 100 * 1.0 / 100.0 < mutate_rate) {
            ans.push_back(mutation(P[i]));
            ans.push_back(mutation(P[j]));
            ans.push_back(edge_mutation(P[i]));
            ans.push_back(edge_mutation(P[j]));
            cnt += 4;
        }
    }

    score = eval_vector(ans);

    for (auto v:score)best = max(best, v);
    if (best >= 1.0 / bestRoute)return P;

    score = robin(score);
    for (int i = popuSize; i != ans.size(); i++)
        pos.push_back(i);
    sort(pos.begin(), pos.end(), [&](int i, int j) {
        return score[i] > score[j];
    });


    vector<vector<int>> temp;
    for (int i = 0; i < popuSize; i++)temp.push_back(ans[pos[i]]);

    getBest(temp);
    return temp;
}

vector<vector<int>> only_EA(vector<vector<int>> P, double mutate_rate) {

    vector<vector<int>> ans;

    vector<double> score = eval_vector(P);
    int n_p = popuSize;
    double tot = 0;
    for (auto e:score)tot += e;
    for (int i = 0; i < n_p; ++i)score[i] /= tot;
    vector<int> pos;
    for (int i = 0; i < n_p; i++)pos.push_back(i);

    sort(pos.begin(), pos.end(), [&](int i, int j) {
        return score[i] > score[j];
    });
    for (int i = 0; i < n_p / 10; i++)ans.push_back(P[pos[i]]);
    score = rank_score(pos);

    int cnt = ans.size();

    while (cnt < 2 * n_p) {

        int i = wheel(score);
        int j = wheel(score);
        while (j == i)j = wheel(score);
//        ans.push_back(P[i]);ans.push_back(P[j]);cnt+=2;
        ans.push_back(order_crossover(P[i], P[j]));
        ans.push_back(order_crossover(P[j], P[i]));
        cnt += 2;

        if (rand() % 100 * 1.0 / 100.0 < mutate_rate) {
            ans.push_back(mutation(P[i]));
            ans.push_back(mutation(P[j]));
            cnt += 2;
            ans.push_back(edge_mutation(P[i]));
            ans.push_back(edge_mutation(P[j]));
            cnt += 2;
        }
    }

    score = eval_vector(ans);
    score = robin(score);
    for (int i = popuSize; i != ans.size(); i++)
        pos.push_back(i);
    sort(pos.begin(), pos.end(), [&](int i, int j) {
        return score[i] > score[j];
    });
    vector<vector<int>> temp;
    for (int i = 0; i < n_p; i++)temp.push_back(ans[pos[i]]);
    return temp;
}


vector<vector<int>> init_p(int p) {
    bestv = 1.0 / 2047483648.0;
    vector<vector<int>> ans;
    vector<int> temp;
    for (int i = 0; i < n; ++i)temp.push_back(i);
    for (int i = 0; i < p; ++i) {
        random_shuffle(temp.begin(), temp.end());
        ans.push_back(temp);
    }
    return ans;
}

void showtime(vector<double> t) {
    double b = 100000, w = 0, m = 0;
    for (auto ele:t) {
        b = min(b, ele);
        w = max(w, ele);
        m += ele;
    }
    cout <<"best search time: "<< b <<"s"<< endl;
    cout << "worst search time: "<< w <<"s"<< endl;
    cout << "average search time: "<<m / t.size() <<"s"<< endl;
}


int main() {
    memset(G, 0, sizeof(G));
    freopen(dir, "r", stdin);
    for (int i = 0; i < n; ++i) {
        int temp;
        scanf("%d %lf %lf", &temp, &x[i], &y[i]);
    }
    fclose(stdin);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            G[i][j] = dist(i,j);
            G[j][i] = G[i][j];
        }
    }
    
    vector<double> log;
    vector<double> time_log;
    int bestCnt = 0;
    vector<double> best_log;
    for (int i = 0; i < experiment_times; i++) {
        
        hybrid_generation = 100;
        clock_t startTime = clock();
        vector<vector<int>> p0 = init_p(popuSize);
        bestv = 1.0 / 2047483648.0;
		
		if(PRE_EA) 
        for (int pre_EA = 0; pre_EA < 100; pre_EA++)
            p0 = only_EA(p0, m_rate);

        getBest(p0);
        if(PRE_EA)
        cout << "shortest length in the initial population: " << 1.0 / bestv<< endl;
        for (int g = 0; g < hybrid_generation; g++) {

            p0 = hybrid_nextGeneration(p0, m_rate);
            getBest(p0);
            cout << "In generation " << g+1 << ", best so far: " << 1.0 / bestv <<  endl;
            Stat(eval_vector(p0));
            if (1.0 / bestv < bestRoute)break;
            else hybrid_generation++;
            if ((clock() - startTime) * 1.0 / CLOCKS_PER_SEC > TIME_LIMIT)break;

        }
        time_log.push_back((clock() - startTime) * 1.0 / CLOCKS_PER_SEC);
        getBest(p0);
        cout << "Final best route length: " << 1 / bestv << endl;
        cout << "Time : " << time_log.back() << endl;
        best_log.push_back(1.0 / bestv);
        if (1.0 / bestv < bestRoute)bestCnt++;
    }
    cout << "BestCnt = " << bestCnt << endl;

    showtime(time_log);

    return 0;
}
