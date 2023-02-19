#include <bits/stdc++.h>
#define ld long double
#define pb push_back

using namespace std;

const ld MaxCost = 1e7, MaxArea = 2000, MinCost = 1e4, MinArea = 10;
const int MaxnR = 20, MaxnBR = 10, MaxswDist = 50, MinnR = 1, MinnBR = 1, MinswDist = 0;

int pos (int x, int y, int z){
    return (x * y + z);
}

void ginv(vector <ld> matr, int n, int m, vector <ld> &d)
{
    // Обернена матриця по Муру-Пенроузу (або псевдообернена)
    int u,v ;
    matr[0] = sqrt(matr[0]);
    for (int i = 1; i < n; ++i) matr[i] /= matr[0];
    for (int i = 1; i < n; ++i)
    {
        for (int j = 0; j < i; ++j)
        {
            matr[pos(i, n, i)] -= matr[pos(j, n, i)] * matr[pos(j, n, i)];
        }
        matr[pos(i, n, i)] = sqrt(matr[pos(i, n, i)]);
        for (int j = i + 1; j < n; ++j)
            {
                for (int k = 0; k < i; ++k)
                    matr[pos(i, n, j)] -= matr[pos(k, n, i)] * matr[pos(k, n, j)];
                    matr[pos(i, n, j)] /= matr[pos(i, n, i)];
            }
    }
    for (int j = 0; j < m; ++j)
    {
        d[j] /= matr[0];
        for (int i = 1; i < n; ++i)
        {
            for (int k = 0; k < i; ++k) d[pos(i, m, j)] -= matr[pos(k, n, i)] * d[pos(k, m, j)];
            d[pos(i, m, j)] /= matr[pos(i, n, i)];
        }
    }
    for (int j = 0; j < m; ++j)
    {
        d[pos(n - 1, m, j)] /= matr[pos(n, n, -1)];
        for (int k = n - 1; k >= 1; --k)
        {
            for (int i = k; i < n; ++i)
            {
                d[pos(k - 1, m, j)] -= matr[pos(k - 1, n, i)] * d[pos(i, m, j)];
            }
            d[pos(k - 1, m, j)] /= matr[pos(k - 1, n, k - 1)];
        }
    }
    return;
}

void leastsq(vector <ld> x, vector <ld> y, int m, int n, vector <ld> &w)
{
    // метод найменших квадратів для побудови лінійної регресії від декількох змінних
    ld p;
    vector <ld> b((m + 1) * (m + 1));
    b[pos(m + 1, m + 1, -1)] = n ;
    for (int j = 0; j < m; ++j)
    {
        p = 0;
        for (int i = 0; i < n; ++i)
            p = p + x[pos(j, n, i)];
            b[pos(m, m + 1, j)] = p;
            b[pos(j, m + 1, m)] = p;
    }
    for (int i = 0; i < m; ++i){
        for (int j = i; j < m; ++j)
        {
            p = 0;
            for (int k = 0; k < n; ++k) p += x[pos(i, n, k)] * x[pos(j, n, k)];
            b[pos(j, m + 1, i)] = p;
            b[pos(i, m + 1, j)] = p;
        }
    }
    w[m] = 0;
    for (int i = 0; i < n; ++i) w[m] += y[i];
    for (int i = 0; i < m; ++i)
    {
        w[i] = 0;
        for (int j = 0; j < n; ++j) w[i] += x[pos(i, n, j)] * y[j];
    }
    ginv(b, m + 1, 1, w);
    return;
}

ld eps = 1e-4;


struct estate {
    ld area, cost; // площа, ціна в $
    int nR, nBR, swDist; // кількість кімнат, кількість санвузлів, відстань до метро
    bool liv; // комерційна чи житлова
};

int check (estate x, bool phase, int p = 5){
    if (x.area < MinArea || x.area > MaxArea) return -1;
    if (x.nR < MinnR || x.nR > MaxnR) return -1;
    if (x.nBR < MinnBR || x.nBR > MaxnBR) return -1;
    if (x.swDist < MinswDist || x.swDist > MaxswDist) return -1;
    if (p < 5) return -1;
    if (phase == 1) return 1;
    if (x.cost < MinCost || x.cost > MaxCost) return -1;
    return 1;
}

void train(vector <estate> v, vector <ld> &w){
    vector <ld> x, y;
    /* створюємо набiр векторiв даних, лiнiйна регресiя не працює на
    лінійно-залежних векторах, тому додамо епсілон помножений на значення та
    випадкове ціле число із відрізку [-5; 5] за нормального розподілу
    */
    for (int i = 0; i < min((int)v.size(), 1000); ++i) x.pb(v[i].area + eps * v[i].area * ((rand() % 11) - 5));
    for (int i = 0; i < min((int)v.size(), 1000); ++i) x.pb(v[i].nR + eps * v[i].nR * ((rand() % 11) - 5));
    for (int i = 0; i < min((int)v.size(), 1000); ++i) x.pb(v[i].nBR + eps * v[i].nBR * ((rand() % 11) - 5));
    for (int i = 0; i < min((int)v.size(), 1000); ++i) x.pb(v[i].swDist + eps * v[i].swDist * ((rand() % 11) - 5));
    for (int i = 0; i < min((int)v.size(), 1000); ++i) y.pb(v[i].cost + eps * v[i].cost * ((rand() % 11) - 5));
    leastsq(x, y, 4, min((int)v.size(), 1000), w);
}

template<typename T>
int input (T &x){
    cin >> x;
    while (cin.fail()){
        cin.clear();
        cin.ignore(numeric_limits<streamsize>::max(),'\n');
        cout << "\n!!! Помилка типу даних! Почнiть спочатку !!!\n\n";
        return -1;
    }
    return 1;
}

void what_price(estate x, vector <ld> w){
    cout << fixed << setprecision(0) << x.area * w[0] + x.nR * w[1] + x.nBR * w[2] + x.swDist * w[3] + w[4] << "$\n";
}

signed main (){
    srand(time(0));
    setlocale(LC_ALL, "ukr");
    /* У моделі прогнозування припущено, що залежність ціни від параметрів квартири є лінійна. За допомогою
    регресії обирається рівняння виду y = X[i] * w[i] + const
    */
    vector <estate> v[2];
    vector <ld> w[2]; w[0].resize(5), w[1].resize(5);
    estate temp;
    cout << "---Training phase---\n";
    cout << '\n';
    cout << "!!! Для прогнозування потрiбно не менше 5 прикладiв типу нерухомостi (комерцiйна чи житлова), який буде прогнозуватися !!!\n\n";
    while (true){
        cout << "Введiть площу нерухомостi у м^2, для завершення введiть 0:\n";
        int valin = input(temp.area);
        if (valin == -1) continue;
        if (temp.area == 0) break;
        cout << "Кiлькiсть кiмнат:\n";
        valin = input(temp.nR);
        if (valin == -1) continue;
        cout << "Кiлькiсть сан-вузлiв:\n";
        valin = input(temp.nBR);
        if (valin == -1) continue;
        cout << "Вiдстань вiд найближчої станцiї у хвилинах пiшки:\n";
        valin = input(temp.swDist);
        if (valin == -1) continue;
        cout << "Тип нерухомостi, введiть 0 для комерцiйної або 1 для житлової:\n";
        valin = input(temp.liv);
        if (valin == -1) continue;
        cout << "Вартiсть нерухомостi у дол.:\n";
        valin = input(temp.cost);
        if (valin == -1) continue;
        int correct = check(temp, 0);
        if (correct == 1) v[temp.liv].pb(temp);
        else if (correct == -1) cout << "\n!!! Було введено некоректнi данi, екземпляр не було додано до списку !!!\n\n";
    }
    bool fl = 0;
    for (int i = 0; i < 2; ++i){
        if (v[i].size()) train(v[i], (w[i])), fl = 1;
    }
    if (!fl) cout << "Недостатньо даних для прогнозування будь-якого типу нерухомості\n";
    cout << "---Forecast phase---\n";
    ld res = 0;
    while (true){
        cout << "Введiть площу нерухомостi у м^2, для завершення введiть 0:\n";
        int valin = input(temp.area);
        if (valin == -1) continue;
        if (temp.area == 0) break;
        cout << "Кiлькiсть кiмнат:\n";
        valin = input(temp.nR);
        if (valin == -1) continue;
        cout << "Кiлькiсть сан-вузлiв:\n";
        valin = input(temp.nBR);
        if (valin == -1) continue;
        cout << "Вiдстань вiд найближчої станцiї у хвилинах пiшки:\n";
        valin = input(temp.swDist);
        if (valin == -1) continue;
        cout << "Тип нерухомостi, введiть 0 для комерцiйної або 1 для житлової:\n";
        valin = input(temp.liv);
        if (valin == -1) continue;
        int correct = check(temp, 1, v[temp.liv].size());
        if (correct == 1) {
            cout << "Прогнозована цiна нерухомостi у дол.:\n";
            what_price(temp, w[temp.liv]);
        }
        else if (correct == -1) cout << "\n!!! Було введено некоректнi данi, для даного екземпляру не можна створити прогноз !!!\n\n";
    }
return 0;
}
