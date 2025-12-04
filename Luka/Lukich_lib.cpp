#include "Lukich_lib.h"

#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <boost/multiprecision/cpp_int.hpp>

using namespace std;
using boost::multiprecision::cpp_int;

// ---------- НСД для великих цілих ----------
cpp_int bigint_gcd(cpp_int a, cpp_int b) {
    while (b != 0) {
        cpp_int t = a % b;
        a = b;
        b = t;
    }
    return a;
}

// ---------- Розмірність представлення по діаграмі Юнга ----------
cpp_int compute_dimension(int N, const vector<int>& row) {
    int r = (int)row.size();
    int totalBoxes = 0;
    for (int x : row) totalBoxes += x;

    vector<cpp_int> num; // чисельники N + j - i
    vector<cpp_int> den; // знаменники (довжини гачків)
    num.reserve(totalBoxes);
    den.reserve(totalBoxes);

    for (int i = 0; i < r; ++i) {          // індекс рядка (0 зверху)
        for (int j = 0; j < row[i]; ++j) { // індекс стовпця (0 зліва)
            //кількість клітинок праворуч у тому ж рядку
            int right = row[i] - j - 1;
            //перебираємо рядки знизу (k > i) й перевіряємо, чи в них є клітинка в стовпці j
            int below = 0;
            for (int k = i + 1; k < r; ++k) {
                if (row[k] > j) ++below;
            }
            int hook = 1 + right + below;

            cpp_int n = cpp_int(N + j - i); // N + (j - i)
            num.push_back(n);
            den.push_back(cpp_int(hook));
        }
    }

    // скорочуємо дроби
    int m = (int)den.size();
    int k = (int)num.size();
    //для кожного знаменника намагаємось поділити його на НСД з чисельниками
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < k; ++j) {
            if (den[i] == 1) break;
            cpp_int g = bigint_gcd(den[i], num[j]);
            if (g > 1) {
                den[i] /= g;
                num[j] /= g;
            }
        }
    }

    cpp_int dim = 1;
    for (const auto& x : num) dim *= x;
    return dim;
}

// ---------- Структура таблиці Юнга з мітками (для продукту) ----------
struct Tableau {
    vector<int> rowLen;             // довжини рядків
    vector<vector<int>> labels;     // матриця тієї ж форми; 0 — клітини з A, >0 — з B
};

// ---------- Інфо про одну форму  λ ----------
struct Info {
    int mult = 0;           // кратність m_λ
    vector<int> shape;      // сама форма діаграми λ
};

// ключ для map по формі діаграми
string shape_key(const vector<int>& row) {
    string s;
    for (size_t i = 0; i < row.size(); ++i) {
        if (i) s.push_back(',');
        s += to_string(row[i]);
    }
    return s;
}

// чи можна додати нову клітинку з міткою label у таблицю T в рядок rowIndex для групи SU(N).
bool can_add_box(const Tableau& T, int N, int rowIndex, int label) {
    int numRows = (int)T.rowLen.size(); // кл-ть зараз
    bool newRow = (rowIndex == numRows); // чи додаємо новий рядок
    //Обчислюємо, у який стовпець потрапить нова клітина:
    //якщо це новий рядок — це 1-ша клітинка, стовпець 1;
    //інакше — в кінець рядка: довжина рядка + 1.
    int col = newRow ? 1 : (T.rowLen[rowIndex] + 1); 

    // форма діаграми: довжина рядка не може перевищувати довжину рядка вище
    if (!newRow && rowIndex > 0 && T.rowLen[rowIndex - 1] < col)
        return false;

    // перевірка стовпця: заборонено повторну мітку в одному стовпці
    int heightBefore = 0;
    for (int i = 0; i < numRows; ++i) {
        if (T.rowLen[i] >= col) {
            ++heightBefore;
            if (T.labels[i][col - 1] == label) return false;
        }
    }

    // обмеження SU(N): висота стовпця <= N
    int heightAfter = heightBefore + 1;
    if (heightAfter > N) return false;

    return true;
}

// ---------- Рекурсивне додавання k клітин з однією міткою (без дублів) ----------
void add_boxes_label_rec(const Tableau& T, int N, int label, int remaining,
                         int minRow, vector<Tableau>& out) {
    if (remaining == 0) {
        out.push_back(T);
        return;
    }

    int numRows = (int)T.rowLen.size();

    // додаємо лише в рядки з iндексом >= minRow (щоб уникнути перестановок)
    for (int r = minRow; r <= numRows; ++r) {
        if (!can_add_box(T, N, r, label)) continue;

        Tableau T2 = T;
        if (r == numRows) {
            // новий рядок
            T2.rowLen.push_back(1);
            T2.labels.push_back(vector<int>(1, label));
        } else {
            T2.rowLen[r]++;
            T2.labels[r].push_back(label);
        }

        // наступну клiтину з тим самим label можна ставити
        // лише в рядок r або нижче
        add_boxes_label_rec(T2, N, label, remaining - 1, r, out);
    }
}

// обгортка, щоб не змiнювати сигнатуру
void add_boxes_label(const Tableau& T, int N, int label, int remaining,
                     vector<Tableau>& out) {
    add_boxes_label_rec(T, N, label, remaining, 0, out);
}

// ---------- Умова Яманоучі (балотна) для послідовності міток ----------
bool is_yamanouchi(const Tableau& T, int numLabels) {
    vector<int> seq;
    // читаємо мітки: рядки зверху вниз, у кожному справа наліво
    int rows = (int)T.rowLen.size();
    for (int i = 0; i < rows; ++i) {
        for (int j = T.rowLen[i] - 1; j >= 0; --j) {
            int lab = T.labels[i][j];
            if (lab > 0) seq.push_back(lab); // тільки клітини з B
        }
    }

    vector<int> cnt(numLabels + 2, 0);

    for (int x : seq) {
        cnt[x]++;
        // на кожному префіксі: #a_i >= #a_{i+1}
        for (int k = 1; k < numLabels; ++k) {
            if (cnt[k] < cnt[k + 1]) return false;
        }
    }
    return true;
}

// ---------- Добуток двох діаграм: A ⊗ B ----------
map<string, Info> tensor_product_two(int N,
                                     const vector<int>& A,
                                     const vector<int>& B)
{
    // Початкова таблиця: тільки діаграма A, усі мітки 0
    Tableau base;
    base.rowLen = A;
    base.labels.resize(A.size());
    for (size_t i = 0; i < A.size(); ++i) {
        base.labels[i].assign(A[i], 0);
    }

    // Послідовно додаємо ряди B: a1, a2, ...
    vector<Tableau> current;
    current.push_back(base);

    int rB = (int)B.size();
    for (int lab = 1; lab <= rB; ++lab) {
        vector<Tableau> next;
        int count = B[lab - 1]; // скільки клітин з міткою lab
        for (const auto& T : current) {
            vector<Tableau> tmp;
            add_boxes_label(T, N, lab, count, tmp);
            next.insert(next.end(), tmp.begin(), tmp.end());
        }
        current.swap(next);
    }

    // Фільтруємо за умовою Яманоучі
    vector<Tableau> valid;
    for (const auto& T : current) {
        if (is_yamanouchi(T, rB)) {
            valid.push_back(T);
        }
    }

    // Групуємо за формою діаграми: мультиплічності
    map<string, Info> result;

    for (const auto& T : valid) {
        string key = shape_key(T.rowLen);
        auto& info = result[key];
        if (info.shape.empty()) info.shape = T.rowLen;
        info.mult += 1;
    }

    return result;
}

int compute() {
    cout << "==============================================\n";
    cout << "   КАЛЬКУЛЯТОР ТЕНЗОРНОГО ДОБУТКУ ДIАГРАМ ЮНГА\n";
    cout << "       ДЛЯ ПРЕДСТАВЛЕНЬ SU(N)\n";
    cout << "==============================================\n\n";

    int N;
    cout << "Введiть N (для SU(N)): ";
    if (!(cin >> N)) return 0;

    int K;
    cout << "Введiть кiлькiсть дiаграм K (>=1): ";
    cin >> K;
    if (K <= 0) {
        cout << "K повинно бути >= 1\n";
        return 0;
    }

    vector<vector<int>> diagrams(K);

    for (int idx = 0; idx < K; ++idx) {
        int r;
        cout << "\nДiаграма #" << (idx + 1) << ":\n";
        cout << "  Введiть кiлькiсть рядкiв r: ";
        cin >> r;
        if (r <= 0) {
            cout << "Помилка: кiлькiсть рядкiв повинна бути додатньою.\n";
            return 0;
        }

        diagrams[idx].resize(r);
        cout << "  Введiть довжини рядкiв (через пробiл, зверху донизу): ";
        for (int i = 0; i < r; ++i) {
            cin >> diagrams[idx][i];
        }

        // --- ПЕРЕВIРКИ ФОРМИ ДIАГРАМИ ---

        for (int i = 0; i < r; ++i) {
            // 1) довжина рядка не може перевищувати кiлькiсть рядкiв
            if (diagrams[idx][i] > r) {
                cout << "Помилка у дiаграмi #" << (idx + 1)
                     << ": довжина рядка " << i
                     << " (" << diagrams[idx][i]
                     << ") не може перевищувати кiлькiсть рядкiв r = "
                     << r << "\n";
                return 0;
            }

            // 2) кожний наступний рядок не довший за попереднiй
            if (i > 0 && diagrams[idx][i] > diagrams[idx][i - 1]) {
                cout << "Помилка у дiаграмi #" << (idx + 1)
                     << ": довжина нижнього рядка " << i
                     << " (" << diagrams[idx][i]
                     << ") не може бути бiльшою за довжину рядка вище "
                     << (i - 1) << " (" << diagrams[idx][i - 1] << ")\n";
                return 0;
            }
        }
    }

    // Якщо всього одна дiаграма — просто виводимо її розмiрнiсть
    if (K == 1) {
        cpp_int dim = compute_dimension(N, diagrams[0]);
        cout << "\nЄдина дiаграма з формою (";
        for (size_t i = 0; i < diagrams[0].size(); ++i) {
            if (i) cout << ", ";
            cout << diagrams[0][i];
        }
        cout << "),  dim = " << dim << "\n";
        return 0;
    }

    // ---------- Послiдовно перемножуємо усi дiаграми ----------

    // Поточний розклад: спочатку це просто D1 з кратнiстю 1
    map<string, Info> current;
    {
        string key = shape_key(diagrams[0]);
        Info& info = current[key];
        info.shape = diagrams[0];
        info.mult = 1;
    }

    // послiдовно множимо на D2, D3, ..., DK
    for (int idx = 1; idx < K; ++idx) {
        map<string, Info> next;

        for (const auto& kv : current) {
            const Info& term = kv.second;      // поточна форма λ iз кратнiстю term.mult
            const vector<int>& shape = term.shape;

            // розклад λ ⊗ D_idx
            map<string, Info> partial = tensor_product_two(N, shape, diagrams[idx]);

            // додаємо в next, враховуючи кратнiсть term.mult
            for (const auto& kv2 : partial) {
                const Info& info2 = kv2.second;   // форма μ з кратнiстю m_μ
                Info& target = next[kv2.first];

                if (target.shape.empty()) target.shape = info2.shape;
                target.mult += term.mult * info2.mult;
            }
        }

        current.swap(next); // тепер current — розклад D1 ⊗ ... ⊗ D_idx
    }

    // ---------- Вивiд фінального розкладу ----------

    cout << "\n----------------------------------------------\n";
    cout << "Результат розкладу тензорного добутку D1 ⊗ D2 ⊗ ... ⊗ DK:\n\n";

    cpp_int sumDims = 0;

    for (const auto& kv : current) {
        const Info& info = kv.second;
        cpp_int dimRep = compute_dimension(N, info.shape);

        cout << "Форма λ = (";
        for (size_t i = 0; i < info.shape.size(); ++i) {
            if (i) cout << ", ";
            cout << info.shape[i];
        }
        cout << "):  кратнiсть m = " << info.mult
             << ",  dim(λ) = " << dimRep << "\n";

        sumDims += dimRep * info.mult;
    }

    // ---------- Перевiрка розмiрностей для всiх дiаграм ----------

    cpp_int dimInput = 1;
    for (int idx = 0; idx < K; ++idx) {
        cpp_int d = compute_dimension(N, diagrams[idx]);
        dimInput *= d;
    }

    cout << "\n----------------------------------------------\n";
    cout << "Перевiрка розмiрностей:\n";
    cout << "  ∏ dim(D_i) = " << dimInput << "\n";
    cout << "  Σ m_λ * dim(λ) = " << sumDims << "\n";

    if (dimInput == sumDims)
        cout << "  ✓ Спiвпадає (розклад коректний)\n";
    else
        cout << "  ⚠ НЕ спiвпадає (можлива помилка або занадто великий приклад)\n";

    return 0;
}
