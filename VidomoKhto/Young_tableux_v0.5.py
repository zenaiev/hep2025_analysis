import copy

while True:
    try:
        user_input_n = input("Симетрія Su(n), введіть n = ")
        user_n = int(user_input_n)
        
        if user_n > 0:
            break
        else:
            print("n - натуральне число\nСпробуте ще раз, n = ")
    except ValueError:
        print("Помилка. n має бути натуральним числом. Ведіть натуральне число")

def get_user_table(n):
    matrix = []
    row_index = 0
    
    print("Вводьте кількість елементів для кожного рядка.")
    print("Натисніть Enter (порожній рядок), щоб завершити.")

    while True:
        if row_index >= n-1:
            print(f"Досягнуто максимально дозволеної кількості рядків в Su({n}). Таблиця завершена.")
            break
        
        try:
            # 1. Читаємо введений рядок
            user_input = input(f"Кількість елементів для рядка {row_index}: ")
            
            # 2. Якщо рядок порожній — виходимо з циклу
            if user_input == "":
                break
                
            # 3. Конвертуємо в число
            num_cols = int(user_input)
            
            if num_cols < 0:
                print("Кількість не може бути від'ємною. Спробуйте ще раз.")
                continue
            
            if matrix and num_cols > len(matrix[-1]):
                print(f"Максимальна дозволена довжина:{len(matrix[-1])}. Спробуйте ще раз.")
                continue

            # 4. Створюємо та додаємо рядок з нулів
            matrix.append([0] * num_cols)
            row_index += 1

        except ValueError:
            print("Помилка: Введіть число або натисніть Enter для завершення.")
    
    return matrix

def reduce_matrix_by_n(n, table):
    """
    Якщо кількість рядків досягає n, видаляє перші елементи всіх рядків,
    поки останній рядок не стане порожнім (і не видалиться).
    Порожні рядки видаляються з таблиці.
    """
    # Поки висота таблиці дорівнює n (або більше, про всяк випадок)
    while len(table) >= n:
        
        # 1. Видаляємо перший елемент (index 0) у кожного рядка
        for row in table:
            if row: # Перевірка, щоб не робити pop з порожнього (хоча такого не має бути при n)
                row.pop(0)
        
        # 2. Очищення: Видаляємо рядки, які стали порожніми
        # Ми перезаписуємо table, залишаючи тільки не порожні рядки
        table[:] = [row for row in table if row]
        
        # Цикл while перевірить знову: якщо після чистки висота все ще n,
        # значить треба зрізати ще один стовпчик.
        
    return table

def append_to_row(n, original_table, row_index, text):
    """
    Додає текст в кінець рядка 'row_index' у таблиці 'table'.
    Якщо row_index вказує на новий рядок (дорівнює довжині таблиці),
    створюється новий рядок.
    """
    if row_index >= n:
        return None
    
    if row_index > 0:
        prev_len = len(original_table[row_index - 1])
        curent_len = len(original_table[row_index]) if row_index < len(original_table) else 0
        
        if curent_len + 1 > prev_len:
            return None
    
    target_col = len(original_table[row_index]) if row_index < len(original_table) else 0
    
    
    for r in range(row_index):
        prev_val = original_table[r][target_col]
        if isinstance(prev_val, str) and prev_val == text:
            return None
    
    new_table=copy.deepcopy(original_table)
    
    current_length = len(new_table)

    if row_index < current_length:
        # Варіант 1: Рядок існує. Просто додаємо елемент в кінець.
        new_table[row_index].append(text)
        
    elif row_index == current_length:
        # Варіант 2: Це "уявний порожній рядок". 
        # Створюємо новий список з цим елементом і додаємо до таблиці.
        new_table.append([text])
        
    else:
        # Варіант 3 (про всяк випадок): Індекс занадто великий (розрив у таблиці)
        print(f"Помилка: Неможливо додати у рядок {row_index}. У таблиці лише {current_length} рядків.")
    
    """
    if len(new_table) >= n:
        new_table = reduce_matrix_by_n(n, new_table)
    """
    
    return new_table

def fill_labels(table):
    new_table = copy.deepcopy(table)

    # Заповнюємо цю копію мітками a1, a2...
    for i in range(len(new_table)):
        label = f"a{i + 1}"
        for j in range(len(new_table[i])):
            # Перезаписуємо значення (наприклад, нулі) на мітки
            new_table[i][j] = label
    return new_table

def check_lattice_condition(table):
    """
    Перевіряє умову решітчастої перестановки (Lattice Word):
    При зчитуванні справа-наліво, зверху-вниз, кількість a(n) 
    ніколи не повинна перевищувати кількість a(n-1).
    """
    # Словник для підрахунку кількості зустрінутих a1, a2, a3...
    # Ключ - номер (int), Значення - кількість (int)
    counts = {}

    # 1. Проходимо по рядках згори донизу (0, 1, 2...)
    for row in table:
        
        # 2. Проходимо по клітинках справа наліво
        for cell in reversed(row):
            
            # Перевіряємо, чи це наші мітки (a1, a2...). Ігноруємо числа (0, 10, 20).
            if isinstance(cell, str) and cell.startswith('a'):
                try:
                    # Витягуємо номер: "a2" -> 2
                    idx = int(cell[1:])
                except ValueError:
                    continue # Якщо раптом там щось дивне, пропускаємо
                
                # Збільшуємо лічильник для цього номера
                counts[idx] = counts.get(idx, 0) + 1
                
                # === ГОЛОВНА ПЕРЕВІРКА ===
                # Якщо це a1 - перевіряти нічого не треба (він найголовніший)
                # Якщо це a2, a3... - перевіряємо відносно попереднього (a(n-1))
                if idx > 1:
                    prev_count = counts.get(idx - 1, 0)
                    
                    # Умова: якщо поточних (напр. a2) стало БІЛЬШЕ ніж попередніх (a1) -> Помилка
                    if counts[idx] > prev_count:
                        return False # Таблиця не валідна

    return True # Якщо пройшли все і не впали - таблиця валідна

def table_expand(n, table_1, table_2):
    # Перевірка, чи є що переміщати
    if not table_2 or not table_2[0]:
        return [], table_2

    # Беремо елемент (правий верхній з вже заповненої таблиці)
    element_to_move = table_2[0][-1]
    
    # --- КРОК 3: Створюємо "залишок" другої таблиці ---
    next_table_2 = copy.deepcopy(table_2)
    next_table_2[0].pop() # Видаляємо елемент, який забрали
    
    if not next_table_2[0]:
        next_table_2.pop(0) # Видаляємо порожній рядок

    # --- КРОК 4: Генеруємо варіанти першої таблиці ---
    list_of_new_table1s = []
    
    for i in range(len(table_1) + 1):
        # append_to_row_new сама робить deepcopy table_1
        variant = append_to_row(n, table_1, i, element_to_move)
        
        if variant is not None:
            if check_lattice_condition(variant):
                list_of_new_table1s.append(variant)

    return list_of_new_table1s, next_table_2 

def generate_all_tableaux(n, start_table1, start0_table2):
    """
    Послідовно приєднує всі елементи з start_table2 до start_table1,
    перебираючи всі можливі валідні варіанти.
    """
    
    start_table2 = fill_labels(start0_table2)
    
    # Поточний список варіантів першої таблиці. Спочатку тут тільки одна таблиця.
    current_variants = [start_table1]
    
    # Поточний стан другої таблиці (джерело елементів)
    current_source_t2 = start_table2

    step_counter = 1

    # Цикл триває, доки джерело (t2) не стане порожнім
    # Ми перевіряємо, чи є там рядки і чи є в першому рядку елементи
    print("\n")
    while current_source_t2 and current_source_t2[0]:
        
        print(f"КРОК {step_counter}: Варіантів у роботі {len(current_variants)}")
        
        next_generation_variants = []
        next_source_t2_ref = None

        # Проходимо по кожному варіанту T1, який ми отримали на попередньому кроці
        for table1_variant in current_variants:
            
            # Викликаємо table_expand для поточної пари (t1, t2)
            new_variants, reduced_t2 = table_expand(n, table1_variant, current_source_t2)
            
            #Фільтрація дублікатів
            for variant in new_variants:
                if variant not in next_generation_variants:
                    next_generation_variants.append(variant)
            
            # Зберігаємо зменшену таблицю T2 для наступного кроку циклу while.
            # Вона однакова для всіх варіантів у цьому циклі, тому просто перезаписуємо.
            next_source_t2_ref = reduced_t2

        # Оновлюємо змінні для наступної ітерації
        current_variants = next_generation_variants
        current_source_t2 = next_source_t2_ref
        step_counter += 1

        # Якщо на якомусь етапі не залишилось валідних варіантів - перериваємо
        if not current_variants:
            print("Валідних варіантів більше немає (глухий кут).")
            break

    return current_variants

def finalize_tables(n, tables):
    f_tables = []
    for t in tables:
        reduced_t = reduce_matrix_by_n(n, t) #!!!!!!!!!!!!!!!!
        #if reduced_t not in f_tables:
        f_tables.append(reduced_t)
    return f_tables

def clean_tables(tables):
    c_tables = []
    for table in tables:
        c_tables.append([[0] * len(row) for row in table])
    return c_tables

def show_table(matrix):
    for i in range(len(matrix)):
        print(f"\ntablle {i + 1}")
        if len(matrix[i]) == 0: print("1")
        else:
            for ii in range(len(matrix[i])):
                for iii in range(len(matrix[i][ii])):
                    if isinstance(matrix[i][ii][iii], int): print("[ ]", end="  ")
                    if isinstance(matrix[i][ii][iii], str): print(f"[{matrix[i][ii][iii]}]", end=" ")
                print("")

def tables_dimensions(n, tables):
    dymesions_of_tables = []
    for table in tables:
        if len(table) == 0:
            dymesions_of_tables.append(1)
            continue
        dymesion_of_table = 1
        dymesion_of_table_numerator = 1
        dymesion_of_table_denominator = 1
        for i in range(len(table[0])):
            j = 0
            for ii in range(len(table)):
                if i < len(table[len(table) - ii - 1]):
                    j += 1
                    slot_in_numerator = n + i - (len(table) - ii - 1)
                    slot_in_denominator = (len(table[len(table) - ii - 1]) - i - 1) + j
                    dymesion_of_table_numerator = dymesion_of_table_numerator * slot_in_numerator
                    dymesion_of_table_denominator = dymesion_of_table_denominator * slot_in_denominator
                else: continue
        dymesion_of_table = dymesion_of_table_numerator / dymesion_of_table_denominator
        dymesions_of_tables.append(dymesion_of_table)
        #print(dymesion_of_table_numerator)
        #print(dymesion_of_table_denominator)
    return dymesions_of_tables

start_tables = []
start_tables.append(get_user_table(user_n))
start_tables.append(get_user_table(user_n))

#table1[0][1] = "a1"

print("Вхідні дані")
show_table(start_tables)

list_of_table1s2 = generate_all_tableaux (user_n, start_tables[0], start_tables[1])
final_tables = finalize_tables(user_n, list_of_table1s2)
cleaned_tables = clean_tables(final_tables)
print("\nОтримуємо")
#show_table(list_of_table1s2)
show_table(final_tables)
#show_table(cleaned_tables)

d_in = tables_dimensions(user_n, start_tables)
d_out = tables_dimensions(user_n, final_tables)

print(f"\nОтримуємо\n{d_in[0]} \u2A02  {d_in[1]} =", end="")
sum = 0
for jjj in range(len(d_out) - 1):
    print(f"  {d_out[jjj]} \u2A01", end="")
    sum = sum + d_out[jjj]
print(f"  {d_out[-1]}")
print(d_in[0] * d_in[1], end=" = ")
print(sum + d_out[-1])
