import copy

def get_user_table():
    matrix = []
    row_index = 0
    
    print("Вводьте кількість елементів для кожного рядка.")
    print("Натисніть Enter (порожній рядок), щоб завершити.")

    while True:
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

def reduce_matrix_by_n(table, n):
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

def append_to_row(original_table, row_index, text):
    """
    Додає текст в кінець рядка 'row_index' у таблиці 'table'.
    Якщо row_index вказує на новий рядок (дорівнює довжині таблиці),
    створюється новий рядок.
    """
    if row_index > 0:
        prev_len = len(original_table[row_index - 1])
        curent_len = len(original_table[row_index]) if row_index < len(original_table) else 0
        
        if curent_len + 1 > prev_len:
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

def table_expand(table_1, table_2):
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
        variant = append_to_row(table_1, i, element_to_move)
        if variant is not None:
            if check_lattice_condition(variant):
                list_of_new_table1s.append(variant)

    return list_of_new_table1s, next_table_2 

def generate_all_tableaux(start_table1, start0_table2):
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
    while current_source_t2 and current_source_t2[0]:
        
        print(f"--- КРОК {step_counter}: Варіантів у роботі {len(current_variants)} ---")
        
        next_generation_variants = []
        next_source_t2_ref = None

        # Проходимо по кожному варіанту T1, який ми отримали на попередньому кроці
        for table1_variant in current_variants:
            
            # Викликаємо table_expand для поточної пари (t1, t2)
            new_variants, reduced_t2 = table_expand(table1_variant, current_source_t2)
            
            # Додаємо знайдені валідні варіанти до списку наступного покоління
            next_generation_variants.extend(new_variants)
            
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

def process_one_generation(current_variants, current_source_t2):
    """
    Приймає список варіантів T1 та поточну T2.
    Застосовує expand_step до КОЖНОГО варіанту з current_variants using current_source_t2.
    
    Повертає:
    1. next_generation_variants (новий, більший список таблиць)
    2. next_source_t2 (оновлена T2, зменшена на 1 елемент)
    """
    next_generation_variants = []
    next_source_t2_ref = None

    # Якщо варіантів немає, то і робити нічого
    if not current_variants:
        return [], current_source_t2

    # Проходимо по кожному варіанту таблиці 1
    for table___1 in current_variants:
        # Викликаємо нашу функцію expand_step
        # Вона повертає список нових варіантів (new_vars) і зменшену T2 (reduced_t2)
        new_vars, reduced_t2 = table_expand(table___1, current_source_t2)
        
        # Додаємо знайдені варіанти до загального списку наступного покоління
        next_generation_variants.extend(new_vars)
        
        # Зберігаємо посилання на зменшену T2.
        # Оскільки ми застосовуємо ТОЙ САМИЙ елемент до всіх таблиць,
        # reduced_t2 буде однаковою для всіх ітерацій циклу.
        # Ми просто беремо останню (або першу).
        next_source_t2_ref = reduced_t2

    return next_generation_variants, next_source_t2_ref

def show_table(name):
    #print(name)
    matrix = name
    row_index = 0
    for row in matrix:
        print(row)

table1 = get_user_table()
table2 = get_user_table()
#table3 = append_to_row(table1, 1, "ho-ho-ho")

#table1[0][1] = "a1"

print("Результати")
print("\ntable1")
show_table(table1)
print("\ntable2")
show_table(table2)

#print("\nДодаємо")
#show_table(table3)

t_table2 = fill_labels(table2)
list_of_table1s1, next_table2 = table_expand(table1, t_table2) #table_multiplication
print("\nnext_table2")
show_table(next_table2)
print("\nlist_of_table1s1")
show_table(list_of_table1s1)

#list_of_table1s2, a = process_one_generation (list_of_table1s1, next_table2)
#print("\nlist_of_table1s2")
#show_table(list_of_table1s2)

list_of_table1s2 = generate_all_tableaux (table1, table2)
print("\nlist_of_table1s2")
show_table(list_of_table1s2)

    # 5. Виводимо результат
#print("\nСтворена таблиця:")
#for row in matrix:
    #print(row)
#print("Результати")
#print("table1", table1)
#print("table2", table2)
