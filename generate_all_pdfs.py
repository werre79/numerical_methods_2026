from fpdf import FPDF
import subprocess
import os

class PDF(FPDF):
    def __init__(self):
        super().__init__()
        font_dir = '/usr/share/fonts/truetype/liberation/'
        if os.path.exists(font_dir + 'LiberationSerif-Regular.ttf'):
            self.add_font('LiberationSerif', '', font_dir + 'LiberationSerif-Regular.ttf', uni=True)
            self.add_font('LiberationSerif', 'B', font_dir + 'LiberationSerif-Bold.ttf', uni=True)
            self.add_font('LiberationSerif', 'I', font_dir + 'LiberationSerif-Italic.ttf', uni=True)
            self.add_font('LiberationMono', '', font_dir + 'LiberationMono-Regular.ttf', uni=True)
        else:
            self.add_font('LiberationSerif', '', 'Arial', uni=True)
            self.add_font('LiberationSerif', 'B', 'Arial', style='B', uni=True)

    def add_title_page(self, title, lab_number):
        self.add_page()
        self.set_font('LiberationSerif', '', 14)
        
        self.cell(0, 7, 'Міністерство освіти і науки України', ln=True, align='C')
        self.cell(0, 7, 'Львівський національний університет імені Івана Франка', ln=True, align='C')
        self.cell(0, 7, 'Факультет електроніки та комп’ютерних технологій', ln=True, align='C')
        self.cell(0, 7, 'Кафедра системного програмування', ln=True, align='C')
        
        self.ln(40)
        
        self.set_font('LiberationSerif', 'B', 16)
        self.cell(0, 8, 'Звіт', ln=True, align='C')
        self.set_font('LiberationSerif', '', 14)
        self.cell(0, 8, f'Про виконання лабораторної роботи №{lab_number}', ln=True, align='C')
        self.cell(0, 8, f'«{title}»', ln=True, align='C')
        
        self.ln(40)
        
        self.set_x(120)
        self.cell(50, 7, 'Виконав:', ln=True)
        self.set_x(120)
        self.cell(50, 7, 'Студент групи ФЕІ-12, Шутяк В. В.', ln=True)
        self.ln(10)
        self.set_x(120)
        self.cell(50, 7, 'Перевірив:', ln=True)
        self.set_x(120)
        self.cell(50, 7, 'ас. Патинко А.М.', ln=True)
        
        self.set_y(-30)
        self.cell(0, 10, 'Львів 2026', ln=True, align='C')

    def chapter_title(self, title):
        self.set_font('LiberationSerif', 'B', 14)
        self.cell(0, 10, title, ln=True, align='L')
        self.ln(2)

    def chapter_body(self, text):
        self.set_font('LiberationSerif', '', 12)
        self.multi_cell(0, 8, text)
        self.ln()

    def add_script_output(self, script_path, cwd):
        result = subprocess.run(['python3', script_path], capture_output=True, text=True, cwd=cwd)
        output = result.stdout.replace('\t', '    ')
        return output

def create_lab8_report():
    pdf = PDF()
    pdf.set_auto_page_break(auto=True, margin=15)
    pdf.add_title_page('Чисельні методи розв’язування нелінійних рівнянь з одним невідомим', '8')
    
    pdf.add_page()
    pdf.chapter_title('Мета:')
    pdf.chapter_body('Вивчити однокрокові та багатокрокові ітераційні методи розв’язку нелінійних рівнянь з одним невідомим. Вивчити чисельні методи розв’язку алгебраїчних рівнянь по методу Ньютона з використанням схеми Горнера та знаходження комплексних коренів по методу Ліна.')
    
    pdf.chapter_title('Хід роботи:')
    text_hid = ('1. Табуляція трансцендентної функції f(x) = sin(x) - x/2 та збереження даних.\n'
                '2. Локалізація коренів та уточнення їх методів простої ітерації, Ньютона, Чебишева, хорд, парабол та зворотної інтерполяції.\n'
                '3. Знаходження дійсного кореня алгебраїчного рівняння третього порядку x^3 - 2x^2 + x - 2 = 0 за допомогою методу Ньютона та схеми Горнера.\n'
                '4. Знаходження комплексних коренів цього рівняння методом Ліна (Баєрстоу).')
    pdf.chapter_body(text_hid)

    pdf.chapter_title('Результати:')
    output = pdf.add_script_output('lab8.py', 'LAB_8')
    pdf.set_font('LiberationMono', '', 10)
    pdf.multi_cell(0, 5, output)
    pdf.ln()

    if os.path.exists('LAB_8/algebraic_plot.png'):
        pdf.image('LAB_8/algebraic_plot.png', w=140)
        pdf.ln()
    
    pdf.chapter_title('Висновки:')
    pdf.chapter_body('Досліджено поведінку однокрокових та багатокрокових ітераційних методів. На практиці перевірено, що методи вищих порядків, такі як метод Чебишева та метод Ньютона, збігаються до точного розв\'язку за значно меншу кількість кроків порівняно з простою ітерацією. Крім того, засвоєно підхід до оптимізації обчислення поліномів за схемою Горнера та виділення комплексних коренів методом Ліна.')
    
    pdf.output('LAB_8/Шутяк_ЧисельніМетоди_лаб8.pdf')

def create_lab9_report():
    pdf = PDF()
    pdf.set_auto_page_break(auto=True, margin=15)
    pdf.add_title_page('Метод Хука-Дживса багатовимірної оптимізації', '9')
    
    pdf.add_page()
    pdf.chapter_title('Мета:')
    pdf.chapter_body('Вивчити використання методу нульового порядку Хука-Дживса для розв’язку системи нелінійних рівнянь шляхом мінімізації спеціально побудованої цільової функції.')
    
    pdf.chapter_title('Хід роботи:')
    text_hid = ('1. Задання цільової функції Розенброка та системи двох нелінійних рівнянь.\n'
                '2. Програмна реалізація методу Хука-Дживса з етапами досліджуючого пошуку та пошуку за зразком.\n'
                '3. Тестування алгоритму на функції Розенброка.\n'
                '4. Розв\'язок системи нелінійних рівнянь шляхом мінімізації цільової функції суми квадратів нев\'язок. Збереження координат траєкторії спуску у файл.\n'
                '5. Візуалізація результатів: побудова графіків рівнянь системи та траєкторії спуску алгоритму Хука-Дживса.')
    pdf.chapter_body(text_hid)

    pdf.chapter_title('Результати:')
    output = pdf.add_script_output('lab9.py', 'LAB_9')
    pdf.set_font('LiberationMono', '', 10)
    pdf.multi_cell(0, 5, output)
    pdf.ln()
    
    pdf.image('LAB_9/plot.png', w=140)
    pdf.ln()

    pdf.chapter_title('Висновки:')
    pdf.chapter_body('Програмно реалізовано алгоритм Хука-Дживса та показано його ефективність на прикладі складної функції Розенброка. Застосувавши метод найменших квадратів для зведення задачі розв\'язку системи нелінійних рівнянь до задачі оптимізації, я успішно знайшов точки перетину двох нелінійних кривих, візуально підтвердивши точність алгоритмічного пошуку та його збіжність до правильного розв\'язку.')
    
    pdf.output('LAB_9/Шутяк_ЧисельніМетоди_лаб9.pdf')

def create_lab10_report():
    pdf = PDF()
    pdf.set_auto_page_break(auto=True, margin=15)
    pdf.add_title_page('Методи Рунге-Кутта та Адамса розв’язання задачі Коші', '10')
    
    pdf.add_page()
    pdf.chapter_title('Мета:')
    pdf.chapter_body('Вивчити метод Рунге-Кутта четвертого порядку та багатокроковий метод прогнозу і корекції Адамса для чисельного розв’язку задачі Коші для звичайних диференційних рівнянь першого порядку.')
    
    pdf.chapter_title('Хід роботи:')
    text_hid = ('1. Аналітичний та чисельний розв\'язок ЗДР y\' = -y + x + 1 з початковою умовою y(0)=1.\n'
                '2. Програмна реалізація методу прогнозу і корекції Адамса 2-го порядку з оцінкою похибки та автоматичним вибором кроку.\n'
                '3. Побудова графіків точної та наближеної (за формулою) локальної похибки, а також графіка зміни розміру кроку.\n'
                '4. Програмна реалізація методу Рунге-Кутта 4-го порядку.\n'
                '5. Оцінка локальної похибки методу Рунге-Кутта за допомогою правила Рунге (порівняння результатів з кроком h та h/2) та реалізація автоматичного вибору кроку.\n'
                '6. Візуалізація графіків похибок та адаптивного кроку для РК4.')
    pdf.chapter_body(text_hid)

    pdf.chapter_title('Результати:')
    output = pdf.add_script_output('lab10.py', 'LAB_10')
    pdf.set_font('LiberationMono', '', 10)
    pdf.multi_cell(0, 5, output)
    pdf.ln()
    
    pdf.image('LAB_10/lab10_part1_adams.png', w=170)
    pdf.ln()
    pdf.image('LAB_10/lab10_part2_rk4.png', w=170)
    pdf.ln()

    pdf.chapter_title('Висновки:')
    pdf.chapter_body('Практично підтверджено ефективність правила Рунге для оцінки локальної похибки без знання точного аналітичного розв\'язку. До того ж, запровадження адаптивної зміни кроку довело свою доцільність — алгоритми успішно зменшували крок на ділянках із швидкою зміною функції і збільшували його там, де крива ставала більш пологою, оптимізуючи таким чином обчислювальне навантаження.')
    
    pdf.output('LAB_10/Шутяк_ЧисельніМетоди_лаб10.pdf')

if __name__ == '__main__':
    create_lab8_report()
    create_lab9_report()
    create_lab10_report()
    print("Всі PDF звіти згенеровано!")
