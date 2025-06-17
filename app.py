from flask import Flask, render_template, request, jsonify, make_response
import io
import re
import base64
from fpdf import FPDF
from fpdf.enums import XPos, YPos
from calculations import calculate_transformer_design

app = Flask(__name__)

# СТВОРЕННЯ PDF

class PDF(FPDF):
    def header(self):
        self.add_font('DejaVu', '', 'DejaVuSansCondensed.ttf')
        self.set_font('DejaVu', '', 15)
        self.cell(0, 10, 'Повний розрахунок силового трансформатора', align='C', new_x=XPos.LMARGIN, new_y=YPos.NEXT)
        self.ln(10)

    def chapter_title(self, title):
        self.set_font('DejaVu', 'B', 12)
        self.multi_cell(0, 8, title, align='L', new_x=XPos.LMARGIN, new_y=YPos.NEXT)
        self.ln(4)

    # Допоміжна функція для сортування
    def _get_sort_key(self, key_string):
        # Шукаємо число в дужках, наприклад (1), (25), (115)
        match = re.search(r'\((\d+)\)', key_string)
        if match:
            return int(match.group(1))
        # Якщо це графік, даємо йому великий номер, щоб він був у кінці
        if '(Графік' in key_string:
            match_plot = re.search(r'\(Графік (\d+)\)', key_string)
            if match_plot:
                return 1000 + int(match_plot.group(1))
        # Якщо номера немає, ставимо в кінець
        return float('inf')

    def chapter_body(self, data):
        self.set_font('DejaVu', '', 10)
        col_width_key = self.w * 0.6
        col_width_val = self.w * 0.3

        # СОРТУЄМО КЛЮЧІ перед виводом
        sorted_keys = sorted(data.keys(), key=self._get_sort_key)

        for key in sorted_keys:
            value = data[key]
            # Пропускаємо графіки, їх обробить інший метод
            if '(Графік' in key:
                continue

            current_y = self.get_y()
            
            self.set_font('DejaVu', 'B', 9)
            self.multi_cell(col_width_key, 5, key, border=0)
            
            self.set_xy(self.l_margin + col_width_key, current_y)
            self.set_font('DejaVu', '', 9)

            if isinstance(value, dict) and 'result' in value:
                result_str = f"{value['result']:.4f}"
                self.multi_cell(col_width_val, 5, result_str, border=0, new_x=XPos.LMARGIN, new_y=YPos.NEXT)
                if value.get('calc_text'):
                    self.set_font('DejaVu', '', 8)
                    self.set_text_color(100, 100, 100)
                    self.multi_cell(0, 4, f"   └─ {value['calc_text']}", new_x=XPos.LMARGIN, new_y=YPos.NEXT)
                    self.set_text_color(0, 0, 0)
            else: # Для даних старого формату (якщо є)
                 self.multi_cell(col_width_val, 5, str(value), border=0, new_x=XPos.LMARGIN, new_y=YPos.NEXT)
            
            self.ln(1)

    def chapter_plots(self, data):
        self.set_font('DejaVu', 'B', 11)
        
        # СОРТУЄМО КЛЮЧІ ГРАФІКІВ перед виводом
        plot_keys = sorted([k for k in data.keys() if '(Графік' in k], key=self._get_sort_key)

        for key in plot_keys:
            value = data[key]
            if isinstance(value, str):
                try:
                    self.ln(10)
                    self.cell(0, 8, key, new_x=XPos.LMARGIN, new_y=YPos.NEXT, align='C')
                    self.ln(2)
                    
                    image_bytes = base64.b64decode(value)
                    img_file = io.BytesIO(image_bytes)
                    
                    img_x = (self.w - 180) / 2
                    self.image(img_file, x=img_x, w=180)
                    img_file.close()
                    
                except Exception:
                    self.set_text_color(255, 0, 0)
                    self.cell(0, 10, f"Не вдалося вставити графік: {key}")
                    self.set_text_color(0, 0, 0)
                    
    def footer(self):
        self.set_y(-15)
        self.set_font('DejaVu', '', 8)
        self.cell(0, 10, f'Сторінка {self.page_no()}', align='C')

@app.route('/generate_pdf', methods=['POST'])
def generate_pdf():
    data = request.get_json().get('results', {})
    if not data or "input_data" not in data:
        return "Invalid data for PDF generation", 400

    pdf = PDF('P', 'mm', 'A4')
    pdf.add_font('DejaVu', '', 'DejaVuSansCondensed.ttf')
    pdf.add_font('DejaVu', 'B', 'DejaVuSansCondensed-Bold.ttf')
    pdf.add_page()

    section_titles = {
        "input_data": "Вхідні дані",
        "section_1": "1. Попереднє визначення взаємопов’язаних параметрів трансформатора.",
        "section_2": "2.1 Розрахунок обмотки низької напруги.",
        "section_3": "2.2 Розрахунок обмотки високої напруги.",
        "section_4": "2.3 Розрахунок регулювальної обмотки.",
        "section_5": "2.4 Втрати в обмотках.",
        "section_6": "2.5 Розрахунок напруги короткого замикання.",
        "section_7": "2.6 Розрахунок механічної стійкості.",
        "section_8": "2.7 Розрахунок температурної стійкості.",
        "section_9": "2.8 Конструкція кістяка трансформатора.",
        "section_10": "2.9 Розрахунок втрат неробочого ходу.",
        "section_11": "2.10 Розрахунок струму неробочого ходу.",
        "section_12": "2.11 Побудова векторної діаграми.",
        "section_13": "2.12 Розрахунок баку.",
        "section_14": "2.13 Тепловий розрахунок.",
        "section_15": "2.14 Визначення маси."
    }
    
    section_order = ['input_data'] + [f'section_{i}' for i in range(1, 16)]

    for section_key in section_order:
        if section_key in data and data[section_key]:
            section_data = data[section_key]
            title = section_titles.get(section_key, "Розділ")

            if section_key != 'input_data':
                pdf.add_page()

            pdf.chapter_title(title)
            pdf.chapter_body(section_data) 
            pdf.chapter_plots(section_data)

    response = make_response(bytes(pdf.output()))
    response.headers['Content-Type'] = 'application/pdf'
    response.headers['Content-Disposition'] = 'attachment; filename=transformer_full_calculation.pdf'
    
    return response

@app.route('/')
def index():
    default_params = {
        'S_nom_kVA': 20000.0, 'U_HV_nom_kV': 15.75, 'U_LV_nom_kV': 6.3,
        'ukz_percent': 7.5, 'P_nl_kW': 13.5, 'P_sc_kW': 127.0, 'i0_percent': 0.75
    }
    return render_template('index.html', defaults=default_params)

@app.route('/calculate', methods=['POST'])
def calculate():
    data = request.get_json()
    params = {key: float(value) for key, value in data.items()}
    results = calculate_transformer_design(**params)
    return jsonify(results)

if __name__ == '__main__':
    app.run(debug=True)