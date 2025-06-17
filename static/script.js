document.addEventListener('DOMContentLoaded', () => {
    // --- БЛОК 1: ВСТАНОВЛЕННЯ ЗНАЧЕНЬ ЗА ЗАМОВЧУВАННЯМ ---
    const defaultValues = {
        S_nom_kVA: 20000.0,
        U_HV_nom_kV: 15.75,
        U_LV_nom_kV: 6.3,
        ukz_percent: 7.5,
        P_nl_kW: 13.5,
        P_sc_kW: 127.0,
        i0_percent: 0.75
    };

    for (const key in defaultValues) {
        const input = document.getElementById(key);
        if (input) {
            input.value = defaultValues[key];
        }
    }

    // --- БЛОК 2: ОТРИМАННЯ ДОСТУПУ ДО ЕЛЕМЕНТІВ СТОРІНКИ ---
    const form = document.getElementById('calc-form');
    const resultsPlaceholder = document.getElementById('results-placeholder');
    const resultsOutput = document.getElementById('results-output');
    const errorContainer = document.getElementById('error-container');
    const downloadBtn = document.getElementById('download-pdf');

    // Глобальна змінна для зберігання результатів для PDF
    let calculationResults = {};

    // --- БЛОК 3: ОБРОБКА ВІДПРАВКИ ФОРМИ ---
    form.addEventListener('submit', async (e) => {
        e.preventDefault(); // Запобігаємо перезавантаженню сторінки
        
        resultsPlaceholder.textContent = 'Виконується розрахунок...';
        resultsOutput.classList.add('hidden');
        errorContainer.classList.add('hidden');

        // Збір даних
        const data = {};
        form.querySelectorAll('input').forEach(input => {
            data[input.id] = input.value;
        });

        // Відправляємо запит на сервер
        const response = await fetch('/calculate', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify(data),
        });

        // Отримуємо та обробляємо відповідь
        const results = await response.json();
        calculationResults = results; // Для PDF
        
        if (results.error) {
            errorContainer.innerHTML = `<p class="error-message">Помилка: ${results.error}</p>`;
            errorContainer.classList.remove('hidden');
            resultsPlaceholder.textContent = 'Виникла помилка під час розрахунків.';
        } else {
            resultsPlaceholder.textContent = 'Результати розрахунків з\'являться тут.';
            displayResults(results);
            resultsOutput.classList.remove('hidden');
        }
    });

    // --- БЛОК 4: ФУНКЦІЯ ВІДОБРАЖЕННЯ РЕЗУЛЬТАТІВ ---
    function displayResults(results) {
    const sectionsContainer = document.getElementById('results-sections');
    sectionsContainer.innerHTML = ''; // Очищуємо попередні результати

    // Словник з назвами секцій
    const sectionTitles = {
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
    };

    // Масив для правильного порядку секцій
    const sectionOrder = ['input_data'].concat(Array.from({length: 15}, (_, i) => `section_${i + 1}`));

    sectionOrder.forEach(key => {
        if (results[key] && Object.keys(results[key]).length > 0) {
            const sectionData = results[key];
            const title = sectionTitles[key] || key.replace(/_/g, ' ').toUpperCase();

            const details = document.createElement('details');
            if (key === 'input_data' || key === 'section_1' || key === 'section_12') {
                details.open = true;
            }
            
            const summary = document.createElement('summary');
            summary.textContent = title;
            details.appendChild(summary);

            const contentContainer = document.createElement('div');
            contentContainer.className = 'section-content';
            let ul = document.createElement('ul');

            // Сортування ключів всередині секції
            const sortedKeys = Object.keys(sectionData).sort((a, b) => {
                const getNum = (str) => {
                    if (str.includes('(Графік')) return 1000 + parseInt(str.replace(/\D/g, ''), 10);
                    const match = str.match(/\((\d+)\)/);
                    return match ? parseInt(match[1], 10) : Infinity;
                };
                return getNum(a) - getNum(b);
            });

            for (const itemKey of sortedKeys) {
                const itemValue = sectionData[itemKey];

                if (itemKey.includes('(Графік')) {
                    const plotDiv = document.createElement('div');
                    plotDiv.className = 'plot-container';
                    const plotTitle = document.createElement('h4');
                    plotTitle.textContent = itemKey;
                    const img = document.createElement('img');
                    img.src = 'data:image/png;base64,' + itemValue;
                    plotDiv.appendChild(plotTitle);
                    plotDiv.appendChild(img);
                    contentContainer.appendChild(plotDiv);
                } else if (itemValue && typeof itemValue === 'object') {
                    const listItem = document.createElement('li');
                    // Перевіряємо, чи є у об'єкта поле 'formula', щоб відобразити повний розрахунок
                    if (itemValue.formula && itemValue.calculation) {
                        listItem.innerHTML = `<strong>${itemKey}</strong>`;
                        const formulaBlock = document.createElement('div');
                        formulaBlock.className = 'formula-block';

                        const formulaGeneral = document.createElement('div');
                        katex.render(itemValue.formula, formulaGeneral, { throwOnError: false, displayMode: true });

                        const formulaCalc = document.createElement('div');
                        katex.render(itemValue.calculation, formulaCalc, { throwOnError: false, displayMode: true });
                        
                        const formulaResult = document.createElement('div');
                        formulaResult.className = 'formula-result';
                        formulaResult.innerHTML = `<strong>Результат: ${Number(itemValue.result).toFixed(4)}</strong>`;

                        formulaBlock.appendChild(formulaGeneral);
                        formulaBlock.appendChild(formulaCalc);
                        formulaBlock.appendChild(formulaResult);
                        listItem.appendChild(formulaBlock);
                    } else if (itemValue.result !== undefined) {
                         // Для вхідних даних та інших простих результатів
                         listItem.innerHTML = `<strong>${itemKey}:</strong> ${Number(itemValue.result).toFixed(4)}`;
                    }
                    if(listItem.innerHTML) {
                        ul.appendChild(listItem);
                    }
                }
            }
            
            if (ul.hasChildNodes()) {
                contentContainer.prepend(ul);
            }
            details.appendChild(contentContainer);
            sectionsContainer.appendChild(details);
        }
    });
}

    // --- БЛОК 5: PDF ---
    downloadBtn.addEventListener('click', async () => {
        const response = await fetch('/generate_pdf', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ results: calculationResults }),
        });

        if (response.ok) {
            const blob = await response.blob();
            const url = window.URL.createObjectURL(blob);
            const a = document.createElement('a');
            a.style.display = 'none';
            a.href = url;
            a.download = "transformer_full_calculation.pdf";
            document.body.appendChild(a);
            a.click();
            window.URL.revokeObjectURL(url);
            a.remove();
        } else {
            alert("Помилка при генерації PDF.");
        }
    });
});