import math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import io
import base64
from scipy.optimize import fsolve

def _generate_plot_base64(plot_function, *args):
    buf = io.BytesIO()
    plot_function(*args)
    plt.savefig(buf, format='png', bbox_inches='tight')
    plt.close()
    buf.seek(0)
    base64_string = base64.b64encode(buf.getvalue()).decode('utf-8')
    buf.close()
    return base64_string

# ==============================================================================
# ФУНКЦІЇ ДЛЯ РОЗРАХУНКУ ДАНИХ  ДЛЯ ГРАФІКІВ
# ==============================================================================
def calculate_G_u(x_val, A_val, k_c_val, k_ya_val):
    if x_val == 0: return 0
    return 0.492 * (10**4) * k_c_val * k_ya_val * (A_val**3) * (x_val**3)
def calculate_G_c(x_val, A1_val, A2_val):
    return (A1_val + A2_val) * x_val
def calculate_G_ya(x_val, B1_val, B2_val):
    if x_val == 0: return 0
    return B1_val * (x_val**3) + B2_val * (x_val**2)
def calculate_P_x(G_c_val, G_u_val, G_ya_val, k_pd_val, p_c_loss_val, k_pu_val, p_ya_loss_val):
    term_Gya_minus_6Gu = G_ya_val - 6 * G_u_val
    return k_pd_val * p_c_loss_val * G_c_val + k_pu_val * G_u_val + k_pd_val * p_ya_loss_val * term_Gya_minus_6Gu
def calculate_Q_x(G_c_val, G_u_val, G_ya_val, k_td_p, k_td_dp, q_c_r, k_tu_val, k_tpl_val, q_ya_r):
    term_Gya_minus_6Gu = G_ya_val - 6 * G_u_val
    return k_td_p * k_td_dp * q_c_r * G_c_val + k_tu_val * G_u_val * k_tpl_val + k_td_p * k_td_dp * q_ya_r * term_Gya_minus_6Gu
def calculate_i_nx_percent(Q_x_val, S_n_VA_val):
    if S_n_VA_val == 0: return 0
    return (100 * Q_x_val) / S_n_VA_val
def calculate_G_obm(x_val, C1_val):
    if x_val == 0: return float('inf')
    return C1_val / (x_val**2)
def calculate_J(P_k_W_val, G_obm_val, k_d_val):
    if G_obm_val == 0 or (k_d_val * P_k_W_val / G_obm_val) < 0: return 0
    return math.sqrt((k_d_val * P_k_W_val) / (2.4 * G_obm_val))
def calculate_sigma_p(x_val, M_val):
    return M_val * (x_val**3)
def calculate_d_core(x_val, A_val):
    return A_val * x_val
def calculate_d12(d_core_val, a_val):
    return a_val * d_core_val
def calculate_l_core_winding(d12_val, beta_val):
    if beta_val == 0: return float('inf')
    return (math.pi * d12_val) / beta_val
def calculate_L_ser(d12_val, A01_val, A02_val, b_val, d_core_val):
    return d12_val + A01_val + A02_val + b_val * d_core_val

# ==============================================================================
# Оголошення ВСІХ 11 функцій для побудови графіків
# ==============================================================================
def plot_1_losses(beta_range, Px_vals):
    plt.figure(figsize=(8, 5)); plt.plot(beta_range, Px_vals); plt.title('Графік 1: Втрати НХ ($P_x$) від $\\beta$'); plt.xlabel('$\\beta$'); plt.ylabel('$P_x$, Вт'); plt.grid(True)
def plot_2_current(beta_range, i_nx_vals):
    plt.figure(figsize=(8, 5)); plt.plot(beta_range, i_nx_vals); plt.title('Графік 2: Струм НХ ($i_{нх}$) від $\\beta$'); plt.xlabel('$\\beta$'); plt.ylabel('$i_{нх}$, %'); plt.grid(True)
def plot_3_density(beta_range, J_vals):
    plt.figure(figsize=(8, 5)); plt.plot(beta_range, J_vals); plt.title('Графік 3: Щільність струму ($J$) від $\\beta$'); plt.xlabel('$\\beta$'); plt.ylabel('$J$, А/мм$^2$'); plt.grid(True)
def plot_4_stress(beta_range, sigma_p_vals):
    plt.figure(figsize=(8, 5)); plt.plot(beta_range, sigma_p_vals); plt.title('Графік 4: Механічні напруження ($\\sigma_p$) від $\\beta$'); plt.xlabel('$\\beta$'); plt.ylabel('$\\sigma_p$, МПа'); plt.grid(True)
def plot_5_d_core(beta_range, d_core_vals):
    plt.figure(figsize=(8, 5)); plt.plot(beta_range, d_core_vals); plt.title('Графік 5: Діаметр стрижня ($d$) від $\\beta$'); plt.xlabel('$\\beta$'); plt.ylabel('$d$, м'); plt.grid(True)
def plot_6_l_core(beta_range, l_core_vals):
    plt.figure(figsize=(8, 5)); plt.plot(beta_range, l_core_vals); plt.title('Графік 6: Висота стрижня/обмотки ($l$) від $\\beta$'); plt.xlabel('$\\beta$'); plt.ylabel('$l$, м'); plt.grid(True)
def plot_7_L_ser(beta_range, L_ser_vals):
    plt.figure(figsize=(8, 5)); plt.plot(beta_range, L_ser_vals); plt.title('Графік 7: Відстань між осями стрижнів ($L_{сер}$) від $\\beta$'); plt.xlabel('$\\beta$'); plt.ylabel('$L_{сер}$, м'); plt.grid(True)
def plot_8_external_V(beta_vals, results_U_V):
    plt.figure(figsize=(8, 5)); [plt.plot(beta_vals, values, label=case) for case, values in results_U_V.items()]; plt.title('Графік 8: Зовнішні характеристики, B'); plt.xlabel('Коефіцієнт навантаження (β)'); plt.ylabel('Напруга на виході (U2), В'); plt.legend(); plt.grid(True); plt.gca().invert_yaxis()
def plot_9_external_percent(beta_vals, results_U_percent):
    plt.figure(figsize=(8, 5)); [plt.plot(beta_vals, values, label=case) for case, values in results_U_percent.items()]; plt.title('Графік 9: Зовнішні характеристики, %'); plt.xlabel('Коефіцієнт навантаження (β)'); plt.ylabel('Напруга на виході (U2), %'); plt.legend(); plt.grid(True)
def plot_10_delta_U_V(beta_vals, results_delta_U_V):
    plt.figure(figsize=(8, 5)); [plt.plot(beta_vals, values, label=case) for case, values in results_delta_U_V.items()]; plt.title('Графік 10: Падіння напруги, В'); plt.xlabel('Коефіцієнт навантаження (β)'); plt.ylabel('Падіння напруги (ΔU), В'); plt.legend(); plt.grid(True)
def plot_11_delta_U_percent(beta_vals, results_delta_U_percent):
    plt.figure(figsize=(8, 5)); [plt.plot(beta_vals, values, label=case) for case, values in results_delta_U_percent.items()]; plt.title('Графік 11: Падіння напруги, %'); plt.xlabel('Коефіцієнт навантаження (β)'); plt.ylabel('Падіння напруги (ΔU), %'); plt.legend(); plt.grid(True)

# ==============================================================================
# Головна функція розрахунку
# ==============================================================================
def calculate_transformer_design(S_nom_kVA, U_HV_nom_kV, U_LV_nom_kV, ukz_percent, P_nl_kW, P_sc_kW, i0_percent):
    try:
        results = {f'section_{i}': {} for i in range(1, 16)}
        results['input_data'] = {}
        results['formulas'] = {}

        results['input_data'] = {
            '(Вхід) Номінальна потужність, кВА': {'result': S_nom_kVA},
            '(Вхід) Номінальна напруга ВН, кВ': {'result': U_HV_nom_kV},
            '(Вхід) Номінальна напруга НН, кВ': {'result': U_LV_nom_kV},
            '(Вхід) Напруга КЗ, %': {'result': ukz_percent},
            '(Вхід) Втрати неробочого ходу, кВт': {'result': P_nl_kW},
            '(Вхід) Втрати КЗ, кВт': {'result': P_sc_kW},
            '(Вхід) Струм неробочого ходу, %': {'result': i0_percent}
        }
        

        S_n_VA = S_nom_kVA * 1000
        U_l_HV_V = U_HV_nom_kV * 1000
        U_l_LV_V = U_LV_nom_kV * 1000
        m_phases = 3
        
        # ================== РОЗДІЛ 1: ПОПЕРЕДНЄ ВИЗНАЧЕННЯ ==================
        sec1 = results['section_1']
        S_st_VA = S_n_VA / m_phases
        sec1['(1) Потужність фази (Sст), ВА'] = {
            'formula': r'S_{ст} = \frac{S_{п}}{m}',
            'calculation': fr'S_{{ст}} = \frac{{{S_n_VA}}}{{{m_phases}}}',
            'calc_text': f'S_st = {S_n_VA} / {m_phases}',
            'result': S_st_VA
        }

        # (2) Номінальний лінійний струм ВН
        I_l_HV_A = S_n_VA / (math.sqrt(3) * U_l_HV_V)
        sec1['(2) Номінальний лінійний струм ВН (Iлв), А'] = {
            'formula': r'I_{лв} = \frac{S_{п}}{\sqrt{3} \cdot U_{лв}}',
            'calculation': fr'I_{{лв}} = \frac{{{S_n_VA:.0f}}}{{\sqrt{{3}} \cdot {U_l_HV_V:.0f}}}',
            'calc_text': f'I_lv = {S_n_VA:.0f} / (sqrt(3) * {U_l_HV_V:.0f})',
            'result': I_l_HV_A
        }

        # (3) Номінальний лінійний струм НН
        I_l_LV_A = S_n_VA / (math.sqrt(3) * U_l_LV_V)
        sec1['(3) Номінальний лінійний струм НН (Iлн), А'] = {
            'formula': r'I_{лн} = \frac{S_{п}}{\sqrt{3} \cdot U_{лн}}',
            'calculation': fr'I_{{лн}} = \frac{{{S_n_VA:.0f}}}{{\sqrt{{3}} \cdot {U_l_LV_V:.0f}}}',
            'calc_text': f'I_ln = {S_n_VA:.0f} / (sqrt(3) * {U_l_LV_V:.0f})',
            'result': I_l_LV_A
        }
        
        # (4) Фазний струм ВН
        I_ph_HV_A = I_l_HV_A / math.sqrt(3)
        sec1['(4) Фазний струм ВН (Iфв) (для D), А'] = {
            'formula': r'I_{фв} = \frac{I_{лв}}{\sqrt{3}}',
            'calculation': fr'I_{{фв}} = \frac{{{I_l_HV_A:.2f}}}{{\sqrt{{3}}}}',
            'calc_text': f'{I_l_HV_A:.2f} / sqrt(3)',
            'result': I_ph_HV_A
        }

        # (5) Фазна напруга ВН
        U_ph_HV_V = U_l_HV_V
        sec1['(5) Фазна напруга ВН (Uфв) (для D), В'] = {
            'formula': r'U_{фв} = U_{лв}',
            'calculation': fr'U_{{фв}} = {U_l_HV_V:.2f}',
            'calc_text': f'{U_l_HV_V:.2f}',
            'result': U_ph_HV_V
        }

        # (6) Фазний струм НН
        I_ph_LV_A = I_l_LV_A
        sec1['(6) Фазний струм НН (Iфн) (для Y), А'] = {
            'formula': r'I_{фн} = I_{лн}',
            'calculation': fr'I_{{фн}} = {I_l_LV_A:.2f}',
            'calc_text': f'{I_l_LV_A:.2f}',
            'result': I_ph_LV_A
        }

        # (7) Фазна напруга НН
        U_ph_LV_V = U_l_LV_V / math.sqrt(3)
        sec1['(7) Фазна напруга НН (Uфн) (для Y), В'] = {
            'formula': r'U_{фн} = \frac{U_{лн}}{\sqrt{3}}',
            'calculation': fr'U_{{фн}} = \frac{{{U_l_LV_V:.2f}}}{{\sqrt{{3}}}}',
            'calc_text': f'{U_l_LV_V:.2f} / sqrt(3)',
            'result': U_ph_LV_V
        }

        A_00, A_01, A_02, A_10, A_20 = 0.022, 0.022, 0.023, 0.050, 0.080
        k_pya, k_zk, k_kz_const, k_d_const, B_c = 1.02, 0.97, 0.927, 0.82, 1.565
        k_formula11, k_p_formula16, f_hz = 0.506, 0.95, 50.0
        a_const_formula17, e_const_formula19, b_const_formula19 = 1.38, 0.41, 0.312
        K_B_const, k_ip_const, sigma_p_MPa_const = 1.9, 1.06, 60.0

        u_akz_percent_calc = (P_sc_kW / S_nom_kVA) * 100
        sec1['(8) Активна складова напруги КЗ (uакз), %'] = {
            'formula': r'u_{а\%} = \frac{P_{кз}}{S_{п}} \cdot 100\%',
            'calculation': fr'u_{{а\%}} = \frac{{{P_sc_kW}}}{{{S_nom_kVA}}} \cdot 100\%',
            'calc_text': f'u_a% = ({P_sc_kW} / {S_nom_kVA}) * 100',
            'result': u_akz_percent_calc
        }

        # (9) Реактивна складова напруги КЗ (uркз), %
        u_rkz_percent_calc = math.sqrt(ukz_percent**2 - u_akz_percent_calc**2)
        sec1['(9) Реактивна складова напруги КЗ (uркз), %'] = {
            'formula': r'u_{р\%} = \sqrt{u_{кз\%}^2 - u_{а\%}^2}',
            'calculation': fr'u_{{р\%}} = \sqrt{{{ukz_percent:.2f}^2 - {u_akz_percent_calc:.2f}^2}}', # <-- Виправлено фігурні дужки
            'calc_text': f'u_r% = sqrt({ukz_percent:.2f}^2 - {u_akz_percent_calc:.2f}^2)', # <-- Додана кома
            'result': u_rkz_percent_calc
        }

        # (10) Коефіцієнт kку
        k_ku_calc = math.sqrt(2) * (100 / ukz_percent) * (1 + math.exp(-math.pi * u_akz_percent_calc / u_rkz_percent_calc))
        sec1['(10) Коефіцієнт kку'] = {
            'formula': r'k_{ку} = \sqrt{2} \cdot \frac{100}{u_{кз\%}} \cdot \left(1 + e^{-\frac{\pi \cdot u_{а\%}}{u_{р\%}}}\right)',
            'calculation': fr'k_{{ку}} = \sqrt{{2}} \cdot \frac{{100}}{{{ukz_percent:.2f}}} \cdot \left(1 + e^{{-\frac{{\pi \cdot {u_akz_percent_calc:.2f}}}{{{u_rkz_percent_calc:.2f}}}}}\right)',
            'calc_text': f'k_ku = sqrt(2) * (100 / {ukz_percent:.2f}) * (1 + exp(-pi * {u_akz_percent_calc:.2f} / {u_rkz_percent_calc:.2f}))', # <-- Додана кома
            'result': k_ku_calc
        }
        S_st_kVA = S_st_VA / 1000
        # (11) Приведена ширина двох обмоток (ac), м
        a_c_m_calc = k_formula11 * ((S_st_VA * (10**-3)) **0.25) * (10**-2)
        sec1['(11) Приведена ширина двох обмоток (ac), м'] = {
            'formula': r'a_c = k_{11} \cdot \left(S_{ст} \cdot 10^{-3}\right)^{0.25} \cdot 10^{-2}',
            'calculation': fr'a_c = {k_formula11} \cdot \left({S_st_kVA} \cdot 10^{{-3}}\right)^{{0.25}} \cdot 10^{{-2}}',
            'calc_text': f'a_c = {k_formula11} * (({S_st_kVA} * 10**-3)**0.25) * 10**-2',
            'result': a_c_m_calc
        }

        # (12) Сума радіальних розмірів обмоток (ap), м
        a_p_m_calc = A_01 + a_c_m_calc
        sec1['(12) Сума радіальних розмірів обмоток (ap), м'] = {
            'formula': r'a_p = A_{01} + a_c',
            'calculation': fr'a_p = {A_01} + {a_c_m_calc:.4f}',
            'calc_text': f'a_p = {A_01} + {a_c_m_calc:.4f}',
            'result': a_p_m_calc
        }

        # (13) Коефіцієнт заповнення сталлю (kc)
        k_c_coeff_calc = k_zk * k_kz_const
        sec1['(13) Коефіцієнт заповнення сталлю (kc)'] = {
            'formula': r'k_c = k_{зк} \cdot k_{кз}',
            'calculation': fr'k_c = {k_zk} \cdot {k_kz_const}',
            'calc_text': f'k_c = {k_zk} * {k_kz_const}',
            'result': k_c_coeff_calc
        }

        # (14) Індукція у ярмі (Bя), Тл
        B_ya_T_calc = B_c / k_pya
        sec1['(14) Індукція у ярмі (Bя), Тл'] = {
            'formula': r'B_я = \frac{B_c}{k_{ря}}',
            'calculation': fr'B_я = \frac{{{B_c}}}{{{k_pya}}}',
            'calc_text': f'B_ya = {B_c} / {k_pya}',
            'result': B_ya_T_calc
        }

        # (15) Індукція в зазорі на косому стику (Bкс), Тл
        B_ks_T_calc = B_c / math.sqrt(2)
        sec1['(15) Індукція в зазорі на косому стику (Bкс), Тл'] = {
            'formula': r'B_{кс} = \frac{B_c}{\sqrt{2}}',
            'calculation': fr'B_{{кс}} = \frac{{{B_c}}}{{\sqrt{{2}}}}',
            'calc_text': f'B_ks = {B_c} / sqrt(2)',
            'result': B_ks_T_calc
        }

        # (16) Коефіцієнт A
        A_coeff_calc = 0.507 * (((S_st_VA * a_p_m_calc * k_p_formula16) / (f_hz * u_rkz_percent_calc * (B_c**2) * (k_c_coeff_calc**2) * (10**3)))**0.25)
        sec1['(16) Коефіцієнт A'] = {
            'formula': r'A = 0.507 \cdot \left( \frac{S_{ст} \cdot a_p \cdot k_p}{f \cdot u_{р\%} \cdot B_c^2 \cdot k_c^2 \cdot 10^3} \right)^{0.25}',
            'calculation': fr'A = 0.507 \cdot \left( \frac{{{S_st_VA} \cdot {a_p_m_calc:.4f} \cdot {k_p_formula16}}}{{{f_hz} \cdot {u_rkz_percent_calc:.2f} \cdot {B_c}^2 \cdot {k_c_coeff_calc:.2f}^2 \cdot 10^3}} \right)^{{0.25}}',
            'calc_text': f'A = 0.507 * (({S_st_VA} * {a_p_m_calc:.4f} * {k_p_formula16}) / ({f_hz} * {u_rkz_percent_calc:.2f} * ({B_c}**2) * ({k_c_coeff_calc:.2f}**2) * 10**3))**0.25',
            'result': A_coeff_calc
        }

        # (17) Коефіцієнт A1, кг
        A1_kg_calc = 5.663 * (10**4) * k_c_coeff_calc * (A_coeff_calc**3) * a_const_formula17
        sec1['(17) Коефіцієнт A1, кг'] = {
            'formula': r'A_1 = 5.663 \cdot 10^4 \cdot k_c \cdot A^3 \cdot a',
            'calculation': fr'A_1 = 5.663 \cdot 10^4 \cdot {k_c_coeff_calc:.2f} \cdot {A_coeff_calc:.4f}^3 \cdot {a_const_formula17}',
            'calc_text': f'A1 = 5.663e4 * {k_c_coeff_calc:.2f} * ({A_coeff_calc:.4f}**3) * {a_const_formula17}',
            'result': A1_kg_calc
        }

        # (18) Коефіцієнт A2, кг
        A2_kg_calc = 3.605 * (10**4) * k_c_coeff_calc * (A_coeff_calc**2) * ((A_10 + A_20) / 2)
        sec1['(18) Коефіцієнт A2, кг'] = {
            'formula': r'A_2 = 3.605 \cdot 10^4 \cdot k_c \cdot A^2 \cdot \frac{A_{10} + A_{20}}{2}',
            'calculation': fr'A_2 = 3.605 \cdot 10^4 \cdot {k_c_coeff_calc:.2f} \cdot {A_coeff_calc:.4f}^2 \cdot \frac{{{A_10} + {A_20}}}{{2}}',
            'calc_text': f'A2 = 3.605e4 * {k_c_coeff_calc:.2f} * ({A_coeff_calc:.4f}**2) * (({A_10} + {A_20}) / 2)',
            'result': A2_kg_calc
        }

        # (19) Коефіцієнт B1, кг
        B1_kg_calc = 2.4 * (10**4) * k_c_coeff_calc * k_pya * (A_coeff_calc**3) * (a_const_formula17 + b_const_formula19 + e_const_formula19)
        sec1['(19) Коефіцієнт B1, кг'] = {
            'formula': r'B_1 = 2.4 \cdot 10^4 \cdot k_c \cdot k_{ря} \cdot A^3 \cdot (a + b + e)',
            'calculation': fr'B_1 = 2.4 \cdot 10^4 \cdot {k_c_coeff_calc:.2f} \cdot {k_pya} \cdot {A_coeff_calc:.4f}^3 \cdot ({a_const_formula17} + {b_const_formula19} + {e_const_formula19})',
            'calc_text': f'B1 = 2.4e4 * {k_c_coeff_calc:.2f} * {k_pya} * ({A_coeff_calc:.4f}**3) * ({a_const_formula17} + {b_const_formula19} + {e_const_formula19})',
            'result': B1_kg_calc
        }

        # (20) Коефіцієнт B2, кг
        B2_kg_calc = 2.4 * (10**4) * k_c_coeff_calc * k_pya * (A_coeff_calc**2) * (A_01 + A_02)
        sec1['(20) Коефіцієнт B2, кг'] = {
            'formula': r'B_2 = 2.4 \cdot 10^4 \cdot k_c \cdot k_{ря} \cdot A^2 \cdot (A_{01} + A_{02})',
            'calculation': fr'B_2 = 2.4 \cdot 10^4 \cdot {k_c_coeff_calc:.2f} \cdot {k_pya} \cdot {A_coeff_calc:.4f}^2 \cdot ({A_01} + {A_02})',
            'calc_text': f'B2 = 2.4e4 * {k_c_coeff_calc:.2f} * {k_pya} * ({A_coeff_calc:.4f}**2) * ({A_01} + {A_02})',
            'result': B2_kg_calc
        }

        # (21) Коефіцієнт C1, кг
        C1_kg_calc = 2.46 * (10**-2) * ((S_n_VA * (a_const_formula17**2)) / (k_d_const * (k_c_coeff_calc**2) * (B_c**2) * u_akz_percent_calc * (A_coeff_calc**2) * (10**3)))
        sec1['(21) Коефіцієнт C1, кг'] = {
            'formula': r'C_1 = 2.46 \cdot 10^{-2} \cdot \frac{S_{п} \cdot a^2}{k_д \cdot k_c^2 \cdot B_c^2 \cdot u_{а\%} \cdot A^2 \cdot 10^3}',
            'calculation': fr'C_1 = 2.46 \cdot 10^{{-2}} \cdot \frac{{{S_n_VA} \cdot {a_const_formula17}^2}}{{{k_d_const} \cdot {k_c_coeff_calc:.2f}^2 \cdot {B_c}^2 \cdot {u_akz_percent_calc:.2f} \cdot {A_coeff_calc:.4f}^2 \cdot 10^3}}',
            'calc_text': f'C1 = 2.46e-2 * (({S_n_VA} * {a_const_formula17}**2) / ({k_d_const} * ({k_c_coeff_calc:.2f}**2) * ({B_c}**2) * {u_akz_percent_calc:.2f} * ({A_coeff_calc:.4f}**2) * 10**3))',
            'result': C1_kg_calc
        }

        # (22) Коефіцієнт M, МПа
        M_MPa_calc = 0.244 * (10**-6) * (k_ku_calc**2) * k_d_const * k_p_formula16 * ((P_sc_kW * (10**3)) / (a_const_formula17 * A_coeff_calc))
        sec1['(22) Коефіцієнт M, МПа'] = {
            'formula': r'M = 0.244 \cdot 10^{-6} \cdot k_{ку}^2 \cdot k_д \cdot k_p \cdot \frac{P_{кз} \cdot 10^3}{a \cdot A}',
            'calculation': fr'M = 0.244 \cdot 10^{{-6}} \cdot {k_ku_calc:.2f}^2 \cdot {k_d_const} \cdot {k_p_formula16} \cdot \frac{{{P_sc_kW} \cdot 10^3}}{{{a_const_formula17} \cdot {A_coeff_calc:.4f}}}',
            'calc_text': f'M = 0.244e-6 * ({k_ku_calc:.2f}**2) * {k_d_const} * {k_p_formula16} * (({P_sc_kW} * 10**3) / ({a_const_formula17} * {A_coeff_calc:.4f}))',
            'result': M_MPa_calc
        }

        # (23) Коефіцієнт B (для рівняння x)
        B_for_eq_calc = (2 * (A2_kg_calc + B2_kg_calc)) / (3 * B1_kg_calc)
        sec1['(23) Коефіцієнт B (для рівняння x)'] = {
            'formula': r'B_x = \frac{2 \cdot (A_2 + B_2)}{3 \cdot B_1}',
            'calculation': fr'B_x = \frac{{2 \cdot ({A2_kg_calc:.2f} + {B2_kg_calc:.2f})}}{{3 \cdot {B1_kg_calc:.2f}}}',
            'calc_text': f'B_x = (2 * ({A2_kg_calc:.2f} + {B2_kg_calc:.2f})) / (3 * {B1_kg_calc:.2f})',
            'result': B_for_eq_calc
        }

        # (24) Коефіцієнт C (для рівняння x)
        C_for_eq_calc = A1_kg_calc / (3 * B1_kg_calc)
        sec1['(24) Коефіцієнт C (для рівняння x)'] = {
            'formula': r'C_x = \frac{A_1}{3 \cdot B_1}',
            'calculation': fr'C_x = \frac{{{A1_kg_calc:.2f}}}{{3 \cdot {B1_kg_calc:.2f}}}',
            'calc_text': f'C_x = {A1_kg_calc:.2f} / (3 * {B1_kg_calc:.2f})',
            'result': C_for_eq_calc
        }

        # (25) Коефіцієнт D (для рівняння x)
        D_for_eq_calc = K_B_const * k_ip_const * (2/3) * (C1_kg_calc / B1_kg_calc)
        sec1['(25) Коефіцієнт D (для рівняння x)'] = {
            'formula': r'D_x = K_B \cdot k_{іп} \cdot \frac{2}{3} \cdot \frac{C_1}{B_1}',
            'calculation': fr'D_x = {K_B_const} \cdot {k_ip_const} \cdot \frac{{2}}{{3}} \cdot \frac{{{C1_kg_calc:.2f}}}{{{B1_kg_calc:.2f}}}',
            'calc_text': f'D_x = {K_B_const} * {k_ip_const} * (2/3) * ({C1_kg_calc:.2f} / {B1_kg_calc:.2f})',
            'result': D_for_eq_calc
        }

        # (27)
        def equation_to_solve(x, B, C, D):
            return x**5 + B * x**4 - C * x - D
        x_solution = fsolve(equation_to_solve, 1.0, args=(B_for_eq_calc, C_for_eq_calc, D_for_eq_calc))[0]
        sec1["(27) Розв'язок рівняння для x"] = {
            'formula': r'x^5 + B \cdot x^4 - C \cdot x - D = 0',
            'calculation': fr'\text{{Розв\'язок рівняння з B={B_for_eq_calc:.3f}, C={C_for_eq_calc:.3f}, D={D_for_eq_calc:.3f}}}',
            'calc_text': f'Solving x^5 + {B_for_eq_calc:.3f}*x^4 - {C_for_eq_calc:.3f}*x - {D_for_eq_calc:.3f} = 0 using fsolve',
            'result': x_solution
        }

        # (28) Значення бета
        beta_val_calc = x_solution**4
        sec1['(28) Значення бета'] = {
            'formula': r'\beta = x^4',
            'calculation': fr'\beta = {x_solution:.4f}^4',
            'calc_text': f'beta = {x_solution:.4f}**4',
            'result': beta_val_calc
        }

        # (29) Перевірка xj (<=)
        limit_xj_rhs_calc = 4.5 * (10**6) * math.sqrt((2.4 * (10**-12) * C1_kg_calc) / (k_d_const * (P_sc_kW * (10**3))))
        sec1['(29) Перевірка xj (<=)'] = {
            'formula': r'x_j \le 4.5 \cdot 10^6 \sqrt{\frac{2.4 \cdot 10^{-12} \cdot C_1}{k_д \cdot P_{кз} \cdot 10^3}}',
            'calculation': fr'x_j \le 4.5 \cdot 10^6 \sqrt{{\frac{{2.4 \cdot 10^{{-12}} \cdot {C1_kg_calc:.2f}}}{{{k_d_const} \cdot {P_sc_kW} \cdot 10^3}}}}',
            'calc_text': f'xj_limit <= 4.5e6 * sqrt((2.4e-12 * {C1_kg_calc:.2f}) / ({k_d_const} * ({P_sc_kW} * 1e3)))',
            'result': limit_xj_rhs_calc
        }

        # (30) Перевірка x_sigma (<=)
        limit_x_sigma_rhs_calc = (sigma_p_MPa_const / M_MPa_calc)**(1/3)
        sec1['(30) Перевірка x_sigma (<=)'] = {
            'formula': r'x_{\sigma} \le \sqrt[3]{\frac{\sigma_{p.доп}}{M}}',
            'calculation': fr'x_{{\sigma}} \le \sqrt[3]{{\frac{{{sigma_p_MPa_const}}}{{{M_MPa_calc:.4f}}}}}',
            'calc_text': f'x_sigma_limit <= ({sigma_p_MPa_const} / {M_MPa_calc:.4f})**(1/3)',
            'result': limit_x_sigma_rhs_calc
        }
        
        #константи
        k_pd = 1.2
        p_c_loss =  0.845
        p_ya_loss = 0.805
        k_pu = 9.74
        k_td_prime = 1.13
        k_td_double_prime = 1.15
        q_c_reac = 1.975
        q_koc = 2200
        q_sz = 29700 
        k_tu = 40.0
        k_tpl = 1.26
        q_ya_r = 1.887
        n_z_prime = 6
        n_z_double_prime = 0
        k_ya = 1
        Pi_c1 = 0.14512




        # (31) Маса кутів
        G_u_final = 0.492 * (10**4) * k_c_coeff_calc * 1.0 * (A_coeff_calc**3) * (x_solution**3)
        sec1['(31) Маса Gу, кг'] = {
            'formula': r'G_у = 0.492 \cdot 10^4 \cdot k_c \cdot A^3 \cdot x^3',
            'calculation': fr'G_у = 0.492 \cdot 10^4 \cdot {k_c_coeff_calc:.3f} \cdot {A_coeff_calc:.4f}^3 \cdot {x_solution:.4f}^3',
            'calc_text': f'G_u = 0.492e4 * {k_c_coeff_calc:.3f} * {A_coeff_calc:.4f}**3 * {x_solution:.4f}**3',
            'result': G_u_final
        }

        # (32) Маса стержнів
        G_c_final = (A1_kg_calc + A2_kg_calc) * x_solution
        sec1['(32) Маса Gc, кг'] = {
            'formula': r'G_c = (A_1 + A_2) \cdot x',
            'calculation': fr'G_c = ({A1_kg_calc:.2f} + {A2_kg_calc:.2f}) \cdot {x_solution:.4f}',
            'calc_text': f'G_c = ({A1_kg_calc:.2f} + {A2_kg_calc:.2f}) * {x_solution:.4f}',
            'result': G_c_final
        }

        # (33) Маса ярма
        G_ya_final = B1_kg_calc * (x_solution**3) + B2_kg_calc * (x_solution**2)
        sec1['(33) Маса ярма Gя, кг'] = {
            'formula': r'G_я = B_1 x^3 + B_2 x^2',
            'calculation': fr'G_я = {B1_kg_calc:.2f} \cdot {x_solution:.4f}^3 + {B2_kg_calc:.2f} \cdot {x_solution:.4f}^2',
            'calc_text': f'G_ya = {B1_kg_calc:.2f} * {x_solution:.4f}**3 + {B2_kg_calc:.2f} * {x_solution:.4f}**2',
            'result': G_ya_final
        }   

       # (34) Площа перерізу стержня Пс, м^2
        P_c = 0.785 * k_c_coeff_calc * (A_coeff_calc**2) * (x_solution**2)
        sec1['(34) Площа перерізу стержня Пс, м^2'] = {
            'formula': r'\Pi_c = \frac{\pi}{4} \cdot k_c \cdot A^2 \cdot x^2',
            'calculation': fr'\Pi_c = \frac{{\pi}}{{4}} \cdot {k_c_coeff_calc:.3f} \cdot {A_coeff_calc:.4f}^2 \cdot {x_solution:.4f}^2',
            'calc_text': f'P_c = (pi/4) * {k_c_coeff_calc:.3f} * {A_coeff_calc:.4f}**2 * {x_solution:.4f}**2',
            'result': P_c
        }

        # (35) Площа косого стику Пкс, м^2
        P_kc = P_c * math.sqrt(2)
        sec1['(35) Площа косого стику Пкс, м^2'] = {
            'formula': r'\Pi_{кс} = \Pi_c \cdot \sqrt{2}',
            'calculation': fr'\Pi_{{кс}} = {P_c:.4f} \cdot \sqrt{{2}}',
            'calc_text': f'P_kc = {P_c:.4f} * sqrt(2)',
            'result': P_kc
        }

        # (36) Втрати холостого ходу Px, Вт
        P_x = k_pd * p_c_loss * (G_c_final + 0.5 * k_pu * G_u_final) * k_pd * p_ya_loss * (G_ya_final - 6 * G_u_final + 0.5 * k_pu * G_u_final)
        sec1['(36) Втрати холостого ходу Px, Вт'] = {
            'formula': r'P_x = [k_{рд} p_{c.в} (G_c + 0.5 k_{ру} G_у)] \cdot [k_{рд} p_{я.в} (G_я - 6 G_у + 0.5 k_{ру} G_у)]',
            'calculation': fr'P_x = [{k_pd} \cdot {p_c_loss} ({G_c_final:.1f} + 0.5 \cdot {k_pu} \cdot {G_u_final:.1f})] \cdot [{k_pd} \cdot {p_ya_loss} ({G_ya_final:.1f} - 6 \cdot {G_u_final:.1f} + 0.5 \cdot {k_pu} \cdot {G_u_final:.1f})]',
            'calc_text': f'P_x = {k_pd}*{p_c_loss}*({G_c_final:.1f}+0.5*{k_pu}*{G_u_final:.1f}) * {k_pd}*{p_ya_loss}*({G_ya_final:.1f}-6*{G_u_final:.1f}+0.5*{k_pu}*{G_u_final:.1f})',
            'result': P_x
        }

        # (37) Реактивна потужність Qx, ВАр
        Q_x = k_td_prime * k_td_double_prime * q_c_reac * (G_c_final + 0.5 * k_tu * G_u_final * k_tpl)\
            + k_td_prime * k_td_double_prime * q_ya_r * (G_ya_final - 6 * G_u_final + 0.5 * k_tu * k_tpl * G_u_final)\
            + k_td_double_prime * q_koc * n_z_prime * 0 + q_sz * k_td_double_prime * n_z_double_prime * 0
        sec1['(46) Реактивна потужність Qx, ВАр'] = {
            'formula': r'Q_x = Q_{ст} + Q_{яр} + Q_{кос} + Q_{сз}',
            'calculation': fr'Q_x = ({k_td_prime} \cdot {k_td_double_prime} \cdot {q_c_reac} \cdot ({G_c_final:.1f} + ...)) + ({k_td_prime} \cdot {k_td_double_prime} \cdot {q_ya_r} \cdot ({G_ya_final:.1f} + ...)) + 0 + 0',
            'calc_text': f'Q_x = {k_td_prime}*{k_td_double_prime}*{q_c_reac}*({G_c_final:.1f}+0.5*{k_tu}*{G_u_final:.1f}*{k_tpl}) + {k_td_prime}*{k_td_double_prime}*{q_ya_r}*({G_ya_final:.1f}-6*{G_u_final:.1f}+0.5*{k_tu}*{k_tpl}*{G_u_final:.1f}) + ...',
            'result': Q_x
        }

        # (38) Загальна маса сталі Gст, кг
        G_st_final = G_c_final + G_ya_final
        sec1['(38) Загальна маса сталі Gст, кг'] = {
            'formula': r'G_{ст} = G_c + G_я',
            'calculation': fr'G_{{ст}} = {G_c_final:.2f} + {G_ya_final:.2f}',
            'calc_text': f'G_st = {G_c_final:.2f} + {G_ya_final:.2f}',
            'result': G_st_final
        }

        # (40) Струм холостого ходу i_хх, %
        i_nx_curr = ((10**3) * Q_x) / (10 * S_n_VA)
        sec1['(40) Струм холостого ходу i_хх, %'] = {
            'formula': r'i_{х\%} = \frac{10^3 \cdot Q_x}{10 \cdot S_{п}}',
            'calculation': fr'i_{{х\%}} = \frac{{10^3 \cdot {Q_x:.2f}}}{{10 \cdot {S_n_VA}}}',
            'calc_text': f'i_xx % = (1000 * {Q_x:.2f}) / (10 * {S_n_VA})',
            'result': i_nx_curr
        }

        # (41) Маса обмоток Gобм, кг
        G_obm_final = C1_kg_calc / (x_solution**2)
        sec1['(39) Маса обмоток Gобм, кг'] = {
            'formula': r'G_{обм} = \frac{C_1}{x^2}',
            'calculation': fr'G_{{обм}} = \frac{{{C1_kg_calc:.2f}}}{{{x_solution:.4f}^2}}',
            'calc_text': f'G_obm = {C1_kg_calc:.2f} / {x_solution:.4f}**2',
            'result': G_obm_final
        }


        # (42) Маса деталей Gд, кг
        G_d_final = 1.1 * 1.03 * G_obm_final
        sec1['(42) Маса деталей Gд, кг'] = {
            'formula': r'G_д = 1.1 \cdot 1.03 \cdot G_{обм}',
            'calculation': fr'G_д = 1.1 \cdot 1.03 \cdot {G_obm_final:.2f}',
            'calc_text': f'G_d = 1.1 * 1.03 * {G_obm_final:.2f}',
            'result': G_d_final
        }       
 
 
        # (43) Загальна маса активних матеріалів Gач, кг
        G_ach_final = K_B_const * G_d_final + G_st_final
        sec1['(43) Загальна маса активних матеріалів Gач, кг'] = {
            'formula': r'G_{ач} = K_B \cdot G_д + G_{ст}',
            'calculation': fr'G_{{ач}} = {K_B_const} \cdot {G_d_final:.2f} + {G_st_final:.2f}',
            'calc_text': f'G_ach = {K_B_const} * {G_d_final:.2f} + {G_st_final:.2f}',
            'result': G_ach_final
        }

        # (44) Густина струму J
        J_calc = (math.sqrt(k_d_const * (P_sc_kW * (10**3)) / (2.4 * G_obm_final))) * (10**-6) # Changed to 10**-6 for А/м^2 -> А/мм^2
        sec1['(44) Густина струму J, А/мм^2'] = {
            'formula': r'J = \sqrt{\frac{k_д \cdot P_{кз} \cdot 10^3}{2.4 \cdot G_{обм}}} \cdot 10^{-6}',
            'calculation': fr'J = \sqrt{{\frac{{{k_d_const} \cdot {P_sc_kW} \cdot 10^3}}{{2.4 \cdot {G_obm_final:.2f}}}}} \cdot 10^{{-6}}',
            'calc_text': f'J = sqrt(({k_d_const} * ({P_sc_kW}*1e3)) / (2.4 * {G_obm_final:.2f})) * 1e-6',
            'result': J_calc
        }

        # (45) Механічне напруження σp
        sigma_p_calc = M_MPa_calc * (x_solution ** 3)
        sec1['(45) Механічне напруження σp, МПа'] = {
            'formula': r'\sigma_p = M \cdot x^3',
            'calculation': fr'\sigma_p = {M_MPa_calc:.4f} \cdot {x_solution:.4f}^3',
            'calc_text': f'sigma_p = {M_MPa_calc:.4f} * {x_solution:.4f}**3',
            'result': sigma_p_calc
        }

        # (46) Діаметр стержня d
        d_calc = A_coeff_calc * x_solution
        sec1['(46) Діаметр стержня d, м'] = {
            'formula': r'd = A \cdot x',
            'calculation': fr'd = {A_coeff_calc:.4f} \cdot {x_solution:.4f}',
            'calc_text': f'd = {A_coeff_calc:.4f} * {x_solution:.4f}',
            'result': d_calc
        }

        # (47) Середній діаметр обмоток d12
        d_12_calc = a_const_formula17 * d_calc
        sec1['(47) Середній діаметр обмоток d12, м'] = {
            'formula': r'd_{12} = a \cdot d',
            'calculation': fr'd_{{12}} = {a_const_formula17} \cdot {d_calc:.4f}',
            'calc_text': f'd_12 = {a_const_formula17} * {d_calc:.4f}',
            'result': d_12_calc
        }

        # (48) Висота обмотки l
        l_calc = math.pi * d_12_calc / beta_val_calc
        sec1['(48) Висота обмотки l, м'] = {
            'formula': r'l = \frac{\pi \cdot d_{12}}{\beta}',
            'calculation': fr'l = \frac{{\pi \cdot {d_12_calc:.4f}}}{{{beta_val_calc:.4f}}}',
            'calc_text': f'l = (pi * {d_12_calc:.4f}) / {beta_val_calc:.4f}',
            'result': l_calc
        }

        # (49) Довжина середньої лінії магнітної індукції
        L_ser_calc = d_12_calc + A_01 + A_02 + b_const_formula19 * d_calc
        sec1['(49) Довжина середньої лінії магнітної індукції Lсер, м'] = {
            'formula': r'L_{сер} = d_{12} + A_{01} + A_{02} + b \cdot d',
            'calculation': fr'L_{{сер}} = {d_12_calc:.4f} + {A_01} + {A_02} + {b_const_formula19} \cdot {d_calc:.4f}',
            'calc_text': f'L_ser = {d_12_calc:.4f} + {A_01} + {A_02} + {b_const_formula19} * {d_calc:.4f}',
            'result': L_ser_calc
        }

        # (50) Перерахований діаметр стержня (d), м
        beta_formula_50 = 1.1
        d_recalc_calc = A_coeff_calc * (beta_formula_50 ** (1/4))
        sec1['(50) Перерахований діаметр стержня (d), м'] = {
            'formula': r'd = A \cdot \sqrt[4]{\beta}',
            'calculation': fr'd = {A_coeff_calc:.4f} \cdot \sqrt[4]{{{beta_formula_50}}}',
            'calc_text': f'd = {A_coeff_calc:.4f} * {beta_formula_50}**(1/4)',
            'result': d_recalc_calc
        }

        # Прийняте значення d для подальших розрахунків
        d_recalc = 0.450
        sec1['(50.a) Прийнятий діаметр стержня (d), м'] = {
            'formula': r'd_{прийн}',
            'calculation': r'\text{Прийняте конструктивне значення}',
            'calc_text': 'Adopted design value',
            'result': d_recalc
        }

        # (51) Перераховане значення бета (β)
        beta_val_recalc = (d_recalc / A_coeff_calc)**4
        sec1['(51) Перераховане значення бета (β)'] = {
            'formula': r'\beta = \left(\frac{d_{прийн}}{A}\right)^4',
            'calculation': fr'\beta = \left(\frac{{{d_recalc}}}{{{A_coeff_calc:.4f}}}\right)^4',
            'calc_text': f'beta = ({d_recalc} / {A_coeff_calc:.4f})**4',
            'result': beta_val_recalc
        }

        # (52) Площа перерізу стержня
        Pi_s_calc = k_zk * Pi_c1
        sec1['(52) Площа перерізу стержня, м^2'] = {
            'formula': r'\Pi_s = k_{зк} \cdot \Pi_{c1}',
            'calculation': fr'\Pi_s = {k_zk} \cdot {Pi_c1:.4f}',
            'calc_text': f'Pi_s = {k_zk} * {Pi_c1:.4f}',
            'result': Pi_s_calc
        }

        # (53) Середній діаметр обмотки
        d2_obm_calc = a_const_formula17 * d_recalc
        sec1['(53) Середній діаметр обмотки d2, м'] = {
            'formula': r'd_2 = a \cdot d_{прийн}',
            'calculation': fr'd_2 = {a_const_formula17} \cdot {d_recalc}',
            'calc_text': f'd2_obm = {a_const_formula17} * {d_recalc}',
            'result': d2_obm_calc
        }

        # (54) Довжина витка обмотки
        l_obm_calc = (math.pi * d2_obm_calc) / beta_val_recalc
        sec1['(54) Довжина витка обмотки l_обм, м'] = {
            'formula': r'l_{обм} = \frac{\pi \cdot d_2}{\beta}',
            'calculation': fr'l_{{обм}} = \frac{{\pi \cdot {d2_obm_calc:.3f}}}{{{beta_val_recalc:.3f}}}',
            'calc_text': f'l_obm = (pi * {d2_obm_calc:.3f}) / {beta_val_recalc:.3f}',
            'result': l_obm_calc
        }

        # ================== РОЗДІЛ 2.1: ОБМОТКА НН ==================

        sec2 = results['section_2']

        # (55) Початкова ЕРС витка (eв), В
        e_v_55 = 4.44 * f_hz * Pi_s_calc * B_c
        sec2['(55) Початкова ЕРС витка (eв), В'] = {
            'formula': r'e_в = 4.44 \cdot f \cdot \Pi_s \cdot B_c',
            'calculation': fr'e_в = 4.44 \cdot {f_hz} \cdot {Pi_s_calc:.4f} \cdot {B_c}',
            'calc_text': f'e_v = 4.44 * {f_hz} * {Pi_s_calc:.4f} * {B_c}',
            'result': e_v_55
        }

        # (56) Попередня кількість витків НН
        N_nn_prelim = U_ph_LV_V / e_v_55
        sec2['(56) Попередня кількість витків НН'] = {
            'formula': r'N_{НН.поп} = \frac{U_{ф.НН}}{e_в}',
            'calculation': fr'N_{{НН.поп}} = \frac{{{U_ph_LV_V}}}{{{e_v_55:.3f}}}',
            'calc_text': f'N_nn_prelim = {U_ph_LV_V} / {e_v_55:.3f}',
            'result': N_nn_prelim
        }

        # Прийнята кількість витків НН
        N_nn_adopted = 74
        sec2['(доп.) Прийнята кількість витків НН'] = {
            'formula': r'N_{НН.прийн}',
            'calculation': r'\text{Прийняте значення}',
            'calc_text': 'Adopted value',
            'result': N_nn_adopted
        }

        # (57) Уточнена ЕРС витка, В
        e_v_57 = U_ph_LV_V / N_nn_adopted
        sec2['(57) Уточнена ЕРС витка, В'] = {
            'formula': r'e_{в.ут} = \frac{U_{ф.НН}}{N_{НН.прийн}}',
            'calculation': fr'e_{{в.ут}} = \frac{{{U_ph_LV_V}}}{{{N_nn_adopted}}}',
            'calc_text': f'e_v_ уточн. = {U_ph_LV_V} / {N_nn_adopted}',
            'result': e_v_57
        }

        # (58) Середня густина струму (Jcp), МА/м^2
        J_cp_58 = 0.746 * k_d_const * (((P_sc_kW * (10**3)) * e_v_57) / (S_n_VA * d2_obm_calc)) * 10
        sec2['(58) Середня густина струму (Jcp), МА/м^2'] = {
            'formula': r'J_{сер} = 0.746 \cdot k_д \cdot \frac{P_{кз} \cdot 10^3 \cdot e_{в.ут}}{S_{ном} \cdot d_2} \cdot 10',
            'calculation': fr'J_{{сер}} = 0.746 \cdot {k_d_const} \cdot \frac{{{P_sc_kW} \cdot 10^3 \cdot {e_v_57:.3f}}}{{{S_n_VA} \cdot {d2_obm_calc:.3f}}} \cdot 10',
            'calc_text': f'J_cp = 0.746 * {k_d_const} * (({P_sc_kW}e3 * {e_v_57:.3f}) / ({S_n_VA} * {d2_obm_calc:.3f})) * 10',
            'result': J_cp_58
        }

        # (59) Попередній переріз витка НН, мм^2
        S_nn_prelim_mm2 = I_ph_LV_A / J_cp_58
        sec2['(59) Попередній переріз витка НН, мм^2'] = {
            'formula': r'S_{НН.поп} = \frac{I_{ф.НН}}{J_{сер}}',
            'calculation': fr'S_{{НН.поп}} = \frac{{{I_ph_LV_A:.2f}}}{{{J_cp_58:.3f}}}',
            'calc_text': f'S_nn_prelim = {I_ph_LV_A:.2f} / {J_cp_58:.3f}',
            'result': S_nn_prelim_mm2
        }

        # (60) Перевірка макс. розміру проводу (b_NN_max), м
        b_NN_max_m = ((2000 * 1) / (1.07 * (J_cp_58 * 1e6)**2 * 1e-8))
        sec2['(60) Перевірка макс. розміру проводу (b_NN_max), м'] = {
            'formula': r'b_{НН.макс} \le \frac{2000}{1.07 \cdot (J_{сер} \cdot 10^6)^2 \cdot 10^{-8}}',
            'calculation': fr'b_{{НН.макс}} \le \frac{{2000}}{{1.07 \cdot ({J_cp_58:.3f} \cdot 10^6)^2 \cdot 10^{{-8}}}}',
            'calc_text': f'b_NN_max <= 2000 / (1.07 * ({J_cp_58:.3f}e6)**2 * 1e-8)',
            'result': b_NN_max_m
        }

        # (61) Фактична площа витка НН, мм^2
        n1, S_nn_prime_mm2 = 20, 29.3 # Прийняті значення
        S_nn_actual_mm2 = S_nn_prime_mm2 * n1
        sec2['(61) Фактична площа витка НН, мм^2'] = {
            'formula': r'S_{НН.факт} = S^{\prime}_{НН} \cdot n_1',
            'calculation': fr'S_{{НН.факт}} = {S_nn_prime_mm2} \cdot {n1}',
            'calc_text': f'S_nn_actual = {S_nn_prime_mm2} * {n1}',
            'result': S_nn_actual_mm2
        }

        # (62) Фактична густина струму (J_NN), МА/м^2
        S_nn_actual_m2 = S_nn_actual_mm2 * 1e-6
        J_NN_A_m2 = I_ph_LV_A / S_nn_actual_m2
        sec2['(62) Фактична густина струму (J_NN), МА/м^2'] = {
            'formula': r'J_{НН} = \frac{I_{ф.НН}}{S_{НН.факт}}',
            'calculation': fr'J_{{НН}} = \frac{{{I_ph_LV_A:.2f}}}{{{S_nn_actual_m2}}}',
            'calc_text': f'J_NN = {I_ph_LV_A:.2f} / {S_nn_actual_m2}',
            'result': J_NN_A_m2 / 1e6
        }

        # (63) Осьовий розмір витка (hNN), м
        h_NN_m = l_obm_calc / (N_nn_adopted + 1) - 0.004
        sec2['(63) Осьовий розмір витка (hNN), м'] = {
            'formula': r'h_{НН} = \frac{l_{обм}}{N_{НН} + 1} - 0.004',
            'calculation': fr'h_{{НН}} = \frac{{{l_obm_calc:.3f}}}{{{N_nn_adopted} + 1}} - 0.004',
            'calc_text': f'h_NN = {l_obm_calc:.3f} / ({N_nn_adopted} + 1) - 0.004',
            'result': h_NN_m
        }

        # (64) Радіальний розмір обмотки (ANN), м
        a1_prime_insulated_mm = 5.00 # Прийняте значення
        a1_prime_m = a1_prime_insulated_mm / 1000
        A_NN_m = a1_prime_m * (n1 / 2 - 1)
        sec2['(64) Радіальний розмір обмотки (ANN), м'] = {
            'formula': r'A_{НН} = a^{\prime}_1 \cdot (\frac{n_1}{2} - 1)',
            'calculation': fr'A_{{НН}} = {a1_prime_m} \cdot (\frac{{{n1}}}{{2}} - 1)',
            'calc_text': f'A_NN = {a1_prime_m} * ({n1}/2 - 1)',
            'result': A_NN_m
        }

        # (65) Горизонтальний ізоляційний канал (h1), м
        sec2['(65) Горизонтальний ізоляційний канал (h1), м'] = {
            'formula': r'h_1',
            'calculation': r'\text{Прийняте значення}',
            'calc_text': 'Adopted value',
            'result': 0.005
        }

        # (66) Осьовий розмір обмотки (lNN), м
        b1_prime_insulated_mm = 7.20 # Прийняте значення
        b1_prime_m = b1_prime_insulated_mm / 1000
        l_NN_m = 2 * b1_prime_m * (N_nn_adopted + 1) + 0.96 * 0.005 * (2 * N_nn_adopted + 1)
        sec2['(66) Осьовий розмір обмотки (lNN), м'] = {
            'formula': r'l_{НН} = 2b^{\prime}_1(N_{НН}+1) + 0.96 \cdot h_1 \cdot (2N_{НН}+1)',
            'calculation': fr'l_{{НН}} = 2 \cdot {b1_prime_m} ({N_nn_adopted}+1) + 0.96 \cdot 0.005 \cdot (2 \cdot {N_nn_adopted}+1)',
            'calc_text': f'l_NN = 2*{b1_prime_m}*({N_nn_adopted}+1) + 0.96*0.005*(2*{N_nn_adopted}+1)',
            'result': l_NN_m
        }

        # (67) Внутрішній діаметр обмотки (dNN), м
        d_NN_m = d_recalc + 2 * A_00
        sec2['(67) Внутрішній діаметр обмотки (dNN), м'] = {
            'formula': r'd_{НН} = d_{прийн} + 2 \cdot A_{00}',
            'calculation': fr'd_{{НН}} = {d_recalc} + 2 \cdot {A_00}',
            'calc_text': f'd_NN = {d_recalc} + 2 * {A_00}',
            'result': d_NN_m
        }

        # (68) Зовнішній діаметр обмотки (DNN), м
        D_NN_m = d_NN_m + 2 * A_NN_m
        sec2['(68) Зовнішній діаметр обмотки (DNN), м'] = {
            'formula': r'D_{НН} = d_{НН} + 2 \cdot A_{НН}',
            'calculation': fr'D_{{НН}} = {d_NN_m:.3f} + 2 \cdot {A_NN_m:.3f}',
            'calc_text': f'D_NN = {d_NN_m:.3f} + 2 * {A_NN_m:.3f}',
            'result': D_NN_m
        }

        # (69) Маса металу обмотки НН (M01), кг
        M01_kg = (28e3 * 3) * (d_NN_m + D_NN_m) * 0.5 * N_nn_adopted * S_nn_actual_m2
        sec2['(69) Маса металу обмотки НН (M01), кг'] = {
            'formula': r'M_{01} = (m \cdot \rho) \cdot \frac{d_{НН} + D_{НН}}{2} \cdot \pi \cdot N_{НН} \cdot S_{НН.факт}',
            'calculation': fr'M_{{01}} = (28 \cdot 10^3 \cdot 3) \cdot \frac{{{d_NN_m:.3f} + {D_NN_m:.3f}}}{{2}} \cdot {N_nn_adopted} \cdot {S_nn_actual_m2}',
            'calc_text': f'M01 = (28e3*3) * ({d_NN_m:.3f} + {D_NN_m:.3f})*0.5 * {N_nn_adopted} * {S_nn_actual_m2}',
            'result': M01_kg
        }

        # (70) Маса обмотки НН з ізоляцією (M'01), кг
        M_prime_01_kg = M01_kg * 1.02
        sec2['(70) Маса обмотки НН з ізоляцією (M\'01), кг'] = {
            'formula': r'M^{\prime}_{01} = M_{01} \cdot 1.02',
            'calculation': fr'M^{{\prime}}_{{01}} = {M01_kg:.2f} \cdot 1.02',
            'calc_text': f'M_prime_01 = {M01_kg:.2f} * 1.02',
            'result': M_prime_01_kg
        }

        # (71) Коефіцієнт (beta_дод1)
        b1_bare_mm = 6.70 # Прийняте значення
        b1_bare_m = b1_bare_mm / 1000
        beta_dod1 = (b1_bare_m * N_nn_adopted * 2 * k_p_formula16) / l_NN_m
        sec2['(71) Коефіцієнт (beta_дод1)'] = {
            'formula': r'\beta_{дод1} = \frac{2 \cdot b_{1.гол} \cdot N_{НН} \cdot k_p}{l_{НН}}',
            'calculation': fr'\beta_{{дод1}} = \frac{{2 \cdot {b1_bare_m} \cdot {N_nn_adopted} \cdot {k_p_formula16}}}{{{l_NN_m:.3f}}}',
            'calc_text': f'beta_dod1 = (2 * {b1_bare_m} * {N_nn_adopted} * {k_p_formula16}) / {l_NN_m:.3f}',
            'result': beta_dod1
        }

        # (72) Коефіцієнт додаткових втрат (кдод1)
        a1_bare_mm = 4.50 # Прийняте значення
        a1_bare_m = a1_bare_mm / 1000
        k_dod1 = 1 + 0.095 * 1e8 * beta_dod1**2 * a1_bare_m**4 * (n1 / 2)**2
        sec2['(72) Коефіцієнт додаткових втрат (кдод1)'] = {
            'formula': r'k_{дод1} = 1 + 0.095 \cdot 10^8 \beta_{дод1}^2 a_{1.гол}^4 (\frac{n_1}{2})^2',
            'calculation': fr'k_{{дод1}} = 1 + 0.095 \cdot 10^8 \cdot {beta_dod1:.3f}^2 \cdot {a1_bare_m}^4 \cdot (\frac{{{n1}}}{{2}})^2',
            'calc_text': f'k_dod1 = 1 + 0.095e8 * {beta_dod1:.3f}**2 * {a1_bare_m}**4 * ({n1}/2)**2',
            'result': k_dod1
        }

        # (73) Щільність теплового потоку (QNN), Вт/м^2
        Q_NN_W_m2 = (107e-10 * J_NN_A_m2 * I_ph_LV_A * k_dod1) / ((b1_prime_m + A_NN_m) * 0.75)
        sec2['(73) Щільність теплового потоку (QNN), Вт/м^2'] = {
            'formula': r'Q_{НН} = \frac{107 \cdot 10^{-10} \cdot J_{НН} \cdot I_{ф.НН} \cdot k_{дод1}}{0.75 \cdot (b^{\prime}_{1} + A_{НН})}',
            'calculation': fr'Q_{{НН}} = \frac{{107 \cdot 10^{{-10}} \cdot {J_NN_A_m2:.1f} \cdot {I_ph_LV_A:.2f} \cdot {k_dod1:.3f}}}{{0.75 \cdot ({b1_prime_m} + {A_NN_m:.3f})}}',
            'calc_text': f'Q_NN = (107e-10 * {J_NN_A_m2:.1f} * {I_ph_LV_A:.2f} * {k_dod1:.3f}) / (({b1_prime_m} + {A_NN_m:.3f}) * 0.75)',
            'result': Q_NN_W_m2
        }
        # ================== РОЗДІЛ 2.2: ОБМОТКА ВН ==================
        sec3 = results['section_3']
        # (74) Напруга на один рівень регулювання (delta_u), В
        tap_step_percent = 1.25
        delta_u_V = U_l_HV_V * (tap_step_percent / 100)
        sec3['(74) Напруга на один рівень регулювання (delta_u), В'] = {
            'formula': r'\Delta U = U_{л.ВН} \cdot \frac{\text{крок}_{\%}}{100}',
            'calculation': fr'\Delta U = {U_l_HV_V} \cdot \frac{{{tap_step_percent}}}{{100}}',
            'calc_text': f'delta_u = {U_l_HV_V} * ({tap_step_percent} / 100)',
            'result': delta_u_V
        }

        # (75) Попередня кількість витків на рівень
        N_p_prelim = delta_u_V / (math.sqrt(3) * e_v_57)
        sec3['(75) Попередня кількість витків на рівень'] = {
            'formula': r'N_{р.поп} = \frac{\Delta U}{\sqrt{3} \cdot e_{в.ут}}',
            'calculation': fr'N_{{р.поп}} = \frac{{{delta_u_V:.2f}}}{{\sqrt{{3}} \cdot {e_v_57:.3f}}}',
            'calc_text': f'N_p_prelim = {delta_u_V:.2f} / (sqrt(3) * {e_v_57:.3f})',
            'result': N_p_prelim
        }

        # Прийнята кількість витків на рівень
        N_p_adopted = 2
        sec3['(доп.) Прийнята кількість витків на рівень'] = {
            'formula': r'N_{р.прийн}',
            'calculation': r'\text{Прийняте значення}',
            'calc_text': 'Adopted value',
            'result': N_p_adopted
        }

        # (76) Уточнена напруга на рівень, В
        delta_u_recalc_V = e_v_57 * N_p_adopted * math.sqrt(3)
        sec3['(76) Уточнена напруга на рівень, В'] = {
            'formula': r'\Delta U_{ут} = e_{в.ут} \cdot N_{р.прийн} \cdot \sqrt{3}',
            'calculation': fr'\Delta U_{{ут}} = {e_v_57:.3f} \cdot {N_p_adopted} \cdot \sqrt{{3}}',
            'calc_text': f'delta_u_recalc = {e_v_57:.3f} * {N_p_adopted} * sqrt(3)',
            'result': delta_u_recalc_V
        }

        # (77) Загальна кількість регулювальних витків
        N_p_prime_logical = N_p_adopted * 8
        sec3['(77) Загальна кількість регулювальних витків'] = {
            'formula': r'N^{\prime}_{р} = N_{р.прийн} \cdot 8',
            'calculation': fr'N^{{\prime}}_{{р}} = {N_p_adopted} \cdot 8',
            'calc_text': f'N_p_prime = {N_p_adopted} * 8',
            'result': N_p_prime_logical
        }

        # (78) Попередня кількість витків ВН
        N_vn_prelim = U_ph_HV_V / e_v_57
        sec3['(78) Попередня кількість витків ВН'] = {
            'formula': r'N_{ВН.поп} = \frac{U_{ф.ВН}}{e_{в.ут}}',
            'calculation': fr'N_{{ВН.поп}} = \frac{{{U_ph_HV_V}}}{{{e_v_57:.3f}}}',
            'calc_text': f'N_vn_prelim = {U_ph_HV_V} / {e_v_57:.3f}',
            'result': N_vn_prelim
        }

        # Прийнята кількість витків ВН
        N_vn_adopted = 320
        sec3['(доп.) Прийнята кількість витків ВН'] = {
            'formula': r'N_{ВН.прийн}',
            'calculation': r'\text{Прийняте значення}',
            'calc_text': 'Adopted value',
            'result': N_vn_adopted
        }

        # (79) Попередня густина струму ВН, МА/м^2
        J_VN_prelim_MA_m2 = 2 * J_cp_58 - (J_NN_A_m2 / 1e6)
        sec3['(79) Попередня густина струму ВН, МА/м^2'] = {
            'formula': r'J_{ВН.поп} = 2 \cdot J_{сер} - J_{НН}',
            'calculation': fr'J_{{ВН.поп}} = 2 \cdot {J_cp_58:.3f} - {J_NN_A_m2/1e6:.3f}',
            'calc_text': f'J_VN_prelim = 2 * {J_cp_58:.3f} - {J_NN_A_m2/1e6:.3f}',
            'result': J_VN_prelim_MA_m2
        }

        # (80) Попередній переріз витка ВН, мм^2
        S_vn_prelim_mm2 = I_ph_HV_A / J_VN_prelim_MA_m2
        sec3['(80) Попередній переріз витка ВН, мм^2'] = {
            'formula': r'S_{ВН.поп} = \frac{I_{ф.ВН}}{J_{ВН.поп}}',
            'calculation': fr'S_{{ВН.поп}} = \frac{{{I_ph_HV_A:.2f}}}{{{J_VN_prelim_MA_m2:.3f}}}',
            'calc_text': f'S_vn_prelim = {I_ph_HV_A:.2f} / {J_VN_prelim_MA_m2:.3f}',
            'result': S_vn_prelim_mm2
        }

        # (81) Фактична площа витка ВН, мм^2
        n2, S_vn_prime_mm2 = 4, 36.6 # Прийняті значення
        S_vn_actual_mm2 = S_vn_prime_mm2 * n2
        sec3['(81) Фактична площа витка ВН, мм^2'] = {
            'formula': r'S_{ВН.факт} = S^{\prime}_{ВН} \cdot n_2',
            'calculation': fr'S_{{ВН.факт}} = {S_vn_prime_mm2} \cdot {n2}',
            'calc_text': f'S_vn_actual = {S_vn_prime_mm2} * {n2}',
            'result': S_vn_actual_mm2
        }

        # (82) Фактична густина струму ВН, МА/м^2
        S_vn_actual_m2 = S_vn_actual_mm2 * 1e-6
        J_VN_A_m2 = I_ph_HV_A / S_vn_actual_m2
        sec3['(82) Фактична густина струму ВН, МА/м^2'] = {
            'formula': r'J_{ВН} = \frac{I_{ф.ВН}}{S_{ВН.факт}}',
            'calculation': fr'J_{{ВН}} = \frac{{{I_ph_HV_A:.2f}}}{{{S_vn_actual_m2}}}',
            'calc_text': f'J_VN = {I_ph_HV_A:.2f} / {S_vn_actual_m2}',
            'result': J_VN_A_m2 / 1e6
        }

        # (83) Попередня кількість котушок
        b2_prime_insulated_mm = 9.00 # Прийняте значення
        b2_prime_m = b2_prime_insulated_mm / 1000
        K2_prelim = l_NN_m / (b2_prime_m + 0.005)
        sec3['(83) Попередня кількість котушок'] = {
            'formula': r'K_{2.поп} = \frac{l_{НН}}{b^{\prime}_2 + 0.005}',
            'calculation': fr'K_{{2.поп}} = \frac{{{l_NN_m:.3f}}}{{{b2_prime_m} + 0.005}}',
            'calc_text': f'K2_prelim = {l_NN_m:.3f} / ({b2_prime_m} + 0.005)',
            'result': K2_prelim
        }

        # Прийнята кількість котушок
        K2_adopted = 130
        sec3['(доп.) Прийнята кількість котушок'] = {
            'formula': r'K_{2.прийн}',
            'calculation': r'\text{Прийняте значення}',
            'calc_text': 'Adopted value',
            'result': K2_adopted
        }

        # (84) Витків в одній котушці
        N_vn_k = N_vn_adopted / K2_adopted
        sec3['(84) Витків в одній котушці'] = {
            'formula': r'N_{ВН/к} = \frac{N_{ВН.прийн}}{K_{2.прийн}}',
            'calculation': fr'N_{{ВН/к}} = \frac{{{N_vn_adopted}}}{{{K2_adopted}}}',
            'calc_text': f'N_vn_k = {N_vn_adopted} / {K2_adopted}',
            'result': N_vn_k
        }

        # Округлене значення витків в котушці
        N_vn_k_rounded = 3
        sec3['(доп.) Округлене число витків в котушці'] = {
            'formula': r'N_{ВН/к.окр}',
            'calculation': r'\text{Округлене значення}',
            'calc_text': 'Rounded value',
            'result': N_vn_k_rounded
        }

        # (85) Радіальний розмір обмотки ВН, м
        a2_prime_insulated_mm = 5.50 # Прийняте значення
        a2_prime_m = a2_prime_insulated_mm / 1000
        A_VN_m = N_vn_k_rounded * a2_prime_m * n2
        sec3['(85) Радіальний розмір обмотки ВН, м'] = {
            'formula': r'A_{ВН} = N_{ВН/к.окр} \cdot a^{\prime}_2 \cdot n_2',
            'calculation': fr'A_{{ВН}} = {N_vn_k_rounded} \cdot {a2_prime_m} \cdot {n2}',
            'calc_text': f'A_VN = {N_vn_k_rounded} * {a2_prime_m} * {n2}',
            'result': A_VN_m
        }

        # (86) Висота обмотки ВН, м
        l_VN_m = l_NN_m
        sec3['(86) Висота обмотки ВН, м'] = {
            'formula': r'l_{ВН} = l_{НН}',
            'calculation': fr'l_{{ВН}} = {l_NN_m:.3f}',
            'calc_text': f'l_VN = {l_NN_m:.3f}',
            'result': l_VN_m
        }

        # (87) Внутрішній діаметр обмотки ВН, м
        d_VN_m = D_NN_m + 2 * A_01
        sec3['(87) Внутрішній діаметр обмотки ВН, м'] = {
            'formula': r'd_{ВН} = D_{НН} + 2 \cdot A_{01}',
            'calculation': fr'd_{{ВН}} = {D_NN_m:.3f} + 2 \cdot {A_01}',
            'calc_text': f'd_VN = {D_NN_m:.3f} + 2 * {A_01}',
            'result': d_VN_m
        }

        # (88) Зовнішній діаметр обмотки ВН, м
        D_VN_m = d_VN_m + 2 * A_VN_m
        sec3['(88) Зовнішній діаметр обмотки ВН, м'] = {
            'formula': r'D_{ВН} = d_{ВН} + 2 \cdot A_{ВН}',
            'calculation': fr'D_{{ВН}} = {d_VN_m:.3f} + 2 \cdot {A_VN_m:.3f}',
            'calc_text': f'D_VN = {d_VN_m:.3f} + 2 * {A_VN_m:.3f}',
            'result': D_VN_m
        }

        # (89) Маса металу обмотки ВН, кг
        M02_kg = (28e3 * 3) * (d_VN_m + D_VN_m) * 0.5 * N_vn_adopted * S_vn_actual_m2
        sec3['(89) Маса металу обмотки ВН, кг'] = {
            'formula': r'M_{02} = (m \cdot \rho) \cdot \frac{d_{ВН} + D_{ВН}}{2} \cdot \pi \cdot N_{ВН} \cdot S_{ВН.факт}',
            'calculation': fr'M_{{02}} = (28 \cdot 10^3 \cdot 3) \cdot \frac{{{d_VN_m:.3f} + {D_VN_m:.3f}}}{{2}} \cdot {N_vn_adopted} \cdot {S_vn_actual_m2}',
            'calc_text': f'M02 = (28e3*3) * ({d_VN_m:.3f} + {D_VN_m:.3f})*0.5 * {N_vn_adopted} * {S_vn_actual_m2}',
            'result': M02_kg
        }

        # (90) Маса обмотки ВН з ізоляцією (M'02), кг
        M_prime_02_kg = M02_kg * 1.015
        sec3['(90) Маса обмотки ВН з ізоляцією (M\'02), кг'] = {
            'formula': r'M^{\prime}_{02} = M_{02} \cdot 1.015',
            'calculation': fr'M^{{\prime}}_{{02}} = {M02_kg:.2f} \cdot 1.015',
            'calc_text': f'M_prime_02 = {M02_kg:.2f} * 1.015',
            'result': M_prime_02_kg
        }

        # (91) Коефіцієнт (beta_дод2)
        b2_bare_mm = 8.50 # Прийняте значення
        b2_bare_m = b2_bare_mm / 1000
        beta_dod2 = (b2_bare_m * K2_adopted * k_p_formula16) / l_VN_m
        sec3['(91) Коефіцієнт (beta_дод2)'] = {
            'formula': r'\beta_{дод2} = \frac{b_{2.гол} \cdot K_2 \cdot k_p}{l_{ВН}}',
            'calculation': fr'\beta_{{дод2}} = \frac{{{b2_bare_m} \cdot {K2_adopted} \cdot {k_p_formula16}}}{{{l_VN_m:.3f}}}',
            'calc_text': f'beta_dod2 = ({b2_bare_m} * {K2_adopted} * {k_p_formula16}) / {l_VN_m:.3f}',
            'result': beta_dod2
        }

        # (92) Коефіцієнт додаткових втрат (кдод2)
        a2_bare_mm = 5.00 # Прийняте значення
        a2_bare_m = a2_bare_mm / 1000
        k_dod2 = 1 + 0.095 * 1e8 * beta_dod2**2 * a2_bare_m**4 * N_vn_k_rounded**2
        sec3['(92) Коефіцієнт додаткових втрат (кдод2)'] = {
            'formula': r'k_{дод2} = 1 + 0.095 \cdot 10^8 \beta_{дод2}^2 a_{2.гол}^4 N_{ВН/к}^2',
            'calculation': fr'k_{{дод2}} = 1 + 0.095 \cdot 10^8 \cdot {beta_dod2:.3f}^2 \cdot {a2_bare_m}^4 \cdot {N_vn_k_rounded}^2',
            'calc_text': f'k_dod2 = 1 + 0.095e8 * {beta_dod2:.3f}**2 * {a2_bare_m}**4 * {N_vn_k_rounded}**2',
            'result': k_dod2
        }

        # (93) Щільність теплового потоку ВН, Вт/м^2
        Q_VN_W_m2 = (107e-10 * J_VN_A_m2 * I_ph_HV_A * k_dod2) / ((b2_prime_m + A_VN_m) * 0.75)
        sec3['(93) Щільність теплового потоку ВН, Вт/м^2'] = {
            'formula': r'Q_{ВН} = \frac{107 \cdot 10^{-10} \cdot J_{ВН} \cdot I_{ф.ВН} \cdot k_{дод2}}{0.75 \cdot (b^{\prime}_{2} + A_{ВН})}',
            'calculation': fr'Q_{{ВН}} = \frac{{107 \cdot 10^{{-10}} \cdot {J_VN_A_m2:.1f} \cdot {I_ph_HV_A:.2f} \cdot {k_dod2:.3f}}}{{0.75 \cdot ({b2_prime_m} + {A_VN_m:.3f})}}',
            'calc_text': f'Q_VN = (107e-10 * {J_VN_A_m2:.1f} * {I_ph_HV_A:.2f} * {k_dod2:.3f}) / (({b2_prime_m} + {A_VN_m:.3f}) * 0.75)',
            'result': Q_VN_W_m2
        }
        # ================== РОЗДІЛ 2.3: РЕГУЛЮВАЛЬНА ОБМОТКА ==================
        sec4 = results['section_4']
        # Прийняті значення для РО
        a2_prime_ro_mm, b2_prime_ro_mm, n3_passes = 6.35, 9.85, 4

        # (94) Радіальний розмір обмотки РО, м
        a2_prime_ro_m = a2_prime_ro_mm / 1000
        A_RO_m = a2_prime_ro_m
        sec4['(94) Радіальний розмір обмотки РО, м'] = {
            'formula': r'A_{РО} = a^{\prime}_{2.РО}',
            'calculation': fr'A_{{РО}} = {a2_prime_ro_m}',
            'calc_text': f'A_RO = {a2_prime_ro_m}',
            'result': A_RO_m
        }

        # (95) Висота обмотки РО, м
        b2_prime_ro_m = b2_prime_ro_mm / 1000
        l_RO_m = n3_passes * b2_prime_ro_m * (N_p_prime_logical + 1)
        sec4['(95) Висота обмотки РО, м'] = {
            'formula': r'l_{РО} = n_3 \cdot b^{\prime}_{2.РО} \cdot (N^{\prime}_р + 1)',
            'calculation': fr'l_{{РО}} = {n3_passes} \cdot {b2_prime_ro_m} \cdot ({N_p_prime_logical} + 1)',
            'calc_text': f'l_RO = {n3_passes} * {b2_prime_ro_m} * ({N_p_prime_logical} + 1)',
            'result': l_RO_m
        }

        # (96) Внутрішній діаметр обмотки РО, м
        d_RO_m = D_VN_m + 2 * A_02
        sec4['(96) Внутрішній діаметр обмотки РО, м'] = {
            'formula': r'd_{РО} = D_{ВН} + 2 \cdot A_{02}',
            'calculation': fr'd_{{РО}} = {D_VN_m:.3f} + 2 \cdot {A_02}',
            'calc_text': f'd_RO = {D_VN_m:.3f} + 2 * {A_02}',
            'result': d_RO_m
        }

        # (97) Зовнішній діаметр обмотки РО, м
        D_RO_m = d_RO_m + 2 * A_RO_m
        sec4['(97) Зовнішній діаметр обмотки РО, м'] = {
            'formula': r'D_{РО} = d_{РО} + 2 \cdot A_{РО}',
            'calculation': fr'D_{{РО}} = {d_RO_m:.3f} + 2 \cdot {A_RO_m}',
            'calc_text': f'D_RO = {d_RO_m:.3f} + 2 * {A_RO_m}',
            'result': D_RO_m
        }

        # (98) Маса металу обмотки РО (M03), кг
        M03_kg = (28e3 * 3) * (d_RO_m + D_RO_m) * 0.5 * N_p_prime_logical * S_vn_actual_m2
        sec4['(98) Маса металу обмотки РО (M03), кг'] = {
            'formula': r'M_{03} = (m \cdot \rho) \cdot \frac{d_{РО} + D_{РО}}{2} \cdot \pi \cdot N^{\prime}_р \cdot S_{ВН.факт}',
            'calculation': fr'M_{{03}} = (28 \cdot 10^3 \cdot 3) \cdot \frac{{{d_RO_m:.3f} + {D_RO_m:.3f}}}{{2}} \cdot {N_p_prime_logical} \cdot {S_vn_actual_m2}',
            'calc_text': f'M03 = (28e3*3) * ({d_RO_m:.3f}+{D_RO_m:.3f})*0.5 * {N_p_prime_logical} * {S_vn_actual_m2}',
            'result': M03_kg
        }

        # (99) Маса обмотки РО з ізоляцією (M'03), кг
        M_prime_03_kg = M03_kg * 1.015
        sec4['(99) Маса обмотки РО з ізоляцією (M\'03), кг'] = {
            'formula': r'M^{\prime}_{03} = M_{03} \cdot 1.015',
            'calculation': fr'M^{{\prime}}_{{03}} = {M03_kg:.2f} \cdot 1.015',
            'calc_text': f'M_prime_03 = {M03_kg:.2f} * 1.015',
            'result': M_prime_03_kg
        }

        # (100) Коефіцієнт (beta_дод3)
        beta_dod3 = (b2_bare_m * 1 * k_p_formula16) / l_VN_m
        sec4['(100) Коефіцієнт (beta_дод3)'] = {
            'formula': r'\beta_{дод3} = \frac{b_{2.гол} \cdot 1 \cdot k_p}{l_{ВН}}',
            'calculation': fr'\beta_{{дод3}} = \frac{{{b2_bare_m} \cdot 1 \cdot {k_p_formula16}}}{{{l_VN_m:.3f}}}',
            'calc_text': f'beta_dod3 = ({b2_bare_m} * 1 * {k_p_formula16}) / {l_VN_m:.3f}',
            'result': beta_dod3
        }

        # (101) Коефіцієнт додаткових втрат (кдод3)
        k_dod3 = 1 + 0.095 * 1e8 * beta_dod3**2 * a2_bare_m**4 * (1**2)
        sec4['(101) Коефіцієнт додаткових втрат (кдод3)'] = {
            'formula': r'k_{дод3} = 1 + 0.095 \cdot 10^8 \beta_{дод3}^2 a_{2.гол}^4 \cdot (1^2)',
            'calculation': fr'k_{{дод3}} = 1 + 0.095 \cdot 10^8 \cdot {beta_dod3:.3f}^2 \cdot {a2_bare_m}^4 \cdot (1^2)',
            'calc_text': f'k_dod3 = 1 + 0.095e8 * {beta_dod3:.3f}**2 * {a2_bare_m}**4 * (1**2)',
            'result': k_dod3
        }

        # (102) Щільність теплового потоку РО, Вт/м^2
        J_RO_A_m2 = J_VN_A_m2
        Q_RO_W_m2 = (107e-10 * J_RO_A_m2 * I_ph_HV_A * 0.25 * k_dod3) / ((b2_prime_ro_m + A_RO_m) * 0.75)
        sec4['(102) Щільність теплового потоку РО, Вт/м^2'] = {
            'formula': r'Q_{РО} = \frac{107 \cdot 10^{-10} \cdot J_{РО} \cdot I_{ф.ВН} \cdot 0.25 \cdot k_{дод3}}{0.75 \cdot (b^{\prime}_{2.РО} + A_{РО})}',
            'calculation': fr'Q_{{РО}} = \frac{{107 \cdot 10^{{-10}} \cdot {J_RO_A_m2:.1f} \cdot {I_ph_HV_A:.2f} \cdot 0.25 \cdot {k_dod3:.3f}}}{{0.75 \cdot ({b2_prime_ro_m} + {A_RO_m})}}',
            'calc_text': f'Q_RO = (107e-10 * {J_RO_A_m2:.1f} * {I_ph_HV_A:.2f} * 0.25 * {k_dod3:.3f}) / (({b2_prime_ro_m} + {A_RO_m}) * 0.75)',
            'result': Q_RO_W_m2
        }

        # ================== РОЗДІЛ 2.4: ВТРАТИ В ОБМОТКАХ ==================
        sec5 = results['section_5']
        # (103) Основні втрати в обмотці НН, Вт
        loss_const = 2.4e-12
        P_NN_W = loss_const * J_NN_A_m2**2 * M01_kg
        sec5['(103) Основні втрати в обмотці НН, Вт'] = {
            'formula': r'P_{НН} = C_{втрат} \cdot J_{НН}^2 \cdot M_{01}',
            'calculation': fr'P_{{НН}} = {loss_const} \cdot ({J_NN_A_m2:.1f})^2 \cdot {M01_kg:.2f}',
            'calc_text': f'P_NN = {loss_const} * {J_NN_A_m2:.1f}**2 * {M01_kg:.2f}',
            'result': P_NN_W
        }

        # (104) Основні втрати в обмотці ВН, Вт
        P_VN_W = loss_const * J_VN_A_m2**2 * M02_kg
        sec5['(104) Основні втрати в обмотці ВН, Вт'] = {
            'formula': r'P_{ВН} = C_{втрат} \cdot J_{ВН}^2 \cdot M_{02}',
            'calculation': fr'P_{{ВН}} = {loss_const} \cdot ({J_VN_A_m2:.1f})^2 \cdot {M02_kg:.2f}',
            'calc_text': f'P_VN = {loss_const} * {J_VN_A_m2:.1f}**2 * {M02_kg:.2f}',
            'result': P_VN_W
        }

        # (105) Основні втрати в обмотці РО, Вт
        P_RO_W = loss_const * J_VN_A_m2**2 * M03_kg
        sec5['(105) Основні втрати в обмотці РО, Вт'] = {
            'formula': r'P_{РО} = C_{втрат} \cdot J_{ВН}^2 \cdot M_{03}',
            'calculation': fr'P_{{РО}} = {loss_const} \cdot ({J_VN_A_m2:.1f})^2 \cdot {M03_kg:.2f}',
            'calc_text': f'P_RO = {loss_const} * {J_VN_A_m2:.1f}**2 * {M03_kg:.2f}',
            'result': P_RO_W
        }

        # (106) Довжина відводів ВН, м
        l_prime_VN_m = 7.5 * l_VN_m
        sec5['(106) Довжина відводів ВН, м'] = {
            'formula': r"l'_{ВН} = 7.5 \cdot l_{ВН}",
            'calculation': fr"l'_{{ВН}} = 7.5 \cdot {l_VN_m:.3f}",
            'calc_text': f"l'_VN = 7.5 * {l_VN_m:.3f}",
            'result': l_prime_VN_m
        }

        # (107) Довжина відводів НН, м
        l_prime_NN_m = 14 * l_NN_m
        sec5['(107) Довжина відводів НН, м'] = {
            'formula': r"l'_{НН} = 14 \cdot l_{НН}",
            'calculation': fr"l'_{{НН}} = 14 \cdot {l_NN_m:.3f}",
            'calc_text': f"l'_NN = 14 * {l_NN_m:.3f}",
            'result': l_prime_NN_m
        }

        # (108) Маса відводів ВН, кг
        gamma_cu = 8900
        M_v1_kg = l_prime_VN_m * S_vn_actual_m2 * gamma_cu
        sec5['(108) Маса відводів ВН, кг'] = {
            'formula': r"M_{в1} = l'_{ВН} \cdot S_{ВН.факт} \cdot \gamma_{cu}",
            'calculation': fr"M_{{в1}} = {l_prime_VN_m:.3f} \cdot {S_vn_actual_m2} \cdot {gamma_cu}",
            'calc_text': f"M_v1 = {l_prime_VN_m:.3f} * {S_vn_actual_m2} * {gamma_cu}",
            'result': M_v1_kg
        }

        # (109) Маса відводів НН, кг
        M_v2_kg = l_prime_NN_m * S_nn_actual_m2 * gamma_cu
        sec5['(109) Маса відводів НН, кг'] = {
            'formula': r"M_{в2} = l'_{НН} \cdot S_{НН.факт} \cdot \gamma_{cu}",
            'calculation': fr"M_{{в2}} = {l_prime_NN_m:.3f} \cdot {S_nn_actual_m2} \cdot {gamma_cu}",
            'calc_text': f"M_v2 = {l_prime_NN_m:.3f} * {S_nn_actual_m2} * {gamma_cu}",
            'result': M_v2_kg
        }

        # (110) Втрати на відводах ВН, Вт
        P_otv_v_W = loss_const * J_VN_A_m2**2 * M_v1_kg * k_dod2
        sec5['(110) Втрати на відводах ВН, Вт'] = {
            'formula': r'P_{від.ВН} = C_{втрат} \cdot J_{ВН}^2 \cdot M_{в1} \cdot k_{дод2}',
            'calculation': fr'P_{{від.ВН}} = {loss_const} \cdot ({J_VN_A_m2:.1f})^2 \cdot {M_v1_kg:.2f} \cdot {k_dod2:.3f}',
            'calc_text': f'P_otv_v = {loss_const} * {J_VN_A_m2:.1f}**2 * {M_v1_kg:.2f} * {k_dod2:.3f}',
            'result': P_otv_v_W
        }

        # (111) Втрати на відводах НН, Вт
        P_otv_n_W = loss_const * J_NN_A_m2**2 * M_v2_kg * k_dod1
        sec5['(111) Втрати на відводах НН, Вт'] = {
            'formula': r'P_{від.НН} = C_{втрат} \cdot J_{НН}^2 \cdot M_{в2} \cdot k_{дод1}',
            'calculation': fr'P_{{від.НН}} = {loss_const} \cdot ({J_NN_A_m2:.1f})^2 \cdot {M_v2_kg:.2f} \cdot {k_dod1:.3f}',
            'calc_text': f'P_otv_n = {loss_const} * {J_NN_A_m2:.1f}**2 * {M_v2_kg:.2f} * {k_dod1:.3f}',
            'result': P_otv_n_W
        }

        # (112) Втрати у баці, Вт
        k_bak = 0.05
        P_bak_W = 10 * k_bak * S_n_VA * 1e-3
        sec5['(112) Втрати у баці, Вт'] = {
            'formula': r'P_{бак} = 10 \cdot k_{бак} \cdot S_{ном} \cdot 10^{-3}',
            'calculation': fr'P_{{бак}} = 10 \cdot {k_bak} \cdot {S_n_VA} \cdot 10^{{-3}}',
            'calc_text': f'P_bak = 10 * {k_bak} * {S_n_VA} * 1e-3',
            'result': P_bak_W
        }

        # (113) Повні втрати КЗ (номінал), Вт
        P_k_prime_W = P_NN_W * k_dod1 + P_VN_W * k_dod2 + P_otv_n_W + P_otv_v_W + P_bak_W
        sec5['(113) Повні втрати КЗ (номінал), Вт'] = {
            'formula': r"P'_{кз} = P_{НН}k_{дод1} + P_{ВН}k_{дод2} + P_{від.НН} + P_{від.ВН} + P_{бак}",
            'calculation': fr"P'_{{кз}} = {P_NN_W:.1f} \cdot {k_dod1:.3f} + {P_VN_W:.1f} \cdot {k_dod2:.3f} + {P_otv_n_W:.1f} + {P_otv_v_W:.1f} + {P_bak_W:.1f}",
            'calc_text': f"P_k' = {P_NN_W:.1f}*{k_dod1:.3f} + {P_VN_W:.1f}*{k_dod2:.3f} + {P_otv_n_W:.1f} + {P_otv_v_W:.1f} + {P_bak_W:.1f}",
            'result': P_k_prime_W
        }

        # (114) Відхилення втрат КЗ від вихідних даних, %
        P_k_initial_W = P_sc_kW * 1000
        deviation_prime_percent = ((P_k_prime_W - P_k_initial_W) / P_k_initial_W) * 100
        sec5['(114) Відхилення втрат КЗ від вихідних даних, %'] = {
            'formula': r"\delta' = \frac{P'_{кз} - P_{кз.поч}}{P_{кз.поч}} \cdot 100\%",
            'calculation': fr"\delta' = \frac{{{P_k_prime_W:.2f} - {P_k_initial_W:.2f}}}{{{P_k_initial_W:.2f}}} \cdot 100\%",
            'calc_text': f"deviation' = (({P_k_prime_W:.2f} - {P_k_initial_W:.2f}) / {P_k_initial_W:.2f}) * 100",
            'result': deviation_prime_percent
        }

        # (115) Повні втрати КЗ (з регулюванням), Вт
        P_k_double_prime_W = P_NN_W * k_dod1 + P_VN_W * k_dod2 + P_RO_W * k_dod3 + P_otv_n_W + P_otv_v_W + P_bak_W
        sec5['(115) Повні втрати КЗ (з регулюванням), Вт'] = {
            'formula': r"P''_{кз} = P_{НН}k_{дод1} + P_{ВН}k_{дод2} + P_{РО}k_{дод3} + P_{від.НН} + P_{від.ВН} + P_{бак}",
            'calculation': fr"P''_{{кз}} = {P_NN_W:.1f} \cdot {k_dod1:.3f} + {P_VN_W:.1f} \cdot {k_dod2:.3f} + {P_RO_W:.1f} \cdot {k_dod3:.3f} + {P_otv_n_W:.1f} + {P_otv_v_W:.1f} + {P_bak_W:.1f}",
            'calc_text': f"P_k'' = {P_NN_W:.1f}*{k_dod1:.3f} + {P_VN_W:.1f}*{k_dod2:.3f} + {P_RO_W:.1f}*{k_dod3:.3f} + {P_otv_n_W:.1f} + {P_otv_v_W:.1f} + {P_bak_W:.1f}",
            'result': P_k_double_prime_W
        }

        # (116) Відхилення втрат КЗ (з регулюванням), %
        deviation_double_prime_percent = ((P_k_double_prime_W - P_k_initial_W) / P_k_initial_W) * 100
        sec5['(116) Відхилення втрат КЗ (з регулюванням), %'] = {
            'formula': r"\delta'' = \frac{P''_{кз} - P_{кз.поч}}{P_{кз.поч}} \cdot 100\%",
            'calculation': fr"\delta'' = \frac{{{P_k_double_prime_W:.2f} - {P_k_initial_W:.2f}}}{{{P_k_initial_W:.2f}}} \cdot 100\%",
            'calc_text': f"deviation'' = (({P_k_double_prime_W:.2f} - {P_k_initial_W:.2f}) / {P_k_initial_W:.2f}) * 100",
            'result': deviation_double_prime_percent
        }
        # ================== РОЗДІЛ 2.5: НАПРУГА КЗ ==================
        sec6 = results['section_6']
        # (117) Активна складова напруги КЗ (uка), %
        u_ka_percent = (P_k_prime_W * 100) / S_n_VA
        sec6['(117) Активна складова напруги КЗ (uка), %'] = {
            'formula': r'u_{ка\%} = \frac{P^{\prime}_{кз} \cdot 100}{S_{ном}}',
            'calculation': fr'u_{{ка\%}} = \frac{{{P_k_prime_W:.2f} \cdot 100}}{{{S_n_VA}}}',
            'calc_text': f'u_ka% = ({P_k_prime_W:.2f} * 100) / {S_n_VA}',
            'result': u_ka_percent
        }

        # (118, 119) Коефіцієнт Роговського (k'p)
        sigma_119 = (A_01 + A_NN_m + A_VN_m) / (math.pi * l_NN_m)
        k_p_prime = 1 - sigma_119 * (1 - math.exp(-1 / sigma_119))
        sec6['(118, 119) Коефіцієнт Роговського (k\'p)'] = {
            'formula': r'\sigma = \frac{A_{01} + A_{НН} + A_{ВН}}{\pi l_{НН}}; \quad k^{\prime}_p = 1 - \sigma(1 - e^{-1/\sigma})',
            'calculation': fr'\sigma = \frac{{{A_01} + {A_NN_m:.3f} + {A_VN_m:.3f}}}{{\pi \cdot {l_NN_m:.3f}}} = {sigma_119:.3f}; \quad k^{{\prime}}_p = 1 - {sigma_119:.3f}(1 - e^{{-1/{sigma_119:.3f}}})',
            'calc_text': f'sigma = ({A_01}+{A_NN_m:.3f}+{A_VN_m:.3f})/(pi*{l_NN_m:.3f})={sigma_119:.3f}; k_p_prime = 1-{sigma_119:.3f}*(1-exp(-1/{sigma_119:.3f}))',
            'result': k_p_prime
        }

        # (120) Ширина приведеного каналу розсіювання (ap), м
        numerator_ap = (D_NN_m + A_01) * A_01 + ((D_NN_m + d_NN_m) / 2) * (A_NN_m / 3) + ((D_VN_m + d_VN_m) / 2) * (A_VN_m / 3)
        a_p_m = numerator_ap / (D_NN_m + A_01)
        sec6['(120) Ширина приведеного каналу розсіювання (ap), м'] = {
            'formula': r'a_p = \frac{(D_{НН} + A_{01})A_{01} + \frac{D_{НН}+d_{НН}}{2}\frac{A_{НН}}{3} + \frac{D_{ВН}+d_{ВН}}{2}\frac{A_{ВН}}{3}}{D_{НН} + A_{01}}',
            'calculation': fr'a_p = \frac{{{numerator_ap:.4f}}}{{{D_NN_m:.3f} + {A_01}}}',
            'calc_text': f'a_p = (({D_NN_m:.3f}+{A_01})*{A_01} + ...)/({D_NN_m:.3f}+{A_01})',
            'result': a_p_m
        }

        # (121) Середній діаметр каналу розсіювання (d12), м
        d12_m = (D_NN_m + d_VN_m) / 2
        sec6['(121) Середній діаметр каналу розсіювання (d12), м'] = {
            'formula': r'd_{12} = \frac{D_{НН} + d_{ВН}}{2}',
            'calculation': fr'd_{{12}} = \frac{{{D_NN_m:.3f} + {d_VN_m:.3f}}}{{2}}',
            'calc_text': f'd12 = ({D_NN_m:.3f} + {d_VN_m:.3f}) / 2',
            'result': d12_m
        }

        # (122) Коефіцієнт (beta1)
        beta1 = (math.pi * d12_m) / l_VN_m
        sec6['(122) Коефіцієнт (beta1)'] = {
            'formula': r'\beta_1 = \frac{\pi \cdot d_{12}}{l_{ВН}}',
            'calculation': fr'\beta_1 = \frac{{\pi \cdot {d12_m:.3f}}}{{{l_VN_m:.3f}}}',
            'calc_text': f'beta1 = (pi * {d12_m:.3f}) / {l_VN_m:.3f}',
            'result': beta1
        }

        # (123) Коефіцієнт розриву (X)
        X_123 = 0.02 / l_VN_m
        sec6['(123) Коефіцієнт розриву (X)'] = {
            'formula': r'X = \frac{0.02}{l_{ВН}}',
            'calculation': fr'X = \frac{{0.02}}{{{l_VN_m:.3f}}}',
            'calc_text': f'X = 0.02 / {l_VN_m:.3f}',
            'result': X_123
        }

        # (124) Коефіцієнт (kq)
        kq = 1 + (l_VN_m * X_123**2) / (3 * a_p_m * k_p_prime)
        sec6['(124) Коефіцієнт (kq)'] = {
            'formula': r'k_q = 1 + \frac{l_{ВН} \cdot X^2}{3 \cdot a_p \cdot k^{\prime}_p}',
            'calculation': fr'k_q = 1 + \frac{{{l_VN_m:.3f} \cdot {X_123:.3f}^2}}{{3 \cdot {a_p_m:.3f} \cdot {k_p_prime:.3f}}}',
            'calc_text': f'kq = 1 + ({l_VN_m:.3f} * {X_123:.3f}**2) / (3 * {a_p_m:.3f} * {k_p_prime:.3f})',
            'result': kq
        }

        # (125) Реактивна складова напруги КЗ (uкр), %
        u_kr_percent = (7.9 * f_hz * beta1 * S_st_kVA * a_p_m * k_p_prime * kq) / (10 * e_v_57**2)
        sec6['(125) Реактивна складова напруги КЗ (uкр), %'] = {
            'formula': r'u_{кр\%} = \frac{7.9 \cdot f \cdot \beta_1 \cdot S_{ст} \cdot a_p \cdot k^{\prime}_p \cdot k_q}{10 \cdot e_в^2}',
            'calculation': fr'u_{{кр\%}} = \frac{{7.9 \cdot {f_hz} \cdot {beta1:.3f} \cdot {S_st_kVA} \cdot {a_p_m:.3f} \cdot {k_p_prime:.3f} \cdot {kq:.3f}}}{{10 \cdot {e_v_57:.3f}^2}}',
            'calc_text': f'u_kr% = (7.9*{f_hz}*{beta1:.3f}*{S_st_kVA}*{a_p_m:.3f}*{k_p_prime:.3f}*{kq:.3f})/(10*{e_v_57:.3f}**2)',
            'result': u_kr_percent
        }

        # (126) Розрахована напруга КЗ (u'k), %
        u_k_prime_percent = math.sqrt(u_ka_percent**2 + u_kr_percent**2)
        sec6['(126) Розрахована напруга КЗ (u\'k), %'] = {
            'formula': r"u'_{к\%} = \sqrt{u_{ка\%}^2 + u_{кр\%}^2}",
            'calculation': fr"u'_{{к\%}} = \sqrt{{{u_ka_percent:.2f}^2 + {u_kr_percent:.2f}^2}}",
            'calc_text': f"u_k'% = sqrt({u_ka_percent:.2f}**2 + {u_kr_percent:.2f}**2)",
            'result': u_k_prime_percent
        }

        # (127) Відхилення напруги КЗ від вихідних даних, %
        deviation_uk_percent = ((u_k_prime_percent - ukz_percent) / ukz_percent) * 100
        sec6['(127) Відхилення напруги КЗ від вихідних даних, %'] = {
            'formula': r"\delta_{uk} = \frac{u'_{к\%} - u_{к\%.поч}}{u_{к\%.поч}} \cdot 100\%",
            'calculation': fr"\delta_{{uk}} = \frac{{{u_k_prime_percent:.2f} - {ukz_percent:.2f}}}{{{ukz_percent:.2f}}} \cdot 100\%",
            'calc_text': f"deviation_uk = (({u_k_prime_percent:.2f} - {ukz_percent:.2f}) / {ukz_percent:.2f}) * 100",
            'result': deviation_uk_percent
        }
        # ================== РОЗДІЛ 2.6: МЕХАНІЧНА СТІЙКІСТЬ ==================
        sec7 = results['section_7']
        # (128) Сталий струм КЗ (IКУ), А
        S_kz_VA = 2500e6
        I_KU_A = 100 * I_ph_HV_A / (u_k_prime_percent * (1 + ((S_n_VA * 100) / (10 * u_k_prime_percent * S_kz_VA))))
        sec7['(128) Сталий струм КЗ (IКУ), А'] = {
            'formula': r'I_{КУ} = \frac{100 \cdot I_{ф.ВН}}{u^{\prime}_{к\%} \left(1 + \frac{S_{ном} \cdot 100}{10 \cdot u^{\prime}_{к\%} \cdot S_{кз}}\right)}',
            'calculation': fr'I_{{КУ}} = \frac{{100 \cdot {I_ph_HV_A:.2f}}}{{{u_k_prime_percent:.2f} \left(1 + \frac{{{S_n_VA} \cdot 100}}{{10 \cdot {u_k_prime_percent:.2f} \cdot {S_kz_VA}}}\right)}}',
            'calc_text': f'I_KU = (100*{I_ph_HV_A:.2f}) / ({u_k_prime_percent:.2f}*(1+(({S_n_VA}*100)/(10*{u_k_prime_percent:.2f}*{S_kz_VA}))))',
            'result': I_KU_A
        }

        # (129) Коефіцієнт (kmax)
        k_max = 1 + math.exp(-math.pi * (u_ka_percent / u_kr_percent))
        sec7['(129) Коефіцієнт (kmax)'] = {
            'formula': r'k_{max} = 1 + e^{-\pi \frac{u_{ка\%}}{u_{кр\%}}}',
            'calculation': fr'k_{{max}} = 1 + e^{{-\pi \frac{{{u_ka_percent:.2f}}}{{{u_kr_percent:.2f}}}}}',
            'calc_text': f'k_max = 1 + exp(-pi * ({u_ka_percent:.2f} / {u_kr_percent:.2f}))',
            'result': k_max
        }

        # (130) Ударний струм КЗ (Іуд), А
        I_ud_A = math.sqrt(2) * k_max * I_KU_A
        sec7['(130) Ударний струм КЗ (Іуд), А'] = {
            'formula': r'I_{уд} = \sqrt{2} \cdot k_{max} \cdot I_{КУ}',
            'calculation': fr'I_{{уд}} = \sqrt{{2}} \cdot {k_max:.3f} \cdot {I_KU_A:.2f}',
            'calc_text': f'I_ud = sqrt(2) * {k_max:.3f} * {I_KU_A:.2f}',
            'result': I_ud_A
        }

        # (131) Радіальна сила (Fp), Н
        F_p_N = 0.628 * (I_ud_A * N_vn_adopted)**2 * beta1 * k_p_prime * 1e-6
        sec7['(131) Радіальна сила (Fp), Н'] = {
            'formula': r'F_p = 0.628 \cdot (I_{уд} \cdot N_{ВН})^2 \cdot \beta_1 \cdot k^{\prime}_p \cdot 10^{-6}',
            'calculation': fr'F_p = 0.628 \cdot ({I_ud_A:.2f} \cdot {N_vn_adopted})^2 \cdot {beta1:.3f} \cdot {k_p_prime:.3f} \cdot 10^{{-6}}',
            'calc_text': f'F_p = 0.628 * ({I_ud_A:.2f} * {N_vn_adopted})**2 * {beta1:.3f} * {k_p_prime:.3f} * 1e-6',
            'result': F_p_N
        }

        # (132) Розтягувальна напруга НН, МПа
        sigma_p_NN_Pa = F_p_N / (2 * math.pi * N_nn_adopted * S_nn_actual_m2)
        sec7['(132) Розтягувальна напруга НН, МПа'] = {
            'formula': r'\sigma_{р.НН} = \frac{F_p}{2\pi \cdot N_{НН} \cdot S_{НН.факт}}',
            'calculation': fr'\sigma_{{р.НН}} = \frac{{{F_p_N:.2f}}}{{2\pi \cdot {N_nn_adopted} \cdot {S_nn_actual_m2}}}',
            'calc_text': f'sigma_p_NN = {F_p_N:.2f} / (2*pi * {N_nn_adopted} * {S_nn_actual_m2})',
            'result': sigma_p_NN_Pa / 1e6
        }

        # (133) Стискальна напруга ВН, МПа
        sigma_c_VN_Pa = F_p_N / (2 * math.pi * N_vn_adopted * S_vn_actual_m2)
        sec7['(133) Стискальна напруга ВН, МПа'] = {
            'formula': r'\sigma_{с.ВН} = \frac{F_p}{2\pi \cdot N_{ВН} \cdot S_{ВН.факт}}',
            'calculation': fr'\sigma_{{с.ВН}} = \frac{{{F_p_N:.2f}}}{{2\pi \cdot {N_vn_adopted} \cdot {S_vn_actual_m2}}}',
            'calc_text': f'sigma_c_VN = {F_p_N:.2f} / (2*pi * {N_vn_adopted} * {S_vn_actual_m2})',
            'result': sigma_c_VN_Pa / 1e6
        }

        # (134) Осьова сила (Fос1), Н
        F_os1_N = F_p_N * a_p_m / (2 * l_VN_m)
        sec7['(134) Осьова сила (Fос1), Н'] = {
            'formula': r'F_{ос1} = \frac{F_p \cdot a_p}{2 \cdot l_{ВН}}',
            'calculation': fr'F_{{ос1}} = \frac{{{F_p_N:.2f} \cdot {a_p_m:.3f}}}{{2 \cdot {l_VN_m:.3f}}}',
            'calc_text': f'F_os1 = ({F_p_N:.2f} * {a_p_m:.3f}) / (2 * {l_VN_m:.3f})',
            'result': F_os1_N
        }

        # (135) Осьова сила (Fос2), Н
        F_os2_N = (F_p_N * (0.02)) / (0.35 * k_p_prime * 4)
        sec7['(135) Осьова сила (Fос2), Н'] = {
            'formula': r'F_{ос2} = \frac{F_p \cdot 0.02}{0.35 \cdot k^{\prime}_p \cdot 4}',
            'calculation': fr'F_{{ос2}} = \frac{{{F_p_N:.2f} \cdot 0.02}}{{0.35 \cdot {k_p_prime:.3f} \cdot 4}}',
            'calc_text': f'F_os2 = ({F_p_N:.2f} * 0.02) / (0.35 * {k_p_prime:.3f} * 4)',
            'result': F_os2_N
        }

        # (136) Зусилля на обмотці НН, Н
        F_NN_N = F_os1_N + F_os2_N
        sec7['(136) Зусилля на обмотці НН, Н'] = {
            'formula': r'F_{НН} = F_{ос1} + F_{ос2}',
            'calculation': fr'F_{{НН}} = {F_os1_N:.2f} + {F_os2_N:.2f}',
            'calc_text': f'F_NN = {F_os1_N:.2f} + {F_os2_N:.2f}',
            'result': F_NN_N
        }

        # (137) Зусилля на обмотці ВН, Н
        F_VN_N = F_os1_N - F_os2_N
        sec7['(137) Зусилля на обмотці ВН, Н'] = {
            'formula': r'F_{ВН} = F_{ос1} - F_{ос2}',
            'calculation': fr'F_{{ВН}} = {F_os1_N:.2f} - {F_os2_N:.2f}',
            'calc_text': f'F_VN = {F_os1_N:.2f} - {F_os2_N:.2f}',
            'result': F_VN_N
        }

        # (138) Напруга стискання НН, МПа
        sigma_szh1_Pa = F_NN_N / (A_NN_m * 20 * 0.05)
        sec7['(138) Напруга стискання НН, МПа'] = {
            'formula': r'\sigma_{сж.НН} = \frac{F_{НН}}{A_{НН} \cdot N_{прокл} \cdot b_{прокл}}',
            'calculation': fr'\sigma_{{сж.НН}} = \frac{{{F_NN_N:.2f}}}{{{A_NN_m:.3f} \cdot 20 \cdot 0.05}}',
            'calc_text': f'sigma_szh1 = {F_NN_N:.2f} / ({A_NN_m:.3f} * 20 * 0.05)',
            'result': sigma_szh1_Pa / 1e6
        }

        # (139) Напруга стискання ВН, МПа
        sigma_szh2_Pa = F_VN_N / (A_VN_m * 20 * 0.05)
        sec7['(139) Напруга стискання ВН, МПа'] = {
            'formula': r'\sigma_{сж.ВН} = \frac{F_{ВН}}{A_{ВН} \cdot N_{прокл} \cdot b_{прокл}}',
            'calculation': fr'\sigma_{{сж.ВН}} = \frac{{{F_VN_N:.2f}}}{{{A_VN_m:.3f} \cdot 20 \cdot 0.05}}',
            'calc_text': f'sigma_szh2 = {F_VN_N:.2f} / ({A_VN_m:.3f} * 20 * 0.05)',
            'result': sigma_szh2_Pa / 1e6
        }

        # ================== РОЗДІЛ 2.7: ТЕМПЕРАТУРНА СТІЙКІСТЬ ==================
        sec8 = results['section_8']
        # (140) Гранична умовна температура обмотки, °С
        T_kz = 4.0
        t_um_C = (670 * T_kz) / (12.5 * (u_k_prime_percent / (J_VN_A_m2/1e6))**2) - T_kz + 90.0
        sec8['(140) Гранична умовна температура обмотки, °С'] = {
            'formula': r'\theta_{ум} = \frac{670 \cdot T_{кз}}{12.5 \left(\frac{u^{\prime}_{к\%}}{J_{ВН}}\right)^2} - T_{кз} + 90',
            'calculation': fr'\theta_{{ум}} = \frac{{670 \cdot {T_kz}}}{{12.5 \left(\frac{{{u_k_prime_percent:.2f}}}{{{J_VN_A_m2/1e6:.3f}}}\right)^2}} - {T_kz} + 90',
            'calc_text': f't_um = (670*{T_kz})/(12.5*({u_k_prime_percent:.2f}/{J_VN_A_m2/1e6:.3f})**2) - {T_kz} + 90',
            'result': t_um_C
        }

        # (141) Час досягнення 250°С, с
        T250_s = 2.5 * (u_k_prime_percent / (J_VN_A_m2/1e6))**2
        sec8['(141) Час досягнення 250°С, с'] = {
            'formula': r'T_{250} = 2.5 \left(\frac{u^{\prime}_{к\%}}{J_{ВН}}\right)^2',
            'calculation': fr'T_{{250}} = 2.5 \left(\frac{{{u_k_prime_percent:.2f}}}{{{J_VN_A_m2/1e6:.3f}}}\right)^2',
            'calc_text': f'T_250 = 2.5 * ({u_k_prime_percent:.2f} / {J_VN_A_m2/1e6:.3f})**2',
            'result': T250_s
        }

        # ================== РОЗДІЛ 2.8: КОНСТРУКЦІЯ КІСТЯКА ==================
        sec9 = results['section_9']
        # (142) Активний переріз стрижня, м^2
        S_core_active_m2 = 1451.2 * k_zk * 1e-4
        sec9['(142) Активний переріз стрижня, м^2'] = {
            'formula': r'S_{с.акт} = 1451.2 \cdot k_{зк} \cdot 10^{-4}',
            'calculation': fr'S_{{с.акт}} = 1451.2 \cdot {k_zk} \cdot 10^{{-4}}',
            'calc_text': f'S_core_active = 1451.2 * {k_zk} * 1e-4',
            'result': S_core_active_m2
        }

        # (143) Активний переріз ярма, м^2
        S_yoke_active_m2 = 1460.2 * k_zk * 1e-4
        sec9['(143) Активний переріз ярма, м^2'] = {
            'formula': r'S_{я.акт} = 1460.2 \cdot k_{зк} \cdot 10^{-4}',
            'calculation': fr'S_{{я.акт}} = 1460.2 \cdot {k_zk} \cdot 10^{{-4}}',
            'calc_text': f'S_yoke_active = 1460.2 * {k_zk} * 1e-4',
            'result': S_yoke_active_m2
        }

        # (144) Площа на косому стику, м^2
        S_kos_m2 = S_core_active_m2 * math.sqrt(2)
        sec9['(144) Площа на косому стику, м^2'] = {
            'formula': r'S_{кос} = S_{с.акт} \cdot \sqrt{2}',
            'calculation': fr'S_{{кос}} = {S_core_active_m2:.4f} \cdot \sqrt{{2}}',
            'calc_text': f'S_kos = {S_core_active_m2:.4f} * sqrt(2)',
            'result': S_kos_m2
        }

        # (145) Об’єм сталі у куті, м^3
        V_st_m3 = k_zk * 55860.0 * 1e-6
        sec9['(145) Об’єм сталі у куті, м^3'] = {
            'formula': r'V_{ст.кута} = k_{зк} \cdot 55860.0 \cdot 10^{-6}',
            'calculation': fr'V_{{ст.кута}} = {k_zk} \cdot 55860.0 \cdot 10^{{-6}}',
            'calc_text': f'V_st = {k_zk} * 55860.0 * 1e-6',
            'result': V_st_m3
        }

        # (146) Довжина стрижня, м
        l_ct_m = l_VN_m + A_10 + A_20
        sec9['(146) Довжина стрижня, м'] = {
            'formula': r'l_{стр} = l_{ВН} + A_{10} + A_{20}',
            'calculation': fr'l_{{стр}} = {l_VN_m:.3f} + {A_10} + {A_20}',
            'calc_text': f'l_ct = {l_VN_m:.3f} + {A_10} + {A_20}',
            'result': l_ct_m
        }

        # (147) Прийнята відстань між осями, м
        L_c_adopted_m = 0.865
        sec9['(147) Прийнята відстань між осями, м'] = {
            'formula': r'L_c',
            'calculation': r'\text{Прийняте значення}',
            'calc_text': 'Adopted value',
            'result': L_c_adopted_m
        }

        # (148) Маса сталі у кутах, кг
        gamma_steel = 7650
        Gy_kg = gamma_steel * V_st_m3
        sec9['(148) Маса сталі у кутах, кг'] = {
            'formula': r'G_у = \gamma_{ст} \cdot V_{ст.кута}',
            'calculation': fr'G_у = {gamma_steel} \cdot {V_st_m3:.6f}',
            'calc_text': f'Gy = {gamma_steel} * {V_st_m3:.6f}',
            'result': Gy_kg
        }

        # (149) Маса сталі ярем, кг
        G_yoke_kg = 4 * S_yoke_active_m2 * L_c_adopted_m * gamma_steel + 2 * Gy_kg
        sec9['(149) Маса сталі ярем, кг'] = {
            'formula': r'G_{я} = 4 \cdot S_{я.акт} \cdot L_c \cdot \gamma_{ст} + 2 \cdot G_у',
            'calculation': fr'G_{{я}} = 4 \cdot {S_yoke_active_m2:.4f} \cdot {L_c_adopted_m} \cdot {gamma_steel} + 2 \cdot {Gy_kg:.2f}',
            'calc_text': f'G_yoke = 4 * {S_yoke_active_m2:.4f} * {L_c_adopted_m} * {gamma_steel} + 2 * {Gy_kg:.2f}',
            'result': G_yoke_kg
        }

        # (150) Маса сталі стрижнів, кг
        G_core_kg = 3 * l_ct_m * S_core_active_m2 * gamma_steel + S_core_active_m2 * gamma_steel * 0.440 - Gy_kg
        sec9['(150) Маса сталі стрижнів, кг'] = {
            'formula': r'G_{с} = 3 l_{стр} S_{с.акт} \gamma_{ст} + S_{с.акт} \gamma_{ст} \cdot 0.440 - G_у',
            'calculation': fr'G_{{с}} = 3 \cdot {l_ct_m:.3f} \cdot {S_core_active_m2:.4f} \cdot {gamma_steel} + {S_core_active_m2:.4f} \cdot {gamma_steel} \cdot 0.440 - {Gy_kg:.2f}',
            'calc_text': f'G_core = 3*{l_ct_m:.3f}*{S_core_active_m2:.4f}*{gamma_steel} + {S_core_active_m2:.4f}*{gamma_steel}*0.440 - {Gy_kg:.2f}',
            'result': G_core_kg
        }

        # (151) Повна маса сталі, кг
        G_steel_total_kg = G_core_kg + G_yoke_kg
        sec9['(151) Повна маса сталі, кг'] = {
            'formula': r'G_{сталі} = G_с + G_я',
            'calculation': fr'G_{{сталі}} = {G_core_kg:.2f} + {G_yoke_kg:.2f}',
            'calc_text': f'G_steel_total = {G_core_kg:.2f} + {G_yoke_kg:.2f}',
            'result': G_steel_total_kg
        }

        # ================== РОЗДІЛ 2.9: ВТРАТИ НЕРОБОЧОГО ХОДУ ==================
        sec10 = results['section_10']
        # (152) Розрахована індукція у стрижні, Тл
        B_c_recalc = e_v_57 / (4.44 * f_hz * S_core_active_m2)
        sec10['(152) Розрахована індукція у стрижні, Тл'] = {
            'formula': r'B_c = \frac{e_{в.ут}}{4.44 \cdot f \cdot S_{с.акт}}',
            'calculation': fr'B_c = \frac{{{e_v_57:.3f}}}{{4.44 \cdot {f_hz} \cdot {S_core_active_m2:.4f}}}',
            'calc_text': f'B_c_recalc = {e_v_57:.3f} / (4.44 * {f_hz} * {S_core_active_m2:.4f})',
            'result': B_c_recalc
        }

        # (153) Розрахована індукція у ярмі, Тл
        B_ya_recalc = e_v_57 / (4.44 * f_hz * S_yoke_active_m2)
        sec10['(153) Розрахована індукція у ярмі, Тл'] = {
            'formula': r'B_я = \frac{e_{в.ут}}{4.44 \cdot f \cdot S_{я.акт}}',
            'calculation': fr'B_я = \frac{{{e_v_57:.3f}}}{{4.44 \cdot {f_hz} \cdot {S_yoke_active_m2:.4f}}}',
            'calc_text': f'B_ya_recalc = {e_v_57:.3f} / (4.44 * {f_hz} * {S_yoke_active_m2:.4f})',
            'result': B_ya_recalc
        }

        # (154) Розрахована індукція на косому стику, Тл
        B_ks_recalc = B_c_recalc / math.sqrt(2)
        sec10['(154) Розрахована індукція на косому стику, Тл'] = {
            'formula': r'B_{кс} = \frac{B_c}{\sqrt{2}}',
            'calculation': fr'B_{{кс}} = \frac{{{B_c_recalc:.3f}}}{{\sqrt{{2}}}}',
            'calc_text': f'B_ks_recalc = {B_c_recalc:.3f} / sqrt(2)',
            'result': B_ks_recalc
        }

        # (155) Розраховані втрати НХ, Вт
        P_nh_prime_W = (1.025 * 1.0 * (0.843*G_core_kg + 0.824*G_yoke_kg - 4*0.824*Gy_kg + 0.5*(0.824+0.843)*8.4*Gy_kg) + 6*S_kos_m2*440.366)*1.0*1.05*1.09
        sec10['(155) Розраховані втрати НХ, Вт'] = {
            'formula': r"P'_{НХ} = (P_{питом} \cdot G_{сталі} + P_{стик}) \cdot k_{обр} \cdot k_{зб} \cdot k_{зап}",
            'calculation': fr"P'_{{НХ}} = (1.025 \cdot (...) + 6 \cdot {S_kos_m2:.4f} \cdot 440.366) \cdot 1.05 \cdot 1.09",
            'calc_text': f"P_nh_prime = (1.025*(0.843*{G_core_kg:.2f} + ...) + 6*{S_kos_m2:.4f}*440.366)*1.05*1.09",
            'result': P_nh_prime_W
        }

        # (156) Відхилення втрат НХ, %
        P_nh_initial_W = P_nl_kW * 1000
        deviation_P_nh_percent = ((P_nh_prime_W - P_nh_initial_W) / P_nh_initial_W) * 100
        sec10['(156) Відхилення втрат НХ, %'] = {
            'formula': r"\delta_{НХ} = \frac{P'_{НХ} - P_{НХ.поч}}{P_{НХ.поч}} \cdot 100\%",
            'calculation': fr"\delta_{{НХ}} = \frac{{{P_nh_prime_W:.2f} - {P_nh_initial_W:.2f}}}{{{P_nh_initial_W:.2f}}} \cdot 100\%",
            'calc_text': f"deviation_P_nh = (({P_nh_prime_W:.2f} - {P_nh_initial_W:.2f}) / {P_nh_initial_W:.2f}) * 100",
            'result': deviation_P_nh_percent
        }
        
        # ================== РОЗДІЛ 2.10: СТРУМ НЕРОБОЧОГО ХОДУ ==================
        sec11 = results['section_11']
        # (157) Намагнічуюча потужність (Qx), ВА
        Qx_VA = (1.225*1.0*(1.798*G_core_kg+1.697*G_yoke_kg-4*1.697*Gy_kg+0.5*(1.798+1.697)*27.95*1.25*Gy_kg)+6*S_kos_m2*2119.5)*1.03*1.06*1.09
        sec11['(157) Намагнічуюча потужність (Qx), ВА'] = {
            'formula': r"Q_x = (q_{питом} \cdot G_{сталі} + Q_{стик}) \cdot k_{...}",
            'calculation': fr"Q_x = (1.225 \cdot (...) + 6 \cdot {S_kos_m2:.4f} \cdot 2119.5) \cdot 1.03 \cdot 1.06 \cdot 1.09",
            'calc_text': f"Qx = (1.225*(1.798*{G_core_kg:.2f} + ...) + 6*{S_kos_m2:.4f}*2119.5)*1.03*1.06*1.09",
            'result': Qx_VA
        }

        # (158) Розрахований струм НХ (i'0), %
        i0_prime_percent = (Qx_VA * 100) / S_n_VA
        sec11['(158) Розрахований струм НХ (i\'0), %'] = {
            'formula': r"i'_{0\%} = \frac{Q_x \cdot 100}{S_{ном}}",
            'calculation': fr"i'_{{0\%}} = \frac{{{Qx_VA:.2f} \cdot 100}}{{{S_n_VA}}}",
            'calc_text': f"i'0% = ({Qx_VA:.2f} * 100) / {S_n_VA}",
            'result': i0_prime_percent
        }

        # (159) Відхилення струму НХ, %
        deviation_i0_percent = ((i0_prime_percent - i0_percent) / i0_percent) * 100
        sec11['(159) Відхилення струму НХ, %'] = {
            'formula': r"\delta_{i0} = \frac{i'_{0\%} - i_{0\%.поч}}{i_{0\%.поч}} \cdot 100\%",
            'calculation': fr"\delta_{{i0}} = \frac{{{i0_prime_percent:.3f} - {i0_percent:.3f}}}{{{i0_percent:.3f}}} \cdot 100\%",
            'calc_text': f"deviation_i0 = (({i0_prime_percent:.3f} - {i0_percent:.3f}) / {i0_percent:.3f}) * 100",
            'result': deviation_i0_percent
        }

        # (160) Активна складова струму НХ, %
        i0a_percent = (P_nh_prime_W * 100) / S_n_VA
        sec11['(160) Активна складова струму НХ, %'] = {
            'formula': r"i_{0а\%} = \frac{P'_{НХ} \cdot 100}{S_{ном}}",
            'calculation': fr"i_{{0а\%}} = \frac{{{P_nh_prime_W:.2f} \cdot 100}}{{{S_n_VA}}}",
            'calc_text': f"i0a% = ({P_nh_prime_W:.2f} * 100) / {S_n_VA}",
            'result': i0a_percent
        }

        # (161) Реактивна складова струму НХ, %
        i0p_percent = math.sqrt(i0_prime_percent**2 - i0a_percent**2)
        sec11['(161) Реактивна складова струму НХ, %'] = {
            'formula': r"i_{0р\%} = \sqrt{i'^{2}_{0\%} - i^2_{0а\%}}",
            'calculation': fr"i_{{0р\%}} = \sqrt{{{i0_prime_percent:.3f}^2 - {i0a_percent:.3f}^2}}",
            'calc_text': f"i0p% = sqrt({i0_prime_percent:.3f}**2 - {i0a_percent:.3f}**2)",
            'result': i0p_percent
        }

        # (162) ККД (номінальний режим), %
        eta_nominal = (1 - (P_k_prime_W + P_nh_prime_W) / (P_k_prime_W + P_nh_prime_W + S_n_VA)) * 100
        sec11['(162) ККД (номінальний режим), %'] = {
            'formula': r"\eta = \left(1 - \frac{P'_{кз} + P'_{НХ}}{P'_{кз} + P'_{НХ} + S_{ном}}\right) \cdot 100\%",
            'calculation': fr"\eta = \left(1 - \frac{{{P_k_prime_W:.2f} + {P_nh_prime_W:.2f}}}{{{P_k_prime_W:.2f} + {P_nh_prime_W:.2f} + {S_n_VA}}}\right) \cdot 100\%",
            'calc_text': f"eta_nom = (1 - ({P_k_prime_W:.2f}+{P_nh_prime_W:.2f})/({P_k_prime_W:.2f}+{P_nh_prime_W:.2f}+{S_n_VA}))*100",
            'result': eta_nominal
        }

        # (163) ККД (з повним регулюванням), %
        eta_full_reg = (1 - (P_k_double_prime_W + P_nh_prime_W) / (P_k_double_prime_W + P_nh_prime_W + S_n_VA)) * 100
        sec11['(163) ККД (з повним регулюванням), %'] = {
            'formula': r"\eta_{рег} = \left(1 - \frac{P''_{кз} + P'_{НХ}}{P''_{кз} + P'_{НХ} + S_{ном}}\right) \cdot 100\%",
            'calculation': fr"\eta_{{рег}} = \left(1 - \frac{{{P_k_double_prime_W:.2f} + {P_nh_prime_W:.2f}}}{{{P_k_double_prime_W:.2f} + {P_nh_prime_W:.2f} + {S_n_VA}}}\right) \cdot 100\%",
            'calc_text': f"eta_reg = (1-({P_k_double_prime_W:.2f}+{P_nh_prime_W:.2f})/({P_k_double_prime_W:.2f}+{P_nh_prime_W:.2f}+{S_n_VA}))*100",
            'result': eta_full_reg
        }

        # ================== РОЗДІЛ 2.11: ВЕКТОРНА ДІАГРАМА  ==================
        sec12 = results['section_12']

        # (165) Коефіцієнт трансформації
        k_t = U_ph_HV_V / U_ph_LV_V
        sec12['(165) Коефіцієнт трансформації'] = {
            'formula': r'k_t = \frac{U_{ф.ВН}}{U_{ф.НН}}',
            'calculation': fr'k_t = \frac{{{U_ph_HV_V}}}{{{U_ph_LV_V}}}',
            'calc_text': f'k_t = {U_ph_HV_V} / {U_ph_LV_V}',
            'result': k_t
        }

        # (166) Напруга КЗ, В
        U_kz_V = U_ph_HV_V * u_k_prime_percent / 100
        sec12['(166) Напруга КЗ, В'] = {
            'formula': r"U_{кз} = U_{ф.ВН} \cdot \frac{u'_{к\%}}{100}",
            'calculation': fr"U_{{кз}} = {U_ph_HV_V} \cdot \frac{{{u_k_prime_percent:.2f}}}{{100}}",
            'calc_text': f"U_kz = {U_ph_HV_V} * {u_k_prime_percent:.2f} / 100",
            'result': U_kz_V
        }

        # (167) Повний опір КЗ, Ом
        I_kz_A = I_ph_HV_A
        Z_k_Ohm = U_kz_V / I_kz_A
        sec12['(167) Повний опір КЗ, Ом'] = {
            'formula': r'Z_k = \frac{U_{кз}}{I_{кз}}',
            'calculation': fr'Z_k = \frac{{{U_kz_V:.2f}}}{{{I_kz_A:.2f}}}',
            'calc_text': f'Z_k = {U_kz_V:.2f} / {I_kz_A:.2f}',
            'result': Z_k_Ohm
        }

        # (168) Активний опір КЗ, Ом
        R_k_Ohm = P_k_prime_W / (m_phases * I_kz_A ** 2)
        sec12['(168) Активний опір КЗ, Ом'] = {
            'formula': r"R_k = \frac{P'_{кз}}{m \cdot I_{кз}^2}",
            'calculation': fr"R_k = \frac{{{P_k_prime_W:.2f}}}{{{m_phases} \cdot {I_kz_A:.2f}^2}}",
            'calc_text': f"R_k = {P_k_prime_W:.2f} / ({m_phases} * {I_kz_A:.2f}**2)",
            'result': R_k_Ohm
        }

        # (169) Реактивний опір КЗ, Ом
        X_k_Ohm = math.sqrt(Z_k_Ohm ** 2 - R_k_Ohm ** 2)
        sec12['(169) Реактивний опір КЗ, Ом'] = {
            'formula': r'X_k = \sqrt{Z_k^2 - R_k^2}',
            'calculation': fr'X_k = \sqrt{{{Z_k_Ohm:.3f}^2 - {R_k_Ohm:.3f}^2}}',
            'calc_text': f'X_k = sqrt({Z_k_Ohm:.3f}**2 - {R_k_Ohm:.3f}**2)',
            'result': X_k_Ohm
        }

        # (170) Повний опір магнітного контуру, Ом
        E1x = U_ph_HV_V
        Zm_Ohm = E1x / i0_prime_percent
        sec12['(170) Повний опір магнітного контуру, Ом'] = {
            'formula': r'Z_m = \frac{E_1}{I_0}',
            'calculation': fr'Z_m = \frac{{{E1x}}}{{{i0_prime_percent:.3f}}}',
            'calc_text': f'Zm = {E1x} / {i0_prime_percent:.3f}',
            'result': Zm_Ohm
        }

        # (171) Активний опір магнітного контуру, Ом
        Rm_Ohm = P_nh_prime_W / (m_phases * i0_prime_percent**2) if i0_prime_percent > 0 else float('inf')
        sec12['(171) Активний опір магнітного контуру, Ом'] = {
            'formula': r'R_m = \frac{P^{\prime}_{НХ}}{m \cdot I_0^2}',
            'calculation': fr'R_m = \frac{{{P_nh_prime_W:.2f}}}{{{m_phases} \cdot {i0_prime_percent:.3f}^2}}',
            'calc_text': f'Rm = {P_nh_prime_W:.2f} / ({m_phases} * {i0_prime_percent:.3f}**2)',
            'result': Rm_Ohm
        }

        # (172) Реактивний опір магнітного контуру, Ом
        Xm_Ohm = math.sqrt(Zm_Ohm**2 - Rm_Ohm**2) if Zm_Ohm**2 > Rm_Ohm**2 else 0
        sec12['(172) Реактивний опір магнітного контуру, Ом'] = {
            'formula': r'X_m = \sqrt{Z_m^2 - R_m^2}',
            'calculation': fr'X_m = \sqrt{{{Zm_Ohm:.2f}^2 - {Rm_Ohm:.2f}^2}}',
            'calc_text': f'Xm = sqrt({Zm_Ohm:.2f}**2 - {Rm_Ohm:.2f}**2)',
            'result': Xm_Ohm
        }

        # (173) Опір R1 = R'2, Ом
        R1_Ohm = R_k_Ohm / 2
        R2_prime_Ohm = R1_Ohm
        sec12['(173) Опори R1 = R\'2, Ом'] = {
            'formula': r"R_1=R'_2 = \frac{R_k}{2}",
            'calculation': fr"R_1=R'_2 = \frac{{{R_k_Ohm:.3f}}}{{2}}",
            'calc_text': f'R1 = R\'2 = {R_k_Ohm:.3f} / 2',
            'result': R1_Ohm
        }

        # (174) Опір X1 = X'2, Ом
        X1_Ohm = X_k_Ohm / 2
        X2_prime_Ohm = X1_Ohm
        sec12['(174) Опори X1 = X\'2, Ом'] = {
            'formula': r"X_1=X'_2 = \frac{X_k}{2}",
            'calculation': fr"X_1=X'_2 = \frac{{{X_k_Ohm:.3f}}}{{2}}",
            'calc_text': f'X1 = X\'2 = {X_k_Ohm:.3f} / 2',
            'result': X1_Ohm
        }

        # (175) Опір Z1 = Z'2, Ом
        Z1_Ohm = Z_k_Ohm / 2
        Z2_prime_Ohm = Z1_Ohm
        sec12['(175) Опори Z1 = Z\'2, Ом'] = {
            'formula': r"Z_1=Z'_2 = \frac{Z_k}{2}",
            'calculation': fr"Z_1=Z'_2 = \frac{{{Z_k_Ohm:.3f}}}{{2}}",
            'calc_text': f'Z1 = Z\'2 = {Z_k_Ohm:.3f} / 2',
            'result': Z1_Ohm
        }



        mF, mi, mU = 0.002, 4.0, 125.0

        # (176) Магнітний потік (Ф), Вб
        Flux_Wb = B_c_recalc * S_core_active_m2
        sec12['(176) Магнітний потік (Ф), Вб'] = {
            'formula': r'\Phi = B_c \cdot S_{с.акт}',
            'calculation': fr'\Phi = {B_c_recalc:.3f} \cdot {S_core_active_m2:.4f}',
            'calc_text': f'Flux = {B_c_recalc:.3f} * {S_core_active_m2:.4f}',
            'result': Flux_Wb
        }

        # (177) Масштабований магнітний потік, мм
        Flux_m = Flux_Wb / mF
        sec12['(177) Масштабований магнітний потік, мм'] = {
            'formula': r'\Phi_m = \frac{\Phi}{m_F}',
            'calculation': fr'\Phi_m = \frac{{{Flux_Wb:.4f}}}{{{mF}}}',
            'calc_text': f'Flux_m = {Flux_Wb:.4f} / {mF}',
            'result': Flux_m
        }

        # (178) Приведений струм НН (I'2), А
        I2_prime_A = I_ph_LV_A / k_t
        sec12['(178) Приведений струм НН (I\'2), А'] = {
            'formula': r"I'_2 = \frac{I_{ф.НН}}{k_t}",
            'calculation': fr"I'_2 = \frac{{{I_ph_LV_A:.2f}}}{{{k_t:.3f}}}",
            'calc_text': f"I'2 = {I_ph_LV_A:.2f} / {k_t:.3f}",
            'result': I2_prime_A
        }

        # (179) Масштабований струм НН, мм
        I2_prime_m = I2_prime_A / mi
        sec12['(179) Масштабований струм НН, мм'] = {
            'formula': r"I'_{2.m} = \frac{I'_2}{m_I}",
            'calculation': fr"I'_{{2.m}} = \frac{{{I2_prime_A:.2f}}}{{{mi}}}",
            'calc_text': f"I'2_m = {I2_prime_A:.2f} / {mi}",
            'result': I2_prime_m
        }

        # (180) Масштабований струм НХ, мм
        sec12['(180) Масштабований струм НХ, мм'] = {
            'formula': r"I_{1xm} = \frac{i'_0}{m_i}",
            'calculation': fr'I_{{1xm}} = \frac{{{i0_prime_percent:.3f}}}{{{mi}}}',
            'calc_text': f'i0_m = \frsc{i0_prime_percent} / {mi}',
            'result': i0_prime_percent / mi
        }

        # (181) Масштабована реактивна складова струму НХ, мм
        sec12['(181) Масштабована реактивна складова струму НХ, мм'] = {
            'formula': r'I_{1xpm} = \frac{i_{0р}}{m_i}',
            'calculation': fr'I_{{0р.m}} = \frac{{{i0p_percent:.3f}}}{{{mi}}}',
            'calc_text': f'i0p_m = \frac({i0p_percent}) / {mi}',
            'result': i0p_percent / mi
        }

        # (182) Масштабована активна складова струму НХ, мм
        sec12['(182) Масштабована активна складова струму НХ, мм'] = {
            'formula': r'I_{1xam} = \frac{i_{0а}}{m_i}',
            'calculation': fr'I_{{1xam}} = \frac{{{i0a_percent:.3f}}}{{{mi}}}',
            'calc_text': f'i0a_m = \frac({i0a_percent}) / {mi}',
            'result': i0a_percent / mi
        }

        # (183) Кут магнітного запізнювання, °
        gamma_deg = math.degrees(math.atan(i0a_percent / i0p_percent)) if i0p_percent > 0 else 0
        sec12['(183) Кут магнітного запізнювання, °'] = {
            'formula': r'\gamma = \arctan\left(\frac{i_{0а\%}}{i_{0р\%}}\right)',
            'calculation': fr'\gamma = \arctan\left(\frac{{{i0a_percent:.3f}}}{{{i0p_percent:.3f}}}\right)',
            'calc_text': f'gamma = degrees(atan({i0a_percent:.3f}/{i0p_percent:.3f}))',
            'result': gamma_deg
        }


        # (184) Масштабована приведена напруга НН, мм
        U2_prime_V = U_ph_LV_V * k_t
        U2_prime_m = U2_prime_V / mU
        sec12['(184) Масштабована приведена напруга НН, мм'] = {
            'formula': r"U'_{2.m} = \frac{U'_{2}}{m_U}",
            'calculation': fr"U'_{{2.m}} = \frac{{{U2_prime_V:.2f}}}{{{mU}}}",
            'calc_text': f"U'2_m = ({U_ph_LV_V}*{k_t:.3f}) / {mU}",
            'result': U2_prime_m
        }

        # (185) Масштабоване активне падіння напруги, мм
        I2R2_drop_m = (I2_prime_A * R2_prime_Ohm) / mU
        sec12['(185) Масштабоване активне падіння напруги, мм'] = {
            'formula': r"\Delta U_{a.m} = \frac{I'_2 R'_2}{m_U}",
            'calculation': fr"\Delta U_{{a.m}} = \frac{{{I2_prime_A:.2f} \cdot {R2_prime_Ohm:.3f}}}{{{mU}}}",
            'calc_text': f"I2R2_m = ({I2_prime_A:.2f} * {R2_prime_Ohm:.3f}) / {mU}",
            'result': I2R2_drop_m
        }

        # (186) Масштабоване реактивне падіння напруги, мм
        I2X2_drop_m = (I2_prime_A * X2_prime_Ohm) / mU
        sec12['(186) Масштабоване реактивне падіння напруги, мм'] = {
            'formula': r"\Delta U_{р.m} = \frac{I'_2 X'_2}{m_U}",
            'calculation': fr"\Delta U_{{р.m}} = \frac{{{I2_prime_A:.2f} \cdot {X2_prime_Ohm:.3f}}}{{{mU}}}",
            'calc_text': f"I2X2_m = ({I2_prime_A:.2f} * {X2_prime_Ohm:.3f}) / {mU}",
            'result': I2X2_drop_m
        }

        # (187) Кут навантаження
        cos_phi2 = 0.98
        phi2_deg = math.degrees(math.acos(cos_phi2))
        sec12[f'(187) Кут навантаження для cos(phi2)={cos_phi2}, °'] = {
            'formula': r'\phi_2 = \arccos(\text{cos}\phi_2)',
            'calculation': fr'\phi_2 = \arccos({cos_phi2})',
            'calc_text': f'phi2 = degrees(acos({cos_phi2}))',
            'result': phi2_deg
        }


        # (189) u1ka, В
        u1ka_V = U_kz_V * u_ka_percent / 100
        sec12['(189) u1ka, В'] = {
            'formula': r'U_{1ка} = U_{кз} \cdot \frac{u_{ка\%}}{100}',
            'calculation': fr'U_{{1ка}} = {U_kz_V:.2f} \cdot \frac{{{u_ka_percent:.2f}}}{{100}}',
            'calc_text': f'u1ka = {U_kz_V:.2f} * {u_ka_percent:.2f} / 100',
            'result': u1ka_V
        }

        # (190) u1kp, В
        u1kp_V = U_kz_V * u_kr_percent / 100
        sec12['(190) u1kp, В'] = {
            'formula': r'U_{1кр} = U_{кз} \cdot \frac{u_{кр\%}}{100}',
            'calculation': fr'U_{{1кр}} = {U_kz_V:.2f} \cdot \frac{{{u_kr_percent:.2f}}}{{100}}',
            'calc_text': f'u1kp = {U_kz_V:.2f} * {u_kr_percent:.2f} / 100',
            'result': u1kp_V
        } 

        
        # ================== РОЗДІЛ 2.12: РОЗРАХУНОК БАКУ ==================
        sec13 = results['section_13']
        # (195) Внутрішній перепад t° НН, °С
        lambda_iz = 0.17
        delta_iz1_mm = 0.5
        Theta_01_C = (Q_NN_W_m2 * delta_iz1_mm * 1e-3) / (2 * lambda_iz)
        sec13['(195) Внутрішній перепад t° НН, °С'] = {
            'formula': r'\Theta_{01} = \frac{Q_{НН} \cdot \delta_{із1} \cdot 10^{-3}}{2 \lambda_{із}}',
            'calculation': fr'\Theta_{{01}} = \frac{{{Q_NN_W_m2:.2f} \cdot {delta_iz1_mm} \cdot 10^{{-3}}}}{{2 \cdot {lambda_iz}}}',
            'calc_text': f'Theta_01 = ({Q_NN_W_m2:.2f} * {delta_iz1_mm} * 1e-3) / (2 * {lambda_iz})',
            'result': Theta_01_C
        }

        # (196) Внутрішній перепад t° ВН, °С
        delta_iz2_mm = 0.5
        Theta_02_C = (Q_VN_W_m2 * delta_iz2_mm * 1e-3) / (2 * lambda_iz)
        sec13['(196) Внутрішній перепад t° ВН, °С'] = {
            'formula': r'\Theta_{02} = \frac{Q_{ВН} \cdot \delta_{із2} \cdot 10^{-3}}{2 \lambda_{із}}',
            'calculation': fr'\Theta_{{02}} = \frac{{{Q_VN_W_m2:.2f} \cdot {delta_iz2_mm} \cdot 10^{{-3}}}}{{2 \cdot {lambda_iz}}}',
            'calc_text': f'Theta_02 = ({Q_VN_W_m2:.2f} * {delta_iz2_mm} * 1e-3) / (2 * {lambda_iz})',
            'result': Theta_02_C
        }

        # (197) Внутрішній перепад t° РО, °С
        delta_iz3_mm = 1.35
        Theta_03_C = (Q_RO_W_m2 * delta_iz3_mm * 1e-3) / (2 * lambda_iz)
        sec13['(197) Внутрішній перепад t° РО, °С'] = {
            'formula': r'\Theta_{03} = \frac{Q_{РО} \cdot \delta_{із3} \cdot 10^{-3}}{2 \lambda_{із}}',
            'calculation': fr'\Theta_{{03}} = \frac{{{Q_RO_W_m2:.2f} \cdot {delta_iz3_mm} \cdot 10^{{-3}}}}{{2 \cdot {lambda_iz}}}',
            'calc_text': f'Theta_03 = ({Q_RO_W_m2:.2f} * {delta_iz3_mm} * 1e-3) / (2 * {lambda_iz})',
            'result': Theta_03_C
        }

        # (198) Перепад t° на поверхні НН, °С
        Theta_om1_C = 1.0 * 1.1 * 1.05 * 0.35 * Q_NN_W_m2**0.6
        sec13['(198) Перепад t° на поверхні НН, °С'] = {
            'formula': r'\Theta_{ом1} = K_1 \cdot K_2 \cdot K_3 \cdot C \cdot Q_{НН}^{0.6}',
            'calculation': fr'\Theta_{{ом1}} = 1.0 \cdot 1.1 \cdot 1.05 \cdot 0.35 \cdot {Q_NN_W_m2:.2f}^{{0.6}}',
            'calc_text': f'Theta_om1 = 1.0 * 1.1 * 1.05 * 0.35 * {Q_NN_W_m2:.2f}**0.6',
            'result': Theta_om1_C
        }

        # (199) Перепад t° на поверхні ВН, °С
        Theta_om2_C = 1.0 * 1.0 * 1.1 * 0.35 * Q_VN_W_m2**0.6
        sec13['(199) Перепад t° на поверхні ВН, °С'] = {
            'formula': r'\Theta_{ом2} = K_1 \cdot K_2 \cdot K_3 \cdot C \cdot Q_{ВН}^{0.6}',
            'calculation': fr'\Theta_{{ом2}} = 1.0 \cdot 1.0 \cdot 1.1 \cdot 0.35 \cdot {Q_VN_W_m2:.2f}^{{0.6}}',
            'calc_text': f'Theta_om2 = 1.0 * 1.0 * 1.1 * 0.35 * {Q_VN_W_m2:.2f}**0.6',
            'result': Theta_om2_C
        }

        # (200) Перепад t° на поверхні РО, °С
        Theta_om3_C = 1.0 * 1.0 * 0.8 * 0.35 * Q_RO_W_m2**0.6
        sec13['(200) Перепад t° на поверхні РО, °С'] = {
            'formula': r'\Theta_{ом3} = K_1 \cdot K_2 \cdot K_3 \cdot C \cdot Q_{РО}^{0.6}',
            'calculation': fr'\Theta_{{ом3}} = 1.0 \cdot 1.0 \cdot 0.8 \cdot 0.35 \cdot {Q_RO_W_m2:.2f}^{{0.6}}',
            'calc_text': f'Theta_om3 = 1.0 * 1.0 * 0.8 * 0.35 * {Q_RO_W_m2:.2f}**0.6',
            'result': Theta_om3_C
        }

        # (201) Повний перепад НН -> масло, °С
        Theta_omsr1_C = Theta_om1_C + Theta_01_C
        sec13['(201) Повний перепад НН -> масло, °С'] = {
            'formula': r'\Theta_{ом.сер1} = \Theta_{ом1} + \Theta_{01}',
            'calculation': fr'\Theta_{{ом.сер1}} = {Theta_om1_C:.2f} + {Theta_01_C:.2f}',
            'calc_text': f'Theta_omsr1 = {Theta_om1_C:.2f} + {Theta_01_C:.2f}',
            'result': Theta_omsr1_C
        }

        # (202) Повний перепад ВН -> масло, °С
        Theta_omsr2_C = Theta_om2_C + Theta_02_C
        sec13['(202) Повний перепад ВН -> масло, °С'] = {
            'formula': r'\Theta_{ом.сер2} = \Theta_{ом2} + \Theta_{02}',
            'calculation': fr'\Theta_{{ом.сер2}} = {Theta_om2_C:.2f} + {Theta_02_C:.2f}',
            'calc_text': f'Theta_omsr2 = {Theta_om2_C:.2f} + {Theta_02_C:.2f}',
            'result': Theta_omsr2_C
        }

        # (203) Повний перепад РО -> масло, °С
        Theta_omsr3_C = Theta_om3_C + Theta_03_C
        sec13['(203) Повний перепад РО -> масло, °С'] = {
            'formula': r'\Theta_{ом.сер3} = \Theta_{ом3} + \Theta_{03}',
            'calculation': fr'\Theta_{{ом.сер3}} = {Theta_om3_C:.2f} + {Theta_03_C:.2f}',
            'calc_text': f'Theta_omsr3 = {Theta_om3_C:.2f} + {Theta_03_C:.2f}',
            'result': Theta_omsr3_C
        }

        # (204) Мінімальна ширина баку, м
        s1, s2, s3, s4, delta1, delta2 = 0.05, 0.05, 0.025, 0.025, 0.02, 0.02
        B_m = D_RO_m + s1 + s2 + s3 + s4 + delta1 + delta2
        sec13['(204) Мінімальна ширина баку, м'] = {
            'formula': r'B_{мін} = D_{РО} + s_1+s_2+s_3+s_4+\delta_1+\delta_2',
            'calculation': fr'B_{{мін}} = {D_RO_m:.3f} + {s1}+{s2}+{s3}+{s4}+{delta1}+{delta2}',
            'calc_text': f'B_min = {D_RO_m:.3f} + {s1}+{s2}+{s3}+{s4}+{delta1}+{delta2}',
            'result': B_m
        }

        # (доп.) Прийнята ширина баку, м
        B_adopted_m = 1.44
        sec13['(доп.) Прийнята ширина баку, м'] = { 'formula': r'B_{прийн}', 'calc_text': 'Adopted value', 'result': B_adopted_m}

        # (205) Довжина баку, м
        A_m = 2 * L_c_adopted_m + B_adopted_m + 0.75
        sec13['(205) Довжина баку, м'] = {
            'formula': r'A = 2L_c + B_{прийн} + 0.75',
            'calculation': fr'A = 2 \cdot {L_c_adopted_m} + {B_adopted_m} + 0.75',
            'calc_text': f'A = 2 * {L_c_adopted_m} + {B_adopted_m} + 0.75',
            'result': A_m
        }

        # (доп.) Прийнята довжина баку, м
        A_adopted_m = 3.94
        sec13['(доп.) Прийнята довжина баку, м'] = { 'formula': r'A_{прийн}', 'calc_text': 'Adopted value', 'result': A_adopted_m}

        # (206) Висота активної частини, м
        n_bya_m = 30 / 1000
        H_ach_m = l_ct_m + 2 * 0.44 + n_bya_m
        sec13['(206) Висота активної частини, м'] = {
            'formula': r'H_{ач} = l_{стр} + 2 \cdot 0.44 + n_{бя}',
            'calculation': fr'H_{{ач}} = {l_ct_m:.3f} + 2 \cdot 0.44 + {n_bya_m}',
            'calc_text': f'H_ach = {l_ct_m:.3f} + 2 * 0.44 + {n_bya_m}',
            'result': H_ach_m
        }

        # (207) Загальна висота баку, м
        H_b_m = H_ach_m + 300 * 1e-3
        sec13['(207) Загальна висота баку, м'] = {
            'formula': r'H_б = H_{ач} + 0.3',
            'calculation': fr'H_б = {H_ach_m:.3f} + 0.3',
            'calc_text': f'H_b = {H_ach_m:.3f} + 0.3',
            'result': H_b_m
        }

        # (доп.) Прийнята висота баку, м
        Hb_adopted_m = 3.14
        sec13['(доп.) Прийнята висота баку, м'] = { 'formula': r'H_{б.прийн}', 'calc_text': 'Adopted value', 'result': Hb_adopted_m}

        # ================== РОЗДІЛ 2.13: ТЕПЛОВИЙ РОЗРАХУНОК ==================
        sec14 = results['section_14']
        # (208) Припустиме підвищення t° масла, °С
        Theta_mv_C = 65 - Theta_omsr1_C
        sec14['(208) Припустиме підвищення t° масла, °С'] = {
            'formula': r'\Theta_{м.в.доп} = 65 - \Theta_{ом.сер1}',
            'calculation': fr'\Theta_{{м.в.доп}} = 65 - {Theta_omsr1_C:.2f}',
            'calc_text': f'Theta_mv = 65 - {Theta_omsr1_C:.2f}',
            'result': Theta_mv_C
        }

        # (209) Попереднє перевищення t° стінки баку, °С
        Theta_bv_C = Theta_mv_C - 6 - 2
        sec14['(209) Попереднє перевищення t° стінки баку, °С'] = {
            'formula': r'\Theta_{б.в.поп} = \Theta_{м.в.доп} - 8',
            'calculation': fr'\Theta_{{б.в.поп}} = {Theta_mv_C:.2f} - 8',
            'calc_text': f'Theta_bv = {Theta_mv_C:.2f} - 8',
            'result': Theta_bv_C
        }

        # (210) Попереднє перевищення t° масла у верхніх шарах, °С
        Theta_mvv_C = 1.2 * Theta_mv_C
        sec14['(210) Попереднє перевищення t° масла у верхніх шарах, °С'] = {
            'formula': r'\Theta_{м.в.в.поп} = 1.2 \cdot \Theta_{м.в.доп}',
            'calculation': fr'\Theta_{{м.в.в.поп}} = 1.2 \cdot {Theta_mv_C:.2f}',
            'calc_text': f'Theta_mvv = 1.2 * {Theta_mv_C:.2f}',
            'result': Theta_mvv_C
        }

        # (211) Площа гладкої стінки баку, м^2
        perimeter = 2 * (A_adopted_m - B_adopted_m) + math.pi * B_adopted_m
        P_glb_m2 = Hb_adopted_m * perimeter * 1.6
        sec14['(211) Площа гладкої стінки баку, м^2'] = {
            'formula': r'P_{гл.б} = H_б \cdot (2(A-B) + \pi B) \cdot 1.6',
            'calculation': fr'P_{{гл.б}} = {Hb_adopted_m:.2f} \cdot (2({A_adopted_m}-{B_adopted_m}) + \pi \cdot {B_adopted_m}) \cdot 1.6',
            'calc_text': f'P_glb = {Hb_adopted_m:.2f} * (2*({A_adopted_m}-{B_adopted_m})+pi*{B_adopted_m}) * 1.6',
            'result': P_glb_m2
        }

        # (212) Площа випромінювання баку, м^2
        P_vip_m2 = P_glb_m2 * 1.75
        sec14['(212) Площа випромінювання баку, м^2'] = {
            'formula': r'P_{випр} = P_{гл.б} \cdot 1.75',
            'calculation': fr'P_{{випр}} = {P_glb_m2:.2f} \cdot 1.75',
            'calc_text': f'P_vip = {P_glb_m2:.2f} * 1.75',
            'result': P_vip_m2
        }

        # (213) Необхідна поверхня конвекції, м^2
        numerator_Pk = 1.05 * (P_nh_prime_W + P_k_prime_W)
        denominator_Pk = 2.5 * Theta_bv_C**1.25
        P_k_m2 = (numerator_Pk / denominator_Pk) - (1.12 * P_vip_m2)
        sec14['(213) Необхідна поверхня конвекції, м^2'] = {
            'formula': r'P_к = \frac{1.05 (P^{\prime}_{НХ} + P^{\prime}_{кз})}{2.5 \cdot \Theta_{б.в}^{1.25}} - 1.12 \cdot P_{випр}',
            'calculation': fr'P_к = \frac{{{numerator_Pk:.2f}}}{{2.5 \cdot {Theta_bv_C:.2f}^{{1.25}}}} - 1.12 \cdot {P_vip_m2:.2f}',
            'calc_text': f'P_k = ({numerator_Pk:.2f} / (2.5*{Theta_bv_C:.2f}**1.25)) - (1.12*{P_vip_m2:.2f})',
            'result': P_k_m2
        }

        # (214) Поверхня конвекції кришки, м^2
        P_kkr_m2 = 0.5 * ((A_adopted_m - B_adopted_m)*(B_adopted_m+0.16) + (math.pi*((B_adopted_m+0.16)**2 / 4)))*1.6
        sec14['(214) Поверхня конвекції кришки, м^2'] = {
            'formula': r'P_{к.кр} = 0.5 \left( (A-B)(B+0.16) + \frac{\pi(B+0.16)^2}{4} \right) \cdot 1.6',
            'calculation': fr'P_{{к.кр}} = 0.5 \cdot (({A_adopted_m}-{B_adopted_m}) \cdot ... ) \cdot 1.6',
            'calc_text': f'P_kkr = 0.5*(({A_adopted_m}-{B_adopted_m})*({B_adopted_m}+0.16) + (pi*(({B_adopted_m}+0.16)**2)/4))*1.6',
            'result': P_kkr_m2
        }

        # (215) Необхідна поверхня конвекції радіаторів, м^2
        P_rad_m2 = P_k_m2 - P_kkr_m2 - P_glb_m2
        sec14['(215) Необхідна поверхня конвекції радіаторів, м^2'] = {
            'formula': r'P_{рад} = P_к - P_{к.кр} - P_{гл.б}',
            'calculation': fr'P_{{рад}} = {P_k_m2:.2f} - {P_kkr_m2:.2f} - {P_glb_m2:.2f}',
            'calc_text': f'P_rad = {P_k_m2:.2f} - {P_kkr_m2:.2f} - {P_glb_m2:.2f}',
            'result': P_rad_m2
        }

        # (216) Приведена поверхня одного радіатору, м^2
        P_kr_m2 = (0.908 + 42.831) * 2.24
        sec14['(216) Приведена поверхня одного радіатору, м^2'] = {
            'formula': r'P_{к.1рад} = (C_1 + C_2) \cdot C_3',
            'calculation': r'P_{{к.1рад}} = (0.908 + 42.831) \cdot 2.24',
            'calc_text': 'P_kr = (0.908 + 42.831) * 2.24',
            'result': P_kr_m2
        }

        # (217) Кількість радіаторів
        n_rad_calc = P_rad_m2 / P_kr_m2
        n_rad_adopted = math.ceil(n_rad_calc)
        sec14['(217) Кількість радіаторів'] = {
            'formula': r'n_{рад} = \lceil \frac{P_{рад}}{P_{к.1рад}} \rceil',
            'calculation': fr'n_{{рад}} = \lceil \frac{{{P_rad_m2:.2f}}}{{{P_kr_m2:.2f}}} \rceil',
            'calc_text': f'n_rad = ceil({P_rad_m2:.2f} / {P_kr_m2:.2f})',
            'result': n_rad_adopted
        }

        # (218) Повна поверхня конвекції, м^2
        P_sum_m2 = P_kkr_m2 + P_glb_m2 + P_kr_m2 * n_rad_adopted
        sec14['(218) Повна поверхня конвекції, м^2'] = {
            'formula': r'P_{сум} = P_{к.кр} + P_{гл.б} + n_{рад} \cdot P_{к.1рад}',
            'calculation': fr'P_{{сум}} = {P_kkr_m2:.2f} + {P_glb_m2:.2f} + {n_rad_adopted} \cdot {P_kr_m2:.2f}',
            'calc_text': f'P_sum = {P_kkr_m2:.2f} + {P_glb_m2:.2f} + {n_rad_adopted} * {P_kr_m2:.2f}',
            'result': P_sum_m2
        }

        # (219) Уточнене перевищення t° стінки баку, °С
        numerator_Theta_bv_recalc = 1.05 * (P_nh_prime_W + P_k_prime_W)
        denominator_Theta_bv_recalc = 2.8 * P_vip_m2 + 2.5 * P_sum_m2
        Theta_bv_recalc_C = (numerator_Theta_bv_recalc / denominator_Theta_bv_recalc)**0.8
        sec14['(219) Уточнене перевищення t° стінки баку, °С'] = {
            'formula': r'\Theta_{б.в.ут} = \left( \frac{1.05 (P^{\prime}_{НХ} + P^{\prime}_{кз})}{2.8 P_{випр} + 2.5 P_{сум}} \right)^{0.8}',
            'calculation': fr'\Theta_{{б.в.ут}} = \left( \frac{{{numerator_Theta_bv_recalc:.2f}}}{{2.8 \cdot {P_vip_m2:.2f} + 2.5 \cdot {P_sum_m2:.2f}}} \right)^{{0.8}}',
            'calc_text': f'Theta_bv_recalc = ({numerator_Theta_bv_recalc:.2f} / (2.8*{P_vip_m2:.2f}+2.5*{P_sum_m2:.2f}))**0.8',
            'result': Theta_bv_recalc_C
        }

        # (220) Перевищення t° масла поблизу стінки, °С
        numerator_Theta_mb = 1.05 * (P_nh_prime_W + P_k_prime_W)
        denominator_Theta_mb = (P_kkr_m2 + P_glb_m2 + 42.831 * n_rad_adopted)
        Theta_mb_C = (numerator_Theta_mb / denominator_Theta_mb)**0.6 * 0.165 * 0.9
        sec14['(220) Перевищення t° масла поблизу стінки, °С'] = {
            'formula': r'\Theta_{м.б} = \left( \frac{1.05(P^{\prime}_{НХ}+P^{\prime}_{кз})}{P_{к.кр}+P_{гл.б}+C \cdot n_{рад}} \right)^{0.6} \cdot 0.165 \cdot 0.9',
            'calculation': fr'\Theta_{{м.б}} = \left( \frac{{{numerator_Theta_mb:.2f}}}{{{denominator_Theta_mb:.2f}}} \right)^{{0.6}} \cdot 0.165 \cdot 0.9',
            'calc_text': f'Theta_mb = ({numerator_Theta_mb:.2f}/{denominator_Theta_mb:.2f})**0.6*0.165*0.9',
            'result': Theta_mb_C
        }

        # (221) Перевищення середньої t° масла над повітрям, °С
        Theta_prime_mvv_C = Theta_mb_C + Theta_bv_recalc_C
        sec14['(221) Перевищення середньої t° масла над повітрям, °С'] = {
            'formula': r"\Theta'_{м.в.в} = \Theta_{м.б} + \Theta_{б.в.ут}",
            'calculation': fr"\Theta'_{{м.в.в}} = {Theta_mb_C:.2f} + {Theta_bv_recalc_C:.2f}",
            'calc_text': f"Theta'_mvv = {Theta_mb_C:.2f} + {Theta_bv_recalc_C:.2f}",
            'result': Theta_prime_mvv_C
        }

        # (222) Перевищення t° масла у верхніх шарах, °С
        Theta_mvv_recalc_C = 1.2 * Theta_prime_mvv_C
        sec14['(222) Перевищення t° масла у верхніх шарах, °С'] = {
            'formula': r"\Theta_{м.в.в.ут} = 1.2 \cdot \Theta'_{м.в.в}",
            'calculation': fr"\Theta_{{м.в.в.ут}} = 1.2 \cdot {Theta_prime_mvv_C:.2f}",
            'calc_text': f"Theta_mvv_recalc = 1.2 * {Theta_prime_mvv_C:.2f}",
            'result': Theta_mvv_recalc_C
        }

        # (223) Повне перевищення t° обмотки НН, °С
        Theta_ov1_C = Theta_prime_mvv_C + Theta_om1_C + Theta_01_C
        sec14['(223) Повне перевищення t° обмотки НН, °С'] = {
            'formula': r"\Theta_{ов1} = \Theta'_{м.в.в} + \Theta_{ом1} + \Theta_{01}",
            'calculation': fr"\Theta_{{ов1}} = {Theta_prime_mvv_C:.2f} + {Theta_om1_C:.2f} + {Theta_01_C:.2f}",
            'calc_text': f"Theta_ov1 = {Theta_prime_mvv_C:.2f} + {Theta_om1_C:.2f} + {Theta_01_C:.2f}",
            'result': Theta_ov1_C
        }

        # (224) Повне перевищення t° обмотки ВН, °С
        Theta_ov2_C = Theta_prime_mvv_C + Theta_om2_C + Theta_02_C
        sec14['(224) Повне перевищення t° обмотки ВН, °С'] = {
            'formula': r"\Theta_{ов2} = \Theta'_{м.в.в} + \Theta_{ом2} + \Theta_{02}",
            'calculation': fr"\Theta_{{ов2}} = {Theta_prime_mvv_C:.2f} + {Theta_om2_C:.2f} + {Theta_02_C:.2f}",
            'calc_text': f"Theta_ov2 = {Theta_prime_mvv_C:.2f} + {Theta_om2_C:.2f} + {Theta_02_C:.2f}",
            'result': Theta_ov2_C
        }

        # (225) Повне перевищення t° обмотки РО, °С
        Theta_ov3_C = Theta_prime_mvv_C + Theta_om3_C + Theta_03_C
        sec14['(225) Повне перевищення t° обмотки РО, °С'] = {
            'formula': r"\Theta_{ов3} = \Theta'_{м.в.в} + \Theta_{ом3} + \Theta_{03}",
            'calculation': fr"\Theta_{{ов3}} = {Theta_prime_mvv_C:.2f} + {Theta_om3_C:.2f} + {Theta_03_C:.2f}",
            'calc_text': f"Theta_ov3 = {Theta_prime_mvv_C:.2f} + {Theta_om3_C:.2f} + {Theta_03_C:.2f}",
            'result': Theta_ov3_C
        }
        
        # ================== РОЗДІЛ 2.14: ВИЗНАЧЕННЯ МАСИ ==================
        sec15 = results['section_15']

        # (226) Об’єм овального баку, м^3
        tank_area_m2 = (A_adopted_m - B_adopted_m) * B_adopted_m + (B_adopted_m**2 * math.pi) / 4
        V_tank_m3 = tank_area_m2 * Hb_adopted_m
        sec15['(226) Об’єм овального баку, м^3'] = {
            'formula': r'V_{бак} = \left( (A-B)B + \frac{\pi B^2}{4} \right) \cdot H_б',
            'calculation': fr'V_{{бак}} = \left( ({A_adopted_m}-{B_adopted_m}){B_adopted_m} + \frac{{\pi \cdot {B_adopted_m}^2}}{{4}} \right) \cdot {Hb_adopted_m}',
            'calc_text': f'V_tank = (({A_adopted_m}-{B_adopted_m})*{B_adopted_m} + ({B_adopted_m}**2*pi)/4) * {Hb_adopted_m}',
            'result': V_tank_m3
        }

        # (227) Маса активної частини, кг
        G_ach_kg = 1.2 * (G_steel_total_kg + M_prime_01_kg + M_prime_02_kg + M_prime_03_kg)
        sec15['(227) Маса активної частини, кг'] = {
            'formula': r"G_{ач} = 1.2 (G_{сталі} + M'_{01} + M'_{02} + M'_{03})",
            'calculation': fr"G_{{ач}} = 1.2 ({G_steel_total_kg:.2f} + {M_prime_01_kg:.2f} + {M_prime_02_kg:.2f} + {M_prime_03_kg:.2f})",
            'calc_text': f"G_ach = 1.2 * ({G_steel_total_kg:.2f}+{M_prime_01_kg:.2f}+{M_prime_02_kg:.2f}+{M_prime_03_kg:.2f})",
            'result': G_ach_kg
        }

        # (228) Об’єм активної частини, м^3
        V_ach_m3 = G_ach_kg / 5700
        sec15['(228) Об’єм активної частини, м^3'] = {
            'formula': r'V_{ач} = \frac{G_{ач}}{\rho_{ач}}',
            'calculation': fr'V_{{ач}} = \frac{{{G_ach_kg:.2f}}}{{5700}}',
            'calc_text': f'V_ach = {G_ach_kg:.2f} / 5700',
            'result': V_ach_m3
        }

        # (229) Об’єм масла у баку, м^3
        V_oil_tank_m3 = V_tank_m3 - V_ach_m3
        sec15['(229) Об’єм масла у баку, м^3'] = {
            'formula': r'V_{масла.бак} = V_{бак} - V_{ач}',
            'calculation': fr'V_{{масла.бак}} = {V_tank_m3:.3f} - {V_ach_m3:.3f}',
            'calc_text': f'V_oil_tank = {V_tank_m3:.3f} - {V_ach_m3:.3f}',
            'result': V_oil_tank_m3
        }

        # (230) Маса масла у баку, кг
        G_oil_tank_kg = 900 * V_oil_tank_m3
        sec15['(230) Маса масла у баку, кг'] = {
            'formula': r'G_{масла.бак} = \rho_{масла} \cdot V_{масла.бак}',
            'calculation': fr'G_{{масла.бак}} = 900 \cdot {V_oil_tank_m3:.3f}',
            'calc_text': f'G_oil_tank = 900 * {V_oil_tank_m3:.3f}',
            'result': G_oil_tank_kg
        }

        # (231) Маса масла у радіаторах, кг
        G_oil_radiators_kg = 497.75 * n_rad_adopted
        sec15['(231) Маса масла у радіаторах, кг'] = {
            'formula': r'G_{масла.рад} = C_{масла} \cdot n_{рад}',
            'calculation': fr'G_{{масла.рад}} = 497.75 \cdot {n_rad_adopted}',
            'calc_text': f'G_oil_radiators = 497.75 * {n_rad_adopted}',
            'result': G_oil_radiators_kg
        }

        # (232) Загальна маса масла, кг
        G_oil_total_kg = G_oil_radiators_kg + G_oil_tank_kg
        sec15['(232) Загальна маса масла, кг'] = {
            'formula': r'G_{масла.заг} = G_{масла.рад} + G_{масла.бак}',
            'calculation': fr'G_{{масла.заг}} = {G_oil_radiators_kg:.2f} + {G_oil_tank_kg:.2f}',
            'calc_text': f'G_oil_total = {G_oil_radiators_kg:.2f} + {G_oil_tank_kg:.2f}',
            'result': G_oil_total_kg
        }

        # (233) Маса сталі бака, кг
        gamma_steel = 7650
        G_tank_steel_kg = 1.1 * (P_glb_m2 + 2 * P_kkr_m2) * gamma_steel * (3 / 1000)
        sec15['(233) Маса сталі бака, кг'] = {
            'formula': r'G_{ст.бака} = 1.1 (P_{гл.б} + 2 P_{к.кр}) \cdot \gamma_{ст} \cdot t_{ст}',
            'calculation': fr'G_{{ст.бака}} = 1.1 ({P_glb_m2:.2f} + 2 \cdot {P_kkr_m2:.2f}) \cdot {gamma_steel} \cdot \frac{{3}}{{1000}}',
            'calc_text': f'G_tank_steel = 1.1*({P_glb_m2:.2f}+2*{P_kkr_m2:.2f})*{gamma_steel}*(3/1000)',
            'result': G_tank_steel_kg
        }

        # (234) Маса сталі радіаторів, кг
        G_radiator_steel_kg = 686.125 * n_rad_adopted
        sec15['(234) Маса сталі радіаторів, кг'] = {
            'formula': r'G_{ст.рад} = C_{сталі} \cdot n_{рад}',
            'calculation': fr'G_{{ст.рад}} = 686.125 \cdot {n_rad_adopted}',
            'calc_text': f'G_radiator_steel = 686.125 * {n_rad_adopted}',
            'result': G_radiator_steel_kg
        }

        # (235) Повна маса трансформатора, кг
        G_transformer_total_kg = 1.2 * (G_tank_steel_kg + G_oil_total_kg + G_ach_kg + G_radiator_steel_kg)
        sec15['(235) Повна маса трансформатора, кг'] = {
            'formula': r'G_{заг} = 1.2 (G_{ст.бака} + G_{масла.заг} + G_{ач} + G_{ст.рад})',
            'calculation': fr'G_{{заг}} = 1.2 ({G_tank_steel_kg:.2f} + {G_oil_total_kg:.2f} + {G_ach_kg:.2f} + {G_radiator_steel_kg:.2f})',
            'calc_text': f'G_total = 1.2*({G_tank_steel_kg:.2f}+{G_oil_total_kg:.2f}+{G_ach_kg:.2f}+{G_radiator_steel_kg:.2f})',
            'result': G_transformer_total_kg
        }
        
        beta_range = np.linspace(0.7, 3.2, 100)
        x_range = beta_range**(1/4)
        Px_vals, i_nx_vals, J_vals, sigma_p_vals, d_core_vals, l_core_vals, L_ser_vals = [],[],[],[],[],[],[]
        for x_curr, beta_curr in zip(x_range, beta_range):
            G_u_curr = calculate_G_u(x_curr, A_coeff_calc, k_c_coeff_calc, k_ya)
            G_c_curr = calculate_G_c(x_curr, A1_kg_calc, A2_kg_calc)
            G_ya_curr = calculate_G_ya(x_curr, B1_kg_calc, B2_kg_calc)
            Px_vals.append(calculate_P_x(G_c_curr, G_u_curr, G_ya_curr, k_pd, p_c_loss, k_pu, p_ya_loss))
            Qx_curr = calculate_Q_x(G_c_curr, G_u_curr, G_ya_curr, k_td_prime, k_td_double_prime, q_c_reac, k_tu, k_tpl, q_ya_r)
            i_nx_vals.append(calculate_i_nx_percent(Qx_curr, S_n_VA))
            G_obm_curr = calculate_G_obm(x_curr, C1_kg_calc)
            J_vals.append(calculate_J(P_sc_kW, G_obm_curr, k_d_const))
            sigma_p_vals.append(calculate_sigma_p(x_curr, M_MPa_calc))
            d_core_curr = calculate_d_core(x_curr, A_coeff_calc)
            d_core_vals.append(d_core_curr)
            d12_curr = calculate_d12(d_core_curr, a_const_formula17)
            l_core_curr = calculate_l_core_winding(d12_curr, beta_curr)
            l_core_vals.append(l_core_curr)
            L_ser_curr = calculate_L_ser(d12_curr, A_01, A_02, b_const_formula19, d_core_curr)
            L_ser_vals.append(L_ser_curr)

        # Додаємо графіки 1-7 до першої секції
        sec1 = results['section_1']
        sec1['(Графік 1) Втрати НХ'] = _generate_plot_base64(plot_1_losses, beta_range, Px_vals)
        sec1['(Графік 2) Струм НХ'] = _generate_plot_base64(plot_2_current, beta_range, i_nx_vals)
        sec1['(Графік 3) Щільність струму'] = _generate_plot_base64(plot_3_density, beta_range, J_vals)
        sec1['(Графік 4) Механічні напруження'] = _generate_plot_base64(plot_4_stress, beta_range, sigma_p_vals)
        sec1['(Графік 5) Діаметр стрижня'] = _generate_plot_base64(plot_5_d_core, beta_range, d_core_vals)
        sec1['(Графік 6) Висота стрижня'] = _generate_plot_base64(plot_6_l_core, beta_range, l_core_vals)
        sec1['(Графік 7) Відстань між осями'] = _generate_plot_base64(plot_7_L_ser, beta_range, L_ser_vals)

        # =============================================================
        # === БЛОК 2: Розрахунок даних для графіків 8-11 (зовнішні характеристики) ===
        # =============================================================
        cos_phi2 = 0.98
        phi2_rad = math.acos(cos_phi2)
        beta_vals_plot = np.linspace(0, 1.25, 100)
        load_cases = {
            "R (Активне)": {"cos_phi": 1.0, "sin_phi": 0.0},
            "RL (Активно-індуктивне)": {"cos_phi": cos_phi2, "sin_phi": math.sin(phi2_rad)},
            "RC (Активно-ємнісне)": {"cos_phi": cos_phi2, "sin_phi": -math.sin(phi2_rad)}
        }
        results_U_V_plot = {case: [] for case in load_cases}
        results_U_percent_plot = {case: [] for case in load_cases}
        results_delta_U_V_plot = {case: [] for case in load_cases}
        results_delta_U_percent_plot = {case: [] for case in load_cases}
        
        for case, angles in load_cases.items():
            for beta in beta_vals_plot:
                delta_u1ka = beta * u1ka_V * angles["cos_phi"]
                delta_u1kp = beta * u1kp_V * angles["sin_phi"]
                delta_U1k_V = delta_u1ka + delta_u1kp
                U2_V = U_ph_HV_V - delta_U1k_V
                results_U_V_plot[case].append(U2_V)
                results_U_percent_plot[case].append(U2_V / U_ph_HV_V * 100)
                results_delta_U_V_plot[case].append(delta_U1k_V)
                results_delta_U_percent_plot[case].append(delta_U1k_V / U_ph_HV_V * 100)
                
        # Додаємо графіки 8-11 до дванадцятої секції
        sec12 = results['section_12']
        sec12['(Графік 8) Зовнішні характеристики, В'] = _generate_plot_base64(plot_8_external_V, beta_vals_plot, results_U_V_plot)
        sec12['(Графік 9) Зовнішні характеристики, %'] = _generate_plot_base64(plot_9_external_percent, beta_vals_plot, results_U_percent_plot)
        sec12['(Графік 10) Падіння напруги, В'] = _generate_plot_base64(plot_10_delta_U_V, beta_vals_plot, results_delta_U_V_plot)
        sec12['(Графік 11) Падіння напруги, %'] = _generate_plot_base64(plot_11_delta_U_percent, beta_vals_plot, results_delta_U_percent_plot)
        
        return results

    except Exception as e:
        import traceback
        print(traceback.format_exc())
        return {"error": f"Помилка під час розрахунку: {str(e)}"}
    
