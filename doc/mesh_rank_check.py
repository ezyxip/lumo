"""
Проверка ранга якобиана MESH-системы (N=3, C=2: total condenser + feed + kettle reboiler)
для варианта A (Q_W задан) и варианта B (Q_W — неизвестная), плюс проверка обусловленности.

Требует: sympy, numpy (pip install sympy numpy)
"""

import sympy as sp
import numpy as np

# ---------- generic smooth stand-in thermo functions (не настоящий Antoine, просто гладкие функции с реалистичным порядком нелинейности) ----------
a1, a2, b1, b2, c1, c2 = 8.1, 7.9, 1500, 1600, 220, 210
p1, p2, q1, q2 = 0.09, 0.11, 5.0, 6.0
P1, P2, Q1, Q2 = 0.03, 0.05, 30.0, 28.0
Dv, Fv, z1v, z2v, Rv = 40.0, 100.0, 0.4, 0.6, 3.0

def Ksym(T, comp):
    a, b, c = (a1, b1, c1) if comp == 1 else (a2, b2, c2)
    return sp.exp(a - b/(T + c))

def hsym(T, x1, x2):
    return x1*(p1*T + q1) + x2*(p2*T + q2)

def Hsym(T, y1, y2):
    return y1*(P1*T + Q1) + y2*(P2*T + Q2)

# ---------- unknowns ----------
x1_1, x2_1, y1_1, y2_1, T1, L1, QD = sp.symbols('x1_1 x2_1 y1_1 y2_1 T1 L1 QD')
x1_2, x2_2, y1_2, y2_2, T2, L2, V2 = sp.symbols('x1_2 x2_2 y1_2 y2_2 T2 L2 V2')
x1_3, x2_3, y1_3, y2_3, T3, W, V3, QW = sp.symbols('x1_3 x2_3 y1_3 y2_3 T3 W V3 QW')
D, F, z1, z2, R = sp.symbols('D F z1 z2 R')  # параметры (фиксированы числом ниже)

QW_fixed_value = 55.0  # для варианта A ("задан")

def build_equations(QW_is_unknown):
    eqs = []
    # Stage 1: total condenser
    eqs.append(V2*y1_2 - (L1 + D)*x1_1)
    eqs.append(V2*y2_2 - (L1 + D)*x2_1)
    eqs.append(y1_1 - Ksym(T1, 1)*x1_1)
    eqs.append(y2_1 - Ksym(T1, 2)*x2_1)
    eqs.append(x1_1 + x2_1 - 1)
    eqs.append(y1_1 + y2_1 - 1)
    eqs.append(V2*Hsym(T2, y1_2, y2_2) - (L1 + D)*hsym(T1, x1_1, x2_1) - QD)

    # Stage 2: feed
    eqs.append(L1*x1_1 + V3*y1_3 + F*z1 - L2*x1_2 - V2*y1_2)
    eqs.append(L1*x2_1 + V3*y2_3 + F*z2 - L2*x2_2 - V2*y2_2)
    eqs.append(y1_2 - Ksym(T2, 1)*x1_2)
    eqs.append(y2_2 - Ksym(T2, 2)*x2_2)
    eqs.append(x1_2 + x2_2 - 1)
    eqs.append(y1_2 + y2_2 - 1)
    hF = hsym(sp.Float(350.0), z1, z2)
    eqs.append(L1*hsym(T1, x1_1, x2_1) + V3*Hsym(T3, y1_3, y2_3) + F*hF
               - L2*hsym(T2, x1_2, x2_2) - V2*Hsym(T2, y1_2, y2_2))

    # Stage 3: kettle reboiler
    eqs.append(L2*x1_2 - W*x1_3 - V3*y1_3)
    eqs.append(L2*x2_2 - W*x2_3 - V3*y2_3)
    eqs.append(y1_3 - Ksym(T3, 1)*x1_3)
    eqs.append(y2_3 - Ksym(T3, 2)*x2_3)
    eqs.append(x1_3 + x2_3 - 1)
    eqs.append(y1_3 + y2_3 - 1)
    QW_term = QW if QW_is_unknown else QW_fixed_value
    eqs.append(L2*hsym(T2, x1_2, x2_2) - W*hsym(T3, x1_3, x2_3) - V3*Hsym(T3, y1_3, y2_3) + QW_term)

    # Closure (Rating): L1 - R*D = 0
    eqs.append(L1 - R*D)
    return eqs

eqs_A = build_equations(QW_is_unknown=False)
unknowns_A = [x1_1, x2_1, y1_1, y2_1, T1, L1, QD,
              x1_2, x2_2, y1_2, y2_2, T2, L2, V2,
              x1_3, x2_3, y1_3, y2_3, T3, W, V3]          # 21 unknowns

eqs_B = build_equations(QW_is_unknown=True)
unknowns_B = unknowns_A + [QW]                            # 22 unknowns

print("Вариант A: #уравнений =", len(eqs_A), " #неизвестных =", len(unknowns_A))
print("Вариант B: #уравнений =", len(eqs_B), " #неизвестных =", len(unknowns_B))

J_A = sp.Matrix(eqs_A).jacobian(unknowns_A)   # 22 x 21
J_B = sp.Matrix(eqs_B).jacobian(unknowns_B)   # 22 x 22

params_subs = {D: Dv, F: Fv, z1: z1v, z2: z2v, R: Rv}
rng = np.random.default_rng(12345)

def random_point():
    """Случайная точка, удовлетворяющая S (Σx=1, Σy=1) на каждой ступени."""
    vals = {}
    for xx1, xx2, yy1, yy2, TT in [
        (x1_1, x2_1, y1_1, y2_1, T1),
        (x1_2, x2_2, y1_2, y2_2, T2),
        (x1_3, x2_3, y1_3, y2_3, T3),
    ]:
        fx = rng.uniform(0.1, 0.9)
        fy = rng.uniform(0.1, 0.9)
        vals[xx1], vals[xx2] = fx, 1 - fx
        vals[yy1], vals[yy2] = fy, 1 - fy
        vals[TT] = rng.uniform(300, 400)
    vals[L1] = rng.uniform(50, 150)
    vals[L2] = rng.uniform(50, 150)
    vals[V2] = rng.uniform(50, 150)
    vals[V3] = rng.uniform(50, 150)
    vals[W]  = rng.uniform(20, 100)
    vals[QD] = rng.uniform(-50, 50)
    vals[QW] = rng.uniform(-50, 50)
    return vals

# ---------- 1) базовая проверка ранга ----------
print("\n--- ранг якобиана в случайных точках (S=0) ---")
for trial in range(3):
    pt = random_point()
    full_subs = {**params_subs, **pt}
    JA_num = np.array(J_A.subs(full_subs)).astype(np.float64)
    JB_num = np.array(J_B.subs(full_subs)).astype(np.float64)
    rankA = np.linalg.matrix_rank(JA_num, tol=1e-8)
    rankB = np.linalg.matrix_rank(JB_num, tol=1e-8)
    print(f"trial {trial}: shape A={JA_num.shape}, rank A={rankA}  |  "
          f"shape B={JB_num.shape}, rank B={rankB}")

# ---------- 2) обусловленность варианта B ----------
print("\n--- обусловленность якобиана варианта B (22x22), 8 точек ---")
for trial in range(8):
    pt = random_point()
    full_subs = {**params_subs, **pt}
    JB_num = np.array(J_B.subs(full_subs)).astype(np.float64)
    sv = np.linalg.svd(JB_num, compute_uv=False)
    print(f"trial {trial}: min sv={sv[-1]:.6g}, max sv={sv[0]:.6g}, "
          f"cond={sv[0]/sv[-1]:.3g}")

# ---------- 3) какая строка "лишняя" в варианте A ----------
print("\n--- вариант A: доминирующие строки в линейной зависимости ---")
pt = random_point()
full_subs = {**params_subs, **pt}
JA_num = np.array(J_A.subs(full_subs)).astype(np.float64)
u, s, vt = np.linalg.svd(JA_num.T)
null_c = u[:, -1]
labels = ["M1_s1","M2_s1","E1_s1","E2_s1","Sx_s1","Sy_s1","H_s1",
          "M1_s2","M2_s2","E1_s2","E2_s2","Sx_s2","Sy_s2","H_s2",
          "M1_s3","M2_s3","E1_s3","E2_s3","Sx_s3","Sy_s3","H_s3",
          "closure"]
order = np.argsort(-np.abs(null_c))
print("smallest singular value of J_A^T (≈0, подтверждает зависимость):", s[-1])
for i in order[:8]:
    print(f"  {labels[i]:8s}  coeff={null_c[i]:+.4f}")

# ---------- 4) физическое масштабирование ----------
print("\n--- проверка: обусловленность после физического scaling ---")
scales = []
for u_ in unknowns_B:
    name = str(u_)
    if name.startswith('x') or name.startswith('y'):
        scales.append(1.0)
    elif name.startswith('T'):
        scales.append(350.0)
    elif name.startswith('L') or name.startswith('V') or name == 'W':
        scales.append(100.0)
    elif name.startswith('Q'):
        scales.append(50.0)
    else:
        scales.append(1.0)
scales = np.array(scales)
eq_scales = np.array([100,100,1,1,1,1,3.5e4]*3 + [100])  # M~100,E~1,S~1,H~3.5e4 на ступень, closure~100

for trial in range(4):
    pt = random_point()
    full_subs = {**params_subs, **pt}
    JB_num = np.array(J_B.subs(full_subs)).astype(np.float64)
    sv_raw = np.linalg.svd(JB_num, compute_uv=False)
    JB_col = JB_num * scales[np.newaxis, :]
    sv_col = np.linalg.svd(JB_col, compute_uv=False)
    JB_rc = (JB_num / eq_scales[:, np.newaxis]) * scales[np.newaxis, :]
    sv_rc = np.linalg.svd(JB_rc, compute_uv=False)
    print(f"trial {trial}: cond(raw)={sv_raw[0]/sv_raw[-1]:.3g}  "
          f"cond(col-scaled)={sv_col[0]/sv_col[-1]:.3g}  "
          f"cond(row+col-scaled)={sv_rc[0]/sv_rc[-1]:.3g}")
