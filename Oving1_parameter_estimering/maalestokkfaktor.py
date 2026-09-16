# ============================================================
# BEREGNING AV MÅLESTOKKFAKTOR I UTM
# ============================================================
#
# Vi bruker formelen for målestokkfaktor i UTM:
#
#                    y_m^2
#   m_k = m_0 + -------------
#                   2 * R^2
#
# der:
#   m_k = målestokkfaktor i måleområdet
#   m_0 = målestokkfaktor langs UTM-sonens sentralmeridian
#         m_0 = 0.9996 for UTM
#   y_m = midlere avstand fra sentralmeridianen / nordaksen
#   R   = jordens radius, her satt til 6 380 000 m
#
#
# I UTM har sentralmeridianen en falsk østverdi på 500 000 m.
# Derfor beregnes y_m fra:
#
#   y_m = E_middel - 500 000
#
#
# Siden måleområdet er lite bruker vi én gjennomsnittlig
# målestokkfaktor for hele området.
#
# E_middel beregnes som gjennomsnittet av Easting-koordinatene
# til P111, P112 og NMB1:
#
#              E_P111 + E_P112 + E_NMB1
#   E_middel = --------------------------
#                          3
#
# ============================================================


# Easting-koordinater for punktene
E_P111 = 599821.802
E_P112 = 599772.346
E_NMB1 = 599809.075


# Konstantverdier
R = 6380000.0     # jordradius i meter
m0 = 0.9996       # målestokkfaktor på sentralmeridianen i UTM


# ============================================================
# 1. BEREGN MIDDEL EASTING
# ============================================================

E_middel = (
    E_P111 +
    E_P112 +
    E_NMB1
) / 3


# ============================================================
# 2. BEREGN y_m
# ============================================================
#
# UTM bruker falsk østverdi E = 500 000 m på sentralmeridianen.
#
# Formel:
#
#   y_m = E_middel - 500 000
#

y_m = E_middel - 500000.0


# ============================================================
# 3. BEREGN MÅLESTOKKFAKTOR
# ============================================================
#
# Formel:
#
#                    y_m^2
#   m_k = m_0 + -------------
#                   2 * R^2
#

m_k = m0 + (y_m**2) / (2 * R**2)


# ============================================================
# RESULTAT
# ============================================================

print("BEREGNING AV MÅLESTOKKFAKTOR I UTM")


print(f"E middel      = {E_middel:.3f} m")
print(f"y_m           = {y_m:.3f} m")
print(f"m0            = {m0:.7f}")
print(f"målestokk m_k = {m_k:.9f}")
