E_P111 = 599821.802
E_P112 = 599772.346
E_NMB1 = 599809.075


R = 6380000.0
m0 = 0.9996


E_middel = (
    E_P111 +
    E_P112 +
    E_NMB1
) / 3


y_m = E_middel - 500000.0


m_k = m0 + (y_m**2) / (2 * R**2)


print("BEREGNING AV MÅLESTOKKFAKTOR I UTM")


print(f"E middel      = {E_middel:.3f} m")
print(f"y_m           = {y_m:.3f} m")
print(f"m0            = {m0:.7f}")
print(f"målestokk m_k = {m_k:.9f}")
