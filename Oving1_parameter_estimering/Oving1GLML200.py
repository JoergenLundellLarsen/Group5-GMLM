from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path

import numpy as np
from math import sqrt, atan2, pi


def beregn_parameterestimering():
    np.set_printoptions(
        precision=6,
        suppress=True,
        linewidth=160
    )


    punkter = {
        "P111": [6615743.203, 599821.802, True],
        "NMB1": [6615708.451, 599809.075, True],


        "P112": [6615726.406, 599772.346, False],
    }


    orienteringer = {
        "P111": 0,
        "P112": 1,
        "NMB1": 2,
    }


    ukjente = {
        "P112": [3, 4],
    }

    addisjonskonstant_col = 5
    antall_ukjente = 6


    retning_obs = [


        ["P111", "P112", 43.1038],
        ["P111", "P112", 43.1014],
        ["P111", "P112", 43.1005],
        ["P111", "P112", 43.0987],


        ["P111", "NMB1", 386.3027],
        ["P111", "NMB1", 386.3060],
        ["P111", "NMB1", 386.3036],
        ["P111", "NMB1", 386.3041],


        ["P112", "P111", 354.8015],
        ["P112", "P111", 354.8037],
        ["P112", "P111", 354.8040],


        ["P112", "NMB1", 4.5993],

        ["P112", "NMB1", 4.5995],
        ["P112", "NMB1", 4.6014],


        ["NMB1", "P112", 171.4201],
        ["NMB1", "P112", 171.4202],
        ["NMB1", "P112", 171.4248],


        ["NMB1", "P111", 264.8323],
        ["NMB1", "P111", 264.8329],
        ["NMB1", "P111", 264.8347],
        ["NMB1", "P111", 264.8340],
    ]


    avstand_obs = [


        ["P111", "P112", 52.3410, 103.0064],
        ["P111", "P112", 52.3405, 103.0069],
        ["P111", "P112", 52.3405, 103.0078],
        ["P111", "P112", 52.3405, 103.0063],


        ["P111", "NMB1", 37.2615, 106.6486],
        ["P111", "NMB1", 37.2615, 106.6486],
        ["P111", "NMB1", 37.2615, 106.6478],
        ["P111", "NMB1", 37.2615, 106.6468],


        ["P112", "P111", 52.3405, 96.9925],
        ["P112", "P111", 52.3405, 96.9956],
        ["P112", "P111", 52.3405, 96.9961],
        ["P112", "P111", 52.3405, 96.9940],


        ["P112", "NMB1", 40.9480, 102.1976],
        ["P112", "NMB1", 40.9480, 102.1883],
        ["P112", "NMB1", 40.9475, 102.1974],
        ["P112", "NMB1", 40.9475, 102.1893],


        ["NMB1", "P112", 40.9480, 97.8054],
        ["NMB1", "P112", 40.9480, 97.8077],
        ["NMB1", "P112", 40.9485, 97.8072],
        ["NMB1", "P112", 40.9485, 97.8062],


        ["NMB1", "P111", 37.2610, 93.3555],
        ["NMB1", "P111", 37.2610, 93.3565],
        ["NMB1", "P111", 37.2610, 93.3572],
        ["NMB1", "P111", 37.2615, 93.3567],
    ]


    sigma_retning = 0.0005

    def sigma_avstand(D):

        return 0.003 + 0.003 * (D / 1000)


    def normaliser_gon(v):

        while v > 200:
            v -= 400
        while v < -200:
            v += 400
        return v


    def beregn_retning(fra_punkt, til_punkt):
        XA, YA = punkter[fra_punkt][0], punkter[fra_punkt][1]
        XB, YB = punkter[til_punkt][0], punkter[til_punkt][1]

        dX = XB - XA
        dY = YB - YA

        r = atan2(dY, dX) * 200 / pi

        if r < 0:
            r += 400

        return r


    def beregn_avstand(fra_punkt, til_punkt):
        XA, YA = punkter[fra_punkt][0], punkter[fra_punkt][1]
        XB, YB = punkter[til_punkt][0], punkter[til_punkt][1]

        dX = XB - XA
        dY = YB - YA

        return sqrt(dX**2 + dY**2)


    orientering_start = {}

    for fra_punkt in orienteringer:
        for obs in retning_obs:
            if obs[0] == fra_punkt:
                til_punkt = obs[1]
                obs_retning = obs[2]

                beregnet = beregn_retning(fra_punkt, til_punkt)

                orientering_start[fra_punkt] = normaliser_gon(
                    beregnet - obs_retning
                )
                break


    def lag_retning_obs(fra_punkt, til_punkt, obs_retning):
        XA, YA = punkter[fra_punkt][0], punkter[fra_punkt][1]
        XB, YB = punkter[til_punkt][0], punkter[til_punkt][1]

        dX = XB - XA
        dY = YB - YA
        D = sqrt(dX**2 + dY**2)

        beregnet = beregn_retning(fra_punkt, til_punkt)
        orientering = orientering_start[fra_punkt]


        l = normaliser_gon(obs_retning + orientering - beregnet)

        A_rad = [0] * antall_ukjente


        A_rad[orienteringer[fra_punkt]] = -1


        if til_punkt in ukjente:
            colX, colY = ukjente[til_punkt]

            A_rad[colX] = -(dY / D**2) * 200 / pi
            A_rad[colY] =  (dX / D**2) * 200 / pi


        if fra_punkt in ukjente:
            colX, colY = ukjente[fra_punkt]

            A_rad[colX] =  (dY / D**2) * 200 / pi
            A_rad[colY] = -(dX / D**2) * 200 / pi


        p = 1

        return A_rad, l, p

    R = 6380000.0
    H = 104.044
    m_k = 0.9997224
    def lag_avstand_obs(fra_punkt, til_punkt, D_skraa, Z_middel):

        XA, YA = punkter[fra_punkt][0], punkter[fra_punkt][1]
        XB, YB = punkter[til_punkt][0], punkter[til_punkt][1]

        dX = XB - XA
        dY = YB - YA

        D = sqrt(dX**2 + dY**2)


        Z_rad = Z_middel * pi / 200


        D_obs_kart = (
            D_skraa
            * np.sin(Z_rad)
            * (1 - H / R)
            * m_k
        )


        l = D_obs_kart - D

        A_rad = [0] * antall_ukjente


        if til_punkt in ukjente:
            colX, colY = ukjente[til_punkt]

            A_rad[colX] = dX / D
            A_rad[colY] = dY / D


        if fra_punkt in ukjente:
            colX, colY = ukjente[fra_punkt]

            A_rad[colX] = -dX / D
            A_rad[colY] = -dY / D


        A_rad[addisjonskonstant_col] = (
            -np.sin(Z_rad)
            * (1 - H / R)
            * m_k
        )

        sigma = sigma_avstand(D)

        p = sigma_retning**2 / sigma**2

        return A_rad, l, p


    A_liste = []
    l_liste = []
    p_liste = []
    navn_liste = []

    for fra_punkt, til_punkt, obs_retning in retning_obs:
        A_rad, l_verdi, p_verdi = lag_retning_obs(
            fra_punkt,
            til_punkt,
            obs_retning
        )

        A_liste.append(A_rad)
        l_liste.append(l_verdi)
        p_liste.append(p_verdi)
        navn_liste.append(f"Retning {fra_punkt}->{til_punkt}")

    for fra_punkt, til_punkt, D_skraa, Z_middel in avstand_obs:
        A_rad, l_verdi, p_verdi = lag_avstand_obs(
            fra_punkt,
            til_punkt,
            D_skraa,
            Z_middel
        )

        A_liste.append(A_rad)
        l_liste.append(l_verdi)
        p_liste.append(p_verdi)
        navn_liste.append(f"Avstand {fra_punkt}->{til_punkt}")

    A = np.array(A_liste, dtype=float)
    l = np.array(l_liste, dtype=float).reshape(-1, 1)
    P = np.diag(p_liste)


    N = A.T @ P @ A
    h = A.T @ P @ l

    dx = np.linalg.solve(N, h)

    v = A @ dx - l

    frihetsgrader = len(l) - antall_ukjente

    s0_hat = sqrt((v.T @ P @ v)[0, 0] / frihetsgrader)

    Qxx = np.linalg.inv(N)

    std_parametre = s0_hat * np.sqrt(np.diag(Qxx)).reshape(-1, 1)


    colN, colE = ukjente["P112"]
    colA = addisjonskonstant_col


    N_start = punkter["P112"][0]
    E_start = punkter["P112"][1]

    dN = dx[colN, 0]
    dE = dx[colE, 0]

    N_ny = N_start + dN
    E_ny = E_start + dE


    a = dx[colA, 0]


    std_N = std_parametre[colN, 0]
    std_E = std_parametre[colE, 0]
    std_a = std_parametre[colA, 0]


    print("ESTIMERING AV ADDISJONSKONSTANT")
    print(f"Observasjoner: {len(l)} ({len(retning_obs)} retninger, {len(avstand_obs)} avstander)")
    print(f"Ukjente: {antall_ukjente}")
    print(f"Frihetsgrader: {frihetsgrader}")

    print("\nAddisjonskonstant")
    print(f"a = {a * 1000:.2f} mm")
    print(f"Standardavvik: {std_a * 1000:.2f} mm")
    print(f"D_sann = D_obs + ({a:.6f} m)")

    print("\nEstimerte koordinater for P112")
    print(f"Nord: {N_ny:.4f} m (standardavvik {std_N * 1000:.2f} mm)")
    print(f"Øst:  {E_ny:.4f} m (standardavvik {std_E * 1000:.2f} mm)")

    print("\nKontroll av utjevningen")
    print(f"Estimert standardavvik for vektenheten: {s0_hat:.8f} gon")
    print(f"Forhold til antatt standardavvik: {s0_hat / sigma_retning:.2f}")
    print(f"Største absolutte retningsresidual: {np.max(np.abs(v[:len(retning_obs), 0])):.6f} gon")
    print(f"Største absolutte avstandsresidual: {np.max(np.abs(v[len(retning_obs):, 0])) * 1000:.2f} mm")


def main():
    resultat_buffer = StringIO()
    with redirect_stdout(resultat_buffer):
        beregn_parameterestimering()

    resultat = resultat_buffer.getvalue()
    output_fil = Path(__file__).with_name("Oving1GLML200_output.txt")
    output_fil.write_text(resultat, encoding="utf-8")

    print(resultat, end="")


if __name__ == "__main__":
    main()
