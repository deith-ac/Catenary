# -*- coding: utf-8 -*-
import math
import matplotlib.pyplot as plt

# Precyzja Pythona dla liczba zmiennoprzecinkowych typu float to 53 bit
# w naszym przypadku niedokładności z tym związane są
# pomijalne. Nawet dla większych spadków i większych rozpiętości otrzymane
# wartości byłyby akceptowalne. Otrzymanie dokładniejszych wyników
# wiązałoby się z dużym nakładem pracy w odpowiednią implementację
# co z punktu widzenia konsturkcji budowlanych nie miałoby sensu.

# Dane projekotwe:
# spad SPAD_B [m]
SPAD_B = 39 - 19 
# rozpiętość SPAN_A [m]
SPAN_A = 674

# STAN 0 - dane:
# obciążenie ciężarem własnym q [N/m]
q_0 = 6.09
# temperatura zawieszenia przewodu [*c]
T_0 = 10
# założona wartość zwisu w środku przęsła [m]
f_0_m = 20

# STAN 1 - dane:
# Obciążenie wiatrem prostopadle do cięgna [N/m]
qw_1 = 0
# Ciężar oblodzenia [N/m]:
q_1k = 3.59
# Obciażenie całkowite bez wiatru [N/m]:
q_tw = q_1k + q_0
# Obciążenie całkowiete z wiatrem [N/m]:
q_t = math.sqrt((q_1k+q_0)**2 + qw_1**2)
# Temperatura T1 [*C]
T_1 = -5
# Różnica temperatur:
T_d = T_1 - T_0

# Does it even work? Apparently so!

# Charakterystyka pręta ALF-8 525mm^2 (520-AL1/67-ST1A)
# Pole przewodu [m^2]
# A = 5.869 * 10 ** (-4)
A = 1.34303 * 10 ** (-4)
# Średni moduł sprężystości [kPa]
# E = 7.6881581189 * 10 ** 10  # = 76.881 GPa
E = 10.6579 * 10 ** 10
# Średni współczynnik odkształcalności [1/K]
# e_T = 1.944 * 10 ** (-5)
e_T = 1.525 * 10 ** (-5)


def calculateM(k):
    """
    Wylicza wielkość mimośrodu m [m] na podstawie parametru
    lini zwisu k [m]

    :param k: parametr lini zwisu [m]
    :return: wartośc mimośrodu [m]
    """
    m = k * math.asinh(SPAD_B / (2 * k * math.sinh(SPAN_A / (2 * k))))
    return m

def calculateMprime(k):
    """
    Wyznacza odległość m' od środka cięgna do miejsca osiągnięcia ekstremum
    przez funkcję zwisu
    :param k:
    :return:
    """
    return k * math.asinh(SPAD_B / SPAN_A)


def calculateK(h, q):
    """
    Wylicza parametr linii zwisania k w [m]

    :param h: siła [N]
    :param q: wartość obciążenia jednorodnego [N/m]
    :return: parametr linii zwisania k w [m]
    """
    return h / q


def calculatefm(h, q):
    """
    wylicza wartość zwisu f.m dla podanej siły h [N], obciążenia q [N/m]
    w środku cięgna

    :param h: siła h w [N]
    :return: zwis fm w [m]
    """
    k = calculateK(h, q)
    m = calculateM(k)

    f_m = (SPAD_B / 2) + k * (
    math.cosh((-0.5 * SPAN_A + m) / k) - math.cosh(m / k))
    return f_m


def findh_s0(h_max, h_min, q):
    """
    Znajduje siłę naciągu metodą numeryczną (wykorzystana metoda bisekcji),
    należy podać granice górną i dolną dla metody bisekcji

    :param h_max: Górna granica dla szukania siły naciągu
    :param h_min: Dolna granica dla szukania siły naciągu
    :param q: całkowite obciążenie kabla [N/m]
    :return: h - siła naciągu, i - ilość potrzebnych iteracji
    """
    i = 1
    h = (h_min + h_max) / 2
    print("Wstępne H = " + str(h))
    f_m = calculatefm(h, q)

    while (math.fabs(f_m - f_0_m) >= 1 * 10 ** -8):
        if f_m < f_0_m:
            h_max = h
        else:
            h_min = h
        # print("iteracja #" + str(i) + " h_max = " + str(h_max) + " h_min = "
        # + str(h_min) + " nowe H: " + str(h) + " f_m = " + str(f_m)
        #      + " docelowe: " + str(f_0_m))
        h = (h_min + h_max) / 2
        f_m = calculatefm(h, q)
        i += 1

    return h, i


def testIfCorrect(h, q):
    f_0 = calculatefm(h, q)
    if math.fabs(f_0 - f_0_m) <= 1 * 10 ** -4:
        return True
    else:
        return False


def calculateInclination(x, h_s0, q):
    """
    Wyznacza nachylenie krzywej na podstawie współrzędnej X oraz
    siły naciągu ze stanu 0

    :param x: współrzędna X od lewej podpory
    :param h_s0: siła naciągu ze stanu 0
    :param q: całkowite obciążenie kabla [N/m]
    :return: nachylenie krzywej do płaszczyzny w danym miejscu
    """
    return math.sinh(x / calculateK(h_s0, q))


def calculateGeoLength(H, q):
    """
    Wylicza długość kabla wykorzystując warunki geometryczne dla danego H w [N]

    :param H: Siła naciągu kabla w [N]
    :param q: całkowite obciążenie kabla [N/m]
    :return: Zwraca geometryczną długość kabla w [m]
    """
    k = calculateK(H, q)
    if SPAD_B != 0:
        # length = SPAN_A * (1 + ((SPAN_A ** 2) * (q_0 ** 2)) / (24 * (H ** 2))
        #                + (SPAD_B ** 2) / (2 * (SPAN_A ** 2)))
        length = math.sqrt(SPAD_B**2 + (2 * k * math.sinh(SPAN_A / (2*k)))**2)
        return length
    elif SPAD_B == 0:
        length = 2 * k * math.sinh(SPAN_A / (2 * k))
        return length


def getBasePoints(m_0):
    """
    Zwraca wspolrzedne podpór.

    :param m_0: mimosród m
    :return: (xa, xb) wsp podpór
    """
    return math.fabs(-0.5 * SPAN_A + m_0), 0.5 * SPAN_A + m_0


def calculateForces(x, H, L, q):
    """
    Zwraca dwie siły w danym stanie, pierwszą jest składowa pionowa V_0(x)
    a drugą jest siła pozioma od wiatru W(x)

    :param x: miejsce dla którego wyznaczone zostaną siły [m]
    :param H: siła naciągu [N]
    :param q: obciążenie całkowite [N/m]
    :return: siła pionowa V, siła pozioma W (od wiatru)
    """
    k_0 = calculateK(H, q)
    #    Nachylenie:
    tg_A = math.sinh(x / k_0)
    #    Pionowa składowa siły naciągu V_0b dla stanu 0 w punkcie x:
    F_0 = H * tg_A  # w [N]
    sV = q_tw / q_t
    sH = qw_1 / q_t

    V_X = q_tw*k_0*math.sinh(x / k_0)
    W_X = 0.5*qw_1*L

    return round(V_X, 2), round(W_X, 2)


def findH_S1(h_max, h_min, l_0, h_0):
    """
    :param h_max: Górna granica poszukiwanej siły naciągu
    :param h_min: Dolna granica poszukiwanej siły naciągu
    :param l_0: Geometryczna długość kabla
    :param h_0: Siła naciągu ze stanu 0
    :param q: Całkowite obciążenie na 1 mb
    :return: Siła naciągu, ilość iteracji, długość kabli
     Siła naciągu, ilość iteracji, długość kabli
    """
    # print("$$$$$$$$$$$$ S1 CALCULATION START $$$$$$$$$$$$$$$$$")

    k_min = calculateK(h_min, q_t)
    l_geo_min = math.sqrt(
        SPAD_B ** 2 + (2 * k_min * math.sinh(SPAN_A / (2 * k_min))) ** 2)
    l_phys_min = l_0 * (1 + e_T * T_d + (h_min - h_0) / (E * A))
    # print("Wstępne min: k_min = %s, l_geo_min = %s, l_phys_min = %s"
    #      % (k_min, l_geo_min, l_phys_min))

    k_max = calculateK(h_max, q_t)
    l_geo_max = math.sqrt(
        SPAD_B ** 2 + (2 * k_max * math.sinh(SPAN_A / (2 * k_max))) ** 2)
    l_phys_max = l_0 * (1 + e_T * T_d + (h_max - h_0) / (E * A))
    # print("Wstępne max: k_max = %s, l_geo_max = %s, l_phys_max = %s"
    #      % (k_max, l_geo_max, l_phys_max))

    l_phys = 1
    l_geo = 5
    i = 1

    if l_geo_min >= l_phys_min and l_phys_max >= l_geo_max:
        while (math.fabs(l_phys - l_geo) >= 1 * 10 ** (-8)):
            h = (h_max + h_min) / 2
            k = calculateK(h, q_t)
            l_geo = math.sqrt(
                SPAD_B ** 2 + (2 * k * math.sinh(SPAN_A / (2 * k))) ** 2)
            l_phys = l_0 * (1 + e_T * T_d + (h - h_0) / (E * A))
            i += 1
            # print("Iteracja #%s - H = %s, l_geo = %s, l_phys = %s"
            #      % (i, h, l_geo, l_phys))
            if l_geo >= l_phys:
                h_min = h
            else:
                h_max = h
        # print("$$$$$$$$$$$$$$$$ CALCULATION FINISHED $$$$$$$$$$$$$$$")
        return h, i, l_geo
    else:
        print("Dla podanych wartości granicznych rozwiązanie nie istnieje")
        return False


def calculateHorizontalForce(qwind, span):
    """
    :param qwind: Obciążenie wiatrem prostopadłe do lini cięgna [N/m]
    :param span: Rozpiętość pozioma cięgna (długość w rzucie "z góry") [m]
    """
    return (qwind * span) / 2


def calculateTensileForce(Hten, H, V):
    """
    :param Hten: Wstępna siła naciągu [N]
    :param H: Pozioma siła na podporze [N]
    :param V: Pionowa siła na podporze [N]
    """
    return math.sqrt(Hten**2 + H**2 + V**2)


def drawCatenary(H, x_A, x_B):
    """
    :param H: Siła naciągu [N]
    :param x_A: Współrzędna X podpory A [m]
    :param x_B: Współrzędna X podpory B [m]
    :param q: Obciążenie całkowite [N/m]
    :return: void
    """
    k = calculateK(H, q_t)
    m = int(round(calculateM(k), 0))


    def catenary(x):
        # Krzywa łańcuchowa:
        f = -k * (math.cosh(x / k) - 1)
        return f

    def drange(start, stop, step):
        r = start
        while r < stop:
            yield r
            r += step

    plotRangeGenerator = drange((-SPAN_A / 2) + m, (SPAN_A / 2) + m, 2)
    plotRange = [i for i in plotRangeGenerator]
    Yaxis = [-catenary(x) for x in plotRange]
    Xaxis = [i for i in plotRange]

    plt.plot(Xaxis, Yaxis)
    plt.plot([Xaxis[0], Xaxis[-1]], [Yaxis[0], Yaxis[-1]], ".r")
    plt.annotate(s="A", xy=(Xaxis[0] - SPAN_A / 60, Yaxis[0] + 5), size=20)
    plt.annotate(s="B", xy=(Xaxis[-1] + SPAN_A / 60, Yaxis[-1] + 5), size=20)
    plt.axis([-SPAN_A / 2 * 1.2, SPAN_A / 2 * 2, -SPAD_B * 0.3, SPAD_B * 2])
    plt.ylabel("Y [m]")
    plt.xlabel("X [m] (rozpietosc)")

    plt.ylim([-40, 70])
    plt.show()


def main():
    # Stan 0
    """
    Wyświetla kluczowe wartości dla stanu 0 i stanu 1 w konsoli
    """
    H_S0, iters = findh_s0(5000000000, 1, q_0)
    x_A, x_B = getBasePoints(calculateM(calculateK(H_S0, q_0)))
    L_0 = calculateGeoLength(H_S0, q_0)

    print("Przyjęte dane geometryczne: \nRozpiętość: %s m\nSpad: %s m"
          "\nZałożona wartość zwisu: %s m" % (SPAN_A, SPAD_B, f_0_m))
    print("\n######### STAN 0 -> T0 [*C] = %s, q_0 = %s [N/m] ###########\n"
          % (T_0, q_0))
    print("Wyznaczona siła naciągu H dla S0: %s [N] w %s iteracji."
          % ((round(H_S0, 2)), iters))
    print("Poprawność obliczeń: %s" % testIfCorrect(H_S0, q_0))
    print("Długość kabla dla S0: %s m" % round(calculateGeoLength(H_S0, q_0), 1))
    print("Nachylenie w podporze B: %s"
          % round(calculateInclination(x_B, H_S0, q_0), 4))
    print("Parametr zwisu dla S0: %s" % round(calculateK(H_S0, q_0), 1))
    print("Mimośród przęsła: %s [m]"
          % round(calculateM(calculateK(H_S0, q_0)), 2))
    print("Współrzędne podpór: A = %s m, B = %s m"
          % (round(x_A, 2), round(x_B, 2)))
    VB, HB = calculateForces(x_B, H_S0, L_0, q_0)
    VA, HA = calculateForces(x_A, H_S0, L_0, q_0)
    print("Siła pionowa podpora A: %s N, W = %s N" % calculateForces(x_A, H_S0, L_0, q_0))
    print("Siła pionowa podpora B: %s N, W = %s N" % calculateForces(x_B, H_S0, L_0, q_0))
    print("Maksymalna siła naciągu B: %s N" % math.sqrt(VB**2 + H_S0**2))
    print("Maksymalna siła naciągu A: %s N" % math.sqrt(VA**2 + H_S0**2))

    #    Stan 1
    print("\n######### STAN 1 -> T1 = %s [*C], T_d = %s [*C],"
          " q_t = %s [N/m] ###########\n" % (T_1, T_d, q_t))

    L_0 = calculateGeoLength(H_S0, q_0)
    # print("Długość L_0 stanu 0: %s m" % L_0)
    H_S1, iters_S1, cableLength = findH_S1(1400000, 10000, L_0, H_S0)
    L_1 = calculateGeoLength(H_S1, q_t)
    x_A, x_B = getBasePoints(calculateM(calculateK(H_S1, q_t)))
    V_A, W_A = calculateForces(x_A, H_S1, L_1, q_tw)
    V_B, W_B = calculateForces(x_B, H_S1, L_1, q_tw)

    if H_S1 is not False:
        print("Wyznaczona siła naciągu H dla S1: %s [N] w %s iteracji."
              % (round(H_S1, 2), iters_S1))
        print("Długość kabla dla S1: %s m" % round(cableLength, 1))
        print("Mimośród przęsła: %s [m]"
              % round(calculateM(calculateK(H_S1, q_t)), 2))
        print("Zwis przewodu w połowie rozpiętości przęsła dla stanu 1: %s [m]"
              % round(calculatefm(H_S1, q_t), 2))
        x_A, x_B = getBasePoints(calculateM(calculateK(H_S1, q_t)))
        print("Współrzędne podpór: A = %s m, B = %s m"
              % (round(x_A, 2), round(x_B, 2)))
        print("Podpora A:")
        print("Siła pionowa: %s N"
              % (V_A))
        print("Siła pozioma: %s N" % W_A)
        print("Maksymalna siła naciągu: %s N"
              % calculateTensileForce(H_S1, W_A, V_A))
        print("Podpora B:")
        print("Siła pionowa: %s N"
              % (V_B))
        print("Siła pozioma: %s N" % W_B)
        print("Maksymalna siła naciągu: %s N"
              % calculateTensileForce(H_S1, W_B, V_B))

        drawCatenary(H_S1, x_A, x_B)


if __name__ == "__main__":
    main()
