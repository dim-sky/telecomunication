from re import I
import math, sys
from typing import Self
import numpy as np
import matplotlib.pyplot as plt




# Σκυλάκος Δημήτριος it21891

# Υπολογισμός Μέσης ισχύς του σήματος
def calculate_E_a_square(input_bit_symbols_array):
    # sigma__square 
    average = sum(input_bit_symbols_array) / len(input_bit_symbols_array)
    temp = np.array(input_bit_symbols_array) - average
    sigma__square = temp**2

    return sigma__square + average

def P_function(t, Ts):
    if 0 < t < Ts/2: 
        return (2/Ts)*t
    elif Ts/2 <= t <= Ts:
        return (2/Ts)*(Ts-t)
    else: 
        return 0



def canlculateXt(input_bit_symbols_array, t, Ts):
    # sinngle_palm = []
    x_t = []
    # print(Ts)
    for k , symbol in enumerate(input_bit_symbols_array):
        # print(k , symbol)
        for single_t in t:
            x_t.append(symbol*P_function(single_t, Ts))
        # x_t.append(sinngle_palm)
        # sinngle_palm = []
    # print(x_t)
    return x_t


def createSymbolArrayFromBits(input_bits, symbols, m):

    symbols_array = []

    # devide string into list of substrings 
    list_of_bits = []
    start = 0
    end = int(m)
    step = int(m)
    while (end <= len(input_bits)):
        list_of_bits.append(input_bits[start:end])
        start = end
        end += step

    # Για κάθε στοιχείο του πίνακα παίρνω τον δεκαδικό αριθμό που αντιστοιχεί σε αυτό, το οποίο μου δηλώνει την θέση του συμβόλου στον πίνακα symbols
    for el in list_of_bits:
        index_of_symbol = int(el, 2)
        symbols_array.append(symbols[index_of_symbol])


    return symbols_array


def createSymbols(m, M):
    b = 1
    symbols = []
    for i in range(int(M)):
        symbol = (2*i - M + 1)*b
        symbols.append(symbol)

    return symbols


def main():
    input_bits='110101100011100010010011'
    M = 16 # 2^4
    m = math.log(M, 2)
    Ts = m*1e-09  # Ρυθμός μετάδοσης = 1 / διάρκεια συμβόλου ==> 1 / 10^9 αν κάθε σύμβολο αναπαρηστούσε 1 bit | επειδή κάθε σύμβολο αναπαρηστά m bits, πολλαπλασιάζω με m
    N = 1000 # Αριθμός δειγμάτων ανά περίοδο


    symbols = createSymbols(m,M)
    print(symbols)
    input_bit_symbols_array = createSymbolArrayFromBits(input_bits, symbols, m)
    print(input_bit_symbols_array)


    # Άξονας του χρόνου για Τs
    t = np.linspace(0, Ts, N )
        

    # κυματομορφή
    x_t = canlculateXt(input_bit_symbols_array, t, Ts)
    # print(x_t[1])

    # Άξονας του χρόνου για k*Τs
    k = len(input_bit_symbols_array)
    t_n = np.linspace(0, Ts, k*N )
    
    
    plt.figure(1)
    plt.plot(t_n, x_t, label="Κυματομορφή Τριγωνικού Παλμού")
    plt.xlabel("Χρόνος (ns)")
    plt.ylabel("Πλάτος")
    plt.title("Κυματομορφή Τριγωνικού Παλμού με Κώδικα Gray")
    plt.grid()
    plt.legend()
   



 #############################################
 ##                                         ##
 ##                 Μέρος Β                 ##
 ##                                         ##
 #############################################

    Nf = 100
    zero_padding = 50

    # θέλω να υπολογίσω ένα μόνο τρίγωνο"-παλμό στο πεδίο του χρόνου.  
    t2 = np.linspace(0, Ts, Nf) # ο άξωνας του χρόνου
    x_t_2 = [P_function(t, Ts) for t in t2] # Οι τιμές του σήματος μου

    # Extend the time axis to include zero padding
    total_points = Nf + 2 * zero_padding
    t2_padded = np.linspace(-zero_padding * (Ts / Nf), Ts + zero_padding * (Ts / Nf), total_points)
    # Add zeros before and after the signal
    x_t_2_padded = np.pad(x_t_2, (zero_padding, zero_padding), mode='constant', constant_values=0)



    # x_t_2 = np.pad(x_t_2, (100, 100), mode='constant', constant_values=0
    plt.figure(4)
    plt.plot(t2_padded,x_t_2_padded, label="στο πεδίο του χρόνου")
    plt.xlabel("t")
    plt.ylabel("")
    plt.title("Το σήμα μου στο πεδίο του χρόνου")
    plt.grid()
    plt.legend()


    # Συχνότητα 
    Fmax = 1/Ts
    
    
    # freq_axis = np.linspace(-Fmax/2, Fmax/2, Nf )
    Nf = len(t2_padded)
    n = np.arange(-Nf/2, Nf/2, 1 )

    Dt = 2*Ts/ Nf
    Df = 1.0 / (Nf * Dt)

    freq_axis = n * Df
    
    # freq_axis = np.linspace(-0.1,0.1, 50)

    # Μαθηματικός Υπολογισμός P(f)
    P_f = (Ts) * (np.sinc(np.pi * freq_axis * Ts))**2
    print(P_f)
    
    
   
    fft_x_t_2 = Dt * np.fft.fftshift(np.fft.fft( np.fft.fftshift(x_t_2_padded)))
    
    
    print(fft_x_t_2.astype(float))

    plt.figure(3)
    plt.plot(t2_padded,P_f, label="Μαθηματικός Υπολογισμός P(f)")
    plt.xlabel("Συχνότητα (Hz)")
    plt.ylabel("")
    plt.title("Το σήμα μου στο πεδίο των συχνοτήτων")
    plt.grid()
    plt.legend()


    # E_a_square = calculate_E_a_square(input_bit_symbols_array)


    plt.figure(2)
    plt.plot(freq_axis,P_f,freq_axis,fft_x_t_2, 'o')
    plt.xlabel("Συχνότητα (Hz)")
    plt.ylabel("")
    plt.title("Το σήμα μου στο πεδίο των συχνοτήτων")
    plt.grid()
    plt.legend()

    plt.show()

if __name__ == "__main__":
    main()




