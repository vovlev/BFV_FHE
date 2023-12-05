from BFV_tools import *

import tkinter as tk

from tkinter import ttk
from tkinter import *

def open_popup(win, text):
    top = Toplevel(win)
    top.geometry("1260x400")
    top.title("Результаты расчетов")
    Label(top, text = text, font=('Mistral 12')).pack(anchor="nw")
    
# Непосредственно калькулятор
def calculate_fhe():
    t = int(t_tf.get())
    n = int(n_tf.get())
    logq = int(logq_tf.get())
    n1 = int(n1_tf.get())
    n2 = int(n2_tf.get())
    argx=str(arx1.get())
    



    # Enter proper parameters below
    #t, n, logq = 16, 1024, 27
    #t, n, logq = 256, 2048, 37
    # t, n, logq = 1024, 4096, 58

    # other necessary parameters (based on n and log(q) determine other parameter)
    q,psi,psiv,w,wv = ParamGen(n,logq) 

    # Determine mu, sigma (for discrete gaussian distribution)
    mu    = 0
    sigma = 0.5 * 3.2

    # Determine T, p (for relinearization and galois keys) based on noise analysis 
    T = 256
    p = q**3 + 1

    # Generate polynomial arithmetic tables
    w_table    = [1]*n
    wv_table   = [1]*n
    psi_table  = [1]*n
    psiv_table = [1]*n
    for i in range(1,n):
        w_table[i]    = ((w_table[i-1]   *w)    % q)
        wv_table[i]   = ((wv_table[i-1]  *wv)   % q)
        psi_table[i]  = ((psi_table[i-1] *psi)  % q)
        psiv_table[i] = ((psiv_table[i-1]*psiv) % q)

    qnp = [w_table,wv_table,psi_table,psiv_table]

    # Generate BFV evaluator
    Evaluator = BFV(n, q, t, mu, sigma, qnp)

    # Generate Keys
    Evaluator.SecretKeyGen()
    Evaluator.PublicKeyGen()
    Evaluator.EvalKeyGenV1(T)
    Evaluator.EvalKeyGenV2(p)

    # print system parameters
    #print(Evaluator)

    op_text="Введенные данные:\n"
    op_text+="n1: {}\n".format(n1)
    op_text+="n2: {}\n".format(n2)
    add_text=op_text
    sub_text=op_text
    mul_text=op_text
    add_text+="n1+n2: {}\n".format(n1+n2)
    sub_text+="n1-n2: {}\n".format(n1-n2)
    mul_text+="n1*n2: {}\n".format(n1*n2)

    # Encode random messages into plaintext polynomials
    op_text1="Данные n1 и n2 закодированы в полиномы m1(x) и m2(x):\n"
    m1 = Evaluator.IntEncode(n1)
    m2 = Evaluator.IntEncode(n2)

    op_text1+="m1(x): {}\n".format(m1)
    op_text1+="m2(x): {}\n".format(m2)

    # Encrypt message
    ct1 = Evaluator.Encryption(m1)
    ct2 = Evaluator.Encryption(m2)

    op_text1+="Полиномы m1 и m2 зашифрованы в шифртексты ct1 и ct2:\n"
    op_text1+="ct1[0]: {}\n".format(ct1[0])
    op_text1+="ct1[1]: {}\n".format(ct1[1])
    op_text1+="ct2[0]: {}\n".format(ct2[0])
    op_text1+="ct2[1]: {}\n".format(ct2[1])
    add_text+=op_text1
    sub_text+=op_text1
    mul_text+=op_text1
    
    # Homomorphic Addition
    ct = Evaluator.HomomorphicAddition(ct1,ct2)
    mt = Evaluator.Decryption(ct)

    nr = Evaluator.IntDecode(mt) 
    ne = (n1+n2) 

    add_text+="Сложение шифртекстов ct_add = Enc(m1) + Enc(m2)\n"
    add_text+="ct_add[0] :{}\n".format(ct[0])
    add_text+="ct_add[1] :{}\n".format(ct[1])
    add_text+="Расшифровка результата в полином ct_dec = Dec(ct_add)\n"
    add_text+="ct_dec    :{}\n".format(mt)
    add_text+="Декодирование полинома в данные ct_dcd = Decode(ct_dec)\n"
    add_text+="ct_dcd    :{}\n".format(nr)

    if nr == ne:
        add_text+="Гомоморфное сложение успешно.\n"
    else:
        add_text+="Гомоморфное сложение неудачно.\n"

    # Homomorphic Subtraction
    ct = Evaluator.HomomorphicSubtraction(ct1,ct2)
    mt = Evaluator.Decryption(ct)

    nr = Evaluator.IntDecode(mt) 
    ne = (n1-n2) 

    sub_text+="Вычитание шифртекстов ct_sub = Enc(m1) - Enc(m2)\n"
    sub_text+="ct_sub[0] :{}\n".format(ct[0])
    sub_text+="ct_sub[1] :{}\n".format(ct[1])
    sub_text+="Расшифровка результата в полином ct_dec = Dec(ct_sub)\n"
    sub_text+="ct_dec    :{}\n".format(mt)
    sub_text+="Декодирование полинома в данные ct_dcd = Decode(ct_dec)\n"
    sub_text+="ct_dcd    :{}\n".format(nr)

    if nr == ne:
        sub_text+="Гомоморфное вычитание успешно.\n"
    else:
        sub_text+="Гомоморфное вычитание неудачно.\n"

    # Multiply two message (no relinearization)
    ct = Evaluator.HomomorphicMultiplication(ct1,ct2)
    mt = Evaluator.DecryptionV2(ct)

    nr = Evaluator.IntDecode(mt) 
    ne = (n1*n2)
    
    mul_text1=mul_text
    mul_text1+="Умножение шифртекстов ct_mul = Enc(m1) * Enc(m2) (без релинеаризации)\n"
    mul_text1+="ct_mul[0] :{}\n".format(ct[0])
    mul_text1+="ct_mul[1] :{}\n".format(ct[1])
    mul_text1+="Расшифровка результата в полином ct_dec = Dec(ct_sub)\n"
    mul_text1+="ct_dec    :{}\n".format(mt)
    mul_text1+="Декодирование полинома в данные ct_dcd = Decode(ct_dec)\n"
    mul_text1+="ct_dcd    :{}\n".format(nr)

    if nr == ne:
        mul_text1+="Гомоморфное умножение успешно.\n"
    else:
        mul_text1+="Гомоморфное умножение неудачно.\n"

    # Multiply two message (relinearization v1)
    ct = Evaluator.HomomorphicMultiplication(ct1,ct2)
    ct = Evaluator.RelinearizationV1(ct)
    mt = Evaluator.Decryption(ct)

    nr = Evaluator.IntDecode(mt)
    ne = (n1*n2)

    mul_text2=mul_text
    mul_text2+="Умножение шифртекстов ct_mul = Enc(m1) * Enc(m2) (с релинеаризацией)\n"
    mul_text2+="ct_mul[0] :{}\n".format(ct[0])
    mul_text2+="ct_mul[1] :{}\n".format(ct[1])
    mul_text2+="Расшифровка результата в полином ct_dec = Dec(ct_sub)\n"
    mul_text2+="ct_dec    :{}\n".format(mt)
    mul_text2+="Декодирование полинома в данные ct_dcd = Decode(ct_dec)\n"
    mul_text2+="ct_dcd    :{}\n".format(nr)
    

    if nr == ne:
        mul_text2+="Гомоморфное умножение успешно.\n"
    else:
        mul_text2+="Гомоморфное умножение неудачно.\n"
        
   
    if argx == "Сложение":
       res=add_text
    if argx == "Вычитание":
       res=sub_text
    if argx == "Умножение":
       res=mul_text1
    if argx == "Умножение (с релинеаризацией)":
       res=mul_text2
    open_popup(window, res)



window = Tk()
window.title('Калькулятор')
window.geometry('480x320')

frame = Frame(
   window,
   padx=100,
   pady=100
)
frame.pack(expand=True)

t_label = Label(
   frame,
   text="Введите t",
)
t_label.grid(row=2, column=1)

n_label = Label(
   frame,
   text="Введите n"
)
n_label.grid(row=3, column=1)

logq_label = Label(
   frame,
   text="Введите log(q)",
)
logq_label.grid(row=4, column=1)

n1_label = Label(
   frame,
   text="Введите 1-ое число",
)
n1_label.grid(row=5, column=1)

n2_label = Label(
   frame,
   text="Введите 2-ое число",
)
n2_label.grid(row=6, column=1)

t_tf = Entry(
   frame,
)
t_tf.grid(row=2, column=2, pady=5)
t_tf.insert(0, '256')
 
n_tf = Entry(
   frame,
)
n_tf.grid(row=3, column=2, pady=5)
n_tf.insert(0, '2048')

logq_tf = Entry(
   frame,
)
logq_tf.grid(row=4, column=2, pady=5)
logq_tf.insert(0, '37')

n1_tf = Entry(
   frame,
)
n1_tf.grid(row=5, column=2, pady=5)
n1_tf.insert(0, '687')

n2_tf = Entry(
   frame,
)
n2_tf.grid(row=6, column=2, pady=5)
n2_tf.insert(0, '705')

cb_label = Label(
   frame,
   text="Выберите тип операции",
)
cb_label.grid(row=7, column=1, pady=5)

arx1 = ttk.Combobox(frame,values=["Сложение", "Вычитание", "Умножение", "Умножение (с релинеаризацией)"], state="readonly")

arx1.grid(row=7, column=2, pady=5)
arx1.current(0)



cal_btn = Button(
   frame,
   text='Вычислить',
   command=calculate_fhe
)
cal_btn.grid(row=14, column=2, pady=5)

window.mainloop()