import math

def f(x):
    return 3*x[0]**2 + x[1]**2 - x[0]*x[1] + x[0]


def gradient(x):
    h = 0.1 ** 6
    result = []
    for i in range(len(x)):
        x_temp = [x[j] if j != i else x[j] + h for j in range(len(x))]
        result.append((f(x_temp) - f(x)) / h)
    return result


def vector_mul_const(c, x):
    return [c * x[i] for i in range(len(x))]


def vector_sum(x1, x2):
    return [x1[i] + x2[i] for i in range(len(x1))]


def norm(x):
    return max([math.fabs(elem) for elem in x])


def get_t(x):
    a, b, eps = -100, 100, 0.1 ** 6
    delta = eps / 2
    gradxk = gradient(x)
    while True:
        x1 = (a + b - delta) / 2
        X_new = [x[i] - x1 * gradxk[i] for i in range(len(x))]
        y1 = (a + b + delta) / 2
        Y_new = [x[i] - y1 * gradxk[i] for i in range(len(x))]
        if f(X_new) < f(Y_new):
            b = y1
        else:
            a = x1
        if math.fabs(b - a) <= 2 * eps:
            x1 = (a + b) / 2
            return x1


N = 10
xk = [1.5, 1.5]
eps1, eps2 = 0.1, 0.15
k = 0
xres = [-1, -1]
gradxk = 0
flag = 0
print("1. Зададим x0=[", xk[0], ",", xk[1], "], eps1= ", eps1, ", eps2= ", eps2, ", M= ", N, sep='')
print("2. Зададим k=", k)

while True:
    gradxk = gradient(xk)
    print("\n3.", k, " Вычислим градиент в x", k, "= (", xk[0], ",", xk[1], "):\n\t grad(f(x", k, "))= (", gradxk[0],
          ",", gradxk[1], ")", sep='')
    normgradxk = norm(gradxk)
    print("4.", k, " Вычислим ||grad(f(x", k, "))||= ", normgradxk, sep='')

    if normgradxk < eps1:
        xres = xk
        print("4.", k, " ||grad(f(x", k, "))||= ", normgradxk, " <", " eps1=", eps1, ", а значит x_res=x", k, "=(",
              xk[0], ",", xk[1], ") f(x_res)=", f(xres), sep='')
        break

    print("4.", k, " ||grad(f(x", k, "))||=", normgradxk, " >", " eps1= ", eps1, ", а значит переходим к шагу 5",
          sep='')
    if k >= N:
        print("5.", k, " Количество итераций больше либо равно M=", N, ", а значит x_res= x", k, " =(", xk[0], ",",
              xk[1],
              ")", sep='')
        xres = xk
        break

    print("5.", k, " Количество итераций меньше M= ", N, ", а значит а значит переходим к шагу 6", sep='')
    t = get_t(xk)
    print("6.", k, " Определим величину шага t*_", k, " из условия, что p(t*_", k, ")= f(x", k, "-t*_", k, " *grad(f(x",
          k, ")))->min:", "\n\t t", k, "= ", t, sep='')
    xk_next = vector_sum(xk, vector_mul_const(-1 * t, gradxk))
    print("7.", k, " Вычислим x", k + 1, "= x", k, "-t*_", k, " *grad(f(x", k, "))= (", xk_next[0], ",", xk_next[1],
          ")",
          sep='')
    norm1 = norm(vector_sum(xk_next, vector_mul_const(-1, xk)))
    norm2 = math.fabs(f(xk_next) - f(xk))
    print("8.", k, " ||x", k + 1, "-x", k, "||= ", norm1, ", |f(x", k + 1, ")-f(x", k, ")|= ", norm2, sep='')

    if norm1 < eps2 and norm2 < eps2:
        flag += 1
        print("Условия выполняются еще раз, т.к. ||x", k + 1, "-x", k, "||= ", norm1, " < ", "eps2= ", eps2,
              ", и |f(x", k + 1, ")-f(x", k, ")|= ", norm2, " < ", "eps2= ", eps2, sep='')
    else:
        print("Условия не выполняются")
        flag = 0

    if flag >= 2:
        xres = xk_next
        print("\nУсловия выполняются как минимум два раза, а значит:\n x_res= x", k + 1, "= (", xres[0], ",", xres[1],
              "), ||grad(f(x", k + 1, "))||= ", norm(gradient(xk_next)), " f(x_res)= ", f(xres), sep='')
        break
    print("Переходим к следующей итерации")

    xk = xk_next
    k += 1