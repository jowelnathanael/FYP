import cProfile
import time
import sage.all as sage
import numpy as np
import itertools as it
import random
import matplotlib.pyplot as plt

def random_code(config):
    q = config["q"]
    n = config["n"]
    k = config["k"]
    F = config["F"]
    gen = config["gen"]
    Q = config["Q"]
    primepower = config["primepower"]
    eta = config["eta"]
    e = config["e"]
    s = config["s"]
    t = config["t"]
    field_dict = config["field_dict"]
    index_dict = config["index_dict"]
    score = config["score"]
    A = config["A"]
    
    m = []
    for i in range(k):
        m += [random.choice(Q)] 
    #m = [1,0,1,0,0]
    #print(m)
    #print("m:",m)
    codeword = encode(m,config)
    return codeword

def random_eval_pts(config):
    
    F = config["F"]
    n = config["n"]
    
    B = F.list().copy()
    B.remove(0)
    A = random.sample(B, n)   
    return A
        

def encode(m,config): # message of length k >> codeword of length n
    q = config["q"]
    n = config["n"]
    k = config["k"]
    F = config["F"]
    gen = config["gen"]
    Q = config["Q"]
    primepower = config["primepower"]
    eta = config["eta"]
    e = config["e"]
    s = config["s"]
    t = config["t"]
    field_dict = config["field_dict"]
    index_dict = config["index_dict"]
    score = config["score"]
    A = config["A"]

    codeword = []
    
    for i in range(n):
        c_i = 0
        for j in range(k):        
            c_i += m[j]*A[i]**j
        codeword += [c_i]
           
    return codeword

def encode2(m,config): # message of length k >> codeword of length n
    q = config["q"]
    n = config["n"]
    k = config["k"]
    F = config["F"]
    gen = config["gen"]
    Q = config["Q"]
    primepower = config["primepower"]
    eta = config["eta"]
    e = config["e"]
    s = config["s"]
    t = config["t"]
    field_dict = config["field_dict"]
    index_dict = config["index_dict"]
    score = config["score"]
    A = config["A"]

    codeword = []
    
    for i in range(n):
        c_i = 0
        for j in range(k):        
            c_i += m[j]*A[i]**j
        codeword += [c_i]
           
    return codeword

def encode_f(f,config): #message of length k in the form of f >> codeword of length n
    q = config["q"]
    n = config["n"]
    k = config["k"]
    F = config["F"]
    gen = config["gen"]
    Q = config["Q"]
    primepower = config["primepower"]
    eta = config["eta"]
    e = config["e"]
    s = config["s"]
    t = config["t"]
    field_dict = config["field_dict"]
    index_dict = config["index_dict"]
    score = config["score"]
    A = config["A"]

    codeword = []
    for i in range(len(A)):
        eval = f(A[i],0)
        codeword += [eval]

    return codeword

def encode_f1(f,config): #message of length k in the form of f >> codeword of length n
    q = config["q"]
    n = config["n"]
    k = config["k"]
    F = config["F"]
    gen = config["gen"]
    Q = config["Q"]
    primepower = config["primepower"]
    eta = config["eta"]
    e = config["e"]
    s = config["s"]
    t = config["t"]
    field_dict = config["field_dict"]
    index_dict = config["index_dict"]
    score = config["score"]
    A = config["A"]
    
    codeword = []
    for i in range(len(A)):
        eval = f(A[i])
        codeword += [eval]

    return codeword

def change_to_field(c,config):
    q = config["q"]
    n = config["n"]
    k = config["k"]
    F = config["F"]
    gen = config["gen"]
    Q = config["Q"]
    primepower = config["primepower"]
    eta = config["eta"]
    e = config["e"]
    s = config["s"]
    t = config["t"]
    field_dict = config["field_dict"]
    index_dict = config["index_dict"]
    score = config["score"]
    A = config["A"]
    
    y = []
    for i in c:
        y += [Q[i]]        
    return y

def change_to_index(c,config):
    q = config["q"]
    n = config["n"]
    k = config["k"]
    F = config["F"]
    gen = config["gen"]
    Q = config["Q"]
    primepower = config["primepower"]
    eta = config["eta"]
    e = config["e"]
    s = config["s"]
    t = config["t"]
    field_dict = config["field_dict"]
    index_dict = config["index_dict"]
    score = config["score"]
    A = config["A"]
       
    y = []
    
    for i in c:
        y += [index_dict[i]]       
    return y

def dist(x,y):
    
    return len([1 for i in range(len(x)) if x[i]!=y[i]]) 


def min_dist_decoder(A,k,y,l):

  for m in it.product(A, repeat = k):
    x = encode(m,A)
    if dist(y,x) <= l :
      print(x)

    
def error_code(c,e):

  y = sage.copy(c)
  positions = random.sample(range(len(c)), e)
  for i in positions:
    while y[i] == c[i]:
      y[i] = F.random_element();

  return y


def pascal_dict(T,config):
    q = config["q"]
    
    pascal  = {}

    for i in range(T+1):

        for j in range(i+1):

              if j == 0 or j == i: pascal[(i,j)] = 1
        
                #else: pascal[(i,j)] = pascal[(i-1,j-1)] +  pascal[(i-1,j)]
              else: pascal[(i,j)] = (pascal[(i-1,j-1)] +  pascal[(i-1,j)]) % q


    return pascal


def eval_dict(D,config):
    
    F = config["F"]

    eval = {}
    for i in F.list():
        for j in range(D+1):
            eval[(i,j)] = i**j
    return eval
        
    
def erasure_decoder(c):
    M = []
    y = []
    for i in c:
        a = i[1]
        b = i[0]
        
        row = []
        for j in range(k):
            row += [a**j]
            
        y += [[b]]
        M += [row]
   

    M   = sage.Matrix(F,M)
    y = sage.Matrix(F,y)    
    sol = M.solve_right(y)    
   
    m = []
    for i in sol:
        m += [i[0]]
        
    c = encode2(m)
    
    return c


# -------------------------
# Berlekamp = Welch for MBW
# -------------------------

# Majority-and-BW
def MBW(c,readBW,config):
    
        return int(c == algo_0(readBW,config))

def algo_0(c,config):
    
    q = config["q"]
    n = config["n"]
    k = config["k"]
    F = config["F"]
    gen = config["gen"]
    Q = config["Q"]
    primepower = config["primepower"]
    eta = config["eta"]
    e = config["e"]
    s = config["s"]
    t = config["t"]
    field_dict = config["field_dict"]
    index_dict = config["index_dict"]
    score = config["score"]
    A = config["A"]
      
    y = c.copy()
    M   = []
    RHS = []
    for i in range(n):
        ai = A[i]
        row = [ ai**j for j in range(e+k) ]
        row += [ - y[i]*ai**j for j in range(e) ]
        M += [row]
        RHS += [ [y[i]*ai**e] ]

    M   = sage.Matrix(F,M)
    RHS = sage.Matrix(F,RHS)
    
    try:
        sol = M.solve_right(RHS)
    except ValueError:
      #print("System is inconsistent!")
      return None


    R = F['X, Y']
    X,Y = R.gens()
    Q = 0

    for j in range(e+k):
        Q += sol[j,0]*X**j
    #print("QX:",Q)

    E = X**e
    for j in range(e):
        E += sol[e+k+j,0]*X**j
    #print("EX:", E)

    if Q%E == 0:
        #print("PX:", Q/E )
        return encode_f(Q/E,config)
    else:
        #print("No solution")
        return None
        
        
# --
# KV
# --
def mult_matrix(reliability,s,config):
    
    q = config["q"]
    n = config["n"]
    k = config["k"]
    F = config["F"]
    gen = config["gen"]
    Q = config["Q"]
    primepower = config["primepower"]
    eta = config["eta"]
    e = config["e"]
    s = config["s"]
    t = config["t"]
    field_dict = config["field_dict"]
    index_dict = config["index_dict"]
    score = config["score"]
    A = config["A"]
    
    a = len(reliability)
    b = len(reliability[0])
    M = np.zeros((a,b),dtype=int) 
    r = np.copy(reliability)
    
    
    
    for i in range(score):
        max_value = -1
        max_index= (0,0)
        
        for j in range(len(r)):
            for l in range(len(r[0])):
                if r[j][l] > max_value:
                    max_value = r[j][l]
                    max_index = (j,l)
        
        #print(max_value, max_index)
        
        
        j,l = max_index[0],max_index[1]
        r[j][l] = reliability[j][l]/(M[j][l]+2)
        old_mult = M[j][l]
        M[j][l] += 1
    
    
    mult = []
    cost = 0
    for i in range(b):
        for j in range(a):
            if M[j][i] > 0:
                mult.append([i,j,M[j][i]])
                cost += (M[j][i])*(M[j][i]+1)*0.5   

    return M,cost,mult


def cost_mult_matrix(M,config):
    q = config["q"]
    n = config["n"]
    k = config["k"]
    F = config["F"]
    gen = config["gen"]
    Q = config["Q"]
    primepower = config["primepower"]
    eta = config["eta"]
    e = config["e"]
    s = config["s"]
    t = config["t"]
    field_dict = config["field_dict"]
    index_dict = config["index_dict"]
    score = config["score"]
    A = config["A"]
    
    cost = 0
    for i in range(len(M)):
        for j in range(len(M[0])):
            cost += (M[i][j])*(M[i][j]+1)
            
    cost = cost/2
    
    return cost

def algo_3_mult(mult,cost,config):
    
    q = config["q"]
    n = config["n"]
    k = config["k"]
    F = config["F"]
    gen = config["gen"]
    Q = config["Q"]
    primepower = config["primepower"]
    eta = config["eta"]
    e = config["e"]
    s = config["s"]
    t = config["t"]
    field_dict = config["field_dict"]
    index_dict = config["index_dict"]
    score = config["score"]
    A = config["A"]
    
    Dy = sage.floor((1+(1+(8*cost/(k-1)))**0.5)/2)-1
    Dx = sage.floor(cost/(Dy+1)+0.5*(k-1)*Dy)

    powers = []
    for i in range(0,Dx+1):
        for j in range(0,Dy+1):
            if i + j*(k-1) <= Dx:
                powers += [[i,j]]

    #print(Dx,Dy,len(powers))

    binomial = pascal_dict(Dx,config)
    eval = eval_dict(Dx,config)
    #print(binomial)
    #print(eval)

    M = []

    for point in mult:
        #print("POINT")
        ai = A[point[0]]
        yi = Q[point[1]]
        m = point[2]
        #print(ai,yi,m)
              
        C = []
        for i in range(0,m):
            for j in range(0,m-i):
                C += [[i,j]]

        #print(C)
        for c in C:
            #print(c)
            row = []
            for power in powers:
                if power[0] >= c[0] and power[1] >= c[1]:
                    row += [binomial[(power[0],c[0])] * binomial[(power[1],c[1])] * eval[(ai,power[0])] * eval[(yi,power[1])]]
                else:
                    row += [0]

            M += [row]
            
    print(M)
    M = sage.Matrix(F,M)
    


    #print("Solving now...")
    sol = M.right_kernel()[1]
    #print(sol)
    #print(powers)

    R = F['X, Y']
    X,Y = R.gens()
    
    QQ = 0
    i = 0
    for power in powers:
        QQ += sol[i]*(X**power[0])*(Y**[power[1]])
        i+=1

    print(QQ)
    factors = QQ.factor()
    print(factors)
    #print(factors)
    L = []
    for f in factors:
        if (Y-f[0]).is_univariate():
            new_code = encode_f((Y-f[0]),config)

        #if dist(new_code,y) <= max(0,e):
            L += [Y-f[0]]
        elif (Y+f[0]).is_univariate():
            new_code = encode_f((Y+f[0]),config)
  #if dist(new_code,y) <= max(0,e):
            L += [Y+f[0]]

    #print("L = ", L)
    for i in L:
        a = encode_f(i,config)
    #print(L)   
    return L

def KV(c,reads,reliability,config):
    print("im here")
    q = config["q"]
    n = config["n"]
    k = config["k"]
    F = config["F"]
    gen = config["gen"]
    Q = config["Q"]
    primepower = config["primepower"]
    eta = config["eta"]
    e = config["e"]
    s = config["s"]
    t = config["t"]
    field_dict = config["field_dict"]
    index_dict = config["index_dict"]
    score = config["score"]
    A = config["A"]
    
    M,cost,mult = mult_matrix(reliability,score,config)
    #print(M,cost,mult)
    L = algo_3_mult(mult,cost,config)
    #print(L)
    
    #c = change_to_field(c)
    
    if len(L) == 0:
        #print("fail",c,L,reads)
        return 0   
    else:
        
        min_dist = float("inf")
        
        for c_l in L:
            c_l = encode_f(c_l,config)
            distance = 0
            for read in reads:
                #read = change_to_field(read)
                distance += dist(c_l,read)
            
            if distance < min_dist:
                min_dist = distance
                min_dist_c = c_l 
                
        if min_dist_c == c:
            return 1
        else:
            #print("fail wrong c",c,L,reads)
            
            return 0 

def assign_fg(G,a,y):
    #print("G,a,y:",G,a,y)
    g0,g1 = G[0] , G[1]
    #print("g0(a,y)=",g0(a,y))
    #print("g1(a,y)=",g1(a,y))
    if g0(a,y) == 0:
        f,g = g1,g0
    
    elif g1(a,y) == 0:
        f,g = g0,g1
    
    elif g0.weighted_degree(1,-1) < g1.weighted_degree(1,-1):
        f,g = g0,g1
    else:
        f,g = g1,g0
   
    if f(a,y) == 0:
        g = f
    
    else:
        g = g - (g(a,y)/f(a,y))*f
    
    
    return f,g

def print_section(section, time):
    print(f"{section}: {time:.4f}s")
        
def chase(c,reliability,config ,print_times=False):
    q = config["q"]
    n = config["n"]
    k = config["k"]
    F = config["F"]
    gen = config["gen"]
    Q = config["Q"]
    primepower = config["primepower"]
    eta = config["eta"]
    e = config["e"]
    s = config["s"]
    t = config["t"]
    field_dict = config["field_dict"]
    index_dict = config["index_dict"]
    score = config["score"]
    A = config["A"]
    
    start = time.time()
    r = reliability.copy().transpose()
    data = np.zeros((n,5),dtype=float)
    i = 0
    for row in r:
        max_row = max(row)
        j = np.where(row == max_row)[0]
        if len(j) > 1:
            data[i][0] = max_row
            data[i][1] = j[0]
            data[i][2] = max_row
            data[i][3] = j[1]
            data[i][4] = 1
                   
        else:
            new_row = np.delete(row,j[0])
            max_new_row = max(new_row)
            h = np.where(row == max_new_row)[0][0]
            #print(max(new_row),h)
            l = max_new_row/max_row

            data[i][0] = max_row
            data[i][1] = j[0]
            data[i][2] = max_new_row
            data[i][3] = h
            data[i][4] = l

        i += 1
                
    #print(data.transpose())    
    gammas = data.transpose()[4]
    gamma_order = []
    
    
    #print("gammas",gammas)
    #print("gamma_order",gamma_order)
    row  = gammas
    while len(gamma_order) != n:
        max_row = max(row)
        j = np.where(gammas == max_row)[0]
        for i in j:
            gamma_order += [i]
            h = np.where(row == max_row)[0][0] 
            row = np.delete(row,h)
    
    #print(gamma_order)
    end = time.time()
    
    if print_times: print_section("Choosing reliabilities", end-start)
    
    start = time.time()
    I = gamma_order[:eta]
    Ibar = gamma_order[eta:n]
    J = Ibar[:k]
    Jbar = Ibar[k:]
    
    #print("J",J)
    #print("J_:",Jbar)
    #print("I:",I)
    
    A_ = []
    A_index = []
    yHD = []
    y2HD = []
    
    #print(k)
    for i in range(k):
        eval_pt = J[i]
        y_pt = int(data[eval_pt][1])
        
        A_index += [eval_pt]
        A_ += [A[eval_pt]]
        yHD += [Q[y_pt]]
        y2HD += [Q[y_pt]]
    
    for i in range(n-eta-k):
        eval_pt = Jbar[i]
        y_pt = int(data[eval_pt][1])
        
        A_index += [eval_pt]
        A_ += [A[eval_pt]]
        yHD += [Q[y_pt]]
        y2HD += [Q[y_pt]]
        
    for i in range(eta):
        eval_pt = I[i]
        y_pt = int(data[eval_pt][1])
        y2_pt = int(data[eval_pt][3])
        
        A_index += [eval_pt]
        A_ += [A[eval_pt]]
        yHD += [Q[y_pt]]
        y2HD += [Q[y2_pt]]
    
    #reencoding
    y_ = []
    for i in range(k):
        y_ += [(A_[i],yHD[i])]
    #print(y_)
    LP = F['x']
    ep = LP.lagrange_polynomial(y_)
    #print(ep)
    psi_ = encode_f1(ep,config)          

#     print("A_index:",A_index)
#     print("A_:", A_)
#     print("yHD:",yHD)
#     print("y2HD:",y2HD)
#     print("psi_",psi_)    
    
    psi = []
    for i in range(n):
        j = A_index[i]
        psi += [psi_[j]]
    
#   print("psi:",psi)
    end = time.time()
    if print_times: print_section("Choosing sets", end-start)
    
    start = time.time()
    # start of algorithm
    S = sage.PolynomialRing(F,'x,z')
    x,z = S.gens()

    v = 1
    for i in range(k):
        v*= (x-A_[i])
    
    #print("v(x):", v.factor())
    g0 = x**0
    g1 = z
    G = [g0,g1]
    for i in range((n-k)-eta):
        j = i + k
        a = A_[j]
        y = (yHD[j] - psi[j])/v(a,0)

        fHD,gHD = assign_fg(G,a,y)
        fHD = (x-a)*fHD
        
        G = [fHD,gHD]
       
    W = dict()
    end = time.time()
    if print_times: print_section("Algo 1", end-start)
    W[(0,0)] = G
    
    start = time.time()

    for j in range(0,eta):
        
        i = n - eta + j
        a = A_[i]
        y = (yHD[i] - psi[i])/v(a,0)
        y2 = (y2HD[i] - psi[i])/v(a,0)
        
        for m in range(2**(j)):
            
            G = W[(m,j)]
            
            fHD,gHD = assign_fg(G,a,y)
            fHD = (x-a)*fHD
            f2HD,g2HD = assign_fg(G,a,y2)
            f2HD = (x-a)*f2HD

            W[(2*m,j+1)] = [fHD,gHD]
            W[(2*m + 1,j+1)] = [f2HD,g2HD]
            
            
    #print(W)
    end = time.time()
    if print_times: print_section("Algo 2", end-start)
    
    start = time.time()
    Q_ = []
    for i in range(2**eta):
        
        w0,w1 = W[(i,eta)]
        #print(w0,w1)
        
        if w0.weighted_degree({z:-1, x:1}) < w1.weighted_degree({z:-1, x:1}):
            Q_.append(w0)
        else: 
            Q_.append(w1)
    end = time.time()
    if print_times: print_section("Choosing min degree", end-start)
        
    start = time.time()
    Q_m = []
    Q_m2= []
    Q_m3= []
    
    Q_ = list(dict.fromkeys(Q_))

    for q in Q_:
        if q !=0:
            
            start = time.time()
            d = q.polynomial(z)
            d_0 = d.list()[0]
            d_1 = d.list()[1]
            end = time.time()
            if print_times: print_section("Factoring by z", end-start)
            
            start = time.time()
            Q_m3.append((d_0*v).quo_rem(d_1))
            #Q_m3.append(((d_0*v)-(d_0*v).reduce([d_1]))/d_1)
            end = time.time()
            if print_times: print_section("Remainder", end-start)
            
    end = time.time()
    if print_times: print_section("Change variable 2", end-start)

    M2 = []
    for q_m3 in Q_m3:
        if (q_m3[1] == 0) and q_m3[1].degree() < k:
            M2 += [q_m3[0]]
    #print(M)
    #print(M2)

    #print(M)
    L = []
    for m in M2:
        c_ = encode_f(m,config)
        L+= [c_]
    
    #print(L)
    LL = []
    for c_ in L:
        co = np.array(c_) + np.array(psi_)
        LL += [list(co)]
    
    for i in LL:
        if c == i:
            #print("yes")
            return 1
    
    return 0


def produceReads(c,N,p,config):
    q = config["q"]
    n = config["n"]
    k = config["k"]
    F = config["F"]
    gen = config["gen"]
    Q = config["Q"]
    primepower = config["primepower"]
    eta = config["eta"]
    e = config["e"]
    s = config["s"]
    t = config["t"]
    field_dict = config["field_dict"]
    index_dict = config["index_dict"]
    score = config["score"]
    A = config["A"]
    error_model = config["error_model"]
    burst_length = config["burst_length"]
   
    start_time = time.time()
    
    reads = []
    count = np.zeros((q, n), dtype=int)
    
 
    for _ in range(N):
        # Introduce errors according to the specified error model
        errored_codeword = introduce_errors(c, error_model, p, burst_length,config)
        reads.append(errored_codeword)

        for j, symbol in enumerate(errored_codeword):
            index = index_dict[symbol]
            count[index][j] += 1

    reliability = count / N
    readBW = [Q[np.argmax(col)] for col in count.T]

    end_time = time.time()
    #print("Produce Reads time:", end_time - start_time)
    
    return reads, readBW, reliability

#     reads = []
#     count = np.zeros((q,n),dtype=int)
    
#     for i in range(N):
#         y = sage.copy(c)
        
#         for j in range(n):
#             coeff = y[j].list() 
            
#             for l in range(t): 
#                 a = coeff[l]
#                 if random.random() > p:
#                      while coeff[l] == a:
#                         coeff[l] = random.randint(0,s-1)
            
#             y[j] = field_dict[tuple(coeff)]
#             #print(coeff)
#             index = index_dict[y[j]]
#             #print(index)
#             count[index][j] += 1
            
#         reads.append(y)
#     reliability = np.zeros((q,n))
    
#     for i in range(len(count)):
#         for j in range(len(count[0])):
#             reliability[i][j] = count[i][j]/N

#     readBW = []
#     count = np.transpose(count)

#     for i in count:        
#         j = np.where(i == max(i))[0][0]
#         #print(j)
#         readBW.append(Q[j]) 
#     end_time = time.time()
  
    #readBW = [Q[np.argmax(col)] for col in count.T]

    #print("Produce Reads:",end_time-start_time)
    
#     reads = []
#     count = np.zeros((q, n), dtype=int)
    
#     for _ in range(N):
#         # Introduce errors according to the specified error model
#         errored_codeword = introduce_errors(c, error_model, error_rate, burst_length)
#         reads.append(errored_codeword)

#         for j, symbol in enumerate(errored_codeword):
#             index = index_dict[symbol]
#             count[index][j] += 1
    
#     reliability = count / N
#     readBW = [Q[np.argmax(col)] for col in count.T]
    
#     end_time = time.time()
#     print("Produce Reads time:", end_time - start_time)
    
    return reads, readBW, reliability

def introduce_burst_and_random_errors(codeword, burst_length, error_rate,config):
    Q = config["Q"]
    errored_codeword = []
    n = len(codeword)
    
    # Determine the starting point for the burst error
    start_point = random.randint(0, n - burst_length)
    end_point = start_point + burst_length

    for i in range(n):
        if start_point <= i < end_point:
            # Apply burst errors within the burst range
            if random.random() < error_rate:
                new_symbol = random.choice([x for x in Q if x != codeword[i]])
                errored_codeword.append(new_symbol)
            else:
                errored_codeword.append(codeword[i])
        else:
            # Apply random errors outside the burst range
            if random.random() < error_rate:
                new_symbol = random.choice([x for x in Q if x != codeword[i]])
                errored_codeword.append(new_symbol)
            else:
                errored_codeword.append(codeword[i])

    return errored_codeword

def introduce_random_errors(codeword, error_rate,config):
    errored_codeword = []
    for symbol in codeword:
        if random.random() < error_rate:
            new_symbol = random.choice(Q)  # Assuming Q is your finite field
            errored_codeword.append(new_symbol)
        else:
            errored_codeword.append(symbol)
    return errored_codeword

def introduce_random_errors_coeff(codeword, error_rate, config):
    
    q = config["q"]
    n = config["n"]
    k = config["k"]
    F = config["F"]
    gen = config["gen"]
    Q = config["Q"]
    primepower = config["primepower"]
    eta = config["eta"]
    e = config["e"]
    s = config["s"]
    t = config["t"]
    field_dict = config["field_dict"]
    index_dict = config["index_dict"]
    score = config["score"]
    A = config["A"]
    error_model = config["error_model"]
    burst_length = config["burst_length"]
    
 
    errored_codeword = []
    for symbol in codeword:
        coeff = symbol.list()

        for idx, c in enumerate(coeff):
            if random.random() > error_rate:
                # Create a list of possible new values excluding the current one
                possible_values = [x for x in range(s) if x != c]
                # Randomly select a new value from the possible ones
                new_value = random.choice(possible_values)

                coeff[idx] = new_value

        # Convert the altered coefficients back to a symbol in the finite field
        new_symbol = field_dict[tuple(coeff)]
        errored_codeword.append(new_symbol)

    return errored_codeword

        
def get_position_dependent_error_rate(base_error_rate, position, n):
    # Example function: linearly increasing error rate based on position
#     return base_error_rate - (position / n) * base_error_rate *0.8
    return base_error_rate - (position / (n)) * base_error_rate*base_error_rate
def introduce_position_dependent_errors(codeword, base_error_rate,config):
    
    Q = config["Q"]
    
    errored_codeword = []
    n = len(codeword)
    for i, symbol in enumerate(codeword):
        # Calculate the error rate for this position
        error_rate = get_position_dependent_error_rate(base_error_rate, i, n)
       
        if random.random() > error_rate:
            # Introduce an error by changing the symbol
            new_symbol = random.choice([x for x in Q if x != symbol])
            errored_codeword.append(new_symbol)
        else:
            errored_codeword.append(symbol)
    return errored_codeword

def introduce_errors(codeword, error_model, error_rate, burst_length,config):
    
    q = config["q"]
    n = config["n"]
    k = config["k"]
    F = config["F"]
    gen = config["gen"]
    Q = config["Q"]
    primepower = config["primepower"]
    eta = config["eta"]
    e = config["e"]
    s = config["s"]
    t = config["t"]
    field_dict = config["field_dict"]
    index_dict = config["index_dict"]
    score = config["score"]
    A = config["A"]
    error_model = config["error_model"]
    burst_length = config["burst_length"]
    
    
    
    
    if error_model == 'random':
        return introduce_random_errors(codeword, error_rate,config)
    elif error_model == 'random_coeff':
        return introduce_random_errors_coeff(codeword, error_rate,config)
    elif error_model == 'position':
         return introduce_position_dependent_errors(codeword, error_rate,config)
    elif error_model == 'burst':
        if burst_length is None:
            raise ValueError("Burst length must be provided for burst errors")
        return introduce_burst_errors(codeword, burst_length, error_rate,config)
    else:
        raise ValueError("Unknown error model")



def round(c,N,p,config):
    
    reads, readBW, reliability = produceReads(c,N,p,config)
    MBW_result = MBW(c,readBW,config)
    #print("MBW Done")
    KV_result = KV(c,reads,reliability,config)
    #KV_result = 0
    chase_result = chase(c,reliability,config)
    
    
    
    return [MBW_result,KV_result,chase_result]

def simulate_decoding(args):
    
    N, p_values, config = args
    
    mbw_success_rates = []
    KV_success_rates = []
    chase_success_rates = []

    for p in p_values:
        mbw_success = 0
        KV_success = 0
        chase_success = 0
        
        iterations = 500  # Number of iterations for each probability value

        for i in range(iterations):
            print(i)
            if i % 20 == 0:
                print("20 iterations done", N,p)

            c = random_code(config)  # Generate a random codeword
            round_result = round(c, N, p,config)
            mbw_result, KV_result, chase_result = round_result[0],round_result[1], round_result[2]
            mbw_success += mbw_result
            KV_success += KV_result
            chase_success += chase_result

        # Calculate success rates
        
        mbw_success_rates.append(mbw_success / iterations)
        KV_success_rates.append(KV_success / iterations)
        chase_success_rates.append(chase_success / iterations)

    return mbw_success_rates, KV_success_rates,chase_success_rates

