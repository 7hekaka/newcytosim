# Some utility methods to analyze collision of microtubules
# Francois Nedelec October 1--15 2021, Sep. 2022
# Copyright Sainsbury Laboratory, Cambridge University, UK

try:
    import os, sys, re, math
except ImportError as e:
    sys.stderr.write("Error loading module: %s\n"%str(e))
    sys.exit()

def get_collision(line, path=''):
    """scan one line of data and returns (angle, outcome)"""
    if len(line) < 2 or line[0] in ('%', '#'):
        return ()
    col = re.split('[\\s]+', line.strip())
    if len(col) == 4:  # Maud's analysis
        try:
            a = float(col[2])
        except:
            return ()
        table = 'UKXZ'
        r = table[int(col[3])]
        return (a, r)
    elif len(col) in {7, 8}:
        if col[0] == 'nan':
            return ()
        try:
            a = float(col[0])
        except:
            return ()
        r = col[-1]
        return (a, r)
    elif len(line) > 1:
        sys.stderr.write(f'Unexpected line in {path}: {line}\n')
    return ()
    
    
def get_collision_csv(line):
    """read a line and returns (angle, outcome)"""
    if len(line) < 2 or line[0] in ('%', '#'):
        return ()
    col = line.split(',')
    if len(col) == 6:  # Maud's analysis
        val = [ s.strip('"') for s in col ]
        if val[3] != 'Angle':
            a = float(val[3]) * math.pi/180
            if a > math.pi/2:
                print(val)
            t = val[4]
            r = '?'
            if t == 'catastrophe':
                r = 'K'
            elif t == 'zippering':
                r = 'Z'
            elif t == 'crossover':
                r = 'X'
            elif t == 'undefined':
                r = 'U'
            else:
                print(val)
            return (a, r)
    elif len(col) == 9:  # Carlos's picker
        a = float(col[6])
        r = col[7]
        return (a, r)
    elif len(col) > 1:
        if col[0] == 'nan':
            return ()
        a = float(col[0])
        r = col[-1]
        return (a, r)
    elif len(line) > 1:
        sys.stderr.write(f'Unexpected ({len(col)}): {line}')
    return ()


def collect_collision(data):
    """read lines in string, and return list of (angle, outcome)"""
    res = []
    for s in data:
        #print('Collision: '+s, end='')
        c = get_collision(s)
        if c:
            res.append(c)
    return res


def read_collision(path):
    """read file, and return list of (angle, outcome)"""
    res = []
    file = open(path, "r")
    for s in file:
        c = get_collision(s.rstrip('\n'), path)
        if c:
            res.append(c)
    file.close()
    return res


def read_collision_csv(path):
    """read file, and return list of (angle, outcome)"""
    res = []
    file = open(path, "r")
    for s in file:
        c = get_collision_csv(s.rstrip('\n'))
        if c:
            res.append(c)
    file.close()
    return res


def fold_angles(data):
    """ mirror-image all angles above pi/2 """
    cut = math.pi * 0.5
    res = [ x for x in data if x[0] < cut ]
    ser = [ (round(math.pi-a,5),c) for a,c in data if cut <= a ]
    #print(f" folded {len(ser)} angles above 90")
    return res + ser


def mirror_angles(data):
    """ mirror-image all angles around pi/2 """
    mirror = [ (round(math.pi-a,5),c) for a in data ]
    return data + mirror


def bin_collisions(data, nbins, sup):
    """
    data = list of ( angle, outcome )
    returns bin centers and array of counts for oucome
    """
    cnt = 0
    A = [0] * nbins
    K = [0] * nbins
    X = [0] * nbins
    Z = [0] * nbins
    B = [0] * nbins
    delta = sup / nbins
    for i in range(nbins):
        inf = delta * i
        sup = delta + inf
        #print("bin %i : [ %6.3f %6.3f ]" %(i, inf, sup))
        sel = [ c for a,c in data if ( inf <= a and a < sup ) ]
        #print(sel)
        A[i] = ( inf + sup ) * 0.5
        K[i] = sel.count('K')
        X[i] = sel.count('X')
        Z[i] = sel.count('Z') + sel.count('Y')
        B[i] = sel.count('B')
        cnt += K[i] + X[i] + Z[i]
    # put anything above into the last bin:
    if 1:
        sup = delta * nbins
        sel = [ c for a,c in data if ( sup <= a ) ]
        i = nbins - 1
        K[i] += sel.count('K')
        X[i] += sel.count('X')
        Z[i] += sel.count('Z') + sel.count('Y')
        B[i] += sel.count('B')
        if 0:
            print(f" events with angle above 90 : [ {sel.count('K')} {sel.count('X')} {sel.count('Z')} ]")
            for e in data:
                if e[0] > sup:
                    print(f"     {e[0]} {e[1]}\n")
    #print("distributed %i data points out of %i in %i bins" %(cnt, len(data), nbins))
    return (A, K, X, Z, B)


def dump_histogram(A, K, X, Z):
    """print bin centers and content"""
    nb = len(A)
    dA = ( A[1] - A[0] ) * 0.5
    print(f' {nb:2} bins:')
    for i in range(nb):
        t = repr(K[i] + X[i] + Z[i])
        print(f' {f:8} in [{A[i]-dA:6.2f} {A[i]+dA:6.2f}]', end='')
        integers = [x==int(x) for x in (K[i], X[i], Z[i])]
        if all(integers):
            print(f'{K[i]:5} {X[i]:5} {Z[i]:5}')
        else:
            print(f'{K[i]:5.2f} {X[i]:5.2f} {Z[i]:5.2f}')
    if 1:
        sK = sum(K)
        sX = sum(X)
        sZ = sum(Z)
        print(f'%% total: {sK:5} {sX:5} {sZ:5}   {sK+sX+sZ:6}')


def print_histogram(A, K, X, Z):
    """print bin centers and content"""
    nb = len(A)
    print('% angle       K     X     Z    total')
    for i in range(nb):
        t = K[i] + X[i] + Z[i]
        print(f' {A[i]:6.4f}   ', end='')
        integers = [x==int(x) for x in (K[i], X[i], Z[i])]
        if all(integers):
            print(f'{K[i]:5} {X[i]:5} {Z[i]:5}   {t:5}')
        else:
            print(f'{K[i]:5.2f} {X[i]:5.2f} {Z[i]:5.2f}   {t:5.2f}')
    if 1:
        sK = sum(K)
        sX = sum(X)
        sZ = sum(Z)
        print(f'%% total: {sK:5} {sX:5} {sZ:5}   {sK+sX+sZ:6}')


def normalize_data(A, K, X, Z):
    """for each bin, return total and number/total"""
    nb = len(A)
    T = [0] * nb
    nK = [0] * nb
    nX = [0] * nb
    nZ = [0] * nb
    for i in range(nb):
        T[i] = K[i] + X[i] + Z[i]
        if T[i] > 0:
            nK[i] = K[i] / T[i]
            nX[i] = X[i] / T[i]
            nZ[i] = Z[i] / T[i]
    return (T, nK, nX, nZ)


def load_histogram(file):
    """read histogram data from file"""
    A = []
    K = []
    X = []
    Z = []
    f = open(file, "r")
    for s in f:
        if len(s) < 1 or s[0] == '#' or s[0] == '%':
            continue
        if s.isspace():
            continue
        col = re.split('[\\s]+', s.strip())
        if len(col) < 4:
            sys.stderr.write(f"Error: histogram `{file}` must have at least 4 columns\n")
            sys.exit(1);
        T = 1.0
        if len(col) > 5:
            T = float(col[5])
        try:
            A.append(float(col[0]))
            K.append(round(T * float(col[1]), 3))
            X.append(round(T * float(col[2]), 3))
            Z.append(round(T * float(col[3]), 3))
        except Exception as e:
            sys.stderr.write(f"Error reading `{file}`: {e!s}\n")
            sys.exit(1);
    f.close()
    return (A, K, X, Z)


def histogram_deviation(K, X, Z, target):
    """Sum squared differences on all histogram columns"""
    tK = target[1]
    tX = target[2]
    tZ = target[3]
    #print_bins(target[0], tK, tX, tZ)
    nbin = len(target[0])
    res = 0
    for i in range(nbin):
        dK = K[i] - tK[i]
        dX = X[i] - tX[i]
        dZ = Z[i] - tZ[i]
        res += (dK)**2 + (dX)**2 + (dZ)**2
    return res
