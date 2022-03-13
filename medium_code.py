import cmath
from math import sqrt, log2
from functools import reduce

from numpy import zeros
from matplotlib import cm
from matplotlib.colors import Normalize


# evaluates the polynomial p at point z
def horner(p, z):
    return reduce(lambda x, y: x * z + y, p)

def escapetime(p, z, R, K):
    k, zk = 1, z
    while k < K and abs(zk) <= R:
        zk = horner(p, zk)
        k += 1
    return k

def drawPPM(filename, rgbfun, n):
    with open(filename, 'a') as _: pass
    with open(filename, 'w') as ppm:
        ppm.writelines(['P3\n', f'{n} {n}\n', '255\n'])
        for i in range(n):
            ppm.write('\n')
            for j in range(n):
                rgb = rgbfun(i, j)
                ppm.write(rgb + '  ')

def drawPPMCircle(n, c, r):
    def rgbcircle(i, j):
        radius = sqrt((i-c)**2 + (j-c)**2)
        return '0 0 0' if radius < r else '255 255 255'
    drawPPM('circle.ppm', rgbcircle, n)

def mapToComplexPlaneCenter(n, c, r, i, j):
    return c + r * complex(2 * j / n - 1, 2 * i / n - 1)

def drawEscapetimeMandelbrot(n, ctr, r, colormap, K):
    p = lambda c: [1, 0, c]
    def rgbMandelbrot(i, j):
        c = mapToComplexPlaneCenter(n, ctr, r, i, j)
        k = escapetime(p(c), 0, 2, K)
        if k == K: return '0 0 0'
        else:
            cmap = colormap(k / K)[:3]
            return ' '.join(str(round(255 * cm)) for cm in cmap)
    drawPPM('escapetime_mandelbrot.ppm', rgbMandelbrot, n)

def radiusJulia(poly, L):
    n = len(poly) - 1
    an = abs(poly[0])
    C = sum(map(abs, poly)) - an
    return max(1, 2 * C / 2, pow(2 * L / an, 1 / (n-1)))

def drawEscapetimeJulia(n, p, colormap, K):
    rp = radiusJulia(p, 1.000001)
    def rgbJulia(i, j):
        z = mapToComplexPlaneCenter(n, 0, rp, i, j)
        k = escapetime(p, z, rp, K)
        if k == K: return '0 0 0'
        else:
            cmap = colormap(k / K)[:3]
            return ' '.join(str(round(255 * cm)) for cm in cmap)
    drawPPM('escapetime_julia.ppm', rgbJulia, n)

def demMandelbrot(c, K, overflow):
    ck, dk = c, 1
    for _ in range(K):
        if max(
            abs(ck.real), abs(ck.imag),
            abs(dk.real), abs(dk.imag)
        ) > overflow: break
        dk = 2 * ck * dk + 1
        ck *= ck
        ck += c
    absck = abs(ck)
    if absck <= 2: return 0
    else:
        absdk = abs(dk)
        if absdk == 0: return -1 # this will probably never happen
        estimate = log2(absck) * absck / absdk
        return -log2(estimate)

def drawDemMandelbrot(n, ctr, r, colormap, K, overflow):
    arr = zeros((n, n), dtype=float)
    
    for i in range(n):
        for j in range(n):
            c = mapToComplexPlaneCenter(n, ctr, r, i, j)
            arr[i,j] = demMandelbrot(c, K, overflow)

    m, M = arr.min(), arr.max()
    arr[arr == 0] = M # 0 only denotes the inner set and it could spoil our normalization
    colortable = colormap(Normalize(m, M)(arr))
            
    def rgbMandelbrot(i, j):
        if arr[i, j] == M: return '0 0 0'
        cmap = colortable[i,j][:3]
        return ' '.join(str(round(255 * cm)) for cm in cmap)

    drawPPM('demMandelbrot.ppm', rgbMandelbrot, n)

# derivative of the given polynomial
def differentiate(poly):
    n, an = len(poly) - 1, poly[0]
    return [(n - i) * an for (i, an) in enumerate(poly[:-1])]

def demJulia(p, dp, z, K, R, overflow):
    zk, dk = z, 1
    for _ in range(K):
        if max(
            abs(zk.real), abs(zk.imag),
            abs(dk.real), abs(dk.imag)
        ) > overflow: break
        dk = horner(dp, zk) * dk
        zk = horner(p, zk)
    abszk = abs(zk)
    if abszk < R: return 0
    else:
        absdk = abs(dk)
        if absdk == 0: return -1
        estimate = log2(abszk) * abszk / absdk
        return -log2(estimate)

def drawDemJulia(n, p, colormap, K, pow_, overflow):
    arr = zeros((n, n), dtype=float)
    dp = differentiate(p)
    r = radiusJulia(p, 1.000001)

    for i in range(n):
        for j in range(n):
            z = mapToComplexPlaneCenter(n, 0, r, i, j)
            arr[i,j] = demJulia(p, dp, z, K, r, overflow)

    m, M = arr.min(), arr.max()
    arr[arr == 0] = M
    normalized = Normalize(m, M)(arr)
    adjusted = pow(normalized, pow_)
    colortable = colormap(adjusted)
            
    def rgbJulia(i, j):
        if arr[i, j] == M: return '0 0 0'
        cmap = colortable[i,j][:3]
        return ' '.join(str(round(255 * cm)) for cm in cmap)

    drawPPM('demJulia.ppm', rgbJulia, n)
    

def inCardioidOrCircle(c):
    s = cmath.sqrt(1 - 4 * c)
    return abs(1 + s) <= 1 or abs(1 - s) <= 1 or abs(1 + c) <= 1/4


if __name__ == '__main__':

    drawDemMandelbrot(
        2000, -0.8, 1.4,
        cm.get_cmap('gist_stern').reversed(),
        100, 10**20
    )
    
