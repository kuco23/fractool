import math, cmath
from functools import reduce
from matplotlib import cm
from matplotlib.colors import Normalize
import numpy as np

I = 0 + 1j

overflow = 10**20
non_escape_count = 3
phisplit = 80
L = 1.000001

quadabs = lambda z: z.real * z.real + z.imag * z.imag

def diferentiate(poly):
    n, an = len(poly) - 1, poly[0]
    return [(n - i) * an for (i, an) in enumerate(poly[:-1])]

def horner(p, z):
    return reduce(lambda x, y: z * x + y, p)

def algorithmJulia(
    poly, radius, alg, iterlim, px
):
    n = len(poly) - 1
    an = abs(poly[0])
    C = sum(map(abs, poly)) - an
    eps = max(1, 2 * C / 2, pow(2 * L / an, 1 / (n-1)))
    eps_quad = pow(eps, 2)

    if alg == 'dem': dpdz = diferentiate(poly)

    def escapetimeCount(z, lim=iterlim):
        zk, count = z, 0
        while count < lim and quadabs(zk) <= eps_quad:
            zk = horner(poly, zk)
            count += 1
        return count

    def demCount(z):
        zk, dk = z, 1
        for _ in range(iterlim):
            if max(
                abs(zk.real) + abs(zk.imag),
                abs(dk.real) + abs(dk.imag)
            ) > overflow: break
            dk = horner(dpdz, zk) * dk
            zk = horner(poly, zk)
        abszk = abs(zk)
        if abszk < eps: return 0
        else:
            absdk = abs(dk)
            if absdk == 0: return -1
            estimate = math.log2(abszk) * abszk / absdk
            return -math.log2(estimate)

    def simulatedRadius():
        maxradius = 0
        r = eps / px
        rquad = r * r
        phi, dphi = 0, 2 * math.pi / phisplit
        while phi < 2 * math.pi:
            dz = r * cmath.exp(I * phi)
            z = px * dz
            while quadabs(z) > rquad:
                count = escapetimeCount(z, non_escape_count)
                if count == non_escape_count: break
                z -= dz
            if quadabs(z) > maxradius:
                maxradius = quadabs(z)
            phi += dphi
        return math.sqrt(maxradius) + 10 * r

    return (
        locals()[alg + 'Count'],
        radius or simulatedRadius()
    )


def algorithmMandelbrot(
    radius, alg, iterlim
):
    def escapetimeCount(c):
        ck, count = complex(0, 0), 0
        while count < iterlim and quadabs(ck) <= 4:
            ck *= ck
            ck += c
            count += 1
        return count

    def demCount(c):
        ck, dk = c, 1
        for _ in range(iterlim):
            if max(
                abs(ck.real) + abs(ck.imag),
                abs(dk.real) + abs(dk.imag)
            ) > overflow: break
            dk = 2 * ck * dk + 1
            ck *= ck
            ck += c
        absck = abs(ck)
        if absck <= 2: return 0
        else:
            absdk = abs(dk)
            if absdk == 0: return -1
            estimate = math.log2(absck) * absck / absdk
            return -math.log2(estimate)

    return (
        locals()[alg + 'Count'],
        radius or 2.2
    )


def drawFractalPPM(
    ppm, poly, center, radius, alg, iterlim,
    px, cmap, order, power, shift, perc
):
    colormap = cm.get_cmap(cmap)
    if order == 'reversed':
        colormap = colormap.reversed()
    
    if poly == 'mandelbrot':
        countAlgo, radius = algorithmMandelbrot(
            radius, alg, iterlim
        )
    else:
        poly = list(map(complex, poly.split()))
        countAlgo, radius = algorithmJulia(
            poly, radius, alg, iterlim, px
        )

    def pointGenerator():
        cx, cy = center.real, center.imag
        ReS, ReT = cx - radius, cx + radius
        ImS, ImT = cy - radius, cy + radius
        dim = dre = 2 * radius / px
        dz = complex(dre, 0)
        zk = complex(ReS, ImS)
        for i in range(px):
            for j in range(px):
                zk += dz
                yield zk
            zk = complex(ReS, zk.imag + dim)
    
    points = pointGenerator()
    pixels = np.empty((px, px), dtype=float)
    for j in range(px):
        for i in range(px):
            pixels[i][j] = countAlgo(next(points))
            
    m, M = pixels.min(), pixels.max()
    if alg == 'dem': pixels[pixels == -1] = m
    pixels[pixels == 0] = m if shift else M
    normed = Normalize(m, M)(pixels)
        
    if perc is not None:
        p, pf, pc = perc
        q = np.percentile(normed, p)
        filt = normed <= q
        normed[filt] = pow(normed[filt], pf)
        filt = ~filt
        normed[filt] = pow(normed[filt], pc)
    elif power != 1:
        normed = pow(normed, power)

    gradient = colormap(normed)
    ppm.writelines(['P3\n', f'{px} {px}\n', '255\n'])
    for j in range(px):
        ppm.write('\n')
        for i in range(px):
            rgb = map(round, 255 * gradient[i][j][:3])
            ppm.write(' '.join(map(str, rgb)) + '  ')
    
          
if __name__ == '__main__':

    from argparse import ArgumentParser

    args = ArgumentParser()
    args.add_argument(
        'poly', type = str,
        metavar = 'polynomial coefficients'
    )
    args.add_argument(
        '-c', metavar=('x', 'y'),
        type=float, nargs=2,
        default=(0, 0),
        help='center point'
    )
    args.add_argument(
        '-r', metavar='radius',
        type=float, default=None,
        help='radius around the center point'
    )
    args.add_argument(
        '-px', metavar='pixels',
        type=int, default=1000,
        help='image pixels (px * px)'
    )
    args.add_argument(
        '-it', metavar='iterations',
        type=int, default=250,
        help='number of iterations to aproximate limit'
    )
    args.add_argument(
        '-alg', metavar='algorithm',
        type=str, default='dem',
        choices = ['escapetime', 'dem'],
        help='algorithm name'
    )
    args.add_argument(
        '-cm', metavar='colormap',
        type=str, default='viridis',
        help='color mapping'
    )
    args.add_argument(
        '-cmp', metavar='colormap power',
        type=float, default=1,
        help='colormap power normalization'
    )
    args.add_argument(
        '-cmo', metavar='colormap order',
        type=str, default='reversed',
        choices = ['normal', 'reversed'],
        help='color mapping order'
    )
    args.add_argument(
        '-csh', metavar='color shift',
        type=bool, default=False,
        help=(
            "color the set interior with "
            "the outer most color"
        )
    )
    args.add_argument(
        '-cpc', metavar=('q', 'p1', 'p2'),
        type=float, default=None, nargs=3,
        help=(
            "q - percentile; "
            "p1 - power applied to x <= q'; "
            "p2 - power applied to x > q'"
        )
    )
    args.add_argument(
        '-fn', metavar='file name',
        type=str, default='fractal',
        help='image file name'
    )
    vals = args.parse_args()
    
    if vals.c != (0, 0) and vals.r is None:
        raise Exception(
            "if center is non-trivial, then "
            "radius must be provided"
        )
    elif (vals.cpc and vals.cmp != 1):
        raise Exception(
            "color map percentile and "
            "colormap power cannot be provided both"
        )
    
    filename = vals.fn + '.ppm'
    open(filename, 'a').close()
    with open(filename, 'w', encoding='utf-8') as ppm:
        drawFractalPPM(
            ppm, vals.poly, complex(*vals.c),
            vals.r, vals.alg, vals.it,
            vals.px, vals.cm, vals.cmo,
            vals.cmp, vals.csh, vals.cpc
        )
