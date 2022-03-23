from sys import argv
import math, cmath
from functools import reduce
from matplotlib import cm
from matplotlib.colors import Normalize
from numpy import empty, percentile, nan
from tqdm import tqdm

I = 0 + 1j
L = 1.000001
OVERFLOW = 10**20
EPS = 0.000001
PHISPLIT = 80
NON_ESCAPE_COUNT = 3

quadabs = lambda z: z.real * z.real + z.imag * z.imag

def differentiate(p):
    n = len(p) - 1
    return [(n - i) * an for (i, an) in enumerate(p[:-1])]

def horner(p, z):
    return reduce(lambda x, y: z * x + y, p)

def juliaRadius(poly, l=L):
    n = len(poly) - 1
    an = abs(poly[0])
    C = sum(map(abs, poly)) - an
    return max(1, 2 * C / 2, pow(2 * l / an, 1 / (n-1)))

def escapetimeMandelbrot(c, K):
    k, ck = 1, c
    while k < K and quadabs(ck) <= 4:
        ck *= ck
        ck += c
        k += 1
    return k if k < K else 0

def escapetimeJulia(z, p, K, R2):
    k, zk = 1, z
    while k < K and quadabs(zk) <= R2:
        zk = horner(p, zk)
        k += 1
    return k if k < K else 0

def demMandelbrot(c, K, overflow=OVERFLOW):
    ck, dk = c, 1
    for _ in range(K):
        if max(
            abs(ck.real), abs(ck.imag), 
            abs(dk.real), abs(dk.imag)
        ) > overflow: break
        dk *= 2 * ck
        dk += 1
        ck *= ck
        ck += c
    absck = abs(ck)
    if absck <= 2: return 0
    absdk = abs(dk)
    if absdk == 0: return nan # rarely happens
    estimate = math.log2(absck) * absck / absdk
    return -math.log2(estimate)

def demJulia(z, p, dp, K, R, overflow=OVERFLOW):
    zk, dk = z, 1
    for _ in range(K):
        if max(
            abs(zk.real) + abs(zk.imag),
            abs(dk.real) + abs(dk.imag)
        ) > overflow: break
        dk = horner(dp, zk) * dk
        zk = horner(p, zk)
    abszk = abs(zk)
    if abszk < R: return 0
    absdk = abs(dk)
    if absdk == 0: return nan # rarely happens
    estimate = math.log2(abszk) * abszk / absdk
    return -math.log2(estimate)

# simulates radius for julia set to find its position
def simulateRadius(
    p, R, px, phisplit=PHISPLIT, 
    non_escape_count=NON_ESCAPE_COUNT
):
    R2 = R * R
    maxradius = 0
    r = R / px
    rquad = r * r
    phi, dphi = 0, 2 * math.pi / phisplit
    while phi < 2 * math.pi:
        dz = r * cmath.exp(I * phi)
        z = px * dz
        while quadabs(z) > rquad:
            count = escapetimeJulia(
                z, p, non_escape_count, R2)
            if count == non_escape_count: break
            z -= dz
        if quadabs(z) > maxradius:
            maxradius = quadabs(z)
        phi += dphi
    return math.sqrt(maxradius) + 10 * r

def drawFractalPPM(
    ppm, algo, center, radius,
    px, colormap, power, cpc, ic,
    pb1, pb2
): 
    def pointGenerator():
        cx, cy = center.real, center.imag
        ReS = cx - radius
        ImS = cy - radius
        dim = dre = 2 * radius / px
        dz = complex(dre, 0)
        zk = complex(ReS, ImS)
        for _ in range(px):
            for _ in range(px):
                zk += dz
                yield zk
            zk = complex(ReS, zk.imag + dim)
    
    points = pointGenerator()
    pixels = empty((px, px), dtype=float)
    for j in range(px):
        for i in range(px):
            pixels[i,j] = algo(next(points))
        pb1.update(px)
            
    m, M = pixels.min(), pixels.max()
    pixels[pixels == nan] = M # none or few points
    pixels[pixels == 0] = m if ic == 'continuous' else M
    normed = Normalize(m, M)(pixels)
    
    if cpc is not None:
        p, pf, pc = cpc
        q = percentile(normed, p)
        filt = normed <= q
        normed[filt] = pow(normed[filt], pf)
        filt = ~filt
        normed[filt] = pow(normed[filt], pc)
    elif power != 1:
        normed = pow(normed, power)

    colors = colormap(normed)

    ppm.writelines(['P3\n', f'{px} {px}\n', '255\n'])
    for j in range(px):
        ppm.write('\n')
        for i in range(px):
            rgb = map(round, 255 * colors[i][j][:3])
            ppm.write(' '.join(map(str, rgb)) + '  ')
        pb2.update(px)
    
          
if __name__ == '__main__':

    from pathlib import Path
    from enum import Enum
    from typing import Tuple
    from typer import Typer, Option

    img_dir = Path('img/')
    img_dir.mkdir(exist_ok=True)

    app = Typer()

    class Algorithm(Enum): 
        escapetime = 'escapetime'
        DEM = 'DEM'
    
    class ColorMapOrder(Enum):
        normal = 'normal'
        reversed = 'reversed'
    
    class InteriorColor(Enum):
        continuous = 'continuous'
        inverted = 'inverted'
    
    def checkConditions(center, radius, cpc, cmp, **kwargs):
        if center != 0 and radius is None:
            raise Exception(
                'if center is non-trivial, '
                'then radius must be provided'
            )
        if cpc and cmp != 1: 
            raise Exception(
                'color map percentile and '
                'colormap power cannot be both provided'
            )
    
    def processAndDraw(
        algo, radius, center, px, fn,
        cmap, cmo, cmp, cpc, ic, 
        **kwargs
    ):
        if fn == '': fn = ' '.join(argv)
        filepath = img_dir / Path(fn + '.ppm')
        filepath.touch()

        colormap = cm.get_cmap(cmap)
        if cmo == 'reversed': 
            colormap = colormap.reversed()
        
        with (
            open(filepath, 'w', encoding='utf-8') as ppm,
            tqdm(total=px*px) as pb1,
            tqdm(total=px*px) as pb2
        ):
            drawFractalPPM(
                ppm, algo, complex(center), radius, px, 
                colormap, cmp, cpc, ic.name, pb1, pb2
            )

    @app.command('julia')
    def julia(
        polynomial: str,
        center: str = Option("0", '-c'),
        radius: float = Option(0, '-r'),
        px: int = Option(1000, '-px'),
        fn: str = Option('', '-fn'),
        it: int = Option(250, '-it'),
        alg: Algorithm = Option('DEM', '-alg'),
        cmap: str = Option('viridis', '-cm'),
        cmo: ColorMapOrder = Option('normal', '-cmo'),
        cmp: float = Option(1, '-cmp'),
        cpc: Tuple[int,float,float] = Option(None, '-cpc'),
        ic: InteriorColor = Option('continuous', '-ic')
    ):
        checkConditions(**locals())

        p = list(map(complex, polynomial.split()))
        R = juliaRadius(p)
        if not radius: 
            radius = simulateRadius(p, R, px)
        
        if alg.name == 'DEM':
            dp = differentiate(p)
            algo = lambda z: demJulia(z, p, dp, it, R)
        if alg.name == 'escapetime':
            R2 = R * R
            algo = lambda z: escapetimeJulia(z, p, it, R2)

        processAndDraw(**locals())

    @app.command('mandelbrot')
    def mandelbrot(
        center: str = Option("-0.8", '-c'),
        radius: float = Option(1.4, '-r'),
        px: int = Option(1000, '-px'),
        fn: str = Option('', '-fn'),
        it: int = Option(250, '-it'),
        alg: Algorithm = Option('DEM', '-alg'),
        cmap: str = Option('viridis', '-cm'),
        cmo: ColorMapOrder = Option('normal', '-cmo'),
        cmp: float = Option(1, '-cmp'),
        cpc: Tuple[int,float,float] = Option(None, '-cpc'),
        ic: InteriorColor = Option('continuous', '-ic')
    ):
        checkConditions(**locals())

        if alg.name == 'DEM': 
            algo = lambda c: demMandelbrot(c, it)
        if alg.name == 'escapetime': 
            algo = lambda c: escapetimeMandelbrot(c, it)
        
        processAndDraw(**locals())

    app()
