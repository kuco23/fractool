import math
from functools import reduce
from matplotlib import cm
from matplotlib.colors import Normalize
import numpy as np
from tqdm import tqdm

L = 1.000001
OVERFLOW = 10**20

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
    if absdk == 0: return np.nan # rarely happens
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
    if absdk == 0: return np.nan # rarely happens
    estimate = math.log2(abszk) * abszk / absdk
    return -math.log2(estimate)


# generates px^2 points on complex plane 
# with radius <= radius and center center
def _complexlattice(center, radius, px):
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

def _algoValues(algo, center, radius, px, pb): 
    points = _complexlattice(center, radius, px)
    pixels = np.empty((px, px), dtype=float)
    for j in range(px):
        for i in range(px):
            pixels[j,i] = algo(next(points))
        pb.update(px)
    return pixels

def _applyColors(values, colormap, cmp, cpc, ic):    
    m, M = values.min(), values.max()
    values[values == np.nan] = M # none or few points
    values[values == 0] = m if ic == 'continuous' else M
    normed = Normalize(m, M)(values)
    
    if cpc is not None:
        p, pf, pc = cpc
        q = np.percentile(normed, p)
        filt = normed <= q
        normed[filt] = pow(normed[filt], pf)
        filt = ~filt
        normed[filt] = pow(normed[filt], pc)
    elif cmp != 1:
        normed = pow(normed, cmp)

    colormat = colormap(normed)
    return (255 * colormat[:,:,:3]).astype(int)


if __name__ == '__main__':

    from sys import argv
    from pathlib import Path
    from enum import Enum
    from typing import Tuple
    from typer import Typer, Option
    from cv2 import imwrite # pip install opencv-python

    img_dir = Path('img/')
    data_dir = Path('data/')
    img_dir.mkdir(exist_ok=True)
    data_dir.mkdir(exist_ok=True)

    app = Typer()

    class Algorithm(Enum): 
        escapetime = 'escapetime'
        DEM = 'DEM'
    
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
    
    def drawFractal(
        algo, radius, center, px, 
        cmap, cmp, cpc, ic, cache,
        filepath, cachepath, **kwargs
    ):
        filepath.touch()
        
        # check if data for image is cached
        if cachepath.exists():
            print('using cached data...')
            values = np.load(cachepath)
        else:
            with tqdm(total=px*px, desc='data') as pb:
                values = _algoValues(algo, center, radius, px, pb)
            if cache:
                cachepath.touch()
                np.save(cachepath, values)

        colormat = _applyColors(values, cmap, cmp, cpc, ic)
        imwrite(str(filepath), colormat[...,::-1])

    def initializeArgs(
        cmap, fn, ext, alg,  
        center, radius, px, it, **kwargs
    ):  
        filename = Path((fn or ' '.join(argv)) + ext)
        if argv[1] == 'mandelbrot': arg = 'mandelbrot'
        if argv[1] == 'julia': arg = f'julia [{argv[2]}]'
        cachename = f'{arg} {alg.name} {center} {radius} {px} {it}.npy' 

        return {
            'center': complex(center),
            'cmap': cm.get_cmap(cmap), 
            'filepath': img_dir / filename, 
            'cachepath': data_dir / cachename
        }

    @app.command('julia')
    def julia(
        polynomial: str,
        center: str = Option("0", '-c'),
        radius: float = Option(0, '-r'),
        px: int = Option(1000, '-px'),
        fn: str = Option('', '-fn'),
        ext: str = Option('.png', '-ext'),
        it: int = Option(250, '-it'),
        alg: Algorithm = Option('DEM', '-alg'),
        cmap: str = Option('inferno', '-cm'),
        cmp: float = Option(1, '-cmp'),
        cpc: Tuple[int,float,float] = Option(None, '-cpc'),
        ic: InteriorColor = Option('continuous', '-ic'),
        cache: bool = False
    ):
        checkConditions(**locals())

        p = list(map(complex, polynomial.split()))
        R = juliaRadius(p)
        
        if not radius:
            radius = R
            print(f'using radius {R}')
        
        if alg.name == 'DEM':
            dp = differentiate(p)
            algo = lambda z: demJulia(z, p, dp, it, R)
        if alg.name == 'escapetime':
            R2 = R * R
            algo = lambda z: escapetimeJulia(z, p, it, R2)

        args = initializeArgs(**locals())
        drawFractal(**{**locals(), **args})

    @app.command('mandelbrot')
    def mandelbrot(
        center: str = Option("-0.8", '-c'),
        radius: float = Option(1.4, '-r'),
        px: int = Option(1000, '-px'),
        fn: str = Option('', '-fn'),
        ext: str = Option('.png', '-ext'),
        it: int = Option(250, '-it'),
        alg: Algorithm = Option('DEM', '-alg'),
        cmap: str = Option('gist_stern_r', '-cm'),
        cmp: float = Option(1, '-cmp'),
        cpc: Tuple[int,float,float] = Option(None, '-cpc'),
        ic: InteriorColor = Option('continuous', '-ic'),
        cache: bool = False
    ):
        checkConditions(**locals())

        if alg.name == 'DEM': 
            algo = lambda c: demMandelbrot(c, it)
        if alg.name == 'escapetime': 
            algo = lambda c: escapetimeMandelbrot(c, it)
        
        args = initializeArgs(**locals())
        drawFractal(**{**locals(), **args})

    app()
