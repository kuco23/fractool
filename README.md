# Fractool

This repo features a CLI app implementing the fractal drawing process via two algorithms - escapetime and DEM. Images are generated in the newly created `img` folder. 
> **Note:** 
> This is for learning purposes, the code is not optimized and is written in Python, which means image rendering is slow.

## CLI app usage

Implementation is in the file `fractool.py`. It enables you to draw two types of fractals - Julia sets and the Mandelbrot set.
It accepts multiple arguments that help with fractal image creation. Two examples:

```bash
python fractool.py julia "1 0 -0.760+0.0698j" \
-px 2000 -it 500 -cm cubehelix_r -cmp 0.8 -ext .png
```

```bash
python fractool.py mandelbrot \
-px 2000 -cmp 0.7 -fn mandelbrot_img -ext .jpg --cache
```

Below is the list of all possible arguments, with explanations:
- **type:** The first argument is either `julia` or `mandelbrot` depending on the fractal you wish to draw.
- **polynomial:** If you chose `julia` as the first argument, you must also specify the associated polynomial. Those are represented as strings, e.g. `"1 2 3"` means the polynomial `1*z^2 + 2*z + 3`.
- **center:** This specifies the center of your image. It defaults to 0 and is set by e.g. `-c "1 + 0.5j"`.
- **radius:** This specifies the radius around the center and is set by e.g. `-r 2`. If center is 0, then the radius can be calculated algorithmically. Else, you need to specify it.
- **algorithm:** You can choose the algorithm used for drawing a chosen fractal as `-alg` followed by either `escapetime` or the default `DEM` (Distance estimation method).
- **iterations:** This is the number of iterations that approximate convergance, used by the chosen algorithm. It can be set by `-it` followed by some positive integer (default is `100`).
- **pixels:** This is the square root of the number of pixels of the generated square image. Set by `-px`, followed by a positive integer (default is `1000`).
- **colormap:** You can specify the colormap by adding `-cm` followed by some `matplotlib` colormap name. The default for julia sets is `inferno_r`, while for the mandelbrot set it's `gist_stern_r`. You can check all available options [here](https://matplotlib.org/stable/tutorials/colors/colormaps.html) or `import cm from matplotlib` and doing `dir(cm)`.
- **invert colormap:** You can color the set interior with the outermost colour (defined by the chosen colormap) by typing `-ic` followed by either `continuous` or `inverted`. Default is `continuous`, but `inverted` is appropriate when using escapetime, so the interior is dark like the area far from the set. In this way the bright border is better visible.
- **colormap power:** Sometimes the border is thin and barely visible, because all iteration values returned by the algorithm are either very close to 0 or very close to 1. This can be fixed by applying roots or powers. It can be done by e.g. `-cmp 2`, which squares all values before applying the colormap. So, it brings numbers close to 1 lower and out of hiding.
- **colormap percentage power:** This is a generalization of colormap power. Powering all the values by the same power may bring those close to 1 lower as well as those already close to 0. This is fixed by differentating the powers applied to some high and low values, which are specified by the given percentage. Use as e.g. `-cpc 91 4 0.25`. See examples for more.
- **file name:** To set the filename of the file use `-fn`. By default files are saved by the command that generated them (trust me, that's a good idea!).
- **extension:** To set the image type, set `-ext` followed by some extension. Default is `.png`.
- **cache:** Drawing a fractal image can take a lot of time, so you might save the general fractal-related info by adding `--cache` to your command. Then you can test different colorings without having to again generate the fractal itself. Whenever you have cached data available (in the `data` folder) for a chosen command, the program recognized this and uses that cached data. The data can be reused if the generating command used the same **type**, **algorithm**, **center**, **radius**, **pixels** and **iterations**.

## Examples

To understand `-cpc` command, see the following examples:

```bash
python fractool.py julia "1 0 -0.7510894579318156+0.11771693494277351j" \
-cm cubehelix
```

```bash
python fractool.py julia "1 0 -0.7510894579318156+0.11771693494277351j" \
-cm cubehelix -cpc 91 4 0.25
```

```bash
python fractool.py julia "1 0 0.1567002004882749+0.6527033090669409j" \
-cm magma -it 1000 -cpc 85 2 0.25
```
