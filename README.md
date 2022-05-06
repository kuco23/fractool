# fractal
This repo features a CLI app implementing the fractal drawing process via two algorithms - escapetime and DEM. The files generated are in the PPM format. To deal with this use [GIMP](https://www.gimp.org/). Images will be generated in the newly created `img` folder.

## CLI app usage
Implementation is in the file `fractal.py`. It allows you to draw two types of fractals - Julia sets and the Mandelbrot set. 
It accepts multiple arguments that help with fractal image creation. Two examples:

```bash
python fractal.py julia "1 0 0.5" 
-cm cubehelix -px 2000 -fn julia_img -ext .png --cached
```
```bash
python fractal.py mandelbrot 
-cm gist_stern_r -it 200 -alg DEM -fn mandelbrot_img -ext .jpg
```

### CLI arguments

- **type:** The first argument is either `julia` or `mandelbrot` depending on which type of fractal you want to draw.
- **polynomial:** If you chose `julia` as the first argument, then you must also specify the associated polynomial. This is specified as a string, e.g. `"1 2 3"` which means the polynomial `1*z^2 + 2*z + 3`.
- **center:** This is the center around which we focus our image. It defaults to 0 and is set by e.g. `-c "1 + 0.5j"`.
- **radius:** This is the radius around the center, set by e.g. `-r 2`. If center is 0, then the radius can be calculated algorithmically. Else, you have to specify it.
- **colormap:** You can specify the colormap by adding `-cm` following by some `matplotlib` colormap name. You can check the available options [here](https://matplotlib.org/stable/tutorials/colors/colormaps.html) or `import cm from matplotlib` and doing `dir(cm)`.
- **algorithm:** You can specify the algorithm that is used for drawing a chosen fractal by using `-alg` following by either `escapetime` or `DEM` (Distance estimation method). 
- **iterations:** This is the number of iterations that approximate convergance, used by the chosen algorithm. This can be set by `-it`.
- **invert colormap:** You can imply that the set interior will be coloured with the outermost colour (defined by the chosen colormap) by doing `-ic` followed by either `continuous` or `inverted`. This is done mostly when using escapetime, so the interior is dark like the area far from the set. In this way the bright border is better visible.
- **colormap power:** Sometimes the border is thin and barely visible, because all iteration values returned by the algorithm are either very close to 0 or very close to 1. This can be fixed by applying roots or powers. It can be done by e.g. `-cmp 2`, which  squares all values before applying the colormap. So, it brings numbers close to 1 lower and out of hiding.
- **colormap percentage power:** This is a generalization of colormap power. Powering all the values by the same power may bring those close to 1 lower as well as those already close to 0. This is fixed by differentating the powers applied to some high and low values, which are specified by the given percentage. Use as e.g. `-cpc 91 4 0.25`. See examples for more.
- **file name:** To set the filename of the file use `-fn`. By default files are saved by the command that generated them (you'll eventually thank me for this).
- **extension:** To set the image type, set `-ext` following by some extension. Default is `.png`.
- **cache:** Drawing a fractal image can take a lot of time, so you might save that raw data by adding `--cached` to your command. Then you can test different colorings without having to again generate the fractal data. Whenever you have cached data available (in the `data` folder) for a chosen command, the program recognized this and uses that cached data.


### Examples

To understand `-cpc` command see the following examples:
```bash
python fractal.py julia "1 0 -0.7510894579318156+0.11771693494277351j" 
-cm cubehelix
```
```bash
python fractal.py julia "1 0 -0.7510894579318156+0.11771693494277351j" 
-cm cubehelix -cpc 91 4 0.25
```
```bash
python fractal.py julia "1 0 0.1567002004882749+0.6527033090669409j" 
-cm magma -it 1000 -cpc 85 2 0.25
```