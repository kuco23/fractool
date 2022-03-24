# fractal
This repo features the code from the [medium post](https://medium.com/p/6ad53bbc8208) and additionally a CLI app implementing the ideas derived there. The files generated are in the PPM format. To deal with this use [GIMP](https://www.gimp.org/). Images will be generated in the newly created `img` folder.

## CLI app usage
Implementation is in the file `fractal.py`. An example usage is

```bash
python fractal.py julia "1 2 3" -cm cubehelix -cmo normal -px 2000
````

Here we wanted to produce a Julia set of a polynomial defined **z*z + 2 + 3** with the [colormap](https://matplotlib.org/stable/tutorials/colors/colormaps.html) cubehelix non-reversed (colormap order is normal), using 2000^2 pixels. If we want to draw the Mandelbrot set, then instead of the polynomial we pass the string `mandelbrot` e.g.

``` bash
python fractal.py mandelbrot -cm gist_stern -it 200
```

Here `-it 200` implies that 200 iterations will be used in the algorithm. There are also two possible algorithms: escapetime and DEM (Distance Estimation Method). The default is DEM, but escapetime can be used as

```bash
python fractal.py mandelbrot -alg escapetime -cm gist_stern -it 100 -ic inverted
```

Here `-ic inverted` implies that the colour of the set interior will be coloured with the outermost colour (defined inside the given colourmap - in this case `gist_stern`). This is done mostly when using escapetime, so the interior is dark like the area far from the set. In this way the bright border is better visible. More advance usage can be seen in this two cases

```bash
python fractal.py julia "1 0 -0.7510894579318156+0.11771693494277351j" 
-cm cubehelix -cmo normal
```

```bash
python fractal.py julia "1 0 -0.7510894579318156+0.11771693494277351j" 
-cm cubehelix -cmo normal -cpc 91 4 0.25
```

Here we use `-cpc` argument to specify that we want to take the upper 91% of values and apply x -> power(x, 4) on them. The other 9% will get applied x -> power(x, 1/4). In this way we can increase the contrast between the points very close to interior and those a little further. Another example is

```bash
python fractal.py julia "1 0 0.1567002004882749+0.6527033090669409j" 
-cm magma -it 1000 -cmo normal -cpc 85 2 0.25
```

Another argument is `-cmp`, which takes a positive real q and applies x -> power(x, q) to all values (less general then `-cpc`).

We can also zoom into specific points by specifying centre with `-c` and the radius around the centre with `-r`. Note that if nothing is passed, the centre point will be 0 and the radius will be determined automatically. To set the filename of the ppm file use `-fn`. By default the files are saved by the command that generated them (you'll eventually thank me for this).

## Cache
Getting fractal iteration data can take a lot of time and computer power. But you may want to test a couple of different color configurations and don't want to generate fractal data each time. To solve this, cli takes an additional argument `--cache`, indicating that you want to save the data, used to generate the (non-colored) fractal. This is then saved in `data` folder and each time you run cli on an argument that was used to generate a cached fractal, it will be found in the data folder and used, so you don't have to wait again for all its iterations to complete.