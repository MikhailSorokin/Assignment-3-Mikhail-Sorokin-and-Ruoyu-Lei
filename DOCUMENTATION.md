Authors - Mikhail Sorokin and Ruoyu Lei

Programming Assignment 2: Image Processing
----------
#### The following methods are implemented:

1. [Brightness](#brightness)
2. [Contrast](#contrast)
3. [Black and White](#black-and-white)
4. [Gaussian Blur](#gaussian-blur)
5. [Channel Extract](#channel-extract)
6. [Median Filter](#median-filter)
7. [Rotate](#rotate)
8. [Scale](#scale)
9. [Crop](#crop)
10. [Fun](#fun)

### This image will be used to have different effects applied on:
![foo](input/Mountain_side.jpg)

# Brightness
Brightness is implemented by using the interpolation formula provided 

>out = (1 - alpha) * in0 + alpha * in1

where alpha is the factor, in0 is zero alpha image, and in1 is the original image.

In this method specifically, the zero alpha image is a pure black image that has the same width and height as the original image.

The following image is produced by running the command:

```
./cmsc427 input/Mountain_side.jpg output/Brightness.jpg -brightness 0.5
```

![foo](output/Brightness.jpg)

The following image is produced by running the command:

```
./cmsc427 input/Mountain_side.jpg output/Brightness_2.jpg -brightness 2
```

![foo](output/Brightness_2.jpg)

# Contrast
Contrast is also implemented by using this interpolation formula:
>out = (1 - alpha) * in0 + alpha * in1

and in this method in0 is a gray image that has r,g,b values of the average luminance of the original image, calculated by
>luminance = 0.30*r + 0.59*g + 0.11*b

Running this command will give us:
```
./cmsc427 input/Mountain_side.jpg output/contrast_-0.5.jpg -contrast -0.5
```

![foo](output/contrast_-0.5.jpg)

Running this command will give us:
```
./cmsc427 input/Mountain_side.jpg output/contrast_2.jpg -contrast 2
```

![foo](output/Contrast_2.jpg)

# Black and White
Black & White method is implemented by replacing each pixel with its luminance.

Running this command will give us:
```
./cmsc427 input/Mountain_side.jpg output/blackandwhite.jpg -blackandwhite
```

![foo](output/blackandwhite.jpg)

# Gaussian Blur
Gaussian Blue is produced using the gaussian formula and using the standard deviation of neighboring points of every pixel in order to produce a blurring effect on the entire image.

Running this command will give us:
```
./cmsc427 input/Mountain_side.jpg output/gaussianblur_04.jpg -gaussian_blur 0.4
```

![foo](output/gaussianblur_04.jpg)

Running this command will give us:
```
./cmsc427 input/Mountain_side.jpg output/gaussianblur_5.jpg -gaussian_blur 5
```

![foo](output/gaussianblur_5.jpg)

# Channel Extract
Channel extract is implemented by keeping one color channel's value as it is while the other twos' are set to zero.

Running this command will give us:
```
./cmsc427 input/Mountain_side.jpg output/chnext_r.jpg -channel_extract 0
```

![foo](output/chnext_r.jpg)

Running this command will give us:
```
./cmsc427 input/Mountain_side.jpg output/chnext_g.jpg -channel_extract 1
```

![foo](output/chnext_g.jpg)

Running this command will give us:
```
./cmsc427 input/Mountain_side.jpg output/chnext_b.jpg -channel_extract 2
```

![foo](output/chnext_b.jpg)

# Median Filter
Median filter is implemented by using the median of values adjacent to a pixel using the "scale width" formula. Once the median of the pixels are computed, then the image turns out to remove some of the "dirty" or "noisy" pixels by filling them in with the median value (assume odd input width according to Piazza!).

Running this command will give us an example of that action in process:
```
./cmsc427 output/Median_Filter_Test.jpg output/MedWidth.jpg -median_filter 5.0
```

# Rotate

Rotation of image is implemented by applying the formula for every (x,y) in the target image:

>x=ucosΘ-vsinΘ

>y=usinΘ+vcosΘ

where (u,v) is the coordinates in the original image, and Θ is the angle in the range [0,360].

There are three -sampling options you can call alone with -scale. They are: 

>0 Point sampling and gaussian.

>1 Bilinear sampling

>2 Gaussian sampling

Point sampling is also called nearest-neighbor sampling. It's the simplest sampling method because it just naively finds the flooring coordinates from the original image when doing reverse-mapping.

This command rotates the image -using point sampling- by 90 degress and gives the output:

```
./cmsc427 input/Mountain_side.jpg output/rotate_90.jpg -rotate 90 -sampling 0
```

![foo](output/rotate_90.jpg)

This command rotates the image -using bilinear filterin- by 180 degress and gives the output:

```
./cmsc427 input/Mountain_side.jpg output/rotate_180.jpg -rotate 180 -sampling 1
```

![foo](output/rotate_180.jpg)

This command scales down the image by half using Gaussian sampling, which has been described in detail in the Gaussian Blur function.

```
./cmsc427 input/Mountain_side.jpg output/rotate_270.jpg -rotate 270 -sampling 2
```

![foo](output/rotate_270.jpg)

# Scale
Scale is implemented using reverse-mapping. We are given a value in the range of [0.05,20]. If it's small than 1 then the image is scaled down. Similarly, the image is scaled up when the value is greater than 1.

The sampling options are discussed above in rotate.

This command scales down the image by half using point sampling.

```
./cmsc427 input/Mountain_side.jpg output/scale_0.jpg -scale 0.5 0.5 -sampling 0
```

![foo](output/scale_0.jpg)

This command scales down the image by half using bilinear sampling.

```
./cmsc427 input/Mountain_side.jpg output/scale_1.jpg -scale 0.5 0.5 -sampling 1
```

![foo](output/scale_1.jpg)

This command scales down the image by half using Gaussian sampling.

```
./cmsc427 input/Mountain_side.jpg output/scale_2.jpg -scale 0.5 0.5 -sampling 2
```

![foo](output/scale_2.jpg)

# Crop
In crop method, we are given 4 arguments: top_left_x, top_left_y, crop_width and crop_height. It is implemented by copying pixels starting from (top_left_x, top_left_y) all the way to (top_left_x + crop_width, top_left_y + crop_height) as a rectangle.

This commands crops the image from (200,200). The width of the cutting window is 400 and the height is 500.

```
./cmsc427 input/Mountain_side.jpg output/crop_400.jpg -crop 200 200 400 500
```

![foo](output/crop_400.jpg)

# Fun
Looks like a wave swept through the image.

```
./cmsc427 input/cballs.jpg output/so_fun.jpg -fun -sampling 0
```

![foo](output/so_fun.jpg)