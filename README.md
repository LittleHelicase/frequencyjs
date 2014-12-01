# frequencyjs

A library for computing the frequencies of signals. It will include 
*Fourier Transformation* (under development) and *Wavelet Transformation* with
different basis functions. It also includes a very simple tone generator.

# Usage

frequencyjs consists of the following packages

 - Generators
 - Transformators

## Generators

To generate a signal simply do

```
var fjs = require("frequencyjs");
fjs.Generator
  .sine({frequency:440})  // there will be different types here
  .create({length: 2200, sampling: 2200})
  .data() // data will give you a list of objects [{t:0,value:0},...]
```

## Transformators

To analyse a signal you can use the DFT on any signal:

```
// sampling = samples per second
var spectrum = fjs.Transform
  .toSpectrum([1,0,-1,0,1,0,-1,0],{sampling: 8, method: 'dft'});

// get the whole spectrum [{frequency: 0, amplitude: 0},...]
var spectrumArray = spectrum.data();

// get the dominant frequency { frequency: 2, amplitude: 1 }
var dominantFrequency = spectrum.dominantFrequency();
```

To switch to FFT you must keep in mind that the FFT only allows for signals with a power of two length. Simply switch
the method to `fft` and you get your FFT.

```
// sampling = samples per second
var spectrum = fjs.Transform
  .toSpectrum([1,0,-1,0,1,0,-1,0],{sampling: 8, method: 'fft'});

// get the whole spectrum [{frequency: 0, amplitude: 0},...]
var spectrumArray = spectrum.data();

// get the dominant frequency { frequency: 2, amplitude: 1 }
var dominantFrequency = spectrum.dominantFrequency();
```


# Acknowledgement

A part of this library is based on [dsp.js](https://github.com/corbanbrook/dsp.js)
which has a very nice and fast implementation of the FFT. It currently uses
a huge portion of the **dsp.js** code.

# License

This library is under the MIT License.
