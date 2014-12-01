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

```javascript
var fjs = require("frequencyjs");
fjs.Generator
  .sine({frequency:440})  // there will be different types here
  .create({length: 2200, sampling: 2200})
  .data() // data will give you a list of objects [{t:0,value:0},...]
```

## Transformators

To analyse a signal you can use the DFT on any signal:

```javascript
var spectrum = fjs.Transform
  .toSpectrum([1,0,-1,0,1,0,-1,0],{sampling: 8, method: 'dft'});
```

If you want to use FFT you must keep in mind that the FFT only works on signals with a power of two length. Simply switch the method to `fft` and you get your FFT.

```javascript
var spectrum = fjs.Transform
  .toSpectrum(signal,{method: 'fft'});
```

The signals created with the Generator above contain information about their sampling rate (`signal.sampling`) and thus it is possible to omit the sampling rate.

## Spectrum

A spectrum calculated with the `toSpectrum` method has the following methods:

 - `data()`: returns a list of Frequencies with their amplitude <br>
  e.g `[{frequency:0, amplitude:0},{...},...]`
  The frequencies are always ascending.
 - `dominantFrequency()`: returns the dominant frequency of the signal as an object containing the frequency and the amplitude <br>
  e.g. `{frequency: 2, amplitude: 1}``

# Acknowledgement

A part of this library is based on [dsp.js](https://github.com/corbanbrook/dsp.js)
which has a very nice and fast implementation of the FFT. It currently uses
a huge portion of the **dsp.js** code.

# License

This library is under the MIT License.
