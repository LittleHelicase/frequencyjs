(function e(t,n,r){function s(o,u){if(!n[o]){if(!t[o]){var a=typeof require=="function"&&require;if(!u&&a)return a(o,!0);if(i)return i(o,!0);throw new Error("Cannot find module '"+o+"'")}var f=n[o]={exports:{}};t[o][0].call(f.exports,function(e){var n=t[o][1][e];return s(n?n:e)},f,f.exports,e,t,n,r)}return n[o].exports}var i=typeof require=="function"&&require;for(var o=0;o<r.length;o++)s(r[o]);return s})({1:[function(require,module,exports){
/* 
 *  DSP.js - a comprehensive digital signal processing  library for javascript
 * 
 *  Created by Corban Brook <corbanbrook@gmail.com> on 2010-01-01.
 *  Copyright 2010 Corban Brook. All rights reserved.
 *
 */

////////////////////////////////////////////////////////////////////////////////
//                                  CONSTANTS                                 //
////////////////////////////////////////////////////////////////////////////////

/**
 * DSP is an object which contains general purpose utility functions and constants
 */
var DSP = {
  // Channels
  LEFT:           0,
  RIGHT:          1,
  MIX:            2,

  // Waveforms
  SINE:           1,
  TRIANGLE:       2,
  SAW:            3,
  SQUARE:         4,

  // Filters
  LOWPASS:        0,
  HIGHPASS:       1,
  BANDPASS:       2,
  NOTCH:          3,

  // Window functions
  BARTLETT:       1,
  BARTLETTHANN:   2,
  BLACKMAN:       3,
  COSINE:         4,
  GAUSS:          5,
  HAMMING:        6,
  HANN:           7,
  LANCZOS:        8,
  RECTANGULAR:    9,
  TRIANGULAR:     10,

  // Loop modes
  OFF:            0,
  FW:             1,
  BW:             2,
  FWBW:           3,

  // Math
  TWO_PI:         2*Math.PI
};

// Setup arrays for platforms which do not support byte arrays
function setupTypedArray(name, fallback) {
  // check if TypedArray exists
  // typeof on Minefield and Chrome return function, typeof on Webkit returns object.
  if (typeof this[name] !== "function" && typeof this[name] !== "object") {
    // nope.. check if WebGLArray exists
    if (typeof this[fallback] === "function" && typeof this[fallback] !== "object") {
      this[name] = this[fallback];
    } else {
      // nope.. set as Native JS array
      this[name] = function(obj) {
        if (obj instanceof Array) {
          return obj;
        } else if (typeof obj === "number") {
          return new Array(obj);
        }
      };
    }
  }
}

setupTypedArray("Float32Array", "WebGLFloatArray");
setupTypedArray("Int32Array",   "WebGLIntArray");
setupTypedArray("Uint16Array",  "WebGLUnsignedShortArray");
setupTypedArray("Uint8Array",   "WebGLUnsignedByteArray");


////////////////////////////////////////////////////////////////////////////////
//                            DSP UTILITY FUNCTIONS                           //
////////////////////////////////////////////////////////////////////////////////

/**
 * Inverts the phase of a signal
 *
 * @param {Array} buffer A sample buffer
 *
 * @returns The inverted sample buffer
 */
DSP.invert = function(buffer) {
  for (var i = 0, len = buffer.length; i < len; i++) {
    buffer[i] *= -1;
  }

  return buffer;
};

/**
 * Converts split-stereo (dual mono) sample buffers into a stereo interleaved sample buffer
 *
 * @param {Array} left  A sample buffer
 * @param {Array} right A sample buffer
 *
 * @returns The stereo interleaved buffer
 */
DSP.interleave = function(left, right) {
  if (left.length !== right.length) {
    throw "Can not interleave. Channel lengths differ.";
  }
 
  var stereoInterleaved = new Float32Array(left.length * 2);
 
  for (var i = 0, len = left.length; i < len; i++) {
    stereoInterleaved[2*i]   = left[i];
    stereoInterleaved[2*i+1] = right[i];
  }
 
  return stereoInterleaved;
};

/**
 * Converts a stereo-interleaved sample buffer into split-stereo (dual mono) sample buffers
 *
 * @param {Array} buffer A stereo-interleaved sample buffer
 *
 * @returns an Array containing left and right channels
 */
DSP.deinterleave = (function() {
  var left, right, mix, deinterleaveChannel = []; 

  deinterleaveChannel[DSP.MIX] = function(buffer) {
    for (var i = 0, len = buffer.length/2; i < len; i++) {
      mix[i] = (buffer[2*i] + buffer[2*i+1]) / 2;
    }
    return mix;
  };

  deinterleaveChannel[DSP.LEFT] = function(buffer) {
    for (var i = 0, len = buffer.length/2; i < len; i++) {
      left[i]  = buffer[2*i];
    }
    return left;
  };

  deinterleaveChannel[DSP.RIGHT] = function(buffer) {
    for (var i = 0, len = buffer.length/2; i < len; i++) {
      right[i]  = buffer[2*i+1];
    }
    return right;
  };

  return function(channel, buffer) { 
    left  = left  || new Float32Array(buffer.length/2);
    right = right || new Float32Array(buffer.length/2);
    mix   = mix   || new Float32Array(buffer.length/2);

    if (buffer.length/2 !== left.length) {
      left  = new Float32Array(buffer.length/2);
      right = new Float32Array(buffer.length/2);
      mix   = new Float32Array(buffer.length/2);
    }

    return deinterleaveChannel[channel](buffer);
  };
}());

/**
 * Separates a channel from a stereo-interleaved sample buffer
 *
 * @param {Array}  buffer A stereo-interleaved sample buffer
 * @param {Number} channel A channel constant (LEFT, RIGHT, MIX)
 *
 * @returns an Array containing a signal mono sample buffer
 */
DSP.getChannel = DSP.deinterleave;

/**
 * Helper method (for Reverb) to mix two (interleaved) samplebuffers. It's possible
 * to negate the second buffer while mixing and to perform a volume correction
 * on the final signal.
 *
 * @param {Array} sampleBuffer1 Array containing Float values or a Float32Array
 * @param {Array} sampleBuffer2 Array containing Float values or a Float32Array
 * @param {Boolean} negate When true inverts/flips the audio signal
 * @param {Number} volumeCorrection When you add multiple sample buffers, use this to tame your signal ;)
 *
 * @returns A new Float32Array interleaved buffer.
 */
DSP.mixSampleBuffers = function(sampleBuffer1, sampleBuffer2, negate, volumeCorrection){
  var outputSamples = new Float32Array(sampleBuffer1);

  for(var i = 0; i<sampleBuffer1.length; i++){
    outputSamples[i] += (negate ? -sampleBuffer2[i] : sampleBuffer2[i]) / volumeCorrection;
  }
 
  return outputSamples;
}; 

// Biquad filter types
DSP.LPF = 0;                // H(s) = 1 / (s^2 + s/Q + 1)
DSP.HPF = 1;                // H(s) = s^2 / (s^2 + s/Q + 1)
DSP.BPF_CONSTANT_SKIRT = 2; // H(s) = s / (s^2 + s/Q + 1)  (constant skirt gain, peak gain = Q)
DSP.BPF_CONSTANT_PEAK = 3;  // H(s) = (s/Q) / (s^2 + s/Q + 1)      (constant 0 dB peak gain)
DSP.NOTCH = 4;              // H(s) = (s^2 + 1) / (s^2 + s/Q + 1)
DSP.APF = 5;                // H(s) = (s^2 - s/Q + 1) / (s^2 + s/Q + 1)
DSP.PEAKING_EQ = 6;         // H(s) = (s^2 + s*(A/Q) + 1) / (s^2 + s/(A*Q) + 1)
DSP.LOW_SHELF = 7;          // H(s) = A * (s^2 + (sqrt(A)/Q)*s + A)/(A*s^2 + (sqrt(A)/Q)*s + 1)
DSP.HIGH_SHELF = 8;         // H(s) = A * (A*s^2 + (sqrt(A)/Q)*s + 1)/(s^2 + (sqrt(A)/Q)*s + A)

// Biquad filter parameter types
DSP.Q = 1;
DSP.BW = 2; // SHARED with BACKWARDS LOOP MODE
DSP.S = 3;

// Find RMS of signal
DSP.RMS = function(buffer) {
  var total = 0;
  
  for (var i = 0, n = buffer.length; i < n; i++) {
    total += buffer[i] * buffer[i];
  }
  
  return Math.sqrt(total / n);
};

// Find Peak of signal
DSP.Peak = function(buffer) {
  var peak = 0;
  
  for (var i = 0, n = buffer.length; i < n; i++) {
    peak = (Math.abs(buffer[i]) > peak) ? Math.abs(buffer[i]) : peak; 
  }
  
  return peak;
};

// Fourier Transform Module used by DFT, FFT, RFFT
function FourierTransform(bufferSize, sampleRate) {
  this.bufferSize = bufferSize;
  this.sampleRate = sampleRate;
  this.bandwidth  = 2 / bufferSize * sampleRate / 2;

  this.spectrum   = new Float32Array(bufferSize/2);
  this.real       = new Float32Array(bufferSize);
  this.imag       = new Float32Array(bufferSize);

  this.peakBand   = 0;
  this.peak       = 0;

  /**
   * Calculates the *middle* frequency of an FFT band.
   *
   * @param {Number} index The index of the FFT band.
   *
   * @returns The middle frequency in Hz.
   */
  this.getBandFrequency = function(index) {
    return this.bandwidth * index + this.bandwidth / 2;
  };

  this.calculateSpectrum = function() {
    var spectrum  = this.spectrum,
        real      = this.real,
        imag      = this.imag,
        bSi       = 2 / this.bufferSize,
        sqrt      = Math.sqrt,
        rval, 
        ival,
        mag;

    for (var i = 0, N = bufferSize/2; i < N; i++) {
      rval = real[i];
      ival = imag[i];
      mag = bSi * sqrt(rval * rval + ival * ival);

      if (mag > this.peak) {
        this.peakBand = i;
        this.peak = mag;
      }

      spectrum[i] = mag;
    }
  };
}

/**
 * DFT is a class for calculating the Discrete Fourier Transform of a signal.
 *
 * @param {Number} bufferSize The size of the sample buffer to be computed
 * @param {Number} sampleRate The sampleRate of the buffer (eg. 44100)
 *
 * @constructor
 */
function DFT(bufferSize, sampleRate) {
  FourierTransform.call(this, bufferSize, sampleRate);

  var N = bufferSize/2 * bufferSize;
  var TWO_PI = 2 * Math.PI;

  this.sinTable = new Float32Array(N);
  this.cosTable = new Float32Array(N);

  for (var i = 0; i < N; i++) {
    this.sinTable[i] = Math.sin(i * TWO_PI / bufferSize);
    this.cosTable[i] = Math.cos(i * TWO_PI / bufferSize);
  }
}

/**
 * Performs a forward transform on the sample buffer.
 * Converts a time domain signal to frequency domain spectra.
 *
 * @param {Array} buffer The sample buffer
 *
 * @returns The frequency spectrum array
 */
DFT.prototype.forward = function(buffer) {
  var real = this.real, 
      imag = this.imag,
      rval,
      ival;

  for (var k = 0; k < this.bufferSize/2; k++) {
    rval = 0.0;
    ival = 0.0;

    for (var n = 0; n < buffer.length; n++) {
      rval += this.cosTable[k*n] * buffer[n];
      ival += this.sinTable[k*n] * buffer[n];
    }

    real[k] = rval;
    imag[k] = ival;
  }

  return this.calculateSpectrum();
};


/**
 * FFT is a class for calculating the Discrete Fourier Transform of a signal
 * with the Fast Fourier Transform algorithm.
 *
 * @param {Number} bufferSize The size of the sample buffer to be computed. Must be power of 2
 * @param {Number} sampleRate The sampleRate of the buffer (eg. 44100)
 *
 * @constructor
 */
function FFT(bufferSize, sampleRate) {
  FourierTransform.call(this, bufferSize, sampleRate);
   
  this.reverseTable = new Uint32Array(bufferSize);

  var limit = 1;
  var bit = bufferSize >> 1;

  var i;

  while (limit < bufferSize) {
    for (i = 0; i < limit; i++) {
      this.reverseTable[i + limit] = this.reverseTable[i] + bit;
    }

    limit = limit << 1;
    bit = bit >> 1;
  }

  this.sinTable = new Float32Array(bufferSize);
  this.cosTable = new Float32Array(bufferSize);

  for (i = 0; i < bufferSize; i++) {
    this.sinTable[i] = Math.sin(-Math.PI/i);
    this.cosTable[i] = Math.cos(-Math.PI/i);
  }
}

/**
 * Performs a forward transform on the sample buffer.
 * Converts a time domain signal to frequency domain spectra.
 *
 * @param {Array} buffer The sample buffer. Buffer Length must be power of 2
 *
 * @returns The frequency spectrum array
 */
FFT.prototype.forward = function(buffer) {
  // Locally scope variables for speed up
  var bufferSize      = this.bufferSize,
      cosTable        = this.cosTable,
      sinTable        = this.sinTable,
      reverseTable    = this.reverseTable,
      real            = this.real,
      imag            = this.imag,
      spectrum        = this.spectrum;

  var k = Math.floor(Math.log(bufferSize) / Math.LN2);

  if (Math.pow(2, k) !== bufferSize) { throw "Invalid buffer size, must be a power of 2."; }
  if (bufferSize !== buffer.length)  { throw "Supplied buffer is not the same size as defined FFT. FFT Size: " + bufferSize + " Buffer Size: " + buffer.length; }

  var halfSize = 1,
      phaseShiftStepReal,
      phaseShiftStepImag,
      currentPhaseShiftReal,
      currentPhaseShiftImag,
      off,
      tr,
      ti,
      tmpReal,
      i;

  for (i = 0; i < bufferSize; i++) {
    real[i] = buffer[reverseTable[i]];
    imag[i] = 0;
  }

  while (halfSize < bufferSize) {
    //phaseShiftStepReal = Math.cos(-Math.PI/halfSize);
    //phaseShiftStepImag = Math.sin(-Math.PI/halfSize);
    phaseShiftStepReal = cosTable[halfSize];
    phaseShiftStepImag = sinTable[halfSize];
    
    currentPhaseShiftReal = 1;
    currentPhaseShiftImag = 0;

    for (var fftStep = 0; fftStep < halfSize; fftStep++) {
      i = fftStep;

      while (i < bufferSize) {
        off = i + halfSize;
        tr = (currentPhaseShiftReal * real[off]) - (currentPhaseShiftImag * imag[off]);
        ti = (currentPhaseShiftReal * imag[off]) + (currentPhaseShiftImag * real[off]);

        real[off] = real[i] - tr;
        imag[off] = imag[i] - ti;
        real[i] += tr;
        imag[i] += ti;

        i += halfSize << 1;
      }

      tmpReal = currentPhaseShiftReal;
      currentPhaseShiftReal = (tmpReal * phaseShiftStepReal) - (currentPhaseShiftImag * phaseShiftStepImag);
      currentPhaseShiftImag = (tmpReal * phaseShiftStepImag) + (currentPhaseShiftImag * phaseShiftStepReal);
    }

    halfSize = halfSize << 1;
  }

  return this.calculateSpectrum();
};

FFT.prototype.inverse = function(real, imag) {
  // Locally scope variables for speed up
  var bufferSize      = this.bufferSize,
      cosTable        = this.cosTable,
      sinTable        = this.sinTable,
      reverseTable    = this.reverseTable,
      spectrum        = this.spectrum;
     
      real = real || this.real;
      imag = imag || this.imag;

  var halfSize = 1,
      phaseShiftStepReal,
      phaseShiftStepImag,
      currentPhaseShiftReal,
      currentPhaseShiftImag,
      off,
      tr,
      ti,
      tmpReal,
      i;

  for (i = 0; i < bufferSize; i++) {
    imag[i] *= -1;
  }

  var revReal = new Float32Array(bufferSize);
  var revImag = new Float32Array(bufferSize);
 
  for (i = 0; i < real.length; i++) {
    revReal[i] = real[reverseTable[i]];
    revImag[i] = imag[reverseTable[i]];
  }
 
  real = revReal;
  imag = revImag;

  while (halfSize < bufferSize) {
    phaseShiftStepReal = cosTable[halfSize];
    phaseShiftStepImag = sinTable[halfSize];
    currentPhaseShiftReal = 1;
    currentPhaseShiftImag = 0;

    for (var fftStep = 0; fftStep < halfSize; fftStep++) {
      i = fftStep;

      while (i < bufferSize) {
        off = i + halfSize;
        tr = (currentPhaseShiftReal * real[off]) - (currentPhaseShiftImag * imag[off]);
        ti = (currentPhaseShiftReal * imag[off]) + (currentPhaseShiftImag * real[off]);

        real[off] = real[i] - tr;
        imag[off] = imag[i] - ti;
        real[i] += tr;
        imag[i] += ti;

        i += halfSize << 1;
      }

      tmpReal = currentPhaseShiftReal;
      currentPhaseShiftReal = (tmpReal * phaseShiftStepReal) - (currentPhaseShiftImag * phaseShiftStepImag);
      currentPhaseShiftImag = (tmpReal * phaseShiftStepImag) + (currentPhaseShiftImag * phaseShiftStepReal);
    }

    halfSize = halfSize << 1;
  }

  var buffer = new Float32Array(bufferSize); // this should be reused instead
  for (i = 0; i < bufferSize; i++) {
    buffer[i] = real[i] / bufferSize;
  }

  return buffer;
};

/**
 * RFFT is a class for calculating the Discrete Fourier Transform of a signal
 * with the Fast Fourier Transform algorithm.
 *
 * This method currently only contains a forward transform but is highly optimized.
 *
 * @param {Number} bufferSize The size of the sample buffer to be computed. Must be power of 2
 * @param {Number} sampleRate The sampleRate of the buffer (eg. 44100)
 *
 * @constructor
 */

// lookup tables don't really gain us any speed, but they do increase
// cache footprint, so don't use them in here

// also we don't use sepearate arrays for real/imaginary parts

// this one a little more than twice as fast as the one in FFT
// however I only did the forward transform

// the rest of this was translated from C, see http://www.jjj.de/fxt/
// this is the real split radix FFT

function RFFT(bufferSize, sampleRate) {
  FourierTransform.call(this, bufferSize, sampleRate);

  this.trans = new Float32Array(bufferSize);

  this.reverseTable = new Uint32Array(bufferSize);

  // don't use a lookup table to do the permute, use this instead
  this.reverseBinPermute = function (dest, source) {
    var bufferSize  = this.bufferSize, 
        halfSize    = bufferSize >>> 1, 
        nm1         = bufferSize - 1, 
        i = 1, r = 0, h;

    dest[0] = source[0];

    do {
      r += halfSize;
      dest[i] = source[r];
      dest[r] = source[i];
      
      i++;

      h = halfSize << 1;
      while (h = h >> 1, !((r ^= h) & h));

      if (r >= i) { 
        dest[i]     = source[r]; 
        dest[r]     = source[i];

        dest[nm1-i] = source[nm1-r]; 
        dest[nm1-r] = source[nm1-i];
      }
      i++;
    } while (i < halfSize);
    dest[nm1] = source[nm1];
  };

  this.generateReverseTable = function () {
    var bufferSize  = this.bufferSize, 
        halfSize    = bufferSize >>> 1, 
        nm1         = bufferSize - 1, 
        i = 1, r = 0, h;

    this.reverseTable[0] = 0;

    do {
      r += halfSize;
      
      this.reverseTable[i] = r;
      this.reverseTable[r] = i;

      i++;

      h = halfSize << 1;
      while (h = h >> 1, !((r ^= h) & h));

      if (r >= i) { 
        this.reverseTable[i] = r;
        this.reverseTable[r] = i;

        this.reverseTable[nm1-i] = nm1-r;
        this.reverseTable[nm1-r] = nm1-i;
      }
      i++;
    } while (i < halfSize);

    this.reverseTable[nm1] = nm1;
  };

  this.generateReverseTable();
}


// Ordering of output:
//
// trans[0]     = re[0] (==zero frequency, purely real)
// trans[1]     = re[1]
//             ...
// trans[n/2-1] = re[n/2-1]
// trans[n/2]   = re[n/2]    (==nyquist frequency, purely real)
//
// trans[n/2+1] = im[n/2-1]
// trans[n/2+2] = im[n/2-2]
//             ...
// trans[n-1]   = im[1] 

RFFT.prototype.forward = function(buffer) {
  var n         = this.bufferSize, 
      spectrum  = this.spectrum,
      x         = this.trans, 
      TWO_PI    = 2*Math.PI,
      sqrt      = Math.sqrt,
      i         = n >>> 1,
      bSi       = 2 / n,
      n2, n4, n8, nn, 
      t1, t2, t3, t4, 
      i1, i2, i3, i4, i5, i6, i7, i8, 
      st1, cc1, ss1, cc3, ss3,
      e, 
      a,
      rval, ival, mag; 

  this.reverseBinPermute(x, buffer);

  /*
  var reverseTable = this.reverseTable;

  for (var k = 0, len = reverseTable.length; k < len; k++) {
    x[k] = buffer[reverseTable[k]];
  }
  */

  for (var ix = 0, id = 4; ix < n; id *= 4) {
    for (var i0 = ix; i0 < n; i0 += id) {
      //sumdiff(x[i0], x[i0+1]); // {a, b}  <--| {a+b, a-b}
      st1 = x[i0] - x[i0+1];
      x[i0] += x[i0+1];
      x[i0+1] = st1;
    } 
    ix = 2*(id-1);
  }

  n2 = 2;
  nn = n >>> 1;

  while((nn = nn >>> 1)) {
    ix = 0;
    n2 = n2 << 1;
    id = n2 << 1;
    n4 = n2 >>> 2;
    n8 = n2 >>> 3;
    do {
      if(n4 !== 1) {
        for(i0 = ix; i0 < n; i0 += id) {
          i1 = i0;
          i2 = i1 + n4;
          i3 = i2 + n4;
          i4 = i3 + n4;
     
          //diffsum3_r(x[i3], x[i4], t1); // {a, b, s} <--| {a, b-a, a+b}
          t1 = x[i3] + x[i4];
          x[i4] -= x[i3];
          //sumdiff3(x[i1], t1, x[i3]);   // {a, b, d} <--| {a+b, b, a-b}
          x[i3] = x[i1] - t1; 
          x[i1] += t1;
     
          i1 += n8;
          i2 += n8;
          i3 += n8;
          i4 += n8;
         
          //sumdiff(x[i3], x[i4], t1, t2); // {s, d}  <--| {a+b, a-b}
          t1 = x[i3] + x[i4];
          t2 = x[i3] - x[i4];
         
          t1 = -t1 * Math.SQRT1_2;
          t2 *= Math.SQRT1_2;
     
          // sumdiff(t1, x[i2], x[i4], x[i3]); // {s, d}  <--| {a+b, a-b}
          st1 = x[i2];
          x[i4] = t1 + st1; 
          x[i3] = t1 - st1;
          
          //sumdiff3(x[i1], t2, x[i2]); // {a, b, d} <--| {a+b, b, a-b}
          x[i2] = x[i1] - t2;
          x[i1] += t2;
        }
      } else {
        for(i0 = ix; i0 < n; i0 += id) {
          i1 = i0;
          i2 = i1 + n4;
          i3 = i2 + n4;
          i4 = i3 + n4;
     
          //diffsum3_r(x[i3], x[i4], t1); // {a, b, s} <--| {a, b-a, a+b}
          t1 = x[i3] + x[i4]; 
          x[i4] -= x[i3];
          
          //sumdiff3(x[i1], t1, x[i3]);   // {a, b, d} <--| {a+b, b, a-b}
          x[i3] = x[i1] - t1; 
          x[i1] += t1;
        }
      }
   
      ix = (id << 1) - n2;
      id = id << 2;
    } while (ix < n);
 
    e = TWO_PI / n2;

    for (var j = 1; j < n8; j++) {
      a = j * e;
      ss1 = Math.sin(a);
      cc1 = Math.cos(a);

      //ss3 = sin(3*a); cc3 = cos(3*a);
      cc3 = 4*cc1*(cc1*cc1-0.75);
      ss3 = 4*ss1*(0.75-ss1*ss1);
   
      ix = 0; id = n2 << 1;
      do {
        for (i0 = ix; i0 < n; i0 += id) {
          i1 = i0 + j;
          i2 = i1 + n4;
          i3 = i2 + n4;
          i4 = i3 + n4;
       
          i5 = i0 + n4 - j;
          i6 = i5 + n4;
          i7 = i6 + n4;
          i8 = i7 + n4;
       
          //cmult(c, s, x, y, &u, &v)
          //cmult(cc1, ss1, x[i7], x[i3], t2, t1); // {u,v} <--| {x*c-y*s, x*s+y*c}
          t2 = x[i7]*cc1 - x[i3]*ss1; 
          t1 = x[i7]*ss1 + x[i3]*cc1;
          
          //cmult(cc3, ss3, x[i8], x[i4], t4, t3);
          t4 = x[i8]*cc3 - x[i4]*ss3; 
          t3 = x[i8]*ss3 + x[i4]*cc3;
       
          //sumdiff(t2, t4);   // {a, b} <--| {a+b, a-b}
          st1 = t2 - t4;
          t2 += t4;
          t4 = st1;
          
          //sumdiff(t2, x[i6], x[i8], x[i3]); // {s, d}  <--| {a+b, a-b}
          //st1 = x[i6]; x[i8] = t2 + st1; x[i3] = t2 - st1;
          x[i8] = t2 + x[i6]; 
          x[i3] = t2 - x[i6];
         
          //sumdiff_r(t1, t3); // {a, b} <--| {a+b, b-a}
          st1 = t3 - t1;
          t1 += t3;
          t3 = st1;
          
          //sumdiff(t3, x[i2], x[i4], x[i7]); // {s, d}  <--| {a+b, a-b}
          //st1 = x[i2]; x[i4] = t3 + st1; x[i7] = t3 - st1;
          x[i4] = t3 + x[i2]; 
          x[i7] = t3 - x[i2];
         
          //sumdiff3(x[i1], t1, x[i6]);   // {a, b, d} <--| {a+b, b, a-b}
          x[i6] = x[i1] - t1; 
          x[i1] += t1;
          
          //diffsum3_r(t4, x[i5], x[i2]); // {a, b, s} <--| {a, b-a, a+b}
          x[i2] = t4 + x[i5]; 
          x[i5] -= t4;
        }
     
        ix = (id << 1) - n2;
        id = id << 2;
   
      } while (ix < n);
    }
  }

  while (--i) {
    rval = x[i];
    ival = x[n-i-1];
    mag = bSi * sqrt(rval * rval + ival * ival);

    if (mag > this.peak) {
      this.peakBand = i;
      this.peak = mag;
    }

    spectrum[i] = mag;
  }

  spectrum[0] = bSi * x[0];

  return spectrum;
};

function Sampler(file, bufferSize, sampleRate, playStart, playEnd, loopStart, loopEnd, loopMode) {
  this.file = file;
  this.bufferSize = bufferSize;
  this.sampleRate = sampleRate;
  this.playStart  = playStart || 0; // 0%
  this.playEnd    = playEnd   || 1; // 100%
  this.loopStart  = loopStart || 0;
  this.loopEnd    = loopEnd   || 1;
  this.loopMode   = loopMode  || DSP.OFF;
  this.loaded     = false;
  this.samples    = [];
  this.signal     = new Float32Array(bufferSize);
  this.frameCount = 0;
  this.envelope   = null;
  this.amplitude  = 1;
  this.rootFrequency = 110; // A2 110
  this.frequency  = 550;
  this.step       = this.frequency / this.rootFrequency;
  this.duration   = 0;
  this.samplesProcessed = 0;
  this.playhead   = 0;
 
  var audio = /* new Audio();*/ document.createElement("AUDIO");
  var self = this;
 
  this.loadSamples = function(event) {
    var buffer = DSP.getChannel(DSP.MIX, event.frameBuffer);
    for ( var i = 0; i < buffer.length; i++) {
      self.samples.push(buffer[i]);
    }
  };
 
  this.loadComplete = function() {
    // convert flexible js array into a fast typed array
    self.samples = new Float32Array(self.samples);
    self.loaded = true;
  };
 
  this.loadMetaData = function() {
    self.duration = audio.duration;
  };
 
  audio.addEventListener("MozAudioAvailable", this.loadSamples, false);
  audio.addEventListener("loadedmetadata", this.loadMetaData, false);
  audio.addEventListener("ended", this.loadComplete, false);
  audio.muted = true;
  audio.src = file;
  audio.play();
}

Sampler.prototype.applyEnvelope = function() {
  this.envelope.process(this.signal);
  return this.signal;
};

Sampler.prototype.generate = function() {
  var frameOffset = this.frameCount * this.bufferSize;
 
  var loopWidth = this.playEnd * this.samples.length - this.playStart * this.samples.length;
  var playStartSamples = this.playStart * this.samples.length; // ie 0.5 -> 50% of the length
  var playEndSamples = this.playEnd * this.samples.length; // ie 0.5 -> 50% of the length
  var offset;

  for ( var i = 0; i < this.bufferSize; i++ ) {
    switch (this.loopMode) {
      case DSP.OFF:
        this.playhead = Math.round(this.samplesProcessed * this.step + playStartSamples);
        if (this.playhead < (this.playEnd * this.samples.length) ) {
          this.signal[i] = this.samples[this.playhead] * this.amplitude;
        } else {
          this.signal[i] = 0;
        }
        break;
     
      case DSP.FW:
        this.playhead = Math.round((this.samplesProcessed * this.step) % loopWidth + playStartSamples);
        if (this.playhead < (this.playEnd * this.samples.length) ) {
          this.signal[i] = this.samples[this.playhead] * this.amplitude;
        }
        break;
       
      case DSP.BW:
        this.playhead = playEndSamples - Math.round((this.samplesProcessed * this.step) % loopWidth);
        if (this.playhead < (this.playEnd * this.samples.length) ) {
          this.signal[i] = this.samples[this.playhead] * this.amplitude;
        }
        break;
       
      case DSP.FWBW:
        if ( Math.floor(this.samplesProcessed * this.step / loopWidth) % 2 === 0 ) {
          this.playhead = Math.round((this.samplesProcessed * this.step) % loopWidth + playStartSamples);
        } else {
          this.playhead = playEndSamples - Math.round((this.samplesProcessed * this.step) % loopWidth);
        }  
        if (this.playhead < (this.playEnd * this.samples.length) ) {
          this.signal[i] = this.samples[this.playhead] * this.amplitude;
        }
        break;
    }
    this.samplesProcessed++;
  }

  this.frameCount++;

  return this.signal;
};

Sampler.prototype.setFreq = function(frequency) {
    var totalProcessed = this.samplesProcessed * this.step;
    this.frequency = frequency;
    this.step = this.frequency / this.rootFrequency;
    this.samplesProcessed = Math.round(totalProcessed/this.step);
};

Sampler.prototype.reset = function() {
  this.samplesProcessed = 0;
  this.playhead = 0;
};

/**
 * Oscillator class for generating and modifying signals
 *
 * @param {Number} type       A waveform constant (eg. DSP.SINE)
 * @param {Number} frequency  Initial frequency of the signal
 * @param {Number} amplitude  Initial amplitude of the signal
 * @param {Number} bufferSize Size of the sample buffer to generate
 * @param {Number} sampleRate The sample rate of the signal
 *
 * @contructor
 */
function Oscillator(type, frequency, amplitude, bufferSize, sampleRate) {
  this.frequency  = frequency;
  this.amplitude  = amplitude;
  this.bufferSize = bufferSize;
  this.sampleRate = sampleRate;
  //this.pulseWidth = pulseWidth;
  this.frameCount = 0;
 
  this.waveTableLength = 2048;

  this.cyclesPerSample = frequency / sampleRate;

  this.signal = new Float32Array(bufferSize);
  this.envelope = null;

  switch(parseInt(type, 10)) {
    case DSP.TRIANGLE:
      this.func = Oscillator.Triangle;
      break;

    case DSP.SAW:
      this.func = Oscillator.Saw;
      break;

    case DSP.SQUARE:
      this.func = Oscillator.Square;
      break;

    default:
    case DSP.SINE:
      this.func = Oscillator.Sine;
      break;
  }

  this.generateWaveTable = function() {
    Oscillator.waveTable[this.func] = new Float32Array(2048);
    var waveTableTime = this.waveTableLength / this.sampleRate;
    var waveTableHz = 1 / waveTableTime;

    for (var i = 0; i < this.waveTableLength; i++) {
      Oscillator.waveTable[this.func][i] = this.func(i * waveTableHz/this.sampleRate);
    }
  };

  if ( typeof Oscillator.waveTable === 'undefined' ) {
    Oscillator.waveTable = {};
  }

  if ( typeof Oscillator.waveTable[this.func] === 'undefined' ) {
    this.generateWaveTable();
  }
 
  this.waveTable = Oscillator.waveTable[this.func];
}

/**
 * Set the amplitude of the signal
 *
 * @param {Number} amplitude The amplitude of the signal (between 0 and 1)
 */
Oscillator.prototype.setAmp = function(amplitude) {
  if (amplitude >= 0 && amplitude <= 1) {
    this.amplitude = amplitude;
  } else {
    throw "Amplitude out of range (0..1).";
  }
};
  
/**
 * Set the frequency of the signal
 *
 * @param {Number} frequency The frequency of the signal
 */  
Oscillator.prototype.setFreq = function(frequency) {
  this.frequency = frequency;
  this.cyclesPerSample = frequency / this.sampleRate;
};
     
// Add an oscillator
Oscillator.prototype.add = function(oscillator) {
  for ( var i = 0; i < this.bufferSize; i++ ) {
    //this.signal[i] += oscillator.valueAt(i);
    this.signal[i] += oscillator.signal[i];
  }
 
  return this.signal;
};
     
// Add a signal to the current generated osc signal
Oscillator.prototype.addSignal = function(signal) {
  for ( var i = 0; i < signal.length; i++ ) {
    if ( i >= this.bufferSize ) {
      break;
    }
    this.signal[i] += signal[i];
   
    /*
    // Constrain amplitude
    if ( this.signal[i] > 1 ) {
      this.signal[i] = 1;
    } else if ( this.signal[i] < -1 ) {
      this.signal[i] = -1;
    }
    */
  }
  return this.signal;
};
     
// Add an envelope to the oscillator
Oscillator.prototype.addEnvelope = function(envelope) {
  this.envelope = envelope;
};

Oscillator.prototype.applyEnvelope = function() {
  this.envelope.process(this.signal);
};
     
Oscillator.prototype.valueAt = function(offset) {
  return this.waveTable[offset % this.waveTableLength];
};
     
Oscillator.prototype.generate = function() {
  var frameOffset = this.frameCount * this.bufferSize;
  var step = this.waveTableLength * this.frequency / this.sampleRate;
  var offset;

  for ( var i = 0; i < this.bufferSize; i++ ) {
    //var step = (frameOffset + i) * this.cyclesPerSample % 1;
    //this.signal[i] = this.func(step) * this.amplitude;
    //this.signal[i] = this.valueAt(Math.round((frameOffset + i) * step)) * this.amplitude;
    offset = Math.round((frameOffset + i) * step);
    this.signal[i] = this.waveTable[offset % this.waveTableLength] * this.amplitude;
  }

  this.frameCount++;

  return this.signal;
};

Oscillator.Sine = function(step) {
  return Math.sin(DSP.TWO_PI * step);
};

Oscillator.Square = function(step) {
  return step < 0.5 ? 1 : -1;
};

Oscillator.Saw = function(step) {
  return 2 * (step - Math.round(step));
};

Oscillator.Triangle = function(step) {
  return 1 - 4 * Math.abs(Math.round(step) - step);
};

Oscillator.Pulse = function(step) {
  // stub
};
 
function ADSR(attackLength, decayLength, sustainLevel, sustainLength, releaseLength, sampleRate) {
  this.sampleRate = sampleRate;
  // Length in seconds
  this.attackLength  = attackLength;
  this.decayLength   = decayLength;
  this.sustainLevel  = sustainLevel;
  this.sustainLength = sustainLength;
  this.releaseLength = releaseLength;
  this.sampleRate    = sampleRate;
 
  // Length in samples
  this.attackSamples  = attackLength  * sampleRate;
  this.decaySamples   = decayLength   * sampleRate;
  this.sustainSamples = sustainLength * sampleRate;
  this.releaseSamples = releaseLength * sampleRate;
 
  // Updates the envelope sample positions
  this.update = function() {
    this.attack         =                this.attackSamples;
    this.decay          = this.attack  + this.decaySamples;
    this.sustain        = this.decay   + this.sustainSamples;
    this.release        = this.sustain + this.releaseSamples;
  };
 
  this.update();
 
  this.samplesProcessed = 0;
}

ADSR.prototype.noteOn = function() {
  this.samplesProcessed = 0;
  this.sustainSamples = this.sustainLength * this.sampleRate;
  this.update();
};

// Send a note off when using a sustain of infinity to let the envelope enter the release phase
ADSR.prototype.noteOff = function() {
  this.sustainSamples = this.samplesProcessed - this.decaySamples;
  this.update();
};

ADSR.prototype.processSample = function(sample) {
  var amplitude = 0;

  if ( this.samplesProcessed <= this.attack ) {
    amplitude = 0 + (1 - 0) * ((this.samplesProcessed - 0) / (this.attack - 0));
  } else if ( this.samplesProcessed > this.attack && this.samplesProcessed <= this.decay ) {
    amplitude = 1 + (this.sustainLevel - 1) * ((this.samplesProcessed - this.attack) / (this.decay - this.attack));
  } else if ( this.samplesProcessed > this.decay && this.samplesProcessed <= this.sustain ) {
    amplitude = this.sustainLevel;
  } else if ( this.samplesProcessed > this.sustain && this.samplesProcessed <= this.release ) {
    amplitude = this.sustainLevel + (0 - this.sustainLevel) * ((this.samplesProcessed - this.sustain) / (this.release - this.sustain));
  }
 
  return sample * amplitude;
};

ADSR.prototype.value = function() {
  var amplitude = 0;

  if ( this.samplesProcessed <= this.attack ) {
    amplitude = 0 + (1 - 0) * ((this.samplesProcessed - 0) / (this.attack - 0));
  } else if ( this.samplesProcessed > this.attack && this.samplesProcessed <= this.decay ) {
    amplitude = 1 + (this.sustainLevel - 1) * ((this.samplesProcessed - this.attack) / (this.decay - this.attack));
  } else if ( this.samplesProcessed > this.decay && this.samplesProcessed <= this.sustain ) {
    amplitude = this.sustainLevel;
  } else if ( this.samplesProcessed > this.sustain && this.samplesProcessed <= this.release ) {
    amplitude = this.sustainLevel + (0 - this.sustainLevel) * ((this.samplesProcessed - this.sustain) / (this.release - this.sustain));
  }
 
  return amplitude;
};
     
ADSR.prototype.process = function(buffer) {
  for ( var i = 0; i < buffer.length; i++ ) {
    buffer[i] *= this.value();

    this.samplesProcessed++;
  }
 
  return buffer;
};
     
     
ADSR.prototype.isActive = function() {
  if ( this.samplesProcessed > this.release || this.samplesProcessed === -1 ) {
    return false;
  } else {
    return true;
  }
};

ADSR.prototype.disable = function() {
  this.samplesProcessed = -1;
};
 
function IIRFilter(type, cutoff, resonance, sampleRate) {
  this.sampleRate = sampleRate;

  switch(type) {
    case DSP.LOWPASS:
    case DSP.LP12:
      this.func = new IIRFilter.LP12(cutoff, resonance, sampleRate);
      break;
  }
}

Object.defineProperty(IIRFilter, 'cutoff', {
  get: function cutoff() {
    return this.func.cutoff;
  }
});

Object.defineProperty(IIRFilter, 'resonance', {
  get: function cutoff() {
    return this.func.resonance;
  }
});


IIRFilter.prototype.set = function(cutoff, resonance) {
  this.func.calcCoeff(cutoff, resonance);
};

IIRFilter.prototype.process = function(buffer) {
  this.func.process(buffer);
};

// Add an envelope to the filter
IIRFilter.prototype.addEnvelope = function(envelope) {
  if ( envelope instanceof ADSR ) {
    this.func.addEnvelope(envelope);
  } else {
    throw "Not an envelope.";
  }
};

IIRFilter.LP12 = function(cutoff, resonance, sampleRate) {
  this.sampleRate = sampleRate;
  this.vibraPos   = 0;
  this.vibraSpeed = 0;
  this.envelope = false;
 
  this.calcCoeff = function(cutoff, resonance) {
    this.w = 2.0 * Math.PI * cutoff / this.sampleRate;
    this.q = 1.0 - this.w / (2.0 * (resonance + 0.5 / (1.0 + this.w)) + this.w - 2.0);
    this.r = this.q * this.q;
    this.c = this.r + 1.0 - 2.0 * Math.cos(this.w) * this.q;
   
    this.cutoff = cutoff;
    this.resonance = resonance;
  };

  this.calcCoeff(cutoff, resonance);

  this.process = function(buffer) {
    for ( var i = 0; i < buffer.length; i++ ) {
      this.vibraSpeed += (buffer[i] - this.vibraPos) * this.c;
      this.vibraPos   += this.vibraSpeed;
      this.vibraSpeed *= this.r;
   
      /*
      var temp = this.vibraPos;
     
      if ( temp > 1.0 ) {
        temp = 1.0;
      } else if ( temp < -1.0 ) {
        temp = -1.0;
      } else if ( temp != temp ) {
        temp = 1;
      }
     
      buffer[i] = temp;
      */

      if (this.envelope) {
        buffer[i] = (buffer[i] * (1 - this.envelope.value())) + (this.vibraPos * this.envelope.value());
        this.envelope.samplesProcessed++;
      } else {
        buffer[i] = this.vibraPos;
      }
    }
  };
}; 

IIRFilter.LP12.prototype.addEnvelope = function(envelope) {
  this.envelope = envelope;
};

function IIRFilter2(type, cutoff, resonance, sampleRate) {
  this.type = type;
  this.cutoff = cutoff;
  this.resonance = resonance;
  this.sampleRate = sampleRate;

  this.f = Float32Array(4);
  this.f[0] = 0.0; // lp
  this.f[1] = 0.0; // hp
  this.f[2] = 0.0; // bp
  this.f[3] = 0.0; // br 
 
  this.calcCoeff = function(cutoff, resonance) {
    this.freq = 2 * Math.sin(Math.PI * Math.min(0.25, cutoff/(this.sampleRate*2)));  
    this.damp = Math.min(2 * (1 - Math.pow(resonance, 0.25)), Math.min(2, 2/this.freq - this.freq * 0.5));
  };

  this.calcCoeff(cutoff, resonance);
}

IIRFilter2.prototype.process = function(buffer) {
  var input, output;
  var f = this.f;

  for ( var i = 0; i < buffer.length; i++ ) {
    input = buffer[i];

    // first pass
    f[3] = input - this.damp * f[2];
    f[0] = f[0] + this.freq * f[2];
    f[1] = f[3] - f[0];
    f[2] = this.freq * f[1] + f[2];
    output = 0.5 * f[this.type];

    // second pass
    f[3] = input - this.damp * f[2];
    f[0] = f[0] + this.freq * f[2];
    f[1] = f[3] - f[0];
    f[2] = this.freq * f[1] + f[2];
    output += 0.5 * f[this.type];

    if (this.envelope) {
      buffer[i] = (buffer[i] * (1 - this.envelope.value())) + (output * this.envelope.value());
      this.envelope.samplesProcessed++;
    } else {
      buffer[i] = output;
    }
  }
};

IIRFilter2.prototype.addEnvelope = function(envelope) {
  if ( envelope instanceof ADSR ) {
    this.envelope = envelope;
  } else {
    throw "This is not an envelope.";
  }
};

IIRFilter2.prototype.set = function(cutoff, resonance) {
  this.calcCoeff(cutoff, resonance);
};



function WindowFunction(type, alpha) {
  this.alpha = alpha;
 
  switch(type) {
    case DSP.BARTLETT:
      this.func = WindowFunction.Bartlett;
      break;
     
    case DSP.BARTLETTHANN:
      this.func = WindowFunction.BartlettHann;
      break;
     
    case DSP.BLACKMAN:
      this.func = WindowFunction.Blackman;
      this.alpha = this.alpha || 0.16;
      break;
   
    case DSP.COSINE:
      this.func = WindowFunction.Cosine;
      break;
     
    case DSP.GAUSS:
      this.func = WindowFunction.Gauss;
      this.alpha = this.alpha || 0.25;
      break;
     
    case DSP.HAMMING:
      this.func = WindowFunction.Hamming;
      break;
     
    case DSP.HANN:
      this.func = WindowFunction.Hann;
      break;
   
    case DSP.LANCZOS:
      this.func = WindowFunction.Lanczoz;
      break;
     
    case DSP.RECTANGULAR:
      this.func = WindowFunction.Rectangular;
      break;
     
    case DSP.TRIANGULAR:
      this.func = WindowFunction.Triangular;
      break;
  }
}

WindowFunction.prototype.process = function(buffer) {
  var length = buffer.length;
  for ( var i = 0; i < length; i++ ) {
    buffer[i] *= this.func(length, i, this.alpha);
  }
  return buffer;
};

WindowFunction.Bartlett = function(length, index) {
  return 2 / (length - 1) * ((length - 1) / 2 - Math.abs(index - (length - 1) / 2));
};

WindowFunction.BartlettHann = function(length, index) {
  return 0.62 - 0.48 * Math.abs(index / (length - 1) - 0.5) - 0.38 * Math.cos(DSP.TWO_PI * index / (length - 1));
};

WindowFunction.Blackman = function(length, index, alpha) {
  var a0 = (1 - alpha) / 2;
  var a1 = 0.5;
  var a2 = alpha / 2;

  return a0 - a1 * Math.cos(DSP.TWO_PI * index / (length - 1)) + a2 * Math.cos(4 * Math.PI * index / (length - 1));
};

WindowFunction.Cosine = function(length, index) {
  return Math.cos(Math.PI * index / (length - 1) - Math.PI / 2);
};

WindowFunction.Gauss = function(length, index, alpha) {
  return Math.pow(Math.E, -0.5 * Math.pow((index - (length - 1) / 2) / (alpha * (length - 1) / 2), 2));
};

WindowFunction.Hamming = function(length, index) {
  return 0.54 - 0.46 * Math.cos(DSP.TWO_PI * index / (length - 1));
};

WindowFunction.Hann = function(length, index) {
  return 0.5 * (1 - Math.cos(DSP.TWO_PI * index / (length - 1)));
};

WindowFunction.Lanczos = function(length, index) {
  var x = 2 * index / (length - 1) - 1;
  return Math.sin(Math.PI * x) / (Math.PI * x);
};

WindowFunction.Rectangular = function(length, index) {
  return 1;
};

WindowFunction.Triangular = function(length, index) {
  return 2 / length * (length / 2 - Math.abs(index - (length - 1) / 2));
};

function sinh (arg) {
  // Returns the hyperbolic sine of the number, defined as (exp(number) - exp(-number))/2 
  //
  // version: 1004.2314
  // discuss at: http://phpjs.org/functions/sinh    // +   original by: Onno Marsman
  // *     example 1: sinh(-0.9834330348825909);
  // *     returns 1: -1.1497971402636502
  return (Math.exp(arg) - Math.exp(-arg))/2;
}

/* 
 *  Biquad filter
 * 
 *  Created by Ricard Marxer <email@ricardmarxer.com> on 2010-05-23.
 *  Copyright 2010 Ricard Marxer. All rights reserved.
 *
 */
// Implementation based on:
// http://www.musicdsp.org/files/Audio-EQ-Cookbook.txt
function Biquad(type, sampleRate) {
  this.Fs = sampleRate;
  this.type = type;  // type of the filter
  this.parameterType = DSP.Q; // type of the parameter

  this.x_1_l = 0;
  this.x_2_l = 0;
  this.y_1_l = 0;
  this.y_2_l = 0;

  this.x_1_r = 0;
  this.x_2_r = 0;
  this.y_1_r = 0;
  this.y_2_r = 0;

  this.b0 = 1;
  this.a0 = 1;

  this.b1 = 0;
  this.a1 = 0;

  this.b2 = 0;
  this.a2 = 0;

  this.b0a0 = this.b0 / this.a0;
  this.b1a0 = this.b1 / this.a0;
  this.b2a0 = this.b2 / this.a0;
  this.a1a0 = this.a1 / this.a0;
  this.a2a0 = this.a2 / this.a0;

  this.f0 = 3000;   // "wherever it's happenin', man."  Center Frequency or
                    // Corner Frequency, or shelf midpoint frequency, depending
                    // on which filter type.  The "significant frequency".

  this.dBgain = 12; // used only for peaking and shelving filters

  this.Q = 1;       // the EE kind of definition, except for peakingEQ in which A*Q is
                    // the classic EE Q.  That adjustment in definition was made so that
                    // a boost of N dB followed by a cut of N dB for identical Q and
                    // f0/Fs results in a precisely flat unity gain filter or "wire".

  this.BW = -3;     // the bandwidth in octaves (between -3 dB frequencies for BPF
                    // and notch or between midpoint (dBgain/2) gain frequencies for
                    // peaking EQ

  this.S = 1;       // a "shelf slope" parameter (for shelving EQ only).  When S = 1,
                    // the shelf slope is as steep as it can be and remain monotonically
                    // increasing or decreasing gain with frequency.  The shelf slope, in
                    // dB/octave, remains proportional to S for all other values for a
                    // fixed f0/Fs and dBgain.

  this.coefficients = function() {
    var b = [this.b0, this.b1, this.b2];
    var a = [this.a0, this.a1, this.a2];
    return {b: b, a:a};
  };

  this.setFilterType = function(type) {
    this.type = type;
    this.recalculateCoefficients();
  };

  this.setSampleRate = function(rate) {
    this.Fs = rate;
    this.recalculateCoefficients();
  };

  this.setQ = function(q) {
    this.parameterType = DSP.Q;
    this.Q = Math.max(Math.min(q, 115.0), 0.001);
    this.recalculateCoefficients();
  };

  this.setBW = function(bw) {
    this.parameterType = DSP.BW;
    this.BW = bw;
    this.recalculateCoefficients();
  };

  this.setS = function(s) {
    this.parameterType = DSP.S;
    this.S = Math.max(Math.min(s, 5.0), 0.0001);
    this.recalculateCoefficients();
  };

  this.setF0 = function(freq) {
    this.f0 = freq;
    this.recalculateCoefficients();
  }; 
 
  this.setDbGain = function(g) {
    this.dBgain = g;
    this.recalculateCoefficients();
  };

  this.recalculateCoefficients = function() {
    var A;
    if (type === DSP.PEAKING_EQ || type === DSP.LOW_SHELF || type === DSP.HIGH_SHELF ) {
      A = Math.pow(10, (this.dBgain/40));  // for peaking and shelving EQ filters only
    } else {
      A  = Math.sqrt( Math.pow(10, (this.dBgain/20)) );   
    }

    var w0 = DSP.TWO_PI * this.f0 / this.Fs;

    var cosw0 = Math.cos(w0);
    var sinw0 = Math.sin(w0);

    var alpha = 0;
   
    switch (this.parameterType) {
      case DSP.Q:
        alpha = sinw0/(2*this.Q);
        break;
           
      case DSP.BW:
        alpha = sinw0 * sinh( Math.LN2/2 * this.BW * w0/sinw0 );
        break;

      case DSP.S:
        alpha = sinw0/2 * Math.sqrt( (A + 1/A)*(1/this.S - 1) + 2 );
        break;
    }

    /**
        FYI: The relationship between bandwidth and Q is
             1/Q = 2*sinh(ln(2)/2*BW*w0/sin(w0))     (digital filter w BLT)
        or   1/Q = 2*sinh(ln(2)/2*BW)             (analog filter prototype)

        The relationship between shelf slope and Q is
             1/Q = sqrt((A + 1/A)*(1/S - 1) + 2)
    */

    var coeff;

    switch (this.type) {
      case DSP.LPF:       // H(s) = 1 / (s^2 + s/Q + 1)
        this.b0 =  (1 - cosw0)/2;
        this.b1 =   1 - cosw0;
        this.b2 =  (1 - cosw0)/2;
        this.a0 =   1 + alpha;
        this.a1 =  -2 * cosw0;
        this.a2 =   1 - alpha;
        break;

      case DSP.HPF:       // H(s) = s^2 / (s^2 + s/Q + 1)
        this.b0 =  (1 + cosw0)/2;
        this.b1 = -(1 + cosw0);
        this.b2 =  (1 + cosw0)/2;
        this.a0 =   1 + alpha;
        this.a1 =  -2 * cosw0;
        this.a2 =   1 - alpha;
        break;

      case DSP.BPF_CONSTANT_SKIRT:       // H(s) = s / (s^2 + s/Q + 1)  (constant skirt gain, peak gain = Q)
        this.b0 =   sinw0/2;
        this.b1 =   0;
        this.b2 =  -sinw0/2;
        this.a0 =   1 + alpha;
        this.a1 =  -2*cosw0;
        this.a2 =   1 - alpha;
        break;

      case DSP.BPF_CONSTANT_PEAK:       // H(s) = (s/Q) / (s^2 + s/Q + 1)      (constant 0 dB peak gain)
        this.b0 =   alpha;
        this.b1 =   0;
        this.b2 =  -alpha;
        this.a0 =   1 + alpha;
        this.a1 =  -2*cosw0;
        this.a2 =   1 - alpha;
        break;

      case DSP.NOTCH:     // H(s) = (s^2 + 1) / (s^2 + s/Q + 1)
        this.b0 =   1;
        this.b1 =  -2*cosw0;
        this.b2 =   1;
        this.a0 =   1 + alpha;
        this.a1 =  -2*cosw0;
        this.a2 =   1 - alpha;
        break;

      case DSP.APF:       // H(s) = (s^2 - s/Q + 1) / (s^2 + s/Q + 1)
        this.b0 =   1 - alpha;
        this.b1 =  -2*cosw0;
        this.b2 =   1 + alpha;
        this.a0 =   1 + alpha;
        this.a1 =  -2*cosw0;
        this.a2 =   1 - alpha;
        break;

      case DSP.PEAKING_EQ:  // H(s) = (s^2 + s*(A/Q) + 1) / (s^2 + s/(A*Q) + 1)
        this.b0 =   1 + alpha*A;
        this.b1 =  -2*cosw0;
        this.b2 =   1 - alpha*A;
        this.a0 =   1 + alpha/A;
        this.a1 =  -2*cosw0;
        this.a2 =   1 - alpha/A;
        break;

      case DSP.LOW_SHELF:   // H(s) = A * (s^2 + (sqrt(A)/Q)*s + A)/(A*s^2 + (sqrt(A)/Q)*s + 1)
        coeff = sinw0 * Math.sqrt( (A^2 + 1)*(1/this.S - 1) + 2*A );
        this.b0 =    A*((A+1) - (A-1)*cosw0 + coeff);
        this.b1 =  2*A*((A-1) - (A+1)*cosw0);
        this.b2 =    A*((A+1) - (A-1)*cosw0 - coeff);
        this.a0 =       (A+1) + (A-1)*cosw0 + coeff;
        this.a1 =   -2*((A-1) + (A+1)*cosw0);
        this.a2 =       (A+1) + (A-1)*cosw0 - coeff;
        break;

      case DSP.HIGH_SHELF:   // H(s) = A * (A*s^2 + (sqrt(A)/Q)*s + 1)/(s^2 + (sqrt(A)/Q)*s + A)
        coeff = sinw0 * Math.sqrt( (A^2 + 1)*(1/this.S - 1) + 2*A );
        this.b0 =    A*((A+1) + (A-1)*cosw0 + coeff);
        this.b1 = -2*A*((A-1) + (A+1)*cosw0);
        this.b2 =    A*((A+1) + (A-1)*cosw0 - coeff);
        this.a0 =       (A+1) - (A-1)*cosw0 + coeff;
        this.a1 =    2*((A-1) - (A+1)*cosw0);
        this.a2 =       (A+1) - (A-1)*cosw0 - coeff;
        break;
    }
   
    this.b0a0 = this.b0/this.a0;
    this.b1a0 = this.b1/this.a0;
    this.b2a0 = this.b2/this.a0;
    this.a1a0 = this.a1/this.a0;
    this.a2a0 = this.a2/this.a0;
  };

  this.process = function(buffer) {
      //y[n] = (b0/a0)*x[n] + (b1/a0)*x[n-1] + (b2/a0)*x[n-2]
      //       - (a1/a0)*y[n-1] - (a2/a0)*y[n-2]

      var len = buffer.length;
      var output = new Float32Array(len);

      for ( var i=0; i<buffer.length; i++ ) {
        output[i] = this.b0a0*buffer[i] + this.b1a0*this.x_1_l + this.b2a0*this.x_2_l - this.a1a0*this.y_1_l - this.a2a0*this.y_2_l;
        this.y_2_l = this.y_1_l;
        this.y_1_l = output[i];
        this.x_2_l = this.x_1_l;
        this.x_1_l = buffer[i];
      }

      return output;
  };

  this.processStereo = function(buffer) {
      //y[n] = (b0/a0)*x[n] + (b1/a0)*x[n-1] + (b2/a0)*x[n-2]
      //       - (a1/a0)*y[n-1] - (a2/a0)*y[n-2]

      var len = buffer.length;
      var output = new Float32Array(len);
     
      for (var i = 0; i < len/2; i++) {
        output[2*i] = this.b0a0*buffer[2*i] + this.b1a0*this.x_1_l + this.b2a0*this.x_2_l - this.a1a0*this.y_1_l - this.a2a0*this.y_2_l;
        this.y_2_l = this.y_1_l;
        this.y_1_l = output[2*i];
        this.x_2_l = this.x_1_l;
        this.x_1_l = buffer[2*i];

        output[2*i+1] = this.b0a0*buffer[2*i+1] + this.b1a0*this.x_1_r + this.b2a0*this.x_2_r - this.a1a0*this.y_1_r - this.a2a0*this.y_2_r;
        this.y_2_r = this.y_1_r;
        this.y_1_r = output[2*i+1];
        this.x_2_r = this.x_1_r;
        this.x_1_r = buffer[2*i+1];
      }

      return output;
  };
}

/* 
 *  Magnitude to decibels
 * 
 *  Created by Ricard Marxer <email@ricardmarxer.com> on 2010-05-23.
 *  Copyright 2010 Ricard Marxer. All rights reserved.
 *
 *  @buffer array of magnitudes to convert to decibels
 *
 *  @returns the array in decibels
 *
 */
DSP.mag2db = function(buffer) {
  var minDb = -120;
  var minMag = Math.pow(10.0, minDb / 20.0);

  var log = Math.log;
  var max = Math.max;
 
  var result = Float32Array(buffer.length);
  for (var i=0; i<buffer.length; i++) {
    result[i] = 20.0*log(max(buffer[i], minMag));
  }

  return result;
};

/* 
 *  Frequency response
 * 
 *  Created by Ricard Marxer <email@ricardmarxer.com> on 2010-05-23.
 *  Copyright 2010 Ricard Marxer. All rights reserved.
 *
 *  Calculates the frequency response at the given points.
 *
 *  @b b coefficients of the filter
 *  @a a coefficients of the filter
 *  @w w points (normally between -PI and PI) where to calculate the frequency response
 *
 *  @returns the frequency response in magnitude
 *
 */
DSP.freqz = function(b, a, w) {
  var i, j;

  if (!w) {
    w = Float32Array(200);
    for (i=0;i<w.length; i++) {
      w[i] = DSP.TWO_PI/w.length * i - Math.PI;
    }
  }

  var result = Float32Array(w.length);
 
  var sqrt = Math.sqrt;
  var cos = Math.cos;
  var sin = Math.sin;
 
  for (i=0; i<w.length; i++) {
    var numerator = {real:0.0, imag:0.0};
    for (j=0; j<b.length; j++) {
      numerator.real += b[j] * cos(-j*w[i]);
      numerator.imag += b[j] * sin(-j*w[i]);
    }

    var denominator = {real:0.0, imag:0.0};
    for (j=0; j<a.length; j++) {
      denominator.real += a[j] * cos(-j*w[i]);
      denominator.imag += a[j] * sin(-j*w[i]);
    }
 
    result[i] =  sqrt(numerator.real*numerator.real + numerator.imag*numerator.imag) / sqrt(denominator.real*denominator.real + denominator.imag*denominator.imag);
  }

  return result;
};

/* 
 *  Graphical Equalizer
 *
 *  Implementation of a graphic equalizer with a configurable bands-per-octave
 *  and minimum and maximum frequencies
 * 
 *  Created by Ricard Marxer <email@ricardmarxer.com> on 2010-05-23.
 *  Copyright 2010 Ricard Marxer. All rights reserved.
 *
 */
function GraphicalEq(sampleRate) {
  this.FS = sampleRate;
  this.minFreq = 40.0;
  this.maxFreq = 16000.0;

  this.bandsPerOctave = 1.0;

  this.filters = [];
  this.freqzs = [];

  this.calculateFreqzs = true;

  this.recalculateFilters = function() {
    var bandCount = Math.round(Math.log(this.maxFreq/this.minFreq) * this.bandsPerOctave/ Math.LN2);

    this.filters = [];
    for (var i=0; i<bandCount; i++) {
      var freq = this.minFreq*(Math.pow(2, i/this.bandsPerOctave));
      var newFilter = new Biquad(DSP.PEAKING_EQ, this.FS);
      newFilter.setDbGain(0);
      newFilter.setBW(1/this.bandsPerOctave);
      newFilter.setF0(freq);
      this.filters[i] = newFilter;
      this.recalculateFreqz(i);
    }
  };

  this.setMinimumFrequency = function(freq) {
    this.minFreq = freq;
    this.recalculateFilters();
  };

  this.setMaximumFrequency = function(freq) {
    this.maxFreq = freq;
    this.recalculateFilters();
  };

  this.setBandsPerOctave = function(bands) {
    this.bandsPerOctave = bands;
    this.recalculateFilters();
  };

  this.setBandGain = function(bandIndex, gain) {
    if (bandIndex < 0 || bandIndex > (this.filters.length-1)) {
      throw "The band index of the graphical equalizer is out of bounds.";
    }

    if (!gain) {
      throw "A gain must be passed.";
    }
   
    this.filters[bandIndex].setDbGain(gain);
    this.recalculateFreqz(bandIndex);
  };
 
  this.recalculateFreqz = function(bandIndex) {
    if (!this.calculateFreqzs) {
      return;
    }

    if (bandIndex < 0 || bandIndex > (this.filters.length-1)) {
      throw "The band index of the graphical equalizer is out of bounds. " + bandIndex + " is out of [" + 0 + ", " + this.filters.length-1 + "]";
    }
       
    if (!this.w) {
      this.w = Float32Array(400);
      for (var i=0; i<this.w.length; i++) {
         this.w[i] = Math.PI/this.w.length * i;
      }
    }
   
    var b = [this.filters[bandIndex].b0, this.filters[bandIndex].b1, this.filters[bandIndex].b2];
    var a = [this.filters[bandIndex].a0, this.filters[bandIndex].a1, this.filters[bandIndex].a2];

    this.freqzs[bandIndex] = DSP.mag2db(DSP.freqz(b, a, this.w));
  };

  this.process = function(buffer) {
    var output = buffer;

    for (var i = 0; i < this.filters.length; i++) {
      output = this.filters[i].process(output);
    }

    return output;
  };

  this.processStereo = function(buffer) {
    var output = buffer;

    for (var i = 0; i < this.filters.length; i++) {
      output = this.filters[i].processStereo(output);
    }

    return output;
  };
}

/**
 * MultiDelay effect by Almer Thie (http://code.almeros.com).
 * Copyright 2010 Almer Thie. All rights reserved.
 * Example: http://code.almeros.com/code-examples/delay-firefox-audio-api/
 *
 * This is a delay that feeds it's own delayed signal back into its circular
 * buffer. Also known as a CombFilter.
 *
 * Compatible with interleaved stereo (or more channel) buffers and
 * non-interleaved mono buffers.
 *
 * @param {Number} maxDelayInSamplesSize Maximum possible delay in samples (size of circular buffer)
 * @param {Number} delayInSamples Initial delay in samples
 * @param {Number} masterVolume Initial master volume. Float value: 0.0 (silence), 1.0 (normal), >1.0 (amplify)
 * @param {Number} delayVolume Initial feedback delay volume. Float value: 0.0 (silence), 1.0 (normal), >1.0 (amplify)
 *
 * @constructor
 */
function MultiDelay(maxDelayInSamplesSize, delayInSamples, masterVolume, delayVolume) {
  this.delayBufferSamples   = new Float32Array(maxDelayInSamplesSize); // The maximum size of delay
  this.delayInputPointer     = delayInSamples;
  this.delayOutputPointer   = 0;
 
  this.delayInSamples   = delayInSamples;
  this.masterVolume     = masterVolume;
  this.delayVolume     = delayVolume;
}

/**
 * Change the delay time in samples.
 *
 * @param {Number} delayInSamples Delay in samples
 */
MultiDelay.prototype.setDelayInSamples = function (delayInSamples) {
  this.delayInSamples = delayInSamples;
 
  this.delayInputPointer = this.delayOutputPointer + delayInSamples;

  if (this.delayInputPointer >= this.delayBufferSamples.length-1) {
    this.delayInputPointer = this.delayInputPointer - this.delayBufferSamples.length; 
  }
};

/**
 * Change the master volume.
 *
 * @param {Number} masterVolume Float value: 0.0 (silence), 1.0 (normal), >1.0 (amplify)
 */
MultiDelay.prototype.setMasterVolume = function(masterVolume) {
  this.masterVolume = masterVolume;
};

/**
 * Change the delay feedback volume.
 *
 * @param {Number} delayVolume Float value: 0.0 (silence), 1.0 (normal), >1.0 (amplify)
 */
MultiDelay.prototype.setDelayVolume = function(delayVolume) {
  this.delayVolume = delayVolume;
};

/**
 * Process a given interleaved or mono non-interleaved float value Array and adds the delayed audio.
 *
 * @param {Array} samples Array containing Float values or a Float32Array
 *
 * @returns A new Float32Array interleaved or mono non-interleaved as was fed to this function.
 */
MultiDelay.prototype.process = function(samples) {
  // NB. Make a copy to put in the output samples to return.
  var outputSamples = new Float32Array(samples.length);

  for (var i=0; i<samples.length; i++) {
    // delayBufferSamples could contain initial NULL's, return silence in that case
    var delaySample = (this.delayBufferSamples[this.delayOutputPointer] === null ? 0.0 : this.delayBufferSamples[this.delayOutputPointer]);
   
    // Mix normal audio data with delayed audio
    var sample = (delaySample * this.delayVolume) + samples[i];
   
    // Add audio data with the delay in the delay buffer
    this.delayBufferSamples[this.delayInputPointer] = sample;
   
    // Return the audio with delay mix
    outputSamples[i] = sample * this.masterVolume;
   
    // Manage circulair delay buffer pointers
    this.delayInputPointer++;
    if (this.delayInputPointer >= this.delayBufferSamples.length-1) {
      this.delayInputPointer = 0;
    }
     
    this.delayOutputPointer++;
    if (this.delayOutputPointer >= this.delayBufferSamples.length-1) {
      this.delayOutputPointer = 0; 
    } 
  }
 
  return outputSamples;
};

/**
 * SingleDelay effect by Almer Thie (http://code.almeros.com).
 * Copyright 2010 Almer Thie. All rights reserved.
 * Example: See usage in Reverb class
 *
 * This is a delay that does NOT feeds it's own delayed signal back into its 
 * circular buffer, neither does it return the original signal. Also known as
 * an AllPassFilter(?).
 *
 * Compatible with interleaved stereo (or more channel) buffers and
 * non-interleaved mono buffers.
 *
 * @param {Number} maxDelayInSamplesSize Maximum possible delay in samples (size of circular buffer)
 * @param {Number} delayInSamples Initial delay in samples
 * @param {Number} delayVolume Initial feedback delay volume. Float value: 0.0 (silence), 1.0 (normal), >1.0 (amplify)
 *
 * @constructor
 */

function SingleDelay(maxDelayInSamplesSize, delayInSamples, delayVolume) {
  this.delayBufferSamples = new Float32Array(maxDelayInSamplesSize); // The maximum size of delay
  this.delayInputPointer  = delayInSamples;
  this.delayOutputPointer = 0;
 
  this.delayInSamples     = delayInSamples;
  this.delayVolume        = delayVolume;
}

/**
 * Change the delay time in samples.
 *
 * @param {Number} delayInSamples Delay in samples
 */
SingleDelay.prototype.setDelayInSamples = function(delayInSamples) {
  this.delayInSamples = delayInSamples;
  this.delayInputPointer = this.delayOutputPointer + delayInSamples;

  if (this.delayInputPointer >= this.delayBufferSamples.length-1) {
    this.delayInputPointer = this.delayInputPointer - this.delayBufferSamples.length; 
  }
};

/**
 * Change the return signal volume.
 *
 * @param {Number} delayVolume Float value: 0.0 (silence), 1.0 (normal), >1.0 (amplify)
 */
SingleDelay.prototype.setDelayVolume = function(delayVolume) {
  this.delayVolume = delayVolume;
};

/**
 * Process a given interleaved or mono non-interleaved float value Array and
 * returns the delayed audio.
 *
 * @param {Array} samples Array containing Float values or a Float32Array
 *
 * @returns A new Float32Array interleaved or mono non-interleaved as was fed to this function.
 */
SingleDelay.prototype.process = function(samples) {
  // NB. Make a copy to put in the output samples to return.
  var outputSamples = new Float32Array(samples.length);

  for (var i=0; i<samples.length; i++) {

    // Add audio data with the delay in the delay buffer
    this.delayBufferSamples[this.delayInputPointer] = samples[i];
   
    // delayBufferSamples could contain initial NULL's, return silence in that case
    var delaySample = this.delayBufferSamples[this.delayOutputPointer];

    // Return the audio with delay mix
    outputSamples[i] = delaySample * this.delayVolume;

    // Manage circulair delay buffer pointers
    this.delayInputPointer++;

    if (this.delayInputPointer >= this.delayBufferSamples.length-1) {
      this.delayInputPointer = 0;
    }
     
    this.delayOutputPointer++;

    if (this.delayOutputPointer >= this.delayBufferSamples.length-1) {
      this.delayOutputPointer = 0; 
    } 
  }
 
  return outputSamples;
};

/**
 * Reverb effect by Almer Thie (http://code.almeros.com).
 * Copyright 2010 Almer Thie. All rights reserved.
 * Example: http://code.almeros.com/code-examples/reverb-firefox-audio-api/
 *
 * This reverb consists of 6 SingleDelays, 6 MultiDelays and an IIRFilter2
 * for each of the two stereo channels.
 *
 * Compatible with interleaved stereo buffers only!
 *
 * @param {Number} maxDelayInSamplesSize Maximum possible delay in samples (size of circular buffers)
 * @param {Number} delayInSamples Initial delay in samples for internal (Single/Multi)delays
 * @param {Number} masterVolume Initial master volume. Float value: 0.0 (silence), 1.0 (normal), >1.0 (amplify)
 * @param {Number} mixVolume Initial reverb signal mix volume. Float value: 0.0 (silence), 1.0 (normal), >1.0 (amplify)
 * @param {Number} delayVolume Initial feedback delay volume for internal (Single/Multi)delays. Float value: 0.0 (silence), 1.0 (normal), >1.0 (amplify)
 * @param {Number} dampFrequency Initial low pass filter frequency. 0 to 44100 (depending on your maximum sampling frequency)
 *
 * @constructor
 */
function Reverb(maxDelayInSamplesSize, delayInSamples, masterVolume, mixVolume, delayVolume, dampFrequency) {
  this.delayInSamples   = delayInSamples;
  this.masterVolume     = masterVolume;
  this.mixVolume       = mixVolume;
  this.delayVolume     = delayVolume;
  this.dampFrequency     = dampFrequency;
 
  this.NR_OF_MULTIDELAYS = 6;
  this.NR_OF_SINGLEDELAYS = 6;
 
  this.LOWPASSL = new IIRFilter2(DSP.LOWPASS, dampFrequency, 0, 44100);
  this.LOWPASSR = new IIRFilter2(DSP.LOWPASS, dampFrequency, 0, 44100);
 
  this.singleDelays = [];
  
  var i, delayMultiply;

  for (i = 0; i < this.NR_OF_SINGLEDELAYS; i++) {
    delayMultiply = 1.0 + (i/7.0); // 1.0, 1.1, 1.2...
    this.singleDelays[i] = new SingleDelay(maxDelayInSamplesSize, Math.round(this.delayInSamples * delayMultiply), this.delayVolume);
  }
 
  this.multiDelays = [];

  for (i = 0; i < this.NR_OF_MULTIDELAYS; i++) {
    delayMultiply = 1.0 + (i/10.0); // 1.0, 1.1, 1.2... 
    this.multiDelays[i] = new MultiDelay(maxDelayInSamplesSize, Math.round(this.delayInSamples * delayMultiply), this.masterVolume, this.delayVolume);
  }
}

/**
 * Change the delay time in samples as a base for all delays.
 *
 * @param {Number} delayInSamples Delay in samples
 */
Reverb.prototype.setDelayInSamples = function (delayInSamples){
  this.delayInSamples = delayInSamples;

  var i, delayMultiply;
 
  for (i = 0; i < this.NR_OF_SINGLEDELAYS; i++) {
    delayMultiply = 1.0 + (i/7.0); // 1.0, 1.1, 1.2...
    this.singleDelays[i].setDelayInSamples( Math.round(this.delayInSamples * delayMultiply) );
  }
   
  for (i = 0; i < this.NR_OF_MULTIDELAYS; i++) {
    delayMultiply = 1.0 + (i/10.0); // 1.0, 1.1, 1.2...
    this.multiDelays[i].setDelayInSamples( Math.round(this.delayInSamples * delayMultiply) );
  }
};

/**
 * Change the master volume.
 *
 * @param {Number} masterVolume Float value: 0.0 (silence), 1.0 (normal), >1.0 (amplify)
 */
Reverb.prototype.setMasterVolume = function (masterVolume){
  this.masterVolume = masterVolume;
};

/**
 * Change the reverb signal mix level.
 *
 * @param {Number} mixVolume Float value: 0.0 (silence), 1.0 (normal), >1.0 (amplify)
 */
Reverb.prototype.setMixVolume = function (mixVolume){
  this.mixVolume = mixVolume;
};

/**
 * Change all delays feedback volume.
 *
 * @param {Number} delayVolume Float value: 0.0 (silence), 1.0 (normal), >1.0 (amplify)
 */
Reverb.prototype.setDelayVolume = function (delayVolume){
  this.delayVolume = delayVolume;
 
  var i;

  for (i = 0; i<this.NR_OF_SINGLEDELAYS; i++) {
    this.singleDelays[i].setDelayVolume(this.delayVolume);
  } 
 
  for (i = 0; i<this.NR_OF_MULTIDELAYS; i++) {
    this.multiDelays[i].setDelayVolume(this.delayVolume);
  } 
};

/**
 * Change the Low Pass filter frequency.
 *
 * @param {Number} dampFrequency low pass filter frequency. 0 to 44100 (depending on your maximum sampling frequency)
 */
Reverb.prototype.setDampFrequency = function (dampFrequency){
  this.dampFrequency = dampFrequency;
 
  this.LOWPASSL.set(dampFrequency, 0);
  this.LOWPASSR.set(dampFrequency, 0); 
};

/**
 * Process a given interleaved float value Array and copies and adds the reverb signal.
 *
 * @param {Array} samples Array containing Float values or a Float32Array
 *
 * @returns A new Float32Array interleaved buffer.
 */
Reverb.prototype.process = function (interleavedSamples){ 
  // NB. Make a copy to put in the output samples to return.
  var outputSamples = new Float32Array(interleavedSamples.length);
 
  // Perform low pass on the input samples to mimick damp
  var leftRightMix = DSP.deinterleave(interleavedSamples);
  this.LOWPASSL.process( leftRightMix[DSP.LEFT] );
  this.LOWPASSR.process( leftRightMix[DSP.RIGHT] ); 
  var filteredSamples = DSP.interleave(leftRightMix[DSP.LEFT], leftRightMix[DSP.RIGHT]);

  var i;

  // Process MultiDelays in parallel
  for (i = 0; i<this.NR_OF_MULTIDELAYS; i++) {
    // Invert the signal of every even multiDelay
    outputSamples = DSP.mixSampleBuffers(outputSamples, this.multiDelays[i].process(filteredSamples), 2%i === 0, this.NR_OF_MULTIDELAYS);
  }
 
  // Process SingleDelays in series
  var singleDelaySamples = new Float32Array(outputSamples.length);
  for (i = 0; i<this.NR_OF_SINGLEDELAYS; i++) {
    // Invert the signal of every even singleDelay
    singleDelaySamples = DSP.mixSampleBuffers(singleDelaySamples, this.singleDelays[i].process(outputSamples), 2%i === 0, 1);
  }

  // Apply the volume of the reverb signal
  for (i = 0; i<singleDelaySamples.length; i++) {
    singleDelaySamples[i] *= this.mixVolume;
  }
 
  // Mix the original signal with the reverb signal
  outputSamples = DSP.mixSampleBuffers(singleDelaySamples, interleavedSamples, 0, 1);

  // Apply the master volume to the complete signal
  for (i = 0; i<outputSamples.length; i++) {
    outputSamples[i] *= this.masterVolume;
  }
   
  return outputSamples;
};


module.exports = {
  Oscillator: Oscillator,
  DFT: DFT,
  FFT: FFT,
  SINE: DSP.SINE
}

},{}],2:[function(require,module,exports){

var Generators = require("./lib/generator.js");
var Transform = require("./lib/transform.js");
var Signal = require("./lib/signal.js");
var Processing = require("./lib/processing.js");

module.exports = {
  Processing: Processing,
  Generators: Generators,
  Transform: Transform,
  toSpectrum: Transform.toSpectrum,
  toSignal: Transform.toSignal,
  signal: Signal,
  sine: Generators.sine,
  sines: Generators.sines,
  convolve: Processing.convolve,
  version: "0.0.3"
}

},{"./lib/generator.js":5,"./lib/processing.js":6,"./lib/signal.js":7,"./lib/transform.js":9}],3:[function(require,module,exports){
/**
  * implementation of the discrete convolution for arbitrary signals
  * doesn't require any periodicity like the fast convolution
  * 
  * preferable when the second signal (filter) is short.
  */

var _ = require("lodash");

var mod = function(m,n){
  return ((m % n) + n) % n;
}

var nonCircular = function(signal1, signal2){
  var len1 = signal1.length;
  var len2 = signal2.length;
  var half = Math.floor(len2 / 2);
  
  var convolved = _.map(_.range(len1), function(idx1){
    var cVal = _.reduce(_.range(len2), function(acc, idx2){
      var curIdx = idx1 - idx2;
      if(0 <= curIdx && curIdx < len1){
        // actual convolution
        return acc + signal1[curIdx] * signal2[idx2];
      }
      // range test. If our index is not within the signal1 
      // we simply return the current value of the convolution
      return acc;
    }, 0);
    return cVal;
  });
  
  return convolved;
}

var circular = function(signal1, signal2){
  var len1 = signal1.length;
  var len2 = signal2.length;
  var half = Math.floor(len2 / 2);
  
  var convolved = _.map(_.range(len1), function(idx1){
    var cVal = _.reduce(_.range(len1), function(acc, idx2){
      var curIdx = idx1 - idx2;
      return acc + signal1[idx2] * signal2[mod(curIdx,len2)];
    }, 0);
    return cVal;
  });
  
  return convolved;
}

var methodForType = {"non-circular" : nonCircular,"circular": circular};

module.exports = function(signal1, signal2, options){
  var convMethod = methodForType[options.type];
  return convMethod(signal1,signal2);
}

},{"lodash":18}],4:[function(require,module,exports){
/**
  * convolution using a spectrum (usually fft) method to calculate the spectrum
  * multiply them and backtransform back into a signal. Only for periodic
  * signals
  */

var transform = require("../transform.js");
var _ = require("lodash");

module.exports = function(signal1, signal2, options){
  var spec1 = transform.toSpectrum(signal1, options.spectrum);
  var spec2 = transform.toSpectrum(signal2, options.spectrum);
  // spec1.len > spec2.len every spec1.len/spec2.len frequency
  // is used. This only works for power of 2 length, but the FFT
  // requires power 2 lengths
  var idxDist = spec1.length / spec2.length;
  
  // multiply only every idxDist frequency, all other multiplications
  // are zero and thus we only need to back transform the shorter spectrum
  var convSpec = _.map(_.range(spec2.legnth), function(idx){
    return spec1[idx*idxDist] * spec2[idx];
  });
  convSpec.sampling = 1;
  return transform.toSignal(convSpec,options.spectrum).create({
    length: signal1.length,
    sampling: signal1.sampling
  });
}

},{"../transform.js":9,"lodash":18}],5:[function(require,module,exports){
/**
  * Generators for different signals (currently only sine waves)
  */

// for now we use dsp.js
var dsp = require("../dsp.js");
var _ = require("lodash");
var Signal = require("./signal.js");

/** A generator can generate a periodic function
 *  using a given function
 */
var generator = function(genFunction, cfg){
  
  var data = []
  return {
    /** creates a signal
     *  call
     *   - create() [uses default: sampling=44100, length=1000]
     *   - create({sampling: 44100, length: 1000})
     */
    create: function(options){
      options = _.defaults(options || {}, {sampling:44100, length:1000});
//      var generator = new dsp.Oscillator(cfg.__dspGenerator, 
//        cfg.frequency, cfg.amplitude, options.length, options.sampling);
//      generator.generate();
      var signalData = _.map(_.range(options.length), function(idx){
        return cfg.amplitude * genFunction(cfg.functionPeriod * cfg.frequency * idx / options.sampling);
      });
//      var signalData = generator.signal;
      var s = Signal(signalData, options);
      s.amplitude = cfg.amplitude;
      s.frequency = cfg.frequency;
      return s;
    }
  }
}

var mergeTwoSignals = function(s1,s2){
  var sigData = _(s1)
    .zip(s2)
    .map(function(p){
      return p[0].value + p[1].value;
    })
    .value();
  return Signal(sigData,{sampling:s1.sampling});
}

var mergeSignals = function(sigs){
  if(_.isArray(sigs)){
    return _.reduce(sigs,mergeTwoSignals);
  } else {
    return _.reduce(arguments,mergeTwoSignals);
  }
}

var combinator = function(generators){
  return {
    create: function(options){
      options = _.defaults(options || {}, {sampling:44100, length:1000});
      return _(generators)
        .invoke("create",options)
        .reduce(mergeTwoSignals);
    }
  }
}

var Signals = require("./utils/signals.js");

var Generators = {
  /** sine wave generator.
   *  call
   *   - sine() [uses default: frequency=440, phase = 0, amplitude: 1]
   *   - sine({frequency:440, phase:0, amplitude: 1})
   *  returns a generator for the sine wave
   */
  sine: function(options){
    options = _.defaults(options || {}, {frequency: 440, phase: 0, amplitude: 1});
    options.functionPeriod = 2*Math.PI;
    return generator(Math.sin, options);
  },
  sines: function(signals){
    var gens = _.map(signals, function(s){
      return Generators.sine(s); });
    return combinator(gens);
  },
  triangle: function(options){
    options = _.defaults(options || {}, {frequency: 440, phase: 0, amplitude: 1});
    options.functionPeriod = 1;
    return generator(Signals.triangle(1), options);
  },
  rectangle: function(options){
    options = _.defaults(options || {}, {frequency: 440, phase: 0, amplitude: 1});
    options.functionPeriod = 1;
    return generator(Signals.rectangle(1), options);
  }
}

module.exports = Generators;

},{"../dsp.js":1,"./signal.js":7,"./utils/signals.js":17,"lodash":18}],6:[function(require,module,exports){
/**
  * (Advanced) Signal processing
  */

var cauchy = require("./convolution/cauchy.js");
var specConvolution = require("./convolution/spectrum.js");
var _ = require("lodash");
var Signal = require("./signal.js");

var defaultOptions = {
  method: "cauchy",
  type: "non-circular",
  spectrum: {
    method: "dft"
    // sampling is not relevant as it gets transformed back again with the same
    // sampling rate
  }
}

module.exports = {
  /** the convolution function assumes two discrete signals
    * that have the same spacing
    */
  convolve: function(signal1, signal2, options){
    options = options || {};
    options = _.defaults(options, defaultOptions);
    
    // ensure signal2 is not longer than signal2
    if(signal1.length < signal2.length){
      var tmp = signal1;
      signal1 = signal2;
      signal2 = tmp;
    }
    
    if(options.method == "cauchy"){
      return cauchy(signal1,signal2,options);
    } else {
      return specConvolution(signal1, signal2,options);
    }
  },
  /** determines if two signals can be considered equal.
    */
  equal: function(signal1, signal2, options){
    options = options || {};
    options.epsilon = options.epsilon || 1E-08;
    var sig1 = Signal(signal1);
    var sig2 = Signal(signal2);
    if(signal1.length != signal2.length) return false;
    var diffSQ = _.reduce(_.range(signal1.length), function(d,idx){
      var diff = (sig1[idx].value-sig2[idx].value);
        return d + diff * diff;
    },0);
    return diffSQ/signal1.length < options.epsilon;
  }
}

},{"./convolution/cauchy.js":3,"./convolution/spectrum.js":4,"./signal.js":7,"lodash":18}],7:[function(require,module,exports){
/**
  * (Advanced) Signal processing
  */
  

var _ = require("lodash");

function defaultSampling(options){
  var defaultSpectrum = 22050;
  if(!options.disableWarnings || !_.contains(options.disableWarnings, "defaultSampling")) {
    //console.warn("switching to default Sampling rate of "+defaultSpectrum+".");
  }
  return defaultSpectrum;
}

var signal = function(sig, options){
  if("values" in sig && "sampling" in sig) return sig;
  options = options || {};
  var sampling = options.sampling || sig.sampling || defaultSampling(options);
  var curT = 0;
  var data = _.map(sig, function(v){
    ret = {t: curT, value: v};
    curT = curT + 1/sampling;
    return ret;
  });
  data.sampling = sampling;
  data.dt = 1/sampling;
  data.values = function(){
    return _.map(data, function(p){ return p.value; });
  };
  return data;
}

module.exports = signal;

},{"lodash":18}],8:[function(require,module,exports){

var _ = require("lodash");

var spectrumCreate = function(spectrumArray, sampling, signalLength){
  var spectrum = _.map(_.range(spectrumArray.length), function(idx){
    // frequency in hertz
    var freq = idx/signalLength * sampling;
    var freqFac = (idx == 0) ? 0.5 : 1;
    return { frequency: freq, amplitude: spectrumArray[idx]*freqFac };
  });
  spectrum.dominantFrequency = function() {
    return _.max(spectrum, function (s) {
      return s.amplitude;
    });
  };
  spectrum.amplitudes = function(){
    return _.map(spectrum, function(v){
      return v.amplitude;
    });
  }
  spectrum.sampling = sampling;
  return spectrum;
};

var timeDependentSpectrumCreate = function(spectrumArray, sampling){
  var levels = Math.floor(Math.log(spectrumArray.length)/Math.log(2));
  var curLevel = levels-1;
  var idx = Math.pow(2,curLevel);
  var spectrum = [];
  var cnt = 0;
  while(curLevel>=0){
    var freq = sampling * idx /* = 2^curLevel */ / spectrumArray.length;
    for(var i=idx; i<idx*2; i++){
      var diffIdx = i - idx;
      var stepSize = spectrumArray.length/(idx * sampling);
      spectrum[i] = {
        frequency: freq,
        amplitude: spectrumArray[i],
        timeStart: diffIdx * stepSize,
        timeEnd:(diffIdx+1)*stepSize,
        time: (diffIdx + 0.5)*stepSize
      };
    }
    curLevel = curLevel - 1;
    idx = idx / 2;
    cnt++;
  }
  spectrum[0] = { frequency: 0, amplitude: spectrumArray[0] };
  spectrum.dominantFrequency = function() {
    return _.max(spectrum, function (s) {
      return s.amplitude;
    });
  };
  spectrum.amplitudes = function(){
    return _.map(spectrum, function(v){
      return v.amplitude;
    });
  }
  spectrum.sampling = sampling;
  return spectrum;
};

var Spectrum = function(spectrumArray, options){
  if("amplitudes" in spectrumArray)
    return spectrumArray;
  if("timeDependent" in spectrumArray && spectrumArray.timeDependent)
    return timeDependentSpectrumCreate(spectrumArray, options.sampling);
  return spectrumCreate(spectrumArray, options.sampling, options.signalLength);
}

module.exports = Spectrum;

},{"lodash":18}],9:[function(require,module,exports){
/**
  * Transforms a signal into its frequencies or vice versa
  */

var Signal = require("./signal.js");

var Transform = require("./transforms/base.js");
module.exports = {
  toSignal: function(spectrum, options){
    options = options || {};
    return Transform.toSignal(spectrum, options);
  },
  toSpectrum: function(signal, options){
    options = options || {};
    options.method = options.method || "dft";
    options.signalLength = signal.length;
    return Transform.toSpectrum(Signal(signal), options);
  }
}

},{"./signal.js":7,"./transforms/base.js":10}],10:[function(require,module,exports){

var Spectrum = require("../spectrum.js");
var Signal = require("../signal.js");
var generator = require("../generator.js");

var Transform = {
  methods: {},
  /** `register` takes an object that contains a name and three functions.
    *  - `forward` and `backward` which applied to a signal/spectrum calculate the
    *   transformation into the corresponding domain
    *  - `prepare` which calls mechanims to speed up calculation like caching
    *   etc.
    */
  register: function(method){
    Transform.methods[method.name] = method;
  },
  toSpectrum: function(signal, options){
    method = Transform.methods[options.method];
    options.signalLength = signal.length;
    options.sampling = options.sampling || signal.sampling || 440;
    return Spectrum(method.forward(signal,options), options);
  },
  toSignal: function(spectrum, options){
    var spec = Spectrum(spectrum, options);
    return generator.sines(spec);
  }
}

// register DFT method
require("./dft.js").register(Transform);
require("./fft.js").register(Transform);
require("./wavelets/haar.js").register(Transform);
require("./wavelets/daubechies.js").register(Transform);

module.exports = Transform;

},{"../generator.js":5,"../signal.js":7,"../spectrum.js":8,"./dft.js":11,"./fft.js":12,"./wavelets/daubechies.js":13,"./wavelets/haar.js":15}],11:[function(require,module,exports){

var dsp = require("../../dsp.js");
var dspCache = {};
var identifier = function(options){
  return options.length + ":" + options.sampling;
};
var prepare = function(options){
  var id = identifier(options);
  if(!(id in dspCache)){
    dspCache[id] = new dsp.DFT(options.length, options.sampling);
  }
  return dspCache[id];
};

module.exports = {
  name: "DFT DSP.JS",
  register: function(Transform){
    Transform.register({
      name: "dft",
      forward: function(signal, options){
        options.length = signal.length;
        var dft = prepare(options);
        dft.forward(signal.values());
        return dft.spectrum;
      },
      backward: function(spectrum, options){
        options.length = spectrum.length;
        var dft = prepare(options);
        dft.backward(signal);
        throw new Error("?? what to return ??");
      }
    });
  }
};

},{"../../dsp.js":1}],12:[function(require,module,exports){

var dsp = require("../../dsp.js");
var isPowerOfTwo = require("../utils/poweroftwo.js");
var dspCache = {};
var identifier = function(options){
  return options.length + ":" + options.sampling;
};
var prepare = function(options){
  var id = identifier(options);
  if(!(id in dspCache)){
    if(!isPowerOfTwo(options.length)){
      throw new Error("FFT only works for a power of two signal length");
    }
    dspCache[id] = new dsp.FFT(options.length, options.sampling);
  }
  return dspCache[id];
};

module.exports = {
  name: "FFT DSP.JS",
  register: function(Transform){
    Transform.register({
      name: "fft",
      forward: function(signal, options){
        options.length = signal.length;
        var fft = prepare(options);
        fft.forward(signal.values());
        return fft.spectrum;
      },
      backward: function(spectrum, options){
        options.length = spectrum.length;
        var fft = prepare(options);
        fft.backward(signal);
        throw new Error("?? what to return ??");
      }
    });
  }
};

},{"../../dsp.js":1,"../utils/poweroftwo.js":16}],13:[function(require,module,exports){

var daubCoeffs = require("./daubechiesCoefficients.js");

var daubechiesPeriodic = function(signal, options){
  var coeffs = daubCoeffs["D"+options.taps];
  var input = signal.values();
  var copy = [];
  var res = [];
  var len = Math.floor(signal.length / 2);
  while(len > 0){
    for(var i=0; i<len; i++){
      var scaling = 0; var wavelet = 0;
      var fac = 1;
      for(var j=0; j<coeffs.length; j++){
        var idx = (2*i + j) % len;
        scaling += coeffs[j] * input[2*i+j]*0.5;
        wavelet += fac * coeffs[j] * input[2*i+j]* 0.5;
        fac = fac * -1;
      }
      copy[i] = scaling;
      res[len + i] = wavelet;
    }
    var tmp = copy;
    copy = input;
    input = tmp;
    len = Math.floor(len / 2);
  }
  res[0] = input[0];
  return res;
}

module.exports = {
  name: "DWT Daubechies",
  register: function(Transform){
    Transform.register({
      name: "daubechies",
      forward: function(signal, options){
        var trans = daubechiesPeriodic(signal, options);
        trans.timeDependent = true;
        return trans;
      },
      backward: function(spectrum, options){
        options.length = spectrum.length;
        var dft = prepare(options);
        dft.backward(signal);
        throw new Error("?? what to return ??");
      }
    });
  }
};

},{"./daubechiesCoefficients.js":14}],14:[function(require,module,exports){
//taken from http://musicdsp.org/showone.php?id=20

var Daub1 = [
7.071067811865475244008443621048490392848359376884740365883398e-01,
7.071067811865475244008443621048490392848359376884740365883398e-01];

var Daub2 = [
  4.829629131445341433748715998644486838169524195042022752011715e-01,
  8.365163037378079055752937809168732034593703883484392934953414e-01,
  2.241438680420133810259727622404003554678835181842717613871683e-01,
  -1.294095225512603811744494188120241641745344506599652569070016e-01];
  
var Daub3 = [
  3.326705529500826159985115891390056300129233992450683597084705e-01,
  8.068915093110925764944936040887134905192973949948236181650920e-01,
  4.598775021184915700951519421476167208081101774314923066433867e-01,
  -1.350110200102545886963899066993744805622198452237811919756862e-01,
  -8.544127388202666169281916918177331153619763898808662976351748e-02,
  3.522629188570953660274066471551002932775838791743161039893406e-02];
 
var Daub4 = [
  2.303778133088965008632911830440708500016152482483092977910968e-01,
  7.148465705529156470899219552739926037076084010993081758450110e-01,
  6.308807679298589078817163383006152202032229226771951174057473e-01,
  -2.798376941685985421141374718007538541198732022449175284003358e-02,
  -1.870348117190930840795706727890814195845441743745800912057770e-01,
  3.084138183556076362721936253495905017031482172003403341821219e-02,
  3.288301166688519973540751354924438866454194113754971259727278e-02,
  -1.059740178506903210488320852402722918109996490637641983484974e-02];

var Daub5 = [
  1.601023979741929144807237480204207336505441246250578327725699e-01,
  6.038292697971896705401193065250621075074221631016986987969283e-01,
  7.243085284377729277280712441022186407687562182320073725767335e-01,
  1.384281459013207315053971463390246973141057911739561022694652e-01,
  -2.422948870663820318625713794746163619914908080626185983913726e-01,
  -3.224486958463837464847975506213492831356498416379847225434268e-02,
  7.757149384004571352313048938860181980623099452012527983210146e-02,
  -6.241490212798274274190519112920192970763557165687607323417435e-03,
  -1.258075199908199946850973993177579294920459162609785020169232e-02,
  3.335725285473771277998183415817355747636524742305315099706428e-03];
  
var Daub6 = [
  1.115407433501094636213239172409234390425395919844216759082360e-01,
  4.946238903984530856772041768778555886377863828962743623531834e-01,
  7.511339080210953506789344984397316855802547833382612009730420e-01,
  3.152503517091976290859896548109263966495199235172945244404163e-01,
  -2.262646939654398200763145006609034656705401539728969940143487e-01,
  -1.297668675672619355622896058765854608452337492235814701599310e-01,
  9.750160558732304910234355253812534233983074749525514279893193e-02,
  2.752286553030572862554083950419321365738758783043454321494202e-02,
  -3.158203931748602956507908069984866905747953237314842337511464e-02,
  5.538422011614961392519183980465012206110262773864964295476524e-04,
  4.777257510945510639635975246820707050230501216581434297593254e-03,
  -1.077301085308479564852621609587200035235233609334419689818580e-03];
 
var Daub7 = [
  7.785205408500917901996352195789374837918305292795568438702937e-02,
  3.965393194819173065390003909368428563587151149333287401110499e-01,
  7.291320908462351199169430703392820517179660611901363782697715e-01,
  4.697822874051931224715911609744517386817913056787359532392529e-01,
  -1.439060039285649754050683622130460017952735705499084834401753e-01,
  -2.240361849938749826381404202332509644757830896773246552665095e-01,
  7.130921926683026475087657050112904822711327451412314659575113e-02,
  8.061260915108307191292248035938190585823820965629489058139218e-02,
  -3.802993693501441357959206160185803585446196938467869898283122e-02,
  -1.657454163066688065410767489170265479204504394820713705239272e-02,
  1.255099855609984061298988603418777957289474046048710038411818e-02,
  4.295779729213665211321291228197322228235350396942409742946366e-04,
  -1.801640704047490915268262912739550962585651469641090625323864e-03,
  3.537137999745202484462958363064254310959060059520040012524275e-04];

var Daub8 = [
  5.441584224310400995500940520299935503599554294733050397729280e-02,
  3.128715909142999706591623755057177219497319740370229185698712e-01,
  6.756307362972898068078007670471831499869115906336364227766759e-01,
  5.853546836542067127712655200450981944303266678053369055707175e-01,
  -1.582910525634930566738054787646630415774471154502826559735335e-02,
  -2.840155429615469265162031323741647324684350124871451793599204e-01,
  4.724845739132827703605900098258949861948011288770074644084096e-04,
  1.287474266204784588570292875097083843022601575556488795577000e-01,
  -1.736930100180754616961614886809598311413086529488394316977315e-02,
  -4.408825393079475150676372323896350189751839190110996472750391e-02,
  1.398102791739828164872293057263345144239559532934347169146368e-02,
  8.746094047405776716382743246475640180402147081140676742686747e-03,
  -4.870352993451574310422181557109824016634978512157003764736208e-03,
  -3.917403733769470462980803573237762675229350073890493724492694e-04,
  6.754494064505693663695475738792991218489630013558432103617077e-04,
  -1.174767841247695337306282316988909444086693950311503927620013e-04];
  
var Daub9 = [
  3.807794736387834658869765887955118448771714496278417476647192e-02,
  2.438346746125903537320415816492844155263611085609231361429088e-01,
  6.048231236901111119030768674342361708959562711896117565333713e-01,
  6.572880780513005380782126390451732140305858669245918854436034e-01,
  1.331973858250075761909549458997955536921780768433661136154346e-01,
  -2.932737832791749088064031952421987310438961628589906825725112e-01,
  -9.684078322297646051350813353769660224825458104599099679471267e-02,
  1.485407493381063801350727175060423024791258577280603060771649e-01,
  3.072568147933337921231740072037882714105805024670744781503060e-02,
  -6.763282906132997367564227482971901592578790871353739900748331e-02,
  2.509471148314519575871897499885543315176271993709633321834164e-04,
  2.236166212367909720537378270269095241855646688308853754721816e-02,
  -4.723204757751397277925707848242465405729514912627938018758526e-03,
  -4.281503682463429834496795002314531876481181811463288374860455e-03,
  1.847646883056226476619129491125677051121081359600318160732515e-03,
  2.303857635231959672052163928245421692940662052463711972260006e-04,
  -2.519631889427101369749886842878606607282181543478028214134265e-04,
  3.934732031627159948068988306589150707782477055517013507359938e-05];

var Daub10 = [
  2.667005790055555358661744877130858277192498290851289932779975e-02,
  1.881768000776914890208929736790939942702546758640393484348595e-01,
  5.272011889317255864817448279595081924981402680840223445318549e-01,
  6.884590394536035657418717825492358539771364042407339537279681e-01,
  2.811723436605774607487269984455892876243888859026150413831543e-01,
  -2.498464243273153794161018979207791000564669737132073715013121e-01,
  -1.959462743773770435042992543190981318766776476382778474396781e-01,
  1.273693403357932600826772332014009770786177480422245995563097e-01,
  9.305736460357235116035228983545273226942917998946925868063974e-02,
  -7.139414716639708714533609307605064767292611983702150917523756e-02,
  -2.945753682187581285828323760141839199388200516064948779769654e-02,
  3.321267405934100173976365318215912897978337413267096043323351e-02,
  3.606553566956169655423291417133403299517350518618994762730612e-03,
  -1.073317548333057504431811410651364448111548781143923213370333e-02,
  1.395351747052901165789318447957707567660542855688552426721117e-03,
  1.992405295185056117158742242640643211762555365514105280067936e-03,
  -6.858566949597116265613709819265714196625043336786920516211903e-04,
  -1.164668551292854509514809710258991891527461854347597362819235e-04,
  9.358867032006959133405013034222854399688456215297276443521873e-05,
  -1.326420289452124481243667531226683305749240960605829756400674e-05];
  
  var Daub11 = [
 1.869429776147108402543572939561975728967774455921958543286692e-02,
 1.440670211506245127951915849361001143023718967556239604318852e-01,
 4.498997643560453347688940373853603677806895378648933474599655e-01,
 6.856867749162005111209386316963097935940204964567703495051589e-01,
 4.119643689479074629259396485710667307430400410187845315697242e-01,
 -1.622752450274903622405827269985511540744264324212130209649667e-01,
 -2.742308468179469612021009452835266628648089521775178221905778e-01,
 6.604358819668319190061457888126302656753142168940791541113457e-02,
 1.498120124663784964066562617044193298588272420267484653796909e-01,
 -4.647995511668418727161722589023744577223260966848260747450320e-02,
 -6.643878569502520527899215536971203191819566896079739622858574e-02,
 3.133509021904607603094798408303144536358105680880031964936445e-02,
 2.084090436018106302294811255656491015157761832734715691126692e-02,
 -1.536482090620159942619811609958822744014326495773000120205848e-02,
 -3.340858873014445606090808617982406101930658359499190845656731e-03,
 4.928417656059041123170739741708273690285547729915802418397458e-03,
 -3.085928588151431651754590726278953307180216605078488581921562e-04,
 -8.930232506662646133900824622648653989879519878620728793133358e-04,
 2.491525235528234988712216872666801088221199302855425381971392e-04,
 5.443907469936847167357856879576832191936678525600793978043688e-05,
 -3.463498418698499554128085159974043214506488048233458035943601e-05,
 4.494274277236510095415648282310130916410497987383753460571741e-06];
 
 var Daub12 = [
1.311225795722951750674609088893328065665510641931325007748280e-02,
1.095662728211851546057045050248905426075680503066774046383657e-01,
3.773551352142126570928212604879206149010941706057526334705839e-01,
6.571987225793070893027611286641169834250203289988412141394281e-01,
5.158864784278156087560326480543032700677693087036090056127647e-01,
-4.476388565377462666762747311540166529284543631505924139071704e-02,
-3.161784537527855368648029353478031098508839032547364389574203e-01,
-2.377925725606972768399754609133225784553366558331741152482612e-02,
1.824786059275796798540436116189241710294771448096302698329011e-01,
5.359569674352150328276276729768332288862665184192705821636342e-03,
-9.643212009650708202650320534322484127430880143045220514346402e-02,
1.084913025582218438089010237748152188661630567603334659322512e-02,
4.154627749508444073927094681906574864513532221388374861287078e-02,
-1.221864906974828071998798266471567712982466093116558175344811e-02,
-1.284082519830068329466034471894728496206109832314097633275225e-02,
6.711499008795509177767027068215672450648112185856456740379455e-03,
2.248607240995237599950865211267234018343199786146177099262010e-03,
-2.179503618627760471598903379584171187840075291860571264980942e-03,
6.545128212509595566500430399327110729111770568897356630714552e-06,
3.886530628209314435897288837795981791917488573420177523436096e-04,
-8.850410920820432420821645961553726598738322151471932808015443e-05,
-2.424154575703078402978915320531719580423778362664282239377532e-05,
1.277695221937976658714046362616620887375960941439428756055353e-05,
-1.529071758068510902712239164522901223197615439660340672602696e-06];

var Daub13 = [
  9.202133538962367972970163475644184667534171916416562386009703e-03,
  8.286124387290277964432027131230466405208113332890135072514277e-02,
  3.119963221604380633960784112214049693946683528967180317160390e-01,
  6.110558511587876528211995136744180562073612676018239438526582e-01,
  5.888895704312189080710395347395333927665986382812836042235573e-01,
  8.698572617964723731023739838087494399231884076619701250882016e-02,
  -3.149729077113886329981698255932282582876888450678789025950306e-01,
  -1.245767307508152589413808336021260180792739295173634719572069e-01,
  1.794760794293398432348450072339369013581966256244133393042881e-01,
  7.294893365677716380902830610477661983325929026879873553627963e-02,
  -1.058076181879343264509667304196464849478860754801236658232360e-01,
  -2.648840647534369463963912248034785726419604844297697016264224e-02,
  5.613947710028342886214501998387331119988378792543100244737056e-02,
  2.379972254059078811465170958554208358094394612051934868475139e-03,
  -2.383142071032364903206403067757739134252922717636226274077298e-02,
  3.923941448797416243316370220815526558824746623451404043918407e-03,
  7.255589401617566194518393300502698898973529679646683695269828e-03,
  -2.761911234656862178014576266098445995350093330501818024966316e-03,
  -1.315673911892298936613835370593643376060412592653652307238124e-03,
  9.323261308672633862226517802548514100918088299801952307991569e-04,
  4.925152512628946192140957387866596210103778299388823500840094e-05,
  -1.651289885565054894616687709238000755898548214659776703347801e-04,
  3.067853757932549346649483228575476236600428217237900563128230e-05,
  1.044193057140813708170714991080596951670706436217328169641474e-05,
  -4.700416479360868325650195165061771321650383582970958556568059e-06,
  5.220035098454864691736424354843176976747052155243557001531901e-07];
  
  var Daub14 = [
 6.461153460087947818166397448622814272327159419201199218101404e-03,
 6.236475884939889832798566758434877428305333693407667164602518e-02,
 2.548502677926213536659077886778286686187042416367137443780084e-01,
 5.543056179408938359926831449851154844078269830951634609683997e-01,
 6.311878491048567795576617135358172348623952456570017289788809e-01,
 2.186706877589065214917475918217517051765774321270432059030273e-01,
 -2.716885522787480414142192476181171094604882465683330814311896e-01,
 -2.180335299932760447555558812702311911975240669470604752747127e-01,
 1.383952138648065910739939690021573713989900463229686119059119e-01,
 1.399890165844607012492943162271163440328221555614326181333683e-01,
 -8.674841156816968904560822066727795382979149539517503657492964e-02,
 -7.154895550404613073584145115173807990958069673129538099990913e-02,
 5.523712625921604411618834060533403397913833632511672157671107e-02,
 2.698140830791291697399031403215193343375766595807274233284349e-02,
 -3.018535154039063518714822623489137573781575406658652624883756e-02,
 -5.615049530356959133218371367691498637457297203925810387698680e-03,
 1.278949326633340896157330705784079299374903861572058313481534e-02,
 -7.462189892683849371817160739181780971958187988813302900435487e-04,
 -3.849638868022187445786349316095551774096818508285700493058915e-03,
 1.061691085606761843032566749388411173033941582147830863893939e-03,
 7.080211542355278586442977697617128983471863464181595371670094e-04,
 -3.868319473129544821076663398057314427328902107842165379901468e-04,
 -4.177724577037259735267979539839258928389726590132730131054323e-05,
 6.875504252697509603873437021628031601890370687651875279882727e-05,
 -1.033720918457077394661407342594814586269272509490744850691443e-05,
 -4.389704901781394115254042561367169829323085360800825718151049e-06,
 1.724994675367812769885712692741798523587894709867356576910717e-06,
 -1.787139968311359076334192938470839343882990309976959446994022e-07];
 
 var Daub15 = [
4.538537361578898881459394910211696346663671243788786997916513e-03,
4.674339489276627189170969334843575776579151700214943513113197e-02,
2.060238639869957315398915009476307219306138505641930902702047e-01,
4.926317717081396236067757074029946372617221565130932402160160e-01,
6.458131403574243581764209120106917996432608287494046181071489e-01,
3.390025354547315276912641143835773918756769491793554669336690e-01,
-1.932041396091454287063990534321471746304090039142863827937754e-01,
-2.888825965669656462484125009822332981311435630435342594971292e-01,
6.528295284877281692283107919869574882039174285596144125965101e-02,
1.901467140071229823484893116586020517959501258174336696878156e-01,
-3.966617655579094448384366751896200668381742820683736805449745e-02,
-1.111209360372316933656710324674058608858623762165914120505657e-01,
3.387714392350768620854817844433523770864744687411265369463195e-02,
5.478055058450761268913790312581879108609415997422768564244845e-02,
-2.576700732843996258594525754269826392203641634825340138396836e-02,
-2.081005016969308167788483424677000162054657951364899040996166e-02,
1.508391802783590236329274460170322736244892823305627716233968e-02,
5.101000360407543169708860185565314724801066527344222055526631e-03,
-6.487734560315744995181683149218690816955845639388826407928967e-03,
-2.417564907616242811667225326300179605229946995814535223329411e-04,
1.943323980382211541764912332541087441011424865579531401452302e-03,
-3.734823541376169920098094213645414611387630968030256625740226e-04,
-3.595652443624688121649620075909808858194202454084090305627480e-04,
1.558964899205997479471658241227108816255567059625495915228603e-04,
2.579269915531893680925862417616855912944042368767340709160119e-05,
-2.813329626604781364755324777078478665791443876293788904267255e-05,
3.362987181737579803124845210420177472134846655864078187186304e-06,
1.811270407940577083768510912285841160577085925337507850590290e-06,
-6.316882325881664421201597299517657654166137915121195510416641e-07,
6.133359913305752029056299460289788601989190450885396512173845e-08];

var Daub16 = [
  3.189220925347738029769547564645958687067086750131428767875878e-03,
  3.490771432367334641030147224023020009218241430503984146140054e-02,
  1.650642834888531178991252730561134811584835002342723240213592e-01,
  4.303127228460038137403925424357684620633970478036986773924646e-01,
  6.373563320837888986319852412996030536498595940814198125967751e-01,
  4.402902568863569000390869163571679288527803035135272578789884e-01,
  -8.975108940248964285718718077442597430659247445582660149624718e-02,
  -3.270633105279177046462905675689119641757228918228812428141723e-01,
  -2.791820813302827668264519595026873204339971219174736041535479e-02,
  2.111906939471042887209680163268837900928491426167679439251042e-01,
  2.734026375271604136485245757201617965429027819507130220231500e-02,
  -1.323883055638103904500474147756493375092287817706027978798549e-01,
  -6.239722752474871765674503394120025865444656311678760990761458e-03,
  7.592423604427631582148498743941422461530405946100943351940313e-02,
  -7.588974368857737638494890864636995796586975144990925400097160e-03,
  -3.688839769173014233352666320894554314718748429706730831064068e-02,
  1.029765964095596941165000580076616900528856265803662208854147e-02,
  1.399376885982873102950451873670329726409840291727868988490100e-02,
  -6.990014563413916670284249536517288338057856199646469078115759e-03,
  -3.644279621498389932169000540933629387055333973353108668841215e-03,
  3.128023381206268831661202559854678767821471906193608117450360e-03,
  4.078969808497128362417470323406095782431952972310546715071397e-04,
  -9.410217493595675889266453953635875407754747216734480509250273e-04,
  1.142415200387223926440228099555662945839684344936472652877091e-04,
  1.747872452253381803801758637660746874986024728615399897971953e-04,
  -6.103596621410935835162369150522212811957259981965919143961722e-05,
  -1.394566898820889345199078311998401982325273569198675335408707e-05,
  1.133660866127625858758848762886536997519471068203753661757843e-05,
  -1.043571342311606501525454737262615404887478930635676471546032e-06,
  -7.363656785451205512099695719725563646585445545841663327433569e-07,
  2.308784086857545866405412732942006121306306735866655525372544e-07,
  -2.109339630100743097000572623603489906836297584591605307745349e-08];
  
  var Daub17 = [
 2.241807001037312853535962677074436914062191880560370733250531e-03,
 2.598539370360604338914864591720788315473944524878241294399948e-02,
 1.312149033078244065775506231859069960144293609259978530067004e-01,
 3.703507241526411504492548190721886449477078876896803823650425e-01,
 6.109966156846228181886678867679372082737093893358726291371783e-01,
 5.183157640569378393254538528085968046216817197718416402439904e-01,
 2.731497040329363500431250719147586480350469818964563003672942e-02,
 -3.283207483639617360909665340725061767581597698151558024679130e-01,
 -1.265997522158827028744679110933825505053966260104086162103728e-01,
 1.973105895650109927854047044781930142551422414135646917122284e-01,
 1.011354891774702721509699856433434802196622545499664876109437e-01,
 -1.268156917782863110948571128662331680384792185915017065732137e-01,
 -5.709141963167692728911239478651382324161160869845347053990144e-02,
 8.110598665416088507965885748555429201024364190954499194020678e-02,
 2.231233617810379595339136059534813756232242114093689244020869e-02,
 -4.692243838926973733300897059211400507138768125498030602878439e-02,
 -3.270955535819293781655360222177494452069525958061609392809275e-03,
 2.273367658394627031845616244788448969906713741338339498024864e-02,
 -3.042989981354637068592482637907206078633395457225096588287881e-03,
 -8.602921520322854831713706413243659917926736284271730611920986e-03,
 2.967996691526094872806485060008038269959463846548378995044195e-03,
 2.301205242153545624302059869038423604241976680189447476064764e-03,
 -1.436845304802976126222890402980384903503674530729935809561434e-03,
 -3.281325194098379713954444017520115075812402442728749700195651e-04,
 4.394654277686436778385677527317841632289249319738892179465910e-04,
 -2.561010956654845882729891210949920221664082061531909655178413e-05,
 -8.204803202453391839095482576282189866136273049636764338689593e-05,
 2.318681379874595084482068205706277572106695174091895338530734e-05,
 6.990600985076751273204549700855378627762758585902057964027481e-06,
 -4.505942477222988194102268206378312129713572600716499944918416e-06,
 3.016549609994557415605207594879939763476168705217646897702706e-07,
 2.957700933316856754979905258816151367870345628924317307354639e-07,
 -8.423948446002680178787071296922877068410310942222799622593133e-08,
 7.267492968561608110879767441409035034158581719789791088892046e-09];
 
 var Daub18 = [
1.576310218440760431540744929939777747670753710991660363684429e-03,
1.928853172414637705921391715829052419954667025288497572236714e-02,
1.035884658224235962241910491937253596470696555220241672976224e-01,
3.146789413370316990571998255652579931786706190489374509491307e-01,
5.718268077666072234818589370900623419393673743130930561295324e-01,
5.718016548886513352891119994065965025668047882818525060759395e-01,
1.472231119699281415750977271081072312557864107355701387801677e-01,
-2.936540407365587442479030994981150723935710729035053239661752e-01,
-2.164809340051429711237678625668271471437937235669492408388692e-01,
1.495339755653777893509301738913667208804816691893765610261943e-01,
1.670813127632574045149318139950134745324205646353988083152250e-01,
-9.233188415084628060429372558659459731431848000144569612074508e-02,
-1.067522466598284855932200581614984861385266404624112083917702e-01,
6.488721621190544281947577955141911463129382116634147846137149e-02,
5.705124773853688412090768846499622260596226120431038524600676e-02,
-4.452614190298232471556143559744653492971477891439833592755034e-02,
-2.373321039586000103275209582665216110197519330713490233071565e-02,
2.667070592647059029987908631672020343207895999936072813363471e-02,
6.262167954305707485236093144497882501990325204745013190268052e-03,
-1.305148094661200177277636447600807169755191054507571666606133e-02,
1.186300338581174657301741592161819084544899417452317405185615e-04,
4.943343605466738130665529516802974834299638313366477765295203e-03,
-1.118732666992497072800658855238650182318060482584970145512687e-03,
-1.340596298336106629517567228251583609823044524685986640323942e-03,
6.284656829651457125619449885420838217551022796301582874349652e-04,
2.135815619103406884039052814341926025873200325996466522543440e-04,
-1.986485523117479485798245416362489554927797880264017876139605e-04,
-1.535917123534724675069770335876717193700472427021513236587288e-07,
3.741237880740038181092208138035393952304292615793985030731363e-05,
-8.520602537446695203919254911655523022437596956226376512305917e-06,
-3.332634478885821888782452033341036827311505907796498439829337e-06,
1.768712983627615455876328730755375176412501359114058815453100e-06,
-7.691632689885176146000152878539598405817397588156525116769908e-08,
-1.176098767028231698450982356561292561347579777695396953528141e-07,
3.068835863045174800935478294933975372450179787894574492930570e-08,
-2.507934454948598267195173183147126731806317144868275819941403e-09];

var Daub19 = [
  1.108669763181710571099154195209715164245299677773435932135455e-03,
  1.428109845076439737439889152950199234745663442163665957870715e-02,
  8.127811326545955065296306784901624839844979971028620366497726e-02,
  2.643884317408967846748100380289426873862377807211920718417385e-01,
  5.244363774646549153360575975484064626044633641048072116393160e-01,
  6.017045491275378948867077135921802620536565639585963293313931e-01,
  2.608949526510388292872456675310528324172673101301907739925213e-01,
  -2.280913942154826463746325776054637207093787237086425909534822e-01,
  -2.858386317558262418545975695028984237217356095588335149922119e-01,
  7.465226970810326636763433111878819005865866149731909656365399e-02,
  2.123497433062784888090608567059824197077074200878839448416908e-01,
  -3.351854190230287868169388418785731506977845075238966819814032e-02,
  -1.427856950387365749779602731626112812998497706152428508627562e-01,
  2.758435062562866875014743520162198655374474596963423080762818e-02,
  8.690675555581223248847645428808443034785208002468192759640352e-02,
  -2.650123625012304089901835843676387361075068017686747808171345e-02,
  -4.567422627723090805645444214295796017938935732115630050880109e-02,
  2.162376740958504713032984257172372354318097067858752542571020e-02,
  1.937554988917612764637094354457999814496885095875825546406963e-02,
  -1.398838867853514163250401235248662521916813867453095836808366e-02,
  -5.866922281012174726584493436054373773814608340808758177372765e-03,
  7.040747367105243153014511207400620109401689897665383078229398e-03,
  7.689543592575483559749139148673955163477947086039406129546422e-04,
  -2.687551800701582003957363855070398636534038920982478290170267e-03,
  3.418086534585957765651657290463808135214214848819517257794031e-04,
  7.358025205054352070260481905397281875183175792779904858189494e-04,
  -2.606761356786280057318315130897522790383939362073563408613547e-04,
  -1.246007917341587753449784408901653990317341413341980904757592e-04,
  8.711270467219922965416862388191128268412933893282083517729443e-05,
  5.105950487073886053049222809934231573687367992106282669389264e-06,
  -1.664017629715494454620677719899198630333675608812018108739144e-05,
  3.010964316296526339695334454725943632645798938162427168851382e-06,
  1.531931476691193063931832381086636031203123032723477463624141e-06,
  -6.862755657769142701883554613486732854452740752771392411758418e-07,
  1.447088298797844542078219863291615420551673574071367834316167e-08,
  4.636937775782604223430857728210948898871748291085962296649320e-08,
  -1.116402067035825816390504769142472586464975799284473682246076e-08,
  8.666848838997619350323013540782124627289742190273059319122840e-10];
  
  var Daub20 = [
 7.799536136668463215861994818889370970510722039232863880031127e-04,
 1.054939462495039832454480973015641498231961468733236691299796e-02,
 6.342378045908151497587346582668785136406523315729666353643372e-02,
 2.199421135513970450080335972537209392121306761010882209298252e-01,
 4.726961853109016963710241465101446230757804141171727845834637e-01,
 6.104932389385938201631515660084201906858628924695448898824748e-01,
 3.615022987393310629195602665268631744967084723079677894136358e-01,
 -1.392120880114838725806970545155530518264944915437808314813582e-01,
 -3.267868004340349674031122837905370666716645587480021744425550e-01,
 -1.672708830907700757517174997304297054003744303620479394006890e-02,
 2.282910508199163229728429126648223086437547237250290835639880e-01,
 3.985024645777120219790581076522174181104027576954427684456660e-02,
 -1.554587507072679559315307870562464374359996091752285157077477e-01,
 -2.471682733861358401587992299169922262915151413349313513685587e-02,
 1.022917191744425578861013681016866083888381385233081516583444e-01,
 5.632246857307435506953246988215209861566800664402785938591145e-03,
 -6.172289962468045973318658334083283558209278762007041823250642e-02,
 5.874681811811826491300679742081997167209743446956901841959711e-03,
 3.229429953076958175885440860617219117564558605035979601073235e-02,
 -8.789324923901561348753650366700695916503030939283830968151332e-03,
 -1.381052613715192007819606423860356590496904285724730356602106e-02,
 6.721627302259456835336850521405425560520025237915708362002910e-03,
 4.420542387045790963058229526673514088808999478115581153468068e-03,
 -3.581494259609622777556169638358238375765194248623891034940330e-03,
 -8.315621728225569192482585199373230956924484221135739973390038e-04,
 1.392559619323136323905254999347967283760544147397530531142397e-03,
 -5.349759843997695051759716377213680036185796059087353172073952e-05,
 -3.851047486992176060650288501475716463266233035937022303649838e-04,
 1.015328897367029050797488785306056522529979267572003990901472e-04,
 6.774280828377729558011184406727978221295796652200819839464354e-05,
 -3.710586183394712864227221271216408416958225264980612822617745e-05,
 -4.376143862183996810373095822528607606900620592585762190542483e-06,
 7.241248287673620102843105877497181565468725757387007139555885e-06,
 -1.011994010018886150340475413756849103197395069431085005709201e-06,
 -6.847079597000556894163334787575159759109091330092963990364192e-07,
 2.633924226270001084129057791994367121555769686616747162262697e-07,
 2.014322023550512694324757845944026047904414136633776958392681e-10,
 -1.814843248299695973210605258227024081458531110762083371310917e-08,
 4.056127055551832766099146230616888024627380574113178257963252e-09,
 -2.998836489619319566407767078372705385732460052685621923178375e-10];
 
 var Daub21 = [
5.488225098526837086776336675992521426750673054588245523834775e-04,
7.776639052354783754338787398088799862510779059555623704879234e-03,
4.924777153817727491399853378340056968104483161598320693657954e-02,
1.813596254403815156260378722764624190931951510708050516519181e-01,
4.196879449393627730946850609089266339973601543036294871772653e-01,
6.015060949350038975629880664020955953066542593896126705346122e-01,
4.445904519276003403643290994523601016151342743089878478478962e-01,
-3.572291961725529045922914178005307189036762547143966578066838e-02,
-3.356640895305295094832978867114363069987575282256098351499731e-01,
-1.123970715684509813515004981340306901641824212464197973490295e-01,
2.115645276808723923846781645238468659430862736248896128529373e-01,
1.152332984396871041993434411681730428103160016594558944687967e-01,
-1.399404249325472249247758764839776903226503657502071670245304e-01,
-8.177594298086382887387303634193790542522570670234556157566786e-02,
9.660039032372422070232189700372539681627783322249829842275517e-02,
4.572340574922879239251202944731235421034828710753381191345186e-02,
-6.497750489373232063332311106008616685748929419452249544690967e-02,
-1.865385920211851534093244412008141266131208093007217139232170e-02,
3.972683542785044175197464400756126818299918992482587866999707e-02,
3.357756390338110842532604766376200760791669954106679933144723e-03,
-2.089205367797907948785235479746212371728219866525211135343707e-02,
2.403470920805434762380632169785689545910525667396313550679652e-03,
8.988824381971911875349463398395464114417817949738911101372312e-03,
-2.891334348588901247375268718015882610844675931117463495551958e-03,
-2.958374038932831280750770228215510959830170264176955719827510e-03,
1.716607040630624138494506282569230126333308533535502799235333e-03,
6.394185005120302146432543767052865436099994387647359452249347e-04,
-6.906711170821016507268939228893784790518270744313525548714065e-04,
-3.196406277680437193708834220804640347636984901270948088339102e-05,
1.936646504165080615323696689856004910579777568504218782029027e-04,
-3.635520250086338309442855006186370752206331429871136596927137e-05,
-3.499665984987447953974079490046597240276268044409625722689849e-05,
1.535482509276049283124233498646050472096482329299719141107128e-05,
2.790330539814487046106169582691767916283793946025922387556917e-06,
-3.090017164545699197158555936852697325985864588418167982685400e-06,
3.166095442367030556603889009833954440058545355777781782000278e-07,
2.992136630464852794401294607536813682771292352506328096125857e-07,
-1.000400879030597332045460600516621971679363965166249211063755e-07,
-2.254014974673330131563184851456825991617915549643308754828159e-09,
7.058033541231121859020947976903904685464512825731230495144226e-09,
-1.471954197650365265189549600816698778213247061389470277337173e-09,
1.038805571023706553035373138760372703492942617518816122570050e-10];

var Daub22 = [
  3.862632314910982158524358900615460368877852009576899680767316e-04,
  5.721854631334539120809783403484493333555361591386208129183833e-03,
  3.806993723641108494769873046391825574447727068953448390456335e-02,
  1.483675408901114285014404448710249837385836373168215616427030e-01,
  3.677286834460374788614690818452372827430535649696462720334897e-01,
  5.784327310095244271421181831735444106385099957908657145590104e-01,
  5.079010906221639018391523325390716836568713192498711562711282e-01,
  7.372450118363015165570139016530653113725172412104955350368114e-02,
  -3.127265804282961918033226222621788537078452535993545440716988e-01,
  -2.005684061048870939324361244042200174132905844868237447130382e-01,
  1.640931881067664818606223226286885712554385317412228836705888e-01,
  1.799731879928913037252154295313083168387840791424988422757762e-01,
  -9.711079840911470969274209179691733251456735137994201552926799e-02,
  -1.317681376866834107513648518146838345477875022352088357523838e-01,
  6.807631439273221556739202147004580559367442550641388181886023e-02,
  8.455737636682607503362813659356786494357635805197410905877078e-02,
  -5.136425429744413245727949984018884707909441768477091944584584e-02,
  -4.653081182750671347875833607846979997825771277976548080904423e-02,
  3.697084662069802057615318892988581825637896696876361343354380e-02,
  2.058670762756536044060249710676656807281671451609632981487139e-02,
  -2.348000134449318868560142854519364987363882333754753819791381e-02,
  -6.213782849364658499069336123807608293122900450508440420104462e-03,
  1.256472521834337406887017835495604463815382993214296088172221e-02,
  3.001373985076435951229129255588255746904937042979316054485183e-04,
  -5.455691986156717076595353163071679107868762395367234726592273e-03,
  1.044260739186025323350755659184734060807432172611689413745029e-03,
  1.827010495657279080112597436850157110235336772062961041154607e-03,
  -7.706909881231196232880372722955519781655769913634565757339739e-04,
  -4.237873998391800799531947768003976978197438302533528661825758e-04,
  3.286094142136787341983758471405935405823323072829619248523697e-04,
  4.345899904532003379046992625575076092823809665933575578710696e-05,
  -9.405223634815760421845190098352673647881298980040512091599943e-05,
  1.137434966212593172736144274866639210339820203135670505287250e-05,
  1.737375695756189356163565074505405906859746605867772002320509e-05,
  -6.166729316467578372152251668422979152169587307212708981768966e-06,
  -1.565179131995160159307426993578204733378112742579926503832095e-06,
  1.295182057318877573889711232345068147800395721925682566394936e-06,
  -8.779879873361286276888117046153049053917243760475816789226764e-08,
  -1.283336228751754417819693932114064887075096030264748079976736e-07,
  3.761228749337362366156711648187743399164239397803629022612862e-08,
  1.680171404922988885554331183691280245962290247654438114807112e-09,
  -2.729623146632976083449327361739104754443221903317745768938846e-09,
  5.335938821667489905169783227036804533253011117886586305435615e-10,
  -3.602113484339554703794807810939301847299106970237814334104274e-11];
  
  var Daub23 = [
 2.719041941282888414192673609703302357098336003920923958924757e-04,
 4.202748893183833538390034372523511472345215563611003407984701e-03,
 2.931000365788411514736204018929480427874317460676079959515131e-02,
 1.205155317839719336306053895611899089004274336891709067958035e-01,
 3.184508138528652363416527748460472152790575031409830417259640e-01,
 5.449311478735204282674240672421984387504149924834544495466793e-01,
 5.510185172419193913452724227212507720514144116478727269717859e-01,
 1.813926253638400136259098302138614937264260737638175539416540e-01,
 -2.613921480306441118856795735210118413900307577511142987337375e-01,
 -2.714020986078430556604069575184718123763697177381058877113471e-01,
 9.212540708241805260646030910734894258577648089100630012130261e-02,
 2.235736582420402317149513960822561717689875252792817094811874e-01,
 -3.303744709428937875006612792463031409461636228731285046551636e-02,
 -1.640113215318759250156057837165276039181451149292112929401186e-01,
 2.028307457564929974897286607551313323418860610791382310375731e-02,
 1.122970436181072886950734465075645977754665593869789965874572e-01,
 -2.112621235622724100704783293549467048999443844657058425212982e-02,
 -7.020739157490110946204219011957565343899895499962369353294028e-02,
 2.176585683449997560776882472168730165799461445156766923497545e-02,
 3.849533252256919901057154320407596073180564628069920893870768e-02,
 -1.852351365015615979794689960740674782817814176166333519597796e-02,
 -1.753710100303584537915846117408613551147985251726558719415169e-02,
 1.275194393152828646243157404474947115052750581861997731041018e-02,
 6.031840650024162816289878206037841640814102314209075233751820e-03,
 -7.075319273706152814194039481466556204493276773483821748740018e-03,
 -1.134865473356251691289337120013286756337393784110786907825400e-03,
 3.122876449818144997419144765125750522437659393621577492535411e-03,
 -2.465014005163512031940473100375377210862560761576109755841161e-04,
 -1.061231228886651321139357625683805642193648671030425010215075e-03,
 3.194204927099011503676530359692366990929679170022583007683112e-04,
 2.567624520078737205563856675376636092314813400664190770435450e-04,
 -1.500218503490340967673163290447832236259277810659068637402668e-04,
 -3.378894834120903434270962452674534330903724108906662510305045e-05,
 4.426071203109246077621875303440935335701832843654692827539837e-05,
 -2.635207889249186237209225933170897825432335273771458456888097e-06,
 -8.347875567854625544366043748844183086765894974439245409223337e-06,
 2.397569546840240057403739507525641239509517148980849889986407e-06,
 8.147574834779447778085443041422881439860288287528356019216814e-07,
 -5.339005405209421154584783682848780965053642859373536945701365e-07,
 1.853091785633965019353699857864654181728710556702529908304185e-08,
 5.417549179539278736503176166323685597634496102979977037271945e-08,
 -1.399935495437998845130909687361847103274208993447892120341999e-08,
 -9.472885901812050535221582074673490573092096712822067564903012e-10,
 1.050446453696543404071105111096438573423068913105255997908040e-09,
 -1.932405111313417542192651899622541612314066389643607507706887e-10,
 1.250203302351040941433216718217504240541423430995137507404787e-11];
 
 var Daub24 = [
1.914358009475513695026138336474115599435172088053846745168462e-04,
3.082081714905494436206199424544404720984720556128685270556458e-03,
2.248233994971641072358415157184825628226776692231940577581580e-02,
9.726223583362519663806545734008355914527504417674578571164300e-02,
2.729089160677263268706137134412557268751671263458895098625356e-01,
5.043710408399249919771876890402814109246866444441814540282099e-01,
5.749392210955419968460807901923407033144945935105622912839838e-01,
2.809855532337118833442626085115402941842959475929278883281409e-01,
-1.872714068851562376981887159775791469060265778441667840307934e-01,
-3.179430789993627375453948489797707550898087789160025182664299e-01,
4.776613684344728187950198323031360866349104994035553200788631e-03,
2.392373887803108551973268291945824822214858134512317715815616e-01,
4.252872964148383258147364472170645232684343235486951540533893e-02,
-1.711753513703468896897638515080572393949165942335556397917666e-01,
-3.877717357792001620177594726199572688446488033750771020190283e-02,
1.210163034692242362312637311149062286659377039046006801523826e-01,
2.098011370914481534980883827326017063121637262728447783605518e-02,
-8.216165420800166702291466006164189460916816748629968198028898e-02,
-4.578436241819221637997516339765068825260159169893967894877272e-03,
5.130162003998087915555334881398688958843078494595140394873884e-02,
-4.944709428125628299815920032649550811877887219282751174798211e-03,
-2.821310709490189098113895361900699228886900995412759197674058e-02,
7.661721881646585897329899904308764405384658404613669817843430e-03,
1.304997087108573583052494067883717533043101857128653233783396e-02,
-6.291435370018187780721843581169343900864298634085743861509767e-03,
-4.746568786323113800477796959513558401732252800905982385017245e-03,
3.736046178282523345179052160810332868725126356493155728625572e-03,
1.153764936839481504858282495202271984454410046682805375157566e-03,
-1.696456818974824394274534636412116243080312601322325642741589e-03,
-4.416184856141520063365958900079406737636243682138363561877750e-05,
5.861270593183109933716735450272894035425792347806515678695765e-04,
-1.181233237969554740613021227756568966806892308457221016257961e-04,
-1.460079817762616838924301818082729036314539476811023255670666e-04,
6.559388639305634085303738560455061974369354538271316071502698e-05,
2.183241460466558363365044032984257709791187640963509380549307e-05,
-2.022888292612697682860859987200455702614855595412267510558659e-05,
1.341157750809114719319937553186023660581084151828593222893663e-08,
3.901100338597702610409014129024223853127911530009766793352492e-06,
-8.980253143938407724149926669980791166378388013293887718404796e-07,
-4.032507756879971624098983247358983425236092110387724315244646e-07,
2.166339653278574639176393978510246335478946697396400359281412e-07,
-5.057645419792500308492508924343248979317507866520688417567606e-10,
-2.255740388176086107368821674947175804005323153443170526520277e-08,
5.157776789671999638950774266313208715015419699643333784626363e-09,
4.748375824256231118094453549799175824526559994333227456737433e-10,
-4.024658644584379774251499574468195118601698713554294941756559e-10,
6.991801157638230974132696433509625934021677793453732225542951e-11,
-4.342782503803710247259037552886749457951053124203814185811297e-12];

var Daub25 = [
  1.348029793470188994578489247159356055370460656508881471268611e-04,
  2.256959591854779520121391049628056149270016860666661928130747e-03,
  1.718674125404015533817186914954848902241194002444696221013131e-02,
  7.803586287213267559750659320481403668422052199257139168386084e-02,
  2.316935078860218199900621518057089104946216881512075361624214e-01,
  4.596834151460945937896973864539659944010260858049947396093277e-01,
  5.816368967460577833534892038757085635755639698734580573323031e-01,
  3.678850748029466984371319740855532278670733841012809062966976e-01,
  -9.717464096463814276130048169040892607068486428294030952842447e-02,
  -3.364730796417461309562110148848845218930261030262170601615289e-01,
  -8.758761458765466140226687673880006154266689569025041229545538e-02,
  2.245378197451017129525176510409543155930843160711989062118482e-01,
  1.181552867199598604563067876819931882639429216001523151773895e-01,
  -1.505602137505796309518094206831433270850173484773520730186277e-01,
  -9.850861528996022153725952822686729410420350758543226219234795e-02,
  1.066338050184779528831274540522414711301747903916268438037723e-01,
  6.675216449401860666895983072443984697329752470942906490126865e-02,
  -7.708411105657419356208567671699032054872853174701595359329826e-02,
  -3.717396286112250887598137324046870459877639250821705817221557e-02,
  5.361790939877949960629041419546536897037332284703545849594129e-02,
  1.554260592910229163981295854603203625062268043511894295387375e-02,
  -3.404232046065334099320628584033729153497903561399447916116575e-02,
  -3.079836794847036661636693963570288706232460663070983852354326e-03,
  1.892280447662762841086581178691039363674755753459524525886478e-02,
  -1.989425782202736494289461896386235348901617760816745484282494e-03,
  -8.860702618046368399013064252456556969199612331833605310278698e-03,
  2.726936258738495739871469244610042793734119359765762028996059e-03,
  3.322707773973191780118197357194829286271392998979276105842863e-03,
  -1.842484290203331280837780430014195744813667655929909114672154e-03,
  -8.999774237462950491085382524008429604309720852269895692000702e-04,
  8.772581936748274843488806190175921376284150686011179612908221e-04,
  1.153212440466300456460181455345639872216326644527860903202733e-04,
  -3.098800990984697989530544245356271119416614147098459162436317e-04,
  3.543714523276059005284289830559259809540337561365927850248007e-05,
  7.904640003965528255137496303166001735463107762646364003487560e-05,
  -2.733048119960041746353244004225286857636045649642652816856524e-05,
  -1.277195293199783804144903848434605690990373526086311486716394e-05,
  8.990661393062588905369930197413951232059323587543226269327396e-06,
  5.232827708153076417963912065899772684403904504491727061662335e-07,
  -1.779201332653634562565948556039009149458987774189389221295909e-06,
  3.212037518862519094895005816661093988294166712919881121802831e-07,
  1.922806790142371601278104244711267420759978799176017569693322e-07,
  -8.656941732278507163388031517930974947984281611717187862530250e-08,
  -2.611598556111770864259843089151782206922842627174274274741722e-09,
  9.279224480081372372250073354726511359667401736947170444723772e-09,
  -1.880415755062155537197782595740975189878162661203102565611681e-09,
  -2.228474910228168899314793352064795957306403503495743572518755e-10,
  1.535901570162657197021927739530721955859277615795931442682785e-10,
  -2.527625163465644811048864286169758128142169484216932624854015e-11,
  1.509692082823910867903367712096001664979004526477422347957324e-12];
  
  var Daub26 = [
 9.493795750710592117802731381148054398461637804818126397577999e-05,
 1.650520233532988247022384885622071050555268137055829216839523e-03,
 1.309755429255850082057770240106799154079932963479202407364818e-02,
 6.227474402514960484193581705107415937690538641013309745983962e-02,
 1.950394387167700994245891508369324694703820522489789125908612e-01,
 4.132929622783563686116108686666547082846741228042232731476147e-01,
 5.736690430342222603195557147853022060758392664086633396520345e-01,
 4.391583117891662321931477565794105633815363384084590559889493e-01,
 1.774076780986685727823533562031556893226571319881417676492595e-03,
 -3.263845936917800216385340830055349953447745005769416287177497e-01,
 -1.748399612893925042664835683606584215248582345438816346170042e-01,
 1.812918323111226960705459766025430918716233584167982942044424e-01,
 1.827554095896723746537533832033286839689931924709760567945595e-01,
 -1.043239002859270439148009137202747658420968144330108510179290e-01,
 -1.479771932752544935782314546369458188243947772922980064071205e-01,
 6.982318611329236513756591683950208955110603212379412334701145e-02,
 1.064824052498086303236593797715344405836015002929319291715777e-01,
 -5.344856168148319149493577269390074213960237013099439431132086e-02,
 -6.865475960403591525454725258715351280947435823354011140858001e-02,
 4.223218579637203541206570902753288247790857760067894456114927e-02,
 3.853571597111186425832144567362328142994885395255438867968781e-02,
 -3.137811036306775484244644776337594435094096964336402798072360e-02,
 -1.776090356835818354094298625884058170354129044259951019182732e-02,
 2.073492017996382475887790073068984224515077665517103399898854e-02,
 5.829580555318887971939315747596613038479561943085291072787359e-03,
 -1.178549790619302893728624468402138072504226527540325463847390e-02,
 -5.287383992626814439198630765217969804966319971038003993984480e-04,
 5.601947239423804853206514239940474788977188460452053462770324e-03,
 -9.390582504738289646165698675070641765810790863514339205205998e-04,
 -2.145530281567620980305401403432221668847980295600748913748902e-03,
 8.383488056543616046381924054554052104937784379435436426690560e-04,
 6.161382204574344193703789012696411561214682388271673214197731e-04,
 -4.319557074261807466712901913481943478521991611607433971794602e-04,
 -1.060574748283803889966150803551837402553866816191659959347053e-04,
 1.574795238607493590547765666590811258087715699737771458390360e-04,
 -5.277795493037868976293566636015627609248847457646525246271036e-06,
 -4.109673996391477816326502438997466532822639385119090230965252e-05,
 1.074221540872195031273584409245060623104931330938273936484593e-05,
 7.000078682964986734859102495210684809643657474253921074934684e-06,
 -3.887400161856795187587790410706550576033603097954065074023128e-06,
 -4.650463220640262639231145944536092973446596027469833860001618e-07,
 7.939210633709952088373459255067360793370284788682979065122810e-07,
 -1.079004237578671411922961583845716126060658213943840375162654e-07,
 -8.904466370168590769052983362721567202750591914741016835071257e-08,
 3.407795621290730008673832107214820587991557116806912418558069e-08,
 2.169328259850323106986222296525930099935873861026310788086221e-09,
 -3.776010478532324328184043667556576385639846460337894963138621e-09,
 6.780047245828636668305808192607091517605349478677442468580825e-10,
 1.002303191046526913509281844136258004034177309673269533418644e-10,
 -5.840408185341171468465492447799819262905317576847426970757700e-11,
 9.130510016371796243923232926650252570239054815939483900056681e-12,
 -5.251871224244435037810503452564279828539007071678724285717464e-13];
 
 var Daub27 = [
6.687131385431931734918880680779563307675740731544063787599480e-05,
1.205531231673213234251999812212394463872002561229330125152073e-03,
9.952588780876619771874091297340545740163119816300838847749336e-03,
4.945259998290488004302995584228917712171023349013386944893643e-02,
1.629220275023933206396286389387812803673796872000118325233533e-01,
3.671102141253898226423388094379126394383458407087000700420400e-01,
5.538498609904800487605460395549044755068663194750017660900436e-01,
4.934061226779989979265447084358038959373468583404767251300717e-01,
1.028408550618229112710739475157388764479351682549490307668477e-01,
-2.897168033145948463175311101489473923261698802610323264603418e-01,
-2.482645819032605667810198368127693701263349361209208170092197e-01,
1.148230195177853576326445213787661879970642975306605349249036e-01,
2.272732884141708265275037216925482827043581894357907763081103e-01,
-3.878641863180231062443346843661817078060143110529946543683356e-02,
-1.780317409590085821070366277249759321269342801053489323888575e-01,
1.579939746024048431173907799261019471878724997312653292884660e-02,
1.311979717171553289711406975836688896451835867594492827800969e-01,
-1.406275155580876537026622167053147161846397735962817855782362e-02,
-9.102290652956591798241345515773322449830692586525337562864481e-02,
1.731101826549371089085675445961947677452358872325373949295769e-02,
5.796940573471798814748840657698008349462526768238833307489106e-02,
-1.851249356199807710545837861298826718763077900221574749342712e-02,
-3.273906663102087145481936428049519742538150452785563039743756e-02,
1.614696692239566682272152627542980896527822528487665111124260e-02,
1.566559564892457873003263983940819950829497022298967052103291e-02,
-1.157718645897628140054089958116866381056430680879332334217267e-02,
-5.862096345462925972966025215266179082657169806555503857975278e-03,
6.856635609684880675273184141746359000591385833807880272568038e-03,
1.342626877303679609082208800217479591902967766971379107017011e-03,
-3.332854469520006162763300141047111065412307706449049389557931e-03,
1.457529625931728587128588244152604734177322144376309490881599e-04,
1.301177450244135139135787970279897042994109161268159963884641e-03,
-3.418351226915427611946547437228006377896519777431057005796358e-04,
-3.879018574101327604369144470124819695479087900682219330965466e-04,
2.019719879690326857104208791272390315160018069955787875123234e-04,
7.660058387068576876674274961751262847965101108848090019821555e-05,
-7.711145517797584208411720507329584053382646435270054267102827e-05,
-3.517483614907445391752737841583832374184046409747387149129674e-06,
2.063442647736885318487206413360228908558806028468062177953960e-05,
-3.901164070638425528170558032557368703418425915665413541985623e-06,
-3.657500908187104997045760131046655906827644494899206692043298e-06,
1.634369624725637835424610743915128591988676092276368687669255e-06,
3.050880686251999094242671997731089918322345713516567387655763e-07,
-3.472468147394389269364673179891460601330730511237974736379548e-07,
3.286558968055159530983261866450459360074591641809187825408848e-08,
4.026255052866908637178682747490340533992340623231336911661711e-08,
-1.321332273990056558848617809101876846857728483295631388083263e-08,
-1.309465606856955151282041809232358209226373823424148862843577e-09,
1.521614984778521740775073159445241799352681846880808663329946e-09,
-2.415526928011130660506395791946234018673860470542996426005750e-10,
-4.374986224293654395069947682013996351823060759948583134078918e-11,
2.213662088067662485181472969374945928903854605356443772873438e-11,
-3.295790122476585807069953975043096139541415768606924980926275e-12,
1.828188352882424933624530026056448539377272017834175009418822e-13];

var Daub28 = [
  4.710807775014051101066545468288837625869263629358873937759173e-05,
  8.794985159843870273564636742144073059158975665525081816488582e-04,
  7.542650377646859177160195786201116927568410621050693986450538e-03,
  3.909260811540534426092083794403768111329778710541126982205076e-02,
  1.351379142536410450770749411679708279921694061092200363031937e-01,
  3.225633612855224257318486139030596702170126503618082416187649e-01,
  5.249982316303355562348293243640252929543774162151269406404636e-01,
  5.305162934414858075256978195354516449402692654391295761050628e-01,
  2.001761440459844380384404537971725815970574972480152145882083e-01,
  -2.304989540475825257279397658067038304888129374484095837624889e-01,
  -3.013278095326417816909366061441334075444383937588485826752087e-01,
  3.285787916338710468450547883547348694255260871071954509422161e-02,
  2.458081513737595535752949960866466132239832334168533456626848e-01,
  3.690688531571127205290633425993077868843846977265847006108551e-02,
  -1.828773307329849166920408764650763092868965221608724574218473e-01,
  -4.683823374455167616514752420549419665215987106243491879971921e-02,
  1.346275679102260877490923315484152662987698625205479167761416e-01,
  3.447863127509970524678534595639646616244376966117385829345554e-02,
  -9.768535580565244174963692133038973587005628990493154911133358e-02,
  -1.734192283130589908795581592406238282930530566316914040035812e-02,
  6.774789550190933956165341752699717255041141690153626336867769e-02,
  3.448018955540951137600471926079622335842207388713342609755316e-03,
  -4.333336861608628393863254980828284403766309203453808666888800e-02,
  4.431732910062988320487418656322338284504389482966303454010563e-03,
  2.468806001015186586264188361362046240243934625858343309818244e-02,
  -6.815549764552309639259447104811254179605050667281644254737890e-03,
  -1.206359196821849005842466619530619474644989878503490321948471e-02,
  5.838816627748944864497370576838809711476027837762897602935327e-03,
  4.784863112454241718009916669120329848973107781600157214960003e-03,
  -3.725461247074254799171427871442937099025589672466088044410521e-03,
  -1.360373845639692436577650137133777929659265166644839235882291e-03,
  1.875998668202795626152766912508562385106168761893900192731562e-03,
  1.415672393140464257573780581396205840941849282748250523509874e-04,
  -7.486749559114629991320679819683227355746847370960399216568306e-04,
  1.154656063658921251969297916771881248142872975490882572741198e-04,
  2.295790982233456202366621544054366855729175050420515776344878e-04,
  -8.903901490044488099517361247378396756893227855233897357882978e-05,
  -4.907713416190250858324783990436748073854807494400738311968278e-05,
  3.641401211050802781223450761733180188911730291497201507086247e-05,
  4.638664981394294654002871426476885751050837817671843706915388e-06,
  -1.004326041333422601781848560432120920634648692782357855473103e-05,
  1.247900317574834146052381692752796047052443265982232422642017e-06,
  1.840363734517769191684379309039277810350620305330900536404818e-06,
  -6.670215479954892588747450458085225880096882699397256774967304e-07,
  -1.757461173209842779903676264971918635870906983281392939812547e-07,
  1.490660013535362170989340065033061951960933954388633507264360e-07,
  -8.262387315626556965966429243600984899650039704831080988658278e-09,
  -1.784138690875710077191713941441263246560738410213624546116655e-08,
  5.044047056383436444631252840057862002264087720676808580373667e-09,
  6.944540328946226952976704718677697525410051405055662575530111e-10,
  -6.077041247229010224760245305596307803830053533836849384680534e-10,
  8.492220011056382105461206077240377024404404638947591299761197e-11,
  1.867367263783390418963879146175452376940453585791428841004699e-11,
  -8.365490471258800799349289794397908900767054085216008197372193e-12,
  1.188850533405901520842321749021089497203940688882364518455403e-12,
  -6.367772354714857335632692092267254266368934590973693820942617e-14];
  
  var Daub29 = [
 3.318966279841524761813546359818075441349169975922439988843475e-05,
 6.409516803044434540833706729120596322083061716935004987374676e-04,
 5.702126517773375434760843998623507494914551464968126455168657e-03,
 3.077358022140837676716707336516751814713312018344719150923618e-02,
 1.113701169517405304762186166370327770191325772342190715118617e-01,
 2.806534559709829376968881262770480606500920092398534229615289e-01,
 4.897588047621993143592705932993573539235839610055331620240518e-01,
 5.513744327583751951223746071670135992466984391233429663886536e-01,
 2.891052383358291634605691113586264061513180158354460952469246e-01,
 -1.540287344599000542466293779503370141731339982919280951230240e-01,
 -3.300409489175880520295083779487012611959310539629627124613719e-01,
 -5.570680007294085781514541931715795784309410235726214400350351e-02,
 2.361052361530259415983110734054626770649468357328362426830433e-01,
 1.124191748731883764769740670535880543076817816861518667898467e-01,
 -1.608779885941877360771615465531852333085159940159968393590303e-01,
 -1.078459499387214201077881957354707913786241153934264316589273e-01,
 1.144722958938182579734135930060053286267822797640393386903440e-01,
 8.322074716244975790297348835032537357891920536002627784941129e-02,
 -8.512549261563550232832311331420804581881235448862834507281486e-02,
 -5.502748952532572320924541450626650067707344725344841099873446e-02,
 6.347916458421186633577789314698972361081611994794140119302163e-02,
 3.053154327270413646637328212093941030592133225231728964047047e-02,
 -4.518798127778834515979704475304405691390090327474972089790857e-02,
 -1.291714255426679462966473962555410660387671182428076570686472e-02,
 2.947043187174764111028122319949903667638786379520519899154373e-02,
 2.648327307678167915542397563479749119673768286990136051577167e-03,
 -1.704122457360668969234196743407615179099529206118693044741086e-02,
 1.737880332720511164430027824345354801611373419264590068097416e-03,
 8.469725493560752287772961661104710791306496373354237126998903e-03,
 -2.550807127789472659145072247724735637183590942511858255354005e-03,
 -3.473798989681100630649790255076233970957721666820195620598374e-03,
 1.877120925723650133179338154344873477230567340668548016358682e-03,
 1.087053942226062966738944397844498417945523630053411148182206e-03,
 -1.000778327085680541055696707760062870925897014530348262794137e-03,
 -2.000711363076779808296301110796026470163110202848894744316755e-04,
 4.111283454742767033424740543004041500054889660665367490129376e-04,
 -2.292018041214499897382298271438084577065170236103859181134525e-05,
 -1.293044840080720609161466939678226852440475312744714379499074e-04,
 3.645026068562774967665464216602750761690984830805534178557146e-05,
 2.913344750169041218495787251929571015775436967652945386217480e-05,
 -1.657328395306616289863396387854880512976861409870690029695161e-05,
 -3.593644804025187638066915189731950450034629392522542962477168e-06,
 4.750609246452552850197117564759363194953518317428400241629683e-06,
 -3.029054592052818286474228294307141792053791695855058563299597e-07,
 -8.975701750636280734511651941681818767895052287332471537510510e-07,
 2.633898386997696553900967704111473475368019612368922599394214e-07,
 9.387197411095863026484410601284876812292554863800653292318725e-08,
 -6.286156922010786166768503252870590953166867739448102804392389e-08,
 1.076591906619196137385201975028785139607670319821266803566785e-09,
 7.768978854770062238895964639391324551611701293594055935346266e-09,
 -1.893995386171984147774611076618946011337498790609031626697228e-09,
 -3.426800863263089001811012278889864200550342566386405676893537e-10,
 2.407099453509342962399811991929330725186626582891090462239366e-10,
 -2.940589250764532582888473974638273664244682541297835986306504e-11,
 -7.832509733627817032356556582819494794884131433810848844709881e-12,
 3.152762413370310423797539876893861621418382024668704492620948e-12,
 -4.285654870068344101898185073376307686875386259541180967347399e-13,
 2.219191311588302960934661700068023727737812918006011019184982e-14];
 
 var Daub30 = [
2.338616172731421471474407279894891960011661146356580425400538e-05,
4.666379504285509336662000111055365140848987563882199035322085e-04,
4.300797165048069510045016757402827408493482974782286966500398e-03,
2.413083267158837895194919987958311943976725005113561262334092e-02,
9.123830406701570679321575555085899708564500191080751595642650e-02,
2.420206709402140994467599658342919512318194032687898436229538e-01,
4.504878218533178366981351802898336415314944375740699506554771e-01,
5.575722329128364304078082520999850413492571645754785374629734e-01,
3.662426833716279793144871151369089533016299234992584741629624e-01,
-6.618367077593731501909741041813726474911212544474895441395148e-02,
-3.329669750208556069196849320598850505877494561268613506392514e-01,
-1.419685133300829310219026267403758254954270602825020111483505e-01,
1.994621215806643032428990062111230223523226088131364328774921e-01,
1.778298732448367361280250921330425046260289700971176750362566e-01,
-1.145582194327077814891518778613672243404957549114393749173137e-01,
-1.572368179599938126878197378886501553251711910617673398124611e-01,
7.277865897036442699893544326605244235248713804556715604416632e-02,
1.227477460450093778691578797698150091624353365248212907325446e-01,
-5.380646545825707676022015051837304300338645984615639237930800e-02,
-8.765869003638366048026572053699028353846982304851342479893827e-02,
4.380166467141773250305407710250135373016604593736480428415303e-02,
5.671236574473569492590636983030617493807140224924978946302257e-02,
-3.567339749675960965780819743176056734137251336781389369397564e-02,
-3.226375891935220815954913483392725682165778426411705216010280e-02,
2.707861959529418272206848318420006522973840949600186710327776e-02,
1.528796076985739546052896626042375110302102640936712142026221e-02,
-1.839974386811734118728169880549148389603890445324127330811811e-02,
-5.296859666131086629169938675330494864053932988161015674773617e-03,
1.091563165830488927536881480211929049886878831313700460017968e-02,
6.196717564977244383592534999284255315694546230739551683085460e-04,
-5.530730148192003288871383856487027893918513053091795443517653e-03,
8.433845866620933982126003584365932145598126087481400294999080e-04,
2.324520094060099304385756339638431339131122661576649123053845e-03,
-8.609276968110423879660725173525347077801305237644122054954659e-04,
-7.678782504380918697963922441514742758516706160788123977340073e-04,
5.050948239033467796256544554086554367969638627715114003635557e-04,
1.724825842351709725545759714374272164367933578194910678479473e-04,
-2.161718301169633804271038862087964094429005266172702380483361e-04,
-8.548305467584070994787824796256108217987765582429940610377190e-06,
6.982008370808327851082027193100914402221658444151889697045071e-05,
-1.339716863293971629296314599448901465078920406443516550195793e-05,
-1.636152478725426488654528710478856195004608401773950511915162e-05,
7.252145535890469015723401169934327900622894130695550273452916e-06,
2.327549098493686509557358103785598216688723737824121617676858e-06,
-2.187267676996166416699555236143059249832615777542412142603694e-06,
1.099474338526203304286307383463498542376432972308342428764576e-08,
4.261662326011572446469849114416378817419458434583398455985144e-07,
-1.000414682354500898864979332965559934104686157639553850670490e-07,
-4.764379965139453357729154748688006975561934425368712852985388e-08,
2.605442754977625431940885841950955928085338672381046225838880e-08,
5.553397861397053982967618072672572206490972606026556946910028e-10,
-3.331105680467578245901976412732595596538702049437802824373020e-09,
6.984862691832182584221096665570313611280449991512869846064780e-10,
1.613622978270904360610418704685783656905979134344922647926295e-10,
-9.461387997276802120884525814092001871993910062127702293573920e-11,
1.000105131393171192746337860330428369495110180346654025287492e-11,
3.239428638532286114355931428908079297696045600279108835760520e-12,
-1.185237592101582328254231496310584611948560976394420324137742e-12,
1.543997570847620046003616417646988780670333040868954794039905e-13,
-7.737942630954405708679963277418806436871098329050829841696327e-15];

var Daub31 = [
  1.648013386456140748122177817418358316441195236228590958603489e-05,
  3.394122037769956699157160165352942212213928231154233571163033e-04,
  3.236884068627721221829662672296912258338131668810067169630813e-03,
  1.885369161298591269159568944275763468999829139547989648553486e-02,
  7.433609301164788697908776495388047669378919816041031344650271e-02,
  2.070128744852353286198055444111916450619762837756134323019573e-01,
  4.091922000374278563928213235836188963704298775635493549519369e-01,
  5.511398409142754983590484577074663132074992263886810324421617e-01,
  4.294688082061372955430413148799008354573408538414331312236645e-01,
  2.716921249736946422305354732634261873401679092095992827198308e-02,
  -3.109551183195075186926560285811004715398678229333522634202008e-01,
  -2.179784855235633521693544507220105631639547435903112747133934e-01,
  1.401782887652732681656253206993073895422881511380152633441096e-01,
  2.249667114737370933697297905066886078307490136415302624018330e-01,
  -4.992634916046823977000579399730138693074543903234092797936484e-02,
  -1.869623608957154494374577196258383009208655076187653847079167e-01,
  1.543698842948893409652995335281236231845293548571166883219023e-02,
  1.450895009319931981518942907854879059128872873116921504156674e-01,
  -8.139832273469236863527708715566588550006680549152344840146851e-03,
  -1.076127733234956326668605511648013952380301953590447106075614e-01,
  1.094129745236496925725237900637802669504835743555466811796369e-02,
  7.535361174328140695528289751109133941376701984419452638686226e-02,
  -1.488002661810482202699555987503429289100801979910046913257306e-02,
  -4.861907546485433003537603385831190109391263542044516048871113e-02,
  1.615417156598591113619453864586701665635869166193865651960591e-02,
  2.804761936675616906861927211659154977049392281479113764697785e-02,
  -1.427627527776351943309800140756746087215016194775579070599004e-02,
  -1.390055293926652880755898888934447671732373519028670201124816e-02,
  1.051763948737184089128633441244991643331033825102031908858652e-02,
  5.516163573310992566561289762241160214476622662764637181816550e-03,
  -6.520852375874612553325469682628530079210293774541131381751695e-03,
  -1.428264223218909891400516038687842292177211292295049238921068e-03,
  3.393066776715931928419358796960612411097347419792355896915546e-03,
  -6.397901106014600492881202314307290077992972755016494062875201e-05,
  -1.459041741985160943114515221598080223845239255190055621901681e-03,
  3.431398296904734438118401084929505912208229684629857530009147e-04,
  4.998816175637222614896912406679513231966722440032799024979502e-04,
  -2.396583469402949615285646688069476140260781708006174912535660e-04,
  -1.243411617250228669409179807383399199879641177993453588807726e-04,
  1.089584350416766882738651833752634206358441308880869184416670e-04,
  1.501335727444532997071651937630983442758297688087711521441229e-05,
  -3.631255157860086164261313773172162991107348698083164489165837e-05,
  4.034520235184278839752741499546098778993926344831736074409765e-06,
  8.795301342692987765440618030678349427367022581211855857458220e-06,
  -3.035142365891509630069007852947057220760887215249503512783023e-06,
  -1.369060230942940782050489751987123955074404782177163471279285e-06,
  9.810015422044371573950976088058064384946146188110905321673802e-07,
  5.327250656974915426977440959783080593776012130063170688309127e-08,
  -1.975925129170206248152121156696590501303803187231928513867046e-07,
  3.616826517331004805247567218405798591329788122337274956172315e-08,
  2.328309713821409644308538888589329921141948539678106680777082e-08,
  -1.061529602150252306500404266150823962402673780484965538270541e-08,
  -6.474311687959861398702581539341954438747926255671605657095807e-10,
  1.408568151025177427076547804944585301332087108125727813194374e-09,
  -2.524043954153353306183643702933218308617979467184848456565837e-10,
  -7.348930032486263904766913919653624379586487437915175106407348e-11,
  3.692108808871129411604189196259677640440919369478263728899602e-11,
  -3.327008967125979929910636246337150851642079794871116041187279e-12,
  -1.324334917243963163878274345609465717294426628053460151843705e-12,
  4.445467096291932163298411852093011459626037560439178917611592e-13,
  -5.559442050579014337641375730083534521513818164827556763756543e-14,
  2.699382879762665647295493928801387173921314576598505507855504e-15];
  
  var Daub32 = [
 1.161463302135014885567464100760659332951431420121048996305591e-05,
 2.466566906380903352739104211274667134470169443886449124673996e-04,
 2.431261919572266100780423071905958127811969678055971488060574e-03,
 1.468104638141913563547809006402194831107662001343421893488086e-02,
 6.025749912033537081745451975527967031851677384078997261920024e-02,
 1.757507836394388988189299915753348505208376399651864661397588e-01,
 3.675096285973496361995340339143234125206079560406868595968025e-01,
 5.343179193409538322901117858552186425529774700290587495921679e-01,
 4.778091637339484033555130814414794130354053753675509287934741e-01,
 1.206305382656178269538098710665261299391507308342013788891222e-01,
 -2.666981814766755535489784087869865024226542605534080371507405e-01,
 -2.774215815584272153338153320303401666681294506143291967655666e-01,
 6.471335480551623831000090095167664918448659157720155321560811e-02,
 2.483106423568801736064852157222867588791898170114101300999760e-01,
 2.466244483969740441701479334808723214802614938081258920635302e-02,
 -1.921023447085468984341365278247990525863123891147783426068990e-01,
 -4.899511718467173853355943225576377418394280156945986899417475e-02,
 1.452320794752866460838830744051944832326998342053148426312341e-01,
 4.440490819993974022640619534046603571086531544468421519143629e-02,
 -1.094561131160893831027722774343269232755171130623890041619420e-01,
 -2.962787250844770491204452379051215505049068645551070779367843e-02,
 8.087414063848395744090831590426327690818854671836423275412813e-02,
 1.410615151610660772869738802931740150275269382463799031013905e-02,
 -5.692631406247843550478416271158537960555270097953330567652364e-02,
 -2.380264464932573834443178362086503847328134994591954135879789e-03,
 3.705145792354468010437633458013030898015496905609424004450953e-02,
 -4.145907660827218781460700428862611061267328108653649653634276e-03,
 -2.166282283639119347634778516947485598599029367518033869601702e-02,
 6.167527310685675112579059689520105004744367282412921739811164e-03,
 1.101740071540688116532806119564345712473051769079712407908648e-02,
 -5.411568257275791208581502410752383050600045942275647685361370e-03,
 -4.649216751184411528658094984504900172989190128905887602541396e-03,
 3.627224640687864960122122984391704782343548385375321260251988e-03,
 1.468955100468467772528811782840480639166582822577191079260543e-03,
 -1.964740555821778254183647540656746450092725858126595984907304e-03,
 -2.211678729579097916278097586914956834196749138610403102772710e-04,
 8.673058518450555343925662389563539890596549655683386287799624e-04,
 -1.024537310607396186949656796812972062290796122915930356634122e-04,
 -3.059654423826911750479261161552574500739091332121504634422577e-04,
 1.053915461739828114700905192091104141076083602686374410146603e-04,
 8.103678329134838389828091896334156224227821362491626044950428e-05,
 -5.259809282684322782648914338377962890245975842272425408122506e-05,
 -1.294045779405512723950480259110995722517019870286295908085366e-05,
 1.824268401980691220603850117995712615809177092802967489081228e-05,
 -6.361781532260254953363913076575914206506177493714496098327288e-07,
 -4.558309576264423135123964145585288808181431652781253437738445e-06,
 1.202889036321620990296134494079846952404216422923750605507047e-06,
 7.560047625595947819392627283726711361273296630256477108501994e-07,
 -4.285970693151457255418342315045357407199066350632593899896712e-07,
 -5.003361868748230293692887222336390314786090450819216035110269e-08,
 8.965966311957728376981484572655177545054433542721057470726361e-08,
 -1.219924359483373093110396748985081720383992859961285213840740e-08,
 -1.104383021722648979552131128575075255513372249283096583736746e-08,
 4.250422311980592983740943309197245384991941251563471671065543e-09,
 4.384387799940474369553236949848427579687147486892033587998023e-10,
 -5.881091462634605628881794361152305108432139465417759716875076e-10,
 8.904723796221605490455387579189371137903330749397374037644960e-11,
 3.263270741332907875981844980104948375955551273115386408552080e-11,
 -1.430918765169202320188022211739750594608742928641485026836608e-11,
 1.075610653501062115165734990153347111902874668945095034791947e-12,
 5.361482229611801638107331379599434078296259332654994508124989e-13,
 -1.663800489433402369889818192962259823988673359967722467427927e-13,
 2.000715303810524954375796020597627467104635766752154321244151e-14,
 -9.421019139535078421314655362291088223782497046057523323473331e-16];
 
 var Daub33 = [
8.186358314175091939858945975190102731733968885547217619434602e-06,
1.791016153702791479424389068736094134247294413108336017758506e-04,
1.822709435164084208084617771787691709255513374281497713580568e-03,
1.139594337458160925830840619716397130445853638888472948832932e-02,
4.861466653171619508385707681587366397164931431125053574327899e-02,
  1.481863131800528081784673514426737436792606299953305691300616e-01,
  3.267181301177075783930752787756046348844272437670999719562429e-01,
  5.093761725149396552227892926384090200953139820961482931291482e-01,
  5.112547705832674655425831875568453973369927971748064975152374e-01,
  2.095823507130554216526494469993023406452629154801126958766008e-01,
  -2.042026223985421049629055102642279430174095014493415546881477e-01,
  -3.159974107665602561905181464284910961862968513875028980451424e-01,
  -1.927833943695275915600583425408664108893845271616240406358226e-02,
  2.454206121192791114179964351253140999836791489738418857473689e-01,
  9.985155868033815698139640215477639365289384281516885362929979e-02,
  -1.714280990518593279308738113273443832545615219650436927029674e-01,
  -1.108441331167107910806084983056783194189909198734302929909672e-01,
  1.219678564037346149389134584371009777591763921148126952722200e-01,
  9.478808805061595889263191779090571160237408179346345390888721e-02,
  -9.114696835133148913093153757138373418923462847746880902676089e-02,
  -7.030248505405615921453280814171665167171986608963193275084895e-02,
  7.019114394099653254998935842432841393915841096633514680190145e-02,
  4.573456189389667743139040427641638967843459421665709740086516e-02,
  -5.347125133582228919431110824663168583260050383336359554980188e-02,
  -2.524858297747649929258392207837724793937727346177294684700378e-02,
  3.868706076024496481748675031852528047303323816250150793091832e-02,
  1.070326582001954942654534968137727769698168853186071888736311e-02,
  -2.572876175473297336123211392278301875687760837710204579628265e-02,
  -2.167758617353607324783298657172830203896433848418061622436727e-03,
  1.531695411585766548347442266431874060229304787191589430967538e-02,
  -1.594288782414604768637856446111392724059836934455189837500244e-03,
  -7.953540387057939240459305406538116220678495240302592677582773e-03,
  2.389062408165908575935815973439728988151836094753689966108405e-03,
  3.480800953405711999411461002429227385937942254778524257436278e-03,
  -1.860718214455795912074482150710567824317228203897000129729967e-03,
  -1.204309257604658876916644980097327372892008586047095719636829e-03,
  1.074380696351291355073899234941719080473877020595209197706651e-03,
  2.727305847336937211749282358350196461733595290569540045817329e-04,
  -4.908329007590351474487792254066540683724948757382104652497458e-04,
  4.393166251766185755059005296958129844094063524324718175254673e-06,
  1.780431898251245351831728023200069586928513661382622116969992e-04,
  -4.160438516273709306234368807933932360567787692918883118883736e-05,
  -4.929564423417301834310231482621574127409950921583062559483686e-05,
  2.423335398816890365621188379922041046073808819182024026589770e-05,
  9.070805757828453800203677464921508178468256685438211818575040e-06,
  -8.866121366757736169176034432364298134186929098274651022820760e-06,
  -3.607516102879771631230351118595069330196155459105589342866625e-07,
  2.288371276141527305481395545993763010565968667577768164201792e-06,
  -4.426923407952870147984002129341809185622768353983550670755106e-07,
  -3.985791291985944076942626511739220753169387460984290019185514e-07,
  1.822443332571053437467128998002798233969112236553215291639303e-07,
  3.377972703730854377516206663481869099376154259897212784144779e-08,
  -3.987838198518880722819502850814936369197384392561970319349663e-08,
  3.672863576838181340505563759379169099717712645283448779390320e-09,
  5.111211857347453839549366593998758891130921028374576213256027e-09,
  -1.671392677251932495173219614104411841891545601521784559793012e-09,
  -2.496402105246193648073519269370197331176405371538404298745013e-10,
  2.426833102305682309891302883361232297664099485514601790344279e-10,
  -3.049574453945863430361296931455141500128170151643206937547928e-11,
  -1.420236859889936792437077844940412749343225644487770840543290e-11,
  5.509414720765524548752673631197714447818740985929081064907524e-12,
  -3.343481218953278765982532722689984725170758193566174566492199e-13,
  -2.152488386833302618520603545685994753329478275805993737095214e-13,
  6.214740247174398315576214699577230693021307854673557214652751e-14,
  -7.196510545363322414033654470779070592316600780697558361083151e-15,
  3.289373678416306368625564108782095644036415401902518812978798e-16];
  
  var Daub34 = [
 5.770510632730285627466067796809329117324708919047900817738025e-06,
 1.299476200679530037833484815390569400369432658207722720405084e-04,
 1.364061390059049998200014449396877439591680435610837369411339e-03,
 8.819889403884978803182764563095879335330977939541630862804757e-03,
 3.904884135178594138905026219591569204043816577941517019631916e-02,
 1.241524821113768081954449898210969172708199672428635378051285e-01,
 2.877650592337145629334256618087718872558560120999651277991839e-01,
 4.784787462793710621468610706120519466268010329031345843336104e-01,
 5.305550996564631773133260223990794445605699030503652382795600e-01,
 2.903663295072749510455945186199530115755664977934564128822650e-01,
 -1.282468421744371672912377747048558427612774932943748628650824e-01,
 -3.315253015083869417715548463087537345035828886426345397256876e-01,
 -1.038919155156404718287260506925867970596448618647006698388596e-01,
 2.169072201874275950610018667099322465619408030256534197819784e-01,
 1.666017504122074437311574334509261366682993700573488534577890e-01,
 -1.273373582238011562843862636988693890108793629966541695807247e-01,
 -1.609249271778668063014799490429649196614628857267382976958607e-01,
 7.799184693794810738265349531832015087096882277333968473726399e-02,
 1.341259602711361284802399913977387999358280900708582462625539e-01,
 -5.448296806413904636632671383140642554265865948686157271017286e-02,
 -1.029475969928140852342073823689090498245496056845473569066667e-01,
 4.357609464963129726428486610925800727137724136370669421246609e-02,
 7.318523543679560555546221335452045680757998947493883124934567e-02,
 -3.701283841786244960356402125554190040750079009127461655784927e-02,
 -4.743855964527776247220681410983851377889756018716427358008296e-02,
 3.073974657395934459931226513844134346305562928466993208164603e-02,
 2.722835075635419610095839895805858855202745897718117731496534e-02,
 -2.367173792282636485046786438094940427456079528043555566867110e-02,
 -1.314398001665716086105827506126287041342680578404007359439612e-02,
 1.640937419986519252112261495537409592363156309874473310057471e-02,
 4.713649260999809905918876125437488856235874027077755004539205e-03,
 -1.004550670836151917439146861146431000364858401181337134891421e-02,
 -6.194748845153872839014356621835501857322345445234809347431098e-04,
 5.334950768759936032170270195983921511565539100791906952901398e-03,
 -7.692127975067836975989490900561029844887285335804349474993607e-04,
 -2.399453943537055863933124827688081952701780599883067560501870e-03,
 8.589959874363661955444898475746536583497522107459291718900058e-04,
 8.751999064078688732610570055224339733760304773327228476255647e-04,
 -5.527355762144197975516415296735124460550632283763688359649888e-04,
 -2.326732140233531635428863212833942245597361085708567528230733e-04,
 2.650772397558057819755811309071002543822145660933016957735937e-04,
 2.660050018453441903046828468025589086403126180798464347801678e-05,
 -9.914697770780134603580350758869378471802751837608461971022567e-05,
 1.353117227249649581251887376414486225127346352042209141315562e-05,
 2.844951419697807376503080001943765930601242225183893658540032e-05,
 -1.057657494257950623848316304755218120233253479317574337409622e-05,
 -5.710826510998303938275050074333400305512451419983646591762318e-06,
 4.169871758547028398316761659984928804362023643629741358799744e-06,
 4.979718101421307748081857636471761057429219265531618602960147e-07,
 -1.116306534817008428597995070751765080383261658112656948526954e-06,
 1.448195708333185127061180618150009526758658641231104901703561e-07,
 2.025990666667859216690536885693725545344933235432307649205497e-07,
 -7.526701740412589411177481797841044281662555785969415398369019e-08,
 -1.990346501531736915866180448337614967570744211158241514589121e-08,
 1.740423332936068076497051274445147160190783847854409836489662e-08,
 -8.665744261368722215864741166245385888818567571145958531936939e-10,
 -2.316501946995482751582294240136010067415084499025753117941001e-09,
 6.446378210323402313101214894500231181606520211579581132442548e-10,
 1.300410318609415248880403259300467720631189120978928377152233e-10,
 -9.904774537632409015479530333979124540183199174591377762845227e-11,
 1.004208735461769864836516428998306778031143650101842361622330e-11,
 6.080125354000167254059025929915591291115751734288584563131636e-12,
 -2.107879108915301546285370395443778864676275235126044599683271e-12,
 9.799451158211597727901178520526388692140586041163624252991805e-14,
 8.579194051799733179793112298652600511486581216528683482143106e-14,
 -2.317083703906408481078257081903089523234020423092175261925515e-14,
 2.587338381935699555813538163144986688834142571207152879144731e-15,
 -1.148944754480590128244815794312606245147888158018823490936280e-16];
 
 var Daub35 = [
4.067934061148559026665247110206084571051201477121972612218005e-06,
9.421469475576740631603027533116630224451049736050903361458759e-05,
1.019122680375098109319314672751485080202557607467199213778085e-03,
6.807292884319132011971333979015625113494050642797397817625326e-03,
3.123628851149071453063391210769353068187088999495893257051179e-02,
1.034044558614783789938787754929279183985553322796063517049140e-01,
2.513073789944933128513251971488905042866779761014740192816902e-01,
4.435927392240354378183910489448494594782039032807956294826105e-01,
5.370084275091661028670690231716974547580034932361053607723887e-01,
3.603456405180473278744458573988718422538114217890792270621563e-01,
-4.388388187393404111343479394097224312100349011932028865098625e-02,
-3.238228649121161212147302807993176715625480327235512530593160e-01,
-1.817869767667278325788350264528191676841493369460849123538616e-01,
1.660413574907809195438433327470947940538097914525298064477785e-01,
2.172992893210892977675493456199559114036326358517672106972956e-01,
-6.526287131067753892154895911331108284007380738865652420304233e-02,
-1.919195892985939528760786800798636198516495957924798820500876e-01,
1.930954466601835091947734585938109944647435243484967057775110e-02,
1.552924803962371144206753760712566993987319378965231186477630e-01,
-4.752680834111350445288110998030979143710864689041902167119118e-03,
-1.205855226433935545076589480704957722635324456812322150437989e-01,
4.734229172641948763293980314992213293971770695480616789828384e-03,
8.991354757072954417865374195261962983644048998218233900481856e-02,
-9.318558949903924837875002823617504227246562152671894579504378e-03,
-6.335603744044346612098887534020545705731671718057964802006671e-02,
1.322854958503655524455929847605110719648746890497356808289302e-02,
4.125469306470509212749750814299126656151504805845417994651417e-02,
-1.436683978422007182104025173214012797788904894291716373493525e-02,
-2.416949780166026740294880681731084091264533168816746227537030e-02,
1.276645671565674419403918018742432714973656598227939824940035e-02,
1.228943600811871086161967625814297050611100200023898377949151e-02,
-9.577797899235709998147309703713518608283233882793489733491642e-03,
-5.085991649233429881797636583578921194675393807761154549733547e-03,
6.137754586740521089596801883631921221145712545042519987641234e-03,
1.428088794070762107355585870669842132609159040625895090070111e-03,
-3.357644380922383229567732565298665639037348585961127075507937e-03,
7.615969435172736546769649923895317451534703066016116257300160e-06,
1.549637469702362975561719246539787717204438637997824935787688e-03,
-3.346692164250854961608526121524596908041109918361306282201310e-04,
-5.864810318991817532175809224131456738367101035694188223408841e-04,
2.648328819961289039302810122699710966048565368047575218693134e-04,
1.700012283661249043584690194716767771204207742625746308522935e-04,
-1.365883072261161602559926714744746422567509177443594045709653e-04,
-2.976995962848509743944225866488519668585242655980656646544319e-05,
5.304143122913310222538317980686374696005605533475685587486683e-05,
-2.437001526827789860990429478540556752694389693432668831073769e-06,
-1.572442077270281693663288966405861215692805972737981986121447e-05,
4.308047861716731191350493437937513220737450410132878032163179e-06,
3.353345862871309889390877168046133657377105681618708355266688e-06,
-1.895929617693153288493891051875444439753318548105998166574535e-06,
-3.903931733287306166657519468494511920760767388397825775326745e-07,
5.302368616904760917074352633915743250769600635829229600812520e-07,
-3.700308378205124537986402644918879149894035910106489082512364e-08,
-9.990396944534900755781728477561240762191443422318249128866740e-08,
3.008188650719066928230268918661718274504955045022550217051301e-08,
1.084902733789934825266560240100449884702749303326571747323086e-08,
-7.458116552893037631192407611262788593505988638365840409367117e-09,
5.897951310384361575470355861162022501172491937837712969865619e-11,
1.030823345485433383811700481488557422005210168069163779730908e-09,
-2.433545573751672936168877250405940817227367937230289801251648e-10,
-6.407938256501889018430608323235974406219193176918284664973727e-11,
4.000536627253744510742788201354093006471710416671002244302586e-11,
-3.125639357108557540598098228678150768528121565391376265627294e-12,
-2.567065476155081449204643852428401530283519685638256074752850e-12,
8.015088533687900921948605418789324826115616416343391081288979e-13,
-2.597954328893848084315198205094389145706680129208998638802995e-14,
-3.397720856796267431956783825659069596940335130100871912329556e-14,
8.624037434720089202680337663692777682810714650060805832406135e-15,
-9.298012529324185420921555664719863501848315099116725184370339e-16,
4.014628712333488654318569164614220308046021091178184654250982e-17];

var Daub36 = [
  2.867925182755946334630479473029238615535511775894262711054705e-06,
  6.826028678546358691748629102209605362240344266505035981791715e-05,
  7.602151099668488285869792677106082100141275054892389379198545e-04,
  5.240297377409884366201603524392995696042174937194435235003941e-03,
  2.489056564482796484885927333959115579403023347044729739255255e-02,
  8.565209259526409083864716995521111486437594750377856524772704e-02,
  2.177569530979008149637945915719999746248969705650625533415876e-01,
  4.064336977082553467407793990250384445903151630768558142125382e-01,
  5.322668952607286914777444748641462027213554723153906901129337e-01,
  4.178753356009697863620634559374236455222275302996931178265919e-01,
  4.397519752934862993862182898358763783110745559238982179690132e-02,
  -2.944210395891145711100715969898758940722458887377844633443675e-01,
  -2.468070369781255270524798278622698446566520718230313889086016e-01,
  9.811420416311477050518401371401568038943437322299913514049728e-02,
  2.465372776089742110529709111809595434656418762898152706621356e-01,
  7.278515095792229009687682299460382878643139026668958884429641e-03,
  -1.993372056086496198603363400094784142714162256792182570541036e-01,
  -4.586140074639271639145126228774831743002971373998329604574394e-02,
  1.541062366276428841776316300420654875883842819413623395358262e-01,
  5.027618007353842862036816972809884096761706036019748316890913e-02,
  -1.188037543101356316801816931383547446073152951044444224449501e-01,
  -3.988085357551317584091699967924044034100374257075864260934102e-02,
  9.115678225801654406336059281306715151058903055370522031843771e-02,
  2.503872144956848989919484296709846860569180993040383621980546e-02,
  -6.820901663681751124880436344265538690580358108714540763125119e-02,
  -1.131910031681742794381808082173695022123056280821611354577883e-02,
  4.851308354780908538616267662315735632292989749013261207046367e-02,
  1.424972661765391603147802607378542396323429657660009755652404e-03,
  -3.198072067763969654470293513742344601172739688274251641873778e-02,
  3.984040198717004857397179486790082321314291366656151213429068e-03,
  1.906359478062535932877576164368198274858108513696832728889209e-02,
  -5.657813245058818380424016973516714570499161434975761798379020e-03,
  -9.990263473281372348001743806489172665465685056975652497503772e-03,
  5.022989106665829004699819220796538830393945994687289792465541e-03,
  4.413484835350575251918616780287775585471012556848037301025999e-03,
  -3.484541445404883311209541395428535732697661971818727286003028e-03,
  -1.503074066296643749549363655363411879858070202740814054964603e-03,
  1.990793771851737270404293245701878186600899439513475823305914e-03,
  2.776812795712026068152384207605140383490242756921936501940389e-04,
  -9.463403823261101964604918059447913047725482130063492242779878e-04,
  8.614565758992702032613879159402330909634737204578606399403107e-05,
  3.693507284967510502620040341882236687749563414433432842567511e-04,
  -1.155118895843527096848376999413102395191976350936666573818799e-04,
  -1.131899468084665671727391922924411467938450743565106978099456e-04,
  6.694741196930590257104231749283786251555566773398199990337698e-05,
  2.375106683660860777161950832380341362257503761490580896617678e-05,
  -2.731390824654337912922346414722045404779935825834384250023192e-05,
  -1.183471059985615942783182762352360917304348034947412986608322e-06,
  8.372218198160788432628056043217491552198857358432112275253310e-06,
  -1.586145782434577495502614631566211839722879492827911790709498e-06,
  -1.870811602859180713762972281154953528056257451900381097476968e-06,
  8.311421279707778528163597405935375886855029592150424544500718e-07,
  2.548423522556577831218519052844387478819866531902854523544709e-07,
  -2.455377658434232699135878286794578515387138194247693201846263e-07,
  2.753249073339512254085076456700241929492720457889076058451072e-09,
  4.799043465450992009934526867650497683545716858606119786327559e-08,
  -1.156093688817008406756913949175208452083765368825442482226093e-08,
  -5.612784343327791397474114357094368557982413895802980814813369e-09,
  3.138841695782424018351567952158415003571380699236147752239001e-09,
  1.090815553713751810964713058800448676068475673611349566405716e-10,
  -4.512545778563249634425200856088490195004077806062978067796020e-10,
  8.962418203859611987065968320295929679774693465791367610044773e-11,
  3.037429098112535221800013609576297196061786927734556635696416e-11,
  -1.599716689261357143200396922409448515398648489795044468046420e-11,
  8.876846287217374213524399682895564055949886050748321818411161e-13,
  1.070969357114017002424433471621197579059927261727846375968378e-12,
  -3.029285026974877268896134589769473854669758797446795757329862e-13,
  5.542263182639804235231685861028995158694397223907295269180336e-15,
  1.338071386299105896025578761458472955294763310766371178363783e-14,
  -3.204628543401749860439316638848579711789176444320134355253750e-15,
  3.339971984818693213132578777712503670014459411167839211495237e-16,
  -1.403274175373190617489823209168013922564353495443487431242610e-17];
  
  var Daub37 = [
 2.022060862498392121815038335333633351464174415618614893795880e-06,
 4.942343750628132004714286117434454499485737947791397867195910e-05,
 5.662418377066724013768394373249439163518654840493603575144737e-04,
 4.024140368257286770702140124893772447952256842478891548092703e-03,
 1.976228615387959153244055502205017461538589475705618414896893e-02,
 7.058482597718160832030361890793007659963483925312132741868671e-02,
 1.873263318620649448028843491747601576761901656888288838192023e-01,
 3.684409724003061409445838616964941132670287724754729425204047e-01,
 5.181670408556228873104519667534437205387109579265718071174178e-01,
 4.622075536616057145505448401528172070050768534504278694229363e-01,
 1.308789632330201726057701201017649601034381070893275586898075e-01,
 -2.461804297610834132869018581145720710365433914584680691693717e-01,
 -2.943759152626617722808219575932673733674290772235644691367427e-01,
 1.967150045235938977077768648740052380288156507222647187301894e-02,
 2.515232543602686933435224095078166291442923992611593827552710e-01,
 8.180602838721862339029076982652411696000045533716726027662147e-02,
 -1.819622917786080007408824256525225216444443143868752611284260e-01,
 -1.084517138233017845554078812341876568514835176341639783558543e-01,
 1.299296469598537527842528895259188653120602318620944502979726e-01,
 1.017802968388141797470948228505865617480048287983176581607964e-01,
 -9.660754061668439030915405045955772715988585374771282291315496e-02,
 -8.233021190655740867404073660920379414988302492018783774702028e-02,
 7.504761994836017933579005072594245435071674452882148228583865e-02,
 5.956741087152995245435589042520108066877114768216272503684398e-02,
 -5.925681563265897095153806724965924334077555174281436189512239e-02,
 -3.825382947938424882011108885090442116802994193611884738133373e-02,
 4.580794415126833246633256156110381805848138158784734496981778e-02,
 2.097280059259754883313769469036393294461497749083921162354229e-02,
 -3.352358406410096994358662875913243067234786296009238949920582e-02,
 -8.833493890410232394064187990625563257107429109130726291528648e-03,
 2.261865154459947356571431658958802912061105608212828675323452e-02,
 1.690472383484423743663952859090705636512807161536954018400081e-03,
 -1.376398196289478433857985486097070339786225136728067000591187e-02,
 1.519305778833399218481261844599507408563295102235964076544334e-03,
 7.387757452855583640107787619408806919082115520707105052944171e-03,
 -2.248053187003824706127276829147166466869908326245810952521710e-03,
 -3.394523276408398601988475786247462646314228994098320665709345e-03,
 1.816871343801423525477184531347879515909226877688306010517914e-03,
 1.263934258117477182626760951047019242187910977671449470318766e-03,
 -1.111484865318630197259018233162929628309920117691177260742614e-03,
 -3.280788470880198419407186455190899535706232295554613820907245e-04,
 5.490532773373631230219769273898345809368332716288071475378651e-04,
 1.534439023195503211083338679106161291342621676983096723309776e-05,
 -2.208944032455493852493630802748509781675182699536797043565515e-04,
 4.336726125945695214852398433524024058216834313839357806404424e-05,
 7.055138782065465075838703109997365141906130284669094131032488e-05,
 -3.098662927619930052417611453170793938796310141219293329658062e-05,
 -1.639162496160583099236044020495877311072716199713679670940295e-05,
 1.354327718416781810683349121150634031343717637827354228989989e-05,
 1.849945003115590390789683032647334516600314304175482456338006e-06,
 -4.309941556597092389020622638271988877959028012481278949268461e-06,
 4.854731396996411681769911684430785681028852413859386141424939e-07,
 1.002121399297177629772998172241869405763288457224082581829033e-06,
 -3.494948603445727645895194867933547164628229076947330682199174e-07,
 -1.509885388671583553484927666148474078148724554849968758642331e-07,
 1.109031232216439389999036327867142640916239658806376290861690e-07,
 5.350657515461434290618742656970344024396382191417247602674540e-09,
 -2.252193836724805775389816424695618411834716065179297102428180e-08,
 4.224485706362419268050011630338101126995607958955688879525896e-09,
 2.793974465953982659829387370821677112004867350709951380622807e-09,
 -1.297205001469435139867686007585972538983682739297235604327668e-09,
 -1.031411129096974965677950646498153071722880698222864687038596e-10,
 1.946164894082315021308714557636277980079559327508927751052218e-10,
 -3.203398244123241367987902201268363088933939831689591684670080e-11,
 -1.398415715537641487959551682557483348661602836709278513081908e-11,
 6.334955440973913249611879065201632922100533284261000819747915e-12,
 -2.096363194234800541614775742755555713279549381264881030843258e-13,
 -4.421612409872105367333572734854401373201808896976552663098518e-13,
 1.138052830921439682522395208295427884729893377395129205716662e-13,
 -4.518889607463726394454509623712773172513778367070839294449849e-16,
 -5.243025691884205832260354503748325334301994904062750850180233e-15,
 1.189012387508252879928637969242590755033933791160383262132698e-15,
 -1.199280335852879554967035114674445327319437557227036460257649e-16,
 4.906615064935203694857690087429901193139905690549533773201453e-18];
 
 var Daub38 = [
1.425776641674131672055420247567865803211784397464191115245081e-06,
3.576251994264023012742569014888876217958307227940126418281357e-05,
4.211702664727116432247014444906469155300573201130549739553848e-04,
3.083088119253751774288740090262741910177322520624582862578292e-03,
1.563724934757215617277490102724080070486270026632620664785632e-02,
5.788994361285925649727664279317241952513246287766481213301801e-02,
1.600719935641106973482800861166599685169395465055048951307626e-01,
3.307757814110146511493637534404611754800768677041577030757306e-01,
4.965911753117180976599171147718708939352414838951726087564419e-01,
4.933560785171007975728485346997317064969513623594359091115804e-01,
2.130505713555785138286743353458562451255624665951160445122307e-01,
-1.828676677083358907975548507946239135218223185041410632924815e-01,
-3.216756378089978628483471725406916361929841940528189059002548e-01,
-6.226650604782432226643360160478765847565862101045597180310490e-02,
2.321259638353531085028708104285994998671615563662858079262996e-01,
1.499851196187170199586403453788927307298226028262603028635758e-01,
-1.417956859730596216710053144522330276392591055375830654519080e-01,
-1.599125651582443618288533214523534937804208844386102639177693e-01,
8.563812155615105741612217814369165313487129645536001850276987e-02,
1.414147340733826800884683119379170594092606174915755283496153e-01,
-5.658645863072738145681787657843320646815509410635114234947902e-02,
-1.147311707107443752394144019458942779715665489230169950201022e-01,
4.309589543304764288137871223616030624246568683595408792078602e-02,
8.720439826203975011910714164154456762073786124233088471855868e-02,
-3.660510340287429567372071039506772372567938710943432838908247e-02,
-6.176620870841315993604736705613246241897497782373337911398117e-02,
3.198987753153780630818381136366859026137035450576631134176875e-02,
4.005498110511594820952087086241114309038577379366732959648548e-02,
-2.689149388089451438550851767715967313417890393287236700072071e-02,
-2.311413402054931680856913553585621248925303865540203357180768e-02,
2.090464525565524340215982365351342094670261491526831672682244e-02,
1.129049727868596484270081487761544232851115891449843967151657e-02,
-1.470188206539868213708986402816605045648481224662435114088245e-02,
-4.131306656031089274123231103326745723188134548520938157995702e-03,
9.214785032197180512031534870181734003522861645903894504302286e-03,
5.625715748403532005741565594881148757066703437214522101740941e-04,
-5.071314509218348093935061417505663002006821323958752649640329e-03,
7.169821821064019257784165364894915621888541496773370435889585e-04,
2.400697781890973183892306914082592143984140550210130139535193e-03,
-8.448626665537775009068937851465856973251363010924003314643612e-04,
-9.424614077227377964015942271780098283910230639908018778588910e-04,
5.810759750532863662020321063678196633409555706981476723988312e-04,
2.817639250380670746018048967535608190123523180612961062603672e-04,
-3.031020460726611993600629020329784682496477106470427787747855e-04,
-4.555682696668420274688683005987764360677217149927938344795290e-05,
1.262043350166170705382346537131817701361522387904917335958705e-04,
-1.155409103833717192628479047983460953381959342642374175822863e-05,
-4.175141648540397797296325065775711309197411926289412468280801e-05,
1.334176149921350382547503457286060922218070031330137601427324e-05,
1.037359184045599795632258335010065103524959844966094870217687e-05,
-6.456730428469619160379910439617575420986972394137121953806236e-06,
-1.550844350118602575853380148525912999401292473185534395740371e-06,
2.149960269939665207789548199790770596890252405076394885606038e-06,
-8.487087586072593071869805266089426629606479876982221840833098e-08,
-5.187733738874144426008474683378542368066310000602823096009187e-07,
1.396377545508355481227961581059961184519872502493462010264633e-07,
8.400351046895965526933587176781279507953080669259318722910523e-08,
-4.884757937459286762082185411608763964041010392101914854918157e-08,
-5.424274800287298511126684174854414928447521710664476410973981e-09,
1.034704539274858480924046490952803937328239537222908159451039e-08,
-1.436329487795135706854539856979275911183628476521636251660849e-09,
-1.349197753983448821850381770889786301246741304307934955997111e-09,
5.261132557357598494535766638772624572100332209198979659077082e-10,
6.732336490189308685740626964182623159759767536724844030164551e-11,
-8.278256522538134727330692938158991115335384611795874767521731e-11,
1.101692934599454551150832622160224231280195362919498540913658e-11,
6.291537317039508581580913620859140835852886308989584198166174e-12,
-2.484789237563642857043361214502760723611468591833262675852242e-12,
2.626496504065252070488282876470525379851429538389481576454618e-14,
1.808661236274530582267084846343959377085922019067808145635263e-13,
-4.249817819571463006966616371554206572863122562744916796556474e-14,
-4.563397162127373109101691643047923747796563449194075621854491e-16,
2.045099676788988907802272564402310095398641092819367167252952e-15,
-4.405307042483461342449027139838301611006835285455050155842865e-16,
4.304596839558790016251867477122791508849697688058169053134463e-17,
-1.716152451088744188732404281737964277713026087224248235541071e-18];


var daubCoeff = {
  D2: Daub1,
  D4: Daub2,
  D6: Daub3,
  D8: Daub4,
  D10: Daub5,
  D12: Daub6,
  D14: Daub7,
  D16: Daub8,
  D18: Daub9,
  D20: Daub10,
  D22: Daub11,
  D24: Daub12,
  D26: Daub13,
  D28: Daub14,
  D30: Daub15,
  D32: Daub16,
  D34: Daub17,
  D36: Daub18,
  D38: Daub19,
  D40: Daub20,
  D42: Daub21,
  D44: Daub22,
  D46: Daub23,
  D48: Daub24,
  D50: Daub25,
  D52: Daub26,
  D54: Daub27,
  D56: Daub28,
  D58: Daub29,
  D60: Daub30,
  D62: Daub31,
  D64: Daub32,
  D66: Daub33,
  D68: Daub34,
  D70: Daub35,
  D72: Daub36,
  D74: Daub37,
  D76: Daub38
};

module.exports = daubCoeff;

},{}],15:[function(require,module,exports){

var haar = function(signal, options){
  var input = signal.values();
  var copy = [];
  var res = [];
  var len = Math.floor(signal.length / 2);
  while(len > 0){
    for(var i=0; i<len; i++){
      var scaling = (input[2*i] + input[2*i+1])*0.5;
      var wavelet = (input[2*i] - input[2*i+1])*0.5;
      copy[i] = scaling;
      res[len + i] = wavelet;
    }
    var tmp = copy;
    copy = input;
    input = tmp;
    len = Math.floor(len / 2);
  }
  res[0] = input[0];
  return res;
}

module.exports = {
  name: "DWT Haar",
  register: function(Transform){
    Transform.register({
      name: "haar",
      forward: function(signal, options){
        var trans = haar(signal);
        trans.timeDependent = true;
        return trans;
      },
      backward: function(spectrum, options){
        options.length = spectrum.length;
        var dft = prepare(options);
        dft.backward(signal);
        throw new Error("?? what to return ??");
      }
    });
  }
};

},{}],16:[function(require,module,exports){

module.exports = function(number){
  // bit check. every power of two has only one 1 in its binary representation
  // e.g. 32 = 100000
  // and a power of 2 minus 1 has all the lower bits set
  // e.g. 31 = 011111
  // the following statement uses a binary AND and it is always false for a 
  // non power of two, and true for a power of two
  // (it's probably the fastes way to do it)
  return ((number!=0) && !(number & (number-1)));
}

},{}],17:[function(require,module,exports){

module.exports = {
  triangle: function(period){
    return function(idx){
      return 1 - (Math.abs(((idx + period*0.25) % (period)) - period*0.5)) / period*4;
    }
  },
  rectangle: function(period){
    return function(idx){
      return (idx % period) < (period*0.5) ? 1 : -1;
    }
  }
}

},{}],18:[function(require,module,exports){
(function (global){
/**
 * @license
 * Lo-Dash 2.4.2 (Custom Build) <https://lodash.com/>
 * Build: `lodash modern -o ./dist/lodash.js`
 * Copyright 2012-2013 The Dojo Foundation <http://dojofoundation.org/>
 * Based on Underscore.js 1.5.2 <http://underscorejs.org/LICENSE>
 * Copyright 2009-2013 Jeremy Ashkenas, DocumentCloud and Investigative Reporters & Editors
 * Available under MIT license <https://lodash.com/license>
 */
;(function() {

  /** Used as a safe reference for `undefined` in pre ES5 environments */
  var undefined;

  /** Used to pool arrays and objects used internally */
  var arrayPool = [],
      objectPool = [];

  /** Used to generate unique IDs */
  var idCounter = 0;

  /** Used to prefix keys to avoid issues with `__proto__` and properties on `Object.prototype` */
  var keyPrefix = +new Date + '';

  /** Used as the size when optimizations are enabled for large arrays */
  var largeArraySize = 75;

  /** Used as the max size of the `arrayPool` and `objectPool` */
  var maxPoolSize = 40;

  /** Used to detect and test whitespace */
  var whitespace = (
    // whitespace
    ' \t\x0B\f\xA0\ufeff' +

    // line terminators
    '\n\r\u2028\u2029' +

    // unicode category "Zs" space separators
    '\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000'
  );

  /** Used to match empty string literals in compiled template source */
  var reEmptyStringLeading = /\b__p \+= '';/g,
      reEmptyStringMiddle = /\b(__p \+=) '' \+/g,
      reEmptyStringTrailing = /(__e\(.*?\)|\b__t\)) \+\n'';/g;

  /**
   * Used to match ES6 template delimiters
   * http://people.mozilla.org/~jorendorff/es6-draft.html#sec-literals-string-literals
   */
  var reEsTemplate = /\$\{([^\\}]*(?:\\.[^\\}]*)*)\}/g;

  /** Used to match regexp flags from their coerced string values */
  var reFlags = /\w*$/;

  /** Used to detected named functions */
  var reFuncName = /^\s*function[ \n\r\t]+\w/;

  /** Used to match "interpolate" template delimiters */
  var reInterpolate = /<%=([\s\S]+?)%>/g;

  /** Used to match leading whitespace and zeros to be removed */
  var reLeadingSpacesAndZeros = RegExp('^[' + whitespace + ']*0+(?=.$)');

  /** Used to ensure capturing order of template delimiters */
  var reNoMatch = /($^)/;

  /** Used to detect functions containing a `this` reference */
  var reThis = /\bthis\b/;

  /** Used to match unescaped characters in compiled string literals */
  var reUnescapedString = /['\n\r\t\u2028\u2029\\]/g;

  /** Used to assign default `context` object properties */
  var contextProps = [
    'Array', 'Boolean', 'Date', 'Function', 'Math', 'Number', 'Object',
    'RegExp', 'String', '_', 'attachEvent', 'clearTimeout', 'isFinite', 'isNaN',
    'parseInt', 'setTimeout'
  ];

  /** Used to make template sourceURLs easier to identify */
  var templateCounter = 0;

  /** `Object#toString` result shortcuts */
  var argsClass = '[object Arguments]',
      arrayClass = '[object Array]',
      boolClass = '[object Boolean]',
      dateClass = '[object Date]',
      funcClass = '[object Function]',
      numberClass = '[object Number]',
      objectClass = '[object Object]',
      regexpClass = '[object RegExp]',
      stringClass = '[object String]';

  /** Used to identify object classifications that `_.clone` supports */
  var cloneableClasses = {};
  cloneableClasses[funcClass] = false;
  cloneableClasses[argsClass] = cloneableClasses[arrayClass] =
  cloneableClasses[boolClass] = cloneableClasses[dateClass] =
  cloneableClasses[numberClass] = cloneableClasses[objectClass] =
  cloneableClasses[regexpClass] = cloneableClasses[stringClass] = true;

  /** Used as an internal `_.debounce` options object */
  var debounceOptions = {
    'leading': false,
    'maxWait': 0,
    'trailing': false
  };

  /** Used as the property descriptor for `__bindData__` */
  var descriptor = {
    'configurable': false,
    'enumerable': false,
    'value': null,
    'writable': false
  };

  /** Used to determine if values are of the language type Object */
  var objectTypes = {
    'boolean': false,
    'function': true,
    'object': true,
    'number': false,
    'string': false,
    'undefined': false
  };

  /** Used to escape characters for inclusion in compiled string literals */
  var stringEscapes = {
    '\\': '\\',
    "'": "'",
    '\n': 'n',
    '\r': 'r',
    '\t': 't',
    '\u2028': 'u2028',
    '\u2029': 'u2029'
  };

  /** Used as a reference to the global object */
  var root = (objectTypes[typeof window] && window) || this;

  /** Detect free variable `exports` */
  var freeExports = objectTypes[typeof exports] && exports && !exports.nodeType && exports;

  /** Detect free variable `module` */
  var freeModule = objectTypes[typeof module] && module && !module.nodeType && module;

  /** Detect the popular CommonJS extension `module.exports` */
  var moduleExports = freeModule && freeModule.exports === freeExports && freeExports;

  /** Detect free variable `global` from Node.js or Browserified code and use it as `root` */
  var freeGlobal = objectTypes[typeof global] && global;
  if (freeGlobal && (freeGlobal.global === freeGlobal || freeGlobal.window === freeGlobal)) {
    root = freeGlobal;
  }

  /*--------------------------------------------------------------------------*/

  /**
   * The base implementation of `_.indexOf` without support for binary searches
   * or `fromIndex` constraints.
   *
   * @private
   * @param {Array} array The array to search.
   * @param {*} value The value to search for.
   * @param {number} [fromIndex=0] The index to search from.
   * @returns {number} Returns the index of the matched value or `-1`.
   */
  function baseIndexOf(array, value, fromIndex) {
    var index = (fromIndex || 0) - 1,
        length = array ? array.length : 0;

    while (++index < length) {
      if (array[index] === value) {
        return index;
      }
    }
    return -1;
  }

  /**
   * An implementation of `_.contains` for cache objects that mimics the return
   * signature of `_.indexOf` by returning `0` if the value is found, else `-1`.
   *
   * @private
   * @param {Object} cache The cache object to inspect.
   * @param {*} value The value to search for.
   * @returns {number} Returns `0` if `value` is found, else `-1`.
   */
  function cacheIndexOf(cache, value) {
    var type = typeof value;
    cache = cache.cache;

    if (type == 'boolean' || value == null) {
      return cache[value] ? 0 : -1;
    }
    if (type != 'number' && type != 'string') {
      type = 'object';
    }
    var key = type == 'number' ? value : keyPrefix + value;
    cache = (cache = cache[type]) && cache[key];

    return type == 'object'
      ? (cache && baseIndexOf(cache, value) > -1 ? 0 : -1)
      : (cache ? 0 : -1);
  }

  /**
   * Adds a given value to the corresponding cache object.
   *
   * @private
   * @param {*} value The value to add to the cache.
   */
  function cachePush(value) {
    var cache = this.cache,
        type = typeof value;

    if (type == 'boolean' || value == null) {
      cache[value] = true;
    } else {
      if (type != 'number' && type != 'string') {
        type = 'object';
      }
      var key = type == 'number' ? value : keyPrefix + value,
          typeCache = cache[type] || (cache[type] = {});

      if (type == 'object') {
        (typeCache[key] || (typeCache[key] = [])).push(value);
      } else {
        typeCache[key] = true;
      }
    }
  }

  /**
   * Used by `_.max` and `_.min` as the default callback when a given
   * collection is a string value.
   *
   * @private
   * @param {string} value The character to inspect.
   * @returns {number} Returns the code unit of given character.
   */
  function charAtCallback(value) {
    return value.charCodeAt(0);
  }

  /**
   * Used by `sortBy` to compare transformed `collection` elements, stable sorting
   * them in ascending order.
   *
   * @private
   * @param {Object} a The object to compare to `b`.
   * @param {Object} b The object to compare to `a`.
   * @returns {number} Returns the sort order indicator of `1` or `-1`.
   */
  function compareAscending(a, b) {
    var ac = a.criteria,
        bc = b.criteria,
        index = -1,
        length = ac.length;

    while (++index < length) {
      var value = ac[index],
          other = bc[index];

      if (value !== other) {
        if (value > other || typeof value == 'undefined') {
          return 1;
        }
        if (value < other || typeof other == 'undefined') {
          return -1;
        }
      }
    }
    // Fixes an `Array#sort` bug in the JS engine embedded in Adobe applications
    // that causes it, under certain circumstances, to return the same value for
    // `a` and `b`. See https://github.com/jashkenas/underscore/pull/1247
    //
    // This also ensures a stable sort in V8 and other engines.
    // See http://code.google.com/p/v8/issues/detail?id=90
    return a.index - b.index;
  }

  /**
   * Creates a cache object to optimize linear searches of large arrays.
   *
   * @private
   * @param {Array} [array=[]] The array to search.
   * @returns {null|Object} Returns the cache object or `null` if caching should not be used.
   */
  function createCache(array) {
    var index = -1,
        length = array.length,
        first = array[0],
        mid = array[(length / 2) | 0],
        last = array[length - 1];

    if (first && typeof first == 'object' &&
        mid && typeof mid == 'object' && last && typeof last == 'object') {
      return false;
    }
    var cache = getObject();
    cache['false'] = cache['null'] = cache['true'] = cache['undefined'] = false;

    var result = getObject();
    result.array = array;
    result.cache = cache;
    result.push = cachePush;

    while (++index < length) {
      result.push(array[index]);
    }
    return result;
  }

  /**
   * Used by `template` to escape characters for inclusion in compiled
   * string literals.
   *
   * @private
   * @param {string} match The matched character to escape.
   * @returns {string} Returns the escaped character.
   */
  function escapeStringChar(match) {
    return '\\' + stringEscapes[match];
  }

  /**
   * Gets an array from the array pool or creates a new one if the pool is empty.
   *
   * @private
   * @returns {Array} The array from the pool.
   */
  function getArray() {
    return arrayPool.pop() || [];
  }

  /**
   * Gets an object from the object pool or creates a new one if the pool is empty.
   *
   * @private
   * @returns {Object} The object from the pool.
   */
  function getObject() {
    return objectPool.pop() || {
      'array': null,
      'cache': null,
      'criteria': null,
      'false': false,
      'index': 0,
      'null': false,
      'number': null,
      'object': null,
      'push': null,
      'string': null,
      'true': false,
      'undefined': false,
      'value': null
    };
  }

  /**
   * Releases the given array back to the array pool.
   *
   * @private
   * @param {Array} [array] The array to release.
   */
  function releaseArray(array) {
    array.length = 0;
    if (arrayPool.length < maxPoolSize) {
      arrayPool.push(array);
    }
  }

  /**
   * Releases the given object back to the object pool.
   *
   * @private
   * @param {Object} [object] The object to release.
   */
  function releaseObject(object) {
    var cache = object.cache;
    if (cache) {
      releaseObject(cache);
    }
    object.array = object.cache = object.criteria = object.object = object.number = object.string = object.value = null;
    if (objectPool.length < maxPoolSize) {
      objectPool.push(object);
    }
  }

  /**
   * Slices the `collection` from the `start` index up to, but not including,
   * the `end` index.
   *
   * Note: This function is used instead of `Array#slice` to support node lists
   * in IE < 9 and to ensure dense arrays are returned.
   *
   * @private
   * @param {Array|Object|string} collection The collection to slice.
   * @param {number} start The start index.
   * @param {number} end The end index.
   * @returns {Array} Returns the new array.
   */
  function slice(array, start, end) {
    start || (start = 0);
    if (typeof end == 'undefined') {
      end = array ? array.length : 0;
    }
    var index = -1,
        length = end - start || 0,
        result = Array(length < 0 ? 0 : length);

    while (++index < length) {
      result[index] = array[start + index];
    }
    return result;
  }

  /*--------------------------------------------------------------------------*/

  /**
   * Create a new `lodash` function using the given context object.
   *
   * @static
   * @memberOf _
   * @category Utilities
   * @param {Object} [context=root] The context object.
   * @returns {Function} Returns the `lodash` function.
   */
  function runInContext(context) {
    // Avoid issues with some ES3 environments that attempt to use values, named
    // after built-in constructors like `Object`, for the creation of literals.
    // ES5 clears this up by stating that literals must use built-in constructors.
    // See http://es5.github.io/#x11.1.5.
    context = context ? _.defaults(root.Object(), context, _.pick(root, contextProps)) : root;

    /** Native constructor references */
    var Array = context.Array,
        Boolean = context.Boolean,
        Date = context.Date,
        Function = context.Function,
        Math = context.Math,
        Number = context.Number,
        Object = context.Object,
        RegExp = context.RegExp,
        String = context.String,
        TypeError = context.TypeError;

    /**
     * Used for `Array` method references.
     *
     * Normally `Array.prototype` would suffice, however, using an array literal
     * avoids issues in Narwhal.
     */
    var arrayRef = [];

    /** Used for native method references */
    var objectProto = Object.prototype;

    /** Used to restore the original `_` reference in `noConflict` */
    var oldDash = context._;

    /** Used to resolve the internal [[Class]] of values */
    var toString = objectProto.toString;

    /** Used to detect if a method is native */
    var reNative = RegExp('^' +
      String(toString)
        .replace(/[.*+?^${}()|[\]\\]/g, '\\$&')
        .replace(/toString| for [^\]]+/g, '.*?') + '$'
    );

    /** Native method shortcuts */
    var ceil = Math.ceil,
        clearTimeout = context.clearTimeout,
        floor = Math.floor,
        fnToString = Function.prototype.toString,
        getPrototypeOf = isNative(getPrototypeOf = Object.getPrototypeOf) && getPrototypeOf,
        hasOwnProperty = objectProto.hasOwnProperty,
        push = arrayRef.push,
        setTimeout = context.setTimeout,
        splice = arrayRef.splice,
        unshift = arrayRef.unshift;

    /** Used to set meta data on functions */
    var defineProperty = (function() {
      // IE 8 only accepts DOM elements
      try {
        var o = {},
            func = isNative(func = Object.defineProperty) && func,
            result = func(o, o, o) && func;
      } catch(e) { }
      return result;
    }());

    /* Native method shortcuts for methods with the same name as other `lodash` methods */
    var nativeCreate = isNative(nativeCreate = Object.create) && nativeCreate,
        nativeIsArray = isNative(nativeIsArray = Array.isArray) && nativeIsArray,
        nativeIsFinite = context.isFinite,
        nativeIsNaN = context.isNaN,
        nativeKeys = isNative(nativeKeys = Object.keys) && nativeKeys,
        nativeMax = Math.max,
        nativeMin = Math.min,
        nativeParseInt = context.parseInt,
        nativeRandom = Math.random;

    /** Used to lookup a built-in constructor by [[Class]] */
    var ctorByClass = {};
    ctorByClass[arrayClass] = Array;
    ctorByClass[boolClass] = Boolean;
    ctorByClass[dateClass] = Date;
    ctorByClass[funcClass] = Function;
    ctorByClass[objectClass] = Object;
    ctorByClass[numberClass] = Number;
    ctorByClass[regexpClass] = RegExp;
    ctorByClass[stringClass] = String;

    /*--------------------------------------------------------------------------*/

    /**
     * Creates a `lodash` object which wraps the given value to enable intuitive
     * method chaining.
     *
     * In addition to Lo-Dash methods, wrappers also have the following `Array` methods:
     * `concat`, `join`, `pop`, `push`, `reverse`, `shift`, `slice`, `sort`, `splice`,
     * and `unshift`
     *
     * Chaining is supported in custom builds as long as the `value` method is
     * implicitly or explicitly included in the build.
     *
     * The chainable wrapper functions are:
     * `after`, `assign`, `bind`, `bindAll`, `bindKey`, `chain`, `compact`,
     * `compose`, `concat`, `countBy`, `create`, `createCallback`, `curry`,
     * `debounce`, `defaults`, `defer`, `delay`, `difference`, `filter`, `flatten`,
     * `forEach`, `forEachRight`, `forIn`, `forInRight`, `forOwn`, `forOwnRight`,
     * `functions`, `groupBy`, `indexBy`, `initial`, `intersection`, `invert`,
     * `invoke`, `keys`, `map`, `max`, `memoize`, `merge`, `min`, `object`, `omit`,
     * `once`, `pairs`, `partial`, `partialRight`, `pick`, `pluck`, `pull`, `push`,
     * `range`, `reject`, `remove`, `rest`, `reverse`, `shuffle`, `slice`, `sort`,
     * `sortBy`, `splice`, `tap`, `throttle`, `times`, `toArray`, `transform`,
     * `union`, `uniq`, `unshift`, `unzip`, `values`, `where`, `without`, `wrap`,
     * and `zip`
     *
     * The non-chainable wrapper functions are:
     * `clone`, `cloneDeep`, `contains`, `escape`, `every`, `find`, `findIndex`,
     * `findKey`, `findLast`, `findLastIndex`, `findLastKey`, `has`, `identity`,
     * `indexOf`, `isArguments`, `isArray`, `isBoolean`, `isDate`, `isElement`,
     * `isEmpty`, `isEqual`, `isFinite`, `isFunction`, `isNaN`, `isNull`, `isNumber`,
     * `isObject`, `isPlainObject`, `isRegExp`, `isString`, `isUndefined`, `join`,
     * `lastIndexOf`, `mixin`, `noConflict`, `parseInt`, `pop`, `random`, `reduce`,
     * `reduceRight`, `result`, `shift`, `size`, `some`, `sortedIndex`, `runInContext`,
     * `template`, `unescape`, `uniqueId`, and `value`
     *
     * The wrapper functions `first` and `last` return wrapped values when `n` is
     * provided, otherwise they return unwrapped values.
     *
     * Explicit chaining can be enabled by using the `_.chain` method.
     *
     * @name _
     * @constructor
     * @category Chaining
     * @param {*} value The value to wrap in a `lodash` instance.
     * @returns {Object} Returns a `lodash` instance.
     * @example
     *
     * var wrapped = _([1, 2, 3]);
     *
     * // returns an unwrapped value
     * wrapped.reduce(function(sum, num) {
     *   return sum + num;
     * });
     * // => 6
     *
     * // returns a wrapped value
     * var squares = wrapped.map(function(num) {
     *   return num * num;
     * });
     *
     * _.isArray(squares);
     * // => false
     *
     * _.isArray(squares.value());
     * // => true
     */
    function lodash(value) {
      // don't wrap if already wrapped, even if wrapped by a different `lodash` constructor
      return (value && typeof value == 'object' && !isArray(value) && hasOwnProperty.call(value, '__wrapped__'))
       ? value
       : new lodashWrapper(value);
    }

    /**
     * A fast path for creating `lodash` wrapper objects.
     *
     * @private
     * @param {*} value The value to wrap in a `lodash` instance.
     * @param {boolean} chainAll A flag to enable chaining for all methods
     * @returns {Object} Returns a `lodash` instance.
     */
    function lodashWrapper(value, chainAll) {
      this.__chain__ = !!chainAll;
      this.__wrapped__ = value;
    }
    // ensure `new lodashWrapper` is an instance of `lodash`
    lodashWrapper.prototype = lodash.prototype;

    /**
     * An object used to flag environments features.
     *
     * @static
     * @memberOf _
     * @type Object
     */
    var support = lodash.support = {};

    /**
     * Detect if functions can be decompiled by `Function#toString`
     * (all but PS3 and older Opera mobile browsers & avoided in Windows 8 apps).
     *
     * @memberOf _.support
     * @type boolean
     */
    support.funcDecomp = !isNative(context.WinRTError) && reThis.test(runInContext);

    /**
     * Detect if `Function#name` is supported (all but IE).
     *
     * @memberOf _.support
     * @type boolean
     */
    support.funcNames = typeof Function.name == 'string';

    /**
     * By default, the template delimiters used by Lo-Dash are similar to those in
     * embedded Ruby (ERB). Change the following template settings to use alternative
     * delimiters.
     *
     * @static
     * @memberOf _
     * @type Object
     */
    lodash.templateSettings = {

      /**
       * Used to detect `data` property values to be HTML-escaped.
       *
       * @memberOf _.templateSettings
       * @type RegExp
       */
      'escape': /<%-([\s\S]+?)%>/g,

      /**
       * Used to detect code to be evaluated.
       *
       * @memberOf _.templateSettings
       * @type RegExp
       */
      'evaluate': /<%([\s\S]+?)%>/g,

      /**
       * Used to detect `data` property values to inject.
       *
       * @memberOf _.templateSettings
       * @type RegExp
       */
      'interpolate': reInterpolate,

      /**
       * Used to reference the data object in the template text.
       *
       * @memberOf _.templateSettings
       * @type string
       */
      'variable': '',

      /**
       * Used to import variables into the compiled template.
       *
       * @memberOf _.templateSettings
       * @type Object
       */
      'imports': {

        /**
         * A reference to the `lodash` function.
         *
         * @memberOf _.templateSettings.imports
         * @type Function
         */
        '_': lodash
      }
    };

    /*--------------------------------------------------------------------------*/

    /**
     * The base implementation of `_.bind` that creates the bound function and
     * sets its meta data.
     *
     * @private
     * @param {Array} bindData The bind data array.
     * @returns {Function} Returns the new bound function.
     */
    function baseBind(bindData) {
      var func = bindData[0],
          partialArgs = bindData[2],
          thisArg = bindData[4];

      function bound() {
        // `Function#bind` spec
        // http://es5.github.io/#x15.3.4.5
        if (partialArgs) {
          // avoid `arguments` object deoptimizations by using `slice` instead
          // of `Array.prototype.slice.call` and not assigning `arguments` to a
          // variable as a ternary expression
          var args = slice(partialArgs);
          push.apply(args, arguments);
        }
        // mimic the constructor's `return` behavior
        // http://es5.github.io/#x13.2.2
        if (this instanceof bound) {
          // ensure `new bound` is an instance of `func`
          var thisBinding = baseCreate(func.prototype),
              result = func.apply(thisBinding, args || arguments);
          return isObject(result) ? result : thisBinding;
        }
        return func.apply(thisArg, args || arguments);
      }
      setBindData(bound, bindData);
      return bound;
    }

    /**
     * The base implementation of `_.clone` without argument juggling or support
     * for `thisArg` binding.
     *
     * @private
     * @param {*} value The value to clone.
     * @param {boolean} [isDeep=false] Specify a deep clone.
     * @param {Function} [callback] The function to customize cloning values.
     * @param {Array} [stackA=[]] Tracks traversed source objects.
     * @param {Array} [stackB=[]] Associates clones with source counterparts.
     * @returns {*} Returns the cloned value.
     */
    function baseClone(value, isDeep, callback, stackA, stackB) {
      if (callback) {
        var result = callback(value);
        if (typeof result != 'undefined') {
          return result;
        }
      }
      // inspect [[Class]]
      var isObj = isObject(value);
      if (isObj) {
        var className = toString.call(value);
        if (!cloneableClasses[className]) {
          return value;
        }
        var ctor = ctorByClass[className];
        switch (className) {
          case boolClass:
          case dateClass:
            return new ctor(+value);

          case numberClass:
          case stringClass:
            return new ctor(value);

          case regexpClass:
            result = ctor(value.source, reFlags.exec(value));
            result.lastIndex = value.lastIndex;
            return result;
        }
      } else {
        return value;
      }
      var isArr = isArray(value);
      if (isDeep) {
        // check for circular references and return corresponding clone
        var initedStack = !stackA;
        stackA || (stackA = getArray());
        stackB || (stackB = getArray());

        var length = stackA.length;
        while (length--) {
          if (stackA[length] == value) {
            return stackB[length];
          }
        }
        result = isArr ? ctor(value.length) : {};
      }
      else {
        result = isArr ? slice(value) : assign({}, value);
      }
      // add array properties assigned by `RegExp#exec`
      if (isArr) {
        if (hasOwnProperty.call(value, 'index')) {
          result.index = value.index;
        }
        if (hasOwnProperty.call(value, 'input')) {
          result.input = value.input;
        }
      }
      // exit for shallow clone
      if (!isDeep) {
        return result;
      }
      // add the source value to the stack of traversed objects
      // and associate it with its clone
      stackA.push(value);
      stackB.push(result);

      // recursively populate clone (susceptible to call stack limits)
      (isArr ? forEach : forOwn)(value, function(objValue, key) {
        result[key] = baseClone(objValue, isDeep, callback, stackA, stackB);
      });

      if (initedStack) {
        releaseArray(stackA);
        releaseArray(stackB);
      }
      return result;
    }

    /**
     * The base implementation of `_.create` without support for assigning
     * properties to the created object.
     *
     * @private
     * @param {Object} prototype The object to inherit from.
     * @returns {Object} Returns the new object.
     */
    function baseCreate(prototype, properties) {
      return isObject(prototype) ? nativeCreate(prototype) : {};
    }
    // fallback for browsers without `Object.create`
    if (!nativeCreate) {
      baseCreate = (function() {
        function Object() {}
        return function(prototype) {
          if (isObject(prototype)) {
            Object.prototype = prototype;
            var result = new Object;
            Object.prototype = null;
          }
          return result || context.Object();
        };
      }());
    }

    /**
     * The base implementation of `_.createCallback` without support for creating
     * "_.pluck" or "_.where" style callbacks.
     *
     * @private
     * @param {*} [func=identity] The value to convert to a callback.
     * @param {*} [thisArg] The `this` binding of the created callback.
     * @param {number} [argCount] The number of arguments the callback accepts.
     * @returns {Function} Returns a callback function.
     */
    function baseCreateCallback(func, thisArg, argCount) {
      if (typeof func != 'function') {
        return identity;
      }
      // exit early for no `thisArg` or already bound by `Function#bind`
      if (typeof thisArg == 'undefined' || !('prototype' in func)) {
        return func;
      }
      var bindData = func.__bindData__;
      if (typeof bindData == 'undefined') {
        if (support.funcNames) {
          bindData = !func.name;
        }
        bindData = bindData || !support.funcDecomp;
        if (!bindData) {
          var source = fnToString.call(func);
          if (!support.funcNames) {
            bindData = !reFuncName.test(source);
          }
          if (!bindData) {
            // checks if `func` references the `this` keyword and stores the result
            bindData = reThis.test(source);
            setBindData(func, bindData);
          }
        }
      }
      // exit early if there are no `this` references or `func` is bound
      if (bindData === false || (bindData !== true && bindData[1] & 1)) {
        return func;
      }
      switch (argCount) {
        case 1: return function(value) {
          return func.call(thisArg, value);
        };
        case 2: return function(a, b) {
          return func.call(thisArg, a, b);
        };
        case 3: return function(value, index, collection) {
          return func.call(thisArg, value, index, collection);
        };
        case 4: return function(accumulator, value, index, collection) {
          return func.call(thisArg, accumulator, value, index, collection);
        };
      }
      return bind(func, thisArg);
    }

    /**
     * The base implementation of `createWrapper` that creates the wrapper and
     * sets its meta data.
     *
     * @private
     * @param {Array} bindData The bind data array.
     * @returns {Function} Returns the new function.
     */
    function baseCreateWrapper(bindData) {
      var func = bindData[0],
          bitmask = bindData[1],
          partialArgs = bindData[2],
          partialRightArgs = bindData[3],
          thisArg = bindData[4],
          arity = bindData[5];

      var isBind = bitmask & 1,
          isBindKey = bitmask & 2,
          isCurry = bitmask & 4,
          isCurryBound = bitmask & 8,
          key = func;

      function bound() {
        var thisBinding = isBind ? thisArg : this;
        if (partialArgs) {
          var args = slice(partialArgs);
          push.apply(args, arguments);
        }
        if (partialRightArgs || isCurry) {
          args || (args = slice(arguments));
          if (partialRightArgs) {
            push.apply(args, partialRightArgs);
          }
          if (isCurry && args.length < arity) {
            bitmask |= 16 & ~32;
            return baseCreateWrapper([func, (isCurryBound ? bitmask : bitmask & ~3), args, null, thisArg, arity]);
          }
        }
        args || (args = arguments);
        if (isBindKey) {
          func = thisBinding[key];
        }
        if (this instanceof bound) {
          thisBinding = baseCreate(func.prototype);
          var result = func.apply(thisBinding, args);
          return isObject(result) ? result : thisBinding;
        }
        return func.apply(thisBinding, args);
      }
      setBindData(bound, bindData);
      return bound;
    }

    /**
     * The base implementation of `_.difference` that accepts a single array
     * of values to exclude.
     *
     * @private
     * @param {Array} array The array to process.
     * @param {Array} [values] The array of values to exclude.
     * @returns {Array} Returns a new array of filtered values.
     */
    function baseDifference(array, values) {
      var index = -1,
          indexOf = getIndexOf(),
          length = array ? array.length : 0,
          isLarge = length >= largeArraySize && indexOf === baseIndexOf,
          result = [];

      if (isLarge) {
        var cache = createCache(values);
        if (cache) {
          indexOf = cacheIndexOf;
          values = cache;
        } else {
          isLarge = false;
        }
      }
      while (++index < length) {
        var value = array[index];
        if (indexOf(values, value) < 0) {
          result.push(value);
        }
      }
      if (isLarge) {
        releaseObject(values);
      }
      return result;
    }

    /**
     * The base implementation of `_.flatten` without support for callback
     * shorthands or `thisArg` binding.
     *
     * @private
     * @param {Array} array The array to flatten.
     * @param {boolean} [isShallow=false] A flag to restrict flattening to a single level.
     * @param {boolean} [isStrict=false] A flag to restrict flattening to arrays and `arguments` objects.
     * @param {number} [fromIndex=0] The index to start from.
     * @returns {Array} Returns a new flattened array.
     */
    function baseFlatten(array, isShallow, isStrict, fromIndex) {
      var index = (fromIndex || 0) - 1,
          length = array ? array.length : 0,
          result = [];

      while (++index < length) {
        var value = array[index];

        if (value && typeof value == 'object' && typeof value.length == 'number'
            && (isArray(value) || isArguments(value))) {
          // recursively flatten arrays (susceptible to call stack limits)
          if (!isShallow) {
            value = baseFlatten(value, isShallow, isStrict);
          }
          var valIndex = -1,
              valLength = value.length,
              resIndex = result.length;

          result.length += valLength;
          while (++valIndex < valLength) {
            result[resIndex++] = value[valIndex];
          }
        } else if (!isStrict) {
          result.push(value);
        }
      }
      return result;
    }

    /**
     * The base implementation of `_.isEqual`, without support for `thisArg` binding,
     * that allows partial "_.where" style comparisons.
     *
     * @private
     * @param {*} a The value to compare.
     * @param {*} b The other value to compare.
     * @param {Function} [callback] The function to customize comparing values.
     * @param {Function} [isWhere=false] A flag to indicate performing partial comparisons.
     * @param {Array} [stackA=[]] Tracks traversed `a` objects.
     * @param {Array} [stackB=[]] Tracks traversed `b` objects.
     * @returns {boolean} Returns `true` if the values are equivalent, else `false`.
     */
    function baseIsEqual(a, b, callback, isWhere, stackA, stackB) {
      // used to indicate that when comparing objects, `a` has at least the properties of `b`
      if (callback) {
        var result = callback(a, b);
        if (typeof result != 'undefined') {
          return !!result;
        }
      }
      // exit early for identical values
      if (a === b) {
        // treat `+0` vs. `-0` as not equal
        return a !== 0 || (1 / a == 1 / b);
      }
      var type = typeof a,
          otherType = typeof b;

      // exit early for unlike primitive values
      if (a === a &&
          !(a && objectTypes[type]) &&
          !(b && objectTypes[otherType])) {
        return false;
      }
      // exit early for `null` and `undefined` avoiding ES3's Function#call behavior
      // http://es5.github.io/#x15.3.4.4
      if (a == null || b == null) {
        return a === b;
      }
      // compare [[Class]] names
      var className = toString.call(a),
          otherClass = toString.call(b);

      if (className == argsClass) {
        className = objectClass;
      }
      if (otherClass == argsClass) {
        otherClass = objectClass;
      }
      if (className != otherClass) {
        return false;
      }
      switch (className) {
        case boolClass:
        case dateClass:
          // coerce dates and booleans to numbers, dates to milliseconds and booleans
          // to `1` or `0` treating invalid dates coerced to `NaN` as not equal
          return +a == +b;

        case numberClass:
          // treat `NaN` vs. `NaN` as equal
          return (a != +a)
            ? b != +b
            // but treat `+0` vs. `-0` as not equal
            : (a == 0 ? (1 / a == 1 / b) : a == +b);

        case regexpClass:
        case stringClass:
          // coerce regexes to strings (http://es5.github.io/#x15.10.6.4)
          // treat string primitives and their corresponding object instances as equal
          return a == String(b);
      }
      var isArr = className == arrayClass;
      if (!isArr) {
        // unwrap any `lodash` wrapped values
        var aWrapped = hasOwnProperty.call(a, '__wrapped__'),
            bWrapped = hasOwnProperty.call(b, '__wrapped__');

        if (aWrapped || bWrapped) {
          return baseIsEqual(aWrapped ? a.__wrapped__ : a, bWrapped ? b.__wrapped__ : b, callback, isWhere, stackA, stackB);
        }
        // exit for functions and DOM nodes
        if (className != objectClass) {
          return false;
        }
        // in older versions of Opera, `arguments` objects have `Array` constructors
        var ctorA = a.constructor,
            ctorB = b.constructor;

        // non `Object` object instances with different constructors are not equal
        if (ctorA != ctorB &&
              !(isFunction(ctorA) && ctorA instanceof ctorA && isFunction(ctorB) && ctorB instanceof ctorB) &&
              ('constructor' in a && 'constructor' in b)
            ) {
          return false;
        }
      }
      // assume cyclic structures are equal
      // the algorithm for detecting cyclic structures is adapted from ES 5.1
      // section 15.12.3, abstract operation `JO` (http://es5.github.io/#x15.12.3)
      var initedStack = !stackA;
      stackA || (stackA = getArray());
      stackB || (stackB = getArray());

      var length = stackA.length;
      while (length--) {
        if (stackA[length] == a) {
          return stackB[length] == b;
        }
      }
      var size = 0;
      result = true;

      // add `a` and `b` to the stack of traversed objects
      stackA.push(a);
      stackB.push(b);

      // recursively compare objects and arrays (susceptible to call stack limits)
      if (isArr) {
        // compare lengths to determine if a deep comparison is necessary
        length = a.length;
        size = b.length;
        result = size == length;

        if (result || isWhere) {
          // deep compare the contents, ignoring non-numeric properties
          while (size--) {
            var index = length,
                value = b[size];

            if (isWhere) {
              while (index--) {
                if ((result = baseIsEqual(a[index], value, callback, isWhere, stackA, stackB))) {
                  break;
                }
              }
            } else if (!(result = baseIsEqual(a[size], value, callback, isWhere, stackA, stackB))) {
              break;
            }
          }
        }
      }
      else {
        // deep compare objects using `forIn`, instead of `forOwn`, to avoid `Object.keys`
        // which, in this case, is more costly
        forIn(b, function(value, key, b) {
          if (hasOwnProperty.call(b, key)) {
            // count the number of properties.
            size++;
            // deep compare each property value.
            return (result = hasOwnProperty.call(a, key) && baseIsEqual(a[key], value, callback, isWhere, stackA, stackB));
          }
        });

        if (result && !isWhere) {
          // ensure both objects have the same number of properties
          forIn(a, function(value, key, a) {
            if (hasOwnProperty.call(a, key)) {
              // `size` will be `-1` if `a` has more properties than `b`
              return (result = --size > -1);
            }
          });
        }
      }
      stackA.pop();
      stackB.pop();

      if (initedStack) {
        releaseArray(stackA);
        releaseArray(stackB);
      }
      return result;
    }

    /**
     * The base implementation of `_.merge` without argument juggling or support
     * for `thisArg` binding.
     *
     * @private
     * @param {Object} object The destination object.
     * @param {Object} source The source object.
     * @param {Function} [callback] The function to customize merging properties.
     * @param {Array} [stackA=[]] Tracks traversed source objects.
     * @param {Array} [stackB=[]] Associates values with source counterparts.
     */
    function baseMerge(object, source, callback, stackA, stackB) {
      (isArray(source) ? forEach : forOwn)(source, function(source, key) {
        var found,
            isArr,
            result = source,
            value = object[key];

        if (source && ((isArr = isArray(source)) || isPlainObject(source))) {
          // avoid merging previously merged cyclic sources
          var stackLength = stackA.length;
          while (stackLength--) {
            if ((found = stackA[stackLength] == source)) {
              value = stackB[stackLength];
              break;
            }
          }
          if (!found) {
            var isShallow;
            if (callback) {
              result = callback(value, source);
              if ((isShallow = typeof result != 'undefined')) {
                value = result;
              }
            }
            if (!isShallow) {
              value = isArr
                ? (isArray(value) ? value : [])
                : (isPlainObject(value) ? value : {});
            }
            // add `source` and associated `value` to the stack of traversed objects
            stackA.push(source);
            stackB.push(value);

            // recursively merge objects and arrays (susceptible to call stack limits)
            if (!isShallow) {
              baseMerge(value, source, callback, stackA, stackB);
            }
          }
        }
        else {
          if (callback) {
            result = callback(value, source);
            if (typeof result == 'undefined') {
              result = source;
            }
          }
          if (typeof result != 'undefined') {
            value = result;
          }
        }
        object[key] = value;
      });
    }

    /**
     * The base implementation of `_.random` without argument juggling or support
     * for returning floating-point numbers.
     *
     * @private
     * @param {number} min The minimum possible value.
     * @param {number} max The maximum possible value.
     * @returns {number} Returns a random number.
     */
    function baseRandom(min, max) {
      return min + floor(nativeRandom() * (max - min + 1));
    }

    /**
     * The base implementation of `_.uniq` without support for callback shorthands
     * or `thisArg` binding.
     *
     * @private
     * @param {Array} array The array to process.
     * @param {boolean} [isSorted=false] A flag to indicate that `array` is sorted.
     * @param {Function} [callback] The function called per iteration.
     * @returns {Array} Returns a duplicate-value-free array.
     */
    function baseUniq(array, isSorted, callback) {
      var index = -1,
          indexOf = getIndexOf(),
          length = array ? array.length : 0,
          result = [];

      var isLarge = !isSorted && length >= largeArraySize && indexOf === baseIndexOf,
          seen = (callback || isLarge) ? getArray() : result;

      if (isLarge) {
        var cache = createCache(seen);
        indexOf = cacheIndexOf;
        seen = cache;
      }
      while (++index < length) {
        var value = array[index],
            computed = callback ? callback(value, index, array) : value;

        if (isSorted
              ? !index || seen[seen.length - 1] !== computed
              : indexOf(seen, computed) < 0
            ) {
          if (callback || isLarge) {
            seen.push(computed);
          }
          result.push(value);
        }
      }
      if (isLarge) {
        releaseArray(seen.array);
        releaseObject(seen);
      } else if (callback) {
        releaseArray(seen);
      }
      return result;
    }

    /**
     * Creates a function that aggregates a collection, creating an object composed
     * of keys generated from the results of running each element of the collection
     * through a callback. The given `setter` function sets the keys and values
     * of the composed object.
     *
     * @private
     * @param {Function} setter The setter function.
     * @returns {Function} Returns the new aggregator function.
     */
    function createAggregator(setter) {
      return function(collection, callback, thisArg) {
        var result = {};
        callback = lodash.createCallback(callback, thisArg, 3);

        var index = -1,
            length = collection ? collection.length : 0;

        if (typeof length == 'number') {
          while (++index < length) {
            var value = collection[index];
            setter(result, value, callback(value, index, collection), collection);
          }
        } else {
          forOwn(collection, function(value, key, collection) {
            setter(result, value, callback(value, key, collection), collection);
          });
        }
        return result;
      };
    }

    /**
     * Creates a function that, when called, either curries or invokes `func`
     * with an optional `this` binding and partially applied arguments.
     *
     * @private
     * @param {Function|string} func The function or method name to reference.
     * @param {number} bitmask The bitmask of method flags to compose.
     *  The bitmask may be composed of the following flags:
     *  1 - `_.bind`
     *  2 - `_.bindKey`
     *  4 - `_.curry`
     *  8 - `_.curry` (bound)
     *  16 - `_.partial`
     *  32 - `_.partialRight`
     * @param {Array} [partialArgs] An array of arguments to prepend to those
     *  provided to the new function.
     * @param {Array} [partialRightArgs] An array of arguments to append to those
     *  provided to the new function.
     * @param {*} [thisArg] The `this` binding of `func`.
     * @param {number} [arity] The arity of `func`.
     * @returns {Function} Returns the new function.
     */
    function createWrapper(func, bitmask, partialArgs, partialRightArgs, thisArg, arity) {
      var isBind = bitmask & 1,
          isBindKey = bitmask & 2,
          isCurry = bitmask & 4,
          isCurryBound = bitmask & 8,
          isPartial = bitmask & 16,
          isPartialRight = bitmask & 32;

      if (!isBindKey && !isFunction(func)) {
        throw new TypeError;
      }
      if (isPartial && !partialArgs.length) {
        bitmask &= ~16;
        isPartial = partialArgs = false;
      }
      if (isPartialRight && !partialRightArgs.length) {
        bitmask &= ~32;
        isPartialRight = partialRightArgs = false;
      }
      var bindData = func && func.__bindData__;
      if (bindData && bindData !== true) {
        // clone `bindData`
        bindData = slice(bindData);
        if (bindData[2]) {
          bindData[2] = slice(bindData[2]);
        }
        if (bindData[3]) {
          bindData[3] = slice(bindData[3]);
        }
        // set `thisBinding` is not previously bound
        if (isBind && !(bindData[1] & 1)) {
          bindData[4] = thisArg;
        }
        // set if previously bound but not currently (subsequent curried functions)
        if (!isBind && bindData[1] & 1) {
          bitmask |= 8;
        }
        // set curried arity if not yet set
        if (isCurry && !(bindData[1] & 4)) {
          bindData[5] = arity;
        }
        // append partial left arguments
        if (isPartial) {
          push.apply(bindData[2] || (bindData[2] = []), partialArgs);
        }
        // append partial right arguments
        if (isPartialRight) {
          unshift.apply(bindData[3] || (bindData[3] = []), partialRightArgs);
        }
        // merge flags
        bindData[1] |= bitmask;
        return createWrapper.apply(null, bindData);
      }
      // fast path for `_.bind`
      var creater = (bitmask == 1 || bitmask === 17) ? baseBind : baseCreateWrapper;
      return creater([func, bitmask, partialArgs, partialRightArgs, thisArg, arity]);
    }

    /**
     * Used by `escape` to convert characters to HTML entities.
     *
     * @private
     * @param {string} match The matched character to escape.
     * @returns {string} Returns the escaped character.
     */
    function escapeHtmlChar(match) {
      return htmlEscapes[match];
    }

    /**
     * Gets the appropriate "indexOf" function. If the `_.indexOf` method is
     * customized, this method returns the custom method, otherwise it returns
     * the `baseIndexOf` function.
     *
     * @private
     * @returns {Function} Returns the "indexOf" function.
     */
    function getIndexOf() {
      var result = (result = lodash.indexOf) === indexOf ? baseIndexOf : result;
      return result;
    }

    /**
     * Checks if `value` is a native function.
     *
     * @private
     * @param {*} value The value to check.
     * @returns {boolean} Returns `true` if the `value` is a native function, else `false`.
     */
    function isNative(value) {
      return typeof value == 'function' && reNative.test(value);
    }

    /**
     * Sets `this` binding data on a given function.
     *
     * @private
     * @param {Function} func The function to set data on.
     * @param {Array} value The data array to set.
     */
    var setBindData = !defineProperty ? noop : function(func, value) {
      descriptor.value = value;
      defineProperty(func, '__bindData__', descriptor);
      descriptor.value = null;
    };

    /**
     * A fallback implementation of `isPlainObject` which checks if a given value
     * is an object created by the `Object` constructor, assuming objects created
     * by the `Object` constructor have no inherited enumerable properties and that
     * there are no `Object.prototype` extensions.
     *
     * @private
     * @param {*} value The value to check.
     * @returns {boolean} Returns `true` if `value` is a plain object, else `false`.
     */
    function shimIsPlainObject(value) {
      var ctor,
          result;

      // avoid non Object objects, `arguments` objects, and DOM elements
      if (!(value && toString.call(value) == objectClass) ||
          (ctor = value.constructor, isFunction(ctor) && !(ctor instanceof ctor))) {
        return false;
      }
      // In most environments an object's own properties are iterated before
      // its inherited properties. If the last iterated property is an object's
      // own property then there are no inherited enumerable properties.
      forIn(value, function(value, key) {
        result = key;
      });
      return typeof result == 'undefined' || hasOwnProperty.call(value, result);
    }

    /**
     * Used by `unescape` to convert HTML entities to characters.
     *
     * @private
     * @param {string} match The matched character to unescape.
     * @returns {string} Returns the unescaped character.
     */
    function unescapeHtmlChar(match) {
      return htmlUnescapes[match];
    }

    /*--------------------------------------------------------------------------*/

    /**
     * Checks if `value` is an `arguments` object.
     *
     * @static
     * @memberOf _
     * @category Objects
     * @param {*} value The value to check.
     * @returns {boolean} Returns `true` if the `value` is an `arguments` object, else `false`.
     * @example
     *
     * (function() { return _.isArguments(arguments); })(1, 2, 3);
     * // => true
     *
     * _.isArguments([1, 2, 3]);
     * // => false
     */
    function isArguments(value) {
      return value && typeof value == 'object' && typeof value.length == 'number' &&
        toString.call(value) == argsClass || false;
    }

    /**
     * Checks if `value` is an array.
     *
     * @static
     * @memberOf _
     * @type Function
     * @category Objects
     * @param {*} value The value to check.
     * @returns {boolean} Returns `true` if the `value` is an array, else `false`.
     * @example
     *
     * (function() { return _.isArray(arguments); })();
     * // => false
     *
     * _.isArray([1, 2, 3]);
     * // => true
     */
    var isArray = nativeIsArray || function(value) {
      return value && typeof value == 'object' && typeof value.length == 'number' &&
        toString.call(value) == arrayClass || false;
    };

    /**
     * A fallback implementation of `Object.keys` which produces an array of the
     * given object's own enumerable property names.
     *
     * @private
     * @type Function
     * @param {Object} object The object to inspect.
     * @returns {Array} Returns an array of property names.
     */
    var shimKeys = function(object) {
      var index, iterable = object, result = [];
      if (!iterable) return result;
      if (!(objectTypes[typeof object])) return result;
        for (index in iterable) {
          if (hasOwnProperty.call(iterable, index)) {
            result.push(index);
          }
        }
      return result
    };

    /**
     * Creates an array composed of the own enumerable property names of an object.
     *
     * @static
     * @memberOf _
     * @category Objects
     * @param {Object} object The object to inspect.
     * @returns {Array} Returns an array of property names.
     * @example
     *
     * _.keys({ 'one': 1, 'two': 2, 'three': 3 });
     * // => ['one', 'two', 'three'] (property order is not guaranteed across environments)
     */
    var keys = !nativeKeys ? shimKeys : function(object) {
      if (!isObject(object)) {
        return [];
      }
      return nativeKeys(object);
    };

    /**
     * Used to convert characters to HTML entities:
     *
     * Though the `>` character is escaped for symmetry, characters like `>` and `/`
     * don't require escaping in HTML and have no special meaning unless they're part
     * of a tag or an unquoted attribute value.
     * http://mathiasbynens.be/notes/ambiguous-ampersands (under "semi-related fun fact")
     */
    var htmlEscapes = {
      '&': '&amp;',
      '<': '&lt;',
      '>': '&gt;',
      '"': '&quot;',
      "'": '&#39;'
    };

    /** Used to convert HTML entities to characters */
    var htmlUnescapes = invert(htmlEscapes);

    /** Used to match HTML entities and HTML characters */
    var reEscapedHtml = RegExp('(' + keys(htmlUnescapes).join('|') + ')', 'g'),
        reUnescapedHtml = RegExp('[' + keys(htmlEscapes).join('') + ']', 'g');

    /*--------------------------------------------------------------------------*/

    /**
     * Assigns own enumerable properties of source object(s) to the destination
     * object. Subsequent sources will overwrite property assignments of previous
     * sources. If a callback is provided it will be executed to produce the
     * assigned values. The callback is bound to `thisArg` and invoked with two
     * arguments; (objectValue, sourceValue).
     *
     * @static
     * @memberOf _
     * @type Function
     * @alias extend
     * @category Objects
     * @param {Object} object The destination object.
     * @param {...Object} [source] The source objects.
     * @param {Function} [callback] The function to customize assigning values.
     * @param {*} [thisArg] The `this` binding of `callback`.
     * @returns {Object} Returns the destination object.
     * @example
     *
     * _.assign({ 'name': 'fred' }, { 'employer': 'slate' });
     * // => { 'name': 'fred', 'employer': 'slate' }
     *
     * var defaults = _.partialRight(_.assign, function(a, b) {
     *   return typeof a == 'undefined' ? b : a;
     * });
     *
     * var object = { 'name': 'barney' };
     * defaults(object, { 'name': 'fred', 'employer': 'slate' });
     * // => { 'name': 'barney', 'employer': 'slate' }
     */
    var assign = function(object, source, guard) {
      var index, iterable = object, result = iterable;
      if (!iterable) return result;
      var args = arguments,
          argsIndex = 0,
          argsLength = typeof guard == 'number' ? 2 : args.length;
      if (argsLength > 3 && typeof args[argsLength - 2] == 'function') {
        var callback = baseCreateCallback(args[--argsLength - 1], args[argsLength--], 2);
      } else if (argsLength > 2 && typeof args[argsLength - 1] == 'function') {
        callback = args[--argsLength];
      }
      while (++argsIndex < argsLength) {
        iterable = args[argsIndex];
        if (iterable && objectTypes[typeof iterable]) {
        var ownIndex = -1,
            ownProps = objectTypes[typeof iterable] && keys(iterable),
            length = ownProps ? ownProps.length : 0;

        while (++ownIndex < length) {
          index = ownProps[ownIndex];
          result[index] = callback ? callback(result[index], iterable[index]) : iterable[index];
        }
        }
      }
      return result
    };

    /**
     * Creates a clone of `value`. If `isDeep` is `true` nested objects will also
     * be cloned, otherwise they will be assigned by reference. If a callback
     * is provided it will be executed to produce the cloned values. If the
     * callback returns `undefined` cloning will be handled by the method instead.
     * The callback is bound to `thisArg` and invoked with one argument; (value).
     *
     * @static
     * @memberOf _
     * @category Objects
     * @param {*} value The value to clone.
     * @param {boolean} [isDeep=false] Specify a deep clone.
     * @param {Function} [callback] The function to customize cloning values.
     * @param {*} [thisArg] The `this` binding of `callback`.
     * @returns {*} Returns the cloned value.
     * @example
     *
     * var characters = [
     *   { 'name': 'barney', 'age': 36 },
     *   { 'name': 'fred',   'age': 40 }
     * ];
     *
     * var shallow = _.clone(characters);
     * shallow[0] === characters[0];
     * // => true
     *
     * var deep = _.clone(characters, true);
     * deep[0] === characters[0];
     * // => false
     *
     * _.mixin({
     *   'clone': _.partialRight(_.clone, function(value) {
     *     return _.isElement(value) ? value.cloneNode(false) : undefined;
     *   })
     * });
     *
     * var clone = _.clone(document.body);
     * clone.childNodes.length;
     * // => 0
     */
    function clone(value, isDeep, callback, thisArg) {
      // allows working with "Collections" methods without using their `index`
      // and `collection` arguments for `isDeep` and `callback`
      if (typeof isDeep != 'boolean' && isDeep != null) {
        thisArg = callback;
        callback = isDeep;
        isDeep = false;
      }
      return baseClone(value, isDeep, typeof callback == 'function' && baseCreateCallback(callback, thisArg, 1));
    }

    /**
     * Creates a deep clone of `value`. If a callback is provided it will be
     * executed to produce the cloned values. If the callback returns `undefined`
     * cloning will be handled by the method instead. The callback is bound to
     * `thisArg` and invoked with one argument; (value).
     *
     * Note: This method is loosely based on the structured clone algorithm. Functions
     * and DOM nodes are **not** cloned. The enumerable properties of `arguments` objects and
     * objects created by constructors other than `Object` are cloned to plain `Object` objects.
     * See http://www.w3.org/TR/html5/infrastructure.html#internal-structured-cloning-algorithm.
     *
     * @static
     * @memberOf _
     * @category Objects
     * @param {*} value The value to deep clone.
     * @param {Function} [callback] The function to customize cloning values.
     * @param {*} [thisArg] The `this` binding of `callback`.
     * @returns {*} Returns the deep cloned value.
     * @example
     *
     * var characters = [
     *   { 'name': 'barney', 'age': 36 },
     *   { 'name': 'fred',   'age': 40 }
     * ];
     *
     * var deep = _.cloneDeep(characters);
     * deep[0] === characters[0];
     * // => false
     *
     * var view = {
     *   'label': 'docs',
     *   'node': element
     * };
     *
     * var clone = _.cloneDeep(view, function(value) {
     *   return _.isElement(value) ? value.cloneNode(true) : undefined;
     * });
     *
     * clone.node == view.node;
     * // => false
     */
    function cloneDeep(value, callback, thisArg) {
      return baseClone(value, true, typeof callback == 'function' && baseCreateCallback(callback, thisArg, 1));
    }

    /**
     * Creates an object that inherits from the given `prototype` object. If a
     * `properties` object is provided its own enumerable properties are assigned
     * to the created object.
     *
     * @static
     * @memberOf _
     * @category Objects
     * @param {Object} prototype The object to inherit from.
     * @param {Object} [properties] The properties to assign to the object.
     * @returns {Object} Returns the new object.
     * @example
     *
     * function Shape() {
     *   this.x = 0;
     *   this.y = 0;
     * }
     *
     * function Circle() {
     *   Shape.call(this);
     * }
     *
     * Circle.prototype = _.create(Shape.prototype, { 'constructor': Circle });
     *
     * var circle = new Circle;
     * circle instanceof Circle;
     * // => true
     *
     * circle instanceof Shape;
     * // => true
     */
    function create(prototype, properties) {
      var result = baseCreate(prototype);
      return properties ? assign(result, properties) : result;
    }

    /**
     * Assigns own enumerable properties of source object(s) to the destination
     * object for all destination properties that resolve to `undefined`. Once a
     * property is set, additional defaults of the same property will be ignored.
     *
     * @static
     * @memberOf _
     * @type Function
     * @category Objects
     * @param {Object} object The destination object.
     * @param {...Object} [source] The source objects.
     * @param- {Object} [guard] Allows working with `_.reduce` without using its
     *  `key` and `object` arguments as sources.
     * @returns {Object} Returns the destination object.
     * @example
     *
     * var object = { 'name': 'barney' };
     * _.defaults(object, { 'name': 'fred', 'employer': 'slate' });
     * // => { 'name': 'barney', 'employer': 'slate' }
     */
    var defaults = function(object, source, guard) {
      var index, iterable = object, result = iterable;
      if (!iterable) return result;
      var args = arguments,
          argsIndex = 0,
          argsLength = typeof guard == 'number' ? 2 : args.length;
      while (++argsIndex < argsLength) {
        iterable = args[argsIndex];
        if (iterable && objectTypes[typeof iterable]) {
        var ownIndex = -1,
            ownProps = objectTypes[typeof iterable] && keys(iterable),
            length = ownProps ? ownProps.length : 0;

        while (++ownIndex < length) {
          index = ownProps[ownIndex];
          if (typeof result[index] == 'undefined') result[index] = iterable[index];
        }
        }
      }
      return result
    };

    /**
     * This method is like `_.findIndex` except that it returns the key of the
     * first element that passes the callback check, instead of the element itself.
     *
     * If a property name is provided for `callback` the created "_.pluck" style
     * callback will return the property value of the given element.
     *
     * If an object is provided for `callback` the created "_.where" style callback
     * will return `true` for elements that have the properties of the given object,
     * else `false`.
     *
     * @static
     * @memberOf _
     * @category Objects
     * @param {Object} object The object to search.
     * @param {Function|Object|string} [callback=identity] The function called per
     *  iteration. If a property name or object is provided it will be used to
     *  create a "_.pluck" or "_.where" style callback, respectively.
     * @param {*} [thisArg] The `this` binding of `callback`.
     * @returns {string|undefined} Returns the key of the found element, else `undefined`.
     * @example
     *
     * var characters = {
     *   'barney': {  'age': 36, 'blocked': false },
     *   'fred': {    'age': 40, 'blocked': true },
     *   'pebbles': { 'age': 1,  'blocked': false }
     * };
     *
     * _.findKey(characters, function(chr) {
     *   return chr.age < 40;
     * });
     * // => 'barney' (property order is not guaranteed across environments)
     *
     * // using "_.where" callback shorthand
     * _.findKey(characters, { 'age': 1 });
     * // => 'pebbles'
     *
     * // using "_.pluck" callback shorthand
     * _.findKey(characters, 'blocked');
     * // => 'fred'
     */
    function findKey(object, callback, thisArg) {
      var result;
      callback = lodash.createCallback(callback, thisArg, 3);
      forOwn(object, function(value, key, object) {
        if (callback(value, key, object)) {
          result = key;
          return false;
        }
      });
      return result;
    }

    /**
     * This method is like `_.findKey` except that it iterates over elements
     * of a `collection` in the opposite order.
     *
     * If a property name is provided for `callback` the created "_.pluck" style
     * callback will return the property value of the given element.
     *
     * If an object is provided for `callback` the created "_.where" style callback
     * will return `true` for elements that have the properties of the given object,
     * else `false`.
     *
     * @static
     * @memberOf _
     * @category Objects
     * @param {Object} object The object to search.
     * @param {Function|Object|string} [callback=identity] The function called per
     *  iteration. If a property name or object is provided it will be used to
     *  create a "_.pluck" or "_.where" style callback, respectively.
     * @param {*} [thisArg] The `this` binding of `callback`.
     * @returns {string|undefined} Returns the key of the found element, else `undefined`.
     * @example
     *
     * var characters = {
     *   'barney': {  'age': 36, 'blocked': true },
     *   'fred': {    'age': 40, 'blocked': false },
     *   'pebbles': { 'age': 1,  'blocked': true }
     * };
     *
     * _.findLastKey(characters, function(chr) {
     *   return chr.age < 40;
     * });
     * // => returns `pebbles`, assuming `_.findKey` returns `barney`
     *
     * // using "_.where" callback shorthand
     * _.findLastKey(characters, { 'age': 40 });
     * // => 'fred'
     *
     * // using "_.pluck" callback shorthand
     * _.findLastKey(characters, 'blocked');
     * // => 'pebbles'
     */
    function findLastKey(object, callback, thisArg) {
      var result;
      callback = lodash.createCallback(callback, thisArg, 3);
      forOwnRight(object, function(value, key, object) {
        if (callback(value, key, object)) {
          result = key;
          return false;
        }
      });
      return result;
    }

    /**
     * Iterates over own and inherited enumerable properties of an object,
     * executing the callback for each property. The callback is bound to `thisArg`
     * and invoked with three arguments; (value, key, object). Callbacks may exit
     * iteration early by explicitly returning `false`.
     *
     * @static
     * @memberOf _
     * @type Function
     * @category Objects
     * @param {Object} object The object to iterate over.
     * @param {Function} [callback=identity] The function called per iteration.
     * @param {*} [thisArg] The `this` binding of `callback`.
     * @returns {Object} Returns `object`.
     * @example
     *
     * function Shape() {
     *   this.x = 0;
     *   this.y = 0;
     * }
     *
     * Shape.prototype.move = function(x, y) {
     *   this.x += x;
     *   this.y += y;
     * };
     *
     * _.forIn(new Shape, function(value, key) {
     *   console.log(key);
     * });
     * // => logs 'x', 'y', and 'move' (property order is not guaranteed across environments)
     */
    var forIn = function(collection, callback, thisArg) {
      var index, iterable = collection, result = iterable;
      if (!iterable) return result;
      if (!objectTypes[typeof iterable]) return result;
      callback = callback && typeof thisArg == 'undefined' ? callback : baseCreateCallback(callback, thisArg, 3);
        for (index in iterable) {
          if (callback(iterable[index], index, collection) === false) return result;
        }
      return result
    };

    /**
     * This method is like `_.forIn` except that it iterates over elements
     * of a `collection` in the opposite order.
     *
     * @static
     * @memberOf _
     * @category Objects
     * @param {Object} object The object to iterate over.
     * @param {Function} [callback=identity] The function called per iteration.
     * @param {*} [thisArg] The `this` binding of `callback`.
     * @returns {Object} Returns `object`.
     * @example
     *
     * function Shape() {
     *   this.x = 0;
     *   this.y = 0;
     * }
     *
     * Shape.prototype.move = function(x, y) {
     *   this.x += x;
     *   this.y += y;
     * };
     *
     * _.forInRight(new Shape, function(value, key) {
     *   console.log(key);
     * });
     * // => logs 'move', 'y', and 'x' assuming `_.forIn ` logs 'x', 'y', and 'move'
     */
    function forInRight(object, callback, thisArg) {
      var pairs = [];

      forIn(object, function(value, key) {
        pairs.push(key, value);
      });

      var length = pairs.length;
      callback = baseCreateCallback(callback, thisArg, 3);
      while (length--) {
        if (callback(pairs[length--], pairs[length], object) === false) {
          break;
        }
      }
      return object;
    }

    /**
     * Iterates over own enumerable properties of an object, executing the callback
     * for each property. The callback is bound to `thisArg` and invoked with three
     * arguments; (value, key, object). Callbacks may exit iteration early by
     * explicitly returning `false`.
     *
     * @static
     * @memberOf _
     * @type Function
     * @category Objects
     * @param {Object} object The object to iterate over.
     * @param {Function} [callback=identity] The function called per iteration.
     * @param {*} [thisArg] The `this` binding of `callback`.
     * @returns {Object} Returns `object`.
     * @example
     *
     * _.forOwn({ '0': 'zero', '1': 'one', 'length': 2 }, function(num, key) {
     *   console.log(key);
     * });
     * // => logs '0', '1', and 'length' (property order is not guaranteed across environments)
     */
    var forOwn = function(collection, callback, thisArg) {
      var index, iterable = collection, result = iterable;
      if (!iterable) return result;
      if (!objectTypes[typeof iterable]) return result;
      callback = callback && typeof thisArg == 'undefined' ? callback : baseCreateCallback(callback, thisArg, 3);
        var ownIndex = -1,
            ownProps = objectTypes[typeof iterable] && keys(iterable),
            length = ownProps ? ownProps.length : 0;

        while (++ownIndex < length) {
          index = ownProps[ownIndex];
          if (callback(iterable[index], index, collection) === false) return result;
        }
      return result
    };

    /**
     * This method is like `_.forOwn` except that it iterates over elements
     * of a `collection` in the opposite order.
     *
     * @static
     * @memberOf _
     * @category Objects
     * @param {Object} object The object to iterate over.
     * @param {Function} [callback=identity] The function called per iteration.
     * @param {*} [thisArg] The `this` binding of `callback`.
     * @returns {Object} Returns `object`.
     * @example
     *
     * _.forOwnRight({ '0': 'zero', '1': 'one', 'length': 2 }, function(num, key) {
     *   console.log(key);
     * });
     * // => logs 'length', '1', and '0' assuming `_.forOwn` logs '0', '1', and 'length'
     */
    function forOwnRight(object, callback, thisArg) {
      var props = keys(object),
          length = props.length;

      callback = baseCreateCallback(callback, thisArg, 3);
      while (length--) {
        var key = props[length];
        if (callback(object[key], key, object) === false) {
          break;
        }
      }
      return object;
    }

    /**
     * Creates a sorted array of property names of all enumerable properties,
     * own and inherited, of `object` that have function values.
     *
     * @static
     * @memberOf _
     * @alias methods
     * @category Objects
     * @param {Object} object The object to inspect.
     * @returns {Array} Returns an array of property names that have function values.
     * @example
     *
     * _.functions(_);
     * // => ['all', 'any', 'bind', 'bindAll', 'clone', 'compact', 'compose', ...]
     */
    function functions(object) {
      var result = [];
      forIn(object, function(value, key) {
        if (isFunction(value)) {
          result.push(key);
        }
      });
      return result.sort();
    }

    /**
     * Checks if the specified property name exists as a direct property of `object`,
     * instead of an inherited property.
     *
     * @static
     * @memberOf _
     * @category Objects
     * @param {Object} object The object to inspect.
     * @param {string} key The name of the property to check.
     * @returns {boolean} Returns `true` if key is a direct property, else `false`.
     * @example
     *
     * _.has({ 'a': 1, 'b': 2, 'c': 3 }, 'b');
     * // => true
     */
    function has(object, key) {
      return object ? hasOwnProperty.call(object, key) : false;
    }

    /**
     * Creates an object composed of the inverted keys and values of the given object.
     *
     * @static
     * @memberOf _
     * @category Objects
     * @param {Object} object The object to invert.
     * @returns {Object} Returns the created inverted object.
     * @example
     *
     * _.invert({ 'first': 'fred', 'second': 'barney' });
     * // => { 'fred': 'first', 'barney': 'second' }
     */
    function invert(object) {
      var index = -1,
          props = keys(object),
          length = props.length,
          result = {};

      while (++index < length) {
        var key = props[index];
        result[object[key]] = key;
      }
      return result;
    }

    /**
     * Checks if `value` is a boolean value.
     *
     * @static
     * @memberOf _
     * @category Objects
     * @param {*} value The value to check.
     * @returns {boolean} Returns `true` if the `value` is a boolean value, else `false`.
     * @example
     *
     * _.isBoolean(null);
     * // => false
     */
    function isBoolean(value) {
      return value === true || value === false ||
        value && typeof value == 'object' && toString.call(value) == boolClass || false;
    }

    /**
     * Checks if `value` is a date.
     *
     * @static
     * @memberOf _
     * @category Objects
     * @param {*} value The value to check.
     * @returns {boolean} Returns `true` if the `value` is a date, else `false`.
     * @example
     *
     * _.isDate(new Date);
     * // => true
     */
    function isDate(value) {
      return value && typeof value == 'object' && toString.call(value) == dateClass || false;
    }

    /**
     * Checks if `value` is a DOM element.
     *
     * @static
     * @memberOf _
     * @category Objects
     * @param {*} value The value to check.
     * @returns {boolean} Returns `true` if the `value` is a DOM element, else `false`.
     * @example
     *
     * _.isElement(document.body);
     * // => true
     */
    function isElement(value) {
      return value && value.nodeType === 1 || false;
    }

    /**
     * Checks if `value` is empty. Arrays, strings, or `arguments` objects with a
     * length of `0` and objects with no own enumerable properties are considered
     * "empty".
     *
     * @static
     * @memberOf _
     * @category Objects
     * @param {Array|Object|string} value The value to inspect.
     * @returns {boolean} Returns `true` if the `value` is empty, else `false`.
     * @example
     *
     * _.isEmpty([1, 2, 3]);
     * // => false
     *
     * _.isEmpty({});
     * // => true
     *
     * _.isEmpty('');
     * // => true
     */
    function isEmpty(value) {
      var result = true;
      if (!value) {
        return result;
      }
      var className = toString.call(value),
          length = value.length;

      if ((className == arrayClass || className == stringClass || className == argsClass ) ||
          (className == objectClass && typeof length == 'number' && isFunction(value.splice))) {
        return !length;
      }
      forOwn(value, function() {
        return (result = false);
      });
      return result;
    }

    /**
     * Performs a deep comparison between two values to determine if they are
     * equivalent to each other. If a callback is provided it will be executed
     * to compare values. If the callback returns `undefined` comparisons will
     * be handled by the method instead. The callback is bound to `thisArg` and
     * invoked with two arguments; (a, b).
     *
     * @static
     * @memberOf _
     * @category Objects
     * @param {*} a The value to compare.
     * @param {*} b The other value to compare.
     * @param {Function} [callback] The function to customize comparing values.
     * @param {*} [thisArg] The `this` binding of `callback`.
     * @returns {boolean} Returns `true` if the values are equivalent, else `false`.
     * @example
     *
     * var object = { 'name': 'fred' };
     * var copy = { 'name': 'fred' };
     *
     * object == copy;
     * // => false
     *
     * _.isEqual(object, copy);
     * // => true
     *
     * var words = ['hello', 'goodbye'];
     * var otherWords = ['hi', 'goodbye'];
     *
     * _.isEqual(words, otherWords, function(a, b) {
     *   var reGreet = /^(?:hello|hi)$/i,
     *       aGreet = _.isString(a) && reGreet.test(a),
     *       bGreet = _.isString(b) && reGreet.test(b);
     *
     *   return (aGreet || bGreet) ? (aGreet == bGreet) : undefined;
     * });
     * // => true
     */
    function isEqual(a, b, callback, thisArg) {
      return baseIsEqual(a, b, typeof callback == 'function' && baseCreateCallback(callback, thisArg, 2));
    }

    /**
     * Checks if `value` is, or can be coerced to, a finite number.
     *
     * Note: This is not the same as native `isFinite` which will return true for
     * booleans and empty strings. See http://es5.github.io/#x15.1.2.5.
     *
     * @static
     * @memberOf _
     * @category Objects
     * @param {*} value The value to check.
     * @returns {boolean} Returns `true` if the `value` is finite, else `false`.
     * @example
     *
     * _.isFinite(-101);
     * // => true
     *
     * _.isFinite('10');
     * // => true
     *
     * _.isFinite(true);
     * // => false
     *
     * _.isFinite('');
     * // => false
     *
     * _.isFinite(Infinity);
     * // => false
     */
    function isFinite(value) {
      return nativeIsFinite(value) && !nativeIsNaN(parseFloat(value));
    }

    /**
     * Checks if `value` is a function.
     *
     * @static
     * @memberOf _
     * @category Objects
     * @param {*} value The value to check.
     * @returns {boolean} Returns `true` if the `value` is a function, else `false`.
     * @example
     *
     * _.isFunction(_);
     * // => true
     */
    function isFunction(value) {
      return typeof value == 'function';
    }

    /**
     * Checks if `value` is the language type of Object.
     * (e.g. arrays, functions, objects, regexes, `new Number(0)`, and `new String('')`)
     *
     * @static
     * @memberOf _
     * @category Objects
     * @param {*} value The value to check.
     * @returns {boolean} Returns `true` if the `value` is an object, else `false`.
     * @example
     *
     * _.isObject({});
     * // => true
     *
     * _.isObject([1, 2, 3]);
     * // => true
     *
     * _.isObject(1);
     * // => false
     */
    function isObject(value) {
      // check if the value is the ECMAScript language type of Object
      // http://es5.github.io/#x8
      // and avoid a V8 bug
      // http://code.google.com/p/v8/issues/detail?id=2291
      return !!(value && objectTypes[typeof value]);
    }

    /**
     * Checks if `value` is `NaN`.
     *
     * Note: This is not the same as native `isNaN` which will return `true` for
     * `undefined` and other non-numeric values. See http://es5.github.io/#x15.1.2.4.
     *
     * @static
     * @memberOf _
     * @category Objects
     * @param {*} value The value to check.
     * @returns {boolean} Returns `true` if the `value` is `NaN`, else `false`.
     * @example
     *
     * _.isNaN(NaN);
     * // => true
     *
     * _.isNaN(new Number(NaN));
     * // => true
     *
     * isNaN(undefined);
     * // => true
     *
     * _.isNaN(undefined);
     * // => false
     */
    function isNaN(value) {
      // `NaN` as a primitive is the only value that is not equal to itself
      // (perform the [[Class]] check first to avoid errors with some host objects in IE)
      return isNumber(value) && value != +value;
    }

    /**
     * Checks if `value` is `null`.
     *
     * @static
     * @memberOf _
     * @category Objects
     * @param {*} value The value to check.
     * @returns {boolean} Returns `true` if the `value` is `null`, else `false`.
     * @example
     *
     * _.isNull(null);
     * // => true
     *
     * _.isNull(undefined);
     * // => false
     */
    function isNull(value) {
      return value === null;
    }

    /**
     * Checks if `value` is a number.
     *
     * Note: `NaN` is considered a number. See http://es5.github.io/#x8.5.
     *
     * @static
     * @memberOf _
     * @category Objects
     * @param {*} value The value to check.
     * @returns {boolean} Returns `true` if the `value` is a number, else `false`.
     * @example
     *
     * _.isNumber(8.4 * 5);
     * // => true
     */
    function isNumber(value) {
      return typeof value == 'number' ||
        value && typeof value == 'object' && toString.call(value) == numberClass || false;
    }

    /**
     * Checks if `value` is an object created by the `Object` constructor.
     *
     * @static
     * @memberOf _
     * @category Objects
     * @param {*} value The value to check.
     * @returns {boolean} Returns `true` if `value` is a plain object, else `false`.
     * @example
     *
     * function Shape() {
     *   this.x = 0;
     *   this.y = 0;
     * }
     *
     * _.isPlainObject(new Shape);
     * // => false
     *
     * _.isPlainObject([1, 2, 3]);
     * // => false
     *
     * _.isPlainObject({ 'x': 0, 'y': 0 });
     * // => true
     */
    var isPlainObject = !getPrototypeOf ? shimIsPlainObject : function(value) {
      if (!(value && toString.call(value) == objectClass)) {
        return false;
      }
      var valueOf = value.valueOf,
          objProto = isNative(valueOf) && (objProto = getPrototypeOf(valueOf)) && getPrototypeOf(objProto);

      return objProto
        ? (value == objProto || getPrototypeOf(value) == objProto)
        : shimIsPlainObject(value);
    };

    /**
     * Checks if `value` is a regular expression.
     *
     * @static
     * @memberOf _
     * @category Objects
     * @param {*} value The value to check.
     * @returns {boolean} Returns `true` if the `value` is a regular expression, else `false`.
     * @example
     *
     * _.isRegExp(/fred/);
     * // => true
     */
    function isRegExp(value) {
      return value && typeof value == 'object' && toString.call(value) == regexpClass || false;
    }

    /**
     * Checks if `value` is a string.
     *
     * @static
     * @memberOf _
     * @category Objects
     * @param {*} value The value to check.
     * @returns {boolean} Returns `true` if the `value` is a string, else `false`.
     * @example
     *
     * _.isString('fred');
     * // => true
     */
    function isString(value) {
      return typeof value == 'string' ||
        value && typeof value == 'object' && toString.call(value) == stringClass || false;
    }

    /**
     * Checks if `value` is `undefined`.
     *
     * @static
     * @memberOf _
     * @category Objects
     * @param {*} value The value to check.
     * @returns {boolean} Returns `true` if the `value` is `undefined`, else `false`.
     * @example
     *
     * _.isUndefined(void 0);
     * // => true
     */
    function isUndefined(value) {
      return typeof value == 'undefined';
    }

    /**
     * Creates an object with the same keys as `object` and values generated by
     * running each own enumerable property of `object` through the callback.
     * The callback is bound to `thisArg` and invoked with three arguments;
     * (value, key, object).
     *
     * If a property name is provided for `callback` the created "_.pluck" style
     * callback will return the property value of the given element.
     *
     * If an object is provided for `callback` the created "_.where" style callback
     * will return `true` for elements that have the properties of the given object,
     * else `false`.
     *
     * @static
     * @memberOf _
     * @category Objects
     * @param {Object} object The object to iterate over.
     * @param {Function|Object|string} [callback=identity] The function called
     *  per iteration. If a property name or object is provided it will be used
     *  to create a "_.pluck" or "_.where" style callback, respectively.
     * @param {*} [thisArg] The `this` binding of `callback`.
     * @returns {Array} Returns a new object with values of the results of each `callback` execution.
     * @example
     *
     * _.mapValues({ 'a': 1, 'b': 2, 'c': 3} , function(num) { return num * 3; });
     * // => { 'a': 3, 'b': 6, 'c': 9 }
     *
     * var characters = {
     *   'fred': { 'name': 'fred', 'age': 40 },
     *   'pebbles': { 'name': 'pebbles', 'age': 1 }
     * };
     *
     * // using "_.pluck" callback shorthand
     * _.mapValues(characters, 'age');
     * // => { 'fred': 40, 'pebbles': 1 }
     */
    function mapValues(object, callback, thisArg) {
      var result = {};
      callback = lodash.createCallback(callback, thisArg, 3);

      forOwn(object, function(value, key, object) {
        result[key] = callback(value, key, object);
      });
      return result;
    }

    /**
     * Recursively merges own enumerable properties of the source object(s), that
     * don't resolve to `undefined` into the destination object. Subsequent sources
     * will overwrite property assignments of previous sources. If a callback is
     * provided it will be executed to produce the merged values of the destination
     * and source properties. If the callback returns `undefined` merging will
     * be handled by the method instead. The callback is bound to `thisArg` and
     * invoked with two arguments; (objectValue, sourceValue).
     *
     * @static
     * @memberOf _
     * @category Objects
     * @param {Object} object The destination object.
     * @param {...Object} [source] The source objects.
     * @param {Function} [callback] The function to customize merging properties.
     * @param {*} [thisArg] The `this` binding of `callback`.
     * @returns {Object} Returns the destination object.
     * @example
     *
     * var names = {
     *   'characters': [
     *     { 'name': 'barney' },
     *     { 'name': 'fred' }
     *   ]
     * };
     *
     * var ages = {
     *   'characters': [
     *     { 'age': 36 },
     *     { 'age': 40 }
     *   ]
     * };
     *
     * _.merge(names, ages);
     * // => { 'characters': [{ 'name': 'barney', 'age': 36 }, { 'name': 'fred', 'age': 40 }] }
     *
     * var food = {
     *   'fruits': ['apple'],
     *   'vegetables': ['beet']
     * };
     *
     * var otherFood = {
     *   'fruits': ['banana'],
     *   'vegetables': ['carrot']
     * };
     *
     * _.merge(food, otherFood, function(a, b) {
     *   return _.isArray(a) ? a.concat(b) : undefined;
     * });
     * // => { 'fruits': ['apple', 'banana'], 'vegetables': ['beet', 'carrot] }
     */
    function merge(object) {
      var args = arguments,
          length = 2;

      if (!isObject(object)) {
        return object;
      }
      // allows working with `_.reduce` and `_.reduceRight` without using
      // their `index` and `collection` arguments
      if (typeof args[2] != 'number') {
        length = args.length;
      }
      if (length > 3 && typeof args[length - 2] == 'function') {
        var callback = baseCreateCallback(args[--length - 1], args[length--], 2);
      } else if (length > 2 && typeof args[length - 1] == 'function') {
        callback = args[--length];
      }
      var sources = slice(arguments, 1, length),
          index = -1,
          stackA = getArray(),
          stackB = getArray();

      while (++index < length) {
        baseMerge(object, sources[index], callback, stackA, stackB);
      }
      releaseArray(stackA);
      releaseArray(stackB);
      return object;
    }

    /**
     * Creates a shallow clone of `object` excluding the specified properties.
     * Property names may be specified as individual arguments or as arrays of
     * property names. If a callback is provided it will be executed for each
     * property of `object` omitting the properties the callback returns truey
     * for. The callback is bound to `thisArg` and invoked with three arguments;
     * (value, key, object).
     *
     * @static
     * @memberOf _
     * @category Objects
     * @param {Object} object The source object.
     * @param {Function|...string|string[]} [callback] The properties to omit or the
     *  function called per iteration.
     * @param {*} [thisArg] The `this` binding of `callback`.
     * @returns {Object} Returns an object without the omitted properties.
     * @example
     *
     * _.omit({ 'name': 'fred', 'age': 40 }, 'age');
     * // => { 'name': 'fred' }
     *
     * _.omit({ 'name': 'fred', 'age': 40 }, function(value) {
     *   return typeof value == 'number';
     * });
     * // => { 'name': 'fred' }
     */
    function omit(object, callback, thisArg) {
      var result = {};
      if (typeof callback != 'function') {
        var props = [];
        forIn(object, function(value, key) {
          props.push(key);
        });
        props = baseDifference(props, baseFlatten(arguments, true, false, 1));

        var index = -1,
            length = props.length;

        while (++index < length) {
          var key = props[index];
          result[key] = object[key];
        }
      } else {
        callback = lodash.createCallback(callback, thisArg, 3);
        forIn(object, function(value, key, object) {
          if (!callback(value, key, object)) {
            result[key] = value;
          }
        });
      }
      return result;
    }

    /**
     * Creates a two dimensional array of an object's key-value pairs,
     * i.e. `[[key1, value1], [key2, value2]]`.
     *
     * @static
     * @memberOf _
     * @category Objects
     * @param {Object} object The object to inspect.
     * @returns {Array} Returns new array of key-value pairs.
     * @example
     *
     * _.pairs({ 'barney': 36, 'fred': 40 });
     * // => [['barney', 36], ['fred', 40]] (property order is not guaranteed across environments)
     */
    function pairs(object) {
      var index = -1,
          props = keys(object),
          length = props.length,
          result = Array(length);

      while (++index < length) {
        var key = props[index];
        result[index] = [key, object[key]];
      }
      return result;
    }

    /**
     * Creates a shallow clone of `object` composed of the specified properties.
     * Property names may be specified as individual arguments or as arrays of
     * property names. If a callback is provided it will be executed for each
     * property of `object` picking the properties the callback returns truey
     * for. The callback is bound to `thisArg` and invoked with three arguments;
     * (value, key, object).
     *
     * @static
     * @memberOf _
     * @category Objects
     * @param {Object} object The source object.
     * @param {Function|...string|string[]} [callback] The function called per
     *  iteration or property names to pick, specified as individual property
     *  names or arrays of property names.
     * @param {*} [thisArg] The `this` binding of `callback`.
     * @returns {Object} Returns an object composed of the picked properties.
     * @example
     *
     * _.pick({ 'name': 'fred', '_userid': 'fred1' }, 'name');
     * // => { 'name': 'fred' }
     *
     * _.pick({ 'name': 'fred', '_userid': 'fred1' }, function(value, key) {
     *   return key.charAt(0) != '_';
     * });
     * // => { 'name': 'fred' }
     */
    function pick(object, callback, thisArg) {
      var result = {};
      if (typeof callback != 'function') {
        var index = -1,
            props = baseFlatten(arguments, true, false, 1),
            length = isObject(object) ? props.length : 0;

        while (++index < length) {
          var key = props[index];
          if (key in object) {
            result[key] = object[key];
          }
        }
      } else {
        callback = lodash.createCallback(callback, thisArg, 3);
        forIn(object, function(value, key, object) {
          if (callback(value, key, object)) {
            result[key] = value;
          }
        });
      }
      return result;
    }

    /**
     * An alternative to `_.reduce` this method transforms `object` to a new
     * `accumulator` object which is the result of running each of its own
     * enumerable properties through a callback, with each callback execution
     * potentially mutating the `accumulator` object. The callback is bound to
     * `thisArg` and invoked with four arguments; (accumulator, value, key, object).
     * Callbacks may exit iteration early by explicitly returning `false`.
     *
     * @static
     * @memberOf _
     * @category Objects
     * @param {Array|Object} object The object to iterate over.
     * @param {Function} [callback=identity] The function called per iteration.
     * @param {*} [accumulator] The custom accumulator value.
     * @param {*} [thisArg] The `this` binding of `callback`.
     * @returns {*} Returns the accumulated value.
     * @example
     *
     * var squares = _.transform([1, 2, 3, 4, 5, 6, 7, 8, 9, 10], function(result, num) {
     *   num *= num;
     *   if (num % 2) {
     *     return result.push(num) < 3;
     *   }
     * });
     * // => [1, 9, 25]
     *
     * var mapped = _.transform({ 'a': 1, 'b': 2, 'c': 3 }, function(result, num, key) {
     *   result[key] = num * 3;
     * });
     * // => { 'a': 3, 'b': 6, 'c': 9 }
     */
    function transform(object, callback, accumulator, thisArg) {
      var isArr = isArray(object);
      if (accumulator == null) {
        if (isArr) {
          accumulator = [];
        } else {
          var ctor = object && object.constructor,
              proto = ctor && ctor.prototype;

          accumulator = baseCreate(proto);
        }
      }
      if (callback) {
        callback = lodash.createCallback(callback, thisArg, 4);
        (isArr ? forEach : forOwn)(object, function(value, index, object) {
          return callback(accumulator, value, index, object);
        });
      }
      return accumulator;
    }

    /**
     * Creates an array composed of the own enumerable property values of `object`.
     *
     * @static
     * @memberOf _
     * @category Objects
     * @param {Object} object The object to inspect.
     * @returns {Array} Returns an array of property values.
     * @example
     *
     * _.values({ 'one': 1, 'two': 2, 'three': 3 });
     * // => [1, 2, 3] (property order is not guaranteed across environments)
     */
    function values(object) {
      var index = -1,
          props = keys(object),
          length = props.length,
          result = Array(length);

      while (++index < length) {
        result[index] = object[props[index]];
      }
      return result;
    }

    /*--------------------------------------------------------------------------*/

    /**
     * Creates an array of elements from the specified indexes, or keys, of the
     * `collection`. Indexes may be specified as individual arguments or as arrays
     * of indexes.
     *
     * @static
     * @memberOf _
     * @category Collections
     * @param {Array|Object|string} collection The collection to iterate over.
     * @param {...(number|number[]|string|string[])} [index] The indexes of `collection`
     *   to retrieve, specified as individual indexes or arrays of indexes.
     * @returns {Array} Returns a new array of elements corresponding to the
     *  provided indexes.
     * @example
     *
     * _.at(['a', 'b', 'c', 'd', 'e'], [0, 2, 4]);
     * // => ['a', 'c', 'e']
     *
     * _.at(['fred', 'barney', 'pebbles'], 0, 2);
     * // => ['fred', 'pebbles']
     */
    function at(collection) {
      var args = arguments,
          index = -1,
          props = baseFlatten(args, true, false, 1),
          length = (args[2] && args[2][args[1]] === collection) ? 1 : props.length,
          result = Array(length);

      while(++index < length) {
        result[index] = collection[props[index]];
      }
      return result;
    }

    /**
     * Checks if a given value is present in a collection using strict equality
     * for comparisons, i.e. `===`. If `fromIndex` is negative, it is used as the
     * offset from the end of the collection.
     *
     * @static
     * @memberOf _
     * @alias include
     * @category Collections
     * @param {Array|Object|string} collection The collection to iterate over.
     * @param {*} target The value to check for.
     * @param {number} [fromIndex=0] The index to search from.
     * @returns {boolean} Returns `true` if the `target` element is found, else `false`.
     * @example
     *
     * _.contains([1, 2, 3], 1);
     * // => true
     *
     * _.contains([1, 2, 3], 1, 2);
     * // => false
     *
     * _.contains({ 'name': 'fred', 'age': 40 }, 'fred');
     * // => true
     *
     * _.contains('pebbles', 'eb');
     * // => true
     */
    function contains(collection, target, fromIndex) {
      var index = -1,
          indexOf = getIndexOf(),
          length = collection ? collection.length : 0,
          result = false;

      fromIndex = (fromIndex < 0 ? nativeMax(0, length + fromIndex) : fromIndex) || 0;
      if (isArray(collection)) {
        result = indexOf(collection, target, fromIndex) > -1;
      } else if (typeof length == 'number') {
        result = (isString(collection) ? collection.indexOf(target, fromIndex) : indexOf(collection, target, fromIndex)) > -1;
      } else {
        forOwn(collection, function(value) {
          if (++index >= fromIndex) {
            return !(result = value === target);
          }
        });
      }
      return result;
    }

    /**
     * Creates an object composed of keys generated from the results of running
     * each element of `collection` through the callback. The corresponding value
     * of each key is the number of times the key was returned by the callback.
     * The callback is bound to `thisArg` and invoked with three arguments;
     * (value, index|key, collection).
     *
     * If a property name is provided for `callback` the created "_.pluck" style
     * callback will return the property value of the given element.
     *
     * If an object is provided for `callback` the created "_.where" style callback
     * will return `true` for elements that have the properties of the given object,
     * else `false`.
     *
     * @static
     * @memberOf _
     * @category Collections
     * @param {Array|Object|string} collection The collection to iterate over.
     * @param {Function|Object|string} [callback=identity] The function called
     *  per iteration. If a property name or object is provided it will be used
     *  to create a "_.pluck" or "_.where" style callback, respectively.
     * @param {*} [thisArg] The `this` binding of `callback`.
     * @returns {Object} Returns the composed aggregate object.
     * @example
     *
     * _.countBy([4.3, 6.1, 6.4], function(num) { return Math.floor(num); });
     * // => { '4': 1, '6': 2 }
     *
     * _.countBy([4.3, 6.1, 6.4], function(num) { return this.floor(num); }, Math);
     * // => { '4': 1, '6': 2 }
     *
     * _.countBy(['one', 'two', 'three'], 'length');
     * // => { '3': 2, '5': 1 }
     */
    var countBy = createAggregator(function(result, value, key) {
      (hasOwnProperty.call(result, key) ? result[key]++ : result[key] = 1);
    });

    /**
     * Checks if the given callback returns truey value for **all** elements of
     * a collection. The callback is bound to `thisArg` and invoked with three
     * arguments; (value, index|key, collection).
     *
     * If a property name is provided for `callback` the created "_.pluck" style
     * callback will return the property value of the given element.
     *
     * If an object is provided for `callback` the created "_.where" style callback
     * will return `true` for elements that have the properties of the given object,
     * else `false`.
     *
     * @static
     * @memberOf _
     * @alias all
     * @category Collections
     * @param {Array|Object|string} collection The collection to iterate over.
     * @param {Function|Object|string} [callback=identity] The function called
     *  per iteration. If a property name or object is provided it will be used
     *  to create a "_.pluck" or "_.where" style callback, respectively.
     * @param {*} [thisArg] The `this` binding of `callback`.
     * @returns {boolean} Returns `true` if all elements passed the callback check,
     *  else `false`.
     * @example
     *
     * _.every([true, 1, null, 'yes']);
     * // => false
     *
     * var characters = [
     *   { 'name': 'barney', 'age': 36 },
     *   { 'name': 'fred',   'age': 40 }
     * ];
     *
     * // using "_.pluck" callback shorthand
     * _.every(characters, 'age');
     * // => true
     *
     * // using "_.where" callback shorthand
     * _.every(characters, { 'age': 36 });
     * // => false
     */
    function every(collection, callback, thisArg) {
      var result = true;
      callback = lodash.createCallback(callback, thisArg, 3);

      var index = -1,
          length = collection ? collection.length : 0;

      if (typeof length == 'number') {
        while (++index < length) {
          if (!(result = !!callback(collection[index], index, collection))) {
            break;
          }
        }
      } else {
        forOwn(collection, function(value, index, collection) {
          return (result = !!callback(value, index, collection));
        });
      }
      return result;
    }

    /**
     * Iterates over elements of a collection, returning an array of all elements
     * the callback returns truey for. The callback is bound to `thisArg` and
     * invoked with three arguments; (value, index|key, collection).
     *
     * If a property name is provided for `callback` the created "_.pluck" style
     * callback will return the property value of the given element.
     *
     * If an object is provided for `callback` the created "_.where" style callback
     * will return `true` for elements that have the properties of the given object,
     * else `false`.
     *
     * @static
     * @memberOf _
     * @alias select
     * @category Collections
     * @param {Array|Object|string} collection The collection to iterate over.
     * @param {Function|Object|string} [callback=identity] The function called
     *  per iteration. If a property name or object is provided it will be used
     *  to create a "_.pluck" or "_.where" style callback, respectively.
     * @param {*} [thisArg] The `this` binding of `callback`.
     * @returns {Array} Returns a new array of elements that passed the callback check.
     * @example
     *
     * var evens = _.filter([1, 2, 3, 4, 5, 6], function(num) { return num % 2 == 0; });
     * // => [2, 4, 6]
     *
     * var characters = [
     *   { 'name': 'barney', 'age': 36, 'blocked': false },
     *   { 'name': 'fred',   'age': 40, 'blocked': true }
     * ];
     *
     * // using "_.pluck" callback shorthand
     * _.filter(characters, 'blocked');
     * // => [{ 'name': 'fred', 'age': 40, 'blocked': true }]
     *
     * // using "_.where" callback shorthand
     * _.filter(characters, { 'age': 36 });
     * // => [{ 'name': 'barney', 'age': 36, 'blocked': false }]
     */
    function filter(collection, callback, thisArg) {
      var result = [];
      callback = lodash.createCallback(callback, thisArg, 3);

      var index = -1,
          length = collection ? collection.length : 0;

      if (typeof length == 'number') {
        while (++index < length) {
          var value = collection[index];
          if (callback(value, index, collection)) {
            result.push(value);
          }
        }
      } else {
        forOwn(collection, function(value, index, collection) {
          if (callback(value, index, collection)) {
            result.push(value);
          }
        });
      }
      return result;
    }

    /**
     * Iterates over elements of a collection, returning the first element that
     * the callback returns truey for. The callback is bound to `thisArg` and
     * invoked with three arguments; (value, index|key, collection).
     *
     * If a property name is provided for `callback` the created "_.pluck" style
     * callback will return the property value of the given element.
     *
     * If an object is provided for `callback` the created "_.where" style callback
     * will return `true` for elements that have the properties of the given object,
     * else `false`.
     *
     * @static
     * @memberOf _
     * @alias detect, findWhere
     * @category Collections
     * @param {Array|Object|string} collection The collection to iterate over.
     * @param {Function|Object|string} [callback=identity] The function called
     *  per iteration. If a property name or object is provided it will be used
     *  to create a "_.pluck" or "_.where" style callback, respectively.
     * @param {*} [thisArg] The `this` binding of `callback`.
     * @returns {*} Returns the found element, else `undefined`.
     * @example
     *
     * var characters = [
     *   { 'name': 'barney',  'age': 36, 'blocked': false },
     *   { 'name': 'fred',    'age': 40, 'blocked': true },
     *   { 'name': 'pebbles', 'age': 1,  'blocked': false }
     * ];
     *
     * _.find(characters, function(chr) {
     *   return chr.age < 40;
     * });
     * // => { 'name': 'barney', 'age': 36, 'blocked': false }
     *
     * // using "_.where" callback shorthand
     * _.find(characters, { 'age': 1 });
     * // =>  { 'name': 'pebbles', 'age': 1, 'blocked': false }
     *
     * // using "_.pluck" callback shorthand
     * _.find(characters, 'blocked');
     * // => { 'name': 'fred', 'age': 40, 'blocked': true }
     */
    function find(collection, callback, thisArg) {
      callback = lodash.createCallback(callback, thisArg, 3);

      var index = -1,
          length = collection ? collection.length : 0;

      if (typeof length == 'number') {
        while (++index < length) {
          var value = collection[index];
          if (callback(value, index, collection)) {
            return value;
          }
        }
      } else {
        var result;
        forOwn(collection, function(value, index, collection) {
          if (callback(value, index, collection)) {
            result = value;
            return false;
          }
        });
        return result;
      }
    }

    /**
     * This method is like `_.find` except that it iterates over elements
     * of a `collection` from right to left.
     *
     * @static
     * @memberOf _
     * @category Collections
     * @param {Array|Object|string} collection The collection to iterate over.
     * @param {Function|Object|string} [callback=identity] The function called
     *  per iteration. If a property name or object is provided it will be used
     *  to create a "_.pluck" or "_.where" style callback, respectively.
     * @param {*} [thisArg] The `this` binding of `callback`.
     * @returns {*} Returns the found element, else `undefined`.
     * @example
     *
     * _.findLast([1, 2, 3, 4], function(num) {
     *   return num % 2 == 1;
     * });
     * // => 3
     */
    function findLast(collection, callback, thisArg) {
      var result;
      callback = lodash.createCallback(callback, thisArg, 3);
      forEachRight(collection, function(value, index, collection) {
        if (callback(value, index, collection)) {
          result = value;
          return false;
        }
      });
      return result;
    }

    /**
     * Iterates over elements of a collection, executing the callback for each
     * element. The callback is bound to `thisArg` and invoked with three arguments;
     * (value, index|key, collection). Callbacks may exit iteration early by
     * explicitly returning `false`.
     *
     * Note: As with other "Collections" methods, objects with a `length` property
     * are iterated like arrays. To avoid this behavior `_.forIn` or `_.forOwn`
     * may be used for object iteration.
     *
     * @static
     * @memberOf _
     * @alias each
     * @category Collections
     * @param {Array|Object|string} collection The collection to iterate over.
     * @param {Function} [callback=identity] The function called per iteration.
     * @param {*} [thisArg] The `this` binding of `callback`.
     * @returns {Array|Object|string} Returns `collection`.
     * @example
     *
     * _([1, 2, 3]).forEach(function(num) { console.log(num); }).join(',');
     * // => logs each number and returns '1,2,3'
     *
     * _.forEach({ 'one': 1, 'two': 2, 'three': 3 }, function(num) { console.log(num); });
     * // => logs each number and returns the object (property order is not guaranteed across environments)
     */
    function forEach(collection, callback, thisArg) {
      var index = -1,
          length = collection ? collection.length : 0;

      callback = callback && typeof thisArg == 'undefined' ? callback : baseCreateCallback(callback, thisArg, 3);
      if (typeof length == 'number') {
        while (++index < length) {
          if (callback(collection[index], index, collection) === false) {
            break;
          }
        }
      } else {
        forOwn(collection, callback);
      }
      return collection;
    }

    /**
     * This method is like `_.forEach` except that it iterates over elements
     * of a `collection` from right to left.
     *
     * @static
     * @memberOf _
     * @alias eachRight
     * @category Collections
     * @param {Array|Object|string} collection The collection to iterate over.
     * @param {Function} [callback=identity] The function called per iteration.
     * @param {*} [thisArg] The `this` binding of `callback`.
     * @returns {Array|Object|string} Returns `collection`.
     * @example
     *
     * _([1, 2, 3]).forEachRight(function(num) { console.log(num); }).join(',');
     * // => logs each number from right to left and returns '3,2,1'
     */
    function forEachRight(collection, callback, thisArg) {
      var length = collection ? collection.length : 0;
      callback = callback && typeof thisArg == 'undefined' ? callback : baseCreateCallback(callback, thisArg, 3);
      if (typeof length == 'number') {
        while (length--) {
          if (callback(collection[length], length, collection) === false) {
            break;
          }
        }
      } else {
        var props = keys(collection);
        length = props.length;
        forOwn(collection, function(value, key, collection) {
          key = props ? props[--length] : --length;
          return callback(collection[key], key, collection);
        });
      }
      return collection;
    }

    /**
     * Creates an object composed of keys generated from the results of running
     * each element of a collection through the callback. The corresponding value
     * of each key is an array of the elements responsible for generating the key.
     * The callback is bound to `thisArg` and invoked with three arguments;
     * (value, index|key, collection).
     *
     * If a property name is provided for `callback` the created "_.pluck" style
     * callback will return the property value of the given element.
     *
     * If an object is provided for `callback` the created "_.where" style callback
     * will return `true` for elements that have the properties of the given object,
     * else `false`
     *
     * @static
     * @memberOf _
     * @category Collections
     * @param {Array|Object|string} collection The collection to iterate over.
     * @param {Function|Object|string} [callback=identity] The function called
     *  per iteration. If a property name or object is provided it will be used
     *  to create a "_.pluck" or "_.where" style callback, respectively.
     * @param {*} [thisArg] The `this` binding of `callback`.
     * @returns {Object} Returns the composed aggregate object.
     * @example
     *
     * _.groupBy([4.2, 6.1, 6.4], function(num) { return Math.floor(num); });
     * // => { '4': [4.2], '6': [6.1, 6.4] }
     *
     * _.groupBy([4.2, 6.1, 6.4], function(num) { return this.floor(num); }, Math);
     * // => { '4': [4.2], '6': [6.1, 6.4] }
     *
     * // using "_.pluck" callback shorthand
     * _.groupBy(['one', 'two', 'three'], 'length');
     * // => { '3': ['one', 'two'], '5': ['three'] }
     */
    var groupBy = createAggregator(function(result, value, key) {
      (hasOwnProperty.call(result, key) ? result[key] : result[key] = []).push(value);
    });

    /**
     * Creates an object composed of keys generated from the results of running
     * each element of the collection through the given callback. The corresponding
     * value of each key is the last element responsible for generating the key.
     * The callback is bound to `thisArg` and invoked with three arguments;
     * (value, index|key, collection).
     *
     * If a property name is provided for `callback` the created "_.pluck" style
     * callback will return the property value of the given element.
     *
     * If an object is provided for `callback` the created "_.where" style callback
     * will return `true` for elements that have the properties of the given object,
     * else `false`.
     *
     * @static
     * @memberOf _
     * @category Collections
     * @param {Array|Object|string} collection The collection to iterate over.
     * @param {Function|Object|string} [callback=identity] The function called
     *  per iteration. If a property name or object is provided it will be used
     *  to create a "_.pluck" or "_.where" style callback, respectively.
     * @param {*} [thisArg] The `this` binding of `callback`.
     * @returns {Object} Returns the composed aggregate object.
     * @example
     *
     * var keys = [
     *   { 'dir': 'left', 'code': 97 },
     *   { 'dir': 'right', 'code': 100 }
     * ];
     *
     * _.indexBy(keys, 'dir');
     * // => { 'left': { 'dir': 'left', 'code': 97 }, 'right': { 'dir': 'right', 'code': 100 } }
     *
     * _.indexBy(keys, function(key) { return String.fromCharCode(key.code); });
     * // => { 'a': { 'dir': 'left', 'code': 97 }, 'd': { 'dir': 'right', 'code': 100 } }
     *
     * _.indexBy(characters, function(key) { this.fromCharCode(key.code); }, String);
     * // => { 'a': { 'dir': 'left', 'code': 97 }, 'd': { 'dir': 'right', 'code': 100 } }
     */
    var indexBy = createAggregator(function(result, value, key) {
      result[key] = value;
    });

    /**
     * Invokes the method named by `methodName` on each element in the `collection`
     * returning an array of the results of each invoked method. Additional arguments
     * will be provided to each invoked method. If `methodName` is a function it
     * will be invoked for, and `this` bound to, each element in the `collection`.
     *
     * @static
     * @memberOf _
     * @category Collections
     * @param {Array|Object|string} collection The collection to iterate over.
     * @param {Function|string} methodName The name of the method to invoke or
     *  the function invoked per iteration.
     * @param {...*} [arg] Arguments to invoke the method with.
     * @returns {Array} Returns a new array of the results of each invoked method.
     * @example
     *
     * _.invoke([[5, 1, 7], [3, 2, 1]], 'sort');
     * // => [[1, 5, 7], [1, 2, 3]]
     *
     * _.invoke([123, 456], String.prototype.split, '');
     * // => [['1', '2', '3'], ['4', '5', '6']]
     */
    function invoke(collection, methodName) {
      var args = slice(arguments, 2),
          index = -1,
          isFunc = typeof methodName == 'function',
          length = collection ? collection.length : 0,
          result = Array(typeof length == 'number' ? length : 0);

      forEach(collection, function(value) {
        result[++index] = (isFunc ? methodName : value[methodName]).apply(value, args);
      });
      return result;
    }

    /**
     * Creates an array of values by running each element in the collection
     * through the callback. The callback is bound to `thisArg` and invoked with
     * three arguments; (value, index|key, collection).
     *
     * If a property name is provided for `callback` the created "_.pluck" style
     * callback will return the property value of the given element.
     *
     * If an object is provided for `callback` the created "_.where" style callback
     * will return `true` for elements that have the properties of the given object,
     * else `false`.
     *
     * @static
     * @memberOf _
     * @alias collect
     * @category Collections
     * @param {Array|Object|string} collection The collection to iterate over.
     * @param {Function|Object|string} [callback=identity] The function called
     *  per iteration. If a property name or object is provided it will be used
     *  to create a "_.pluck" or "_.where" style callback, respectively.
     * @param {*} [thisArg] The `this` binding of `callback`.
     * @returns {Array} Returns a new array of the results of each `callback` execution.
     * @example
     *
     * _.map([1, 2, 3], function(num) { return num * 3; });
     * // => [3, 6, 9]
     *
     * _.map({ 'one': 1, 'two': 2, 'three': 3 }, function(num) { return num * 3; });
     * // => [3, 6, 9] (property order is not guaranteed across environments)
     *
     * var characters = [
     *   { 'name': 'barney', 'age': 36 },
     *   { 'name': 'fred',   'age': 40 }
     * ];
     *
     * // using "_.pluck" callback shorthand
     * _.map(characters, 'name');
     * // => ['barney', 'fred']
     */
    function map(collection, callback, thisArg) {
      var index = -1,
          length = collection ? collection.length : 0;

      callback = lodash.createCallback(callback, thisArg, 3);
      if (typeof length == 'number') {
        var result = Array(length);
        while (++index < length) {
          result[index] = callback(collection[index], index, collection);
        }
      } else {
        result = [];
        forOwn(collection, function(value, key, collection) {
          result[++index] = callback(value, key, collection);
        });
      }
      return result;
    }

    /**
     * Retrieves the maximum value of a collection. If the collection is empty or
     * falsey `-Infinity` is returned. If a callback is provided it will be executed
     * for each value in the collection to generate the criterion by which the value
     * is ranked. The callback is bound to `thisArg` and invoked with three
     * arguments; (value, index, collection).
     *
     * If a property name is provided for `callback` the created "_.pluck" style
     * callback will return the property value of the given element.
     *
     * If an object is provided for `callback` the created "_.where" style callback
     * will return `true` for elements that have the properties of the given object,
     * else `false`.
     *
     * @static
     * @memberOf _
     * @category Collections
     * @param {Array|Object|string} collection The collection to iterate over.
     * @param {Function|Object|string} [callback=identity] The function called
     *  per iteration. If a property name or object is provided it will be used
     *  to create a "_.pluck" or "_.where" style callback, respectively.
     * @param {*} [thisArg] The `this` binding of `callback`.
     * @returns {*} Returns the maximum value.
     * @example
     *
     * _.max([4, 2, 8, 6]);
     * // => 8
     *
     * var characters = [
     *   { 'name': 'barney', 'age': 36 },
     *   { 'name': 'fred',   'age': 40 }
     * ];
     *
     * _.max(characters, function(chr) { return chr.age; });
     * // => { 'name': 'fred', 'age': 40 };
     *
     * // using "_.pluck" callback shorthand
     * _.max(characters, 'age');
     * // => { 'name': 'fred', 'age': 40 };
     */
    function max(collection, callback, thisArg) {
      var computed = -Infinity,
          result = computed;

      // allows working with functions like `_.map` without using
      // their `index` argument as a callback
      if (typeof callback != 'function' && thisArg && thisArg[callback] === collection) {
        callback = null;
      }
      if (callback == null && isArray(collection)) {
        var index = -1,
            length = collection.length;

        while (++index < length) {
          var value = collection[index];
          if (value > result) {
            result = value;
          }
        }
      } else {
        callback = (callback == null && isString(collection))
          ? charAtCallback
          : lodash.createCallback(callback, thisArg, 3);

        forEach(collection, function(value, index, collection) {
          var current = callback(value, index, collection);
          if (current > computed) {
            computed = current;
            result = value;
          }
        });
      }
      return result;
    }

    /**
     * Retrieves the minimum value of a collection. If the collection is empty or
     * falsey `Infinity` is returned. If a callback is provided it will be executed
     * for each value in the collection to generate the criterion by which the value
     * is ranked. The callback is bound to `thisArg` and invoked with three
     * arguments; (value, index, collection).
     *
     * If a property name is provided for `callback` the created "_.pluck" style
     * callback will return the property value of the given element.
     *
     * If an object is provided for `callback` the created "_.where" style callback
     * will return `true` for elements that have the properties of the given object,
     * else `false`.
     *
     * @static
     * @memberOf _
     * @category Collections
     * @param {Array|Object|string} collection The collection to iterate over.
     * @param {Function|Object|string} [callback=identity] The function called
     *  per iteration. If a property name or object is provided it will be used
     *  to create a "_.pluck" or "_.where" style callback, respectively.
     * @param {*} [thisArg] The `this` binding of `callback`.
     * @returns {*} Returns the minimum value.
     * @example
     *
     * _.min([4, 2, 8, 6]);
     * // => 2
     *
     * var characters = [
     *   { 'name': 'barney', 'age': 36 },
     *   { 'name': 'fred',   'age': 40 }
     * ];
     *
     * _.min(characters, function(chr) { return chr.age; });
     * // => { 'name': 'barney', 'age': 36 };
     *
     * // using "_.pluck" callback shorthand
     * _.min(characters, 'age');
     * // => { 'name': 'barney', 'age': 36 };
     */
    function min(collection, callback, thisArg) {
      var computed = Infinity,
          result = computed;

      // allows working with functions like `_.map` without using
      // their `index` argument as a callback
      if (typeof callback != 'function' && thisArg && thisArg[callback] === collection) {
        callback = null;
      }
      if (callback == null && isArray(collection)) {
        var index = -1,
            length = collection.length;

        while (++index < length) {
          var value = collection[index];
          if (value < result) {
            result = value;
          }
        }
      } else {
        callback = (callback == null && isString(collection))
          ? charAtCallback
          : lodash.createCallback(callback, thisArg, 3);

        forEach(collection, function(value, index, collection) {
          var current = callback(value, index, collection);
          if (current < computed) {
            computed = current;
            result = value;
          }
        });
      }
      return result;
    }

    /**
     * Retrieves the value of a specified property from all elements in the collection.
     *
     * @static
     * @memberOf _
     * @type Function
     * @category Collections
     * @param {Array|Object|string} collection The collection to iterate over.
     * @param {string} property The name of the property to pluck.
     * @returns {Array} Returns a new array of property values.
     * @example
     *
     * var characters = [
     *   { 'name': 'barney', 'age': 36 },
     *   { 'name': 'fred',   'age': 40 }
     * ];
     *
     * _.pluck(characters, 'name');
     * // => ['barney', 'fred']
     */
    var pluck = map;

    /**
     * Reduces a collection to a value which is the accumulated result of running
     * each element in the collection through the callback, where each successive
     * callback execution consumes the return value of the previous execution. If
     * `accumulator` is not provided the first element of the collection will be
     * used as the initial `accumulator` value. The callback is bound to `thisArg`
     * and invoked with four arguments; (accumulator, value, index|key, collection).
     *
     * @static
     * @memberOf _
     * @alias foldl, inject
     * @category Collections
     * @param {Array|Object|string} collection The collection to iterate over.
     * @param {Function} [callback=identity] The function called per iteration.
     * @param {*} [accumulator] Initial value of the accumulator.
     * @param {*} [thisArg] The `this` binding of `callback`.
     * @returns {*} Returns the accumulated value.
     * @example
     *
     * var sum = _.reduce([1, 2, 3], function(sum, num) {
     *   return sum + num;
     * });
     * // => 6
     *
     * var mapped = _.reduce({ 'a': 1, 'b': 2, 'c': 3 }, function(result, num, key) {
     *   result[key] = num * 3;
     *   return result;
     * }, {});
     * // => { 'a': 3, 'b': 6, 'c': 9 }
     */
    function reduce(collection, callback, accumulator, thisArg) {
      if (!collection) return accumulator;
      var noaccum = arguments.length < 3;
      callback = lodash.createCallback(callback, thisArg, 4);

      var index = -1,
          length = collection.length;

      if (typeof length == 'number') {
        if (noaccum) {
          accumulator = collection[++index];
        }
        while (++index < length) {
          accumulator = callback(accumulator, collection[index], index, collection);
        }
      } else {
        forOwn(collection, function(value, index, collection) {
          accumulator = noaccum
            ? (noaccum = false, value)
            : callback(accumulator, value, index, collection)
        });
      }
      return accumulator;
    }

    /**
     * This method is like `_.reduce` except that it iterates over elements
     * of a `collection` from right to left.
     *
     * @static
     * @memberOf _
     * @alias foldr
     * @category Collections
     * @param {Array|Object|string} collection The collection to iterate over.
     * @param {Function} [callback=identity] The function called per iteration.
     * @param {*} [accumulator] Initial value of the accumulator.
     * @param {*} [thisArg] The `this` binding of `callback`.
     * @returns {*} Returns the accumulated value.
     * @example
     *
     * var list = [[0, 1], [2, 3], [4, 5]];
     * var flat = _.reduceRight(list, function(a, b) { return a.concat(b); }, []);
     * // => [4, 5, 2, 3, 0, 1]
     */
    function reduceRight(collection, callback, accumulator, thisArg) {
      var noaccum = arguments.length < 3;
      callback = lodash.createCallback(callback, thisArg, 4);
      forEachRight(collection, function(value, index, collection) {
        accumulator = noaccum
          ? (noaccum = false, value)
          : callback(accumulator, value, index, collection);
      });
      return accumulator;
    }

    /**
     * The opposite of `_.filter` this method returns the elements of a
     * collection that the callback does **not** return truey for.
     *
     * If a property name is provided for `callback` the created "_.pluck" style
     * callback will return the property value of the given element.
     *
     * If an object is provided for `callback` the created "_.where" style callback
     * will return `true` for elements that have the properties of the given object,
     * else `false`.
     *
     * @static
     * @memberOf _
     * @category Collections
     * @param {Array|Object|string} collection The collection to iterate over.
     * @param {Function|Object|string} [callback=identity] The function called
     *  per iteration. If a property name or object is provided it will be used
     *  to create a "_.pluck" or "_.where" style callback, respectively.
     * @param {*} [thisArg] The `this` binding of `callback`.
     * @returns {Array} Returns a new array of elements that failed the callback check.
     * @example
     *
     * var odds = _.reject([1, 2, 3, 4, 5, 6], function(num) { return num % 2 == 0; });
     * // => [1, 3, 5]
     *
     * var characters = [
     *   { 'name': 'barney', 'age': 36, 'blocked': false },
     *   { 'name': 'fred',   'age': 40, 'blocked': true }
     * ];
     *
     * // using "_.pluck" callback shorthand
     * _.reject(characters, 'blocked');
     * // => [{ 'name': 'barney', 'age': 36, 'blocked': false }]
     *
     * // using "_.where" callback shorthand
     * _.reject(characters, { 'age': 36 });
     * // => [{ 'name': 'fred', 'age': 40, 'blocked': true }]
     */
    function reject(collection, callback, thisArg) {
      callback = lodash.createCallback(callback, thisArg, 3);
      return filter(collection, function(value, index, collection) {
        return !callback(value, index, collection);
      });
    }

    /**
     * Retrieves a random element or `n` random elements from a collection.
     *
     * @static
     * @memberOf _
     * @category Collections
     * @param {Array|Object|string} collection The collection to sample.
     * @param {number} [n] The number of elements to sample.
     * @param- {Object} [guard] Allows working with functions like `_.map`
     *  without using their `index` arguments as `n`.
     * @returns {Array} Returns the random sample(s) of `collection`.
     * @example
     *
     * _.sample([1, 2, 3, 4]);
     * // => 2
     *
     * _.sample([1, 2, 3, 4], 2);
     * // => [3, 1]
     */
    function sample(collection, n, guard) {
      if (collection && typeof collection.length != 'number') {
        collection = values(collection);
      }
      if (n == null || guard) {
        return collection ? collection[baseRandom(0, collection.length - 1)] : undefined;
      }
      var result = shuffle(collection);
      result.length = nativeMin(nativeMax(0, n), result.length);
      return result;
    }

    /**
     * Creates an array of shuffled values, using a version of the Fisher-Yates
     * shuffle. See http://en.wikipedia.org/wiki/Fisher-Yates_shuffle.
     *
     * @static
     * @memberOf _
     * @category Collections
     * @param {Array|Object|string} collection The collection to shuffle.
     * @returns {Array} Returns a new shuffled collection.
     * @example
     *
     * _.shuffle([1, 2, 3, 4, 5, 6]);
     * // => [4, 1, 6, 3, 5, 2]
     */
    function shuffle(collection) {
      var index = -1,
          length = collection ? collection.length : 0,
          result = Array(typeof length == 'number' ? length : 0);

      forEach(collection, function(value) {
        var rand = baseRandom(0, ++index);
        result[index] = result[rand];
        result[rand] = value;
      });
      return result;
    }

    /**
     * Gets the size of the `collection` by returning `collection.length` for arrays
     * and array-like objects or the number of own enumerable properties for objects.
     *
     * @static
     * @memberOf _
     * @category Collections
     * @param {Array|Object|string} collection The collection to inspect.
     * @returns {number} Returns `collection.length` or number of own enumerable properties.
     * @example
     *
     * _.size([1, 2]);
     * // => 2
     *
     * _.size({ 'one': 1, 'two': 2, 'three': 3 });
     * // => 3
     *
     * _.size('pebbles');
     * // => 7
     */
    function size(collection) {
      var length = collection ? collection.length : 0;
      return typeof length == 'number' ? length : keys(collection).length;
    }

    /**
     * Checks if the callback returns a truey value for **any** element of a
     * collection. The function returns as soon as it finds a passing value and
     * does not iterate over the entire collection. The callback is bound to
     * `thisArg` and invoked with three arguments; (value, index|key, collection).
     *
     * If a property name is provided for `callback` the created "_.pluck" style
     * callback will return the property value of the given element.
     *
     * If an object is provided for `callback` the created "_.where" style callback
     * will return `true` for elements that have the properties of the given object,
     * else `false`.
     *
     * @static
     * @memberOf _
     * @alias any
     * @category Collections
     * @param {Array|Object|string} collection The collection to iterate over.
     * @param {Function|Object|string} [callback=identity] The function called
     *  per iteration. If a property name or object is provided it will be used
     *  to create a "_.pluck" or "_.where" style callback, respectively.
     * @param {*} [thisArg] The `this` binding of `callback`.
     * @returns {boolean} Returns `true` if any element passed the callback check,
     *  else `false`.
     * @example
     *
     * _.some([null, 0, 'yes', false], Boolean);
     * // => true
     *
     * var characters = [
     *   { 'name': 'barney', 'age': 36, 'blocked': false },
     *   { 'name': 'fred',   'age': 40, 'blocked': true }
     * ];
     *
     * // using "_.pluck" callback shorthand
     * _.some(characters, 'blocked');
     * // => true
     *
     * // using "_.where" callback shorthand
     * _.some(characters, { 'age': 1 });
     * // => false
     */
    function some(collection, callback, thisArg) {
      var result;
      callback = lodash.createCallback(callback, thisArg, 3);

      var index = -1,
          length = collection ? collection.length : 0;

      if (typeof length == 'number') {
        while (++index < length) {
          if ((result = callback(collection[index], index, collection))) {
            break;
          }
        }
      } else {
        forOwn(collection, function(value, index, collection) {
          return !(result = callback(value, index, collection));
        });
      }
      return !!result;
    }

    /**
     * Creates an array of elements, sorted in ascending order by the results of
     * running each element in a collection through the callback. This method
     * performs a stable sort, that is, it will preserve the original sort order
     * of equal elements. The callback is bound to `thisArg` and invoked with
     * three arguments; (value, index|key, collection).
     *
     * If a property name is provided for `callback` the created "_.pluck" style
     * callback will return the property value of the given element.
     *
     * If an array of property names is provided for `callback` the collection
     * will be sorted by each property value.
     *
     * If an object is provided for `callback` the created "_.where" style callback
     * will return `true` for elements that have the properties of the given object,
     * else `false`.
     *
     * @static
     * @memberOf _
     * @category Collections
     * @param {Array|Object|string} collection The collection to iterate over.
     * @param {Array|Function|Object|string} [callback=identity] The function called
     *  per iteration. If a property name or object is provided it will be used
     *  to create a "_.pluck" or "_.where" style callback, respectively.
     * @param {*} [thisArg] The `this` binding of `callback`.
     * @returns {Array} Returns a new array of sorted elements.
     * @example
     *
     * _.sortBy([1, 2, 3], function(num) { return Math.sin(num); });
     * // => [3, 1, 2]
     *
     * _.sortBy([1, 2, 3], function(num) { return this.sin(num); }, Math);
     * // => [3, 1, 2]
     *
     * var characters = [
     *   { 'name': 'barney',  'age': 36 },
     *   { 'name': 'fred',    'age': 40 },
     *   { 'name': 'barney',  'age': 26 },
     *   { 'name': 'fred',    'age': 30 }
     * ];
     *
     * // using "_.pluck" callback shorthand
     * _.map(_.sortBy(characters, 'age'), _.values);
     * // => [['barney', 26], ['fred', 30], ['barney', 36], ['fred', 40]]
     *
     * // sorting by multiple properties
     * _.map(_.sortBy(characters, ['name', 'age']), _.values);
     * // = > [['barney', 26], ['barney', 36], ['fred', 30], ['fred', 40]]
     */
    function sortBy(collection, callback, thisArg) {
      var index = -1,
          isArr = isArray(callback),
          length = collection ? collection.length : 0,
          result = Array(typeof length == 'number' ? length : 0);

      if (!isArr) {
        callback = lodash.createCallback(callback, thisArg, 3);
      }
      forEach(collection, function(value, key, collection) {
        var object = result[++index] = getObject();
        if (isArr) {
          object.criteria = map(callback, function(key) { return value[key]; });
        } else {
          (object.criteria = getArray())[0] = callback(value, key, collection);
        }
        object.index = index;
        object.value = value;
      });

      length = result.length;
      result.sort(compareAscending);
      while (length--) {
        var object = result[length];
        result[length] = object.value;
        if (!isArr) {
          releaseArray(object.criteria);
        }
        releaseObject(object);
      }
      return result;
    }

    /**
     * Converts the `collection` to an array.
     *
     * @static
     * @memberOf _
     * @category Collections
     * @param {Array|Object|string} collection The collection to convert.
     * @returns {Array} Returns the new converted array.
     * @example
     *
     * (function() { return _.toArray(arguments).slice(1); })(1, 2, 3, 4);
     * // => [2, 3, 4]
     */
    function toArray(collection) {
      if (collection && typeof collection.length == 'number') {
        return slice(collection);
      }
      return values(collection);
    }

    /**
     * Performs a deep comparison of each element in a `collection` to the given
     * `properties` object, returning an array of all elements that have equivalent
     * property values.
     *
     * @static
     * @memberOf _
     * @type Function
     * @category Collections
     * @param {Array|Object|string} collection The collection to iterate over.
     * @param {Object} props The object of property values to filter by.
     * @returns {Array} Returns a new array of elements that have the given properties.
     * @example
     *
     * var characters = [
     *   { 'name': 'barney', 'age': 36, 'pets': ['hoppy'] },
     *   { 'name': 'fred',   'age': 40, 'pets': ['baby puss', 'dino'] }
     * ];
     *
     * _.where(characters, { 'age': 36 });
     * // => [{ 'name': 'barney', 'age': 36, 'pets': ['hoppy'] }]
     *
     * _.where(characters, { 'pets': ['dino'] });
     * // => [{ 'name': 'fred', 'age': 40, 'pets': ['baby puss', 'dino'] }]
     */
    var where = filter;

    /*--------------------------------------------------------------------------*/

    /**
     * Creates an array with all falsey values removed. The values `false`, `null`,
     * `0`, `""`, `undefined`, and `NaN` are all falsey.
     *
     * @static
     * @memberOf _
     * @category Arrays
     * @param {Array} array The array to compact.
     * @returns {Array} Returns a new array of filtered values.
     * @example
     *
     * _.compact([0, 1, false, 2, '', 3]);
     * // => [1, 2, 3]
     */
    function compact(array) {
      var index = -1,
          length = array ? array.length : 0,
          result = [];

      while (++index < length) {
        var value = array[index];
        if (value) {
          result.push(value);
        }
      }
      return result;
    }

    /**
     * Creates an array excluding all values of the provided arrays using strict
     * equality for comparisons, i.e. `===`.
     *
     * @static
     * @memberOf _
     * @category Arrays
     * @param {Array} array The array to process.
     * @param {...Array} [values] The arrays of values to exclude.
     * @returns {Array} Returns a new array of filtered values.
     * @example
     *
     * _.difference([1, 2, 3, 4, 5], [5, 2, 10]);
     * // => [1, 3, 4]
     */
    function difference(array) {
      return baseDifference(array, baseFlatten(arguments, true, true, 1));
    }

    /**
     * This method is like `_.find` except that it returns the index of the first
     * element that passes the callback check, instead of the element itself.
     *
     * If a property name is provided for `callback` the created "_.pluck" style
     * callback will return the property value of the given element.
     *
     * If an object is provided for `callback` the created "_.where" style callback
     * will return `true` for elements that have the properties of the given object,
     * else `false`.
     *
     * @static
     * @memberOf _
     * @category Arrays
     * @param {Array} array The array to search.
     * @param {Function|Object|string} [callback=identity] The function called
     *  per iteration. If a property name or object is provided it will be used
     *  to create a "_.pluck" or "_.where" style callback, respectively.
     * @param {*} [thisArg] The `this` binding of `callback`.
     * @returns {number} Returns the index of the found element, else `-1`.
     * @example
     *
     * var characters = [
     *   { 'name': 'barney',  'age': 36, 'blocked': false },
     *   { 'name': 'fred',    'age': 40, 'blocked': true },
     *   { 'name': 'pebbles', 'age': 1,  'blocked': false }
     * ];
     *
     * _.findIndex(characters, function(chr) {
     *   return chr.age < 20;
     * });
     * // => 2
     *
     * // using "_.where" callback shorthand
     * _.findIndex(characters, { 'age': 36 });
     * // => 0
     *
     * // using "_.pluck" callback shorthand
     * _.findIndex(characters, 'blocked');
     * // => 1
     */
    function findIndex(array, callback, thisArg) {
      var index = -1,
          length = array ? array.length : 0;

      callback = lodash.createCallback(callback, thisArg, 3);
      while (++index < length) {
        if (callback(array[index], index, array)) {
          return index;
        }
      }
      return -1;
    }

    /**
     * This method is like `_.findIndex` except that it iterates over elements
     * of a `collection` from right to left.
     *
     * If a property name is provided for `callback` the created "_.pluck" style
     * callback will return the property value of the given element.
     *
     * If an object is provided for `callback` the created "_.where" style callback
     * will return `true` for elements that have the properties of the given object,
     * else `false`.
     *
     * @static
     * @memberOf _
     * @category Arrays
     * @param {Array} array The array to search.
     * @param {Function|Object|string} [callback=identity] The function called
     *  per iteration. If a property name or object is provided it will be used
     *  to create a "_.pluck" or "_.where" style callback, respectively.
     * @param {*} [thisArg] The `this` binding of `callback`.
     * @returns {number} Returns the index of the found element, else `-1`.
     * @example
     *
     * var characters = [
     *   { 'name': 'barney',  'age': 36, 'blocked': true },
     *   { 'name': 'fred',    'age': 40, 'blocked': false },
     *   { 'name': 'pebbles', 'age': 1,  'blocked': true }
     * ];
     *
     * _.findLastIndex(characters, function(chr) {
     *   return chr.age > 30;
     * });
     * // => 1
     *
     * // using "_.where" callback shorthand
     * _.findLastIndex(characters, { 'age': 36 });
     * // => 0
     *
     * // using "_.pluck" callback shorthand
     * _.findLastIndex(characters, 'blocked');
     * // => 2
     */
    function findLastIndex(array, callback, thisArg) {
      var length = array ? array.length : 0;
      callback = lodash.createCallback(callback, thisArg, 3);
      while (length--) {
        if (callback(array[length], length, array)) {
          return length;
        }
      }
      return -1;
    }

    /**
     * Gets the first element or first `n` elements of an array. If a callback
     * is provided elements at the beginning of the array are returned as long
     * as the callback returns truey. The callback is bound to `thisArg` and
     * invoked with three arguments; (value, index, array).
     *
     * If a property name is provided for `callback` the created "_.pluck" style
     * callback will return the property value of the given element.
     *
     * If an object is provided for `callback` the created "_.where" style callback
     * will return `true` for elements that have the properties of the given object,
     * else `false`.
     *
     * @static
     * @memberOf _
     * @alias head, take
     * @category Arrays
     * @param {Array} array The array to query.
     * @param {Function|Object|number|string} [callback] The function called
     *  per element or the number of elements to return. If a property name or
     *  object is provided it will be used to create a "_.pluck" or "_.where"
     *  style callback, respectively.
     * @param {*} [thisArg] The `this` binding of `callback`.
     * @returns {*} Returns the first element(s) of `array`.
     * @example
     *
     * _.first([1, 2, 3]);
     * // => 1
     *
     * _.first([1, 2, 3], 2);
     * // => [1, 2]
     *
     * _.first([1, 2, 3], function(num) {
     *   return num < 3;
     * });
     * // => [1, 2]
     *
     * var characters = [
     *   { 'name': 'barney',  'blocked': true,  'employer': 'slate' },
     *   { 'name': 'fred',    'blocked': false, 'employer': 'slate' },
     *   { 'name': 'pebbles', 'blocked': true,  'employer': 'na' }
     * ];
     *
     * // using "_.pluck" callback shorthand
     * _.first(characters, 'blocked');
     * // => [{ 'name': 'barney', 'blocked': true, 'employer': 'slate' }]
     *
     * // using "_.where" callback shorthand
     * _.pluck(_.first(characters, { 'employer': 'slate' }), 'name');
     * // => ['barney', 'fred']
     */
    function first(array, callback, thisArg) {
      var n = 0,
          length = array ? array.length : 0;

      if (typeof callback != 'number' && callback != null) {
        var index = -1;
        callback = lodash.createCallback(callback, thisArg, 3);
        while (++index < length && callback(array[index], index, array)) {
          n++;
        }
      } else {
        n = callback;
        if (n == null || thisArg) {
          return array ? array[0] : undefined;
        }
      }
      return slice(array, 0, nativeMin(nativeMax(0, n), length));
    }

    /**
     * Flattens a nested array (the nesting can be to any depth). If `isShallow`
     * is truey, the array will only be flattened a single level. If a callback
     * is provided each element of the array is passed through the callback before
     * flattening. The callback is bound to `thisArg` and invoked with three
     * arguments; (value, index, array).
     *
     * If a property name is provided for `callback` the created "_.pluck" style
     * callback will return the property value of the given element.
     *
     * If an object is provided for `callback` the created "_.where" style callback
     * will return `true` for elements that have the properties of the given object,
     * else `false`.
     *
     * @static
     * @memberOf _
     * @category Arrays
     * @param {Array} array The array to flatten.
     * @param {boolean} [isShallow=false] A flag to restrict flattening to a single level.
     * @param {Function|Object|string} [callback=identity] The function called
     *  per iteration. If a property name or object is provided it will be used
     *  to create a "_.pluck" or "_.where" style callback, respectively.
     * @param {*} [thisArg] The `this` binding of `callback`.
     * @returns {Array} Returns a new flattened array.
     * @example
     *
     * _.flatten([1, [2], [3, [[4]]]]);
     * // => [1, 2, 3, 4];
     *
     * _.flatten([1, [2], [3, [[4]]]], true);
     * // => [1, 2, 3, [[4]]];
     *
     * var characters = [
     *   { 'name': 'barney', 'age': 30, 'pets': ['hoppy'] },
     *   { 'name': 'fred',   'age': 40, 'pets': ['baby puss', 'dino'] }
     * ];
     *
     * // using "_.pluck" callback shorthand
     * _.flatten(characters, 'pets');
     * // => ['hoppy', 'baby puss', 'dino']
     */
    function flatten(array, isShallow, callback, thisArg) {
      // juggle arguments
      if (typeof isShallow != 'boolean' && isShallow != null) {
        thisArg = callback;
        callback = (typeof isShallow != 'function' && thisArg && thisArg[isShallow] === array) ? null : isShallow;
        isShallow = false;
      }
      if (callback != null) {
        array = map(array, callback, thisArg);
      }
      return baseFlatten(array, isShallow);
    }

    /**
     * Gets the index at which the first occurrence of `value` is found using
     * strict equality for comparisons, i.e. `===`. If the array is already sorted
     * providing `true` for `fromIndex` will run a faster binary search.
     *
     * @static
     * @memberOf _
     * @category Arrays
     * @param {Array} array The array to search.
     * @param {*} value The value to search for.
     * @param {boolean|number} [fromIndex=0] The index to search from or `true`
     *  to perform a binary search on a sorted array.
     * @returns {number} Returns the index of the matched value or `-1`.
     * @example
     *
     * _.indexOf([1, 2, 3, 1, 2, 3], 2);
     * // => 1
     *
     * _.indexOf([1, 2, 3, 1, 2, 3], 2, 3);
     * // => 4
     *
     * _.indexOf([1, 1, 2, 2, 3, 3], 2, true);
     * // => 2
     */
    function indexOf(array, value, fromIndex) {
      if (typeof fromIndex == 'number') {
        var length = array ? array.length : 0;
        fromIndex = (fromIndex < 0 ? nativeMax(0, length + fromIndex) : fromIndex || 0);
      } else if (fromIndex) {
        var index = sortedIndex(array, value);
        return array[index] === value ? index : -1;
      }
      return baseIndexOf(array, value, fromIndex);
    }

    /**
     * Gets all but the last element or last `n` elements of an array. If a
     * callback is provided elements at the end of the array are excluded from
     * the result as long as the callback returns truey. The callback is bound
     * to `thisArg` and invoked with three arguments; (value, index, array).
     *
     * If a property name is provided for `callback` the created "_.pluck" style
     * callback will return the property value of the given element.
     *
     * If an object is provided for `callback` the created "_.where" style callback
     * will return `true` for elements that have the properties of the given object,
     * else `false`.
     *
     * @static
     * @memberOf _
     * @category Arrays
     * @param {Array} array The array to query.
     * @param {Function|Object|number|string} [callback=1] The function called
     *  per element or the number of elements to exclude. If a property name or
     *  object is provided it will be used to create a "_.pluck" or "_.where"
     *  style callback, respectively.
     * @param {*} [thisArg] The `this` binding of `callback`.
     * @returns {Array} Returns a slice of `array`.
     * @example
     *
     * _.initial([1, 2, 3]);
     * // => [1, 2]
     *
     * _.initial([1, 2, 3], 2);
     * // => [1]
     *
     * _.initial([1, 2, 3], function(num) {
     *   return num > 1;
     * });
     * // => [1]
     *
     * var characters = [
     *   { 'name': 'barney',  'blocked': false, 'employer': 'slate' },
     *   { 'name': 'fred',    'blocked': true,  'employer': 'slate' },
     *   { 'name': 'pebbles', 'blocked': true,  'employer': 'na' }
     * ];
     *
     * // using "_.pluck" callback shorthand
     * _.initial(characters, 'blocked');
     * // => [{ 'name': 'barney',  'blocked': false, 'employer': 'slate' }]
     *
     * // using "_.where" callback shorthand
     * _.pluck(_.initial(characters, { 'employer': 'na' }), 'name');
     * // => ['barney', 'fred']
     */
    function initial(array, callback, thisArg) {
      var n = 0,
          length = array ? array.length : 0;

      if (typeof callback != 'number' && callback != null) {
        var index = length;
        callback = lodash.createCallback(callback, thisArg, 3);
        while (index-- && callback(array[index], index, array)) {
          n++;
        }
      } else {
        n = (callback == null || thisArg) ? 1 : callback || n;
      }
      return slice(array, 0, nativeMin(nativeMax(0, length - n), length));
    }

    /**
     * Creates an array of unique values present in all provided arrays using
     * strict equality for comparisons, i.e. `===`.
     *
     * @static
     * @memberOf _
     * @category Arrays
     * @param {...Array} [array] The arrays to inspect.
     * @returns {Array} Returns an array of shared values.
     * @example
     *
     * _.intersection([1, 2, 3], [5, 2, 1, 4], [2, 1]);
     * // => [1, 2]
     */
    function intersection() {
      var args = [],
          argsIndex = -1,
          argsLength = arguments.length,
          caches = getArray(),
          indexOf = getIndexOf(),
          trustIndexOf = indexOf === baseIndexOf,
          seen = getArray();

      while (++argsIndex < argsLength) {
        var value = arguments[argsIndex];
        if (isArray(value) || isArguments(value)) {
          args.push(value);
          caches.push(trustIndexOf && value.length >= largeArraySize &&
            createCache(argsIndex ? args[argsIndex] : seen));
        }
      }
      var array = args[0],
          index = -1,
          length = array ? array.length : 0,
          result = [];

      outer:
      while (++index < length) {
        var cache = caches[0];
        value = array[index];

        if ((cache ? cacheIndexOf(cache, value) : indexOf(seen, value)) < 0) {
          argsIndex = argsLength;
          (cache || seen).push(value);
          while (--argsIndex) {
            cache = caches[argsIndex];
            if ((cache ? cacheIndexOf(cache, value) : indexOf(args[argsIndex], value)) < 0) {
              continue outer;
            }
          }
          result.push(value);
        }
      }
      while (argsLength--) {
        cache = caches[argsLength];
        if (cache) {
          releaseObject(cache);
        }
      }
      releaseArray(caches);
      releaseArray(seen);
      return result;
    }

    /**
     * Gets the last element or last `n` elements of an array. If a callback is
     * provided elements at the end of the array are returned as long as the
     * callback returns truey. The callback is bound to `thisArg` and invoked
     * with three arguments; (value, index, array).
     *
     * If a property name is provided for `callback` the created "_.pluck" style
     * callback will return the property value of the given element.
     *
     * If an object is provided for `callback` the created "_.where" style callback
     * will return `true` for elements that have the properties of the given object,
     * else `false`.
     *
     * @static
     * @memberOf _
     * @category Arrays
     * @param {Array} array The array to query.
     * @param {Function|Object|number|string} [callback] The function called
     *  per element or the number of elements to return. If a property name or
     *  object is provided it will be used to create a "_.pluck" or "_.where"
     *  style callback, respectively.
     * @param {*} [thisArg] The `this` binding of `callback`.
     * @returns {*} Returns the last element(s) of `array`.
     * @example
     *
     * _.last([1, 2, 3]);
     * // => 3
     *
     * _.last([1, 2, 3], 2);
     * // => [2, 3]
     *
     * _.last([1, 2, 3], function(num) {
     *   return num > 1;
     * });
     * // => [2, 3]
     *
     * var characters = [
     *   { 'name': 'barney',  'blocked': false, 'employer': 'slate' },
     *   { 'name': 'fred',    'blocked': true,  'employer': 'slate' },
     *   { 'name': 'pebbles', 'blocked': true,  'employer': 'na' }
     * ];
     *
     * // using "_.pluck" callback shorthand
     * _.pluck(_.last(characters, 'blocked'), 'name');
     * // => ['fred', 'pebbles']
     *
     * // using "_.where" callback shorthand
     * _.last(characters, { 'employer': 'na' });
     * // => [{ 'name': 'pebbles', 'blocked': true, 'employer': 'na' }]
     */
    function last(array, callback, thisArg) {
      var n = 0,
          length = array ? array.length : 0;

      if (typeof callback != 'number' && callback != null) {
        var index = length;
        callback = lodash.createCallback(callback, thisArg, 3);
        while (index-- && callback(array[index], index, array)) {
          n++;
        }
      } else {
        n = callback;
        if (n == null || thisArg) {
          return array ? array[length - 1] : undefined;
        }
      }
      return slice(array, nativeMax(0, length - n));
    }

    /**
     * Gets the index at which the last occurrence of `value` is found using strict
     * equality for comparisons, i.e. `===`. If `fromIndex` is negative, it is used
     * as the offset from the end of the collection.
     *
     * If a property name is provided for `callback` the created "_.pluck" style
     * callback will return the property value of the given element.
     *
     * If an object is provided for `callback` the created "_.where" style callback
     * will return `true` for elements that have the properties of the given object,
     * else `false`.
     *
     * @static
     * @memberOf _
     * @category Arrays
     * @param {Array} array The array to search.
     * @param {*} value The value to search for.
     * @param {number} [fromIndex=array.length-1] The index to search from.
     * @returns {number} Returns the index of the matched value or `-1`.
     * @example
     *
     * _.lastIndexOf([1, 2, 3, 1, 2, 3], 2);
     * // => 4
     *
     * _.lastIndexOf([1, 2, 3, 1, 2, 3], 2, 3);
     * // => 1
     */
    function lastIndexOf(array, value, fromIndex) {
      var index = array ? array.length : 0;
      if (typeof fromIndex == 'number') {
        index = (fromIndex < 0 ? nativeMax(0, index + fromIndex) : nativeMin(fromIndex, index - 1)) + 1;
      }
      while (index--) {
        if (array[index] === value) {
          return index;
        }
      }
      return -1;
    }

    /**
     * Removes all provided values from the given array using strict equality for
     * comparisons, i.e. `===`.
     *
     * @static
     * @memberOf _
     * @category Arrays
     * @param {Array} array The array to modify.
     * @param {...*} [value] The values to remove.
     * @returns {Array} Returns `array`.
     * @example
     *
     * var array = [1, 2, 3, 1, 2, 3];
     * _.pull(array, 2, 3);
     * console.log(array);
     * // => [1, 1]
     */
    function pull(array) {
      var args = arguments,
          argsIndex = 0,
          argsLength = args.length,
          length = array ? array.length : 0;

      while (++argsIndex < argsLength) {
        var index = -1,
            value = args[argsIndex];
        while (++index < length) {
          if (array[index] === value) {
            splice.call(array, index--, 1);
            length--;
          }
        }
      }
      return array;
    }

    /**
     * Creates an array of numbers (positive and/or negative) progressing from
     * `start` up to but not including `end`. If `start` is less than `stop` a
     * zero-length range is created unless a negative `step` is specified.
     *
     * @static
     * @memberOf _
     * @category Arrays
     * @param {number} [start=0] The start of the range.
     * @param {number} end The end of the range.
     * @param {number} [step=1] The value to increment or decrement by.
     * @returns {Array} Returns a new range array.
     * @example
     *
     * _.range(4);
     * // => [0, 1, 2, 3]
     *
     * _.range(1, 5);
     * // => [1, 2, 3, 4]
     *
     * _.range(0, 20, 5);
     * // => [0, 5, 10, 15]
     *
     * _.range(0, -4, -1);
     * // => [0, -1, -2, -3]
     *
     * _.range(1, 4, 0);
     * // => [1, 1, 1]
     *
     * _.range(0);
     * // => []
     */
    function range(start, end, step) {
      start = +start || 0;
      step = typeof step == 'number' ? step : (+step || 1);

      if (end == null) {
        end = start;
        start = 0;
      }
      // use `Array(length)` so engines like Chakra and V8 avoid slower modes
      // http://youtu.be/XAqIpGU8ZZk#t=17m25s
      var index = -1,
          length = nativeMax(0, ceil((end - start) / (step || 1))),
          result = Array(length);

      while (++index < length) {
        result[index] = start;
        start += step;
      }
      return result;
    }

    /**
     * Removes all elements from an array that the callback returns truey for
     * and returns an array of removed elements. The callback is bound to `thisArg`
     * and invoked with three arguments; (value, index, array).
     *
     * If a property name is provided for `callback` the created "_.pluck" style
     * callback will return the property value of the given element.
     *
     * If an object is provided for `callback` the created "_.where" style callback
     * will return `true` for elements that have the properties of the given object,
     * else `false`.
     *
     * @static
     * @memberOf _
     * @category Arrays
     * @param {Array} array The array to modify.
     * @param {Function|Object|string} [callback=identity] The function called
     *  per iteration. If a property name or object is provided it will be used
     *  to create a "_.pluck" or "_.where" style callback, respectively.
     * @param {*} [thisArg] The `this` binding of `callback`.
     * @returns {Array} Returns a new array of removed elements.
     * @example
     *
     * var array = [1, 2, 3, 4, 5, 6];
     * var evens = _.remove(array, function(num) { return num % 2 == 0; });
     *
     * console.log(array);
     * // => [1, 3, 5]
     *
     * console.log(evens);
     * // => [2, 4, 6]
     */
    function remove(array, callback, thisArg) {
      var index = -1,
          length = array ? array.length : 0,
          result = [];

      callback = lodash.createCallback(callback, thisArg, 3);
      while (++index < length) {
        var value = array[index];
        if (callback(value, index, array)) {
          result.push(value);
          splice.call(array, index--, 1);
          length--;
        }
      }
      return result;
    }

    /**
     * The opposite of `_.initial` this method gets all but the first element or
     * first `n` elements of an array. If a callback function is provided elements
     * at the beginning of the array are excluded from the result as long as the
     * callback returns truey. The callback is bound to `thisArg` and invoked
     * with three arguments; (value, index, array).
     *
     * If a property name is provided for `callback` the created "_.pluck" style
     * callback will return the property value of the given element.
     *
     * If an object is provided for `callback` the created "_.where" style callback
     * will return `true` for elements that have the properties of the given object,
     * else `false`.
     *
     * @static
     * @memberOf _
     * @alias drop, tail
     * @category Arrays
     * @param {Array} array The array to query.
     * @param {Function|Object|number|string} [callback=1] The function called
     *  per element or the number of elements to exclude. If a property name or
     *  object is provided it will be used to create a "_.pluck" or "_.where"
     *  style callback, respectively.
     * @param {*} [thisArg] The `this` binding of `callback`.
     * @returns {Array} Returns a slice of `array`.
     * @example
     *
     * _.rest([1, 2, 3]);
     * // => [2, 3]
     *
     * _.rest([1, 2, 3], 2);
     * // => [3]
     *
     * _.rest([1, 2, 3], function(num) {
     *   return num < 3;
     * });
     * // => [3]
     *
     * var characters = [
     *   { 'name': 'barney',  'blocked': true,  'employer': 'slate' },
     *   { 'name': 'fred',    'blocked': false,  'employer': 'slate' },
     *   { 'name': 'pebbles', 'blocked': true, 'employer': 'na' }
     * ];
     *
     * // using "_.pluck" callback shorthand
     * _.pluck(_.rest(characters, 'blocked'), 'name');
     * // => ['fred', 'pebbles']
     *
     * // using "_.where" callback shorthand
     * _.rest(characters, { 'employer': 'slate' });
     * // => [{ 'name': 'pebbles', 'blocked': true, 'employer': 'na' }]
     */
    function rest(array, callback, thisArg) {
      if (typeof callback != 'number' && callback != null) {
        var n = 0,
            index = -1,
            length = array ? array.length : 0;

        callback = lodash.createCallback(callback, thisArg, 3);
        while (++index < length && callback(array[index], index, array)) {
          n++;
        }
      } else {
        n = (callback == null || thisArg) ? 1 : nativeMax(0, callback);
      }
      return slice(array, n);
    }

    /**
     * Uses a binary search to determine the smallest index at which a value
     * should be inserted into a given sorted array in order to maintain the sort
     * order of the array. If a callback is provided it will be executed for
     * `value` and each element of `array` to compute their sort ranking. The
     * callback is bound to `thisArg` and invoked with one argument; (value).
     *
     * If a property name is provided for `callback` the created "_.pluck" style
     * callback will return the property value of the given element.
     *
     * If an object is provided for `callback` the created "_.where" style callback
     * will return `true` for elements that have the properties of the given object,
     * else `false`.
     *
     * @static
     * @memberOf _
     * @category Arrays
     * @param {Array} array The array to inspect.
     * @param {*} value The value to evaluate.
     * @param {Function|Object|string} [callback=identity] The function called
     *  per iteration. If a property name or object is provided it will be used
     *  to create a "_.pluck" or "_.where" style callback, respectively.
     * @param {*} [thisArg] The `this` binding of `callback`.
     * @returns {number} Returns the index at which `value` should be inserted
     *  into `array`.
     * @example
     *
     * _.sortedIndex([20, 30, 50], 40);
     * // => 2
     *
     * // using "_.pluck" callback shorthand
     * _.sortedIndex([{ 'x': 20 }, { 'x': 30 }, { 'x': 50 }], { 'x': 40 }, 'x');
     * // => 2
     *
     * var dict = {
     *   'wordToNumber': { 'twenty': 20, 'thirty': 30, 'fourty': 40, 'fifty': 50 }
     * };
     *
     * _.sortedIndex(['twenty', 'thirty', 'fifty'], 'fourty', function(word) {
     *   return dict.wordToNumber[word];
     * });
     * // => 2
     *
     * _.sortedIndex(['twenty', 'thirty', 'fifty'], 'fourty', function(word) {
     *   return this.wordToNumber[word];
     * }, dict);
     * // => 2
     */
    function sortedIndex(array, value, callback, thisArg) {
      var low = 0,
          high = array ? array.length : low;

      // explicitly reference `identity` for better inlining in Firefox
      callback = callback ? lodash.createCallback(callback, thisArg, 1) : identity;
      value = callback(value);

      while (low < high) {
        var mid = (low + high) >>> 1;
        (callback(array[mid]) < value)
          ? low = mid + 1
          : high = mid;
      }
      return low;
    }

    /**
     * Creates an array of unique values, in order, of the provided arrays using
     * strict equality for comparisons, i.e. `===`.
     *
     * @static
     * @memberOf _
     * @category Arrays
     * @param {...Array} [array] The arrays to inspect.
     * @returns {Array} Returns an array of combined values.
     * @example
     *
     * _.union([1, 2, 3], [5, 2, 1, 4], [2, 1]);
     * // => [1, 2, 3, 5, 4]
     */
    function union() {
      return baseUniq(baseFlatten(arguments, true, true));
    }

    /**
     * Creates a duplicate-value-free version of an array using strict equality
     * for comparisons, i.e. `===`. If the array is sorted, providing
     * `true` for `isSorted` will use a faster algorithm. If a callback is provided
     * each element of `array` is passed through the callback before uniqueness
     * is computed. The callback is bound to `thisArg` and invoked with three
     * arguments; (value, index, array).
     *
     * If a property name is provided for `callback` the created "_.pluck" style
     * callback will return the property value of the given element.
     *
     * If an object is provided for `callback` the created "_.where" style callback
     * will return `true` for elements that have the properties of the given object,
     * else `false`.
     *
     * @static
     * @memberOf _
     * @alias unique
     * @category Arrays
     * @param {Array} array The array to process.
     * @param {boolean} [isSorted=false] A flag to indicate that `array` is sorted.
     * @param {Function|Object|string} [callback=identity] The function called
     *  per iteration. If a property name or object is provided it will be used
     *  to create a "_.pluck" or "_.where" style callback, respectively.
     * @param {*} [thisArg] The `this` binding of `callback`.
     * @returns {Array} Returns a duplicate-value-free array.
     * @example
     *
     * _.uniq([1, 2, 1, 3, 1]);
     * // => [1, 2, 3]
     *
     * _.uniq([1, 1, 2, 2, 3], true);
     * // => [1, 2, 3]
     *
     * _.uniq(['A', 'b', 'C', 'a', 'B', 'c'], function(letter) { return letter.toLowerCase(); });
     * // => ['A', 'b', 'C']
     *
     * _.uniq([1, 2.5, 3, 1.5, 2, 3.5], function(num) { return this.floor(num); }, Math);
     * // => [1, 2.5, 3]
     *
     * // using "_.pluck" callback shorthand
     * _.uniq([{ 'x': 1 }, { 'x': 2 }, { 'x': 1 }], 'x');
     * // => [{ 'x': 1 }, { 'x': 2 }]
     */
    function uniq(array, isSorted, callback, thisArg) {
      // juggle arguments
      if (typeof isSorted != 'boolean' && isSorted != null) {
        thisArg = callback;
        callback = (typeof isSorted != 'function' && thisArg && thisArg[isSorted] === array) ? null : isSorted;
        isSorted = false;
      }
      if (callback != null) {
        callback = lodash.createCallback(callback, thisArg, 3);
      }
      return baseUniq(array, isSorted, callback);
    }

    /**
     * Creates an array excluding all provided values using strict equality for
     * comparisons, i.e. `===`.
     *
     * @static
     * @memberOf _
     * @category Arrays
     * @param {Array} array The array to filter.
     * @param {...*} [value] The values to exclude.
     * @returns {Array} Returns a new array of filtered values.
     * @example
     *
     * _.without([1, 2, 1, 0, 3, 1, 4], 0, 1);
     * // => [2, 3, 4]
     */
    function without(array) {
      return baseDifference(array, slice(arguments, 1));
    }

    /**
     * Creates an array that is the symmetric difference of the provided arrays.
     * See http://en.wikipedia.org/wiki/Symmetric_difference.
     *
     * @static
     * @memberOf _
     * @category Arrays
     * @param {...Array} [array] The arrays to inspect.
     * @returns {Array} Returns an array of values.
     * @example
     *
     * _.xor([1, 2, 3], [5, 2, 1, 4]);
     * // => [3, 5, 4]
     *
     * _.xor([1, 2, 5], [2, 3, 5], [3, 4, 5]);
     * // => [1, 4, 5]
     */
    function xor() {
      var index = -1,
          length = arguments.length;

      while (++index < length) {
        var array = arguments[index];
        if (isArray(array) || isArguments(array)) {
          var result = result
            ? baseUniq(baseDifference(result, array).concat(baseDifference(array, result)))
            : array;
        }
      }
      return result || [];
    }

    /**
     * Creates an array of grouped elements, the first of which contains the first
     * elements of the given arrays, the second of which contains the second
     * elements of the given arrays, and so on.
     *
     * @static
     * @memberOf _
     * @alias unzip
     * @category Arrays
     * @param {...Array} [array] Arrays to process.
     * @returns {Array} Returns a new array of grouped elements.
     * @example
     *
     * _.zip(['fred', 'barney'], [30, 40], [true, false]);
     * // => [['fred', 30, true], ['barney', 40, false]]
     */
    function zip() {
      var array = arguments.length > 1 ? arguments : arguments[0],
          index = -1,
          length = array ? max(pluck(array, 'length')) : 0,
          result = Array(length < 0 ? 0 : length);

      while (++index < length) {
        result[index] = pluck(array, index);
      }
      return result;
    }

    /**
     * Creates an object composed from arrays of `keys` and `values`. Provide
     * either a single two dimensional array, i.e. `[[key1, value1], [key2, value2]]`
     * or two arrays, one of `keys` and one of corresponding `values`.
     *
     * @static
     * @memberOf _
     * @alias object
     * @category Arrays
     * @param {Array} keys The array of keys.
     * @param {Array} [values=[]] The array of values.
     * @returns {Object} Returns an object composed of the given keys and
     *  corresponding values.
     * @example
     *
     * _.zipObject(['fred', 'barney'], [30, 40]);
     * // => { 'fred': 30, 'barney': 40 }
     */
    function zipObject(keys, values) {
      var index = -1,
          length = keys ? keys.length : 0,
          result = {};

      if (!values && length && !isArray(keys[0])) {
        values = [];
      }
      while (++index < length) {
        var key = keys[index];
        if (values) {
          result[key] = values[index];
        } else if (key) {
          result[key[0]] = key[1];
        }
      }
      return result;
    }

    /*--------------------------------------------------------------------------*/

    /**
     * Creates a function that executes `func`, with  the `this` binding and
     * arguments of the created function, only after being called `n` times.
     *
     * @static
     * @memberOf _
     * @category Functions
     * @param {number} n The number of times the function must be called before
     *  `func` is executed.
     * @param {Function} func The function to restrict.
     * @returns {Function} Returns the new restricted function.
     * @example
     *
     * var saves = ['profile', 'settings'];
     *
     * var done = _.after(saves.length, function() {
     *   console.log('Done saving!');
     * });
     *
     * _.forEach(saves, function(type) {
     *   asyncSave({ 'type': type, 'complete': done });
     * });
     * // => logs 'Done saving!', after all saves have completed
     */
    function after(n, func) {
      if (!isFunction(func)) {
        throw new TypeError;
      }
      return function() {
        if (--n < 1) {
          return func.apply(this, arguments);
        }
      };
    }

    /**
     * Creates a function that, when called, invokes `func` with the `this`
     * binding of `thisArg` and prepends any additional `bind` arguments to those
     * provided to the bound function.
     *
     * @static
     * @memberOf _
     * @category Functions
     * @param {Function} func The function to bind.
     * @param {*} [thisArg] The `this` binding of `func`.
     * @param {...*} [arg] Arguments to be partially applied.
     * @returns {Function} Returns the new bound function.
     * @example
     *
     * var func = function(greeting) {
     *   return greeting + ' ' + this.name;
     * };
     *
     * func = _.bind(func, { 'name': 'fred' }, 'hi');
     * func();
     * // => 'hi fred'
     */
    function bind(func, thisArg) {
      return arguments.length > 2
        ? createWrapper(func, 17, slice(arguments, 2), null, thisArg)
        : createWrapper(func, 1, null, null, thisArg);
    }

    /**
     * Binds methods of an object to the object itself, overwriting the existing
     * method. Method names may be specified as individual arguments or as arrays
     * of method names. If no method names are provided all the function properties
     * of `object` will be bound.
     *
     * @static
     * @memberOf _
     * @category Functions
     * @param {Object} object The object to bind and assign the bound methods to.
     * @param {...string} [methodName] The object method names to
     *  bind, specified as individual method names or arrays of method names.
     * @returns {Object} Returns `object`.
     * @example
     *
     * var view = {
     *   'label': 'docs',
     *   'onClick': function() { console.log('clicked ' + this.label); }
     * };
     *
     * _.bindAll(view);
     * jQuery('#docs').on('click', view.onClick);
     * // => logs 'clicked docs', when the button is clicked
     */
    function bindAll(object) {
      var funcs = arguments.length > 1 ? baseFlatten(arguments, true, false, 1) : functions(object),
          index = -1,
          length = funcs.length;

      while (++index < length) {
        var key = funcs[index];
        object[key] = createWrapper(object[key], 1, null, null, object);
      }
      return object;
    }

    /**
     * Creates a function that, when called, invokes the method at `object[key]`
     * and prepends any additional `bindKey` arguments to those provided to the bound
     * function. This method differs from `_.bind` by allowing bound functions to
     * reference methods that will be redefined or don't yet exist.
     * See http://michaux.ca/articles/lazy-function-definition-pattern.
     *
     * @static
     * @memberOf _
     * @category Functions
     * @param {Object} object The object the method belongs to.
     * @param {string} key The key of the method.
     * @param {...*} [arg] Arguments to be partially applied.
     * @returns {Function} Returns the new bound function.
     * @example
     *
     * var object = {
     *   'name': 'fred',
     *   'greet': function(greeting) {
     *     return greeting + ' ' + this.name;
     *   }
     * };
     *
     * var func = _.bindKey(object, 'greet', 'hi');
     * func();
     * // => 'hi fred'
     *
     * object.greet = function(greeting) {
     *   return greeting + 'ya ' + this.name + '!';
     * };
     *
     * func();
     * // => 'hiya fred!'
     */
    function bindKey(object, key) {
      return arguments.length > 2
        ? createWrapper(key, 19, slice(arguments, 2), null, object)
        : createWrapper(key, 3, null, null, object);
    }

    /**
     * Creates a function that is the composition of the provided functions,
     * where each function consumes the return value of the function that follows.
     * For example, composing the functions `f()`, `g()`, and `h()` produces `f(g(h()))`.
     * Each function is executed with the `this` binding of the composed function.
     *
     * @static
     * @memberOf _
     * @category Functions
     * @param {...Function} [func] Functions to compose.
     * @returns {Function} Returns the new composed function.
     * @example
     *
     * var realNameMap = {
     *   'pebbles': 'penelope'
     * };
     *
     * var format = function(name) {
     *   name = realNameMap[name.toLowerCase()] || name;
     *   return name.charAt(0).toUpperCase() + name.slice(1).toLowerCase();
     * };
     *
     * var greet = function(formatted) {
     *   return 'Hiya ' + formatted + '!';
     * };
     *
     * var welcome = _.compose(greet, format);
     * welcome('pebbles');
     * // => 'Hiya Penelope!'
     */
    function compose() {
      var funcs = arguments,
          length = funcs.length;

      while (length--) {
        if (!isFunction(funcs[length])) {
          throw new TypeError;
        }
      }
      return function() {
        var args = arguments,
            length = funcs.length;

        while (length--) {
          args = [funcs[length].apply(this, args)];
        }
        return args[0];
      };
    }

    /**
     * Creates a function which accepts one or more arguments of `func` that when
     * invoked either executes `func` returning its result, if all `func` arguments
     * have been provided, or returns a function that accepts one or more of the
     * remaining `func` arguments, and so on. The arity of `func` can be specified
     * if `func.length` is not sufficient.
     *
     * @static
     * @memberOf _
     * @category Functions
     * @param {Function} func The function to curry.
     * @param {number} [arity=func.length] The arity of `func`.
     * @returns {Function} Returns the new curried function.
     * @example
     *
     * var curried = _.curry(function(a, b, c) {
     *   console.log(a + b + c);
     * });
     *
     * curried(1)(2)(3);
     * // => 6
     *
     * curried(1, 2)(3);
     * // => 6
     *
     * curried(1, 2, 3);
     * // => 6
     */
    function curry(func, arity) {
      arity = typeof arity == 'number' ? arity : (+arity || func.length);
      return createWrapper(func, 4, null, null, null, arity);
    }

    /**
     * Creates a function that will delay the execution of `func` until after
     * `wait` milliseconds have elapsed since the last time it was invoked.
     * Provide an options object to indicate that `func` should be invoked on
     * the leading and/or trailing edge of the `wait` timeout. Subsequent calls
     * to the debounced function will return the result of the last `func` call.
     *
     * Note: If `leading` and `trailing` options are `true` `func` will be called
     * on the trailing edge of the timeout only if the the debounced function is
     * invoked more than once during the `wait` timeout.
     *
     * @static
     * @memberOf _
     * @category Functions
     * @param {Function} func The function to debounce.
     * @param {number} wait The number of milliseconds to delay.
     * @param {Object} [options] The options object.
     * @param {boolean} [options.leading=false] Specify execution on the leading edge of the timeout.
     * @param {number} [options.maxWait] The maximum time `func` is allowed to be delayed before it's called.
     * @param {boolean} [options.trailing=true] Specify execution on the trailing edge of the timeout.
     * @returns {Function} Returns the new debounced function.
     * @example
     *
     * // avoid costly calculations while the window size is in flux
     * var lazyLayout = _.debounce(calculateLayout, 150);
     * jQuery(window).on('resize', lazyLayout);
     *
     * // execute `sendMail` when the click event is fired, debouncing subsequent calls
     * jQuery('#postbox').on('click', _.debounce(sendMail, 300, {
     *   'leading': true,
     *   'trailing': false
     * });
     *
     * // ensure `batchLog` is executed once after 1 second of debounced calls
     * var source = new EventSource('/stream');
     * source.addEventListener('message', _.debounce(batchLog, 250, {
     *   'maxWait': 1000
     * }, false);
     */
    function debounce(func, wait, options) {
      var args,
          maxTimeoutId,
          result,
          stamp,
          thisArg,
          timeoutId,
          trailingCall,
          lastCalled = 0,
          maxWait = false,
          trailing = true;

      if (!isFunction(func)) {
        throw new TypeError;
      }
      wait = nativeMax(0, wait) || 0;
      if (options === true) {
        var leading = true;
        trailing = false;
      } else if (isObject(options)) {
        leading = options.leading;
        maxWait = 'maxWait' in options && (nativeMax(wait, options.maxWait) || 0);
        trailing = 'trailing' in options ? options.trailing : trailing;
      }
      var delayed = function() {
        var remaining = wait - (now() - stamp);
        if (remaining <= 0) {
          if (maxTimeoutId) {
            clearTimeout(maxTimeoutId);
          }
          var isCalled = trailingCall;
          maxTimeoutId = timeoutId = trailingCall = undefined;
          if (isCalled) {
            lastCalled = now();
            result = func.apply(thisArg, args);
            if (!timeoutId && !maxTimeoutId) {
              args = thisArg = null;
            }
          }
        } else {
          timeoutId = setTimeout(delayed, remaining);
        }
      };

      var maxDelayed = function() {
        if (timeoutId) {
          clearTimeout(timeoutId);
        }
        maxTimeoutId = timeoutId = trailingCall = undefined;
        if (trailing || (maxWait !== wait)) {
          lastCalled = now();
          result = func.apply(thisArg, args);
          if (!timeoutId && !maxTimeoutId) {
            args = thisArg = null;
          }
        }
      };

      return function() {
        args = arguments;
        stamp = now();
        thisArg = this;
        trailingCall = trailing && (timeoutId || !leading);

        if (maxWait === false) {
          var leadingCall = leading && !timeoutId;
        } else {
          if (!maxTimeoutId && !leading) {
            lastCalled = stamp;
          }
          var remaining = maxWait - (stamp - lastCalled),
              isCalled = remaining <= 0;

          if (isCalled) {
            if (maxTimeoutId) {
              maxTimeoutId = clearTimeout(maxTimeoutId);
            }
            lastCalled = stamp;
            result = func.apply(thisArg, args);
          }
          else if (!maxTimeoutId) {
            maxTimeoutId = setTimeout(maxDelayed, remaining);
          }
        }
        if (isCalled && timeoutId) {
          timeoutId = clearTimeout(timeoutId);
        }
        else if (!timeoutId && wait !== maxWait) {
          timeoutId = setTimeout(delayed, wait);
        }
        if (leadingCall) {
          isCalled = true;
          result = func.apply(thisArg, args);
        }
        if (isCalled && !timeoutId && !maxTimeoutId) {
          args = thisArg = null;
        }
        return result;
      };
    }

    /**
     * Defers executing the `func` function until the current call stack has cleared.
     * Additional arguments will be provided to `func` when it is invoked.
     *
     * @static
     * @memberOf _
     * @category Functions
     * @param {Function} func The function to defer.
     * @param {...*} [arg] Arguments to invoke the function with.
     * @returns {number} Returns the timer id.
     * @example
     *
     * _.defer(function(text) { console.log(text); }, 'deferred');
     * // logs 'deferred' after one or more milliseconds
     */
    function defer(func) {
      if (!isFunction(func)) {
        throw new TypeError;
      }
      var args = slice(arguments, 1);
      return setTimeout(function() { func.apply(undefined, args); }, 1);
    }

    /**
     * Executes the `func` function after `wait` milliseconds. Additional arguments
     * will be provided to `func` when it is invoked.
     *
     * @static
     * @memberOf _
     * @category Functions
     * @param {Function} func The function to delay.
     * @param {number} wait The number of milliseconds to delay execution.
     * @param {...*} [arg] Arguments to invoke the function with.
     * @returns {number} Returns the timer id.
     * @example
     *
     * _.delay(function(text) { console.log(text); }, 1000, 'later');
     * // => logs 'later' after one second
     */
    function delay(func, wait) {
      if (!isFunction(func)) {
        throw new TypeError;
      }
      var args = slice(arguments, 2);
      return setTimeout(function() { func.apply(undefined, args); }, wait);
    }

    /**
     * Creates a function that memoizes the result of `func`. If `resolver` is
     * provided it will be used to determine the cache key for storing the result
     * based on the arguments provided to the memoized function. By default, the
     * first argument provided to the memoized function is used as the cache key.
     * The `func` is executed with the `this` binding of the memoized function.
     * The result cache is exposed as the `cache` property on the memoized function.
     *
     * @static
     * @memberOf _
     * @category Functions
     * @param {Function} func The function to have its output memoized.
     * @param {Function} [resolver] A function used to resolve the cache key.
     * @returns {Function} Returns the new memoizing function.
     * @example
     *
     * var fibonacci = _.memoize(function(n) {
     *   return n < 2 ? n : fibonacci(n - 1) + fibonacci(n - 2);
     * });
     *
     * fibonacci(9)
     * // => 34
     *
     * var data = {
     *   'fred': { 'name': 'fred', 'age': 40 },
     *   'pebbles': { 'name': 'pebbles', 'age': 1 }
     * };
     *
     * // modifying the result cache
     * var get = _.memoize(function(name) { return data[name]; }, _.identity);
     * get('pebbles');
     * // => { 'name': 'pebbles', 'age': 1 }
     *
     * get.cache.pebbles.name = 'penelope';
     * get('pebbles');
     * // => { 'name': 'penelope', 'age': 1 }
     */
    function memoize(func, resolver) {
      if (!isFunction(func)) {
        throw new TypeError;
      }
      var memoized = function() {
        var cache = memoized.cache,
            key = resolver ? resolver.apply(this, arguments) : keyPrefix + arguments[0];

        return hasOwnProperty.call(cache, key)
          ? cache[key]
          : (cache[key] = func.apply(this, arguments));
      }
      memoized.cache = {};
      return memoized;
    }

    /**
     * Creates a function that is restricted to execute `func` once. Repeat calls to
     * the function will return the value of the first call. The `func` is executed
     * with the `this` binding of the created function.
     *
     * @static
     * @memberOf _
     * @category Functions
     * @param {Function} func The function to restrict.
     * @returns {Function} Returns the new restricted function.
     * @example
     *
     * var initialize = _.once(createApplication);
     * initialize();
     * initialize();
     * // `initialize` executes `createApplication` once
     */
    function once(func) {
      var ran,
          result;

      if (!isFunction(func)) {
        throw new TypeError;
      }
      return function() {
        if (ran) {
          return result;
        }
        ran = true;
        result = func.apply(this, arguments);

        // clear the `func` variable so the function may be garbage collected
        func = null;
        return result;
      };
    }

    /**
     * Creates a function that, when called, invokes `func` with any additional
     * `partial` arguments prepended to those provided to the new function. This
     * method is similar to `_.bind` except it does **not** alter the `this` binding.
     *
     * @static
     * @memberOf _
     * @category Functions
     * @param {Function} func The function to partially apply arguments to.
     * @param {...*} [arg] Arguments to be partially applied.
     * @returns {Function} Returns the new partially applied function.
     * @example
     *
     * var greet = function(greeting, name) { return greeting + ' ' + name; };
     * var hi = _.partial(greet, 'hi');
     * hi('fred');
     * // => 'hi fred'
     */
    function partial(func) {
      return createWrapper(func, 16, slice(arguments, 1));
    }

    /**
     * This method is like `_.partial` except that `partial` arguments are
     * appended to those provided to the new function.
     *
     * @static
     * @memberOf _
     * @category Functions
     * @param {Function} func The function to partially apply arguments to.
     * @param {...*} [arg] Arguments to be partially applied.
     * @returns {Function} Returns the new partially applied function.
     * @example
     *
     * var defaultsDeep = _.partialRight(_.merge, _.defaults);
     *
     * var options = {
     *   'variable': 'data',
     *   'imports': { 'jq': $ }
     * };
     *
     * defaultsDeep(options, _.templateSettings);
     *
     * options.variable
     * // => 'data'
     *
     * options.imports
     * // => { '_': _, 'jq': $ }
     */
    function partialRight(func) {
      return createWrapper(func, 32, null, slice(arguments, 1));
    }

    /**
     * Creates a function that, when executed, will only call the `func` function
     * at most once per every `wait` milliseconds. Provide an options object to
     * indicate that `func` should be invoked on the leading and/or trailing edge
     * of the `wait` timeout. Subsequent calls to the throttled function will
     * return the result of the last `func` call.
     *
     * Note: If `leading` and `trailing` options are `true` `func` will be called
     * on the trailing edge of the timeout only if the the throttled function is
     * invoked more than once during the `wait` timeout.
     *
     * @static
     * @memberOf _
     * @category Functions
     * @param {Function} func The function to throttle.
     * @param {number} wait The number of milliseconds to throttle executions to.
     * @param {Object} [options] The options object.
     * @param {boolean} [options.leading=true] Specify execution on the leading edge of the timeout.
     * @param {boolean} [options.trailing=true] Specify execution on the trailing edge of the timeout.
     * @returns {Function} Returns the new throttled function.
     * @example
     *
     * // avoid excessively updating the position while scrolling
     * var throttled = _.throttle(updatePosition, 100);
     * jQuery(window).on('scroll', throttled);
     *
     * // execute `renewToken` when the click event is fired, but not more than once every 5 minutes
     * jQuery('.interactive').on('click', _.throttle(renewToken, 300000, {
     *   'trailing': false
     * }));
     */
    function throttle(func, wait, options) {
      var leading = true,
          trailing = true;

      if (!isFunction(func)) {
        throw new TypeError;
      }
      if (options === false) {
        leading = false;
      } else if (isObject(options)) {
        leading = 'leading' in options ? options.leading : leading;
        trailing = 'trailing' in options ? options.trailing : trailing;
      }
      debounceOptions.leading = leading;
      debounceOptions.maxWait = wait;
      debounceOptions.trailing = trailing;

      return debounce(func, wait, debounceOptions);
    }

    /**
     * Creates a function that provides `value` to the wrapper function as its
     * first argument. Additional arguments provided to the function are appended
     * to those provided to the wrapper function. The wrapper is executed with
     * the `this` binding of the created function.
     *
     * @static
     * @memberOf _
     * @category Functions
     * @param {*} value The value to wrap.
     * @param {Function} wrapper The wrapper function.
     * @returns {Function} Returns the new function.
     * @example
     *
     * var p = _.wrap(_.escape, function(func, text) {
     *   return '<p>' + func(text) + '</p>';
     * });
     *
     * p('Fred, Wilma, & Pebbles');
     * // => '<p>Fred, Wilma, &amp; Pebbles</p>'
     */
    function wrap(value, wrapper) {
      return createWrapper(wrapper, 16, [value]);
    }

    /*--------------------------------------------------------------------------*/

    /**
     * Creates a function that returns `value`.
     *
     * @static
     * @memberOf _
     * @category Utilities
     * @param {*} value The value to return from the new function.
     * @returns {Function} Returns the new function.
     * @example
     *
     * var object = { 'name': 'fred' };
     * var getter = _.constant(object);
     * getter() === object;
     * // => true
     */
    function constant(value) {
      return function() {
        return value;
      };
    }

    /**
     * Produces a callback bound to an optional `thisArg`. If `func` is a property
     * name the created callback will return the property value for a given element.
     * If `func` is an object the created callback will return `true` for elements
     * that contain the equivalent object properties, otherwise it will return `false`.
     *
     * @static
     * @memberOf _
     * @category Utilities
     * @param {*} [func=identity] The value to convert to a callback.
     * @param {*} [thisArg] The `this` binding of the created callback.
     * @param {number} [argCount] The number of arguments the callback accepts.
     * @returns {Function} Returns a callback function.
     * @example
     *
     * var characters = [
     *   { 'name': 'barney', 'age': 36 },
     *   { 'name': 'fred',   'age': 40 }
     * ];
     *
     * // wrap to create custom callback shorthands
     * _.createCallback = _.wrap(_.createCallback, function(func, callback, thisArg) {
     *   var match = /^(.+?)__([gl]t)(.+)$/.exec(callback);
     *   return !match ? func(callback, thisArg) : function(object) {
     *     return match[2] == 'gt' ? object[match[1]] > match[3] : object[match[1]] < match[3];
     *   };
     * });
     *
     * _.filter(characters, 'age__gt38');
     * // => [{ 'name': 'fred', 'age': 40 }]
     */
    function createCallback(func, thisArg, argCount) {
      var type = typeof func;
      if (func == null || type == 'function') {
        return baseCreateCallback(func, thisArg, argCount);
      }
      // handle "_.pluck" style callback shorthands
      if (type != 'object') {
        return property(func);
      }
      var props = keys(func),
          key = props[0],
          a = func[key];

      // handle "_.where" style callback shorthands
      if (props.length == 1 && a === a && !isObject(a)) {
        // fast path the common case of providing an object with a single
        // property containing a primitive value
        return function(object) {
          var b = object[key];
          return a === b && (a !== 0 || (1 / a == 1 / b));
        };
      }
      return function(object) {
        var length = props.length,
            result = false;

        while (length--) {
          if (!(result = baseIsEqual(object[props[length]], func[props[length]], null, true))) {
            break;
          }
        }
        return result;
      };
    }

    /**
     * Converts the characters `&`, `<`, `>`, `"`, and `'` in `string` to their
     * corresponding HTML entities.
     *
     * @static
     * @memberOf _
     * @category Utilities
     * @param {string} string The string to escape.
     * @returns {string} Returns the escaped string.
     * @example
     *
     * _.escape('Fred, Wilma, & Pebbles');
     * // => 'Fred, Wilma, &amp; Pebbles'
     */
    function escape(string) {
      return string == null ? '' : String(string).replace(reUnescapedHtml, escapeHtmlChar);
    }

    /**
     * This method returns the first argument provided to it.
     *
     * @static
     * @memberOf _
     * @category Utilities
     * @param {*} value Any value.
     * @returns {*} Returns `value`.
     * @example
     *
     * var object = { 'name': 'fred' };
     * _.identity(object) === object;
     * // => true
     */
    function identity(value) {
      return value;
    }

    /**
     * Adds function properties of a source object to the destination object.
     * If `object` is a function methods will be added to its prototype as well.
     *
     * @static
     * @memberOf _
     * @category Utilities
     * @param {Function|Object} [object=lodash] object The destination object.
     * @param {Object} source The object of functions to add.
     * @param {Object} [options] The options object.
     * @param {boolean} [options.chain=true] Specify whether the functions added are chainable.
     * @example
     *
     * function capitalize(string) {
     *   return string.charAt(0).toUpperCase() + string.slice(1).toLowerCase();
     * }
     *
     * _.mixin({ 'capitalize': capitalize });
     * _.capitalize('fred');
     * // => 'Fred'
     *
     * _('fred').capitalize().value();
     * // => 'Fred'
     *
     * _.mixin({ 'capitalize': capitalize }, { 'chain': false });
     * _('fred').capitalize();
     * // => 'Fred'
     */
    function mixin(object, source, options) {
      var chain = true,
          methodNames = source && functions(source);

      if (!source || (!options && !methodNames.length)) {
        if (options == null) {
          options = source;
        }
        ctor = lodashWrapper;
        source = object;
        object = lodash;
        methodNames = functions(source);
      }
      if (options === false) {
        chain = false;
      } else if (isObject(options) && 'chain' in options) {
        chain = options.chain;
      }
      var ctor = object,
          isFunc = isFunction(ctor);

      forEach(methodNames, function(methodName) {
        var func = object[methodName] = source[methodName];
        if (isFunc) {
          ctor.prototype[methodName] = function() {
            var chainAll = this.__chain__,
                value = this.__wrapped__,
                args = [value];

            push.apply(args, arguments);
            var result = func.apply(object, args);
            if (chain || chainAll) {
              if (value === result && isObject(result)) {
                return this;
              }
              result = new ctor(result);
              result.__chain__ = chainAll;
            }
            return result;
          };
        }
      });
    }

    /**
     * Reverts the '_' variable to its previous value and returns a reference to
     * the `lodash` function.
     *
     * @static
     * @memberOf _
     * @category Utilities
     * @returns {Function} Returns the `lodash` function.
     * @example
     *
     * var lodash = _.noConflict();
     */
    function noConflict() {
      context._ = oldDash;
      return this;
    }

    /**
     * A no-operation function.
     *
     * @static
     * @memberOf _
     * @category Utilities
     * @example
     *
     * var object = { 'name': 'fred' };
     * _.noop(object) === undefined;
     * // => true
     */
    function noop() {
      // no operation performed
    }

    /**
     * Gets the number of milliseconds that have elapsed since the Unix epoch
     * (1 January 1970 00:00:00 UTC).
     *
     * @static
     * @memberOf _
     * @category Utilities
     * @example
     *
     * var stamp = _.now();
     * _.defer(function() { console.log(_.now() - stamp); });
     * // => logs the number of milliseconds it took for the deferred function to be called
     */
    var now = isNative(now = Date.now) && now || function() {
      return new Date().getTime();
    };

    /**
     * Converts the given value into an integer of the specified radix.
     * If `radix` is `undefined` or `0` a `radix` of `10` is used unless the
     * `value` is a hexadecimal, in which case a `radix` of `16` is used.
     *
     * Note: This method avoids differences in native ES3 and ES5 `parseInt`
     * implementations. See http://es5.github.io/#E.
     *
     * @static
     * @memberOf _
     * @category Utilities
     * @param {string} value The value to parse.
     * @param {number} [radix] The radix used to interpret the value to parse.
     * @returns {number} Returns the new integer value.
     * @example
     *
     * _.parseInt('08');
     * // => 8
     */
    var parseInt = nativeParseInt(whitespace + '08') == 8 ? nativeParseInt : function(value, radix) {
      // Firefox < 21 and Opera < 15 follow the ES3 specified implementation of `parseInt`
      return nativeParseInt(isString(value) ? value.replace(reLeadingSpacesAndZeros, '') : value, radix || 0);
    };

    /**
     * Creates a "_.pluck" style function, which returns the `key` value of a
     * given object.
     *
     * @static
     * @memberOf _
     * @category Utilities
     * @param {string} key The name of the property to retrieve.
     * @returns {Function} Returns the new function.
     * @example
     *
     * var characters = [
     *   { 'name': 'fred',   'age': 40 },
     *   { 'name': 'barney', 'age': 36 }
     * ];
     *
     * var getName = _.property('name');
     *
     * _.map(characters, getName);
     * // => ['barney', 'fred']
     *
     * _.sortBy(characters, getName);
     * // => [{ 'name': 'barney', 'age': 36 }, { 'name': 'fred',   'age': 40 }]
     */
    function property(key) {
      return function(object) {
        return object[key];
      };
    }

    /**
     * Produces a random number between `min` and `max` (inclusive). If only one
     * argument is provided a number between `0` and the given number will be
     * returned. If `floating` is truey or either `min` or `max` are floats a
     * floating-point number will be returned instead of an integer.
     *
     * @static
     * @memberOf _
     * @category Utilities
     * @param {number} [min=0] The minimum possible value.
     * @param {number} [max=1] The maximum possible value.
     * @param {boolean} [floating=false] Specify returning a floating-point number.
     * @returns {number} Returns a random number.
     * @example
     *
     * _.random(0, 5);
     * // => an integer between 0 and 5
     *
     * _.random(5);
     * // => also an integer between 0 and 5
     *
     * _.random(5, true);
     * // => a floating-point number between 0 and 5
     *
     * _.random(1.2, 5.2);
     * // => a floating-point number between 1.2 and 5.2
     */
    function random(min, max, floating) {
      var noMin = min == null,
          noMax = max == null;

      if (floating == null) {
        if (typeof min == 'boolean' && noMax) {
          floating = min;
          min = 1;
        }
        else if (!noMax && typeof max == 'boolean') {
          floating = max;
          noMax = true;
        }
      }
      if (noMin && noMax) {
        max = 1;
      }
      min = +min || 0;
      if (noMax) {
        max = min;
        min = 0;
      } else {
        max = +max || 0;
      }
      if (floating || min % 1 || max % 1) {
        var rand = nativeRandom();
        return nativeMin(min + (rand * (max - min + parseFloat('1e-' + ((rand +'').length - 1)))), max);
      }
      return baseRandom(min, max);
    }

    /**
     * Resolves the value of property `key` on `object`. If `key` is a function
     * it will be invoked with the `this` binding of `object` and its result returned,
     * else the property value is returned. If `object` is falsey then `undefined`
     * is returned.
     *
     * @static
     * @memberOf _
     * @category Utilities
     * @param {Object} object The object to inspect.
     * @param {string} key The name of the property to resolve.
     * @returns {*} Returns the resolved value.
     * @example
     *
     * var object = {
     *   'cheese': 'crumpets',
     *   'stuff': function() {
     *     return 'nonsense';
     *   }
     * };
     *
     * _.result(object, 'cheese');
     * // => 'crumpets'
     *
     * _.result(object, 'stuff');
     * // => 'nonsense'
     */
    function result(object, key) {
      if (object) {
        var value = object[key];
        return isFunction(value) ? object[key]() : value;
      }
    }

    /**
     * A micro-templating method that handles arbitrary delimiters, preserves
     * whitespace, and correctly escapes quotes within interpolated code.
     *
     * Note: In the development build, `_.template` utilizes sourceURLs for easier
     * debugging. See http://www.html5rocks.com/en/tutorials/developertools/sourcemaps/#toc-sourceurl
     *
     * For more information on precompiling templates see:
     * https://lodash.com/custom-builds
     *
     * For more information on Chrome extension sandboxes see:
     * http://developer.chrome.com/stable/extensions/sandboxingEval.html
     *
     * @static
     * @memberOf _
     * @category Utilities
     * @param {string} text The template text.
     * @param {Object} data The data object used to populate the text.
     * @param {Object} [options] The options object.
     * @param {RegExp} [options.escape] The "escape" delimiter.
     * @param {RegExp} [options.evaluate] The "evaluate" delimiter.
     * @param {Object} [options.imports] An object to import into the template as local variables.
     * @param {RegExp} [options.interpolate] The "interpolate" delimiter.
     * @param {string} [sourceURL] The sourceURL of the template's compiled source.
     * @param {string} [variable] The data object variable name.
     * @returns {Function|string} Returns a compiled function when no `data` object
     *  is given, else it returns the interpolated text.
     * @example
     *
     * // using the "interpolate" delimiter to create a compiled template
     * var compiled = _.template('hello <%= name %>');
     * compiled({ 'name': 'fred' });
     * // => 'hello fred'
     *
     * // using the "escape" delimiter to escape HTML in data property values
     * _.template('<b><%- value %></b>', { 'value': '<script>' });
     * // => '<b>&lt;script&gt;</b>'
     *
     * // using the "evaluate" delimiter to generate HTML
     * var list = '<% _.forEach(people, function(name) { %><li><%- name %></li><% }); %>';
     * _.template(list, { 'people': ['fred', 'barney'] });
     * // => '<li>fred</li><li>barney</li>'
     *
     * // using the ES6 delimiter as an alternative to the default "interpolate" delimiter
     * _.template('hello ${ name }', { 'name': 'pebbles' });
     * // => 'hello pebbles'
     *
     * // using the internal `print` function in "evaluate" delimiters
     * _.template('<% print("hello " + name); %>!', { 'name': 'barney' });
     * // => 'hello barney!'
     *
     * // using a custom template delimiters
     * _.templateSettings = {
     *   'interpolate': /{{([\s\S]+?)}}/g
     * };
     *
     * _.template('hello {{ name }}!', { 'name': 'mustache' });
     * // => 'hello mustache!'
     *
     * // using the `imports` option to import jQuery
     * var list = '<% jq.each(people, function(name) { %><li><%- name %></li><% }); %>';
     * _.template(list, { 'people': ['fred', 'barney'] }, { 'imports': { 'jq': jQuery } });
     * // => '<li>fred</li><li>barney</li>'
     *
     * // using the `sourceURL` option to specify a custom sourceURL for the template
     * var compiled = _.template('hello <%= name %>', null, { 'sourceURL': '/basic/greeting.jst' });
     * compiled(data);
     * // => find the source of "greeting.jst" under the Sources tab or Resources panel of the web inspector
     *
     * // using the `variable` option to ensure a with-statement isn't used in the compiled template
     * var compiled = _.template('hi <%= data.name %>!', null, { 'variable': 'data' });
     * compiled.source;
     * // => function(data) {
     *   var __t, __p = '', __e = _.escape;
     *   __p += 'hi ' + ((__t = ( data.name )) == null ? '' : __t) + '!';
     *   return __p;
     * }
     *
     * // using the `source` property to inline compiled templates for meaningful
     * // line numbers in error messages and a stack trace
     * fs.writeFileSync(path.join(cwd, 'jst.js'), '\
     *   var JST = {\
     *     "main": ' + _.template(mainText).source + '\
     *   };\
     * ');
     */
    function template(text, data, options) {
      // based on John Resig's `tmpl` implementation
      // http://ejohn.org/blog/javascript-micro-templating/
      // and Laura Doktorova's doT.js
      // https://github.com/olado/doT
      var settings = lodash.templateSettings;
      text = String(text || '');

      // avoid missing dependencies when `iteratorTemplate` is not defined
      options = defaults({}, options, settings);

      var imports = defaults({}, options.imports, settings.imports),
          importsKeys = keys(imports),
          importsValues = values(imports);

      var isEvaluating,
          index = 0,
          interpolate = options.interpolate || reNoMatch,
          source = "__p += '";

      // compile the regexp to match each delimiter
      var reDelimiters = RegExp(
        (options.escape || reNoMatch).source + '|' +
        interpolate.source + '|' +
        (interpolate === reInterpolate ? reEsTemplate : reNoMatch).source + '|' +
        (options.evaluate || reNoMatch).source + '|$'
      , 'g');

      text.replace(reDelimiters, function(match, escapeValue, interpolateValue, esTemplateValue, evaluateValue, offset) {
        interpolateValue || (interpolateValue = esTemplateValue);

        // escape characters that cannot be included in string literals
        source += text.slice(index, offset).replace(reUnescapedString, escapeStringChar);

        // replace delimiters with snippets
        if (escapeValue) {
          source += "' +\n__e(" + escapeValue + ") +\n'";
        }
        if (evaluateValue) {
          isEvaluating = true;
          source += "';\n" + evaluateValue + ";\n__p += '";
        }
        if (interpolateValue) {
          source += "' +\n((__t = (" + interpolateValue + ")) == null ? '' : __t) +\n'";
        }
        index = offset + match.length;

        // the JS engine embedded in Adobe products requires returning the `match`
        // string in order to produce the correct `offset` value
        return match;
      });

      source += "';\n";

      // if `variable` is not specified, wrap a with-statement around the generated
      // code to add the data object to the top of the scope chain
      var variable = options.variable,
          hasVariable = variable;

      if (!hasVariable) {
        variable = 'obj';
        source = 'with (' + variable + ') {\n' + source + '\n}\n';
      }
      // cleanup code by stripping empty strings
      source = (isEvaluating ? source.replace(reEmptyStringLeading, '') : source)
        .replace(reEmptyStringMiddle, '$1')
        .replace(reEmptyStringTrailing, '$1;');

      // frame code as the function body
      source = 'function(' + variable + ') {\n' +
        (hasVariable ? '' : variable + ' || (' + variable + ' = {});\n') +
        "var __t, __p = '', __e = _.escape" +
        (isEvaluating
          ? ', __j = Array.prototype.join;\n' +
            "function print() { __p += __j.call(arguments, '') }\n"
          : ';\n'
        ) +
        source +
        'return __p\n}';

      // Use a sourceURL for easier debugging.
      // http://www.html5rocks.com/en/tutorials/developertools/sourcemaps/#toc-sourceurl
      var sourceURL = '\n/*\n//# sourceURL=' + (options.sourceURL || '/lodash/template/source[' + (templateCounter++) + ']') + '\n*/';

      try {
        var result = Function(importsKeys, 'return ' + source + sourceURL).apply(undefined, importsValues);
      } catch(e) {
        e.source = source;
        throw e;
      }
      if (data) {
        return result(data);
      }
      // provide the compiled function's source by its `toString` method, in
      // supported environments, or the `source` property as a convenience for
      // inlining compiled templates during the build process
      result.source = source;
      return result;
    }

    /**
     * Executes the callback `n` times, returning an array of the results
     * of each callback execution. The callback is bound to `thisArg` and invoked
     * with one argument; (index).
     *
     * @static
     * @memberOf _
     * @category Utilities
     * @param {number} n The number of times to execute the callback.
     * @param {Function} callback The function called per iteration.
     * @param {*} [thisArg] The `this` binding of `callback`.
     * @returns {Array} Returns an array of the results of each `callback` execution.
     * @example
     *
     * var diceRolls = _.times(3, _.partial(_.random, 1, 6));
     * // => [3, 6, 4]
     *
     * _.times(3, function(n) { mage.castSpell(n); });
     * // => calls `mage.castSpell(n)` three times, passing `n` of `0`, `1`, and `2` respectively
     *
     * _.times(3, function(n) { this.cast(n); }, mage);
     * // => also calls `mage.castSpell(n)` three times
     */
    function times(n, callback, thisArg) {
      n = (n = +n) > -1 ? n : 0;
      var index = -1,
          result = Array(n);

      callback = baseCreateCallback(callback, thisArg, 1);
      while (++index < n) {
        result[index] = callback(index);
      }
      return result;
    }

    /**
     * The inverse of `_.escape` this method converts the HTML entities
     * `&amp;`, `&lt;`, `&gt;`, `&quot;`, and `&#39;` in `string` to their
     * corresponding characters.
     *
     * @static
     * @memberOf _
     * @category Utilities
     * @param {string} string The string to unescape.
     * @returns {string} Returns the unescaped string.
     * @example
     *
     * _.unescape('Fred, Barney &amp; Pebbles');
     * // => 'Fred, Barney & Pebbles'
     */
    function unescape(string) {
      return string == null ? '' : String(string).replace(reEscapedHtml, unescapeHtmlChar);
    }

    /**
     * Generates a unique ID. If `prefix` is provided the ID will be appended to it.
     *
     * @static
     * @memberOf _
     * @category Utilities
     * @param {string} [prefix] The value to prefix the ID with.
     * @returns {string} Returns the unique ID.
     * @example
     *
     * _.uniqueId('contact_');
     * // => 'contact_104'
     *
     * _.uniqueId();
     * // => '105'
     */
    function uniqueId(prefix) {
      var id = ++idCounter;
      return String(prefix == null ? '' : prefix) + id;
    }

    /*--------------------------------------------------------------------------*/

    /**
     * Creates a `lodash` object that wraps the given value with explicit
     * method chaining enabled.
     *
     * @static
     * @memberOf _
     * @category Chaining
     * @param {*} value The value to wrap.
     * @returns {Object} Returns the wrapper object.
     * @example
     *
     * var characters = [
     *   { 'name': 'barney',  'age': 36 },
     *   { 'name': 'fred',    'age': 40 },
     *   { 'name': 'pebbles', 'age': 1 }
     * ];
     *
     * var youngest = _.chain(characters)
     *     .sortBy('age')
     *     .map(function(chr) { return chr.name + ' is ' + chr.age; })
     *     .first()
     *     .value();
     * // => 'pebbles is 1'
     */
    function chain(value) {
      value = new lodashWrapper(value);
      value.__chain__ = true;
      return value;
    }

    /**
     * Invokes `interceptor` with the `value` as the first argument and then
     * returns `value`. The purpose of this method is to "tap into" a method
     * chain in order to perform operations on intermediate results within
     * the chain.
     *
     * @static
     * @memberOf _
     * @category Chaining
     * @param {*} value The value to provide to `interceptor`.
     * @param {Function} interceptor The function to invoke.
     * @returns {*} Returns `value`.
     * @example
     *
     * _([1, 2, 3, 4])
     *  .tap(function(array) { array.pop(); })
     *  .reverse()
     *  .value();
     * // => [3, 2, 1]
     */
    function tap(value, interceptor) {
      interceptor(value);
      return value;
    }

    /**
     * Enables explicit method chaining on the wrapper object.
     *
     * @name chain
     * @memberOf _
     * @category Chaining
     * @returns {*} Returns the wrapper object.
     * @example
     *
     * var characters = [
     *   { 'name': 'barney', 'age': 36 },
     *   { 'name': 'fred',   'age': 40 }
     * ];
     *
     * // without explicit chaining
     * _(characters).first();
     * // => { 'name': 'barney', 'age': 36 }
     *
     * // with explicit chaining
     * _(characters).chain()
     *   .first()
     *   .pick('age')
     *   .value();
     * // => { 'age': 36 }
     */
    function wrapperChain() {
      this.__chain__ = true;
      return this;
    }

    /**
     * Produces the `toString` result of the wrapped value.
     *
     * @name toString
     * @memberOf _
     * @category Chaining
     * @returns {string} Returns the string result.
     * @example
     *
     * _([1, 2, 3]).toString();
     * // => '1,2,3'
     */
    function wrapperToString() {
      return String(this.__wrapped__);
    }

    /**
     * Extracts the wrapped value.
     *
     * @name valueOf
     * @memberOf _
     * @alias value
     * @category Chaining
     * @returns {*} Returns the wrapped value.
     * @example
     *
     * _([1, 2, 3]).valueOf();
     * // => [1, 2, 3]
     */
    function wrapperValueOf() {
      return this.__wrapped__;
    }

    /*--------------------------------------------------------------------------*/

    // add functions that return wrapped values when chaining
    lodash.after = after;
    lodash.assign = assign;
    lodash.at = at;
    lodash.bind = bind;
    lodash.bindAll = bindAll;
    lodash.bindKey = bindKey;
    lodash.chain = chain;
    lodash.compact = compact;
    lodash.compose = compose;
    lodash.constant = constant;
    lodash.countBy = countBy;
    lodash.create = create;
    lodash.createCallback = createCallback;
    lodash.curry = curry;
    lodash.debounce = debounce;
    lodash.defaults = defaults;
    lodash.defer = defer;
    lodash.delay = delay;
    lodash.difference = difference;
    lodash.filter = filter;
    lodash.flatten = flatten;
    lodash.forEach = forEach;
    lodash.forEachRight = forEachRight;
    lodash.forIn = forIn;
    lodash.forInRight = forInRight;
    lodash.forOwn = forOwn;
    lodash.forOwnRight = forOwnRight;
    lodash.functions = functions;
    lodash.groupBy = groupBy;
    lodash.indexBy = indexBy;
    lodash.initial = initial;
    lodash.intersection = intersection;
    lodash.invert = invert;
    lodash.invoke = invoke;
    lodash.keys = keys;
    lodash.map = map;
    lodash.mapValues = mapValues;
    lodash.max = max;
    lodash.memoize = memoize;
    lodash.merge = merge;
    lodash.min = min;
    lodash.omit = omit;
    lodash.once = once;
    lodash.pairs = pairs;
    lodash.partial = partial;
    lodash.partialRight = partialRight;
    lodash.pick = pick;
    lodash.pluck = pluck;
    lodash.property = property;
    lodash.pull = pull;
    lodash.range = range;
    lodash.reject = reject;
    lodash.remove = remove;
    lodash.rest = rest;
    lodash.shuffle = shuffle;
    lodash.sortBy = sortBy;
    lodash.tap = tap;
    lodash.throttle = throttle;
    lodash.times = times;
    lodash.toArray = toArray;
    lodash.transform = transform;
    lodash.union = union;
    lodash.uniq = uniq;
    lodash.values = values;
    lodash.where = where;
    lodash.without = without;
    lodash.wrap = wrap;
    lodash.xor = xor;
    lodash.zip = zip;
    lodash.zipObject = zipObject;

    // add aliases
    lodash.collect = map;
    lodash.drop = rest;
    lodash.each = forEach;
    lodash.eachRight = forEachRight;
    lodash.extend = assign;
    lodash.methods = functions;
    lodash.object = zipObject;
    lodash.select = filter;
    lodash.tail = rest;
    lodash.unique = uniq;
    lodash.unzip = zip;

    // add functions to `lodash.prototype`
    mixin(lodash);

    /*--------------------------------------------------------------------------*/

    // add functions that return unwrapped values when chaining
    lodash.clone = clone;
    lodash.cloneDeep = cloneDeep;
    lodash.contains = contains;
    lodash.escape = escape;
    lodash.every = every;
    lodash.find = find;
    lodash.findIndex = findIndex;
    lodash.findKey = findKey;
    lodash.findLast = findLast;
    lodash.findLastIndex = findLastIndex;
    lodash.findLastKey = findLastKey;
    lodash.has = has;
    lodash.identity = identity;
    lodash.indexOf = indexOf;
    lodash.isArguments = isArguments;
    lodash.isArray = isArray;
    lodash.isBoolean = isBoolean;
    lodash.isDate = isDate;
    lodash.isElement = isElement;
    lodash.isEmpty = isEmpty;
    lodash.isEqual = isEqual;
    lodash.isFinite = isFinite;
    lodash.isFunction = isFunction;
    lodash.isNaN = isNaN;
    lodash.isNull = isNull;
    lodash.isNumber = isNumber;
    lodash.isObject = isObject;
    lodash.isPlainObject = isPlainObject;
    lodash.isRegExp = isRegExp;
    lodash.isString = isString;
    lodash.isUndefined = isUndefined;
    lodash.lastIndexOf = lastIndexOf;
    lodash.mixin = mixin;
    lodash.noConflict = noConflict;
    lodash.noop = noop;
    lodash.now = now;
    lodash.parseInt = parseInt;
    lodash.random = random;
    lodash.reduce = reduce;
    lodash.reduceRight = reduceRight;
    lodash.result = result;
    lodash.runInContext = runInContext;
    lodash.size = size;
    lodash.some = some;
    lodash.sortedIndex = sortedIndex;
    lodash.template = template;
    lodash.unescape = unescape;
    lodash.uniqueId = uniqueId;

    // add aliases
    lodash.all = every;
    lodash.any = some;
    lodash.detect = find;
    lodash.findWhere = find;
    lodash.foldl = reduce;
    lodash.foldr = reduceRight;
    lodash.include = contains;
    lodash.inject = reduce;

    mixin(function() {
      var source = {}
      forOwn(lodash, function(func, methodName) {
        if (!lodash.prototype[methodName]) {
          source[methodName] = func;
        }
      });
      return source;
    }(), false);

    /*--------------------------------------------------------------------------*/

    // add functions capable of returning wrapped and unwrapped values when chaining
    lodash.first = first;
    lodash.last = last;
    lodash.sample = sample;

    // add aliases
    lodash.take = first;
    lodash.head = first;

    forOwn(lodash, function(func, methodName) {
      var callbackable = methodName !== 'sample';
      if (!lodash.prototype[methodName]) {
        lodash.prototype[methodName]= function(n, guard) {
          var chainAll = this.__chain__,
              result = func(this.__wrapped__, n, guard);

          return !chainAll && (n == null || (guard && !(callbackable && typeof n == 'function')))
            ? result
            : new lodashWrapper(result, chainAll);
        };
      }
    });

    /*--------------------------------------------------------------------------*/

    /**
     * The semantic version number.
     *
     * @static
     * @memberOf _
     * @type string
     */
    lodash.VERSION = '2.4.2';

    // add "Chaining" functions to the wrapper
    lodash.prototype.chain = wrapperChain;
    lodash.prototype.toString = wrapperToString;
    lodash.prototype.value = wrapperValueOf;
    lodash.prototype.valueOf = wrapperValueOf;

    // add `Array` functions that return unwrapped values
    forEach(['join', 'pop', 'shift'], function(methodName) {
      var func = arrayRef[methodName];
      lodash.prototype[methodName] = function() {
        var chainAll = this.__chain__,
            result = func.apply(this.__wrapped__, arguments);

        return chainAll
          ? new lodashWrapper(result, chainAll)
          : result;
      };
    });

    // add `Array` functions that return the existing wrapped value
    forEach(['push', 'reverse', 'sort', 'unshift'], function(methodName) {
      var func = arrayRef[methodName];
      lodash.prototype[methodName] = function() {
        func.apply(this.__wrapped__, arguments);
        return this;
      };
    });

    // add `Array` functions that return new wrapped values
    forEach(['concat', 'slice', 'splice'], function(methodName) {
      var func = arrayRef[methodName];
      lodash.prototype[methodName] = function() {
        return new lodashWrapper(func.apply(this.__wrapped__, arguments), this.__chain__);
      };
    });

    return lodash;
  }

  /*--------------------------------------------------------------------------*/

  // expose Lo-Dash
  var _ = runInContext();

  // some AMD build optimizers like r.js check for condition patterns like the following:
  if (typeof define == 'function' && typeof define.amd == 'object' && define.amd) {
    // Expose Lo-Dash to the global object even when an AMD loader is present in
    // case Lo-Dash is loaded with a RequireJS shim config.
    // See http://requirejs.org/docs/api.html#config-shim
    root._ = _;

    // define as an anonymous module so, through path mapping, it can be
    // referenced as the "underscore" module
    define(function() {
      return _;
    });
  }
  // check for `exports` after `define` in case a build optimizer adds an `exports` object
  else if (freeExports && freeModule) {
    // in Node.js or RingoJS
    if (moduleExports) {
      (freeModule.exports = _)._ = _;
    }
    // in Narwhal or Rhino -require
    else {
      freeExports._ = _;
    }
  }
  else {
    // in a browser or Rhino
    root._ = _;
  }
}.call(this));

}).call(this,typeof self !== "undefined" ? self : typeof window !== "undefined" ? window : {})
},{}]},{},[2])