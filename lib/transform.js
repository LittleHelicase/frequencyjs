/**
  * Transforms a signal into its frequencies or vice versa
  */

var dsp = require("../dsp.js");
var _ = require("lodash");
var generator = require("./generator.js");
var Signal = require("./signal.js");

function isPowerOfTwo(number){
  // bit check. every power of two has only one 1 in its binary representation
  // e.g. 32 = 100000
  // and a power of 2 minus 1 has all the lower bits set
  // e.g. 31 = 011111
  // the following statement uses a binary AND and it is always false for a 
  // non power of two, and true for a power of two
  // (it's probably the fastes way to do it)
  return ((number!=0) && !(number & (number-1)));
}

function getMethod(methodID, sampling, length) {
  if (methodID === "dft") {
    return new dsp.DFT(length, sampling);
  }
  if (methodID === "fft") {
    if(!isPowerOfTwo(length)){
      throw new Error("FFT only works for a power of two signal length");
    }
    return new dsp.FFT(length, sampling);
  }
  throw new Error("Method " + methodID + " not supported");
}

var Spectrum = function(spectrumArray, sampling, signalLength){
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
}

module.exports = {
  toSignal: function(spectrum, options){
    options = options || {};
    var sampling = options.sampling || spectrum.sampling || defaultSampling(options);
    var methodID = options.method || "dft";
    var method = getMethod(methodID, sampling, spectrum.length);
    if(spectrum.length == 0) return null;
    return generator.sines(spectrum);
  },
  toSpectrum: function(signal, options){
    options = options || {};
    var methodID = options.method || "dft";
    var sig = Signal(signal, options);
    var method = getMethod(methodID, sig.sampling, sig.length);
    if(sig.length == 0) return null;
    method.forward(sig.values());
    return Spectrum(method.spectrum, sig.sampling, sig.length);
  }
}
