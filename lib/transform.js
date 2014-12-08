/**
  * Transforms a signal into its frequencies or vice versa
  */

var dsp = require("../dsp.js");
var _ = require("lodash");

function defaultSampling(options){
  var defaultSpectrum = 22050;
  if(!options.disableWarnings || !_.contains(options.disableWarnings, "defaultSampling")) {
    console.warn("switching to default Sampling rate of "+defaultSpectrum+".");
  }
  return defaultSpectrum;
}

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
  return {
    data: function(){
      return spectrum;
    },
    dominantFrequency: function() {
      return _.max(spectrum, function (s) {
        return s.amplitude;
      });
    }
  }
}

module.exports = {
  toSignal: function(spectrum, options){
    var sampling = options.sampling || signal.sampling || defaultSampling(options);
    var methodID = options.method || "dft";
    var method = getMethod(methodID, sampling, signal.length);
    var spectrumData = spectrum;
    if(spectrum.length == 0) return null;
    if(_.isObject(signal[0])){
      signalData = _.maps(spectrum, function(v){ v.value; });
    }
    method.inverse(spectrumData);
    return Signal(method.signal);
  },
  toSpectrum: function(signal, options){
    var sampling = options.sampling || signal.sampling || defaultSampling(options);
    var methodID = options.method || "dft";
    var method = getMethod(methodID, sampling, signal.length);
    var signalData = signal;
    if(signal.length == 0) return null;
    if(_.isObject(signal[0])){
      signalData = _.map(signal, function(v){ return v.value; });
    }
    method.forward(signalData);
    return Spectrum(method.spectrum, sampling, signal.length);
  }
}
