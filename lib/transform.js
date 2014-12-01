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

function getMethod(methodID, sampling, length) {
  if (methodID === "dft") {
    return new dsp.DFT(length, sampling);
  }
  if (methodID === "fft") {
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
    throw new Error("Not implemented yet..");
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
    return new Spectrum(method.spectrum, sampling, signal.length);
  }
}
