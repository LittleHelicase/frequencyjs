
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

var Spectrum = require("../spectrum.js");

module.exports = {
  name: "FFT DSP.JS",
  register: function(Transform){
    Transform.register({
      name: "fft",
      forward: function(signal, options){
        options.length = signal.length;
        var fft = prepare(options);
        fft.forward(signal.values());
        return Spectrum(fft.spectrum,options);
      },
      backward: function(spectrum, options){
        return generator.sines(spectrum);
      }
    });
  }
};
