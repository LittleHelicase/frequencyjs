
var dsp = require("../../dsp.js");
var generator = require("../generator.js");
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

var Spectrum = require("../spectrum.js");

module.exports = {
  name: "DFT DSP.JS",
  register: function(Transform){
    Transform.register({
      name: "dft",
      forward: function(signal, options){
        options.length = signal.length;
        var dft = prepare(options);
        dft.forward(signal.values());
        return Spectrum(dft.spectrum, options);
      },
      backward: function(spectrum, options){
        return generator.sines(spectrum);
      }
    });
  }
};
