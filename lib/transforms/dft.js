
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
