
var Spectrum = require("../spectrum.js");
var Signal = require("../signal.js");
var generator = require("../generator.js");

Transform = {
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

module.exports = Transform;
