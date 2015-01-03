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
