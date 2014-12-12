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
  return transform.toSignal(convSpec,options.spectrum);
}
