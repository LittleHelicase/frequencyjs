/**
  * convolution using a spectrum (usually fft) method to calculate the spectrum
  * multiply them and backtransform back into a signal. Only for periodic
  * signals
  */

var transform = require("../transform.js");

module.exports = function(signal1, signals2, options){
  var spec1 = transform.toSpectrum(signal1, options.spectrum);
  var spec2 = transform.toSpectrum(signal2, options.spectrum);
  // spec1.len > spec2.len , but the shorter spectrum has only 
  // the high frequencies => multiply from end to start
//  var convSpec = _.map(_.range());
}
