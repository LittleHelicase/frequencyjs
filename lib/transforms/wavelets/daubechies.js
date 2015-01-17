
var daubCoeffs = require("./daubechiesCoefficients.js");

var daubechiesPeriodic = function(signal, options){
  var coeffs = daubCoeffs["D"+options.taps];
  var input = signal.values();
  var copy = [];
  var res = [];
  var len = Math.floor(signal.length / 2);
  while(len >= options.taps){
    for(var i=0; i<len; i++){
      var scaling = 0; var wavelet = 0;
      var fac = 1;
      for(var j=0; j<coeffs.length; j++){
        var idx = (2*i + j) % len;
        scaling += coeffs[j] * input[2*i+j]*0.5;
        wavelet += fac * coeffs[j] * input[2*i+j]* 0.5;
        fac = fac * -1;
      }
      copy[i] = scaling;
      res[len + i] = wavelet;
    }
    var tmp = copy;
    copy = input;
    input = tmp;
    len = Math.floor(len / 2);
  }
  res[0] = input[0];
  return res;
}

var daubechiesReversePeriodic = function(spectrum, options){
  var taps = options.taps;
  var coeffs = daubCoeffs["D" + taps];
  var ampls = spectrum.slice();
  var copy = ampls.slice();
  var len = taps;
  while(len < spectrum.length){
    for(var i=0; i<len; i++){
      var f1 = 0; var f2 = 0;
      for(var j=0; j<taps; j++){
        f1 += ampls[i] + spectrum[i + len];
        f2 += ampls[i] - spectrum[i + len];
      }
      copy[2*i] = f1;
      copy[2*i + 1] = f2;
    }
    var tmp = ampls;
    ampls = copy;
    copy = tmp;
    len = len * 2;
  }
  return ampls;
}

var Spectrum = require("../../spectrum.js");

module.exports = {
  name: "DWT Daubechies",
  register: function(Transform){
    Transform.register({
      name: "daubechies",
      forward: function(signal, options){
        options = options || {};
        options.taps = options.taps || 2;
        var trans = daubechiesPeriodic(signal, options);
        trans.timeDependent = true;
        return Spectrum(trans,options);
      },
      backward: function(spectrum, options){
        options = options || {};
        options.taps = options.taps || 2;
        var rev = daubechiesReversePeriodic(spectrum, options);
        return Signal(rev, 
          {
            sampling: 2*spectrum[spectrum.length-1].frequency
          });
      }
    });
  }
};
