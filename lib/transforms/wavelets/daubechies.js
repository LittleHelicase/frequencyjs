
var daubCoeffs = require("./daubechiesCoefficients.js");

var daubechiesPeriodic = function(signal, options){
  var coeffs = daubCoeffs["D"+options.taps];
  var input = signal.values();
  var copy = [];
  var res = [];
  var len = Math.floor(signal.length / 2);
  while(len > 0){
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

module.exports = {
  name: "DWT Daubechies",
  register: function(Transform){
    Transform.register({
      name: "daubechies",
      forward: function(signal, options){
        var trans = daubechiesPeriodic(signal, options);
        trans.timeDependent = true;
        return trans;
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
