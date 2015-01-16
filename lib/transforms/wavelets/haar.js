
var haar = function(signal){
  var input = signal;
  var copy = [];
  var res = [];
  var len = Math.floor(signal.length / 2);
  while(len > 0){
    for(var i=0; i<len; i++){
      var scaling = (input[2*i] + input[2*i+1])*0.5;
      var wavelet = (input[2*i] - input[2*i+1])*0.5;
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

var reverse_haar = function(spectrum){
  var ampls = spectrum.slice();
  var copy = ampls.slice();
  var len = 1;
  while(len < spectrum.length){
    for(var i=0; i<len; i++){
      var f1 = ampls[i] + spectrum[i + len];
      var f2 = ampls[i] - spectrum[i + len];
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

var Signal = require("../../signal.js");
var Spectrum = require("../../spectrum.js");

module.exports = {
  name: "DWT Haar",
  register: function(Transform){
    Transform.register({
      name: "haar",
      forward: function(signal){
        var trans = haar(signal.values());
        return Spectrum(trans,{sampling:signal.sampling, timeDependent: true});
      },
      backward: function(spectrum){
        return Signal(reverse_haar(spectrum.amplitudes()),
          {
            sampling: 2*spectrum[spectrum.length-1].frequency
          });
      }
    });
  }
};
