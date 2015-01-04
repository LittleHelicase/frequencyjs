
var haar = function(signal, options){
  var res = signal.values();
  var copy = [];
  var len = Math.floor(signal.length / 2);
  while(len > 0){
    for(var i=0; i<len; i++){
      var scaling = (res[2*i] + res[2*i+1]);
      var wavelet = (res[2*i] - res[2*i+1]);
      copy[i] = scaling;
      copy[len + i] = wavelet;
    }
    var tmp = copy;
    copy = res;
    res = tmp;
    len = Math.floor(len / 2);
  }
  return res;
}

module.exports = {
  name: "DWT Haar",
  register: function(Transform){
    Transform.register({
      name: "haar",
      forward: function(signal, options){
        var trans = haar(signal);
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
