
var haar = function(signal, options){
  var input = signal.values();
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
