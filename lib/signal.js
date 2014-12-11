/**
  * (Advanced) Signal processing
  */
  

var _ = require("lodash");

function defaultSampling(options){
  var defaultSpectrum = 22050;
  if(!options.disableWarnings || !_.contains(options.disableWarnings, "defaultSampling")) {
    console.warn("switching to default Sampling rate of "+defaultSpectrum+".");
  }
  return defaultSpectrum;
}

var signal = function(sig, options){
  if("values" in sig && "sampling" in sig) return sig;
  options = options || {};
  var sampling = options.sampling || sig.sampling || defaultSampling(options);
  var curT = 0;
  var data = _.map(sig, function(v){
    ret = {t: curT, value: v};
    curT = curT + 1/sampling;
    return ret;
  });
  data.sampling = sampling;
  data.dt = 1/sampling;
  data.values = function(){
    return _.map(data, function(p){ return p.value; });
  };
  return data;
}

module.exports = signal;
