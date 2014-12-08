/**
  * (Advanced) Signal processing
  */


var cauchy = require("./convolution/cauchy.js");
var specConvolution = require("./convolution/spectrum.js");
var _ = require("lodash");

var defaultOptions = {
  method: "cauchy",
  type: "non-circular",
  spectrum: {
    method: "dft"
    // sampling is not relevant as it gets transformed back again with the same
    // sampling rate
  }
}

module.exports = {
  /** the convolution function assumes two discrete signals
    * that have the same spacing
    */
  convolve: function(signal1, signal2, options){
    options = options || {};
    options = _.defaults(options, defaultOptions);
    
    // ensure signal2 is not longer than signal2
    if(signal1.length < signal2.length){
      var tmp = signal1;
      signal1 = signal2;
      signal2 = tmp;
    }
    
    if(options.method == "cauchy"){
      return cauchy(signal1,signal2,options);
    } else {
      return specConvolution(signal1, signal2,options);
    }
  },
  /** determines if two signals can be considered equal.
    */
  equal: function(signal1, signal2, options){
    options = options || {};
    options.epsilon = options.epsilon || 1E-08;
    if(signal1.length != signal2.length) return false;
    var diffSQ = _.reduce(_.range(signal1.length), function(d,idx){
      diff = (signal1[idx]-signal2[idx]);
      return d + diff * diff;
    },0);
    return diffSQ/signal1.length < options.epsilon;
  }
}
