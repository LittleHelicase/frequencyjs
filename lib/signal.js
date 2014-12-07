/**
  * (Advanced) Signal processing
  */


var cauchy = require("./convolution/cauchy.js");
var _ = require("lodash");

var defaultOptions = {
  method: "cauchy",
  type: "non-circular"
}

module.exports = {
  /** the convolution function assumes two discrete signals
    * that have the same spacing
    */
  convolve: function(signal1, signal2, options){
    options = options || {};
    options = _.defaults(options, defaultOptions);
    
    if(options.method == "cauchy"){
      return cauchy(signal1,signal2,options);
    } else {
      throw new Error("fast convolution not implemented yet");
      return [];
    }
  }
}
