/**
  * (Advanced) Signal processing
  */

var assert = require("assert");
var _ = require("lodash");

module.exports = {
  /** the convolution function assumes two discrete signals
    * that have the same spacing
    */
  convolve: function(signal1, signal2){
    var len1 = signal1.length;
    var len2 = signal2.length;
    var half = Math.floor(len2 / 2);
    
    var convolved = _.map(_.range(len1), function(idx1){
      var cVal = _.reduce(_.range(len2), function(acc, idx2){
        var curIdx = idx1 + idx2 - half;
        // range test. If our index is not within the signal1 
        // we simply return the current value of the convolution
        if(curIdx < 0 || curIdx >= len1){
          return acc;
        }
        // actual convolution
        return acc + signal1[curIdx] * signal2[idx2];
      }, 0);
      return cVal;
    });
    
    return convolved;
  }
}
