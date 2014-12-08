/**
  * implementation of the discrete convolution for arbitrary signals
  * doesn't require any periodicity like the fast convolution
  * 
  * preferable when the second signal (filter) is short.
  */

var _ = require("lodash");

var mod = function(m,n){
  return ((m % n) + n) % n;
}

var nonCircular = function(signal1, signal2){
  var len1 = signal1.length;
  var len2 = signal2.length;
  var half = Math.floor(len2 / 2);
  
  var convolved = _.map(_.range(len1), function(idx1){
    var cVal = _.reduce(_.range(len2), function(acc, idx2){
      var curIdx = idx1 - idx2;
      if(0 <= curIdx && curIdx < len1){
        // actual convolution
        return acc + signal1[curIdx] * signal2[idx2];
      }
      // range test. If our index is not within the signal1 
      // we simply return the current value of the convolution
      return acc;
    }, 0);
    return cVal;
  });
  
  return convolved;
}

var circular = function(signal1, signal2){
  var len1 = signal1.length;
  var len2 = signal2.length;
  var half = Math.floor(len2 / 2);
  
  var convolved = _.map(_.range(len1), function(idx1){
    var cVal = _.reduce(_.range(len1), function(acc, idx2){
      var curIdx = idx1 - idx2;
      return acc + signal1[idx2] * signal2[mod(curIdx,len2)];
    }, 0);
    return cVal;
  });
  
  return convolved;
}

var methodForType = {"non-circular" : nonCircular,"circular": circular};

module.exports = function(signal1, signal2, options){
  var convMethod = methodForType[options.type];
  return convMethod(signal1,signal2);
}
