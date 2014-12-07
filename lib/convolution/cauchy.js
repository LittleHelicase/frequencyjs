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

var circularIndex = function(i1,i2,n){
  return mod(i1 - i2,n);
}

var index = function(i1,i2,n){
  var i = i1-i2;
  return (0<=i&&i<n) ? i : null;
}

var indexerForType = {"non-circular" : index,"circular": circularIndex};

module.exports = function(signal1, signal2, options){
  var len1 = signal1.length;
  var len2 = signal2.length;
  var half = Math.floor(len2 / 2);
  var idx = indexerForType[options.type];
  
  var convolved = _.map(_.range(len1), function(idx1){
    var cVal = _.reduce(_.range(len2), function(acc, idx2){
      var curIdx = idx(idx1,idx2,len1);
      // range test. If our index is not within the signal1 
      // we simply return the current value of the convolution
      if(curIdx == null){
        return acc;
      }
      // actual convolution
      return acc + signal1[curIdx] * signal2[idx2];
    }, 0);
    return cVal;
  });
  
  return convolved;
}
