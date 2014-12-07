
var chai = require("chai");
chai.should();

var Signal = require("../lib/signal.js");
var _ = require("lodash");

describe("Signal Convolution", function(){
  it("dirac impulse is the identity", function(){
    var sig1 = _.range(100);
    var sig2 = [1];
    var conv = Signal.convolve(sig1,sig2);
    conv.should.deep.equal(sig1);
  });
  it("constant sequence", function(){
    // sig = [1,1,1,1...]
    var sig = _.map(_.range(32), function(idx){
      return 1
    });
    var conv = Signal.convolve(sig,sig,{type:"circular"});
    _.each(conv,function(v){
      v.should.equal(sig.length);
    })
  });
});
