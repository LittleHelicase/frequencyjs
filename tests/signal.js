
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
  it("alternating sequences", function(){
    // sig = [1,0,-1,0,...]
    var sig = _.map(_.range(32), function(idx){
      return Math.round(Math.cos(Math.PI / 2 * idx));
    });
    var conv = Signal.convolve(sig,sig);
    conv.should.deep.equal(sig);
  });
});
