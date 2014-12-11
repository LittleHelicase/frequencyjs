
var chai = require("chai");
chai.should();

var Processing = require("../lib/processing.js");
var _ = require("lodash");

describe("Signal Convolution", function(){
  it("dirac impulse is the identity", function(){
    var sig1 = _.range(100);
    var sig2 = [1];
    var conv = Processing.convolve(sig1,sig2);
    conv.should.deep.equal(sig1);
  });
  it("should convolve constant sequences correctly", function(){
    // sig = [1,1,1,1...]
    var sig = _.map(_.range(32), function(idx){
      return 1
    });
    var conv = Processing.convolve(sig,sig,{type:"circular"});
    _.each(conv,function(v){
      v.should.equal(sig.length);
    })
  });
  it("should be able to use periodic signals", function(){
    var sig = _.map(_.range(32), function(){ return 1; });
    var pSig = [1];
    var sig_sig = Processing.convolve(sig,sig,{type:"circular"});
    var sig_pSig = Processing.convolve(sig,pSig,{type:"circular"});
    sig_sig.should.deep.equal(sig_pSig);
  })
});
