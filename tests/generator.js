
var chai = require("chai");
chai.should();

var generator = require("../lib/generator.js");
var signal = require("../lib/signal.js");
var _ = require("lodash");

describe("Signal generator", function(){
  describe("Sine generator", function(){
    it("should start at origin", function(){
      var signal = generator
        .sine()
        .create()
        .values();
      signal[0].should.equal(0);
    });
    it("should create a signal of specified length", function(){
      var signal = generator
        .sine()
        .create({length: 100})
        .data();
      signal.length.should.equal(100);
    });
    it("should use the given options", function(){
      var signal = generator
        .sine({frequency: 340, amplitude: 3})
        .create({sampling: 2200})
        .data();
      signal.frequency.should.equal(340);
      signal.amplitude.should.equal(3);
      signal.sampling.should.equal(2200);
    });
    it("should equal a normal sine", function(){
      var sig = generator
        .sine({frequency: 1, amplitude: 1})
        .create({sampling: 10,length:10})
        .values();
      var sig2 = _.map(_.range(10), function(idx){ return Math.sin(Math.PI * idx / 5); });
      var areEqual = signal.equal(sig,sig2,{epsilon:1e-5});
      areEqual.should.be.true;
    });
    it("should be able to generate multiple sines", function(){
      var sig = generator
        .sines([{frequency: 1, amplitude: 1},
          {frequency: 1, amplitude: 1},
          {frequency: 1, amplitude: 1}])
        .create({sampling: 10, length: 10})
        .values();
      var sig2 = _.map(_.range(10), function(idx){ return 3*Math.sin(Math.PI * idx / 5); });
      var areEqual = signal.equal(sig,sig2,{epsilon:1e-5});
      areEqual.should.be.true;
    });
  });
});
