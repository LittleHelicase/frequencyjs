
var chai = require("chai");
chai.should();

var generator = require("../lib/generator.js");

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
    })
  });
});
