
var chai = require("chai");
var _ = require("lodash");
chai.should();

var Generator = require("../lib/generator.js");
var Transform = require("../lib/transform.js");
var Processing = require("../lib/processing.js");

describe("Discrete Fourier Transform", function(){
  it("should detect the frequency of the sine generator", function(){
    // generate signal of 440 hertz
    var signal = Generator
        .sine({frequency: 440})
        .create({length: 100, sampling: 4400});

    // calculate the dominant frequency of the signal using dft
    var domFreq = Transform
      .toSpectrum(signal,{method: "dft"})
      .dominantFrequency();
    domFreq.frequency.should.equal(440);
  });
  it("should detect different frequencies", function(){
    for(var i=1; i<4; i++) {
      cosFunc = _.map(_.range(8), function (idx) {
        return Math.cos(i * 2 * Math.PI / 8 * idx);
      });

      // calculate the dominant frequency of the signal using dft
      var domFreq = Transform
          .toSpectrum(cosFunc,{method: "dft", sampling:8})
          .dominantFrequency();
      domFreq.frequency.should.equal(i)
    }
  });
});
describe("Backtransfrom via Generators", function(){
  it("transform and backtransform should cancle out", function(){
    var signal = Generator.sine({frequency: 32}).create({length:128,sampling:128});
    var spec = Transform.toSpectrum(signal);
    var baseSignal = Transform.toSignal(spec).create({length:128,sampling:128});
    var areEqual = Processing.equal(signal, baseSignal);
    areEqual.should.be.true;
  });
})
