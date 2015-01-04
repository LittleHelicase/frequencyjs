
var chai = require("chai");
var _ = require("lodash");
chai.should();

var Generator = require("../lib/generator.js");
var Transform = require("../lib/transform.js");
var Processing = require("../lib/processing.js");
var Signals = require("../lib/utils/signals.js");

describe("Haar Wavelet Transform", function(){
  it("should detect different (rectangle) frequencies", function(){
    //for(var i=1; i<4; i++) {
    for(var i=1; i<4; i++) {
      var freq = Math.pow(2,i-1);
      var rectFunc = Generator.rectangle({frequency:freq})
        .create({sampling:8,length:8})
        .values();
      // calculate the dominant frequency of the signal using dft
      var domFreq = Transform
      .toSpectrum(rectFunc,{method: "haar", sampling:8})
      .dominantFrequency();
      domFreq.frequency.should.equal(Math.pow(2,i - 1))
    }
  });
  it("should transform an impuls correctly", function(){
    // example taken from http://en.wikipedia.org/wiki/Discrete_wavelet_transform#Comparison_with_Fourier_transform
    var impuls = [1,0,0,0];
    var spec = Transform
      .toSpectrum(impuls,{method:"haar", sampling: 4})
      .amplitudes();
    var areEqual = Processing.equal(spec, [0.25,0.25,0.5,0]);
    areEqual.should.be.true;
  })
});
