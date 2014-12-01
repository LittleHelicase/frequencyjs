
var chai = require("chai");
var _ = require("lodash");
chai.should();

var Generator = require("../lib/generator.js");
var Transform = require("../lib/transform.js");

describe("Fast Fourier Transform", function(){
    it("should detect the frequency of the sine generator", function(){
        // generate signal of 32 hertz (only powers of two are used in fft)
        var signal = Generator
            .sine({frequency: 32})
            .create({length: 128, sampling: 128})
            .data();
        var domFreq = Transform
            .toSpectrum(signal,{method: "fft"})
            .dominantFrequency();
        domFreq.frequency.should.equal(32);
    });
    it("should detect different frequencies", function(){
        for(var i=1; i<4; i++) {
            cosFunc = _.map(_.range(8), function (idx) {
                return Math.cos(i * 2 * Math.PI / 8 * idx);
            });
            var domFreq = Transform
                .toSpectrum(cosFunc,{method: "fft", sampling:8})
                .dominantFrequency();
            domFreq.frequency.should.equal(i)
        }
    });
});
