
var _ = require("lodash");

var spectrumCreate = function(spectrumArray, sampling, signalLength){
  var spectrum = _.map(_.range(spectrumArray.length), function(idx){
    // frequency in hertz
    var freq = idx/signalLength * sampling;
    var freqFac = (idx == 0) ? 0.5 : 1;
    return { frequency: freq, amplitude: spectrumArray[idx]*freqFac };
  });
  spectrum.dominantFrequency = function() {
    return _.max(spectrum, function (s) {
      return s.amplitude;
    });
  };
  spectrum.amplitudes = function(){
    return _.map(spectrum, function(v){
      return v.amplitude;
    });
  }
  spectrum.sampling = sampling;
  return spectrum;
};

var timeDependentSpectrumCreate = function(spectrumArray, sampling){
  var levels = Math.floor(Math.log(spectrumArray.length)/Math.log(2));
  var curLevel = levels-1;
  var idx = Math.pow(2,curLevel);
  var spectrum = [];
  var cnt = 0;
  while(curLevel>=0){
    var freq = sampling * idx /* = 2^curLevel */ / spectrumArray.length;
    for(var i=idx; i<idx*2; i++){
      var diffIdx = i - idx;
      var stepSize = spectrumArray.length/(idx * sampling);
      spectrum[i] = {
        frequency: freq,
        amplitude: spectrumArray[i],
        timeStart: diffIdx * stepSize,
        timeEnd:(diffIdx+1)*stepSize,
        time: (diffIdx + 0.5)*stepSize
      };
    }
    curLevel = curLevel - 1;
    idx = idx / 2;
    cnt++;
  }
  spectrum[0] = { frequency: 0, amplitude: spectrumArray[0] };
  spectrum.dominantFrequency = function() {
    return _.max(spectrum, function (s) {
      return s.amplitude;
    });
  };
  spectrum.amplitudes = function(){
    return _.map(spectrum, function(v){
      return v.amplitude;
    });
  }
  spectrum.sampling = sampling;
  return spectrum;
};

var Spectrum = function(spectrumArray, options){
  if("amplitudes" in spectrumArray)
    return spectrumArray;
  if("timeDependent" in options && options.timeDependent)
    return timeDependentSpectrumCreate(spectrumArray, options.sampling);
  return spectrumCreate(spectrumArray, options.sampling, options.signalLength);
}

module.exports = Spectrum;
