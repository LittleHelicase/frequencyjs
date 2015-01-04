
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
  var idx = 0;
  var spectrum = [];
  var cnt = 0;
  while(2*idx < spectrumArray.length){
    var freq = (idx == 0) ? 0 : sampling / Math.pow(2,cnt);
    for(var i=idx; i<idx*2+1; i++){
      var diffIdx = i - idx;
      var stepSize = (spectrumArray.length/(idx+1))/sampling;
      spectrum.push({
        frequency: freq,
        amplitude: spectrumArray[i],
        timeStart: diffIdx * stepSize,
        timeEnd:(diffIdx+1)*stepSize,
        time: (diffIdx + 0.5)*stepSize
      });
    }
    idx = idx*2+1;
    cnt++;
  }
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
  if("timeDependent" in spectrumArray && spectrumArray.timeDependent)
    return timeDependentSpectrumCreate(spectrumArray, options.sampling);
  return spectrumCreate(spectrumArray, options.sampling, options.signalLength);
}

module.exports = Spectrum;
