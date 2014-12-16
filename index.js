
var Generators = require("./lib/generator.js");
var Transform = require("./lib/transform.js");
var Signal = require("./lib/signal.js");
var Processing = require("./lib/processing.js");

module.exports = {
  Processing: Processing,
  Generators: Generators,
  Transform: Transform,
  toSpectrum: Transform.toSpectrum,
  toSignal: Transform.toSignal,
  signal: Signal,
  sine: Generators.sine,
  sines: Generators.sines,
  convolve: Processing.convolve,
  version: "0.0.3"
}
