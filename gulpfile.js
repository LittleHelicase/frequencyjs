var gulp = require("gulp");
var mocha = require('gulp-mocha');
var browserify = require('gulp-browserify');
var rename = require("gulp-rename");

var testFiles = [
  'tests/*.js'
];

gulp.task('test', function() {
  // Be sure to return the stream
  return gulp.src(testFiles)
    .pipe(mocha())
    .on('error', function(err) {
      // Make sure failed tests cause gulp to exit non-zero
      throw err;
    });
});

gulp.task('default', function() {
  gulp.src("index.js")
    .pipe(browserify({}))
    .pipe(rename("fjs.js"))
    .pipe(gulp.dest("./dist/"));
});
