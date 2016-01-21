var gulp = require('gulp');
var source = require('vinyl-source-stream');
var browserify = require('browserify');
var watchify = require('watchify');
var browserifyShader = require("browserify-shader")
var tsify = require('tsify');
var uglify = require('uglifyify');

gulp.task('.bower.install', function () {
    var bower = require('gulp-bower');
    return bower();
});

gulp.task('.clean.bower.lib', function (cb) {
    var del = require('del');
    del(['lib/'], cb);
});

gulp.task('clean.npm', function (cb) {
    var del = require('del');
    del(['node_modules/'], cb);
});

gulp.task('watch', function() {
    var bundler = watchify(browserify({debug: true})
        .add('alternator.ts')
        .plugin(tsify)
        .transform(browserifyShader));

    bundler.on('update', rebundle)
 
    function rebundle () {
        return bundler.bundle()
          .pipe(source('bundle.js'))
          .pipe(gulp.dest('.'))
    }
     
    return rebundle();
});

gulp.task('build.debug', function() {
    var bundler = browserify({debug: true})
        .add('alternator.ts')
        .plugin(tsify)
        .transform(browserifyShader);

    return bundler.bundle()
        .pipe(source('bundle.js'))
        .pipe(gulp.dest('.'));
});

gulp.task('build.release', function() {
    var bundler = browserify()
        .add('alternator.ts')
        .plugin(tsify)
        .transform(browserifyShader)
        .transform(uglify);

    return bundler.bundle()
        .pipe(source('bundle.js'))
        .pipe(gulp.dest('.'));
});