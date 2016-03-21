module.exports = function(grunt) {
  "use strict";
  require( 'load-grunt-tasks' )( grunt );

  grunt.initConfig({
    jscs: {
      main: [ 'Gruntfile.js',
              'bin/server.js',
      ],
      options: {
        config: '.jscsrc'
      }
    },
    jshint: {
      all: [
        'Gruntfile.js',
        'bin/server.js'
      ],
      options: {
        jshintrc: '.jshintrc'
      }
    },
    watch: {
      js: {
        files:[
          'Gruntfile.js',
          '.jshintrc',
          '.jscsrc',
          'bin/**/*.js',
        ],
        tasks: [
          'test'
        ]
      }
    }
  });

  grunt.registerTask('lint', ['jshint', 'jscs']);
  grunt.registerTask('test', ['lint']);
  grunt.registerTask('default', ['test']);
};
