module.exports = function(grunt) {
  "use strict";
  require( 'load-grunt-tasks' )( grunt );

  grunt.initConfig({

    clean: {
      api_docs: [ 'docs/api/**' ]
    },

    jsdoc: {
      src: ['lib/**/*.js', 'bin/*.js'],
      options: {
        destination: 'docs/api'
      }
    },

    jsonlint: {
      pkg: {
        src: [
          'package.json',
          'test/server/data/fixtures/*.json'
        ]
      }
    },

    jscs: {
      main: [
        'bin/*.js',
        'lib/**/*.js'
      ],
      options: {
        config: '.jscsrc'
      }
    },

    jshint: {
      all: [
        'Gruntfile.js',
        'bin/*.js',
        'lib/**/*.js',
        'test/**/*.js',
      ],
      options: {
        jshintrc: '.jshintrc'
      }
    },

    jasmine_nodejs: {
      // task specific (default) options
      options: {
        specNameSuffix: "spec.js", // also accepts an array
        helperNameSuffix: "helper.js",
        useHelpers: false,
        random: false,
        seed: null,
        defaultTimeout: 5000,
        stopOnFailure: false,
        traceFatal: 2,
        // configure one or more built-in reporters
        reporters: {
          console: {
            colors: true,        // (0|false)|(1|true)|2
            cleanStack: 0,       // (0|false)|(1|true)|2|3
            verbosity: 4,        // (0|false)|1|2|3|(4|true)
            listStyle: "indent", // "flat"|"indent"
            activity: false
          },
        },
        // add custom Jasmine reporter(s)
        customReporters: []
      },
      'server_tests': {
        // target specific options
        //options: {
        //    useHelpers: true
        //},
        // spec files
        specs: [
          "test/server/**"
        ]
        //helpers: [
        //    "test/helpers/**"
        //]
      },
      'client_tests': {
        specs: [ "test/client/*.js" ]
      }
    },

    watch: {
      js: {
        files: [
          'Gruntfile.js',
          '.jshintrc',
          '.jscsrc',
          'bin/**/*.js',
          'lib/**/*.js',
          'test/client/*.js'
        ],
        tasks: [
          'lint', 'test'
        ]
      }
    }
  });

  grunt.registerTask('lint',    ['jshint', 'jscs', 'jsonlint']);
  grunt.registerTask('jasmine', ['jasmine_nodejs']);
  grunt.registerTask('test',    ['lint', 'jasmine']);
  grunt.registerTask('doc',     ['clean:api_docs', 'jsdoc']);
  grunt.registerTask('default', ['test']);
};
