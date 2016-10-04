/* globals describe, it, expect */

"use strict";
const assert  = require('assert');
const os      = require('os');
const path    = require('path');
const config  = require('../lib/config.js');

describe('Building options', function() {
  it('Must define options before retrieval', function() {
    expect( () => {config.provide();} ).toThrowError(assert.AssertionError,
      'Options is undefined');
  });
  it('If there is a parameter, it must be a function', function() {
    expect( () => {config.provide( 1 );} ).toThrowError(assert.AssertionError,
      'parameter must be a function');
  });
  it('Function must return an object', function() {
    expect( () => {config.provide( () => {return 1;} );} ).toThrowError(assert.AssertionError,
      'parameter must return an object');
    expect( () => {config.provide( () => {return {};} );} ).not.toThrow();
  });
  it('Can get options passed via function', function() {
    let options;
    expect( () => {options = config.provide( () => {return {test: 'pass'};} );} ).not.toThrow();
    expect( options.get('test') === 'pass' ).toBe(true);
  });
  it('Options passed via function will overwrite defaults', function() {
    let options;
    expect( () => {options = config.provide( () => {return {};} );} ).not.toThrow();
    expect( options.get('mongourl') === 'mongodb://localhost:27017/imetacache' ).toBe(true); // Provided in config.js
    expect( () => {options = config.provide( () => {return {mongourl: 'newmongourl'};} );} ).not.toThrow();
    expect( options.get('mongourl') === 'newmongourl' ).toBe(true);
  });
  it('Configs can be passed from a json file', function() {
    let options;
    expect( () => {options = config.provide( () => {
      return { configfile: path.resolve(__dirname, 'server', 'data', 'testConfig.json') };
    });} ).not.toThrow();
    expect( options.get('testConfig') ).toBe(true);
  });
  it('Configs from function will overwrite those from json file', function() {
    let options;
    expect( () => {options = config.provide( () => {
      return {
        configfile: path.resolve(__dirname, 'server', 'data', 'testConfig.json'),
        testConfig: false
      };
    });}).not.toThrow();
    expect( options.get('testConfig') ).toBe(false);
  });
  it('tempdir, port unspecified, defaults created', function() {
    let options;
    expect( () => {options = config.provide( () => {
      return {};
    });}).not.toThrow();
    expect( options.get('tempdir').startsWith(path.join(os.tmpdir(), 'npg_ranger_')) ).toBe(true);
    expect( options.get('port').startsWith(path.join(os.tmpdir(), 'npg_ranger_')) ).toBe(true);
    expect( options.get('port').endsWith('npg_ranger.sock') ).toBe(true);
  });
  it('tempdir specified, port default', function() {
    let options;
    expect( () => {options = config.provide( () => {
      return {tempdir: config.tempFilePath('npg_ranger_config_test_')};
    });}).not.toThrow();
    expect( options.get('tempdir').startsWith(path.join(os.tmpdir(), 'npg_ranger_config_test')) ).toBe(true);
    expect( options.get('port').startsWith(path.join(os.tmpdir(), 'npg_ranger_config_test'))).toBe(true);
    expect( options.get('port').endsWith('npg_ranger.sock') ).toBe(true);
  });
  it('tempdir default, port specified', function() {
    let options;
    expect( () => {options = config.provide( () => {
      return {port: '45678'};
    });}).not.toThrow();
    expect( options.get('tempdir').startsWith(path.join(os.tmpdir(), 'npg_ranger_')) ).toBe(true);
    expect( options.get('port') === 45678 ).toBe(true);
  });
  it('tempdir specified, port specified', function() {
    let options;
    expect( () => {options = config.provide( () => {
      return {
        tempdir: config.tempFilePath('npg_ranger_config_test_'),
        port:    '45678'
      };
    });}).not.toThrow();
    expect( options.get('tempdir').startsWith(path.join(os.tmpdir(), 'npg_ranger_config_test_')) ).toBe(true);
    expect( options.get('port') === 45678 ).toBe(true);
  });
});

describe('Creating temp file path', function() {
  it('Without prefix', function() {
    let temppath = config.tempFilePath();
    expect( temppath.startsWith(os.tmpdir()) ).toBe(true);
    expect( temppath ).toMatch(/\/\d{8}_\d{8}$/);
  });
  it('With prefix', function() {
    let temppath = config.tempFilePath('npg_ranger_config_test_');
    expect( temppath.startsWith(path.join(os.tmpdir(), 'npg_ranger_config_test_')) ).toBe(true);
    expect( temppath ).toMatch(/\/npg_ranger_config_test_\d{8}_\d{8}$/);
  });
});

describe('Listing config options', function() {
  it("Options listing", function() {
    config.provide( () => {return {mongourl: 'mymongourl',
                                   hostname: 'myhost',
                                   tempdir:  '/tmp/mydir',
                                   port:     9999,
                                   debug:    true,
                                   help:     true
                                  };} );
    let expected = ['anyorigin=undefined',
                    'configfile=undefined',
                    'debug=true',
                    'hostname=myhost',
                    'mongourl=mymongourl',
                    'multiref=undefined',
                    'port=9999',
                    'references=undefined',
                    'skipauth=undefined',
                    'tempdir=/tmp/mydir',
                    'timeout=3'];
    let o = config.logOpts();
    expect(o).not.toMatch(/help=/);
    expect(o).toMatch(expected.join(', '));
  });
});

describe('Validating CORS options', function() {

  it('anyorigin and originlist options cannot be set together', function() {
    let dir = config.provide().get('tempdir');
    expect( () => { config.provide( () => {return {mongourl: 'mymongourl',
                                   hostname: 'myhost',
                                   tempdir:  '/tmp/yourdir',
                                   port:     9999,
                                   debug:    true,
                                   help:     true,
                                   anyorigin: 'true'
                    };});
		  }).toThrowError("'anyorigin' should be a boolean type value");
    expect(config.provide().get('tempdir')).toBe(dir, 'old value retained');
  });

  it('anyorigin and originlist options cannot be set together', function() {
    let port = config.provide().get('port');
    expect( () => { config.provide( () => {return {mongourl: 'mymongourl',
                                   hostname: 'myhost',
                                   tempdir:  '/tmp/mydir',
                                   port:     port - 1,
                                   debug:    true,
                                   help:     true,
                                   anyorigin: true,
                                   originlist: 'url1,url2'
                    };});
		  }).toThrowError(
                    "'anyorigin' and 'originlist' options cannot both be set");
    expect(config.provide().get('port')).toBe(port, 'old value retained');
  });

  it('anyorigin is incompatible with authorization', function() {
    expect( () => { config.provide( () => {return {mongourl: 'mymongourl',
                                   hostname: 'myhost',
                                   tempdir:  '/tmp/mydir',
                                   port:     9999,
                                   debug:    true,
                                   help:     true,
                                   anyorigin: true,
                                   skipauth: false
                    };});
		  }).toThrowError(
                    "'anyorigin' option cannot be set if authorization is performed");
  });

  it('validation for originlist', function() {

     expect( () => { config.provide( () => {return {mongourl: 'mymongourl',
                                    hostname: 'myhost',
                                    tempdir:  '/tmp/mydir',
                                    port:     9999,
                                    debug:    true,
                                    help:     true,
                                    anyorigin: false,
     				   originlist: ['some', 'urls']
                    };});
     		  })
       .toThrowError("'originlist' should be a comma-separated string");

    expect( () => { config.provide( () => {return {mongourl: 'mymongourl',
                                   hostname: 'myhost',
                                   tempdir:  '/tmp/mydir',
                                   port:     9999,
                                   debug:    true,
                                   help:     true,
                                   anyorigin: false,
    				   originlist: '*'
                    };});
    		  }).toThrowError(
                    "'originlist' string should not contain the wild card character *");

    expect( () => { config.provide( () => {return {mongourl: 'mymongourl',
                                   hostname: 'myhost',
                                   tempdir:  '/tmp/mydir',
    				   port:     9999,
                                   debug:    true,
                                   help:     true,
                                   anyorigin: false,
    				   originlist: 'myurl'
                    };});
    		  }).toThrowError(/Protocol is absent/);

    expect( () => { config.provide( () => {return {mongourl: 'mymongourl',
                                   hostname: 'myhost',
                                   tempdir:  '/tmp/mydir',
    				   port:     9999,
                                   debug:    true,
                                   help:     true,
                                   anyorigin: false,
    				   originlist: 'myurl'
                    };});
    		  }).toThrowError(/Protocol is absent/);

    expect( () => { config.provide( () => {return {mongourl: 'mymongourl',
                                   hostname: 'myhost',
                                   tempdir:  '/tmp/mydir',
    				   port:     9999,
                                   debug:    true,
                                   help:     true,
                                   anyorigin: false,
    				   originlist: 'https://myurl'
                    };});
    		  }).toThrowError(/URL protocol should match server protocol/);

    expect( () => { config.provide( () => {return {mongourl: 'mymongourl',
                                   hostname: 'myhost',
                                   tempdir:  '/tmp/mydir',
    				   port:     9999,
                                   debug:    true,
                                   help:     true,
                                   anyorigin: false,
    				   originlist: 'http://'
                    };});
    		  }).toThrowError(/Server host is absent/);

    expect( () => { config.provide( () => {return {mongourl: 'mymongourl',
                                   hostname: 'myhost',
                                   tempdir:  '/tmp/mydir',
    				   port:     9999,
                                   debug:    true,
                                   help:     true,
                                   anyorigin: false,
    				   originlist: 'http://server.com/foo'
                    };});
    		  }).toThrowError(/Path cannot be present/);

    expect( () => { config.provide( () => {return {mongourl: 'mymongourl',
                                   hostname: 'myhost',
                                   tempdir:  '/tmp/mydir',
    				   port:     9999,
                                   debug:    true,
                                   help:     true,
                                   anyorigin: false,
    				   originlist: 'http://server.com/'
                    };});
    		  }).not.toThrow();

    expect( () => { config.provide( () => {return {mongourl: 'mymongourl',
                                   hostname: 'myhost',
                                   tempdir:  '/tmp/mydir',
    				   port:     9999,
                                   debug:    true,
                                   help:     true,
                                   anyorigin:  false,
    				   originlist: 'http://server.com/,,http://server.org'
                    };});
    		  }).toThrowError(/Empty string in 'originlist'/);

    expect( () => { config.provide( () => {return {mongourl: 'mymongourl',
                                   hostname: 'myhost',
                                   tempdir:  '/tmp/mydir',
    				   port:     9999,
                                   debug:    true,
                                   help:     true,
                                   anyorigin: false,
    				   originlist: 'http://server.com:8080?foo=2'
                    };});
    		  }).toThrowError(/Search string cannot be present/);

   expect( () => { config.provide( () => {return {mongourl: 'mymongourl',
                                   hostname: 'myhost',
                                   tempdir:  '/tmp/mydir',
    				   port:     9999,
                                   debug:    true,
                                   help:     true,
                                   anyorigin: false,
    				   originlist: 'http://server.com:8080?'
                    };});
    		  }).toThrowError(/Search string cannot be present/);

    expect( () => { config.provide( () => {return {mongourl: 'mymongourl',
                                   hostname: 'myhost',
                                   tempdir:  '/tmp/mydir',
    				   port:     9999,
                                   debug:    true,
                                   help:     true,
                                   anyorigin: false,
    				   originlist: 'http://localhost:9999#mttag'
                    };});
    		  }).toThrowError(/Hash tag cannot be present/);

    expect( () => { config.provide( () => {return {mongourl: 'mymongourl',
                                   hostname: 'myhost',
                                   tempdir:  '/tmp/mydir',
    				   port:     9999,
                                   debug:    true,
                                   help:     true,
                                   anyorigin: false,
    				   originlist: 'http://server.com:9999#'
                    };});
    		  }).toThrowError(/Hash tag cannot be present/);
  });

  it('Setting originlist to an array of URL strings', function() {
    config.provide( () => {return {mongourl: 'mymongourl',
                                   hostname: 'myhost',
                                   tempdir:  '/tmp/mydir',
                                   port:     9999,
                                   debug:    true,
                                   help:     true,
                                   originlist: 'http://my.com'
                                  };});
    expect(config.provide().get('originlist').join()).toEqual(
      'http://my.com', 'one url in the array');

    let expected = 'http://my.com:80,http://your.org';
    config.provide( () => {return {mongourl: 'mymongourl',
                                   hostname: 'myhost',
                                   tempdir:  '/tmp/mydir',
                                   port:     9999,
                                   debug:    true,
                                   help:     true,
                                   originlist: 'http://my.com:80,http://your.org'
                                  };});
    expect(config.provide().get('originlist').join()).toEqual(expected, 'two urls in the array');

    config.provide( () => {return {mongourl: 'mymongourl',
                                   hostname: 'myhost',
                                   tempdir:  '/tmp/mydir',
                                   port:     9999,
                                   debug:    true,
                                   help:     true,
                                   originlist: 'http://my.com:80, http://your.org'
                                  };});
    expect(config.provide().get('originlist').join()).toEqual(
      expected, 'spaces between strings are allowed');
  });
});

