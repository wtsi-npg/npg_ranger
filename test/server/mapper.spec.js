/* globals describe, it, fail, expect, beforeAll, afterAll */

"use strict";

const assert      = require('assert');
const child       = require('child_process');
const MongoClient = require('mongodb').MongoClient;
const tmp         = require('tmp');
const fse         = require('fs-extra');

const DataMapper  = require('../../lib/server/mapper.js');
const config      = require('../../lib/config.js');

const BASE_PORT  = 1400;
const PORT_RANGE = 200;
const PORT       = Math.floor(Math.random() * PORT_RANGE) + BASE_PORT;
const FIXTURES   = 'test/server/data/fixtures/fileinfo.json';

/* ***************************************************************************
 * Test cases (sample accession number, file path, details)
 * ***************************************************************************
 * XYZ120923 /seq/10000/10000_1#58_phix.bam  mqc 0
 *           /seq/10000/10000_1#63.bam       not a target
 *           /seq/10000/10000_2#24.bam       alignment 0
 *           /seq/10000/10000_2#22.bam       OK
 *
 * XYZ238967 /seq/10000/10000_4#43.bam  OK
 *           /seq/10000/10000_5#74.bam  access_control_group_id - empty string
 *           /seq/10000/10000_6#83.bam  access_control_group_id - no key
 *           /seq/10000/10000_7#92.bam  access_control_group_id - zero
 * ***************************************************************************
 */

describe('Static function', function() {
  it('getFilterNames', function() {
    let filterNames = DataMapper.getFilterNames();
    let expectedFilterNames = [
      'target', 'target_not',
      'alignment', 'alignment_not',
      'manual_qc', 'manual_qc_not',
      'alignment_filter', 'alignment_filter_not',
      'alt_target', 'alt_target_not',
      'alt_process', 'alt_process_not'
    ];
    expect(filterNames.length).toBe(12);
    expectedFilterNames.forEach( (filter) => {
      expect(filterNames).toContain(filter);
    });
  });
});

describe('Data info retrieval', function() {
  var tmpobj = tmp.dirSync({ prefix: 'npg_ranger_test_' });
  var tmp_dir = tmpobj.name;
  console.log(`MONGO data directory: ${tmp_dir}`);
  var db_name = 'npg_ranger_test';
  var url = `mongodb://localhost:${PORT}/${db_name}`;

  var compareFiles = function(a, b) {
    if (a.file > b.file) {
      return 1;
    }
    if (a.file < b.file) {
      return -1;
    }
    return 0;
  };

  beforeAll( () => {
    let command = `mongod -f test/server/data/mongodb_conf.yml --port ${PORT} --dbpath ${tmp_dir} --pidfilepath ${tmp_dir}/mpid --logpath ${tmp_dir}/dbserver.log`;
    console.log(`\nCommand to start MONGO DB daemon: ${command}`);
    let out = child.execSync(command);
    console.log(`Started MONGO DB daemon: ${out}`);
    command = `mongoimport --port ${PORT} --db ${db_name} --collection fileinfo --jsonArray --file ${FIXTURES}`;
    out = child.execSync(command);
    console.log(`Loaded data to MONGO DB: ${out}`);
    config.provide(() => { return {multiref: true}; });
  });


  describe('Reference mismatch correctly handled', function() {
    beforeAll( () => {
      config.provide(() => { return {}; });
    });

    afterAll( () => {
      config.provide(() => { return {multiref: true}; });
    });

    it('Do not allow mismatching references when multiref falsy', function(done) {
      MongoClient.connect(url, function(err, db) {
        assert.equal(err, null);
        var dm = new DataMapper(db);
        dm.once('data', () => {
          fail();
          done();
        });
        dm.once('nodata', (reason) => {
          expect(reason).toBe('Not all references match for sample accession XYZ238967');
          done();
        });
        dm.getFileInfo({accession: "XYZ238967"}, 'localhost');
      });
    });

    it('Fail if no reference found', function(done) {
      MongoClient.connect(url, function(err, db) {
        assert.equal(err, null);
        var dm = new DataMapper(db);
        dm.once('data', () => {
          fail();
          done();
        });
        dm.once('nodata', (reason) => {
          expect(reason).toBe('No reference for 10000_1#58_phix.bam');
          done();
        });
        dm.getFileInfo({name: "10000_1#58_phix.bam"}, 'localhost');
      });
    });
  });

  describe('Test filters', function() {
    it('invalid combination', function(done) {
      MongoClient.connect(url, function(err, db) {
        assert.equal(err, null);
        let dm = new DataMapper(db);
        dm.once('nodata', (reason) => {
          // Expect failure because both target and target_not are specified,
          // which is not allowed
          expect(reason).toBe('Invalid query');
          done();
        });
        dm.on('data', (data) => {
          fail('Should have received no data, but instead received ' + data.toString());
          done();
        });
        dm.getFileInfo({accession: 'JKL987654', target: '1', target_not: '0'}, 'localhost');
      });
    });

    it('all defaults', function(done) {
      MongoClient.connect(url, function(err, db) {
        assert.equal(err, null);
        let dm = new DataMapper(db);
        dm.once('nodata', (reason) => {
          fail(reason);
          done();
        });
        dm.on('data', (data) => {
          expect(data).toEqual(
            [{file: 'target1 mqc1 align1', accessGroup : '2136',
              reference : '/test/reference/path.fa'}]
          );
          done();
        });
        dm.getFileInfo({accession: 'JKL987654'}, 'localhost');
      });
    });

    it('target=0', function(done) {
      MongoClient.connect(url, function(err, db) {
        assert.equal(err, null);
        let dm = new DataMapper(db);
        dm.once('nodata', (reason) => {
          fail(reason);
          done();
        });
        dm.once('data', (data) => {
          expect(data.length).toBe(2);
          expect(data).toContain(
            {file: 'target0 mqc1 align1', accessGroup : '2136',
              reference : '/test/reference/path.fa'}
          );
          expect(data).toContain(
            {file: 'target0 mqc1 align1 afphix', accessGroup : '2136',
              reference : '/test/reference/path.fa'}
          );
          done();
        });
        dm.getFileInfo({accession: 'JKL987654', target: '0'}, 'localhost');
      });
    });

    it('target=', function(done) {
      MongoClient.connect(url, function(err, db) {
        assert.equal(err, null);
        let dm = new DataMapper(db);
        dm.once('nodata', (reason) => {
          fail(reason);
          done();
        });
        dm.once('data', (data) => {
          expect(data.length).toBe(4);
          expect(data).toContain(
            {file: 'target0 mqc1 align1', accessGroup : '2136',
              reference : '/test/reference/path.fa'}
          );
          expect(data).toContain(
            {file: 'target1 mqc1 align1', accessGroup : '2136',
              reference : '/test/reference/path.fa'}
          );
          expect(data).toContain(
            {file: 'targetu mqc1 align1', accessGroup : '2136',
              reference : '/test/reference/path.fa'}
          );
          expect(data).toContain(
            {file: 'target0 mqc1 align1 afphix', accessGroup : '2136',
              reference : '/test/reference/path.fa'}
          );
          done();
        });
        dm.getFileInfo({accession: 'JKL987654', target: ''}, 'localhost');
      });
    });

    it('target=undef', function(done) {
      MongoClient.connect(url, function(err, db) {
        assert.equal(err, null);
        let dm = new DataMapper(db);
        dm.once('nodata', (reason) => {
          fail(reason);
          done();
        });
        dm.once('data', (data) => {
          expect(data.length).toBe(1);
          expect(data).toContain(
            {file: 'targetu mqc1 align1', accessGroup : '2136',
              reference : '/test/reference/path.fa'}
          );
          done();
        });
        dm.getFileInfo({accession: 'JKL987654', target: 'undef'}, 'localhost');
      });
    });

    it('target_not=1', function(done) {
      MongoClient.connect(url, function(err, db) {
        assert.equal(err, null);
        let dm = new DataMapper(db);
        dm.once('nodata', (reason) => {
          fail(reason);
          done();
        });
        dm.once('data', (data) => {
          expect(data.length).toBe(3);
          expect(data).toContain(
            {file: 'targetu mqc1 align1', accessGroup : '2136',
              reference : '/test/reference/path.fa'}
          );
          expect(data).toContain(
            {file: 'target0 mqc1 align1', accessGroup : '2136',
              reference : '/test/reference/path.fa'}
          );
          expect(data).toContain(
            {file: 'target0 mqc1 align1 afphix', accessGroup : '2136',
              reference : '/test/reference/path.fa'}
          );
          done();
        });
        dm.getFileInfo({accession: 'JKL987654', target_not: '1'}, 'localhost');
      });
    });

    it('target_not=undef', function(done) {
      MongoClient.connect(url, function(err, db) {
        assert.equal(err, null);
        let dm = new DataMapper(db);
        dm.once('nodata', (reason) => {
          fail(reason);
          done();
        });
        dm.once('data', (data) => {
          expect(data.length).toBe(3);
          expect(data).toContain(
            {file: 'target1 mqc1 align1', accessGroup : '2136',
              reference : '/test/reference/path.fa'}
          );
          expect(data).toContain(
            {file: 'target0 mqc1 align1', accessGroup : '2136',
              reference : '/test/reference/path.fa'}
          );
          expect(data).toContain(
            {file: 'target0 mqc1 align1 afphix', accessGroup : '2136',
              reference : '/test/reference/path.fa'}
          );
          done();
        });
        dm.getFileInfo({accession: 'JKL987654', target_not: 'undef'}, 'localhost');
      });
    });

    it('target_not=', function(done) {
      MongoClient.connect(url, function(err, db) {
        assert.equal(err, null);
        let dm = new DataMapper(db);
        dm.once('nodata', (reason) => {
          fail(reason);
          done();
        });
        dm.once('data', (data) => {
          expect(data.length).toBe(4);
          expect(data).toContain(
            {file: 'target0 mqc1 align1', accessGroup : '2136',
              reference : '/test/reference/path.fa'}
          );
          expect(data).toContain(
            {file: 'target1 mqc1 align1', accessGroup : '2136',
              reference : '/test/reference/path.fa'}
          );
          expect(data).toContain(
            {file: 'targetu mqc1 align1', accessGroup : '2136',
              reference : '/test/reference/path.fa'}
          );
          expect(data).toContain(
            {file: 'target0 mqc1 align1 afphix', accessGroup : '2136',
              reference : '/test/reference/path.fa'}
          );
          done();
        });
        dm.getFileInfo({accession: 'JKL987654', target_not: ''}, 'localhost');
      });
    });

    it('alignment_filter=1', function(done) {
      MongoClient.connect(url, function(err, db) {
        assert.equal(err, null);
        let dm = new DataMapper(db);
        dm.once('nodata', (reason) => {
          // Expect no data because nothing has alignment_filter=1
          expect(reason).toBe('No files for sample accession JKL987654');
          done();
        });
        dm.once('data', (data) => {
          fail('should return no data but returned ' + data.toString());
          done();
        });
        dm.getFileInfo({accession: 'JKL987654', alignment_filter: '1'}, 'localhost');
      });
    });

    it('alignment_filter=phix', function(done) {
      MongoClient.connect(url, function(err, db) {
        assert.equal(err, null);
        let dm = new DataMapper(db);
        dm.once('nodata', (reason) => {
          // Expect no data because only file with alignment_filter=phix
          // also has target=1, which must be explicitly defined
          expect(reason).toBe('No files for sample accession JKL987654');
          done();
        });
        dm.once('data', (data) => {
          fail('should return no data but returned ' + data.toString());
          done();
        });
        dm.getFileInfo({accession: 'JKL987654', alignment_filter: 'phix'}, 'localhost');
      });
    });

    it('target=0, alignment_filter=phix', function(done) {
      MongoClient.connect(url, function(err, db) {
        assert.equal(err, null);
        let dm = new DataMapper(db);
        dm.once('nodata', (reason) => {
          fail(reason);
          done();
        });
        dm.once('data', (data) => {
          expect(data.length).toBe(1);
          expect(data).toEqual(
            [{file: 'target0 mqc1 align1 afphix', accessGroup : '2136',
              reference : '/test/reference/path.fa'}]
          );
          done();
        });
        dm.getFileInfo({accession: 'JKL987654', target: '0', alignment_filter: 'phix'}, 'localhost');
      });
    });

    it('target=0, alignment_filter_not=phix', function(done) {
      MongoClient.connect(url, function(err, db) {
        assert.equal(err, null);
        let dm = new DataMapper(db);
        dm.once('nodata', (reason) => {
          fail(reason);
          done();
        });
        dm.once('data', (data) => {
          expect(data.length).toBe(1);
          expect(data).toEqual(
            [{file: 'target0 mqc1 align1', accessGroup : '2136',
              reference : '/test/reference/path.fa'}]
          );
          done();
        });
        dm.getFileInfo({accession: 'JKL987654', target: '0', alignment_filter_not: 'phix'}, 'localhost');
      });
    });

    it('target=0, alignment_filter=', function(done) {
      MongoClient.connect(url, function(err, db) {
        assert.equal(err, null);
        let dm = new DataMapper(db);
        dm.once('nodata', (reason) => {
          fail(reason);
          done();
        });
        dm.once('data', (data) => {
          expect(data.length).toBe(2);
          expect(data).toContain(
            {file: 'target0 mqc1 align1', accessGroup : '2136',
              reference : '/test/reference/path.fa'}
          );
          expect(data).toContain(
            {file: 'target0 mqc1 align1 afphix', accessGroup : '2136',
              reference : '/test/reference/path.fa'}
          );
          done();
        });
        dm.getFileInfo({accession: 'JKL987654', target: '0', alignment_filter: ''}, 'localhost');
      });
    });


  });


  it('Input validation', function() {
    expect( () => {new DataMapper();} ).toThrowError(ReferenceError,
      'Database handle is required');
    expect( () => {
      let dm = new DataMapper({});
      dm.getFileInfo();
    }).toThrowError(ReferenceError, 'Query object is required');
    expect( () => {
      let dm = new DataMapper({});
      dm.getFileInfo({});
    }).toThrowError(ReferenceError, 'Host name is required');
    expect( () => {
      let dm = new DataMapper({});
      dm.getFileInfo({fruit: "apple"}, 'localhost');
    }).toThrowError(Error,
      'Sample accession number or file name should be given');
  });

  it('No data received - unknown accession number given', function(done) {
    MongoClient.connect(url, function(err, db) {
      assert.equal(err, null);
      var dm = new DataMapper(db);
      dm.on('nodata', (reason) => {
        expect(reason).toBe('No files for sample accession KRT1234');
        dm.removeAllListeners();
        done();
      });
      dm.getFileInfo({accession: "KRT1234"}, 'localhost');
    });
  });

  it('No data received - unknown file name', function(done) {
    MongoClient.connect(url, function(err, db) {
      assert.equal(err, null);
      var dm = new DataMapper(db);
      dm.on('nodata', (reason) => {
        expect(reason).toBe('No files for 1234_2.bam');
        done();
      });
      dm.getFileInfo({name: "1234_2.bam"}, 'localhost');
    });
  });

  it('Filtering on flags - some data left', function(done) {
    MongoClient.connect(url, function(err, db) {
      assert.equal(err, null);
      var dm = new DataMapper(db);
      dm.on('data', (data) => {
        expect(data).toEqual(
          [{file: 'irods:/seq/10000/10000_2#22.bam', accessGroup: '2574',
            reference: '/Homo_sapiens/1000Genomes_hs37d5/all/fasta/hs37d5.fa'}]);
        done();
      });
      dm.getFileInfo({accession: "XYZ120923"}, 'localhost');
    });
  });

  it('No Filtering on flags when querying by name', function(done) {
    MongoClient.connect(url, function(err, db) {
      assert.equal(err, null);
      var dm = new DataMapper(db);
      dm.on('data', (data) => {
        expect(data).toEqual(
          [{file: 'irods:/seq/10000/10000_1#63.bam', accessGroup: '2136',
            reference: '/Homo_sapiens/1000Genomes_hs37d5/all/fasta/hs37d5.fa'}]);
        done();
      });
      dm.getFileInfo({name: "10000_1#63.bam"}, 'localhost');
    });
  });

  it('Query by name', function(done) {
    MongoClient.connect(url, function(err, db) {
      assert.equal(err, null);
      var dm = new DataMapper(db);
      dm.on('data', (data) => {
        expect(data).toEqual(
          [{file: 'irods:/seq/10000/10000_2#22.bam', accessGroup: '2574',
            reference: '/Homo_sapiens/1000Genomes_hs37d5/all/fasta/hs37d5.fa'}]);
        done();
      });
      dm.getFileInfo({name: "10000_2#22.bam"}, 'localhost');
    });
  });

  it('No data received - unknown directory', function(done) {
    MongoClient.connect(url, function(err, db) {
      assert.equal(err, null);
      var dm = new DataMapper(db);
      dm.on('nodata', (reason) => {
        expect(reason).toBe('No files for 10000_2#22.bam in /tmp/files');
        done();
      });
      dm.getFileInfo({name: "10000_2#22.bam", directory: "/tmp/files"}, 'localhost');
    });
  });

  it('Query by name and directory', function(done) {
    MongoClient.connect(url, function(err, db) {
      assert.equal(err, null);
      var dm = new DataMapper(db);
      dm.on('data', (data) => {
        expect(data).toEqual([{
          file: 'irods:/seq/10000/10000_2#22.bam', accessGroup: '2574',
          reference: '/Homo_sapiens/1000Genomes_hs37d5/all/fasta/hs37d5.fa'
        }]);
        done();
      });
      dm.getFileInfo({name: "10000_2#22.bam", directory: "/seq/10000"}, 'localhost');
    });
  });

  it('Query by name, localisation by host', function(done) {
    MongoClient.connect(url, function(err, db) {
      assert.equal(err, null);
      var dm = new DataMapper(db);
      dm.on('data', (data) => {
        expect(data).toEqual(
          [{ file: '/irods-seq-i10-bc/seq/10000/10000_2#22.bam', accessGroup: '2574',
             reference: '/Homo_sapiens/1000Genomes_hs37d5/all/fasta/hs37d5.fa'}]);
        done();
      });
      dm.getFileInfo({name: "10000_2#22.bam"}, 'irods-seq-i10');
    });
  });

  it('Query by accession number, multiple results', function(done) {
    MongoClient.connect(url, function(err, db) {
      assert.equal(err, null);
      var dm = new DataMapper(db);
      dm.on('data', (data) => {
        let d = [
          {
            file: 'irods:/seq/10000/10000_4#43.bam', accessGroup: '2586',
            reference: '/Homo_sapiens/1000Genomes_hs37d5/all/fasta/hs37d5.fa'
          },
          {
            file: 'irods:/seq/10000/10000_5#74.bam', accessGroup: '',
            reference: '/Caenorhabditis_elegans/101019/all/fasta/C_elegans_101019.fasta'
          },
          {
            file: 'irods:/seq/10000/10000_7#92.bam', accessGroup: '0',
            reference: '/Homo_sapiens/1000Genomes_hs37d5/all/fasta/hs37d5.fa'
          },
        ];
        data.sort(compareFiles);
        expect(data).toEqual(d);
        done();
      });
      dm.getFileInfo({accession: "XYZ238967"}, 'localhost');
    });
  });

  it('Query by accession number, multiple results, localisation by host', function(done) {
    MongoClient.connect(url, function(err, db) {
      assert.equal(err, null);
      var dm = new DataMapper(db);
      dm.on('data', (data) => {
        let d = [
          {file: '/irods-seq-sr04-ddn-gc10-30-31-32/seq/10000/10000_4#43.bam', accessGroup: '2586', reference: '/Homo_sapiens/1000Genomes_hs37d5/all/fasta/hs37d5.fa'},
          {file: '/irods-seq-sr04-ddn-gc10-30-31-32/seq/10000/10000_5#74.bam', accessGroup: '', reference: '/Caenorhabditis_elegans/101019/all/fasta/C_elegans_101019.fasta'},
          {file: '/irods-seq-sr04-ddn-gc10-30-31-32/seq/10000/10000_7#92.bam', accessGroup: '0', reference: '/Homo_sapiens/1000Genomes_hs37d5/all/fasta/hs37d5.fa'},
          {file: '/irods-seq-sr04-ddn-gc10-30-31-32/seq/10000/10000_8#97.bam', accessGroup: '2574', reference: '/Homo_sapiens/1000Genomes_hs37d5/all/fasta/hs37d5.fa'}
        ];
        data.sort(compareFiles);
        expect(data).toEqual(d);
        done();
      });
      dm.getFileInfo({accession: "XYZ238967"}, 'irods-seq-sr04');
    });
  });

  it('Some files are local, some not', function(done) {
    MongoClient.connect(url, function(err, db) {
      assert.equal(err, null);
      var dm = new DataMapper(db);
      dm.on('data', (data) => {
        let d = [
          {file: 'irods:/seq/10000/10000_4#43.bam', accessGroup: '2586',
           reference: '/Homo_sapiens/1000Genomes_hs37d5/all/fasta/hs37d5.fa'},
          {file: 'irods:/seq/10000/10000_5#74.bam', accessGroup: '',
           reference: '/Caenorhabditis_elegans/101019/all/fasta/C_elegans_101019.fasta'},
          {file: '/irods-seq-i10-bc/seq/10000/10000_7#92.bam', accessGroup: '0',
           reference: '/Homo_sapiens/1000Genomes_hs37d5/all/fasta/hs37d5.fa'},
        ];
        d.sort(compareFiles);
        data.sort(compareFiles);
        expect(data).toEqual(d);
        done();
      });
      dm.getFileInfo({accession: "XYZ238967"}, 'irods-seq-i10');
    });
  });

  it('Filtering out empty paths, no results remaining', function(done) {
    MongoClient.connect(url, function(err, db) {
      assert.equal(err, null);
      var dm = new DataMapper(db);
      dm.on('nodata', (reason) => {
        expect(reason).toBe('No files for 10000_8#97.bam');
        done();
      });
      dm.getFileInfo({name: "10000_8#97.bam"}, 'irods-seq-i10');
    });
  });

  afterAll(() => {
    child.execSync(`mongo 'mongodb://localhost:${PORT}/admin' --eval 'db.shutdownServer()'`);
    console.log('\nMONGODB server has been shut down');
    fse.remove(tmp_dir, (err) => {if (err) {console.log(`Error removing ${tmp_dir}: ${err}`);}});
    config.provide(() => { return {}; });
  });
});
