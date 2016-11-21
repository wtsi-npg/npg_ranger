/* globals describe, it, expect */

"use strict";

const accUtils = require('../../lib/uiclient/accUtils');

describe('testing buildIndexPath', () => {
  it('returns expected index path', () => {
    let expected = {
      'references/Homo_sapiens/1000Genomes/all/fasta/human.fasta':
        'references/Homo_sapiens/1000Genomes/all/samtools/human.fasta.fai',
      'references/Human/1000Genomes/all/fasta/human.fasta':
        'references/Human/1000Genomes/all/samtools/human.fasta.fai',
      'references/references/wrong/path.fasta': null,
      'references/Homo_sapiens/1000Genomes/all/human.fasta': null,
      '1000Genomes/all/fasta/human.fasta': null
    };
    Object.keys(expected).forEach((referencePath) => {
      expect(accUtils.buildIndexPath(referencePath)).toEqual(
        expected[referencePath]
      );
    });
  });
});

describe('testing getReferenceInfo', () => {
  it('validates parameters', () => {
    expect( () => {
      accUtils.getReferenceInfo();
    }).toThrowError(/referencePath is required/);
  });

  it('validates reference path structure', () => {
    [
      'unexpected/format/for/reference.fast',
      ' ',
      'singleword'
    ].forEach( ( referencePath ) => {
      expect( () => {
        accUtils.getReferenceInfo( referencePath );
      }).toThrowError(/Unexpected reference path structure/);
    });
  });

  it('returns expected reference info', () => {
    let expected = {
      'references/Homo_sapiens/1000Genomes/all/fasta/human.fasta': {
        speciesName: 'Homo sapiens',
        referenceName: '1000Genomes',
        twoBitPath: 'references/Homo_sapiens/1000Genomes/all/blat/human.fasta.2bit'
      },
      'references/Human/1000Genomes/all/fasta/human.fasta': {
        speciesName: 'Human',
        referenceName: '1000Genomes',
        twoBitPath: 'references/Human/1000Genomes/all/blat/human.fasta.2bit'
      },
      'references/Mouse/CGP_NCBI37/all/fasta/mm_ref_NCBI37_1.fasta': {
        speciesName: 'Mouse',
        referenceName: 'CGP NCBI37',
        twoBitPath: 'references/Mouse/CGP_NCBI37/all/blat/mm_ref_NCBI37_1.fasta.2bit'
      },
      'references/Mus_musculus/CGP_NCBI37/all/fasta/mm_ref_NCBI37_1.fasta': {
        speciesName: 'Mus musculus',
        referenceName: 'CGP NCBI37',
        twoBitPath: 'references/Mus_musculus/CGP_NCBI37/all/blat/mm_ref_NCBI37_1.fasta.2bit'
      }
    };
    Object.keys(expected).forEach((referencePath) => {
      expect(accUtils.getReferenceInfo(referencePath)).toEqual(
        expected[referencePath]
      );
    });
  });
});

describe('testing getRange', () => {
  it('validates parameters', () => {
    expect(() => {
      accUtils.getRange();
    }).toThrowError(/chrData is required/);
    expect(() => {
      accUtils.getRange([]);
    }).toThrowError(/chr is required/);
    expect(() => {
      accUtils.getRange([], 'chr1');
    }).toThrowError(/rangeSize is required/);
  });

  it('can find the expected ranges', () => {
    let chrData = [
      ['1', '500'],
      ['chr2', '2000'],
      ['3', '5000'],
      ['4', '50'],
      ['5', '0']
    ];

    expect(accUtils.getRange(chrData, '1', 100)).toEqual(
      { from: 200, to: 300 }
    );

    expect(accUtils.getRange(chrData, 'chr2', 100)).toEqual(
      { from: 950, to: 1050 }
    );

    expect(accUtils.getRange(chrData, '4', 100)).toEqual(
      { from: 0, to: 50 }
    );

    expect(accUtils.getRange(chrData, '5', 100)).toEqual(
      { from: 0, to: 0 }
    );
  });
});

describe('testing parseIndexData', () => {
  it('validates parameters', () => {
    expect(() => {
      accUtils.parseIndexData();
    }).toThrowError(/indexData is required/);
  });

  it('validates parameter type', () =>  {
    expect(() => {
      accUtils.parseIndexData(1);
    }).toThrowError(/must be a string/);
  });

  let indexData = `1	249	52	60	61
2	243	253404903	60	61
X	155	2929051733	60	61
Y	593	3086910193	60	61
GL000207.1	426	3147290265	60	61
GL000192.1	547	3152949898	60	61`;
  let chrs     = '1 2 X Y GL000207.1 GL000192.1'.split(' ');
  let chrSizes = '249 243 155 593 426 547'.split(' ');

  it('can parse index data', () => {
    let parsedData = accUtils.parseIndexData( indexData );
    expect(parsedData.length).toBe(6);
    parsedData.forEach(( row, index ) => {
      expect(row.length).toBe(2);
      expect(row[0]).toEqual(chrs[index]);
      expect(row[1]).toEqual(chrSizes[index]);
    });
  });
});
