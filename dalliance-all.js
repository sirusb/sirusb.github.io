/* -*- mode: javascript; c-basic-offset: 4; indent-tabs-mode: nil -*- */

// 
// Dalliance Genome Explorer
// (c) Thomas Down 2006-2011
//
// bam.js: indexed binary alignments
//

"use strict";

var BAM_MAGIC = 21840194;
var BAI_MAGIC = 21578050;

function BamFile() {
}

function Vob(b, o) {
    this.block = b;
    this.offset = o;
}

Vob.prototype.toString = function() {
    return '' + this.block + ':' + this.offset;
}

function Chunk(minv, maxv) {
    this.minv = minv; this.maxv = maxv;
}

function makeBam(data, bai, callback) {
    var bam = new BamFile();
    bam.data = data;
    bam.bai = bai;

    bam.bai.fetch(function(header) {   // Do we really need to fetch the whole thing? :-(
        if (!header) {
            return callback(null, "Couldn't access BAI");
        }

        var uncba = new Uint8Array(header);
        var baiMagic = readInt(uncba, 0);
        if (baiMagic != BAI_MAGIC) {
            return callback(null, 'Not a BAI file, magic=0x' + baiMagic.toString(16));
        }

        var nref = readInt(uncba, 4);

        bam.indices = [];

        var p = 8;
        var minBlockIndex = 1000000000;
        for (var ref = 0; ref < nref; ++ref) {
            var blockStart = p;
            var nbin = readInt(uncba, p); p += 4;
            for (var b = 0; b < nbin; ++b) {
                var bin = readInt(uncba, p);
                var nchnk = readInt(uncba, p+4);
                p += 8 + (nchnk * 16);
            }
            var nintv = readInt(uncba, p); p += 4;
            
            var q = p;
            for (var i = 0; i < nintv; ++i) {
                var v = readVob(uncba, q); q += 8;
                if (v) {
                    var bi = v.block;
                    if (v.offset > 0)
                        bi += 65536;

                    if (bi < minBlockIndex)
                        minBlockIndex = bi;
                    break;
                }
            }
            p += (nintv * 8);


            if (nbin > 0) {
                bam.indices[ref] = new Uint8Array(header, blockStart, p - blockStart);
            }                     
        }

        bam.data.slice(0, minBlockIndex).fetch(function(r) {
            if (!r) {
                return callback(null, "Couldn't access BAM");
            }
            
            var unc = unbgzf(r, r.byteLength);
            var uncba = new Uint8Array(unc);

            var magic = readInt(uncba, 0);
            if (magic != BAM_MAGIC) {
                return callback(null, "Not a BAM file, magic=0x" + magic.toString(16));
            }
            var headLen = readInt(uncba, 4);
            var header = '';
            for (var i = 0; i < headLen; ++i) {
                header += String.fromCharCode(uncba[i + 8]);
            }

            var nRef = readInt(uncba, headLen + 8);
            var p = headLen + 12;

            bam.chrToIndex = {};
            bam.indexToChr = [];
            for (var i = 0; i < nRef; ++i) {
                var lName = readInt(uncba, p);
                var name = '';
                for (var j = 0; j < lName-1; ++j) {
                    name += String.fromCharCode(uncba[p + 4 + j]);
                }
                var lRef = readInt(uncba, p + lName + 4);
                bam.chrToIndex[name] = i;
                if (name.indexOf('chr') == 0) {
                    bam.chrToIndex[name.substring(3)] = i;
                } else {
                    bam.chrToIndex['chr' + name] = i;
                }
                bam.indexToChr.push(name);

                p = p + 8 + lName;
            }

            if (bam.indices) {
                return callback(bam);
            }
        });
    });
}



BamFile.prototype.blocksForRange = function(refId, min, max) {
    var index = this.indices[refId];
    if (!index) {
        return [];
    }

    var intBinsL = reg2bins(min, max);
    var intBins = [];
    for (var i = 0; i < intBinsL.length; ++i) {
        intBins[intBinsL[i]] = true;
    }
    var leafChunks = [], otherChunks = [];

    var nbin = readInt(index, 0);
    var p = 4;
    for (var b = 0; b < nbin; ++b) {
        var bin = readInt(index, p);
        var nchnk = readInt(index, p+4);
//        dlog('bin=' + bin + '; nchnk=' + nchnk);
        p += 8;
        if (intBins[bin]) {
            for (var c = 0; c < nchnk; ++c) {
                var cs = readVob(index, p);
                var ce = readVob(index, p + 8);
                (bin < 4681 ? otherChunks : leafChunks).push(new Chunk(cs, ce));
                p += 16;
            }
        } else {
            p +=  (nchnk * 16);
        }
    }
//    dlog('leafChunks = ' + miniJSONify(leafChunks));
//    dlog('otherChunks = ' + miniJSONify(otherChunks));

    var nintv = readInt(index, p);
    var lowest = null;
    var minLin = Math.min(min>>14, nintv - 1), maxLin = Math.min(max>>14, nintv - 1);
    for (var i = minLin; i <= maxLin; ++i) {
        var lb =  readVob(index, p + 4 + (i * 8));
        if (!lb) {
            continue;
        }
        if (!lowest || lb.block < lowest.block || lb.offset < lowest.offset) {
            lowest = lb;
        }
    }
    // dlog('Lowest LB = ' + lowest);
    
    var prunedOtherChunks = [];
    if (lowest != null) {
        for (var i = 0; i < otherChunks.length; ++i) {
            var chnk = otherChunks[i];
            if (chnk.maxv.block >= lowest.block && chnk.maxv.offset >= lowest.offset) {
                prunedOtherChunks.push(chnk);
            }
        }
    }
    // dlog('prunedOtherChunks = ' + miniJSONify(prunedOtherChunks));
    otherChunks = prunedOtherChunks;

    var intChunks = [];
    for (var i = 0; i < otherChunks.length; ++i) {
        intChunks.push(otherChunks[i]);
    }
    for (var i = 0; i < leafChunks.length; ++i) {
        intChunks.push(leafChunks[i]);
    }

    intChunks.sort(function(c0, c1) {
        var dif = c0.minv.block - c1.minv.block;
        if (dif != 0) {
            return dif;
        } else {
            return c0.minv.offset - c1.minv.offset;
        }
    });
    var mergedChunks = [];
    if (intChunks.length > 0) {
        var cur = intChunks[0];
        for (var i = 1; i < intChunks.length; ++i) {
            var nc = intChunks[i];
            if (nc.minv.block == cur.maxv.block /* && nc.minv.offset == cur.maxv.offset */) { // no point splitting mid-block
                cur = new Chunk(cur.minv, nc.maxv);
            } else {
                mergedChunks.push(cur);
                cur = nc;
            }
        }
        mergedChunks.push(cur);
    }
    // dlog('mergedChunks = ' + miniJSONify(mergedChunks));

    return mergedChunks;
}

BamFile.prototype.fetch = function(chr, min, max, callback, opts) {
    var thisB = this;
    opts = opts || {};

    var chrId = this.chrToIndex[chr];
    var chunks;
    if (chrId === undefined) {
        chunks = [];
    } else {
        chunks = this.blocksForRange(chrId, min, max);
        if (!chunks) {
            callback(null, 'Error in index fetch');
        }
    }
    
    var records = [];
    var index = 0;
    var data;

    function tramp() {
        if (index >= chunks.length) {
            return callback(records);
        } else if (!data) {
            // dlog('fetching ' + index);
            var c = chunks[index];
            var fetchMin = c.minv.block;
            var fetchMax = c.maxv.block + (1<<16); // *sigh*
            thisB.data.slice(fetchMin, fetchMax - fetchMin).fetch(function(r) {
                data = unbgzf(r, c.maxv.block - c.minv.block + 1);
                return tramp();
            });
        } else {
            var ba = new Uint8Array(data);
            thisB.readBamRecords(ba, chunks[index].minv.offset, records, min, max, chrId, opts);
            data = null;
            ++index;
            return tramp();
        }
    }
    tramp();
}

var SEQRET_DECODER = ['=', 'A', 'C', 'x', 'G', 'x', 'x', 'x', 'T', 'x', 'x', 'x', 'x', 'x', 'x', 'N'];
var CIGAR_DECODER = ['M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X', '?', '?', '?', '?', '?', '?', '?'];

function BamRecord() {
}

BamFile.prototype.readBamRecords = function(ba, offset, sink, min, max, chrId, opts) {
    while (true) {
        var blockSize = readInt(ba, offset);
        var blockEnd = offset + blockSize + 4;
        if (blockEnd >= ba.length) {
            return sink;
        }

        var record = new BamRecord();

        var refID = readInt(ba, offset + 4);
        var pos = readInt(ba, offset + 8);
        
        var bmn = readInt(ba, offset + 12);
        var bin = (bmn & 0xffff0000) >> 16;
        var mq = (bmn & 0xff00) >> 8;
        var nl = bmn & 0xff;

        var flag_nc = readInt(ba, offset + 16);
        var flag = (flag_nc & 0xffff0000) >> 16;
        var nc = flag_nc & 0xffff;
    
        var lseq = readInt(ba, offset + 20);
        
        var nextRef  = readInt(ba, offset + 24);
        var nextPos = readInt(ba, offset + 28);
        
        var tlen = readInt(ba, offset + 32);
    
        record.segment = this.indexToChr[refID];
        record.flag = flag;
        record.pos = pos;
        record.mq = mq;
        if (opts.light)
            record.seqLength = lseq;

        if (!opts.light) {
            var readName = '';
            for (var j = 0; j < nl-1; ++j) {
                readName += String.fromCharCode(ba[offset + 36 + j]);
            }
            record.readName = readName;
        
            var p = offset + 36 + nl;

            var cigar = '';
            for (var c = 0; c < nc; ++c) {
                var cigop = readInt(ba, p);
                cigar = cigar + (cigop>>4) + CIGAR_DECODER[cigop & 0xf];
                p += 4;
            }
            record.cigar = cigar;
        
            var seq = '';
            var seqBytes = (lseq + 1) >> 1;
            for (var j = 0; j < seqBytes; ++j) {
                var sb = ba[p + j];
                seq += SEQRET_DECODER[(sb & 0xf0) >> 4];
                seq += SEQRET_DECODER[(sb & 0x0f)];
            }
            p += seqBytes;
            record.seq = seq;

            var qseq = '';
            for (var j = 0; j < lseq; ++j) {
                qseq += String.fromCharCode(ba[p + j] + 33);
            }
            p += lseq;
            record.quals = qseq;

            while (p < blockEnd) {
                var tag = String.fromCharCode(ba[p], ba[p + 1]);
                var type = String.fromCharCode(ba[p + 2]);
                var value;

                if (type == 'A') {
                    value = String.fromCharCode(ba[p + 3]);
                    p += 4;
                } else if (type == 'i' || type == 'I') {
                    value = readInt(ba, p + 3);
                    p += 7;
                } else if (type == 'c' || type == 'C') {
                    value = ba[p + 3];
                    p += 4;
                } else if (type == 's' || type == 'S') {
                    value = readShort(ba, p + 3);
                    p += 5;
                } else if (type == 'f') {
                    value = readFloat(ba, p + 3);
                    p += 7;
                } else if (type == 'Z' || type == 'H') {
                    p += 3;
                    value = '';
                    for (;;) {
                        var cc = ba[p++];
                        if (cc == 0) {
                            break;
                        } else {
                            value += String.fromCharCode(cc);
                        }
                    }
                } else if (type == 'B') {
                    var atype = String.fromCharCode(ba[p + 3]);
                    var alen = readInt(ba, p + 4);
                    var elen;
                    var reader;
                    if (atype == 'i' || atype == 'I' || atype == 'f') {
                        elen = 4;
                        if (atype == 'f')
                            reader = readFloat;
                        else
                            reader = readInt;
                    } else if (atype == 's' || atype == 'S') {
                        elen = 2;
                        reader = readShort;
                    } else if (atype == 'c' || atype == 'C') {
                        elen = 1;
                        reader = readByte;
                    } else {
                        throw 'Unknown array type ' + atype;
                    }

                    p += 8;
                    value = [];
                    for (var i = 0; i < alen; ++i) {
                        value.push(reader(ba, p));
                        p += elen;
                    }
                } else {
                    throw 'Unknown type '+ type;
                }
                record[tag] = value;
            }
        }

        if (!min || record.pos <= max && record.pos + lseq >= min) {
            if (chrId === undefined || refID == chrId) {
                sink.push(record);
            }
        }
        offset = blockEnd;
    }

    // Exits via top of loop.
};

(function(global) {
    var convertBuffer = new ArrayBuffer(8);
    var ba = new Uint8Array(convertBuffer);
    var fa = new Float32Array(convertBuffer);


    global.readFloat = function(buf, offset) {
        ba[0] = buf[offset];
        ba[1] = buf[offset+1];
        ba[2] = buf[offset+2];
        ba[3] = buf[offset+3];
        return fa[0];
    };
 }(this));

function readInt64(ba, offset) {
    return (ba[offset + 7] << 24) | (ba[offset + 6] << 16) | (ba[offset + 5] << 8) | (ba[offset + 4]);
}

function readInt(ba, offset) {
    return (ba[offset + 3] << 24) | (ba[offset + 2] << 16) | (ba[offset + 1] << 8) | (ba[offset]);
}

function readIntBE(ba, offset) {
    return (ba[offset] << 24) | (ba[offset + 1] << 16) | (ba[offset + 2] << 8) | (ba[offset + 3]);
}

function readShort(ba, offset) {
    return (ba[offset + 1] << 8) | (ba[offset]);
}

function readByte(ba, offset) {
    return ba[offset];
}

function readVob(ba, offset) {
    var block = ((ba[offset+6] & 0xff) * 0x100000000) + ((ba[offset+5] & 0xff) * 0x1000000) + ((ba[offset+4] & 0xff) * 0x10000) + ((ba[offset+3] & 0xff) * 0x100) + ((ba[offset+2] & 0xff));
    var bint = (ba[offset+1] << 8) | (ba[offset]);
    if (block == 0 && bint == 0) {
        return null;  // Should only happen in the linear index?
    } else {
        return new Vob(block, bint);
    }
}

function unbgzf(data, lim) {
    lim = Math.min(lim || 1, data.byteLength - 50);
    var oBlockList = [];
    var ptr = [0];
    var totalSize = 0;

    while (ptr[0] < lim) {
        var ba = new Uint8Array(data, ptr[0], 12); // FIXME is this enough for all credible BGZF block headers?
        var xlen = (ba[11] << 8) | (ba[10]);
        // dlog('xlen[' + (ptr[0]) +']=' + xlen);
        var unc = jszlib_inflate_buffer(data, 12 + xlen + ptr[0], Math.min(65536, data.byteLength - 12 - xlen - ptr[0]), ptr);
        ptr[0] += 8;
        totalSize += unc.byteLength;
        oBlockList.push(unc);
    }

    if (oBlockList.length == 1) {
        return oBlockList[0];
    } else {
        var out = new Uint8Array(totalSize);
        var cursor = 0;
        for (var i = 0; i < oBlockList.length; ++i) {
            var b = new Uint8Array(oBlockList[i]);
            arrayCopy(b, 0, out, cursor, b.length);
            cursor += b.length;
        }
        return out.buffer;
    }
}

//
// Binning (transliterated from SAM1.3 spec)
//

/* calculate bin given an alignment covering [beg,end) (zero-based, half-close-half-open) */
function reg2bin(beg, end)
{
    --end;
    if (beg>>14 == end>>14) return ((1<<15)-1)/7 + (beg>>14);
    if (beg>>17 == end>>17) return ((1<<12)-1)/7 + (beg>>17);
    if (beg>>20 == end>>20) return ((1<<9)-1)/7 + (beg>>20);
    if (beg>>23 == end>>23) return ((1<<6)-1)/7 + (beg>>23);
    if (beg>>26 == end>>26) return ((1<<3)-1)/7 + (beg>>26);
    return 0;
}

/* calculate the list of bins that may overlap with region [beg,end) (zero-based) */
var MAX_BIN = (((1<<18)-1)/7);
function reg2bins(beg, end) 
{
    var i = 0, k, list = [];
    --end;
    list.push(0);
    for (k = 1 + (beg>>26); k <= 1 + (end>>26); ++k) list.push(k);
    for (k = 9 + (beg>>23); k <= 9 + (end>>23); ++k) list.push(k);
    for (k = 73 + (beg>>20); k <= 73 + (end>>20); ++k) list.push(k);
    for (k = 585 + (beg>>17); k <= 585 + (end>>17); ++k) list.push(k);
    for (k = 4681 + (beg>>14); k <= 4681 + (end>>14); ++k) list.push(k);
    return list;
}/* -*- mode: javascript; c-basic-offset: 4; indent-tabs-mode: nil -*- */

// 
// Dalliance Genome Explorer
// (c) Thomas Down 2006-2010
//
// bigwig.js: indexed binary WIG (and BED) files
//

"use strict";

var BIG_WIG_MAGIC = 0x888FFC26;
var BIG_WIG_MAGIC_BE = 0x26FC8F88;
var BIG_BED_MAGIC = 0x8789F2EB;
var BIG_BED_MAGIC_BE = 0xEBF28987;

var BIG_WIG_TYPE_GRAPH = 1;
var BIG_WIG_TYPE_VSTEP = 2;
var BIG_WIG_TYPE_FSTEP = 3;
  
var M1 = 256;
var M2 = 256*256;
var M3 = 256*256*256;
var M4 = 256*256*256*256;

var BED_COLOR_REGEXP = new RegExp("^[0-9]+,[0-9]+,[0-9]+");

function bwg_readOffset(ba, o) {
    var offset = ba[o] + ba[o+1]*M1 + ba[o+2]*M2 + ba[o+3]*M3 + ba[o+4]*M4;
    return offset;
}

function BigWig() {
}

BigWig.prototype.readChromTree = function(callback) {
    var thisB = this;
    this.chromsToIDs = {};
    this.idsToChroms = {};
    this.maxID = 0;

    var udo = this.unzoomedDataOffset;
    var eb = (udo - this.chromTreeOffset) & 3;
    udo = udo + 4 - eb;

    this.data.slice(this.chromTreeOffset, udo - this.chromTreeOffset).fetch(function(bpt) {
        var ba = new Uint8Array(bpt);
        var sa = new Int16Array(bpt);
        var la = new Int32Array(bpt);
        var bptMagic = la[0];
        var blockSize = la[1];
        var keySize = la[2];
        var valSize = la[3];
        var itemCount = bwg_readOffset(ba, 16);
        var rootNodeOffset = 32;

        var bptReadNode = function(offset) {
            var nodeType = ba[offset];
            var cnt = sa[(offset/2) + 1];
            offset += 4;
            for (var n = 0; n < cnt; ++n) {
                if (nodeType == 0) {
                    offset += keySize;
                    var childOffset = bwg_readOffset(ba, offset);
                    offset += 8;
                    childOffset -= thisB.chromTreeOffset;
                    bptReadNode(childOffset);
                } else {
                    var key = '';
                    for (var ki = 0; ki < keySize; ++ki) {
                        var charCode = ba[offset++];
                        if (charCode != 0) {
                            key += String.fromCharCode(charCode);
                        }
                    }
                    var chromId = (ba[offset+3]<<24) | (ba[offset+2]<<16) | (ba[offset+1]<<8) | (ba[offset+0]);
                    var chromSize = (ba[offset + 7]<<24) | (ba[offset+6]<<16) | (ba[offset+5]<<8) | (ba[offset+4]);
                    offset += 8;

                    thisB.chromsToIDs[key] = chromId;
                    if (key.indexOf('chr') == 0) {
                        thisB.chromsToIDs[key.substr(3)] = chromId;
                    }
                    thisB.idsToChroms[chromId] = key;
                    thisB.maxID = Math.max(thisB.maxID, chromId);
                }
            }
        };
        bptReadNode(rootNodeOffset);

        callback(thisB);
    });
}

function BigWigView(bwg, cirTreeOffset, cirTreeLength, isSummary) {
    this.bwg = bwg;
    this.cirTreeOffset = cirTreeOffset;
    this.cirTreeLength = cirTreeLength;
    this.isSummary = isSummary;
}



BigWigView.prototype.readWigData = function(chrName, min, max, callback) {
    var chr = this.bwg.chromsToIDs[chrName];
    if (chr === undefined) {
        // Not an error because some .bwgs won't have data for all chromosomes.

        // dlog("Couldn't find chr " + chrName);
        // dlog('Chroms=' + miniJSONify(this.bwg.chromsToIDs));
        return callback([]);
    } else {
        this.readWigDataById(chr, min, max, callback);
    }
}

BigWigView.prototype.readWigDataById = function(chr, min, max, callback) {
    var thisB = this;
    if (!this.cirHeader) {
        this.bwg.data.slice(this.cirTreeOffset, 48).fetch(function(result) {
            thisB.cirHeader = result;
            var la = new Int32Array(thisB.cirHeader);
            thisB.cirBlockSize = la[1];
            thisB.readWigDataById(chr, min, max, callback);
        });
        return;
    }

    var blocksToFetch = [];
    var outstanding = 0;

    var beforeBWG = Date.now();

    var filter = function(chromId, fmin, fmax, toks) {
        return ((chr < 0 || chromId == chr) && fmin <= max && fmax >= min);
    }

    var cirFobRecur = function(offset, level) {
        if (thisB.bwg.instrument)
            console.log('level=' + level + '; offset=' + offset + '; time=' + (Date.now()|0));

        outstanding += offset.length;

        if (offset.length == 1 && offset[0] - thisB.cirTreeOffset == 48 && thisB.cachedCirRoot) {
            cirFobRecur2(thisB.cachedCirRoot, 0, level);
            --outstanding;
            if (outstanding == 0) {
                thisB.fetchFeatures(filter, blocksToFetch, callback);
            }
            return;
        }

        var maxCirBlockSpan = 4 +  (thisB.cirBlockSize * 32);   // Upper bound on size, based on a completely full leaf node.
        var spans;
        for (var i = 0; i < offset.length; ++i) {
            var blockSpan = new Range(offset[i], offset[i] + maxCirBlockSpan);
            spans = spans ? union(spans, blockSpan) : blockSpan;
        }
        
        var fetchRanges = spans.ranges();
        for (var r = 0; r < fetchRanges.length; ++r) {
            var fr = fetchRanges[r];
            cirFobStartFetch(offset, fr, level);
        }
    }

    var cirFobStartFetch = function(offset, fr, level, attempts) {
        var length = fr.max() - fr.min();
        thisB.bwg.data.slice(fr.min(), fr.max() - fr.min()).fetch(function(resultBuffer) {
            for (var i = 0; i < offset.length; ++i) {
                if (fr.contains(offset[i])) {
                    cirFobRecur2(resultBuffer, offset[i] - fr.min(), level);

                    if (offset[i] - thisB.cirTreeOffset == 48 && offset[i] - fr.min() == 0)
                        thisB.cachedCirRoot = resultBuffer;

                    --outstanding;
                    if (outstanding == 0) {
                        thisB.fetchFeatures(filter, blocksToFetch, callback);
                    }
                }
            }
        });
    }

    var cirFobRecur2 = function(cirBlockData, offset, level) {
        var ba = new Uint8Array(cirBlockData);
        var sa = new Int16Array(cirBlockData);
        var la = new Int32Array(cirBlockData);

        var isLeaf = ba[offset];
        var cnt = sa[offset/2 + 1];
        offset += 4;

        if (isLeaf != 0) {
            for (var i = 0; i < cnt; ++i) {
                var lo = offset/4;
                var startChrom = la[lo];
                var startBase = la[lo + 1];
                var endChrom = la[lo + 2];
                var endBase = la[lo + 3];
                var blockOffset = bwg_readOffset(ba, offset+16);
                var blockSize = bwg_readOffset(ba, offset+24);
                if (((chr < 0 || startChrom < chr) || (startChrom == chr && startBase <= max)) &&
                    ((chr < 0 || endChrom   > chr) || (endChrom == chr && endBase >= min)))
                {
                    blocksToFetch.push({offset: blockOffset, size: blockSize});
                }
                offset += 32;
            }
        } else {
            var recurOffsets = [];
            for (var i = 0; i < cnt; ++i) {
                var lo = offset/4;
                var startChrom = la[lo];
                var startBase = la[lo + 1];
                var endChrom = la[lo + 2];
                var endBase = la[lo + 3];
                var blockOffset = bwg_readOffset(ba, offset+16);
                if ((chr < 0 || startChrom < chr || (startChrom == chr && startBase <= max)) &&
                    (chr < 0 || endChrom   > chr || (endChrom == chr && endBase >= min)))
                {
                    recurOffsets.push(blockOffset);
                }
                offset += 24;
            }
            if (recurOffsets.length > 0) {
                cirFobRecur(recurOffsets, level + 1);
            }
        }
    };

    cirFobRecur([thisB.cirTreeOffset + 48], 1);
}


BigWigView.prototype.fetchFeatures = function(filter, blocksToFetch, callback) {
    var thisB = this;

    blocksToFetch.sort(function(b0, b1) {
        return (b0.offset|0) - (b1.offset|0);
    });

    if (blocksToFetch.length == 0) {
        callback([]);
    } else {
        var features = [];
        var createFeature = function(chr, fmin, fmax, opts) {
            if (!opts) {
                opts = {};
            }
        
            var f = new DASFeature();
            f._chromId = chr;
            f.segment = thisB.bwg.idsToChroms[chr];
            f.min = fmin;
            f.max = fmax;
            f.type = 'bigwig';
            
            for (var k in opts) {
                f[k] = opts[k];
            }
            
            features.push(f);
        };

        var tramp = function() {
            if (blocksToFetch.length == 0) {
                var afterBWG = Date.now();
                // dlog('BWG fetch took ' + (afterBWG - beforeBWG) + 'ms');
                callback(features);
                return;  // just in case...
            } else {
                var block = blocksToFetch[0];
                if (block.data) {
                    thisB.parseFeatures(block.data, createFeature, filter);
                    blocksToFetch.splice(0, 1);
                    tramp();
                } else {
                    var fetchStart = block.offset;
                    var fetchSize = block.size;
                    var bi = 1;
                    while (bi < blocksToFetch.length && blocksToFetch[bi].offset == (fetchStart + fetchSize)) {
                        fetchSize += blocksToFetch[bi].size;
                        ++bi;
                    }

                    thisB.bwg.data.slice(fetchStart, fetchSize).fetch(function(result) {
                        var offset = 0;
                        var bi = 0;
                        while (offset < fetchSize) {
                            var fb = blocksToFetch[bi];
                        
                            var data;
                            if (thisB.bwg.uncompressBufSize > 0) {
                                data = jszlib_inflate_buffer(result, offset + 2, fb.size - 2);
                            } else {
                                var tmp = new Uint8Array(fb.size);    // FIXME is this really the best we can do?
                                arrayCopy(new Uint8Array(result, offset, fb.size), 0, tmp, 0, fb.size);
                                data = tmp.buffer;
                            }
                            fb.data = data;
                            
                            offset += fb.size;
                            ++bi;
                        }
                        tramp();
                    });
                }
            }
        }
        tramp();
    }
}

BigWigView.prototype.parseFeatures = function(data, createFeature, filter) {
    var ba = new Uint8Array(data);

    if (this.isSummary) {
        var sa = new Int16Array(data);
        var la = new Int32Array(data);
        var fa = new Float32Array(data);

        var itemCount = data.byteLength/32;
        for (var i = 0; i < itemCount; ++i) {
            var chromId =   la[(i*8)];
            var start =     la[(i*8)+1];
            var end =       la[(i*8)+2];
            var validCnt =  la[(i*8)+3];
            var minVal    = fa[(i*8)+4];
            var maxVal    = fa[(i*8)+5];
            var sumData   = fa[(i*8)+6];
            var sumSqData = fa[(i*8)+7];
            
            if (filter(chromId, start + 1, end)) {
                var summaryOpts = {type: 'bigwig', score: sumData/validCnt, maxScore: maxVal};
                if (this.bwg.type == 'bigbed') {
                    summaryOpts.type = 'density';
                }
                createFeature(chromId, start + 1, end, summaryOpts);
            }
        }
    } else if (this.bwg.type == 'bigwig') {
        var sa = new Int16Array(data);
        var la = new Int32Array(data);
        var fa = new Float32Array(data);

        var chromId = la[0];
        var blockStart = la[1];
        var blockEnd = la[2];
        var itemStep = la[3];
        var itemSpan = la[4];
        var blockType = ba[20];
        var itemCount = sa[11];
        
        if (blockType == BIG_WIG_TYPE_FSTEP) {
            for (var i = 0; i < itemCount; ++i) {
                var score = fa[i + 6];
                var fmin = blockStart + (i*itemStep) + 1, fmax = blockStart + (i*itemStep) + itemSpan;
                if (filter(chromId, fmin, fmax))
                    createFeature(chromId, fmin, fmax, {score: score});
            }
        } else if (blockType == BIG_WIG_TYPE_VSTEP) {
            for (var i = 0; i < itemCount; ++i) {
                var start = la[(i*2) + 6] + 1;
                var end = start + itemSpan - 1;
                var score = fa[(i*2) + 7];
                if (filter(chromId, start, end))
                    createFeature(chromId, start, end, {score: score});
            }
        } else if (blockType == BIG_WIG_TYPE_GRAPH) {
            for (var i = 0; i < itemCount; ++i) {
                var start = la[(i*3) + 6] + 1;
                var end   = la[(i*3) + 7];
                var score = fa[(i*3) + 8];
                if (start > end) {
                    start = end;
                }
                if (filter(chromId, start, end))
                    createFeature(chromId, start, end, {score: score});
            }
        } else {
            console.log('Currently not handling bwgType=' + blockType);
        }
    } else if (this.bwg.type == 'bigbed') {
        var offset = 0;
        var dfc = this.bwg.definedFieldCount;
        var schema = this.bwg.schema;

        while (offset < ba.length) {
            var chromId = (ba[offset+3]<<24) | (ba[offset+2]<<16) | (ba[offset+1]<<8) | (ba[offset+0]);
            var start = (ba[offset+7]<<24) | (ba[offset+6]<<16) | (ba[offset+5]<<8) | (ba[offset+4]);
            var end = (ba[offset+11]<<24) | (ba[offset+10]<<16) | (ba[offset+9]<<8) | (ba[offset+8]);
            offset += 12;
            var rest = '';
            while (true) {
                var ch = ba[offset++];
                if (ch != 0) {
                    rest += String.fromCharCode(ch);
                } else {
                    break;
                }
            }

            var featureOpts = {};
            
            var bedColumns;
            if (rest.length > 0) {
                bedColumns = rest.split('\t');
            } else {
                bedColumns = [];
            }
            if (bedColumns.length > 0 && dfc > 3) {
                featureOpts.label = bedColumns[0];
            }
            if (bedColumns.length > 1 && dfc > 4) {
                featureOpts.score = stringToInt(bedColumns[1]);
            }
            if (bedColumns.length > 2 && dfc > 5) {
                featureOpts.orientation = bedColumns[2];
            }
            if (bedColumns.length > 5 && dfc > 8) {
                var color = bedColumns[5];
                if (BED_COLOR_REGEXP.test(color)) {
                    featureOpts.itemRgb = 'rgb(' + color + ')';
                }
            }

            if (bedColumns.length > dfc-3 && schema) {
                for (var col = dfc - 3; col < bedColumns.length; ++col) {
                    featureOpts[schema.fields[col+3].name] = bedColumns[col];
                }
            }

            if (filter(chromId, start + 1, end, bedColumns)) {
                if (dfc < 12) {
                    createFeature(chromId, start + 1, end, featureOpts);
                } else {
                    var thickStart = bedColumns[3]|0;
                    var thickEnd   = bedColumns[4]|0;
                    var blockCount = bedColumns[6]|0;
                    var blockSizes = bedColumns[7].split(',');
                    var blockStarts = bedColumns[8].split(',');
                    
                    featureOpts.type = 'bb-transcript'
                    var grp = new DASGroup();
                    for (var k in featureOpts) {
                        grp[k] = featureOpts[k];
                    }
                    grp.id = bedColumns[0];
                    grp.segment = this.bwg.idsToChroms[chromId];
                    grp.min = start + 1;
                    grp.max = end;
                    grp.notes = [];
                    featureOpts.groups = [grp];

                    if (bedColumns.length > 9) {
                        var geneId = bedColumns[9];
                        var geneName = geneId;
                        if (bedColumns.length > 10) {
                            geneName = bedColumns[10];
                        }
                        var gg = shallowCopy(grp);
                        gg.id = geneId;
                        gg.label = geneName;
                        gg.type = 'gene';
                        featureOpts.groups.push(gg);
                    }

                    var spanList = [];
                    for (var b = 0; b < blockCount; ++b) {
                        var bmin = (blockStarts[b]|0) + start;
                        var bmax = bmin + (blockSizes[b]|0);
                        var span = new Range(bmin, bmax);
                        spanList.push(span);
                    }
                    var spans = union(spanList);
                    
                    var tsList = spans.ranges();
                    for (var s = 0; s < tsList.length; ++s) {
                        var ts = tsList[s];
                        createFeature(chromId, ts.min() + 1, ts.max(), featureOpts);
                    }

                    if (thickEnd > thickStart) {
                        var tl = intersection(spans, new Range(thickStart, thickEnd));
                        if (tl) {
                            featureOpts.type = 'bb-translation';
                            var tlList = tl.ranges();
                            for (var s = 0; s < tlList.length; ++s) {
                                var ts = tlList[s];
                                createFeature(chromId, ts.min() + 1, ts.max(), featureOpts);
                            }
                        }
                    }
                }
            }
        }
    } else {
        dlog("Don't know what to do with " + this.bwg.type);
    }
}

//
// nasty cut/paste, should roll back in!
//

BigWigView.prototype.getFirstAdjacent = function(chrName, pos, dir, callback) {
    var chr = this.bwg.chromsToIDs[chrName];
    if (chr === undefined) {
        // Not an error because some .bwgs won't have data for all chromosomes.

        // dlog("Couldn't find chr " + chrName);
        // dlog('Chroms=' + miniJSONify(this.bwg.chromsToIDs));
        return callback([]);
    } else {
        this.getFirstAdjacentById(chr, pos, dir, callback);
    }
}

BigWigView.prototype.getFirstAdjacentById = function(chr, pos, dir, callback) {
    var thisB = this;
    if (!this.cirHeader) {
        this.bwg.data.slice(this.cirTreeOffset, 48).fetch(function(result) {
            thisB.cirHeader = result;
            var la = new Int32Array(thisB.cirHeader);
            thisB.cirBlockSize = la[1];
            thisB.getFirstAdjacentById(chr, pos, dir, callback);
        });
        return;
    }

    var blockToFetch = null;
    var bestBlockChr = -1;
    var bestBlockOffset = -1;

    var outstanding = 0;

    var beforeBWG = Date.now();

    var cirFobRecur = function(offset, level) {
        outstanding += offset.length;

        var maxCirBlockSpan = 4 +  (thisB.cirBlockSize * 32);   // Upper bound on size, based on a completely full leaf node.
        var spans;
        for (var i = 0; i < offset.length; ++i) {
            var blockSpan = new Range(offset[i], offset[i] + maxCirBlockSpan);
            spans = spans ? union(spans, blockSpan) : blockSpan;
        }
        
        var fetchRanges = spans.ranges();
        for (var r = 0; r < fetchRanges.length; ++r) {
            var fr = fetchRanges[r];
            cirFobStartFetch(offset, fr, level);
        }
    }

    var cirFobStartFetch = function(offset, fr, level, attempts) {
        var length = fr.max() - fr.min();
        thisB.bwg.data.slice(fr.min(), fr.max() - fr.min()).fetch(function(resultBuffer) {
            for (var i = 0; i < offset.length; ++i) {
                if (fr.contains(offset[i])) {
                    cirFobRecur2(resultBuffer, offset[i] - fr.min(), level);
                    --outstanding;
                    if (outstanding == 0) {
                        if (!blockToFetch) {
                            if (dir > 0 && (chr != 0 || pos > 0)) {
                                return thisB.getFirstAdjacentById(0, 0, dir, callback);
                            } else if (dir < 0 && (chr != thisB.bwg.maxID || pos < 1000000000)) {
                                return thisB.getFirstAdjacentById(thisB.bwg.maxID, 1000000000, dir, callback);
                            }
                            return callback([]);
                        }

                        thisB.fetchFeatures(function(chrx, fmin, fmax, toks) {
                            return (dir < 0 && (chrx < chr || fmax < pos)) || (dir > 0 && (chrx > chr || fmin > pos));
                        }, [blockToFetch], function(features) {
                            var bestFeature = null;
                            var bestChr = -1;
                            var bestPos = -1;
                            for (var fi = 0; fi < features.length; ++fi) {
                                var f = features[fi];
                                var chrx = f._chromId, fmin = f.min, fmax = f.max;
                                if (bestFeature == null || ((dir < 0) && (chrx > bestChr || fmax > bestPos)) || ((dir > 0) && (chrx < bestChr || fmin < bestPos))) {
                                    bestFeature = f;
                                    bestPos = (dir < 0) ? fmax : fmin;
                                    bestChr = chrx;
                                }
                            }

                            if (bestFeature != null) 
                                return callback([bestFeature]);
                            else
                                return callback([]);
                        });
                    }
                }
            }
        });
    }

    var cirFobRecur2 = function(cirBlockData, offset, level) {
        var ba = new Uint8Array(cirBlockData);
        var sa = new Int16Array(cirBlockData);
        var la = new Int32Array(cirBlockData);

        var isLeaf = ba[offset];
        var cnt = sa[offset/2 + 1];
        offset += 4;

        if (isLeaf != 0) {
            for (var i = 0; i < cnt; ++i) {
                var lo = offset/4;
                var startChrom = la[lo];
                var startBase = la[lo + 1];
                var endChrom = la[lo + 2];
                var endBase = la[lo + 3];
                var blockOffset = bwg_readOffset(ba, offset+16);
                var blockSize = bwg_readOffset(ba, offset+24);
                if ((dir < 0 && ((startChrom < chr || (startChrom == chr && startBase <= pos)))) ||
                    (dir > 0 && ((endChrom > chr || (endChrom == chr && endBase >= pos)))))
                {
                    // console.log('Got an interesting block: startBase=' + startChrom + ':' + startBase + '; endBase=' + endChrom + ':' + endBase + '; offset=' + blockOffset + '; size=' + blockSize);
                    if (/_random/.exec(thisB.bwg.idsToChroms[startChrom])) {
                        // dlog('skipping random: ' + thisB.bwg.idsToChroms[startChrom]);
                    } else if (blockToFetch == null || ((dir < 0) && (endChrom > bestBlockChr || (endChrom == bestBlockChr && endBase > bestBlockOffset)) ||
                                                 (dir > 0) && (startChrom < bestBlockChr || (startChrom == bestBlockChr && startBase < bestBlockOffset))))
                    {
                        //                        dlog('best is: startBase=' + startChrom + ':' + startBase + '; endBase=' + endChrom + ':' + endBase + '; offset=' + blockOffset + '; size=' + blockSize);
                        blockToFetch = {offset: blockOffset, size: blockSize};
                        bestBlockOffset = (dir < 0) ? endBase : startBase;
                        bestBlockChr = (dir < 0) ? endChrom : startChrom;
                    }
                }
                offset += 32;
            }
        } else {
            var bestRecur = -1;
            var bestPos = -1;
            var bestChr = -1;
            for (var i = 0; i < cnt; ++i) {
                var lo = offset/4;
                var startChrom = la[lo];
                var startBase = la[lo + 1];
                var endChrom = la[lo + 2];
                var endBase = la[lo + 3];
                var blockOffset = (la[lo + 4]<<32) | (la[lo + 5]);
                if ((dir < 0 && ((startChrom < chr || (startChrom == chr && startBase <= pos)) &&
                                 (endChrom   >= chr))) ||
                     (dir > 0 && ((endChrom > chr || (endChrom == chr && endBase >= pos)) &&
                                  (startChrom <= chr))))
                {
                    if (bestRecur < 0 || endBase > bestPos) {
                        bestRecur = blockOffset;
                        bestPos = (dir < 0) ? endBase : startBase;
                        bestChr = (dir < 0) ? endChrom : startChrom;
                    }
                }
                offset += 24;
            }
            if (bestRecur >= 0) {
                cirFobRecur([bestRecur], level + 1);
            }
        }
    };
    

    cirFobRecur([thisB.cirTreeOffset + 48], 1);
}

BigWig.prototype.readWigData = function(chrName, min, max, callback) {
    this.getUnzoomedView().readWigData(chrName, min, max, callback);
}

BigWig.prototype.getUnzoomedView = function() {
    if (!this.unzoomedView) {
        var cirLen = 4000;
        var nzl = this.zoomLevels[0];
        if (nzl) {
            cirLen = this.zoomLevels[0].dataOffset - this.unzoomedIndexOffset;
        }
        this.unzoomedView = new BigWigView(this, this.unzoomedIndexOffset, cirLen, false);
    }
    return this.unzoomedView;
}

BigWig.prototype.getZoomedView = function(z) {
    var zh = this.zoomLevels[z];
    if (!zh.view) {
        zh.view = new BigWigView(this, zh.indexOffset, /* this.zoomLevels[z + 1].dataOffset - zh.indexOffset */ 4000, true);
    }
    return zh.view;
}


function makeBwgFromURL(url, callback, creds) {
    makeBwg(new URLFetchable(url, {credentials: creds}), callback, url);
}

function makeBwgFromFile(file, callback) {
    makeBwg(new BlobFetchable(file), callback, 'file');
}

function makeBwg(data, callback, name) {
    var bwg = new BigWig();
    bwg.data = data;
    bwg.name = name;
    bwg.data.slice(0, 512).salted().fetch(function(result) {
        if (!result) {
            return callback(null, "Couldn't fetch file");
        }

        var header = result;
        var ba = new Uint8Array(header);
        var sa = new Int16Array(header);
        var la = new Int32Array(header);
        var magic = ba[0] + (M1 * ba[1]) + (M2 * ba[2]) + (M3 * ba[3]);
        if (magic == BIG_WIG_MAGIC) {
            bwg.type = 'bigwig';
        } else if (magic == BIG_BED_MAGIC) {
            bwg.type = 'bigbed';
        } else if (magic == BIG_WIG_MAGIC_BE || magic == BIG_BED_MAGIC_BE) {
            callback(null, "Currently don't support big-endian BBI files");
        } else {
            callback(null, "Not a supported format, magic=0x" + magic.toString(16));
        }
        // console.log('magic okay');

        bwg.version = sa[2];             // 4
        bwg.numZoomLevels = sa[3];       // 6
        bwg.chromTreeOffset = bwg_readOffset(ba, 8);
        bwg.unzoomedDataOffset = bwg_readOffset(ba, 16);
        bwg.unzoomedIndexOffset = bwg_readOffset(ba, 24);
        bwg.fieldCount = sa[16];         // 32
        bwg.definedFieldCount = sa[17];  // 34
        bwg.asOffset = bwg_readOffset(ba, 36);
        bwg.totalSummaryOffset = bwg_readOffset(ba, 44);
        bwg.uncompressBufSize = la[13];  // 52
        bwg.extHeaderOffset = bwg_readOffset(ba, 56);
        
        // console.log('bwgVersion: ' + bwg.version);
        // dlog('bigType: ' + bwg.type);
        // dlog('chromTree at: ' + bwg.chromTreeOffset);
        // dlog('uncompress: ' + bwg.uncompressBufSize);
        // dlog('data at: ' + bwg.unzoomedDataOffset);
        // dlog('index at: ' + bwg.unzoomedIndexOffset);
        // dlog('field count: ' + bwg.fieldCount);
        // dlog('defined count: ' + bwg.definedFieldCount);

        bwg.zoomLevels = [];
        for (var zl = 0; zl < bwg.numZoomLevels; ++zl) {
            var zlReduction = la[zl*6 + 16]
            var zlData = bwg_readOffset(ba, zl*24 + 72);
            var zlIndex = bwg_readOffset(ba, zl*24 + 80);
            bwg.zoomLevels.push({reduction: zlReduction, dataOffset: zlData, indexOffset: zlIndex});
        }

        bwg.readChromTree(function() {
            bwg.getAutoSQL(function(as) {
                bwg.schema = as;
                return callback(bwg);
            });
        });
    });
}


BigWig.prototype._tsFetch = function(zoom, chr, min, max, callback) {
    var bwg = this;
    if (zoom >= this.zoomLevels.length - 1) {
        if (!this.topLevelReductionCache) {
            this.getZoomedView(this.zoomLevels.length - 1).readWigDataById(-1, 0, 300000000, function(feats) {
                bwg.topLevelReductionCache = feats;
                return bwg._tsFetch(zoom, chr, min, max, callback);
            });
        } else {
            var f = [];
            var c = this.topLevelReductionCache;
            for (var fi = 0; fi < c.length; ++fi) {
                if (c[fi]._chromId == chr) {
                    f.push(c[fi]);
                }
            }
            return callback(f);
        }
    } else {
        var view;
        if (zoom < 0) {
            view = this.getUnzoomedView();
        } else {
            view = this.getZoomedView(zoom);
        }
        return view.readWigDataById(chr, min, max, callback);
    }
}

BigWig.prototype.thresholdSearch = function(chrName, referencePoint, dir, threshold, callback) {
    dir = (dir<0) ? -1 : 1;
    var bwg = this;
    var initialChr = this.chromsToIDs[chrName];
    var candidates = [{chrOrd: 0, chr: initialChr, zoom: bwg.zoomLevels.length - 4, min: 0, max: 300000000, fromRef: true}]
    for (var i = 1; i <= this.maxID + 1; ++i) {
        var chrId = (initialChr + (dir*i)) % (this.maxID + 1);
        if (chrId < 0) 
            chrId += (this.maxID + 1);
        candidates.push({chrOrd: i, chr: chrId, zoom: bwg.zoomLevels.length - 1, min: 0, max: 300000000})
    }
       
    function fbThresholdSearchRecur() {
    	if (candidates.length == 0) {
    	    return callback(null);
    	}
    	candidates.sort(function(c1, c2) {
    	    var d = c1.zoom - c2.zoom;
    	    if (d != 0)
    		    return d;

            d = c1.chrOrd - c2.chrOrd;
            if (d != 0)
                return d;
    	    else
    		    return c1.min - c2.min * dir;
    	});

	    var candidate = candidates.splice(0, 1)[0];
        bwg._tsFetch(candidate.zoom, candidate.chr, candidate.min, candidate.max, function(feats) {
            var rp = dir > 0 ? 0 : 300000000;
            if (candidate.fromRef)
                rp = referencePoint;
            
            for (var fi = 0; fi < feats.length; ++fi) {
    	        var f = feats[fi];
                var score;
                if (f.maxScore != undefined)
                    score = f.maxScore;
                else
                    score = f.score;

                if (dir > 0) {
    	            if (score > threshold) {
        		        if (candidate.zoom < 0) {
        		            if (f.min > rp)
                                return callback(f);
        		        } else if (f.max > rp) {
        		            candidates.push({chr: candidate.chr, chrOrd: candidate.chrOrd, zoom: candidate.zoom - 2, min: f.min, max: f.max, fromRef: candidate.fromRef});
        		        }
                    }
                } else {
                    if (score > threshold) {
            		    if (candidate.zoom < 0) {
                	        if (f.max < rp)
                			    return callback(f);
                        } else if (f.min < rp) {
                            candidates.push({chr: candidate.chr, chrOrd: candidate.chrOrd, zoom: candidate.zoom - 2, min: f.min, max: f.max, fromRef: candidate.fromRef});
                        }
    	            }
                }
    	    }
            fbThresholdSearchRecur();
        });
    }
    
    fbThresholdSearchRecur();
}

BigWig.prototype.getAutoSQL = function(callback) {
    var thisB = this;
    if (!this.asOffset)
        return callback(null);


    this.data.slice(this.asOffset, 2048).fetch(function(result) {
        var ba = new Uint8Array(result);
        var s = '';
        for (var i = 0; i < ba.length; ++i) {
            if (ba[i] == 0)
                break;
            s += String.fromCharCode(ba[i]);
        }
        
        /* 
         * Quick'n'dirty attempt to parse autoSql format.
         * See: http://www.linuxjournal.com/files/linuxjournal.com/linuxjournal/articles/059/5949/5949l2.html
         */

        var header_re = /(\w+)\s+(\w+)\s+("([^"]+)")?\s+\(\s*/;
        var field_re = /([\w\[\]]+)\s+(\w+)\s*;\s*("([^"]+)")?\s*/g;

        var headerMatch = header_re.exec(s);
        if (headerMatch) {
            var as = {
                declType: headerMatch[1],
                name: headerMatch[2],
                comment: headerMatch[4],

                fields: []
            };

            s = s.substring(headerMatch[0]);
            for (var m = field_re.exec(s); m != null; m = field_re.exec(s)) {
                as.fields.push({type: m[1],
                             name: m[2],
                             comment: m[4]});
            }

            return callback(as);
        }
    });
}

BigWig.prototype.getExtraIndices = function(callback) {
    var thisB = this;
    if (this.version < 4 || this.extHeaderOffset == 0 || this.type != 'bigbed') {
        return callback(null);
    } else {
        this.data.slice(this.extHeaderOffset, 64).fetch(function(result) {
            if (!result) {
                return callback(null, "Couldn't fetch extension header");
            }

            var ba = new Uint8Array(result);
            var sa = new Int16Array(result);
            var la = new Int32Array(result);
            
            var extHeaderSize = sa[0];
            var extraIndexCount = sa[1];
            var extraIndexListOffset = bwg_readOffset(ba, 4);

            if (extraIndexCount == 0) {
                return callback(null);
            }

            // FIXME 20byte records only make sense for single-field indices.
            // Right now, these seem to be the only things around, but the format
            // is actually more general.
            thisB.data.slice(extraIndexListOffset, extraIndexCount * 20).fetch(function(eil) {
                if (!eil) {
                    return callback(null, "Couldn't fetch index info");
                }

                var ba = new Uint8Array(eil);
                var sa = new Int16Array(eil);
                var la = new Int32Array(eil);

                var indices = [];
                for (var ii = 0; ii < extraIndexCount; ++ii) {
                    var eiType = sa[ii*10];
                    var eiFieldCount = sa[ii*10 + 1];
                    var eiOffset = bwg_readOffset(ba, ii*20 + 4);
                    var eiField = sa[ii*10 + 8]
                    var index = new BBIExtraIndex(thisB, eiType, eiFieldCount, eiOffset, eiField);
                    indices.push(index);
                }
                callback(indices);
            });
        });
    }
}

function BBIExtraIndex(bbi, type, fieldCount, offset, field) {
    this.bbi = bbi;
    this.type = type;
    this.fieldCount = fieldCount;
    this.offset = offset;
    this.field = field;
}

BBIExtraIndex.prototype.lookup = function(name, callback) {
    var thisB = this;

    this.bbi.data.slice(this.offset, 32).fetch(function(bpt) {
        var ba = new Uint8Array(bpt);
        var sa = new Int16Array(bpt);
        var la = new Int32Array(bpt);
        var bptMagic = la[0];
        var blockSize = la[1];
        var keySize = la[2];
        var valSize = la[3];
        var itemCount = bwg_readOffset(ba, 16);
        var rootNodeOffset = 32;

        function bptReadNode(nodeOffset) {
            thisB.bbi.data.slice(nodeOffset, 4 + (blockSize * (keySize + valSize))).fetch(function(node) {
                var ba = new Uint8Array(node);
                var sa = new Uint16Array(node);
                var la = new Uint32Array(node);

                var nodeType = ba[0];
                var cnt = sa[1];

                var offset = 4;
                if (nodeType == 0) {
                    var lastChildOffset = null;
                    for (var n = 0; n < cnt; ++n) {
                        var key = '';
                        for (var ki = 0; ki < keySize; ++ki) {
                            var charCode = ba[offset++];
                            if (charCode != 0) {
                                key += String.fromCharCode(charCode);
                            }
                        }

                        var childOffset = bwg_readOffset(ba, offset);
                        offset += 8;
                        
                        if (name.localeCompare(key) < 0 && lastChildOffset) {
                            bptReadNode(lastChildOffset);
                            return;
                        }
                        lastChildOffset = childOffset;
                    }
                    bptReadNode(lastChildOffset);
                } else {
                    for (var n = 0; n < cnt; ++n) {
                        var key = '';
                        for (var ki = 0; ki < keySize; ++ki) {
                            var charCode = ba[offset++];
                            if (charCode != 0) {
                                key += String.fromCharCode(charCode);
                            }
                        }
                        
                        // Specific for EI case.
                        if (key == name) {
                            var start = bwg_readOffset(ba, offset);
                            var length = readInt(ba, offset + 8);

                            return thisB.bbi.getUnzoomedView().fetchFeatures(
                                function(chr, min, max, toks) {
                                    if (toks && toks.length > thisB.field - 3)
                                        return toks[thisB.field - 3] == name;
                                }, 
                                [{offset: start, size: length}], 
                                callback);
                        }
                        offset += valSize;
                    }
                    return callback([]);
                }
            });
        }

        bptReadNode(thisB.offset + rootNodeOffset);
    });
}
/* -*- mode: javascript; c-basic-offset: 4; indent-tabs-mode: nil -*- */

// 
// Dalliance Genome Explorer
// (c) Thomas Down 2006-2011
//
// bin.js general binary data support
//

"use strict";

function BlobFetchable(b) {
    this.blob = b;
}

BlobFetchable.prototype.slice = function(start, length) {
    var b;

    if (this.blob.slice) {
        if (length) {
            b = this.blob.slice(start, start + length);
        } else {
            b = this.blob.slice(start);
        }
    } else {
        if (length) {
            b = this.blob.webkitSlice(start, start + length);
        } else {
            b = this.blob.webkitSlice(start);
        }
    }
    return new BlobFetchable(b);
}

BlobFetchable.prototype.salted = function() {return this;}

if (typeof(FileReader) !== 'undefined') {
    // console.log('defining async BlobFetchable.fetch');

    BlobFetchable.prototype.fetch = function(callback) {
        var reader = new FileReader();
        reader.onloadend = function(ev) {
            callback(bstringToBuffer(reader.result));
        };
        reader.readAsBinaryString(this.blob);
    }

} else {
    // if (console && console.log)
    //    console.log('defining sync BlobFetchable.fetch');

    BlobFetchable.prototype.fetch = function(callback) {
        var reader = new FileReaderSync();
        try {
            var res = reader.readAsArrayBuffer(this.blob);
            callback(res);
        } catch (e) {
            callback(null, e);
        }
    }
}

function URLFetchable(url, start, end, opts) {
    if (!opts) {
        if (typeof start === 'object') {
            opts = start;
            start = undefined;
        } else {
            opts = {};
        }
    }

    this.url = url;
    this.start = start || 0;
    if (end) {
        this.end = end;
    }
    this.opts = opts;
}

URLFetchable.prototype.slice = function(s, l) {
    if (s < 0) {
        throw 'Bad slice ' + s;
    }

    var ns = this.start, ne = this.end;
    if (ns && s) {
        ns = ns + s;
    } else {
        ns = s || ns;
    }
    if (l && ns) {
        ne = ns + l - 1;
    } else {
        ne = ne || l - 1;
    }
    return new URLFetchable(this.url, ns, ne, this.opts);
}

var seed=0;
var isSafari = navigator.userAgent.indexOf('Safari') >= 0 && navigator.userAgent.indexOf('Chrome') < 0 ;

URLFetchable.prototype.fetchAsText = function(callback) {
    var req = new XMLHttpRequest();
    var length;
    var url = this.url;
    if (isSafari || this.opts.salt) {
        url = saltURL(url);
        url = url + '?salt=' + b64_sha1('' + Date.now() + ',' + (++seed));
    }
    req.open('GET', url, true);

    if (this.end) {
        if (this.end - this.start > 100000000) {
            throw 'Monster fetch!';
        }
        req.setRequestHeader('Range', 'bytes=' + this.start + '-' + this.end);
        length = this.end - this.start + 1;
    }

    req.onreadystatechange = function() {
        if (req.readyState == 4) {
            if (req.status == 200 || req.status == 206) {
                return callback(req.responseText);
            } else {
                return callback(null);
            }
        }
    };
    if (this.opts.credentials) {
        req.withCredentials = true;
    }
    req.send('');
}

URLFetchable.prototype.salted = function() {
    var o = shallowCopy(this.opts);
    o.salt = true;
    return new URLFetchable(this.url, this.start, this.end, o);
}

URLFetchable.prototype.fetch = function(callback, attempt, truncatedLength) {
    var thisB = this;

    attempt = attempt || 1;
    if (attempt > 3) {
        return callback(null);
    }

    var req = new XMLHttpRequest();
    var length;
    var url = this.url;
    if (isSafari || this.opts.salt) {
        url = url + '?salt=' + b64_sha1('' + Date.now() + ',' + (++seed));
    }
    req.open('GET', url, true);
    req.overrideMimeType('text/plain; charset=x-user-defined');
    if (this.end) {
        if (this.end - this.start > 100000000) {
            throw 'Monster fetch!';
        }
        req.setRequestHeader('Range', 'bytes=' + this.start + '-' + this.end);
        length = this.end - this.start + 1;
    }
    req.responseType = 'arraybuffer';
    req.onreadystatechange = function() {
        if (req.readyState == 4) {
            if (req.status == 200 || req.status == 206) {
                if (req.response) {
                    var bl = req.response.byteLength;
                    if (length && length != bl && (!truncatedLength || bl != truncatedLength)) {
                        return thisB.fetch(callback, attempt + 1, bl);
                    } else {
                        return callback(req.response);
                    }
                } else if (req.mozResponseArrayBuffer) {
                    return callback(req.mozResponseArrayBuffer);
                } else {
                    var r = req.responseText;
                    if (length && length != r.length && (!truncatedLength || r.length != truncatedLength)) {
                        return thisB.fetch(callback, attempt + 1, r.length);
                    } else {
                        return callback(bstringToBuffer(req.responseText));
                    }
                }
            } else {
                return thisB.fetch(callback, attempt + 1);
            }
        }
    };
    if (this.opts.credentials) {
        req.withCredentials = true;
    }
    req.send('');
}

function bstringToBuffer(result) {
    if (!result) {
        return null;
    }

    var ba = new Uint8Array(result.length);
    for (var i = 0; i < ba.length; ++i) {
        ba[i] = result.charCodeAt(i);
    }
    return ba.buffer;
}

/* -*- mode: javascript; c-basic-offset: 4; indent-tabs-mode: nil -*- */

// 
// Dalliance Genome Explorer
// (c) Thomas Down 2006-2011
//
// cbrowser.js: canvas browser container
//

"use strict";

var NS_SVG = 'http://www.w3.org/2000/svg';
var NS_HTML = 'http://www.w3.org/1999/xhtml';
var NS_XLINK = 'http://www.w3.org/1999/xlink';

function Region(chr, min, max) {
    this.min = min;
    this.max = max;
    this.chr = chr;
}

function Browser(opts) {
    if (!opts) {
        opts = {};
    }

    this.uiPrefix = '//www.biodalliance.org/release-0.12/';

    this.sources = [];
    this.tiers = [];

    this.featureListeners = [];
    this.featureHoverListeners = [];
    this.viewListeners = [];
    this.regionSelectListeners = [];
    this.tierListeners = [];
    this.tierSelectionListeners = [];
    this.tierSelectionWrapListeners = [];

    this.cookieKey = 'browser';
    this.registry = 'http://www.dasregistry.org/das/sources';
    this.coordSystem = {
        speciesName: 'Human',
        taxon: 9606,
        auth: 'NCBI',
        version: '36',
        ucscName: 'hg18'
    };
    this.chains = {};

    this.pageName = 'svgHolder'
    this.maxExtra = 2.5;
    this.minExtra = 0.5;
    this.zoomFactor = 1.0;
    this.zoomMin = 10.0;
    this.zoomMax;       // Allow configuration for compatibility, but otherwise clobber.
    this.origin = 0;
    this.targetQuantRes = 1.0;
    this.featurePanelWidth = 750;
    this.zoomBase = 100;
    this.zoomExpt = 30.0; // Back to being fixed....
    this.zoomSliderValue = 100;
    this.entryPoints = null;
    this.currentSeqMax = -1; // init once EPs are fetched.

    this.highlights = [];
    this.selectedTiers = [1];

    this.maxViewWidth = 500000;

    // Options.
    
    this.reverseScrolling = false;
    this.rulerLocation = 'center';
    this.defaultHighlightFill = 'red';
    this.defaultHighlightAlpha = 0.3;
    this.exportHighlights = true;
    this.exportRuler = true;

    // Visual config.

    // this.tierBackgroundColors = ["rgb(245,245,245)", "rgb(230,230,250)" /* 'white' */];
    this.tierBackgroundColors = ["rgb(245,245,245)", 'white'];
    this.minTierHeight = 20;
    this.noDefaultLabels = false;

    // Registry

    this.availableSources = new Observed();
    this.defaultSources = [];
    this.mappableSources = {};

    this.hubs = [];
    this.hubObjects = [];

    this.sourceCache = new SourceCache();
    
    this.retina = true;

    this.useFetchWorkers = true;
    this.maxWorkers = 0;

    this.assemblyNamePrimary = true;
    this.assemblyNameUcsc = true;

    this.initListeners = [];


    for (var k in opts) {
        this[k] = opts[k];
    }

    this.workerPath = '$$js/fetchworker.js';

    var thisB = this;
    window.addEventListener('load', function(ev) {thisB.realInit();}, false);
}

Browser.prototype.resolveURL = function(url) {
    return url.replace('$$', this.uiPrefix);
}

Browser.prototype.realInit = function() {
    if (this.wasInitialized) {
        console.log('Attemping to call realInit on an already-initialized Dalliance instance');
        return;
    }

    this.wasInitialized = true;

    this.defaultChr = this.chr;
    this.defaultStart = this.viewStart;
    this.defaultEnd = this.viewEnd;
    this.defaultSources = [];
    for (var i = 0; i < this.sources.length; ++i) {
        this.defaultSources.push(this.sources[i]);
    }

    if (this.restoreStatus) {
        this.statusRestored = this.restoreStatus();
    }

    var helpPopup;
    var thisB = this;
    this.browserHolderHolder = document.getElementById(this.pageName);
    this.browserHolderHolder.classList.add('dalliance-injection-point');
    this.browserHolder = makeElement('div', null, {className: 'dalliance dalliance-root', tabIndex: -1});
    removeChildren(this.browserHolderHolder);
    this.browserHolderHolder.appendChild(this.browserHolder);
    this.svgHolder = makeElement('div', null, {className: 'main-holder'});

    this.initUI(this.browserHolder, this.svgHolder);

    this.pinnedTierHolder = makeElement('div', null, {className: 'tier-holder tier-holder-pinned'});
    this.tierHolder = makeElement('div', this.makeLoader(24), {className: 'tier-holder tier-holder-rest'});

    this.tierHolderHolder = makeElement('div', [this.pinnedTierHolder, this.tierHolder], {className: 'tier-holder-holder'});
    this.svgHolder.appendChild(this.tierHolderHolder);

    this.bhtmlRoot = makeElement('div');
    if (!this.disablePoweredBy) {
        this.bhtmlRoot.appendChild(makeElement('span', ['Powered by ', makeElement('a', 'Dalliance', {href: 'http://www.biodalliance.org/'}), ' ' + VERSION], {className: 'powered-by'}));
    }
    this.browserHolder.appendChild(this.bhtmlRoot);
    
    window.addEventListener('resize', function(ev) {
        thisB.resizeViewer();
    }, false);

    this.ruler = makeElement('div', null, {className: 'guideline'})
    this.ruler2 = makeElement('div', null, {className: 'guideline'}, {backgroundColor: 'gray', opacity: '0.5', zIndex: 899});
    this.tierHolderHolder.appendChild(this.ruler);
    this.tierHolderHolder.appendChild(this.ruler2);

    if (this.maxWorkers > 0) {
        try {
            var fws = [];
            for (var fi = 0; fi < this.maxWorkers; ++fi)
                fws.push(new FetchWorker());
            this.fetchWorkers = fws;
            this.nextWorker = 0;
        } catch (ex) {
            console.log(ex);
        }
    }

    if (window.getComputedStyle(this.browserHolderHolder).display != 'none') {
        setTimeout(function() {thisB.realInit2()}, 1);
    } else {
        var pollInterval = setInterval(function() {
            if (window.getComputedStyle(thisB.browserHolderHolder).display != 'none') {
                clearInterval(pollInterval);
                thisB.realInit2();
            } 
        }, 300);
    }
}

Browser.prototype.realInit2 = function() {
    var thisB = this;

    // Remove the loader icon, if needed
    removeChildren(this.tierHolder);
    removeChildren(this.pinnedTierHolder);

    this.featurePanelWidth = this.tierHolder.getBoundingClientRect().width | 0;
    this.scale = this.featurePanelWidth / (this.viewEnd - this.viewStart);
    if (!this.zoomMax) {
        this.zoomMax = this.zoomExpt * Math.log(this.maxViewWidth / this.zoomBase);
        this.zoomMin = this.zoomExpt * Math.log(this.featurePanelWidth / 10 / this.zoomBase);
    }
    this.zoomSliderValue = this.zoomExpt * Math.log((this.viewEnd - this.viewStart + 1) / this.zoomBase);

    // Event handlers

    this.tierHolderHolder.addEventListener('mousewheel', function(ev) {
        ev.stopPropagation(); ev.preventDefault();

        if (ev.wheelDeltaX) {
            var delta = ev.wheelDeltaX/5;
            if (!thisB.reverseScrolling) {
                delta = -delta;
            }
            thisB.move(delta);
        }

        if (ev.wheelDeltaY) {
            var delta = ev.wheelDeltaY;
            if (thisB.reverseScrolling) {
                delta = -delta;
            }
            thisB.tierHolder.scrollTop += delta;
        }
    }, false); 

    this.tierHolderHolder.addEventListener('MozMousePixelScroll', function(ev) {
        ev.stopPropagation(); ev.preventDefault();
        if (ev.axis == 1) {
            if (ev.detail != 0) {
                var delta = ev.detail/4;
                if (thisB.reverseScrolling) {
                    delta = -delta;
                }
                thisB.move(delta);
            }
        } else {
            var delta = ev.detail;
            if (!thisB.reverseScrolling) {
              delta = -delta;
            }

            thisB.tierHolder.scrollTop += delta;
        }
    }, false); 

    this.tierHolderHolder.addEventListener('touchstart', function(ev) {return thisB.touchStartHandler(ev)}, false);
    this.tierHolderHolder.addEventListener('touchmove', function(ev) {return thisB.touchMoveHandler(ev)}, false);
    this.tierHolderHolder.addEventListener('touchend', function(ev) {return thisB.touchEndHandler(ev)}, false);
    this.tierHolderHolder.addEventListener('touchcancel', function(ev) {return thisB.touchCancelHandler(ev)}, false);

    var keyHandler = function(ev) {
        // console.log('cbkh: ' + ev.keyCode);
        if (ev.keyCode == 13) { // enter
            var layoutsChanged = false;
            for (var ti = 0; ti < thisB.tiers.length; ++ti) {
                var t = thisB.tiers[ti];
                if (t.wantedLayoutHeight && t.wantedLayoutHeight != t.layoutHeight) {
                    t.layoutHeight = t.wantedLayoutHeight;
                    t.clipTier();
                    layoutsChanged = true;
                }
            }
            if (layoutsChanged) {
                thisB.arrangeTiers();
            }
        } else if (ev.keyCode == 32 || ev.charCode == 32) { // space
            if (!thisB.isSnapZooming) {
                thisB.isSnapZooming = true;
                var newZoom = (thisB.savedZoom || 0.0) + thisB.zoomMin;
                thisB.savedZoom = thisB.zoomSliderValue - thisB.zoomMin;
                thisB.zoomSliderValue = newZoom;
                thisB.zoom(Math.exp((1.0 * newZoom) / thisB.zoomExpt));
            } else {
                thisB.isSnapZooming = false;
                var newZoom = (thisB.savedZoom || 20.0) + thisB.zoomMin;
                thisB.savedZoom = thisB.zoomSliderValue - thisB.zoomMin;
                thisB.zoomSliderValue = newZoom;
                thisB.zoom(Math.exp((1.0 * newZoom) / thisB.zoomExpt));
            }
            thisB.snapZoomLockout = true;
            ev.stopPropagation(); ev.preventDefault();      
        } else if (ev.keyCode == 39) { // right arrow
            ev.stopPropagation(); ev.preventDefault();
            thisB.scrollArrowKey(ev, -1);
        } else if (ev.keyCode == 37) { // left arrow
            ev.stopPropagation(); ev.preventDefault();
            thisB.scrollArrowKey(ev, 1);
        } else if (ev.keyCode == 38 || ev.keyCode == 87) { // up arrow | w
            ev.stopPropagation(); ev.preventDefault();

            if (ev.shiftKey) {
                var st = thisB.getSelectedTier();
                if (st < 0) return;
                var tt = thisB.tiers[st];
                var ch = tt.forceHeight || tt.subtiers[0].height;
                if (ch >= 40) {
                    tt.mergeConfig({height: ch-10});
                }
            } else if (ev.ctrlKey || ev.metaKey) {
                var st = thisB.getSelectedTier();
                if (st < 0) return;
                var tt = thisB.tiers[st];
  
                if (tt.quantLeapThreshold) {
                    var th = tt.subtiers[0].height;
                    var tq = tt.subtiers[0].quant;
                    if (!tq)
                        return;

                    var qmin = 1.0 * tq.min;
                    var qmax = 1.0 * tq.max;

                    var qscale = (qmax - qmin) / th;
                    tt.mergeConfig({quantLeapThreshold: qmin + ((Math.round((tt.quantLeapThreshold - qmin)/qscale)|0)+1)*qscale});

                    tt.notify('Threshold: ' + formatQuantLabel(tt.quantLeapThreshold));
                }                
            } else if (ev.altKey) {
                var cnt = thisB.selectedTiers.length;
                if (cnt == 0)
                    return;

                var st = thisB.selectedTiers[0];
                var contiguous = true;
                var mt = [];
                for (var si = 0; si < thisB.selectedTiers.length; ++si) {
                    mt.push(thisB.tiers[thisB.selectedTiers[si]]);
                    if (si > 0 && thisB.selectedTiers[si] - thisB.selectedTiers[si - 1] != 1)
                        contiguous = false;
                }

                if (contiguous && st <= 0)
                    return;

                for (var si = thisB.selectedTiers.length - 1; si >= 0; --si)
                    thisB.tiers.splice(thisB.selectedTiers[si], 1);

                thisB.selectedTiers.splice(0, cnt);

                var ip = contiguous ? st - 1 : st;
                for (var si = 0; si < mt.length; ++si) {
                    thisB.tiers.splice(ip+si, 0, mt[si]);
                    thisB.selectedTiers.push(ip + si);
                }

                thisB.markSelectedTiers();
                thisB.notifyTierSelection();
                thisB.reorderTiers();
                thisB.notifyTier();
            } else {
                var st = thisB.getSelectedTier();
                if (st > 0) {
                    thisB.setSelectedTier(st - 1);
                    var nst = thisB.tiers[thisB.getSelectedTier()];
                    var top = nst.row.offsetTop, bottom = top + nst.row.offsetHeight;
                    if (top < thisB.tierHolder.scrollTop || bottom > thisB.tierHolder.scrollTop + thisB.tierHolder.offsetHeight) {
                        thisB.tierHolder.scrollTop = top;
                    }
                } else {
                    thisB.notifyTierSelectionWrap(-1);
                }
            }
        } else if (ev.keyCode == 40 || ev.keyCode == 83) { // down arrow | s
            ev.stopPropagation(); ev.preventDefault();

            if (ev.shiftKey) {
                var st = thisB.getSelectedTier();
                if (st < 0) return;
                var tt = thisB.tiers[st];
                var ch = tt.forceHeight || tt.subtiers[0].height;
                tt.mergeConfig({height: ch+10});
            } else if (ev.ctrlKey || ev.metaKey) {
                var st = thisB.getSelectedTier();
                if (st < 0) return;
                var tt = thisB.tiers[st];

                if (tt.quantLeapThreshold) {
                    var th = tt.subtiers[0].height;
                    var tq = tt.subtiers[0].quant;
                    if (!tq)
                        return;

                    var qmin = 1.0 * tq.min;
                    var qmax = 1.0 * tq.max;
                    var qscale = (qmax - qmin) / th;

                    var it = Math.round((tt.quantLeapThreshold - qmin)/qscale)|0;
                    if (it > 1) {
                        tt.mergeConfig({quantLeapThreshold: qmin + (it-1)*qscale});
                        tt.notify('Threshold: ' + formatQuantLabel(tt.quantLeapThreshold));
                    }
                }
            } else if (ev.altKey) {
                var cnt = thisB.selectedTiers.length;
                if (cnt == 0)
                    return;

                var st = thisB.selectedTiers[0];
                var discontig = 0;
                var mt = [];
                for (var si = 0; si < thisB.selectedTiers.length; ++si) {
                    mt.push(thisB.tiers[thisB.selectedTiers[si]]);
                    if (si > 0)
                        discontig += (thisB.selectedTiers[si] - thisB.selectedTiers[si - 1] - 1);
                }
                contiguous = discontig == 0;

                if (contiguous && st + cnt >= thisB.tiers.length)
                    return;

                for (var si = thisB.selectedTiers.length - 1; si >= 0; --si)
                    thisB.tiers.splice(thisB.selectedTiers[si], 1);

                thisB.selectedTiers.splice(0, cnt);

                var ip = contiguous ? st + 1 : st + discontig;
                for (var si = 0; si < mt.length; ++si) {
                    thisB.tiers.splice(ip+si, 0, mt[si]);
                    thisB.selectedTiers.push(ip + si);
                }

                thisB.markSelectedTiers();
                thisB.notifyTierSelection();
                thisB.reorderTiers();
                thisB.notifyTier();
            } else {
                var st = thisB.getSelectedTier();
                if (st < thisB.tiers.length -1) {
                    thisB.setSelectedTier(st + 1);
                    var nst = thisB.tiers[thisB.getSelectedTier()];
                    var top = nst.row.offsetTop, bottom = top + nst.row.offsetHeight;
                    if (top < thisB.tierHolder.scrollTop || bottom > thisB.tierHolder.scrollTop + thisB.tierHolder.offsetHeight) {
                        thisB.tierHolder.scrollTop = bottom - thisB.tierHolder.offsetHeight;
                    }
                }
            }
        } else if (ev.keyCode == 187 || ev.keyCode == 61) { // +
            ev.stopPropagation(); ev.preventDefault();
            thisB.zoomStep(-10);
        } else if (ev.keyCode == 189 || ev.keyCode == 173) { // -
            ev.stopPropagation(); ev.preventDefault();
            thisB.zoomStep(10);
        } else if (ev.keyCode == 73 || ev.keyCode == 105) { // i
            ev.stopPropagation(); ev.preventDefault();
            var st = thisB.getSelectedTier();
            if (st < 0) return;
            var t = thisB.tiers[st];
            if (!t.infoVisible) {
                t.infoElement.style.display = 'block';
                t.updateHeight();
                t.infoVisible = true;
            } else {
                t.infoElement.style.display = 'none';
                t.updateHeight();
                t.infoVisible = false;
            }
        } else if (ev.keyCode == 84 || ev.keyCode == 116) { // t
            var bumpStatus;
            if( ev.shiftKey ){
                ev.stopPropagation(); ev.preventDefault();
                for (var ti = 0; ti < thisB.tiers.length; ++ti) {
                    var t = thisB.tiers[ti];
                    if (t.dasSource.collapseSuperGroups) {
                        if (bumpStatus === undefined) {
                            bumpStatus = !t.bumped;
                        }
                        t.bumped = bumpStatus;
                        t.layoutWasDone = false;
                        t.draw();
                        t.updateLabel();
                    }
                }
            } else if (!ev.ctrlKey && !ev.metaKey) {
                ev.stopPropagation(); ev.preventDefault();
                var st = thisB.getSelectedTier();
                if (st < 0) return;
                var t = thisB.tiers[st];
                if (t.dasSource.collapseSuperGroups) {
                    if (bumpStatus === undefined) {
                        bumpStatus = !t.bumped;
                    }
                    t.bumped = bumpStatus;
                    t.layoutWasDone = false;
                    t.draw();
                    t.updateLabel();
                }
            }
        } else if (ev.keyCode == 77 || ev.keyCode == 109) { // m
            ev.stopPropagation(); ev.preventDefault();
            if ((ev.ctrlKey || ev.metaKey) && thisB.selectedTiers.length > 1) {
                thisB.mergeSelectedTiers();
            }
        } else if (ev.keyCode == 68 || ev.keyCode == 100) { // d
            ev.stopPropagation(); ev.preventDefault();
            if (ev.ctrlKey || ev.metaKey) {
                var st = thisB.getSelectedTier();
                if (st < 0) return;
                thisB.addTier(thisB.tiers[st].dasSource);
            }
        } else if (ev.keyCode == 80 || ev.keyCode == 112) { // p
            if (ev.ctrlKey || ev.metaKey) {
                // Need to be careful because order of tiers could change
                // once we start updating pinning.
                var tt = [];
                for (var st = 0; st < thisB.selectedTiers.length; ++st) {
                    tt.push(thisB.tiers[thisB.selectedTiers[st]]);
                }
                for (var ti = 0; ti < tt.length; ++ti) {
                    tt[ti].mergeConfig({pinned: !tt[ti].pinned});
                }
            }
        } else {
            // console.log('key: ' + ev.keyCode + '; char: ' + ev.charCode);
        }
    };
    var keyUpHandler = function(ev) {
        thisB.snapZoomLockout = false;
    }

    this.browserHolder.addEventListener('focus', function(ev) {
        thisB.browserHolder.addEventListener('keydown', keyHandler, false);
    }, false);
    this.browserHolder.addEventListener('blur', function(ev) {
        thisB.browserHolder.removeEventListener('keydown', keyHandler, false);
    }, false);

    // Popup support (does this really belong here? FIXME)
    this.hPopupHolder = makeElement('div');
    this.hPopupHolder.style['font-family'] = 'helvetica';
    this.hPopupHolder.style['font-size'] = '12pt';
    this.hPopupHolder.classList.add('dalliance');
    document.body.appendChild(this.hPopupHolder);

    for (var t = 0; t < this.sources.length; ++t) {
        var source = this.sources[t];
        var config = {};
        if (this.restoredConfigs) {
            config = this.restoredConfigs[t];
        }

        if (!source.disabled) {
            this.makeTier(source, config);
        }
    }

    thisB.arrangeTiers();
    thisB.refresh();
    thisB.setSelectedTier(1);

    thisB.positionRuler();


    var ss = this.getSequenceSource();
    if (ss) {
        ss.getSeqInfo(this.chr, function(si) {
            thisB.currentSeqMax = si.length;
        });
    }

    this.queryRegistry();
    for (var m in this.chains) {
        this.queryRegistry(m, true);
    }

    if (this.hubs) {
        for (var hi = 0; hi < this.hubs.length; ++hi) {
            var hc = this.hubs[hi];
            if (typeof hc == 'string') {
                hc = {url: hc};
            };

            (function(hc) {
                connectTrackHub(hc.url, function(hub, err) {
                    if (err) {
                        console.log(err);
                    } else {
                        var tdb;
                        if (hc.genome)
                            tdb = hub.genomes[hc.genome];
                        else 
                            tdb = hub.genomes[thisB.coordSystem.ucscName];

                        if (tdb) {
                            if (hc.mapping) 
                                tdb.mapping = hc.mapping;
                            thisB.hubObjects.push(tdb);
                        }
                    }
                }, hc);
            })(hc);
        }
    }

    if (this.fullScreen) {
        this.setFullScreenHeight();
    }

    if (!this.statusRestored && this.storeStatus) {
        this.storeStatus();
    }

    // Ping any init listeners.
    for (var ii = 0; ii < this.initListeners.length; ++ii) {
        try {
            this.initListeners[ii].call(this);
        } catch (e) {
            console.log(e);
        }
    }
}

// 
// iOS touch support

Browser.prototype.touchStartHandler = function(ev) {
    ev.stopPropagation(); ev.preventDefault();
    
    this.touchOriginX = ev.touches[0].pageX;
    this.touchOriginY = ev.touches[0].pageY;
    if (ev.touches.length == 2) {
        var sep = Math.abs(ev.touches[0].pageX - ev.touches[1].pageX);
        this.zooming = true;
        this.zoomLastSep = this.zoomInitialSep = sep;
        this.zoomInitialScale = this.scale;
    }
}

Browser.prototype.touchMoveHandler = function(ev) {
    ev.stopPropagation(); ev.preventDefault();
    
    if (ev.touches.length == 1) {
        var touchX = ev.touches[0].pageX;
        var touchY = ev.touches[0].pageY;
        if (this.touchOriginX && touchX != this.touchOriginX) {
            this.move(touchX - this.touchOriginX);
        }
        if (this.touchOriginY && touchY != this.touchOriginY) {
            this.tierHolder.scrollTop -= (touchY - this.touchOriginY);
        }
        this.touchOriginX = touchX;
        this.touchOriginY = touchY;
    } else if (this.zooming && ev.touches.length == 2) {
        var sep = Math.abs(ev.touches[0].pageX - ev.touches[1].pageX);
        if (sep != this.zoomLastSep) {
            var cp = (ev.touches[0].pageX + ev.touches[1].pageX)/2;
            var scp = this.viewStart + (cp/this.scale)|0
            this.scale = this.zoomInitialScale * (sep/this.zoomInitialSep);
            this.viewStart = scp - (cp/this.scale)|0;
            for (var i = 0; i < this.tiers.length; ++i) {
	           this.tiers[i].draw();
            }
        }
        this.zoomLastSep = sep;
    }
}

Browser.prototype.touchEndHandler = function(ev) {
    ev.stopPropagation(); ev.preventDefault();
}

Browser.prototype.touchCancelHandler = function(ev) {
}


Browser.prototype.makeTier = function(source, config) {
    try {
        return this.realMakeTier(source, config);
    } catch (e) {
        console.log(e.stack || e);
    }
}

Browser.prototype.realMakeTier = function(source, config) {
    var thisB = this;
    var background = null;
    if (this.tierBackgroundColors) {
        background = this.tierBackgroundColors[this.tiers.length % this.tierBackgroundColors.length];
    }

    var tier = new DasTier(this, source, config, background);
    tier.oorigin = this.viewStart

    var isDragging = false;
    var dragOrigin, dragMoveOrigin;
    var hoverTimeout;

    var featureLookup = function(rx, ry) {
        var st = tier.subtiers;
        if (!st) {
            return;
        }

        var sti = 0;
        ry -= MIN_PADDING;
        while (sti < st.length && ry > st[sti].height && sti < (st.length - 1)) {
            ry = ry - st[sti].height - MIN_PADDING;
            ++sti;
        }
        if (sti >= st.length) {
            return;
        }

        var glyphs = st[sti].glyphs;
        var viewCenter = (thisB.viewStart + thisB.viewEnd)/2;
        var offset = (tier.glyphCacheOrigin - thisB.viewStart)*thisB.scale;
        rx -= offset;
       
        return glyphLookup(glyphs, rx, ry);
    }

    var dragMoveHandler = function(ev) {
        ev.preventDefault(); ev.stopPropagation();
        var rx = ev.clientX;
        if (tier.dasSource.tier_type !== 'sequence' && rx != dragMoveOrigin) {
            thisB.move((rx - dragMoveOrigin));
            dragMoveOrigin = rx;
        }
        thisB.isDragging = true;
    }

    var dragUpHandler = function(ev) {
        window.removeEventListener('mousemove', dragMoveHandler, true);
        window.removeEventListener('mouseup', dragUpHandler, true);
    }
        

    tier.viewport.addEventListener('mousedown', function(ev) {
        thisB.browserHolder.focus();
        ev.preventDefault();
        var br = tier.row.getBoundingClientRect();
        var rx = ev.clientX, ry = ev.clientY;

        window.addEventListener('mousemove', dragMoveHandler, true);
        window.addEventListener('mouseup', dragUpHandler, true);
        dragOrigin = dragMoveOrigin = rx;
        thisB.isDragging = false; // Not dragging until a movement event arrives.
    }, false);

    tier.viewport.addEventListener('mousemove', function(ev) {
        var br = tier.row.getBoundingClientRect();
        var rx = ev.clientX - br.left, ry = ev.clientY - br.top;

        if (hoverTimeout) {
            clearTimeout(hoverTimeout);
        }

        if (isDragging) {
            // if (tier.dasSource.tier_type !== 'sequence' && rx != dragMoveOrigin) {
            //    thisB.move((rx - dragMoveOrigin));
            //    dragMoveOrigin = rx;
            // }
        } else {
            hoverTimeout = setTimeout(function() {
                var hit = featureLookup(rx, ry);
                if (hit && hit.length > 0) {
                    thisB.notifyFeatureHover(ev, hit[hit.length - 1], hit, tier);
                }
            }, 1000);
        }
    });

    var doubleClickTimeout = null;
    tier.viewport.addEventListener('mouseup', function(ev) {
        var br = tier.row.getBoundingClientRect();
        var rx = ev.clientX - br.left, ry = ev.clientY - br.top;

        var hit = featureLookup(rx, ry);
        if (hit && hit.length > 0 && !thisB.isDragging) {
            if (doubleClickTimeout) {
                clearTimeout(doubleClickTimeout);
                doubleClickTimeout = null;
                thisB.featureDoubleClick(hit, rx, ry);
            } else {
                doubleClickTimeout = setTimeout(function() {
                    doubleClickTimeout = null;
                    thisB.notifyFeature(ev, hit[hit.length-1], hit, tier);
                }, 500);
            }
        }

        if (thisB.isDragging && rx != dragOrigin && tier.dasSource.tier_type === 'sequence') {
            var a = thisB.viewStart + (rx/thisB.scale);
            var b = thisB.viewStart + (dragOrigin/thisB.scale);

            var min, max;
            if (a < b) {
                min = a|0; max = b|0;
            } else {
                min = b|0; max = a|0;
            }

            thisB.notifyRegionSelect(thisB.chr, min, max);
        }
        thisB.isDragging = false;
    }, false);

    tier.viewport.addEventListener('mouseout', function(ev) {
        isDragging = false;
    });

    tier.removeButton.addEventListener('click', function(ev) {
        ev.stopPropagation(); ev.preventDefault();
        for (var ti = 0; ti < thisB.tiers.length; ++ti) {
            if (thisB.tiers[ti] === tier) {
                thisB.removeTier({index: ti});
                break;
            }
        }
    }, false);
    tier.nameButton.addEventListener('click', function(ev) {
        ev.stopPropagation(); ev.preventDefault();

        if (ev.shiftKey) {
            var hitTier = -1;
            for (var ti = 0; ti < thisB.tiers.length; ++ti) {
                if (thisB.tiers[ti] === tier) {
                    hitTier = ti;
                    break;
                }
            }
            if (hitTier >= 0) {
                var i = thisB.selectedTiers.indexOf(hitTier);
                if (i >= 0) {
                    thisB.selectedTiers.splice(i, 1);
                } else {
                    thisB.selectedTiers.push(hitTier);
                    thisB.selectedTiers.sort();
                }
                thisB.markSelectedTiers();
                thisB.notifyTierSelection();

                if (thisB.selectedTiers.length > 0) {
                    thisB.browserHolder.focus();
                } else {
                    thisB.notifyTierSelectionWrap(-1);
                }
            }
        } else {
            for (var ti = 0; ti < thisB.tiers.length; ++ti) {
                if (thisB.tiers[ti] === tier) {
                    thisB.browserHolder.focus();
                    if (thisB.selectedTiers.length != 1 || thisB.selectedTiers[0] != ti) {
                        thisB.setSelectedTier(ti);
                        return;
                    }
                }
            }

            if (!tier.infoVisible) {
                tier.infoElement.style.display = 'block';
                tier.updateHeight();
                tier.infoVisible = true;
            } else {
                tier.infoElement.style.display = 'none';
                tier.updateHeight();
                tier.infoVisible = false;
            }
        }
    }, false);
    tier.bumpButton.addEventListener('click', function(ev) {
        ev.stopPropagation(); ev.preventDefault();
        var bumpStatus;
        var t = tier;
        if (t.dasSource.collapseSuperGroups) {
            
            if (bumpStatus === undefined) {
                bumpStatus = !t.bumped;
            }
            t.bumped = bumpStatus;
            t.layoutWasDone = false;
            t.draw();
            
            t.updateLabel();
        }
    }, false);

    
    var dragLabel;
    var dragTierHolder;
    var dragTierHolderScrollLimit;
    var tierOrdinal;
    var yAtLastReorder;
    var tiersWereReordered = false;

    var labelDragHandler = function(ev) {
        var label = tier.label;

        ev.stopPropagation(); ev.preventDefault();
        if (!dragLabel) {
            if (tier.pinned) {
                dragTierHolder = thisB.pinnedTierHolder;
            } else {
                dragTierHolder = thisB.tierHolder;
            }
            dragTierHolderScrollLimit = dragTierHolder.scrollHeight - dragTierHolder.offsetHeight;

            dragLabel = label.cloneNode(true);
            dragLabel.style.cursor = 'pointer';
            dragTierHolder.appendChild(dragLabel);
            label.style.visibility = 'hidden';

            for (var ti = 0; ti < thisB.tiers.length; ++ti) {
                if (thisB.tiers[ti] === tier) {
                    tierOrdinal = ti;
                    break;
                }
            }

            yAtLastReorder = ev.clientY;
        }
        
        var holderBCR = dragTierHolder.getBoundingClientRect();
        dragLabel.style.left = (label.getBoundingClientRect().left - holderBCR.left) + 'px'; 
        dragLabel.style.top = (ev.clientY - holderBCR.top + dragTierHolder.scrollTop - 10) + 'px';

        var pty = ev.clientY - holderBCR.top + dragTierHolder.scrollTop;
        for (var ti = 0; ti < thisB.tiers.length; ++ti) {
            var tt = thisB.tiers[ti];
            if (tt.pinned ^ tier.pinned)
                continue; 

            var ttr = tt.row.getBoundingClientRect();
            pty -= (ttr.bottom - ttr.top);
            if (pty < 0) {
                if (ti < tierOrdinal && ev.clientY < yAtLastReorder || ti > tierOrdinal && ev.clientY > yAtLastReorder) {
                    thisB.withPreservedSelection(function() {
                        thisB.tiers.splice(tierOrdinal, 1);
                        thisB.tiers.splice(ti, 0, tier);
                    });

                    tierOrdinal = ti;
                    yAtLastReorder = ev.clientY;
                    thisB.reorderTiers();
                    dragTierHolder.appendChild(dragLabel); // Because reorderTiers removes all children.
                    tiersWereReordered = true;
                }
                break;
            }
        }

        if (dragLabel.offsetTop < dragTierHolder.scrollTop) {
            dragTierHolder.scrollTop -= (dragTierHolder.scrollTop - dragLabel.offsetTop);
        } else if ((dragLabel.offsetTop + dragLabel.offsetHeight) > (dragTierHolder.scrollTop + dragTierHolder.offsetHeight)) {
            dragTierHolder.scrollTop = Math.min(dragTierHolder.scrollTop + 
                                                   (dragLabel.offsetTop + dragLabel.offsetHeight) - 
                                                   (dragTierHolder.scrollTop + dragTierHolder.offsetHeight),
                                                dragTierHolderScrollLimit);
        }
    };

    var labelReleaseHandler = function(ev) {
        var label = tier.label;

        ev.stopPropagation(); ev.preventDefault();
        if (dragLabel) {
            dragLabel.style.cursor = 'auto';
            dragTierHolder.removeChild(dragLabel);
            dragLabel = null;
            label.style.visibility = null;
        }
        document.removeEventListener('mousemove', labelDragHandler, false);
        document.removeEventListener('mouseup', labelReleaseHandler, false);

        if (tiersWereReordered) {
            for (var ti = 0; ti < thisB.tiers.length; ++ti) {
                if (thisB.tiers[ti] == tier) {
                    thisB.setSelectedTier(ti);
                    break;
                }
            }
            thisB.notifyTier();
        }
    };

    tier.label.addEventListener('mousedown', function(ev) {
        ev.stopPropagation(); ev.preventDefault();
        tiersWereReordered = false;
        document.addEventListener('mousemove', labelDragHandler, false);
        document.addEventListener('mouseup', labelReleaseHandler, false);
    }, false);

    this.tiers.push(tier);  // NB this currently tells any extant knownSpace about the new tier.
    
    tier.init(); // fetches stylesheet
    tier.currentlyHeight = 50;
    this.updateHeight();
    tier.updateLabel();

    if (tier.featureSource && tier.featureSource.addActivityListener) {
        tier.featureSource.addActivityListener(function(busy) {
            if (busy > 0) {
                tier.loaderButton.style.display = 'inline-block';
            } else {
                tier.loaderButton.style.display = 'none';
            }
            thisB.pingActivity();
        });
    }

    tier._updateFromConfig();
    this.reorderTiers();

    return tier;
}

Browser.prototype.reorderTiers = function() {
    removeChildren(this.tierHolder);
    removeChildren(this.pinnedTierHolder);
    var hasPinned = false;
    var pinnedTiers = [], unpinnedTiers = [];
    for (var i = 0; i < this.tiers.length; ++i) {
        var t = this.tiers[i];
        if (t.pinned) {
            pinnedTiers.push(t);
            this.pinnedTierHolder.appendChild(this.tiers[i].row);
            hasPinned = true;
        } else {
            unpinnedTiers.push(t);
            this.tierHolder.appendChild(this.tiers[i].row);
        }
    }

    this.withPreservedSelection(function() {
        this.tiers.splice(0, this.tiers.length);
        for (var i = 0; i < pinnedTiers.length; ++i)
            this.tiers.push(pinnedTiers[i]);
        for (var i = 0; i < unpinnedTiers.length; ++i)
            this.tiers.push(unpinnedTiers[i]);
    });

    if (hasPinned)
        this.pinnedTierHolder.classList.add('tier-holder-pinned-full');
    else
        this.pinnedTierHolder.classList.remove('tier-holder-pinned-full');

    this.arrangeTiers();
}

Browser.prototype.withPreservedSelection = function(f) {
    var st = [];
    for (var xi = 0; xi < this.selectedTiers.length; ++xi) {
        st.push(this.tiers[this.selectedTiers[xi]]);
    }

    f.call(this);

    this.selectedTiers = [];
    for (var sti = 0; sti < this.tiers.length; ++sti) {
        if (st.indexOf(this.tiers[sti]) >= 0)
            this.selectedTiers.push(sti);
    }
}

Browser.prototype.refreshTier = function(tier) {
    if (this.knownSpace) {
        this.knownSpace.invalidate(tier);
    }
}

Browser.prototype.arrangeTiers = function() {
    var arrangedTiers = [];
    for (var ti = 0; ti < this.tiers.length; ++ti) {
        var t = this.tiers[ti];
        if (t.pinned) {
            arrangedTiers.push(t);
        }
    }
    for (var ti = 0; ti < this.tiers.length; ++ti) {
        var t = this.tiers[ti];
        if (!t.pinned) {
            arrangedTiers.push(t);
        }
    }

    if (this.tierBackgroundColors) {
        for (var ti = 0; ti < arrangedTiers.length; ++ti) {
            var t = arrangedTiers[ti];
            t.background = this.tierBackgroundColors[ti % this.tierBackgroundColors.length];
        }
    }
}

Browser.prototype.refresh = function() {
    this.notifyLocation();
    var width = (this.viewEnd - this.viewStart) + 1;
    var minExtraW = (100.0/this.scale)|0;
    var maxExtraW = (1000.0/this.scale)|0;

    var newOrigin = (this.viewStart + this.viewEnd) / 2;
    var oh = newOrigin - this.origin;
    this.origin = newOrigin;
    this.scaleAtLastRedraw = this.scale;
    for (var t = 0; t < this.tiers.length; ++t) {
        var od = oh;
        if (this.tiers[t].originHaxx) {
            od += this.tiers[t].originHaxx;
        }
        this.tiers[t].originHaxx = od;
    }

    var scaledQuantRes = this.targetQuantRes / this.scale;

    var innerDrawnStart = Math.max(1, (this.viewStart|0) - minExtraW);
    var innerDrawnEnd = Math.min((this.viewEnd|0) + minExtraW, ((this.currentSeqMax|0) > 0 ? (this.currentSeqMax|0) : 1000000000))
    var outerDrawnStart = Math.max(1, (this.viewStart|0) - maxExtraW);
    var outerDrawnEnd = Math.min((this.viewEnd|0) + maxExtraW, ((this.currentSeqMax|0) > 0 ? (this.currentSeqMax|0) : 1000000000));

    if (!this.knownSpace || this.knownSpace.chr !== this.chr) {
        var ss = this.getSequenceSource();
        this.knownSpace = new KnownSpace(this.tiers, this.chr, outerDrawnStart, outerDrawnEnd, scaledQuantRes, ss);
    }
    
    var seg = this.knownSpace.bestCacheOverlapping(this.chr, innerDrawnStart, innerDrawnEnd);
    if (seg && seg.min <= innerDrawnStart && seg.max >= innerDrawnEnd) {
        this.drawnStart = Math.max(seg.min, outerDrawnStart);
        this.drawnEnd = Math.min(seg.max, outerDrawnEnd);
    } else {
        this.drawnStart = outerDrawnStart;
        this.drawnEnd = outerDrawnEnd;
    }
    
    this.knownSpace.viewFeatures(this.chr, this.drawnStart, this.drawnEnd, scaledQuantRes);
    this.drawOverlays();
}

function setSources(msh, availableSources, maybeMapping) {
    if (maybeMapping) {
        for (var s = 0; s < availableSources.length; ++s) {
            availableSources[s].mapping = maybeMapping;
        }
    }
    msh.set(availableSources);
}

Browser.prototype.queryRegistry = function(maybeMapping, tryCache) {
    var thisB = this;
    var coords, msh;
    if (maybeMapping) {
        coords = this.chains[maybeMapping].coords;
        if (!thisB.mappableSources[maybeMapping]) {
            thisB.mappableSources[maybeMapping] = new Observed();
        }
        msh = thisB.mappableSources[maybeMapping];
    } else {
        coords = this.coordSystem;
        msh = this.availableSources;
    }
    var cacheHash = hex_sha1(miniJSONify(coords));
    if (tryCache) {
        var cacheTime = localStorage['dalliance.registry.' + cacheHash + '.last_queried'];
        if (cacheTime) {
            try {
                setSources(msh, JSON.parse(localStorage['dalliance.registry.' + cacheHash + '.sources']), maybeMapping);
                var cacheAge = (Date.now()|0) - (cacheTime|0);
                if (cacheAge < (12 * 60 * 60 * 1000)) {
                    return;
                }
            } catch (rex) {
                console.log('Bad registry cache: ' + rex);
            }
        }
    }
            
    new DASRegistry(this.registry).sources(function(sources) {
        var availableSources = [];
        for (var s = 0; s < sources.length; ++s) {
            var source = sources[s];
            if (!source.coords || source.coords.length == 0) {
                continue;
            }
            var scoords = source.coords[0];
            if (scoords.taxon != coords.taxon || scoords.auth != coords.auth || scoords.version != coords.version) {
                continue;
            }   
            availableSources.push(source);
        }

        localStorage['dalliance.registry.' + cacheHash + '.sources'] = JSON.stringify(availableSources);
        localStorage['dalliance.registry.' + cacheHash + '.last_queried'] = '' + Date.now();
        
        setSources(msh, availableSources, maybeMapping);
    }, function(error) {
        // msh.set(null);
    }, coords);
}

//
// Navigation
//

Browser.prototype.move = function(pos)
{
    var wid = this.viewEnd - this.viewStart;
    this.viewStart -= pos / this.scale;
    this.viewEnd = this.viewStart + wid;
    if (this.currentSeqMax > 0 && this.viewEnd > this.currentSeqMax) {
        this.viewEnd = this.currentSeqMax;
        this.viewStart = this.viewEnd - wid;
    }
    if (this.viewStart < 1) {
        this.viewStart = 1;
        this.viewEnd = this.viewStart + wid;
    }
    this.notifyLocation();
    
    var viewCenter = (this.viewStart + this.viewEnd)/2;
    
    for (var i = 0; i < this.tiers.length; ++i) {
        var offset = (this.viewStart - this.tiers[i].norigin)*this.scale;
	this.tiers[i].viewport.style.left = '' + ((-offset|0) - 1000) + 'px';
        var ooffset = (this.viewStart - this.tiers[i].oorigin)*this.scale;
        this.tiers[i].overlay.style.left = '' + ((-ooffset|0) - 1000) + 'px';
    }

    this.spaceCheck();
}

Browser.prototype.zoomStep = function(delta) {
    var oz = 1.0 * this.zoomSliderValue;
    var nz = oz + delta;
    if (nz < this.zoomMin) {
        nz= this.zoomMin;
    }
    if (nz > this.zoomMax) {
        nz = this.zoomMax;
    }

    if (nz != oz) {
        this.zoomSliderValue = nz; // FIXME maybe ought to set inside zoom!
        this.zoom(Math.exp((1.0 * nz) / this.zoomExpt));
    }
}

Browser.prototype.zoom = function(factor) {
    this.zoomFactor = factor;
    var viewCenter = Math.round((this.viewStart + this.viewEnd) / 2.0)|0;
    this.viewStart = viewCenter - this.zoomBase * this.zoomFactor / 2;
    this.viewEnd = viewCenter + this.zoomBase * this.zoomFactor / 2;
    if (this.currentSeqMax > 0 && (this.viewEnd > this.currentSeqMax + 5)) {
        var len = this.viewEnd - this.viewStart + 1;
        this.viewEnd = this.currentSeqMax;
        this.viewStart = this.viewEnd - len + 1;
    }
    if (this.viewStart < 1) {
        var len = this.viewEnd - this.viewStart + 1;
        this.viewStart = 1;
        this.viewEnd = this.viewStart + len - 1;
    }
    this.scale = this.featurePanelWidth / (this.viewEnd - this.viewStart)
    var width = this.viewEnd - this.viewStart + 1;
    
    var scaleRat = (this.scale / this.scaleAtLastRedraw);

    this.refresh();
}

Browser.prototype.spaceCheck = function(dontRefresh) {
    if (!this.knownSpace || this.knownSpace.chr !== this.chr) {
        this.refresh();
        return;
    } 

    var width = ((this.viewEnd - this.viewStart)|0) + 1;
    var minExtraW = (100.0/this.scale)|0;
    var maxExtraW = (1000.0/this.scale)|0;

    if ((this.drawnStart|0) > Math.max(1, ((this.viewStart|0) - minExtraW)|0)  || (this.drawnEnd|0) < Math.min((this.viewEnd|0) + minExtraW, ((this.currentSeqMax|0) > 0 ? (this.currentSeqMax|0) : 1000000000)))  {
        this.refresh();
    }
}

Browser.prototype.resizeViewer = function(skipRefresh) {
    var width = this.tierHolder.getBoundingClientRect().width | 0;
    if (width == 0)
        return;

    var oldFPW = Math.max(this.featurePanelWidth, 300); // Can get silly values stored
                                                        // when the browser is hidden.
    this.featurePanelWidth = width|0;

    if (oldFPW != this.featurePanelWidth) {
        this.zoomMax = this.zoomExpt * Math.log(this.maxViewWidth / this.zoomBase);
        this.zoomMin = this.zoomExpt * Math.log(this.featurePanelWidth / 10 / this.zoomBase);
        this.zoomSliderValue = this.zoomExpt * Math.log((this.viewEnd - this.viewStart + 1) / this.zoomBase);

        var viewWidth = this.viewEnd - this.viewStart;
        var nve = this.viewStart + (viewWidth * this.featurePanelWidth) / oldFPW;

        this.viewEnd = nve;

        var wid = this.viewEnd - this.viewStart + 1;
        if (this.currentSeqMax > 0 && this.viewEnd > this.currentSeqMax) {
            this.viewEnd = this.currentSeqMax;
            this.viewStart = this.viewEnd - wid + 1;
        }
        if (this.viewStart < 1) {
            this.viewStart = 1;
            this.viewEnd = this.viewStart + wid - 1;
        }

        this.positionRuler();

        if (!skipRefresh) {
            this.spaceCheck();
        }
        this.notifyLocation();
    }

    if (this.fullScreen) {
        this.setFullScreenHeight();
    }
}

Browser.prototype.setFullScreenHeight = function() {
    var rest = document.body.offsetHeight - this.browserHolder.offsetHeight;
    this.browserHolder.style.maxHeight = Math.max(300, window.innerHeight - rest - 20) + 'px'
}

Browser.prototype.addTier = function(conf) {
    conf = shallowCopy(conf);
    conf.disabled = false;
    
    var tier = this.makeTier(conf);
    this.markSelectedTiers();
    this.positionRuler();
    this.notifyTier();
    return tier;
}

function sourceDataURI(conf) {
    if (conf.uri) {
        return conf.uri;
    } else if (conf.blob) {
        return 'file:' + conf.blob.name;
    } else if (conf.bwgBlob) {
        return 'file:' + conf.bwgBlob.name;
    } else if (conf.bamBlob) {
        return 'file:' + conf.bamBlob.name;
    }

    return conf.bwgURI || conf.bamURI || conf.jbURI || conf.twoBitURI || 'http://www.biodalliance.org/magic/no_uri';
}

function sourceStyleURI(conf) {
    if (conf.stylesheet_uri)
        return conf.stylesheet_uri;
    else if (conf.tier_type == 'sequence')
        return 'http://www.biodalliance.org/magic/sequence'
    else
        return sourceDataURI(conf);
}

function sourcesAreEqual(a, b) {
    if (sourceDataURI(a) != sourceDataURI(b) ||
        sourceStyleURI(a) != sourceStyleURI(b))
        return false;

    if (a.mapping != b.mapping)
        return false;

    if (a.tier_type != b.tier_type)
        return false;

    if (a.overlay) {
        if (!b.overlay || b.overlay.length != a.overlay.length)
            return false;
        for (var oi = 0; oi < a.overlay.length; ++oi) {
            if (!sourcesAreEqual(a.overlay[oi], b.overlay[oi]))
                return false;
        }
    } else {
        if (b.overlay)
            return false;
    }

    return true;
}

Browser.prototype.removeTier = function(conf, force) {
    var target = -1;

    if (typeof conf.index !== 'undefined' && conf.index >=0 && conf.index < this.tiers.length) {
        target = conf.index;
    } else {
        for (var ti = 0; ti < this.tiers.length; ++ti) {
            var ts = this.tiers[ti].dasSource;
            
            if (sourcesAreEqual(conf, ts)) {
                target = ti; break;
            }
        }
    }

    if (target < 0) {
        throw "Couldn't find requested tier";
    }

    this.tiers.splice(target, 1);

    var nst = [];
    for (var sti = 0; sti < this.selectedTiers.length; ++sti) {
        var st = this.selectedTiers[sti];
        if (st < target) {
            nst.push(st);
        } else if (st > target) {
            nst.push(st - 1);
        }
    }
    this.selectedTiers = nst;
    this.markSelectedTiers();

    this.reorderTiers();
    this.notifyTier();
}

Browser.prototype.getSequenceSource = function() {
    if (this._sequenceSource === undefined)
        this._sequenceSource = this._getSequenceSource();
    return this._sequenceSource;
}

Browser.prototype._getSequenceSource = function() {
    for (var ti = 0; ti < this.tiers.length; ++ti) {
        if (this.tiers[ti].sequenceSource) {
            return this.tiers[ti].sequenceSource;
        }
    }

    for (var si = 0; si < this.defaultSources.length; ++si) {
        var s = this.defaultSources[si];
        if (s.provides_entrypoints || s.tier_type == 'sequence' || s.twoBitURI) {
            if (s.twoBitURI) {
                return new TwoBitSequenceSource(s);
            } else {
                return new DASSequenceSource(s);
            }
        }
    }
}

Browser.prototype.setLocation = function(newChr, newMin, newMax, callback) {
    if (!callback) {
        callback = function(err) {
            if (err) {
                throw err;
            }
        }
    }
    var thisB = this;

    if (!newChr || newChr == this.chr) {
        return this._setLocation(null, newMin, newMax, null, callback);
    } else {
        var ss = this.getSequenceSource();
        if (!ss) {
            return callback('Need a sequence source');
        }

        ss.getSeqInfo(newChr, function(si) {
            if (!si) {
                var altChr;
                if (newChr.indexOf('chr') == 0) {
                    altChr = newChr.substr(3);
                } else {
                    altChr = 'chr' + newChr;
                }
                ss.getSeqInfo(altChr, function(si2) {
                    if (!si2) {
                        return callback("Couldn't find sequence '" + newChr + "'");
                    } else {
                        return thisB._setLocation(altChr, newMin, newMax, si2, callback);
                    }
                });
            } else {
                return thisB._setLocation(newChr, newMin, newMax, si, callback);
            }
        });
    }
}

Browser.prototype._setLocation = function(newChr, newMin, newMax, newChrInfo, callback) {
    var chrChanged = false;
    if (newChr) {
        if (newChr.indexOf('chr') == 0)
            newChr = newChr.substring(3);

        if (this.chr != newChr)
            chrChanged = true;
        this.chr = newChr;
        this.currentSeqMax = newChrInfo.length;
    }

    newMin|=0; newMax|=0;
    var newWidth = Math.max(10, newMax-newMin+1);
    if (newMin < 1) {
        newMin = 1; newMax = newMin + newWidth - 1;
    }
    if (newMax > this.currentSeqMax) {
        newMax = this.currentSeqMax;
        newMin = Math.max(1, newMax - newWidth + 1);
    }

    this.viewStart = newMin;
    this.viewEnd = newMax;
    var newScale = Math.max(this.featurePanelWidth, 50) / (this.viewEnd - this.viewStart);
    var oldScale = this.scale;
    var scaleChanged = (Math.abs(newScale - oldScale)) > 0.0001;
    this.scale = newScale;

    var newZS, oldZS;
    oldZS = this.zoomSliderValue;
    this.zoomSliderValue = newZS = this.zoomExpt * Math.log((this.viewEnd - this.viewStart + 1) / this.zoomBase);
    
    if (scaleChanged || chrChanged) {
        for (var i = 0; i < this.tiers.length; ++i) {
            this.tiers[i].viewport.style.left = '5000px';
            this.tiers[i].overlay.style.left = '5000px';
        }

        this.refresh();

        if (this.savedZoom) {
            newZS -= this.zoomMin;
            oldZS -= this.zoomMin;
            var difToActive = newZS - oldZS;
            var difToSaved = newZS - this.savedZoom;
            if (Math.abs(difToActive) > Math.abs(difToSaved)) {
                this.isSnapZooming = !this.isSnapZooming;
                this.savedZoom = oldZS;
            }
        } else {
            this.isSnapZooming = false;
            this.savedZoom = null;
        }
    } else {
        var viewCenter = (this.viewStart + this.viewEnd)/2;
    
        for (var i = 0; i < this.tiers.length; ++i) {
            var offset = (this.viewStart - this.tiers[i].norigin)*this.scale;
	        this.tiers[i].viewport.style.left = '' + ((-offset|0) - 1000) + 'px';
            var ooffset = (this.viewStart - this.tiers[i].oorigin)*this.scale;
            this.tiers[i].overlay.style.left = '' + ((-ooffset|0) - 1000) + 'px';
        }
    }
    this.notifyLocation();

    this.spaceCheck();
    if (this.instrumentActivity)
        this.activityStartTime = Date.now()|0;
    return callback();
}

Browser.prototype.pingActivity = function() {
    if (!this.instrumentActivity || !this.activityStartTime)
        return;

    var activity = 0;
    for (var ti = 0; ti < this.tiers.length; ++ti) {
        if (this.tiers[ti].loaderButton.style.display !== 'none')
            ++activity;
    }

    if (activity == 0) {
        var now = Date.now()|0;
        console.log('Loading took ' + (now-this.activityStartTime) + 'ms');
        this.activityStartTime = null;
    }
}

Browser.prototype.addInitListener = function(handler) {
    this.initListeners.push(handler);
}

Browser.prototype.addFeatureListener = function(handler, opts) {
    opts = opts || {};
    this.featureListeners.push(handler);
}

Browser.prototype.notifyFeature = function(ev, feature, hit, tier) {
  for (var fli = 0; fli < this.featureListeners.length; ++fli) {
      try {
          if (this.featureListeners[fli](ev, feature, hit, tier))
            return;
      } catch (ex) {
          console.log(ex.stack);
      }
  }
}

Browser.prototype.addFeatureHoverListener = function(handler, opts) {
    opts = opts || {};
    this.featureHoverListeners.push(handler);
}

Browser.prototype.notifyFeatureHover = function(ev, feature, hit, tier) {
    for (var fli = 0; fli < this.featureHoverListeners.length; ++fli) {
        try {
            this.featureHoverListeners[fli](ev, feature, hit, tier);
        } catch (ex) {
            console.log(ex.stack);
        }
    }
}

Browser.prototype.addViewListener = function(handler, opts) {
    opts = opts || {};
    this.viewListeners.push(handler);
}

Browser.prototype.notifyLocation = function() {
    for (var lli = 0; lli < this.viewListeners.length; ++lli) {
        try {
            this.viewListeners[lli](this.chr, this.viewStart|0, this.viewEnd|0, this.zoomSliderValue, {current: this.zoomSliderValue, min: this.zoomMin, max: this.zoomMax});
        } catch (ex) {
            console.log(ex.stack);
        }
    }
}

Browser.prototype.addTierListener = function(handler) {
    this.tierListeners.push(handler);
}

Browser.prototype.notifyTier = function() {
    for (var tli = 0; tli < this.tierListeners.length; ++tli) {
        try {
            this.tierListeners[tli]();
        } catch (ex) {
            console.log(ex.stack);
        }
    }
}

Browser.prototype.addRegionSelectListener = function(handler) {
    this.regionSelectListeners.push(handler);
}

Browser.prototype.notifyRegionSelect = function(chr, min, max) {
    for (var rli = 0; rli < this.regionSelectListeners.length; ++rli) {
        try {
            this.regionSelectListeners[rli](chr, min, max);
        } catch (ex) {
            console.log(ex.stack);
        }
    }
}


Browser.prototype.highlightRegion = function(chr, min, max) {
    var thisB = this;
    
    if (chr == this.chr) {
        return this._highlightRegion(chr, min, max);
    }

    var ss = this.getSequenceSource();
    if (!ss) {
        throw 'Need a sequence source';
    }

    ss.getSeqInfo(chr, function(si) {
        if (!si) {
            var altChr;
            if (chr.indexOf('chr') == 0) {
                altChr = chr.substr(3);
            } else {
                altChr = 'chr' + chr;
            }
            ss.getSeqInfo(altChr, function(si2) {
                if (!si2) {
                    // Fail silently.
                } else {
                    return thisB._highlightRegion(altChr, min, max);
                }
            });
        } else {
            return thisB._highlightRegion(chr, min, max);
        }
    });
}

Browser.prototype._highlightRegion = function(chr, min, max) {
    for (var hi = 0; hi < this.highlights.length; ++hi) {
        var h = this.highlights[hi];
        if (h.chr == chr && h.min == min && h.max == max)
            return;
    }

    this.highlights.push(new Region(chr, min, max));
    var visStart = this.viewStart - (1000/this.scale);
    var visEnd = this.viewEnd + (1000/this.scale);
    if ((chr == this.chr || chr == ('chr'+this.chr)) && min < visEnd && max > visStart) {
        this.drawOverlays();
    }

    this.notifyLocation();
}

Browser.prototype.clearHighlights = function() {
    this.highlights = [];
    this.drawOverlays();
    this.notifyLocation();
}

Browser.prototype.drawOverlays = function() {
    for (var ti = 0; ti < this.tiers.length; ++ti) {
        this.tiers[ti].drawOverlay();
    }
}

Browser.prototype.featuresInRegion = function(chr, min, max) {
    var features = [];
    if (chr !== this.chr) {
        return [];
    }

    for (var ti = 0; ti < this.tiers.length; ++ti) {
        var fl = this.tiers[ti].currentFeatures || [];
        for (var fi = 0; fi < fl.length; ++fi) {
            var f = fl[fi];
            if (f.min <= max && f.max >= min) {
                features.push(f);
            }
        }
    }
    return features;
}


Browser.prototype.getSelectedTier = function() {
    if (this.selectedTiers.length > 0) 
        return this.selectedTiers[0];
    else
        return -1;
}

Browser.prototype.setSelectedTier = function(t) {
    if (t == null) {
        this.selectedTiers = [];
    } else {
        this.selectedTiers = [t];
    }
    this.markSelectedTiers();
    this.notifyTierSelection();
}

Browser.prototype.markSelectedTiers = function() {
    for (var ti = 0; ti < this.tiers.length; ++ti) {
        var button = this.tiers[ti].nameButton;

        if (this.selectedTiers.indexOf(ti) >= 0) {
            button.classList.add('active');
        } else {
            button.classList.remove('active');
        }
    }
    if (this.selectedTiers.length > 0) {
        this.browserHolder.focus();
    }
}

Browser.prototype.addTierSelectionListener = function(f) {
    this.tierSelectionListeners.push(f);
}

Browser.prototype.notifyTierSelection = function() {
    for (var fli = 0; fli < this.tierSelectionListeners.length; ++fli) {
        try {
            this.tierSelectionListeners[fli](this.selectedTiers);
        } catch (ex) {
            console.log(ex.stack);
        }
    }
}

Browser.prototype.addTierSelectionWrapListener = function(f) {
    this.tierSelectionWrapListeners.push(f);
}

Browser.prototype.notifyTierSelectionWrap = function(i) {
    for (var fli = 0; fli < this.tierSelectionWrapListeners.length; ++fli) {
        try {
            this.tierSelectionWrapListeners[fli](i);
        } catch (ex) {
            console.log(ex.stack);
        }
    }
}

Browser.prototype.positionRuler = function() {
    var display = 'none';
    var left = '';
    var right = '';

    if (this.rulerLocation == 'center') {
        display = 'block';
        left = '' + ((this.featurePanelWidth/2)|0) + 'px';
    } else if (this.rulerLocation == 'left') {
        display = 'block';
        left = '0px';
    } else if (this.rulerLocation == 'right') {
        display = 'block';
        right = '0px'
    } else {
        display = 'none';
    }

    this.ruler.style.display = display;
    this.ruler.style.left = left;
    this.ruler.style.right = right;

    this.ruler2.style.display = this.rulerLocation == 'center' ? 'none' : 'block';
    this.ruler2.style.left = '' + ((this.featurePanelWidth/2)|0) + 'px';

    for (var ti = 0; ti < this.tiers.length; ++ti) {
        var tier = this.tiers[ti];
        var q = tier.quantOverlay;

        var quant;
        if (tier.subtiers && tier.subtiers.length > 0)
            quant = tier.subtiers[0].quant;

        if (q) {
            q.style.display = quant ? display : 'none';
            q.style.left = left;
            q.style.right = right;
        }
    }
}

Browser.prototype.featureDoubleClick = function(hit, rx, ry) {
    if (!hit || hit.length == 0)
        return;

    var f = hit[hit.length - 1];

    if (!f.min || !f.max) {
        return;
    }

    var fstart = (((f.min|0) - (this.viewStart|0)) * this.scale);
    var fwidth = (((f.max - f.min) + 1) * this.scale);
    
    var newMid = (((f.min|0) + (f.max|0)))/2;
    if (fwidth > 10) {
        var frac = (1.0 * (rx - fstart)) / fwidth;
        if (frac < 0.3) {
            newMid = (f.min|0);
        } else  if (frac > 0.7) {
            newMid = (f.max|0) + 1;
        }
    }

    var width = this.viewEnd - this.viewStart;
    this.setLocation(null, newMid - (width/2), newMid + (width/2));
}

Browser.prototype.zoomForScale = function(scale) {
    var ssScale;
    if (scale > 0.2) {
        ssScale = 'high';
    } else if (scale > 0.01) {
        ssScale = 'medium';
    } else  {
        ssScale = 'low';
    }
    return ssScale;
}

Browser.prototype.zoomForCurrentScale = function() {
    return this.zoomForScale(this.scale);
}

Browser.prototype.updateHeight = function() {
    var tierTotal = 0;
    for (var ti = 0; ti < this.tiers.length; ++ti) 
        tierTotal += (this.tiers[ti].currentHeight || 30);
    this.ruler.style.height = '' + tierTotal + 'px';
    this.ruler2.style.height = '' + tierTotal + 'px';
    // this.svgHolder.style.maxHeight = '' + Math.max(tierTotal, 500) + 'px';
}

Browser.prototype.scrollArrowKey = function(ev, dir) {
    if (this.reverseKeyScrolling)
        dir = -dir;
    
    if (ev.ctrlKey || ev.metaKey) {
        var fedge = false;
        if(ev.shiftKey){
            fedge = true;
        }

        this.leap(dir, fedge);
    } else {
        this.move(ev.shiftKey ? 100*dir : 25*dir);
    }
}

Browser.prototype.leap = function(dir, fedge) {
    var thisB = this;
    var pos=((thisB.viewStart + thisB.viewEnd + 1)/2)|0;

    var st = thisB.getSelectedTier();
    if (st < 0) return;
    var tier = thisB.tiers[st];

    if (tier && ((tier.featureSource && sourceAdapterIsCapable(tier.featureSource, 'quantLeap') && typeof(tier.quantLeapThreshold) == 'number')
                 || (tier.featureSource && sourceAdapterIsCapable(tier.featureSource, 'leap')))) {
        tier.findNextFeature(
              thisB.chr,
              pos,
              -dir,
              fedge,
              function(nxt) {
                  if (nxt) {
                      var nmin = nxt.min;
                      var nmax = nxt.max;
                      if (fedge) { 
                        if (dir > 0) {
                          if (nmin>pos+1) {
                              nmax=nmin;
                          } else {
                              nmax++;
                              nmin=nmax
                          }
                        } else {
                            if (nmax<pos-1) {
                                nmax++;
                                nmin=nmax;
                            } else {
                                nmax=nmin;
                            }
                        } 
                      }
                      var wid = thisB.viewEnd - thisB.viewStart + 1;
                      if(parseFloat(wid/2) == parseInt(wid/2)){wid--;}
                      var newStart = (nmin + nmax - wid)/2 + 1;
                      var newEnd = newStart + wid - 1;
                      var pos2=pos;
                      thisB.setLocation(nxt.segment, newStart, newEnd);
                  } else {
                      alert('no next feature'); // FIXME better reporting would be nice!
                  }
              });
    } else {
        this.move(100*dir);
    }
}

function glyphLookup(glyphs, rx, ry, matches) {
    matches = matches || [];

    for (var gi = glyphs.length - 1; gi >= 0; --gi) {
        var g = glyphs[gi];
        if (!g.notSelectable && g.min() <= rx && g.max() >= rx) {
            if (g.minY) {
                if (ry < g.minY() || ry > g.maxY())
                    continue;
            }

            if (g.feature) {
                matches.push(g.feature);
            } else if (g.group) {
                matches.push(g.group);
            }
    
            if (g.glyphs) {
                return glyphLookup(g.glyphs, rx, ry, matches);
            } else if (g.glyph) {
                return glyphLookup([g.glyph], rx, ry, matches);
            } else {
                return matches;
            }
        }
    }
    return matches;
}

Browser.prototype.nameForCoordSystem = function(cs) {
    var primary = null, ucsc = null;
    if (this.assemblyNamePrimary) {
        primary = '' + cs.auth;
        if (typeof(cs.version) !== 'undefined')
            primary += cs.version;
    }
    if (this.assemblyNameUcsc) {
        ucsc = cs.ucscName;
    }
    if (primary != null && ucsc != null)
        return primary + '/' + ucsc;
    else 
        return primary || ucsc || 'unknown';
}

Browser.prototype.makeLoader = function(size) {
    size = size || 16;
    var retina = window.devicePixelRatio > 1;
    if (size < 20) {
        return makeElement('img', null, {src: this.uiPrefix + 'img/spinner_' + (retina ? 16 : 32) + '.gif', width: '16', height: '16'});
    } else {
        return makeElement('img', null, {src: this.uiPrefix + 'img/spinner_' + (retina ? 24 : 48) + '.gif', width: '24', height: '24'});
    }
}

Browser.prototype.getWorker = function() {
    if (!this.useFetchWorkers || !this.fetchWorkers || this.fetchWorkers.length==0)
        return null;

    if (this.nextWorker >= this.fetchWorkers.length)
        this.nextWorker = 0;
    return this.fetchWorkers[this.nextWorker++];
}

function FetchWorker() {
    var thisB = this;
    this.worker = new Worker('js/fetchworker.js');
    this.tagSeed = 0;
    this.callbacks = {};

    this.worker.onmessage = function(ev) {
        var cb = thisB.callbacks[ev.data.tag];
        if (cb) {
            cb(ev.data.result, ev.data.error);
            delete thisB.callbacks[ev.data.tag];
        }
    }
}

FetchWorker.prototype.postCommand = function(cmd, callback, transfer) {
    var tag = 'x' + (++this.tagSeed);
    cmd.tag = tag;
    this.callbacks[tag] = callback;
    this.worker.postMessage(cmd, transfer);
}
/* -*- mode: javascript; c-basic-offset: 4; indent-tabs-mode: nil -*- */

// 
// Dalliance Genome Explorer
// (c) Thomas Down 2006-2011
//
// feature-popup.js
//

"use strict";

var TAGVAL_NOTE_RE = new RegExp('^([A-Za-z]+)=(.+)');

Browser.prototype.addFeatureInfoPlugin = function(handler) {
    if (!this.featureInfoPlugins) {
        this.featureInfoPlugins = [];
    }
    this.featureInfoPlugins.push(handler);
}

function FeatureInfo(hit, feature, group) {
    var name = pick(group.type, feature.type);
    var fid = pick(group.label, feature.label, group.id, feature.id);
    if (fid && fid.indexOf('__dazzle') != 0) {
        name = name + ': ' + fid;
    }

    this.hit = hit;
    this.feature = feature;
    this.group = group;
    this.title = name;
    this.sections = [];
}

FeatureInfo.prototype.setTitle = function(t) {
    this.title = t;
}

FeatureInfo.prototype.add = function(label, info) {
    if (typeof info === 'string') {
        info = makeElement('span', info);
    }
    this.sections.push({label: label, info: info});
}

Browser.prototype.featurePopup = function(ev, __ignored_feature, hit, tier) {
    var hi = hit.length;
    var feature = --hi >= 0 ? hit[hi] : {};
    var group = --hi >= 0 ? hit[hi] : {};

    var featureInfo = new FeatureInfo(hit, feature, group);
    featureInfo.tier = tier;
    var fips = this.featureInfoPlugins || [];
    for (var fipi = 0; fipi < fips.length; ++fipi) {
        try {
            fips[fipi](feature, featureInfo);
        } catch (e) {
            console.log(e.stack || e);
        }
    }
    fips = tier.featureInfoPlugins || [];
    for (fipi = 0; fipi < fips.length; ++fipi) {
        try {
            fips[fipi](feature, featureInfo);
        } catch (e) {
            console.log(e.stack || e);
        }
    }

    this.removeAllPopups();

    var table = makeElement('table', null, {className: 'table table-striped table-condensed'});
    table.style.width = '100%';
    table.style.margin = '0px';

    var idx = 0;
    if (feature.method) {
        var row = makeElement('tr', [
            makeElement('th', 'Method'),
            makeElement('td', feature.method)
        ]);
        table.appendChild(row);
        ++idx;
    }
    {
        var loc;
        if (group.segment) {
            loc = group;
        } else {
            loc = feature;
        }
        var row = makeElement('tr', [
            makeElement('th', 'Location'),
            makeElement('td', loc.segment + ':' + loc.min + '-' + loc.max, {}, {minWidth: '200px'})
        ]);
        table.appendChild(row);
        ++idx;
    }
    if (feature.score !== undefined && feature.score !== null && feature.score != '-') {
        var row = makeElement('tr', [
            makeElement('th', 'Score'),
            makeElement('td', '' + feature.score)
        ]);
        table.appendChild(row);
        ++idx;
    }
    {
        var links = maybeConcat(group.links, feature.links);
        if (links && links.length > 0) {
            var row = makeElement('tr', [
                makeElement('th', 'Links'),
                makeElement('td', links.map(function(l) {
                    return makeElement('div', makeElement('a', l.desc, {href: l.uri, target: '_new'}));
                }))
            ]);
            table.appendChild(row);
            ++idx;
        }
    }
    {
        var notes = maybeConcat(group.notes, feature.notes);
        for (var ni = 0; ni < notes.length; ++ni) {
            var k = 'Note';
            var v = notes[ni];
            var m = v.match(TAGVAL_NOTE_RE);
            if (m) {
                k = m[1];
                v = m[2];
            }

            var row = makeElement('tr', [
                makeElement('th', k),
                makeElement('td', v)
            ]);
            table.appendChild(row);
            ++idx;
        }
    }

    for (var fisi = 0; fisi < featureInfo.sections.length; ++fisi) {
        var section = featureInfo.sections[fisi];
        table.appendChild(makeElement('tr', [
            makeElement('th', section.label),
            makeElement('td', section.info)]));
    }        

    this.popit(ev, featureInfo.title || 'Feature', table, {width: 450});
}
/* -*- mode: javascript; c-basic-offset: 4; indent-tabs-mode: nil -*- */

// 
// Dalliance Genome Explorer
// (c) Thomas Down 2006-2010
//
// chainset.js: liftover support
//

"use strict";

function Chainset(uri, srcTag, destTag, coords) {
    this.uri = uri;
    this.srcTag = srcTag;
    this.destTag = destTag;
    this.coords = coords;
    this.chainsBySrc = {};
    this.chainsByDest = {};
    this.postFetchQueues = {};
};

(function(scope) {

var CIGAR_REGEXP = new RegExp('([0-9]*)([MIDS])', 'g');

function parseCigar(cigar)
{
    var cigops = [];
    var match;
    while ((match = CIGAR_REGEXP.exec(cigar)) != null) {
        var count = match[1];
        if (count.length == 0) {
            count = 1;
        }
        cigops.push({cnt: count|0, op: match[2]});
    }
    return cigops;
}

scope.parseCigar = parseCigar;

}(this));

Chainset.prototype.fetchChainsTo = function(chr) {
    var thisCS = this;
    new DASSource(this.uri).alignments(chr, {}, function(aligns) {
        if (!thisCS.chainsByDest[chr]) {
            thisCS.chainsByDest[chr] = []; // prevent re-fetching.
        }

        for (var ai = 0; ai < aligns.length; ++ai) {
            var aln = aligns[ai];
            for (var bi = 0; bi < aln.blocks.length; ++bi) {
                var block = aln.blocks[bi];
                var srcSeg, destSeg;
                for (var si = 0; si < block.segments.length; ++si) {
                    var seg = block.segments[si];
                    var obj = aln.objects[seg.object];
                    if (obj.dbSource === thisCS.srcTag) {
                        srcSeg = seg;
                    } else if (obj.dbSource === thisCS.destTag) {
                        destSeg = seg;
                    }
                }
                if (srcSeg && destSeg) {
                    var chain = {
                        srcChr:     aln.objects[srcSeg.object].accession,
                        srcMin:     srcSeg.min|0,
                        srcMax:     srcSeg.max|0,
                        srcOri:     srcSeg.strand,
                        destChr:    aln.objects[destSeg.object].accession,
                        destMin:    destSeg.min|0,
                        destMax:    destSeg.max|0,
                        destOri:    destSeg.strand,
                        blocks:     []
                    }

                    var srcops = parseCigar(srcSeg.cigar), destops = parseCigar(destSeg.cigar);

                    var srcOffset = 0, destOffset = 0;
                    var srci = 0, desti = 0;
                    while (srci < srcops.length && desti < destops.length) {
                        if (srcops[srci].op == 'M' && destops[desti].op == 'M') {
                            var blockLen = Math.min(srcops[srci].cnt, destops[desti].cnt);
                            chain.blocks.push([srcOffset, destOffset, blockLen]);
                            if (srcops[srci].cnt == blockLen) {
                                ++srci;
                            } else {
                                srcops[srci].cnt -= blockLen;
                            }
                            if (destops[desti].cnt == blockLen) {
                                ++desti;
                            } else {
                                destops[desti] -= blockLen;
                            }
                            srcOffset += blockLen;
                            destOffset += blockLen;
                        } else if (srcops[srci].op == 'I') {
                            destOffset += srcops[srci++].cnt;
                        } else if (destops[desti].op == 'I') {
                            srcOffset += destops[desti++].cnt;
                        }
                    }

                    pusho(thisCS.chainsBySrc, chain.srcChr, chain);
                    pusho(thisCS.chainsByDest, chain.destChr, chain);
                }
            }
        }

        if (thisCS.postFetchQueues[chr]) {
            var pfq = thisCS.postFetchQueues[chr];
            for (var i = 0; i < pfq.length; ++i) {
                pfq[i]();
            }
            thisCS.postFetchQueues[chr] = null;
        }
    });
}

Chainset.prototype.mapPoint = function(chr, pos) {
    var chains = this.chainsBySrc[chr] || [];
    for (var ci = 0; ci < chains.length; ++ci) {
        var c = chains[ci];
        if (pos >= c.srcMin && pos <= c.srcMax) {
            var cpos;
            if (c.srcOri == '-') {
                cpos = c.srcMax - pos;
            } else {
                cpos = pos - c.srcMin;
            }
            var blocks = c.blocks;
            for (var bi = 0; bi < blocks.length; ++bi) {
                var b = blocks[bi];
                var bSrc = b[0];
                var bDest = b[1];
                var bSize = b[2];
                if (cpos >= bSrc && cpos <= (bSrc + bSize)) {
                    var apos = cpos - bSrc;

                    var dpos;
                    if (c.destOri == '-') {
                        dpos = c.destMax - bDest - apos;
                    } else {
                        dpos = apos + bDest + c.destMin;
                    }
                    return {seq: c.destChr, pos: dpos, flipped: (c.srcOri != c.destOri)}
                }
            }
        }
    }
    return null;
}

Chainset.prototype.unmapPoint = function(chr, pos) {
    var chains = this.chainsByDest[chr] || [];
    for (var ci = 0; ci < chains.length; ++ci) {
        var c = chains[ci];
        if (pos >= c.destMin && pos <= c.destMax) {
            var cpos;
            if (c.srcOri == '-') {
                cpos = c.destMax - pos;
            } else {
                cpos = pos - c.destMin;
            }    
            
            var blocks = c.blocks;
            for (var bi = 0; bi < blocks.length; ++bi) {
                var b = blocks[bi];
                var bSrc = b[0];
                var bDest = b[1];
                var bSize = b[2];

                if (cpos >= bDest && cpos <= (bDest + bSize)) {
                    var apos = cpos - bDest;

                    var dpos = apos + bSrc + c.srcMin;
                    var dpos;
                    if (c.destOri == '-') {
                        dpos = c.srcMax - bSrc - apos;
                    } else {
                        dpos = apos + bSrc + c.srcMin;
                    }
                    return {seq: c.srcChr, pos: dpos, flipped: (c.srcOri != c.destOri)}
                }
            }
            // return null;
        }
    }
    return null;
}

Chainset.prototype.sourceBlocksForRange = function(chr, min, max, callback) {
    if (!this.chainsByDest[chr]) {
        var fetchNeeded = !this.postFetchQueues[chr];
        var thisCS = this;
        pusho(this.postFetchQueues, chr, function() {
            thisCS.sourceBlocksForRange(chr, min, max, callback);
        });
        if (fetchNeeded) {
            this.fetchChainsTo(chr);
        }
    } else {
        var srcBlocks = [];
        var chains = this.chainsByDest[chr] || [];
        for (var ci = 0; ci < chains.length; ++ci) {
            var c = chains[ci];
            if (min <= c.destMax && max >= c.destMin) {
                var cmin, cmax;
                if (c.srcOri == '-') {
                    cmin = c.destMax - max;
                    cmax = c.destMax - min;
                } else {
                    cmin = min - c.destMin;
                    cmax = max - c.destMin;
                }

                var blocks = c.blocks;
                for (var bi = 0; bi < blocks.length; ++bi) {
                    var b = blocks[bi];
                    var bSrc = b[0];
                    var bDest = b[1];
                    var bSize = b[2];

                    if (cmax >= bDest && cmin <= (bDest + bSize)) {
                        var amin = Math.max(cmin, bDest) - bDest;
                        var amax = Math.min(cmax, bDest + bSize) - bDest;

                        if (c.destOri == '-') {
                            srcBlocks.push(new DASSegment(c.srcChr, c.srcMax - bSrc - amax, c.srcMax - bSrc - amin));
                        } else {
                            srcBlocks.push(new DASSegment(c.srcChr, c.srcMin + amin + bSrc, c.srcMin + amax + bSrc));
                        }
                    }
                }
            }
        }
        callback(srcBlocks);
    }
}
/* -*- mode: javascript; c-basic-offset: 4; indent-tabs-mode: nil -*- */

// 
// Dalliance Genome Explorer
// (c) Thomas Down 2006-2010
//
// color.js
//

"use strict";

function DColour(red, green, blue, name) {
    this.red = red|0;
    this.green = green|0;
    this.blue = blue|0;
    if (name) {
        this.name = name;
    }
}

DColour.prototype.toSvgString = function() {
    if (!this.name) {
        this.name = "rgb(" + this.red + "," + this.green + "," + this.blue + ")";
    }

    return this.name;
}

function hex2(x) {
    var y = '00' + x.toString(16);
    return y.substring(y.length - 2);
}

DColour.prototype.toHexString = function() {
    return '#' + hex2(this.red) + hex2(this.green) + hex2(this.blue);
}

var palette = {
    red: new DColour(255, 0, 0, 'red'),
    green: new DColour(0, 255, 0, 'green'),
    blue: new DColour(0, 0, 255, 'blue'),
    yellow: new DColour(255, 255, 0, 'yellow'),
    white: new DColour(255, 255, 255, 'white'),
    black: new DColour(0, 0, 0, 'black'),
    gray: new DColour(180, 180, 180, 'gray'),
    grey: new DColour(180, 180, 180, 'grey')
};

var COLOR_RE = new RegExp('^#([0-9A-Fa-f]{2})([0-9A-Fa-f]{2})([0-9A-Fa-f]{2})$');
var CSS_COLOR_RE = /rgb\(([0-9]+),([0-9]+),([0-9]+)\)/

function dasColourForName(name) {
    var c = palette[name];
    if (!c) {
        var match = COLOR_RE.exec(name);
        if (match) {
            c = new DColour(('0x' + match[1])|0, ('0x' + match[2])|0, ('0x' + match[3])|0, name);
            palette[name] = c;
        } else {
    	    match = CSS_COLOR_RE.exec(name);
    	    if (match) {
        		c = new DColour(match[1]|0, match[2]|0, match[3]|0, name);
        		palette[name] = c;
	       } else {
		      console.log("couldn't handle color: " + name);
		      c = palette.black;
		      palette[name] = c;
	       }
        }
    }
    return c;
}

function makeColourSteps(steps, stops, colours) {
    var dcolours = [];
    for (var ci = 0; ci < colours.length; ++ci) {
        dcolours.push(dasColourForName(colours[ci]));
    }

    var grad = [];
  STEP_LOOP:
    for (var si = 0; si < steps; ++si) {
        var rs = (1.0 * si) / (steps-1);
        var score = stops[0] + (stops[stops.length -1] - stops[0]) * rs;
        for (var i = 0; i < stops.length - 1; ++i) {
            if (score >= stops[i] && score <= stops[i+1]) {
                var frac = (score - stops[i]) / (stops[i+1] - stops[i]);
                var ca = dcolours[i];
                var cb = dcolours[i+1];

                var fill = new DColour(
                    ((ca.red * (1.0 - frac)) + (cb.red * frac))|0,
                    ((ca.green * (1.0 - frac)) + (cb.green * frac))|0,
                    ((ca.blue * (1.0 - frac)) + (cb.blue * frac))|0
                ).toSvgString();
                grad.push(fill);

                continue STEP_LOOP;
            }
        }
        throw 'Bad step';
    }

    return grad;
}

function makeGradient(steps, color1, color2, color3) {
    if (color3) {
        return makeColourSteps(steps, [0, 0.5, 1], [color1, color2, color3]);
    } else {
        return makeColourSteps(steps, [0, 1], [color1, color2]);
    }
}
/* -*- mode: javascript; c-basic-offset: 4; indent-tabs-mode: nil -*- */

// 
// Dalliance Genome Explorer
// (c) Thomas Down 2006-2010
//
// das.js: queries and low-level data model.
//

"use strict";

var dasLibErrorHandler = function(errMsg) {
    alert(errMsg);
}
var dasLibRequestQueue = new Array();



function DASSegment(name, start, end, description) {
    this.name = name;
    this.start = start;
    this.end = end;
    this.description = description;
}
DASSegment.prototype.toString = function() {
    return this.name + ':' + this.start + '..' + this.end;
};
DASSegment.prototype.isBounded = function() {
    return this.start && this.end;
}
DASSegment.prototype.toDASQuery = function() {
    var q = 'segment=' + this.name;
    if (this.start && this.end) {
        q += (':' + this.start + ',' + this.end);
    }
    return q;
}


function DASSource(a1, a2) {
    var options;
    if (typeof a1 == 'string') {
        this.uri = a1;
        options = a2 || {};
    } else {
        options = a1 || {};
    }
    for (var k in options) {
        if (typeof(options[k]) != 'function') {
            this[k] = options[k];
        }
    }


    if (!this.coords) {
        this.coords = [];
    }
    if (!this.props) {
        this.props = {};
    }

    this.dasBaseURI = this.uri;
    if (this.dasBaseURI && this.dasBaseURI.substr(this.uri.length - 1) != '/') {
        this.dasBaseURI = this.dasBaseURI + '/';
    }
}

function DASCoords() {
}

function coordsMatch(c1, c2) {
    return c1.taxon == c2.taxon && c1.auth == c2.auth && c1.version == c2.version;
}

//
// DAS 1.6 entry_points command
//

DASSource.prototype.entryPoints = function(callback) {
    var dasURI = this.dasBaseURI + 'entry_points';
    this.doCrossDomainRequest(dasURI, function(responseXML) {
            if (!responseXML) {
                return callback([]);
            }

                var entryPoints = new Array();
                
                var segs = responseXML.getElementsByTagName('SEGMENT');
                for (var i = 0; i < segs.length; ++i) {
                    var seg = segs[i];
                    var segId = seg.getAttribute('id');
                    
                    var segSize = seg.getAttribute('size');
                    var segMin, segMax;
                    if (segSize) {
                        segMin = 1; segMax = segSize|0;
                    } else {
                        segMin = seg.getAttribute('start');
                        if (segMin) {
                            segMin |= 0;
                        }
                        segMax = seg.getAttribute('stop');
                        if (segMax) {
                            segMax |= 0;
                        }
                    }
                    var segDesc = null;
                    if (seg.firstChild) {
                        segDesc = seg.firstChild.nodeValue;
                    }
                    entryPoints.push(new DASSegment(segId, segMin, segMax, segDesc));
                }          
               callback(entryPoints);
    });         
}

//
// DAS 1.6 sequence command
// Do we need an option to fall back to the dna command?
//

function DASSequence(name, start, end, alpha, seq) {
    this.name = name;
    this.start = start;
    this.end = end;
    this.alphabet = alpha;
    this.seq = seq;
}

DASSource.prototype.sequence = function(segment, callback) {
    var dasURI = this.dasBaseURI + 'sequence?' + segment.toDASQuery();
    this.doCrossDomainRequest(dasURI, function(responseXML) {
        if (!responseXML) {
            callback([]);
            return;
        } else {
                var seqs = new Array();
                
                var segs = responseXML.getElementsByTagName('SEQUENCE');
                for (var i = 0; i < segs.length; ++i) {
                    var seg = segs[i];
                    var segId = seg.getAttribute('id');
                    var segMin = seg.getAttribute('start');
                    var segMax = seg.getAttribute('stop');
                    var segAlpha = 'DNA';
                    var segSeq = null;
                    if (seg.firstChild) {
                        var rawSeq = seg.firstChild.nodeValue;
                        segSeq = '';
                        var idx = 0;
                        while (true) {
                            var space = rawSeq.indexOf('\n', idx);
                            if (space >= 0) {
                                segSeq += rawSeq.substring(idx, space);
                                idx = space + 1;
                            } else {
                                segSeq += rawSeq.substring(idx);
                                break;
                            }
                        }
                    }
                    seqs.push(new DASSequence(segId, segMin, segMax, segAlpha, segSeq));
                }
                
                callback(seqs);
        }
    });
}

//
// DAS 1.6 features command
//

function DASFeature() {
    // We initialize these in the parser...
}

function DASGroup() {
    // We initialize these in the parser, too...
}

function DASLink(desc, uri) {
    this.desc = desc;
    this.uri = uri;
}

DASSource.prototype.features = function(segment, options, callback) {
    options = options || {};
    var thisB = this;

    var dasURI;
    if (this.features_uri) {
        dasURI = this.features_uri;
    } else {
        var filters = [];

        if (segment) {
            filters.push(segment.toDASQuery());
        } else if (options.group) {
            var g = options.group;
            if (typeof g == 'string') {
                filters.push('group_id=' + g);
            } else {
                for (var gi = 0; gi < g.length; ++gi) {
                    filters.push('group_id=' + g[gi]);
                }
            }
        }

        if (options.adjacent) {
            var adj = options.adjacent;
            if (typeof adj == 'string') {
                adj = [adj];
            }
            for (var ai = 0; ai < adj.length; ++ai) {
                filters.push('adjacent=' + adj[ai]);
            }
        }

        if (options.type) {
            if (typeof options.type == 'string') {
                filters.push('type=' + options.type);
            } else {
                for (var ti = 0; ti < options.type.length; ++ti) {
                    filters.push('type=' + options.type[ti]);
                }
            }
        }
        
        if (options.maxbins) {
            filters.push('maxbins=' + options.maxbins);
        }
        
        if (filters.length > 0) {
            dasURI = this.dasBaseURI + 'features?' + filters.join(';');
        } else {
            callback([], 'No filters specified');
        }
    } 
   

    this.doCrossDomainRequest(dasURI, function(responseXML, req) {
        if (!responseXML) {
            var msg;
            if (req.status == 0) {
                msg = 'server may not support CORS';
            } else {
                msg = 'status=' + req.status;
            }
            callback([], 'Failed request: ' + msg);
            return;
        }
/*      if (req) {
            var caps = req.getResponseHeader('X-DAS-Capabilties');
            if (caps) {
                alert(caps);
            }
        } */

        var features = new Array();
        var segmentMap = {};

        var segs = responseXML.getElementsByTagName('SEGMENT');
        for (var si = 0; si < segs.length; ++si) {
            var segmentXML = segs[si];
            var segmentID = segmentXML.getAttribute('id');
            segmentMap[segmentID] = {
                min: segmentXML.getAttribute('start'),
                max: segmentXML.getAttribute('stop')
            };
            
            var featureXMLs = segmentXML.getElementsByTagName('FEATURE');
            for (var i = 0; i < featureXMLs.length; ++i) {
                var feature = featureXMLs[i];
                var dasFeature = new DASFeature();
                
                dasFeature.segment = segmentID;
                dasFeature.id = feature.getAttribute('id');
                dasFeature.label = feature.getAttribute('label');


/*
                var childNodes = feature.childNodes;
                for (var c = 0; c < childNodes.length; ++c) {
                    var cn = childNodes[c];
                    if (cn.nodeType == Node.ELEMENT_NODE) {
                        var key = cn.tagName;
                        //var val = null;
                        //if (cn.firstChild) {
                        //   val = cn.firstChild.nodeValue;
                        //}
                        dasFeature[key] = 'x';
                    }
                } */


                var spos = elementValue(feature, "START");
                var epos = elementValue(feature, "END");
                if ((spos|0) > (epos|0)) {
                    dasFeature.min = epos|0;
                    dasFeature.max = spos|0;
                } else {
                    dasFeature.min = spos|0;
                    dasFeature.max = epos|0;
                }
                {
                    var tec = feature.getElementsByTagName('TYPE');
                    if (tec.length > 0) {
                        var te = tec[0];
                        if (te.firstChild) {
                            dasFeature.type = te.firstChild.nodeValue;
                        }
                        dasFeature.typeId = te.getAttribute('id');
                        dasFeature.typeCv = te.getAttribute('cvId');
                    }
                }
                dasFeature.type = elementValue(feature, "TYPE");
                if (!dasFeature.type && dasFeature.typeId) {
                    dasFeature.type = dasFeature.typeId; // FIXME?
                }
                
                dasFeature.method = elementValue(feature, "METHOD");
                {
                    var ori = elementValue(feature, "ORIENTATION");
                    if (!ori) {
                        ori = '0';
                    }
                    dasFeature.orientation = ori;
                }
                dasFeature.score = elementValue(feature, "SCORE");
                dasFeature.links = dasLinksOf(feature);
                dasFeature.notes = dasNotesOf(feature);
                
                var groups = feature.getElementsByTagName("GROUP");
                for (var gi  = 0; gi < groups.length; ++gi) {
                    var groupXML = groups[gi];
                    var dasGroup = new DASGroup();
                    dasGroup.type = groupXML.getAttribute('type');
                    dasGroup.id = groupXML.getAttribute('id');
                    dasGroup.links = dasLinksOf(groupXML);
                    dasGroup.notes = dasNotesOf(groupXML);
                    if (!dasFeature.groups) {
                        dasFeature.groups = new Array(dasGroup);
                    } else {
                        dasFeature.groups.push(dasGroup);
                    }
                }

                // Magic notes.  Check with TAD before changing this.
                if (dasFeature.notes) {
                    for (var ni = 0; ni < dasFeature.notes.length; ++ni) {
                        var n = dasFeature.notes[ni];
                        if (n.indexOf('Genename=') == 0) {
                            var gg = new DASGroup();
                            gg.type='gene';
                            gg.id = n.substring(9);
                            if (!dasFeature.groups) {
                                dasFeature.groups = new Array(gg);
                            } else {
                                dasFeature.groups.push(gg);
                            }
                        }
                    }
                }
                
                {
                    var pec = feature.getElementsByTagName('PART');
                    if (pec.length > 0) {
                        var parts = [];
                        for (var pi = 0; pi < pec.length; ++pi) {
                            parts.push(pec[pi].getAttribute('id'));
                        }
                        dasFeature.parts = parts;
                    }
                }
                {
                    var pec = feature.getElementsByTagName('PARENT');
                    if (pec.length > 0) {
                        var parents = [];
                        for (var pi = 0; pi < pec.length; ++pi) {
                            parents.push(pec[pi].getAttribute('id'));
                        }
                        dasFeature.parents = parents;
                    }
                }
                
                features.push(dasFeature);
            }
        }
                
        callback(features, undefined, segmentMap);
    },
    function (err) {
        callback([], err);
    });
}

function DASAlignment(type) {
    this.type = type;
    this.objects = {};
    this.blocks = [];
}

DASSource.prototype.alignments = function(segment, options, callback) {
    var dasURI = this.dasBaseURI + 'alignment?query=' + segment;
    this.doCrossDomainRequest(dasURI, function(responseXML) {
        if (!responseXML) {
            callback([], 'Failed request ' + dasURI);
            return;
        }

        var alignments = [];
        var aliXMLs = responseXML.getElementsByTagName('alignment');
        for (var ai = 0; ai < aliXMLs.length; ++ai) {
            var aliXML = aliXMLs[ai];
            var ali = new DASAlignment(aliXML.getAttribute('alignType'));
            var objXMLs = aliXML.getElementsByTagName('alignObject');
            for (var oi = 0; oi < objXMLs.length; ++oi) {
                var objXML = objXMLs[oi];
                var obj = {
                    id:          objXML.getAttribute('intObjectId'),
                    accession:   objXML.getAttribute('dbAccessionId'),
                    version:     objXML.getAttribute('objectVersion'),
                    dbSource:    objXML.getAttribute('dbSource'),
                    dbVersion:   objXML.getAttribute('dbVersion')
                };
                ali.objects[obj.id] = obj;
            }
            
            var blockXMLs = aliXML.getElementsByTagName('block');
            for (var bi = 0; bi < blockXMLs.length; ++bi) {
                var blockXML = blockXMLs[bi];
                var block = {
                    order:      blockXML.getAttribute('blockOrder'),
                    segments:   []
                };
                var segXMLs = blockXML.getElementsByTagName('segment');
                for (var si = 0; si < segXMLs.length; ++si) {
                    var segXML = segXMLs[si];
                    var seg = {
                        object:      segXML.getAttribute('intObjectId'),
                        min:         segXML.getAttribute('start'),
                        max:         segXML.getAttribute('end'),
                        strand:      segXML.getAttribute('strand'),
                        cigar:       elementValue(segXML, 'cigar')
                    };
                    block.segments.push(seg);
                }
                ali.blocks.push(block);
            }       
                    
            alignments.push(ali);
        }
        callback(alignments);
    });
}


function DASStylesheet() {
/*
    this.highZoomStyles = new Object();
    this.mediumZoomStyles = new Object();
    this.lowZoomStyles = new Object();
*/

    this.styles = [];
}

DASStylesheet.prototype.pushStyle = function(filters, zoom, style) {
    /*

    if (!zoom) {
        this.highZoomStyles[type] = style;
        this.mediumZoomStyles[type] = style;
        this.lowZoomStyles[type] = style;
    } else if (zoom == 'high') {
        this.highZoomStyles[type] = style;
    } else if (zoom == 'medium') {
        this.mediumZoomStyles[type] = style;
    } else if (zoom == 'low') {
        this.lowZoomStyles[type] = style;
    }

    */

    if (!filters) {
        filters = {type: 'default'};
    }
    var styleHolder = shallowCopy(filters);
    if (zoom) {
        styleHolder.zoom = zoom;
    }
    styleHolder.style = style;
    this.styles.push(styleHolder);
}

function DASStyle() {
}

function parseGradient(grad) {
    var steps = grad.getAttribute('steps');
    if (steps) {
        steps = steps|0;
    } else {
        steps = 50;
    }


    var stops = [];
    var colors = [];
    var se = grad.getElementsByTagName('STOP');
    for (var si = 0; si < se.length; ++si) {
        var stop = se[si];
        stops.push(1.0 * stop.getAttribute('score'));
        colors.push(stop.firstChild.nodeValue);
    }

    return makeColourSteps(steps, stops, colors);
}

DASSource.prototype.stylesheet = function(successCB, failureCB) {
    var dasURI, creds = this.credentials;
    if (this.stylesheet_uri) {
        dasURI = this.stylesheet_uri;
        creds = false;
    } else {
        dasURI = this.dasBaseURI + 'stylesheet';
    }

    doCrossDomainRequest(dasURI, function(responseXML) {
        if (!responseXML) {
            if (failureCB) {
                failureCB();
            } 
            return;
        }
        var stylesheet = new DASStylesheet();
        var typeXMLs = responseXML.getElementsByTagName('TYPE');
        for (var i = 0; i < typeXMLs.length; ++i) {
            var typeStyle = typeXMLs[i];
            
            var filter = {};
            filter.type = typeStyle.getAttribute('id'); // Am I right in thinking that this makes DASSTYLE XML invalid?  Ugh.
            filter.label = typeStyle.getAttribute('label');
            filter.method = typeStyle.getAttribute('method');
            var glyphXMLs = typeStyle.getElementsByTagName('GLYPH');
            for (var gi = 0; gi < glyphXMLs.length; ++gi) {
                var glyphXML = glyphXMLs[gi];
                var zoom = glyphXML.getAttribute('zoom');
                var glyph = childElementOf(glyphXML);
                var style = new DASStyle();
                style.glyph = glyph.localName;
                var child = glyph.firstChild;
        
                while (child) {
                    if (child.nodeType == Node.ELEMENT_NODE) {
                        // alert(child.localName);
                        if (child.localName == 'BGGRAD') {
                            style[child.localName] = parseGradient(child);
                        } else {      
                            style[child.localName] = child.firstChild.nodeValue;
                        }
                    }
                    child = child.nextSibling;
                }
                stylesheet.pushStyle(filter, zoom, style);
            }
        }
        successCB(stylesheet);
    }, creds);
}

//
// sources command
// 

function DASRegistry(uri, opts)
{
    opts = opts || {};
    this.uri = uri;
    this.opts = opts;   
}

DASRegistry.prototype.sources = function(callback, failure, opts)
{
    if (!opts) {
        opts = {};
    }

    var filters = [];
    if (opts.taxon) {
        filters.push('organism=' + opts.taxon);
    }
    if (opts.auth) {
        filters.push('authority=' + opts.auth);
    }
    if (opts.version) {
        filters.push('version=' + opts.version);
    }
    var quri = this.uri;
    if (filters.length > 0) {
        quri = quri + '?' + filters.join('&');   // '&' as a separator to hack around dasregistry.org bug.
    }

    doCrossDomainRequest(quri, function(responseXML) {
        if (!responseXML && failure) {
            failure();
            return;
        }

        var sources = [];       
        var sourceXMLs = responseXML.getElementsByTagName('SOURCE');
        for (var si = 0; si < sourceXMLs.length; ++si) {
            var sourceXML = sourceXMLs[si];
            var versionXMLs = sourceXML.getElementsByTagName('VERSION');
            if (versionXMLs.length < 1) {
                continue;
            }
            var versionXML = versionXMLs[0];

            var coordXMLs = versionXML.getElementsByTagName('COORDINATES');
            var coords = [];
            for (var ci = 0; ci < coordXMLs.length; ++ci) {
                var coordXML = coordXMLs[ci];
                var coord = new DASCoords();
                coord.auth = coordXML.getAttribute('authority');
                coord.taxon = coordXML.getAttribute('taxid');
                coord.version = coordXML.getAttribute('version');
                coords.push(coord);
            }
            
            var caps = [];
            var capXMLs = versionXML.getElementsByTagName('CAPABILITY');
            var uri;
            for (var ci = 0; ci < capXMLs.length; ++ci) {
                var capXML = capXMLs[ci];
                
                caps.push(capXML.getAttribute('type'));

                if (capXML.getAttribute('type') == 'das1:features') {
                    var fep = capXML.getAttribute('query_uri');
                    uri = fep.substring(0, fep.length - ('features'.length));
                }
            }

            var props = {};
            var propXMLs = versionXML.getElementsByTagName('PROP');
            for (var pi = 0; pi < propXMLs.length; ++pi) {
                pusho(props, propXMLs[pi].getAttribute('name'), propXMLs[pi].getAttribute('value'));
            }
            
            if (uri) {
                var source = new DASSource(uri, {
                    source_uri: sourceXML.getAttribute('uri'),
                    name:  sourceXML.getAttribute('title'),
                    desc:  sourceXML.getAttribute('description'),
                    coords: coords,
                    props: props,
                    capabilities: caps
                });
                sources.push(source);
            }
        }
        
        callback(sources);
    });
}


//
// Utility functions
//

function elementValue(element, tag)
{
    var children = element.getElementsByTagName(tag);
    if (children.length > 0 && children[0].firstChild) {
        var c = children[0];
        if (c.childNodes.length == 1) {
            return c.firstChild.nodeValue;
        } else {
            var s = '';
            for (var ni = 0; ni < c.childNodes.length; ++ni) {
                s += c.childNodes[ni].nodeValue;
            }
            return s;
        }

    } else {
        return null;
    }
}

function childElementOf(element)
{
    if (element.hasChildNodes()) {
        var child = element.firstChild;
        do {
            if (child.nodeType == Node.ELEMENT_NODE) {
                return child;
            } 
            child = child.nextSibling;
        } while (child != null);
    }
    return null;
}


function dasLinksOf(element)
{
    var links = new Array();
    var maybeLinkChilden = element.getElementsByTagName('LINK');
    for (var ci = 0; ci < maybeLinkChilden.length; ++ci) {
        var linkXML = maybeLinkChilden[ci];
        if (linkXML.parentNode == element) {
            links.push(new DASLink(linkXML.firstChild ? linkXML.firstChild.nodeValue : 'Unknown', linkXML.getAttribute('href')));
        }
    }
    
    return links;
}

function dasNotesOf(element)
{
    var notes = [];
    var maybeNotes = element.getElementsByTagName('NOTE');
    for (var ni = 0; ni < maybeNotes.length; ++ni) {
        if (maybeNotes[ni].firstChild) {
            notes.push(maybeNotes[ni].firstChild.nodeValue);
        }
    }
    return notes;
}

function doCrossDomainRequest(url, handler, credentials, custAuth) {
    // TODO: explicit error handlers?

    if (window.XDomainRequest) {
        var req = new XDomainRequest();
        req.onload = function() {
            var dom = new ActiveXObject("Microsoft.XMLDOM");
            dom.async = false;
            dom.loadXML(req.responseText);
            handler(dom);
        }
        req.open("get", url);
        req.send('');
    } else {
        var reqStart = Date.now();
        var req = new XMLHttpRequest();

        req.onreadystatechange = function() {
            if (req.readyState == 4) {
              if (req.status >= 200 || req.status == 0) {
                  handler(req.responseXML, req);
              }
            }
        };
        req.open("get", url, true);
        if (credentials) {
            req.withCredentials = true;
        }
        if (custAuth) {
            req.setRequestHeader('X-DAS-Authorisation', custAuth);
        }
        req.setRequestHeader('Accept', 'application/xml,*/*');
        req.send('');
    }
}

DASSource.prototype.doCrossDomainRequest = function(url, handler, errHandler) {
    var custAuth;
    if (this.xUser) {
        custAuth = 'Basic ' + btoa(this.xUser + ':' + this.xPass);
    }

    try {
        return doCrossDomainRequest(url, handler, this.credentials, custAuth);
    } catch (err) {
        if (errHandler) {
            errHandler(err);
        } else {
            throw err;
        }
    }
}
/* -*- mode: javascript; c-basic-offset: 4; indent-tabs-mode: nil -*- */

// 
// Dalliance Genome Explorer
// (c) Thomas Down 2006-2010
//
// domui.js: SVG UI components
//

"use strict";

Browser.prototype.removeAllPopups = function() {
    removeChildren(this.hPopupHolder);
    removeChildren(this.popupHolder);
}

Browser.prototype.makeTooltip = function(ele, text)
{
    var isin = false;
    var thisB = this;
    var timer = null;
    var outlistener;
    outlistener = function(ev) {
        isin = false;
        if (timer) {
            clearTimeout(timer);
            timer = null;
        }
        ele.removeEventListener('mouseout', outlistener, false);
    };

    var setup = function(ev) {
        var mx = ev.clientX + window.scrollX, my = ev.clientY + window.scrollY;
        if (!timer) {
            timer = setTimeout(function() {
                var ttt;
                if (typeof(text) === 'function') {
                    ttt = text();
                } else {
                    ttt = text;
                }

                var popup = makeElement('div',
                    [makeElement('div', null, {className: 'tooltip-arrow'}),
                     makeElement('div', ttt, {className: 'tooltip-inner'})], 
                    {className: 'tooltip bottom in'}, {
                    display: 'block',
                    top: '' + (my + 20) + 'px',
                    left: '' + Math.max(mx - 30, 20) + 'px'
                });
                thisB.hPopupHolder.appendChild(popup);
                var moveHandler;
                moveHandler = function(ev) {
                    try {
                        thisB.hPopupHolder.removeChild(popup);
                    } catch (e) {
                        // May have been removed by other code which clears the popup layer.
                    }
                    window.removeEventListener('mousemove', moveHandler, false);
                    if (isin) {
                        if (ele.offsetParent == null) {
                            // dlog('Null parent...');
                        } else {
                            setup(ev);
                        }
                    }
                }
                window.addEventListener('mousemove', moveHandler, false);
                timer = null;
            }, 1000);
        }
    };

    ele.addEventListener('mouseover', function(ev) {
        isin = true
        ele.addEventListener('mouseout', outlistener, false);
        setup(ev);
    }, false);
    ele.addEventListener('DOMNodeRemovedFromDocument', function(ev) {
        isin = false;
        if (timer) {
            clearTimeout(timer);
            timer = null;
        }
    }, false);
}

Browser.prototype.popit = function(ev, name, ele, opts)
{
    var thisB = this;
    if (!opts) 
        opts = {};
    if (!ev) 
        ev = {};

    var width = opts.width || 200;

    var mx, my;

    if (ev.clientX) {
        var mx =  ev.clientX, my = ev.clientY;
    } else {
        mx = 500; my= 50;
    }
    mx +=  document.documentElement.scrollLeft || document.body.scrollLeft;
    my +=  document.documentElement.scrollTop || document.body.scrollTop;
    var winWidth = window.innerWidth;

    var top = my;
    var left = Math.min(mx - (width/2) - 4, (winWidth - width - 30));

    var popup = makeElement('div');
    popup.className = 'popover fade ' + (ev.clientX ? 'bottom ' : '') + 'in';
    popup.style.display = 'block';
    popup.style.position = 'absolute';
    popup.style.top = '' + top + 'px';
    popup.style.left = '' + left + 'px';
    popup.style.width = width + 'px';
    if (width > 276) {
        // HACK Bootstrappification...
        popup.style.maxWidth = width + 'px';
    }

    popup.appendChild(makeElement('div', null, {className: 'arrow'}));

    if (name) {
        var closeButton = makeElement('button', '', {className: 'close'});
        closeButton.innerHTML = '&times;'

        closeButton.addEventListener('mouseover', function(ev) {
            closeButton.style.color = 'red';
        }, false);
        closeButton.addEventListener('mouseout', function(ev) {
            closeButton.style.color = 'black';
        }, false);
        closeButton.addEventListener('click', function(ev) {
            ev.preventDefault(); ev.stopPropagation();
            thisB.removeAllPopups();
        }, false);
        var tbar = makeElement('h4', [makeElement('span', name, null, {maxWidth: '200px'}), closeButton], {/*className: 'popover-title' */}, {paddingLeft: '10px', paddingRight: '10px'});

        var dragOX, dragOY;
        var moveHandler, upHandler;
        moveHandler = function(ev) {
            ev.stopPropagation(); ev.preventDefault();
            left = left + (ev.clientX - dragOX);
            if (left < 8) {
                left = 8;
            } if (left > (winWidth - width - 32)) {
                left = (winWidth - width - 26);
            }
            top = top + (ev.clientY - dragOY);
            top = Math.max(10, top);
            popup.style.top = '' + top + 'px';
            popup.style.left = '' + Math.min(left, (winWidth - width - 10)) + 'px';
            dragOX = ev.clientX; dragOY = ev.clientY;
        }
        upHandler = function(ev) {
            ev.stopPropagation(); ev.preventDefault();
            window.removeEventListener('mousemove', moveHandler, false);
            window.removeEventListener('mouseup', upHandler, false);
        }
        tbar.addEventListener('mousedown', function(ev) {
            ev.preventDefault(); ev.stopPropagation();
            dragOX = ev.clientX; dragOY = ev.clientY;
            window.addEventListener('mousemove', moveHandler, false);
            window.addEventListener('mouseup', upHandler, false);
        }, false);
                              

        popup.appendChild(tbar);
    }

    popup.appendChild(makeElement('div', ele, {className: 'popover-content'}, {
        padding: '0px'
    }));
    this.hPopupHolder.appendChild(popup);

    var popupHandle = {
        node: popup,
        displayed: true
    };
    popup.addEventListener('DOMNodeRemoved', function(ev) {
        if (ev.target == popup) {
            popupHandle.displayed = false;
        }
    }, false);
    return popupHandle;
}

function makeTreeTableSection(title, content, visible) {
    var ttButton = makeElement('i');
    function update() {
        if (visible) {
            ttButton.className = 'fa fa-caret-down';
            content.style.display = 'table';
        } else {
            ttButton.className = 'fa fa-caret-right';
            content.style.display = 'none';
        }
    }
    update();

    ttButton.addEventListener('click', function(ev) {
        ev.preventDefault(); ev.stopPropagation();
        visible = !visible;
        update();
    }, false);

    var heading = makeElement('h6', [ttButton, ' ', title], {}, {display: 'block', background: 'gray', color: 'white', width: '100%', padding: '5px 2px', margin: '0px'});
    return makeElement('div', [heading, content], {});
}

function dlog(msg) {
    console.log(msg);
}
// 
// Dalliance Genome Explorer
// (c) Thomas Down 2006-2010
//
// feature-draw.js: new feature-tier renderer
//

"use strict";

var MIN_PADDING = 3;
var DEFAULT_SUBTIER_MAX = 100;

function isDasBooleanTrue(s) {
    s = ('' + s).toLowerCase();
    return s==='yes' || s==='true';
}

function isDasBooleanNotFalse(s) {
    if (!s)
        return false;

    s = ('' + s).toLowerCase();
    return s!=='no' || s!=='false';
}

function SubTier() {
    this.glyphs = [];
    this.height = 0;
    this.quant = null;
}

SubTier.prototype.indexFor = function(glyph) {
    var gmin = glyph.min();
    var lb = 0, ub = this.glyphs.length;
    while (ub > lb) {
        var mid = ((lb + ub)/2)|0;
        if (mid >= this.glyphs.length)
            return this.glyphs.length;
        var mg = this.glyphs[mid];
        if (gmin < mg.min()) {
            ub = mid;
        } else {
            lb = mid + 1;
        }
    }
    return ub;
}

SubTier.prototype.add = function(glyph) {
    var ind = this.indexFor(glyph);
    this.glyphs.splice(ind, 0, glyph);
    this.height = Math.max(this.height, glyph.height());
    if (glyph.quant && this.quant == null) {
        this.quant = glyph.quant;
    }
}

SubTier.prototype.hasSpaceFor = function(glyph) {
    var ind = this.indexFor(glyph);
    if (ind > 0 && this.glyphs[ind-1].max() >= glyph.min())
        return false;
    if (ind < this.glyphs.length && this.glyphs[ind].min() <= glyph.max())
        return false;

    return true;
}

var GLOBAL_GC;

function drawFeatureTier(tier)
{
    var start = Date.now()|0;
    GLOBAL_GC = tier.viewport.getContext('2d'); // Should only be used for metrics.
    sortFeatures(tier);

    var glyphs = [];
    var specials = false;

    // group by style
    var gbsFeatures = {};
    var gbsStyles = {};

    for (var uft in tier.ungroupedFeatures) {
        var ufl = tier.ungroupedFeatures[uft];
        
        for (var pgid = 0; pgid < ufl.length; ++pgid) {
            var f = ufl[pgid];
            if (f.parts) {  // FIXME shouldn't really be needed
                continue;
            }

            var style = tier.styleForFeature(f);
            if (!style)
                continue;

            if (style.glyph == 'LINEPLOT') {
                pusho(gbsFeatures, style.id, f);
                gbsStyles[style.id] = style;
            } else {
                var g = glyphForFeature(f, 0, style, tier);
                if (g)
                    glyphs.push(g);
            }
        }
    }

    for (var gbs in gbsFeatures) {
        var gf = gbsFeatures[gbs];
        var style = gbsStyles[gbs];
        if (style.glyph == 'LINEPLOT') {
            glyphs.push(makeLineGlyph(gf, style, tier));
            specials = true;
        }
    }

    // Merge supergroups    

    if (tier.dasSource.collapseSuperGroups && !tier.bumped) {
        for (var sg in tier.superGroups) {
            var sgg = tier.superGroups[sg];
            tier.groups[sg] = shallowCopy(tier.groups[sg]);
            tier.groups[sg].type = tier.groups[sgg[0]].type;   // HACK to make styling easier in DAS1.6
            var featsByType = {};

            var sgMin = 10000000000, sgMax = -10000000000;
            var sgSeg = null;
            for (var g = 0; g < sgg.length; ++g) {
                var gf = tier.groupedFeatures[sgg[g]];
                if (!gf)
                    continue;

                for (var fi = 0; fi < gf.length; ++fi) {
                    var f = gf[fi];
                    pusho(featsByType, f.type, f);
                    sgMin = Math.min(f.min, sgMin);
                    sgMax = Math.max(f.max, sgMax);
                    if (f.segment && !sgSeg)
                        sgSeg = f.segment;
                }

                if (tier.groups[sg] && !tier.groups[sg].links || tier.groups[sg].links.length == 0) {
                   tier.groups[sg].links = tier.groups[sgg[0]].links;
                }

                delete tier.groupedFeatures[sgg[g]];  // 'cos we don't want to render the unmerged version.
            }

            tier.groups[sg].max = sgMax;
            tier.groups[sg].min = sgMin;
            tier.groups[sg].segment = sgSeg;

            for (var t in featsByType) {
                var feats = featsByType[t];
                var template = feats[0];
                var loc = null;
                for (var fi = 0; fi < feats.length; ++fi) {
                    var f = feats[fi];
                    var fl = new Range(f.min, f.max);
                    if (!loc) {
                        loc = fl;
                    } else {
                        loc = union(loc, fl);
                    }
                }
                var mergedRanges = loc.ranges();
                for (var si = 0; si < mergedRanges.length; ++si) {
                    var r = mergedRanges[si];

                    // begin coverage-counting
                    var posCoverage = ((r.max()|0) - (r.min()|0) + 1) * sgg.length;
                    var actCoverage = 0;
                    for (var fi = 0; fi < feats.length; ++fi) {
                        var f = feats[fi];
                        if ((f.min|0) <= r.max() && (f.max|0) >= r.min()) {
                            var umin = Math.max(f.min|0, r.min());
                            var umax = Math.min(f.max|0, r.max());
                            actCoverage += (umax - umin + 1);
                        }
                    }
                    var visualWeight = ((1.0 * actCoverage) / posCoverage);
                    // end coverage-counting

                    var newf = new DASFeature();
                    for (var k in template) {
                        newf[k] = template[k];
                    }
                    newf.min = r.min();
                    newf.max = r.max();
                    if (newf.label && sgg.length > 1) {
                        newf.label += ' (' + sgg.length + ' vars)';
                    }
                    newf.visualWeight = ((1.0 * actCoverage) / posCoverage);
                    pusho(tier.groupedFeatures, sg, newf);
                    // supergroups are already in tier.groups.
                }
            }

            delete tier.superGroups[sg]; // Do we want this?
        }       
    }

    // Glyphify groups.

    var gl = new Array();
    for (var gid in tier.groupedFeatures) {
        gl.push(gid);
    }
    gl.sort(function(g1, g2) {
        var d = tier.groupedFeatures[g1][0].score - tier.groupedFeatures[g2][0].score;
        if (d > 0) {
            return -1;
        } else if (d == 0) {
            return 0;
        } else {
            return 1;
        }
    });

    var groupGlyphs = {};
    for (var gx = 0; gx < gl.length; ++gx) {
        var gid = gl[gx];
        var g = glyphsForGroup(tier.groupedFeatures[gid], 0, tier.groups[gid], tier,
                               (tier.dasSource.collapseSuperGroups && !tier.bumped) ? 'collapsed_gene' : 'tent');
        if (g) {
            g.group = tier.groups[gid];
            groupGlyphs[gid] = g;
        }
    }

    for (var sg in tier.superGroups) {
        var sgg = tier.superGroups[sg];
        var sgGlyphs = [];
        var sgMin = 10000000000;
        var sgMax = -10000000000;
        for (var sgi = 0; sgi < sgg.length; ++sgi) {
            var gg = groupGlyphs[sgg[sgi]];
            groupGlyphs[sgg[sgi]] = null;
            if (gg) {
                sgGlyphs.push(gg);
                sgMin = Math.min(sgMin, gg.min());
                sgMax = Math.max(sgMax, gg.max());
            }
        }
        for (var sgi = 0; sgi < sgGlyphs.length; ++sgi) {
            var gg = sgGlyphs[sgi];
            glyphs.push(new PaddedGlyph(gg, sgMin, sgMax));
        }
    }
    for (var g in groupGlyphs) {
        var gg = groupGlyphs[g];
        if (gg) {
            glyphs.push(gg);
        }
    }

    // Bumping

    var unbumpedST = new SubTier();
    var bumpedSTs = [];
    var hasBumpedFeatures = false;
    var subtierMax = tier.dasSource.subtierMax || DEFAULT_SUBTIER_MAX;
    
  GLYPH_LOOP:
    for (var i = 0; i < glyphs.length; ++i) {
        var g = glyphs[i];
        if (g.bump) {
            hasBumpedFeatures = true;
        }
        if (g.bump && (tier.bumped || tier.dasSource.collapseSuperGroups)) {       // kind-of nasty.  supergroup collapsing is different from "normal" unbumping
            for (var sti = 0; sti < bumpedSTs.length;  ++sti) {
                var st = bumpedSTs[sti];
                if (st.hasSpaceFor(g)) {
                    st.add(g);
                    continue GLYPH_LOOP;
                }
            }
            if (bumpedSTs.length >= subtierMax) {
                // tier.status = 'Too many overlapping features, truncating at ' + subtierMax;
            } else {
                var st = new SubTier();
                st.add(g);
                bumpedSTs.push(st);
            }
        } else {
            unbumpedST.add(g);
        }
    }

    if (unbumpedST.glyphs.length > 0) {
        bumpedSTs = [unbumpedST].concat(bumpedSTs);
    }

    for (var sti = 0; sti < bumpedSTs.length; ++sti) {
        var st = bumpedSTs[sti];
        if (st.quant) {
            st.glyphs.unshift(new GridGlyph(st.height));
        }
    }

    for (var sti = 0; sti < bumpedSTs.length; ++sti) {
        var st = bumpedSTs[sti];
        st.glyphs.sort(function (g1, g2) {
            var z1 = g1.zindex || 0;
            var z2 = g2.zindex || 0;
            return z1 - z2;
        });
    }

    tier.subtiers = bumpedSTs;
    tier.glyphCacheOrigin = tier.browser.viewStart;

    var end = Date.now()|0;
    // console.log('dft took ' + (end-start) + 'ms');
}

function formatQuantLabel(v) {
    var t = '' + v;
    var dot = t.indexOf('.');
    if (dot < 0) {
        return t;
    } else {
        var dotThreshold = 2;
        if (t.substring(0, 1) == '-') {
            ++dotThreshold;
        }

        if (dot >= dotThreshold) {
            return t.substring(0, dot);
        } else {
            return t.substring(0, dot + 2);
        }
    }
}

DasTier.prototype.paint = function() {
    var start = Date.now()|0;
    var retina = this.browser.retina && window.devicePixelRatio > 1;

    var subtiers = this.subtiers;
    if (!subtiers) {
	   return;
    }

    var desiredWidth = this.browser.featurePanelWidth + 2000;
    if (retina) {
        desiredWidth *= 2;
    }
    var fpw = this.viewport.width|0;
    if (fpw < desiredWidth - 50) {
        this.viewport.width = fpw = desiredWidth;
    }

    var lh = MIN_PADDING;
    for (var s = 0; s < subtiers.length; ++s) {
        lh = lh + subtiers[s].height + MIN_PADDING;
    }
    lh += 6
    lh = Math.max(lh, this.browser.minTierHeight);

    var canvasHeight = lh;
    if (retina) {
        canvasHeight *= 2;
    }

    if (canvasHeight != this.viewport.height) {
        this.viewport.height = canvasHeight;
    }
    
    var tierHeight = Math.max(lh, this.browser.minTierHeight);
    this.viewport.style.left = '-1000px';
    this.viewport.style.width = retina ? ('' + (fpw/2) + 'px') : ('' + fpw + 'px');
    this.viewport.style.height = '' + lh + 'px';
    this.layoutHeight =  Math.max(lh, this.browser.minTierHeight);

    this.updateHeight();
    this.drawOverlay();
    this.norigin = this.browser.viewStart;

    var gc = this.viewport.getContext('2d');
    gc.clearRect(0, 0, fpw, canvasHeight);

    gc.save();
    if (retina) {
        gc.scale(2, 2);
    }

    if (this.background) {
        gc.fillStyle = this.background;

        if (this.knownCoverage) {
            var knownRanges = this.knownCoverage.ranges();
            for (var ri = 0; ri < knownRanges.length; ++ri) {
                var r = knownRanges[ri];
                var knownMin = (r.min() - this.browser.viewStart) * this.browser.scale + 1000;
                var knownMax = (r.max() - this.browser.viewStart) * this.browser.scale + 1000;
                gc.fillRect(knownMin, 0, knownMax - knownMin, lh);
            }
        }
    }

    var offset = ((this.glyphCacheOrigin - this.browser.viewStart)*this.browser.scale)+1000;
    gc.translate(offset, MIN_PADDING);
   
    for (var s = 0; s < subtiers.length; ++s) {
        var quant = null;
        var glyphs = subtiers[s].glyphs;
        for (var i = 0; i < glyphs.length; ++i) {
            var glyph = glyphs[i];
            if (glyph.min() < fpw-offset && glyph.max() > -offset) { 
                var glyph = glyphs[i];
                glyph.draw(gc);
                if (glyph.quant) {
                    quant = glyph.quant;
                }
            }
        }
        gc.translate(0, subtiers[s].height + MIN_PADDING);
    }
    gc.restore();

    if (quant && this.quantLeapThreshold && this.featureSource && sourceAdapterIsCapable(this.featureSource, 'quantLeap')) {
        var ry = 3 + subtiers[0].height * (1.0 - ((this.quantLeapThreshold - quant.min) / (quant.max - quant.min)));

        gc.save();
        if (retina)
            gc.scale(2, 2);
        gc.strokeStyle = 'red';
        gc.lineWidth = 0.3;
        gc.beginPath();
        gc.moveTo(0, ry);
        gc.lineTo(5000, ry);
        gc.stroke();
        gc.restore();
    }

    this.paintQuant();

    var end = Date.now()|0;
    // console.log('paint took ' + (end-start) + 'ms');
}

DasTier.prototype.paintQuant = function() {
    if (!this.quantOverlay)
        return;

    var retina = this.browser.retina && window.devicePixelRatio > 1;

    var quant;
    if (this.subtiers && this.subtiers.length > 0)
        quant = this.subtiers[0].quant;

    if (quant) {
        var h = this.subtiers[0].height;
        var w = 50;
        this.quantOverlay.height = this.viewport.height;
        this.quantOverlay.width = retina ? w*2 : w;
        this.quantOverlay.style.height = '' + (retina ? this.quantOverlay.height/2 : this.quantOverlay.height) + 'px';
        this.quantOverlay.style.width = '' + w + 'px';
        this.quantOverlay.style.display = 'block';
        var ctx = this.quantOverlay.getContext('2d');
        if (retina)
            ctx.scale(2, 2);

        var numTics = 2;
        if (h > 40) {
            numTics = 1 + ((h/20) | 0);
        }
        var ticSpacing = (h+MIN_PADDING*2) / (numTics - 1);
        var ticInterval = (quant.max - quant.min) / (numTics - 1);

        ctx.fillStyle = 'white'
        ctx.globalAlpha = 0.6;
        if (this.browser.rulerLocation == 'right') {
            ctx.fillRect(w-30, 0, 30, h + MIN_PADDING*2);
        } else {
            ctx.fillRect(0, 0, 30, h + MIN_PADDING*2);
        }
        ctx.globalAlpha = 1.0;

        ctx.strokeStyle = 'black';
        ctx.lineWidth = 1;
        ctx.beginPath();

        if (this.browser.rulerLocation == 'right') {
            ctx.moveTo(w - 8, MIN_PADDING);
            ctx.lineTo(w, MIN_PADDING);
            ctx.lineTo(w, h + MIN_PADDING);
            ctx.lineTo(w - 8, h + MIN_PADDING);
            for (var t = 1; t < numTics-1; ++t) {
                var ty = t*ticSpacing;
                ctx.moveTo(w, ty);
                ctx.lineTo(w - 5, ty);
            }
        } else {
            ctx.moveTo(8, MIN_PADDING);
            ctx.lineTo(0, MIN_PADDING);
            ctx.lineTo(0, h+MIN_PADDING);
            ctx.lineTo(8, h+MIN_PADDING);
            for (var t = 1; t < numTics-1; ++t) {
                var ty = t*ticSpacing;
                ctx.moveTo(0, ty);
                ctx.lineTo(5, ty);
            }
        }
        ctx.stroke();

        ctx.fillStyle = 'black';

        if (this.browser.rulerLocation == 'right') {
            ctx.textAlign = 'right';
            ctx.fillText(formatQuantLabel(quant.max), w-9, 8);
            ctx.fillText(formatQuantLabel(quant.min), w-9, h + MIN_PADDING);
            for (var t = 1; t < numTics-1; ++t) {
                var ty = t*ticSpacing;
                ctx.fillText(formatQuantLabel((1.0*quant.max) - (t*ticInterval)), w - 9, ty + 3);
            }
        } else {
            ctx.textAlign = 'left';
            ctx.fillText(formatQuantLabel(quant.max), 9, 8);
            ctx.fillText(formatQuantLabel(quant.min), 9, h + MIN_PADDING);
            for (var t = 1; t < numTics-1; ++t) {
                var ty = t*ticSpacing;
                ctx.fillText(formatQuantLabel((1.0*quant.max) - (t*ticInterval)), 9, ty + 3);
            }
        }
    } else {
        this.quantOverlay.style.display = 'none';
    }
}

function glyphsForGroup(features, y, groupElement, tier, connectorType) {
    var gstyle = tier.styleForFeature(groupElement);
    var label;
    var labelWanted = false;

    var glyphs = [];
    var strand = null;
    for (var i = 0; i < features.length; ++i) {
        var f = features[i];
        if (f.orientation && strand==null) {
            strand = f.orientation;
        }
         if (!label && f.label) {
            label = f.label;
        }

        var style = tier.styleForFeature(f);
        if (!style) {
            continue;
        }
        if (f.parts) {  // FIXME shouldn't really be needed
            continue;
        }
        if (isDasBooleanTrue(style.LABEL))
            labelWanted = true;

        var g = glyphForFeature(f, 0, style, tier, null, true);
        if (g) {
            glyphs.push(g);
        }
    }

    if (glyphs.length == 0)
        return null;
    
    var connector = 'flat';
    if (tier.dasSource.collapseSuperGroups && !tier.bumped) {
        if (strand === '+') {
            connector = 'collapsed+';
        } else if (strand === '-') {
            connector = 'collapsed-';
        }
    } else {
        if (strand === '+') {
            connector = 'hat+';
        } else if (strand === '-') {
            connector = 'hat-';
        }
    }

    var labelText = null;
    if ((label && labelWanted) || (gstyle && (isDasBooleanTrue(gstyle.LABEL) || isDasBooleanTrue(gstyle.LABELS)))) {  // HACK, LABELS should work.
        labelText = groupElement.label || label;
        var sg = tier.groupsToSupers[groupElement.id];
        if (sg && tier.superGroups[sg]) {    // workaround case where group and supergroup IDs match.
            //if (groupElement.id != tier.superGroups[sg][0]) {
            //    dg.label = null;
            // }
        }
    }

    var gg = new GroupGlyph(glyphs, connector);
    if (labelText) {
        if (strand === '+') {
            labelText = '>' + labelText;
        } else if (strand === '-') {
            labelText = '<' + labelText;
        }
        gg = new LabelledGlyph(gg, labelText, false);
    }
    gg.bump = true;
    return gg;
}

function glyphForFeature(feature, y, style, tier, forceHeight, noLabel)
{
    var scale = tier.browser.scale, origin = tier.browser.viewStart;
    var gtype = style.glyph || 'BOX';
    var glyph;

    var min = feature.min;
    var max = feature.max;
    var type = feature.type;
    var strand = feature.orientation;
    var score = feature.score;
    var label = feature.label || feature.id;

    var minPos = (min - origin) * scale;
    var rawMaxPos = ((max - origin + 1) * scale);
    var maxPos = Math.max(rawMaxPos, minPos + 1);

    var height = tier.forceHeight || style.HEIGHT || forceHeight || 12;
    var requiredHeight = height = 1.0 * height;
    var bump = style.BUMP && isDasBooleanTrue(style.BUMP);

    var gg, quant;

    if (gtype === 'CROSS' || gtype === 'EX' || gtype === 'TRIANGLE' || gtype === 'DOT' || gtype === 'SQUARE' || gtype === 'STAR' || gtype === 'PLIMSOLL') {
        var stroke = style.FGCOLOR || 'black';
        var fill = style.BGCOLOR || 'none';
        var outline = style.STROKECOLOR;

        if (style.BGITEM && feature.itemRgb) {
            stroke = feature.itemRgb;
        } else if (isDasBooleanTrue(style.COLOR_BY_SCORE2)) {
            var grad = style.BGGRAD || style._gradient;
            if (!grad) {
                grad = makeGradient(50, style.COLOR1, style.COLOR2, style.COLOR3);
                style._gradient = grad;
            }

            var sc2 = feature.score2;
            if (sc2 != undefined || !stroke) {
                sc2 = sc2 || 0;

                var smin2 = style.MIN2 ? (1.0 * style.MIN2) : 0.0;
                var smax2 = style.MAX2 ? (1.0 * style.MAX2) : 1.0;
                var relScore2 = ((1.0 * sc2) - smin2) / (smax2-smin2);

                var step = (relScore2*grad.length)|0;
                if (step < 0) step = 0;
                if (step >= grad.length) step = grad.length - 1;
                stroke = grad[step];
            }
        }



        var height = tier.forceHeight || style.HEIGHT || forceHeight || 12;
        requiredHeight = height = 1.0 * height;

        var size = style.SIZE || height;
        if (style.RSIZE) {
            size = (1.0 * style.RSIZE) * height;
        }

        if (style.STROKETHRESHOLD) {
            if (size < (1.0 * style.STROKETHRESHOLD))
                outline = null;
        }
        
        size = 1.0 * size;

        var mid = (minPos + maxPos)/2;
        var hh = size/2;

        var mark;
        var bMinPos = minPos, bMaxPos = maxPos;

        if (gtype === 'EX') {
            gg = new ExGlyph(mid, size, stroke);
        } else if (gtype === 'TRIANGLE') {
            var dir = style.DIRECTION || 'N';
            var width = style.LINEWIDTH || size;
            gg = new TriangleGlyph(mid, size, dir, width, stroke, outline);
        } else if (gtype === 'DOT') {
            gg = new DotGlyph(mid, size, stroke, outline);
        } else if (gtype === 'PLIMSOLL') {
            gg = new PlimsollGlyph(mid, size, 0.2 * size, stroke, outline);
        } else if (gtype === 'SQUARE') {
            gg = new BoxGlyph(mid - hh, 0, size, size, stroke, outline);
        } else if (gtype === 'STAR') {
            var points = 5;
            if (style.POINTS) 
                points = style.POINTS | 0;
            gg = new StarGlyph(mid, hh, points, stroke, outline);
        } else {
            gg = new CrossGlyph(mid, size, stroke);
        }

        if (fill && fill != 'none' && (maxPos - minPos) > 5) {
            var bgg = new BoxGlyph(minPos, 0, (maxPos - minPos), size, fill);
            gg = new GroupGlyph([bgg, gg]);
        }

        if (isDasBooleanTrue(style.SCATTER)) {
            var smin = tier.quantMin(style);
            var smax = tier.quantMax(style);

            if (!smax) {
                if (smin < 0) {
                    smax = 0;
                } else {
                    smax = 10;
                }
            }
            if (!smin) {
                smin = 0;
            }

            var relScore = ((1.0 * score) - smin) / (smax-smin);
            var relOrigin = (-1.0 * smin) / (smax - smin);

            if (relScore < 0.0 || relScore > 1.0) {
                // Glyph is out of bounds.
                // Should we allow for "partially showing" glyphs?

                return null;
            } else {
                if (relScore >= relOrigin) {
                    height = Math.max(1, (relScore - relOrigin) * requiredHeight);
                    y = y + ((1.0 - relOrigin) * requiredHeight) - height;
                } else {
                    height = Math.max(1, (relScore - relOrigin) * requiredHeight);
                    y = y + ((1.0 - relOrigin) * requiredHeight);
                }
                
                quant = {min: smin, max: smax};

                var heightFudge = 0;
                var featureLabel;
                if (typeof(feature.forceLabel) !== 'undefined')
                    featureLabel = feature.forceLabel;
                else
                    featureLabel = style.LABEL;

                if (isDasBooleanNotFalse(featureLabel) && label && !noLabel) {
                    gg = new LabelledGlyph(gg, label, true, null, featureLabel == 'above' ? 'above' : 'below');
                    if (featureLabel == 'above') {
                        heightFudge = gg.textHeight + 2;
                    }
                    noLabel = true;
                }
                gg = new TranslatedGlyph(gg, 0, y - hh - heightFudge, requiredHeight);
            }
        }
    } else if (gtype === 'HISTOGRAM' || gtype === 'GRADIENT' && score !== 'undefined') {
        var smin = tier.quantMin(style);
        var smax = tier.quantMax(style);

        if (!smax) {
            if (smin < 0) {
                smax = 0;
            } else {
                smax = 10;
            }
        }
        if (!smin) {
            smin = 0;
        }

        if ((1.0 * score) < (1.0 *smin)) {
            score = smin;
        }
        if ((1.0 * score) > (1.0 * smax)) {
            score = smax;
        }
        var relScore = ((1.0 * score) - smin) / (smax-smin);
        var relOrigin = (-1.0 * smin) / (smax - smin);

        if (gtype === 'HISTOGRAM') {
            if (relScore >= relOrigin) {
                height = Math.max(1, (relScore - relOrigin) * requiredHeight);
                y = y + ((1.0 - relOrigin) * requiredHeight) - height;
            } else {
                height = Math.max(1, (relOrigin - relScore) * requiredHeight);
                y = y + ((1.0 - relOrigin) * requiredHeight);
            }
            quant = {min: smin, max: smax};
        }

        var stroke = style.FGCOLOR || null;
        var fill = style.BGCOLOR || style.COLOR1 || 'green';
        if (style.BGITEM && feature.itemRgb)
            fill = feature.itemRgb;
        var alpha = style.ALPHA ? (1.0 * style.ALPHA) : null;

        if (style.BGGRAD) {
            var grad = style.BGGRAD;
            var step = (relScore*grad.length)|0;
            if (step < 0) step = 0;
            if (step >= grad.length) step = grad.length - 1;
            fill = grad[step];
        }
        if (style.COLOR2) {
            var grad = style._gradient;
            if (!grad) {
                grad = makeGradient(50, style.COLOR1, style.COLOR2, style.COLOR3);
                style._gradient = grad;
            }

            var step = (relScore*grad.length)|0;
            if (step < 0) step = 0;
            if (step >= grad.length) step = grad.length - 1;
            fill = grad[step];
        }

        gg = new BoxGlyph(minPos, y, (maxPos - minPos), height, fill, stroke, alpha);
    } else if (gtype === 'HIDDEN') {
        gg = new PaddedGlyph(null, minPos, maxPos);
        noLabel = true;
    } else if (gtype === 'ARROW') {
        var color = style.FGCOLOR || 'purple';
        var parallel = isDasBooleanTrue(style.PARALLEL);
        var sw = isDasBooleanTrue(style.SOUTHWEST);
        var ne = isDasBooleanTrue(style.NORTHEAST);
        gg = new ArrowGlyph(minPos, maxPos, height, color, parallel, sw, ne);
    } else if (gtype === 'ANCHORED_ARROW') {
        var stroke = style.FGCOLOR || 'none';
        var fill = style.BGCOLOR || 'green';
        gg = new AArrowGlyph(minPos, maxPos, height, fill, stroke, strand);
        gg.bump = true;
    } else if (gtype === 'SPAN') {
        var stroke = style.FGCOLOR || 'black';
        gg = new SpanGlyph(minPos, maxPos, height, stroke);
    } else if (gtype === 'LINE') {
        var stroke = style.FGCOLOR || 'black';
        var lineStyle = style.STYLE || 'solid';
        gg = new LineGlyph(minPos, maxPos, height, lineStyle, strand, stroke);
    } else if (gtype === 'PRIMERS') {
        var stroke = style.FGCOLOR || 'black';
        var fill = style.BGCOLOR || 'red';
        gg = new PrimersGlyph(minPos, maxPos, height, fill, stroke);
    } else if (gtype === 'TEXT') {
        var string = style.STRING || 'text';
        var fill = style.FGCOLOR || 'black';
        gg = new TextGlyph(minPos, maxPos, height, fill, string);
    } else if (gtype === 'TOOMANY') {
        var stroke = style.FGCOLOR || 'gray';
        var fill = style.BGCOLOR || 'orange';
        gg = new TooManyGlyph(minPos, maxPos, height, fill, stroke);
    } else if (gtype === 'POINT') {
        var height = tier.forceHeight || style.HEIGHT || 30;
        var smin = tier.quantMin(style);
        var smax = tier.quantMax(style);
        var yscale = ((1.0 * height) / (smax - smin));
        var relScore = ((1.0 * score) - smin) / (smax-smin);
        var sc = ((score - (1.0*smin)) * yscale)|0;
        quant = {min: smin, max: smax};

        var fill = style.FGCOLOR || style.COLOR1 || 'black';
        if (style.COLOR2) {
            var grad = style._gradient;
            if (!grad) {
                grad = makeGradient(50, style.COLOR1, style.COLOR2, style.COLOR3);
                style._gradient = grad;
            }

            var step = (relScore*grad.length)|0;
            if (step < 0) step = 0;
            if (step >= grad.length) step = grad.length - 1;
            fill = grad[step];
        } 

        gg = new PointGlyph((minPos + maxPos)/2, height-sc, height, fill);
    } else if (gtype === '__SEQUENCE') {
        var rawseq = feature.seq;
        var seq = rawseq;
        var rawquals = feature.quals;
        var quals = rawquals;
        var insertionLabels = isDasBooleanTrue(style.__INSERTIONS);

        var indels = [];
        if (feature.cigar) {
            var ops = parseCigar(feature.cigar);
            seq = ''
            quals = '';
            var cursor = 0;
            for (var ci = 0; ci < ops.length; ++ci) {
                var co = ops[ci];
                if (co.op == 'M') {
                    seq += rawseq.substr(cursor, co.cnt);
                    quals += rawquals.substr(cursor, co.cnt);
                    cursor += co.cnt;
                } else if (co.op == 'D') {
                    for (var oi = 0; oi < co.cnt; ++oi) {
                        seq += '-';
                        quals += 'Z';
                    }
                } else if (co.op == 'I') {
                    var inseq =  rawseq.substr(cursor, co.cnt);
                    var ig = new TriangleGlyph(minPos + (seq.length*scale), 5, 'S', 5, 'red');
                    if (insertionLabels)
                        ig = new LabelledGlyph(ig, inseq, false, 'center', 'above', '7px sans-serif');
                    ig.feature = {label: 'Insertion: ' + inseq, type: 'insertion', method: 'insertion'};
                    indels.push(ig);

                    cursor += co.cnt;
                } else if (co.op == 'S') {
                    cursor += co.cnt;
                } else {
                    console.log('unknown cigop' + co.op);
                }
            }
        }

        var refSeq = null;
        if (tier.currentSequence) {
            var csStart = tier.currentSequence.start|0;
            var csEnd = tier.currentSequence.end|0;
            if (csStart <= max && csEnd >= min) {
                var sfMin = Math.max(min, csStart);
                var sfMax = Math.min(max, csEnd);

                refSeq = tier.currentSequence.seq.substr(sfMin - csStart, sfMax - sfMin + 1);
                while (min < sfMin) {
                    refSeq = 'N' + refSeq;
                    sfMin--;
                }
                while (max > sfMax) {
                    refSeq = refSeq + 'N';
                    sfMax++;
                }
            }
        }
        gg = new SequenceGlyph(minPos, maxPos, height, seq, refSeq, style.__SEQCOLOR, quals);
        if (insertionLabels)
            gg = new TranslatedGlyph(gg, 0, 7);
        if (indels.length > 0) {
            indels.splice(0, 0, gg);
            gg = new GroupGlyph(indels);
        }
    } else if (gtype === '__INSERTION') {
        var ig = new TriangleGlyph(minPos, 5, 'S', 5, 'red');
        gg = new LabelledGlyph(ig, feature.insertion || feature.altAlleles[0], false, 'center', 'above', '7px sans-serif');
        if ((maxPos - minPos) > 1) {
            var fill = style.BGCOLOR || style.COLOR1 || 'green';
            var bg = new BoxGlyph(minPos, 5, (maxPos - minPos), height, fill, stroke);
            gg = new GroupGlyph([bg, gg]);
        }
    } else if (gtype === '__NONE') {
        return null;
    } else /* default to BOX */ {
        var stroke = style.FGCOLOR || null;
        var fill = style.BGCOLOR || style.COLOR1 || 'green';
        if (style.BGITEM && feature.itemRgb)
            fill = feature.itemRgb;
        gg = new BoxGlyph(minPos, 0, (maxPos - minPos), height, fill, stroke);
        // gg.bump = true;
    }

    if ((isDasBooleanTrue(style.LABEL) || feature.forceLabel) && label && !noLabel) {
        gg = new LabelledGlyph(gg, label, false);
    }

    if (bump) {
        gg.bump = true;
    }

    gg.feature = feature;
    if (quant) {
        gg.quant = quant;
    }

    if (style.ZINDEX) {
        gg.zindex = style.ZINDEX | 0;
    }

    return gg;
}

DasTier.prototype.styleForFeature = function(f) {
    var ssScale = this.browser.zoomForCurrentScale();

    if (!this.stylesheet) {
        return null;
    }

    var maybe = null;
    var ss = this.stylesheet.styles;
    for (var si = 0; si < ss.length; ++si) {
        var sh = ss[si];
        if (sh.zoom && sh.zoom != ssScale) {
            continue;
        }

        if (sh.orientation) {
            if (sh.orientation != f.orientation) {
                continue;
            }
        }

        var labelRE = sh._labelRE;
        if (!labelRE || !labelRE.test) {
            labelRE = new RegExp('^' + sh.label + '$');
            sh._labelRE = labelRE;
        }
        if (sh.label && !(labelRE.test(f.label))) {
            continue;
        }
        var methodRE = sh._methodRE;
        if (!methodRE || !methodRE.test) {
            methodRE = new RegExp('^' + sh.method + '$');
            sh._methodRE = methodRE;
        }
        if (sh.method && !(methodRE.test(f.method))) {
            continue;
        }
        if (sh.type) {
            if (sh.type == 'default') {
                if (!maybe) {
                    maybe = sh.style;
                }
                continue;
            } else {
                var typeRE = sh._typeRE;
                if (!typeRE || !typeRE.test) {
                    typeRE = new RegExp('^' + sh.type + '$');
                    sh._typeRE = typeRE;
                }
                if (!typeRE.test(f.type)) 
                    continue;
            }
        }
        return sh.style;
    }
    return maybe;
}

function makeLineGlyph(features, style, tier) {
    var origin = tier.browser.viewStart, scale = tier.browser.scale;
    var height = tier.forceHeight || style.HEIGHT || 30;
    var min = tier.quantMin(style);
    var max = tier.quantMax(style);
    var yscale = ((1.0 * height) / (max - min));
    var width = style.LINEWIDTH || 1;
    var color = style.FGCOLOR || style.COLOR1 || 'black';

    var points = [];
    for (var fi = 0; fi < features.length; ++fi) {
        var f = features[fi];

        var px = ((((f.min|0) + (f.max|0)) / 2) - origin) * scale;
        var sc = ((f.score - (1.0*min)) * yscale)|0;
        var py = (height - sc);  // FIXME y???
        points.push(px);
        points.push(py);
    }
    var lgg = new LineGraphGlyph(points, color, height);
    lgg.quant = {min: min, max: max};

    if (style.ZINDEX) 
        lgg.zindex = style.ZINDEX|0;

    return lgg;
}

DasTier.prototype.quantMin = function(style) {
    if (this.forceMinDynamic) {
        return this.currentFeaturesMinScore || 0;
    } else if (typeof(this.forceMin) === 'number') {
        return this.forceMin;
    } else {
        return style.MIN || this.currentFeaturesMinScore || 0;
    }
}

DasTier.prototype.quantMax = function(style) {
    if (this.forceMaxDynamic) {
        return this.currentFeaturesMaxScore || 0;
    } else if (typeof(this.forceMax) === 'number') {
        return this.forceMax;
    } else {
        return style.MAX || this.currentFeaturesMaxScore || 0;
    }
}// 
// Dalliance Genome Explorer
// (c) Thomas Down 2006-2012
//
// sequence-draw.js: renderers for sequence-related data
//

"use strict";

var MIN_TILE = 100;
var rulerTileColors = ['black', 'white'];
var baseColors = {A: 'green', C: 'blue', G: 'black', T: 'red'};
var steps = [1,2,5];


function tileSizeForScale(scale, min)
{
    if (!min) {
        min = MIN_TILE;
    }

    function ts(p) {
        return steps[p % steps.length] * Math.pow(10, (p / steps.length)|0);
    }
    var pow = steps.length;
    while (scale * ts(pow) < min) {
        ++pow;
    }
    return ts(pow);
}

function drawSeqTier(tier, seq)
{
	var retina = tier.browser.retina && window.devicePixelRatio > 1;
    var scale = tier.browser.scale, knownStart = tier.browser.viewStart - (1000/scale)|0, knownEnd = tier.browser.viewEnd + (2000/scale), currentSeqMax = tier.browser.currentSeqMax;

	var desiredWidth = tier.browser.featurePanelWidth + 2000;
    if (retina) {
        desiredWidth *= 2;
    }
    var fpw = tier.viewport.width|0; // this.browser.featurePanelWidth;
    if (fpw < desiredWidth - 50) {
        tier.viewport.width = fpw = desiredWidth;
    }

    var height = 50;
    if (seq && seq.seq) {
		height += 25;
    }

    var canvasHeight = height;
    if (retina) 
    	canvasHeight *= 2;

    tier.viewport.height = canvasHeight;
    tier.viewport.style.height = '' + height + 'px';
    tier.viewport.style.width = retina ? ('' + (fpw/2) + 'px') : ('' + fpw + 'px');
    tier.layoutHeight = height;
    tier.updateHeight();

    var gc = tier.viewport.getContext('2d');
    if (tier.background) {
        gc.fillStyle = tier.background;
        gc.fillRect(0, 0, fpw, tier.viewport.height);
    }
    if (retina) {
        gc.scale(2, 2);
    }
    gc.translate(1000,0);

    var seqTierMax = knownEnd;
    if (currentSeqMax > 0 && currentSeqMax < knownEnd) {
        seqTierMax = currentSeqMax;
    }
    var tile = tileSizeForScale(scale);
    var pos = Math.max(0, ((knownStart / tile)|0) * tile);
    
    var origin = tier.browser.viewStart;

    while (pos <= seqTierMax) {
		gc.fillStyle = ((pos / tile) % 2 == 0) ? 'white' : 'black';
		gc.strokeStyle = 'black';
		gc.fillRect((pos - origin) * scale,
			    8,
			    tile*scale,
			    3);
		gc.strokeRect((pos - origin) * scale,
			      8,
			      tile*scale,
			      3);

		gc.fillStyle = 'black';
		gc.fillText(formatLongInt(pos), ((pos - origin) * scale), 22);
		

		pos += tile;
    }

    if (seq && seq.seq) {
		for (var p = knownStart; p <= knownEnd; ++p) {
		    if (p >= seq.start && p <= seq.end) {
				var base = seq.seq.substr(p - seq.start, 1).toUpperCase();
				var color = baseColors[base];
				if (!color) {
		            color = 'gray';
				}

				gc.fillStyle = color;

				if (scale >= 8) {
				    gc.fillText(base, (p - origin) * scale, 52);
				} else {
				    gc.fillRect((p - origin) * scale, 42, scale, 12); 
				}
		    }
		}
    } 

    tier.norigin = tier.browser.viewStart;
    tier.viewport.style.left = '-1000px';
}

function svgSeqTier(tier, seq) {
    var scale = tier.browser.scale, knownStart = tier.browser.viewStart - (1000/scale)|0, knownEnd = tier.browser.viewEnd + (2000/scale), currentSeqMax = tier.browser.currentSeqMax;

    var fpw = tier.viewport.width|0; 

    var seqTierMax = knownEnd;
    if (currentSeqMax > 0 && currentSeqMax < knownEnd) {
        seqTierMax = currentSeqMax;
    }
    var tile = tileSizeForScale(scale);
    var pos = Math.max(0, ((knownStart / tile)|0) * tile);
    
    var origin = tier.browser.viewStart;

    var  g = makeElementNS(NS_SVG, 'g', [], {fontSize: '8pt'}); 
    while (pos <= seqTierMax) {
	g.appendChild(
	    makeElementNS(
		NS_SVG, 'rect',
		null,
		{x: (pos-origin)*scale,
		 y: 8,
		 width: tile*scale,
		 height: 3,
		 fill: ((pos / tile) % 2 == 0) ? 'white' : 'black',
		 stroke: 'black'}));

	g.appendChild(
	    makeElementNS(
		NS_SVG, 'text',
		formatLongInt(pos),
		{x: (pos-origin)*scale,
		 y: 28,
		 fill: 'black', stroke: 'none'}));
	
	pos += tile;
    }

    if (seq && seq.seq) {
	for (var p = knownStart; p <= knownEnd; ++p) {
	    if (p >= seq.start && p <= seq.end) {
		var base = seq.seq.substr(p - seq.start, 1).toUpperCase();
		var color = baseColors[base];
		if (!color) {
                    color = 'gray';
		}


		if (scale >= 8) {
		    // gc.fillText(base, (p - origin) * scale, 12);
		    g.appendChild(
			makeElementNS(NS_SVG, 'text', base, {
			    x: (p-origin)*scale,
			    y: 52,
			    fill: color}));
		} else {
		    g.appendChild(
			makeElementNS(NS_SVG, 'rect', null, {
			    x: (p - origin)*scale,
			    y: 42,
			    width: scale,
			    height: 12,
	                    fill: color}));

		}
	    }
	}
    } 

    return g;
}
"use strict";

function sortFeatures(tier)
{
    var dmin = tier.browser.drawnStart, dmax = tier.browser.drawnEnd;
    var ungroupedFeatures = {};
    var groupedFeatures = {};
    var drawnGroupedFeatures = {};
    var groupMins = {}, groupMaxes = {};
    var groups = {};
    var superGroups = {};
    var groupsToSupers = {};
    var nonPositional = [];
    var minScore, maxScore;
    var fbid;

    var init_fbid = function() {
        fbid = {};
        for (var fi = 0; fi < tier.currentFeatures.length; ++fi) {
            var f = tier.currentFeatures[fi];
            if (f.id) {
                fbid[f.id] = f;
            }
        }
    };
    
    var superParentsOf = function(f) {
        // FIXME: should recur.
        var spids = [];
        if (f.parents) {
            for (var pi = 0; pi < f.parents.length; ++pi) {
                var pid = f.parents[pi];
                var p = fbid[pid];
                if (!p) {
                    continue;
                }
                // alert(p.type + ':' + p.typeCv);
                if (p.typeCv == 'SO:0000704') {
                    pushnew(spids, pid);
                }
            }
        }
        return spids;
    }

    for (var fi = 0; fi < tier.currentFeatures.length; ++fi) {
        var f = tier.currentFeatures[fi];
        if (f.parts) {
            continue;
        }

        var drawn = f.min <= dmax && f.max >= dmin;

        if (!f.min || !f.max) {
            nonPositional.push(f);
            continue;
        }

        if (f.score && f.score != '.' && f.score != '-') {
            var sc = 1.0 * f.score;
            if (!minScore || sc < minScore) {
                minScore = sc;
            }
            if (!maxScore || sc > maxScore) {
                maxScore = sc;
            }
        }

        var fGroups = [];
        var fSuperGroup = null;
        if (f.groups) {
            for (var gi = 0; gi < f.groups.length; ++gi) {
                var g = f.groups[gi];
                var gid = g.id;
                if (g.type == 'gene') {
                    // Like a super-grouper...
                    fSuperGroup = gid; 
                    groups[gid] = g;
                } else if (g.type == 'translation') {
                    // have to ignore this to get sensible results from bj-e :-(.
                } else {
                    pusho(groupedFeatures, gid, f);
                    groups[gid] = g;
                    fGroups.push(gid);

                    var ogm = groupMins[gid];
                    if (!ogm || f.min < ogm)
                        groupMins[gid] = f.min;

                    ogm = groupMaxes[gid];
                    if (!ogm || f.max > ogm)
                        groupMaxes[gid] = f.max;
                }
            }
        }

        if (f.parents) {
            if (!fbid) {
                init_fbid();
            }
            for (var pi = 0; pi < f.parents.length; ++pi) {
                var pid = f.parents[pi];
                var p = fbid[pid];
                if (!p) {
                    // alert("couldn't find " + pid);
                    continue;
                }
                if (!p.parts) {
                    p.parts = [f];
                }
                pushnewo(groupedFeatures, pid, p);
                pusho(groupedFeatures, pid, f);
                
                if (!groups[pid]) {
                    groups[pid] = {
                        type: p.type,
                        id: p.id,
                        label: p.label || p.id
                    };
                }
                fGroups.push(pid);

                var ogm = groupMins[pid];
                if (!ogm || f.min < ogm)
                    groupMins[pid] = f.min;

                ogm = groupMaxes[pid];
                if (!ogm || f.max > ogm)
                    groupMaxes[pid] = f.max;

                var sgs = superParentsOf(p);
                if (sgs.length > 0) {
                    fSuperGroup = sgs[0];
                    var sp = fbid[sgs[0]];
                    groups[sgs[0]] = {
                        type: sp.type,
                        id: sp.id,
                        label: sp.label || sp.id
                    };
                    if (!tier.dasSource.collapseSuperGroups) {
                        tier.dasSource.collapseSuperGroups = true;
                    }
                }
            }   
        }

        if (fGroups.length == 0) {
            if (drawn)
                pusho(ungroupedFeatures, f.type, f);
        } else if (fSuperGroup) {
            for (var g = 0; g < fGroups.length; ++g) {
                var gid = fGroups[g];
                pushnewo(superGroups, fSuperGroup, gid);
                groupsToSupers[gid] = fSuperGroup;
            } 
        }       
    }

    for (var gid in groupedFeatures) {
        var group = groups[gid];
        if (typeof(group.min) !== 'number') 
            group.min = groupMins[gid];
        if (typeof(group.max) !== 'number') 
            group.max = groupMaxes[gid];

        if (groupMaxes[gid] >= dmin && groupMins[gid] <= dmax)
            drawnGroupedFeatures[gid] = groupedFeatures[gid];
    }

    tier.ungroupedFeatures = ungroupedFeatures;
    tier.groupedFeatures = drawnGroupedFeatures;
    tier.groups = groups;
    tier.superGroups = superGroups;
    tier.groupsToSupers = groupsToSupers;

    if (minScore) {
        if (minScore > 0) {
            minScore = 0;
        } else if (maxScore < 0) {
            maxScore = 0;
        }
        tier.currentFeaturesMinScore = minScore;
        tier.currentFeaturesMaxScore = maxScore;
    }
}
/* -*- mode: javascript; c-basic-offset: 4; indent-tabs-mode: nil -*- */

// 
// Dalliance Genome Explorer
// (c) Thomas Down 2006-2011
//
// feature-popup.js
//

"use strict";

var TAGVAL_NOTE_RE = new RegExp('^([A-Za-z]+)=(.+)');

Browser.prototype.addFeatureInfoPlugin = function(handler) {
    if (!this.featureInfoPlugins) {
        this.featureInfoPlugins = [];
    }
    this.featureInfoPlugins.push(handler);
}

function FeatureInfo(hit, feature, group) {
    var name = pick(group.type, feature.type);
    var fid = pick(group.label, feature.label, group.id, feature.id);
    if (fid && fid.indexOf('__dazzle') != 0) {
        name = name + ': ' + fid;
    }

    this.hit = hit;
    this.feature = feature;
    this.group = group;
    this.title = name;
    this.sections = [];
}

FeatureInfo.prototype.setTitle = function(t) {
    this.title = t;
}

FeatureInfo.prototype.add = function(label, info) {
    if (typeof info === 'string') {
        info = makeElement('span', info);
    }
    this.sections.push({label: label, info: info});
}

Browser.prototype.featurePopup = function(ev, __ignored_feature, hit, tier) {
    var hi = hit.length;
    var feature = --hi >= 0 ? hit[hi] : {};
    var group = --hi >= 0 ? hit[hi] : {};

    var featureInfo = new FeatureInfo(hit, feature, group);
    featureInfo.tier = tier;
    var fips = this.featureInfoPlugins || [];
    for (var fipi = 0; fipi < fips.length; ++fipi) {
        try {
            fips[fipi](feature, featureInfo);
        } catch (e) {
            console.log(e.stack || e);
        }
    }
    fips = tier.featureInfoPlugins || [];
    for (fipi = 0; fipi < fips.length; ++fipi) {
        try {
            fips[fipi](feature, featureInfo);
        } catch (e) {
            console.log(e.stack || e);
        }
    }

    this.removeAllPopups();

    var table = makeElement('table', null, {className: 'table table-striped table-condensed'});
    table.style.width = '100%';
    table.style.margin = '0px';

    var idx = 0;
    if (feature.method) {
        var row = makeElement('tr', [
            makeElement('th', 'Method'),
            makeElement('td', feature.method)
        ]);
        table.appendChild(row);
        ++idx;
    }
    {
        var loc;
        if (group.segment) {
            loc = group;
        } else {
            loc = feature;
        }
        var row = makeElement('tr', [
            makeElement('th', 'Location'),
            makeElement('td', loc.segment + ':' + loc.min + '-' + loc.max, {}, {minWidth: '200px'})
        ]);
        table.appendChild(row);
        ++idx;
    }
    if (feature.score !== undefined && feature.score !== null && feature.score != '-') {
        var row = makeElement('tr', [
            makeElement('th', 'Score'),
            makeElement('td', '' + feature.score)
        ]);
        table.appendChild(row);
        ++idx;
    }
    {
        var links = maybeConcat(group.links, feature.links);
        if (links && links.length > 0) {
            var row = makeElement('tr', [
                makeElement('th', 'Links'),
                makeElement('td', links.map(function(l) {
                    return makeElement('div', makeElement('a', l.desc, {href: l.uri, target: '_new'}));
                }))
            ]);
            table.appendChild(row);
            ++idx;
        }
    }
    {
        var notes = maybeConcat(group.notes, feature.notes);
        for (var ni = 0; ni < notes.length; ++ni) {
            var k = 'Note';
            var v = notes[ni];
            var m = v.match(TAGVAL_NOTE_RE);
            if (m) {
                k = m[1];
                v = m[2];
            }

            var row = makeElement('tr', [
                makeElement('th', k),
                makeElement('td', v)
            ]);
            table.appendChild(row);
            ++idx;
        }
    }

    for (var fisi = 0; fisi < featureInfo.sections.length; ++fisi) {
        var section = featureInfo.sections[fisi];
        table.appendChild(makeElement('tr', [
            makeElement('th', section.label),
            makeElement('td', section.info)]));
    }        

    this.popit(ev, featureInfo.title || 'Feature', table, {width: 450});
}
/* -*- mode: javascript; c-basic-offset: 4; indent-tabs-mode: nil -*- */

// 
// Dalliance Genome Explorer
// (c) Thomas Down 2006-2010
//
// karyoscape.js
//

"use strict";

function Karyoscape(browser, dsn)
{
    this.browser = browser; // for tooltips.
    this.dsn = dsn;
    this.svg = makeElementNS(NS_SVG, 'g');
    this.width = 250;
}

Karyoscape.prototype.update = function(chr, start, end) {
    this.start = start;
    this.end = end;
    if (!this.chr || chr != this.chr) {
        this.chr = chr;
        removeChildren(this.svg);

        var kscape = this;
        this.dsn.features(
            new DASSegment(chr),
            {type: 'karyotype'},
            function(karyos, err, segmentMap) {
                if (segmentMap && segmentMap[chr] && segmentMap[chr].max) {
                    kscape.chrLen = segmentMap[chr].max;
                } else {
                    kscape.chrLen = null;
                }
                kscape.karyos = karyos || [];
                kscape.redraw();
            }
        );
    } else {
        this.setThumb();
    }
}

var karyo_palette = {
    gneg: 'white',
    gpos25: 'rgb(200,200,200)',
    gpos33: 'rgb(180,180,180)',
    gpos50: 'rgb(128,128,128)',
    gpos66: 'rgb(100,100,100)',
    gpos75: 'rgb(64,64,64)',
    gpos100: 'rgb(0,0,0)',
    gpos: 'rgb(0,0,0)',
    gvar: 'rgb(100,100,100)',
    acen: 'rgb(100,100,100)',
    stalk: 'rgb(100,100,100)'
};

Karyoscape.prototype.redraw = function() {
    removeChildren(this.svg);
    this.karyos = this.karyos.sort(function(k1, k2) {
        return (k1.min|0) - (k2.min|0);
    });
    if (this.karyos.length > 0) {
        if (!this.chrLen) {
            this.chrLen = this.karyos[this.karyos.length - 1].max;
        }
    } else {
        if (!this.chrLen) {
            alert('Warning: insufficient data to set up spatial navigator');
            this.chrLen = 200000000;
        } 
        this.karyos.push({
            min: 1,
            max: this.chrLen,
            label: 'gneg'
        });
    }
    var bandspans = null;
    for (var i = 0; i < this.karyos.length; ++i) {
        var k = this.karyos[i];
        var bmin = ((1.0 * k.min) / this.chrLen) * this.width;
        var bmax = ((1.0 * k.max) / this.chrLen) * this.width;
        var col = karyo_palette[k.label];
        if (!col) {
            // alert("don't understand " + k.label);
        } else {
            if (bmax > bmin) {
                var band = makeElementNS(NS_SVG, 'rect', null, {
                    x: bmin,
                    y: (k.label == 'stalk' || k.label == 'acen' ? 5 : 0),
                    width: (bmax - bmin),
                    height: (k.label == 'stalk' || k.label == 'acen'? 5 : 15),
                    stroke: 'none',
                    fill: col
                });
                if (k.label.substring(0, 1) == 'g') {
                    var br = new Range(k.min, k.max);
                    if (bandspans == null) {
                        bandspans = br;
                    } else {
                        bandspans = union(bandspans, br);
                    }
                }
                this.browser.makeTooltip(band, k.id);
                this.svg.appendChild(band);
            }
        }
    }

    if (bandspans) {
        var r = bandspans.ranges();

        var pathopsT = 'M 0 10 L 0 0';
        var pathopsB = 'M 0 5 L 0 15';
        
        var curx = 0;
        for (var ri = 0; ri < r.length; ++ri) {
            var rr = r[ri];
            var bmin = ((1.0 * rr.min()) / this.chrLen) * this.width;
            var bmax = ((1.0 * rr.max()) / this.chrLen) * this.width;
            if ((bmin - curx > 0.75)) {
                pathopsT += ' M ' + bmin + ' 0';
                pathopsB += ' M ' + bmin + ' 15';
            }
            pathopsT +=  ' L ' + bmax + ' 0';
            pathopsB +=  ' L ' + bmax + ' 15';
            curx = bmax;
        }
        if ((this.width - curx) > 0.75) {
            pathopsT += ' M ' + this.width + ' 0';
            pathopsB += ' M ' + this.width + ' 15';
        } else {
            pathopsT += ' L ' + this.width + ' 0';
            pathopsB += ' L ' + this.width + ' 15';
        }
        pathopsT +=  ' L ' + this.width + ' 10';
        pathopsB +=  ' L ' + this.width + ' 5';
        this.svg.appendChild(makeElementNS(NS_SVG, 'path', null, {
            d: pathopsT + ' ' + pathopsB,
            stroke: 'black',
            strokeWidth: 2,
            fill: 'none'
        }));
    }

    this.thumb = makeElementNS(NS_SVG, 'rect', null, {
        x: 50, y: -5, width: 8, height: 25,
        fill: 'blue', fillOpacity: 0.5, stroke: 'none'
    });
    this.svg.appendChild(this.thumb);
    this.setThumb();

    var thisKaryo = this;
    var sliderDeltaX;

    var moveHandler = function(ev) {
        ev.stopPropagation(); ev.preventDefault();
        var sliderX = Math.max(-4, Math.min(ev.clientX + sliderDeltaX, thisKaryo.width - 4));
        thisKaryo.thumb.setAttribute('x', sliderX);
//      if (thisSlider.onchange) {
//          thisSlider.onchange(value, false);
//      }
    }
    var upHandler = function(ev) {
        ev.stopPropagation(); ev.preventDefault();
        if (thisKaryo.onchange) {
            thisKaryo.onchange((1.0 * ((thisKaryo.thumb.getAttribute('x')|0) + 4)) / thisKaryo.width, true);
        }
        document.removeEventListener('mousemove', moveHandler, true);
        document.removeEventListener('mouseup', upHandler, true);
    }

    this.thumb.addEventListener('mousedown', function(ev) {
        ev.stopPropagation(); ev.preventDefault();
        sliderDeltaX = thisKaryo.thumb.getAttribute('x') - ev.clientX;
        document.addEventListener('mousemove', moveHandler, true);
        document.addEventListener('mouseup', upHandler, true);
    }, false);
}

Karyoscape.prototype.setThumb = function() {
    var pos = ((this.start|0) + (this.end|0)) / 2
    var gpos = ((1.0 * pos)/this.chrLen) * this.width;
    if (this.thumb) {
        this.thumb.setAttribute('x', gpos - 4);
    }
}
            

/* -*- mode: javascript; c-basic-offset: 4; indent-tabs-mode: nil -*- */

// 
// Dalliance Genome Explorer
// (c) Thomas Down 2006-2013
//
// kspace.js
//

"use strict";

function FetchPool() {
    this.reqs = [];
    this.awaitedFeatures = {};
}

FetchPool.prototype.addRequest = function(xhr) {
    this.reqs.push(xhr);
}

FetchPool.prototype.abortAll = function() {
    for (var i = 0; i < this.reqs.length; ++i) {
        this.reqs[i].abort();
    }
}

function KSCacheBaton(chr, min, max, scale, features, status, coverage) {
    this.chr = chr;
    this.min = min;
    this.max = max;
    this.coverage = coverage;
    this.scale = scale;
    this.features = features || [];
    this.status = status;
}

KSCacheBaton.prototype.toString = function() {
    return this.chr + ":" + this.min + ".." + this.max + ";scale=" + this.scale;
}

function KnownSpace(tierMap, chr, min, max, scale, seqSource) {
    this.tierMap = tierMap;
    this.chr = chr;
    this.min = min;
    this.max = max;
    this.scale = scale;
    this.seqSource = seqSource || new DummySequenceSource();
    this.viewCount = 0;

    this.featureCache = {};
    this.latestViews = {};
}

KnownSpace.prototype.bestCacheOverlapping = function(chr, min, max) {
    var baton = this.featureCache[this.tierMap[0]];
    if (baton) {
        return baton;
    } else {
        return null;
    }
}

KnownSpace.prototype.viewFeatures = function(chr, min, max, scale) {
    // dlog('viewFeatures(' + chr + ', ' + min + ', ' + max + ', ' + scale +')');
    if (scale != scale) {
        throw "viewFeatures called with silly scale";
    }

    if (chr != this.chr) {
        throw "Can't extend Known Space to a new chromosome";
    }
    if (min < 1) {
        min = 1;
    }

    this.min = min;
    this.max = max;
    this.scale = scale;

    if (this.pool) {
        this.pool.abortAll();
    }
    this.pool = new FetchPool();
    this.awaitedSeq = new Awaited();
    this.seqWasFetched = false;
    this.viewCount++;
    
    this.startFetchesForTiers(this.tierMap);
}
    
function filterFeatures(features, min, max) {
    var ff = [];
    var featuresByGroup = {};

    for (var fi = 0; fi < features.length; ++fi) {
        var f = features[fi];
        if (!f.min || !f.max) {
            ff.push(f);
        } else if (f.groups && f.groups.length > 0) {
            pusho(featuresByGroup, f.groups[0].id, f);
        } else if (f.min <= max && f.max >= min) {
            ff.push(f);
        }
    }

    for (var gid in featuresByGroup) {
        var gf = featuresByGroup[gid];
        var gmin = 100000000000, gmax = -100000000000;
        for (var fi = 0; fi < gf.length; ++fi) {
            var f = gf[fi];
            gmin = Math.min(gmin, f.min);
            gmax = Math.max(gmax, f.max);
        }
        if (gmin <= max || gmax >= min) {
            for (var fi = 0; fi < gf.length; ++fi) {
                ff.push(gf[fi]);
            }
        }
    }

    return ff;
}

KnownSpace.prototype.invalidate = function(tier) {
    if (!this.pool) {
        return;
    }

    this.featureCache[tier] = null;
    this.startFetchesForTiers([tier]);
}

KnownSpace.prototype.startFetchesForTiers = function(tiers) {
    var thisB = this;

    var awaitedSeq = this.awaitedSeq;
    var needSeq = false;

    var gex;

    for (var t = 0; t < tiers.length; ++t) {
        try {
            if (this.startFetchesFor(tiers[t], awaitedSeq)) {
                needSeq = true;
            }
        } catch (ex) {
            tiers[t].updateStatus(ex);
            console.log('Error fetching tier source');
            console.log(ex);
            gex = ex;
        }
    }

    if (needSeq && !this.seqWasFetched) {
        this.seqWasFetched = true;
        // dlog('needSeq ' + this.chr + ':' + this.min + '..' + this.max);
        var smin = this.min, smax = this.max;

        if (this.cs) {
            if (this.cs.start <= smin && this.cs.end >= smax) {
                var cachedSeq;
                if (this.cs.start == smin && this.cs.end == smax) {
                    cachedSeq = this.cs;
                } else {
                    cachedSeq = new DASSequence(this.cs.name, smin, smax, this.cs.alphabet, 
                                                this.cs.seq.substring(smin - this.cs.start, smax + 1 - this.cs.start));
                }
                return awaitedSeq.provide(cachedSeq);
            }
        }
        
        this.seqSource.fetch(this.chr, smin, smax, this.pool, function(err, seq) {
            if (seq) {
                if (!thisB.cs || (smin <= thisB.cs.start && smax >= thisB.cs.end) || 
                    (smin >= thisB.cs.end) || (smax <= thisB.cs.start) || 
                    ((smax - smin) > (thisB.cs.end - thisB.cs.start))) 
                {
                    thisB.cs = seq;
                }
                awaitedSeq.provide(seq);
            } else {
                dlog('Noseq: ' + miniJSONify(err));
                awaitedSeq.provide(null);
            }
        });
    } 

    if (gex)
        throw gex;
}

KnownSpace.prototype.startFetchesFor = function(tier, awaitedSeq) {
    var thisB = this;

    var viewID = this.viewCount;
    var source = tier.getSource() || new DummyFeatureSource();
    var needsSeq = tier.needsSequence(this.scale);
    var baton = thisB.featureCache[tier];
    var wantedTypes = tier.getDesiredTypes(this.scale);
    var chr = this.chr, min = this.min, max = this.max;


    if (wantedTypes === undefined) {
//         dlog('skipping because wantedTypes is undef');
        return false;
    }
    if (baton) {
//      dlog('considering cached features: ' + baton);
    }
    if (baton && baton.chr === this.chr && baton.min <= min && baton.max >= max) {
        var cachedFeatures = baton.features;
        if (baton.min < min || baton.max > max) {
            cachedFeatures = filterFeatures(cachedFeatures, min, max);
        }
        
        // dlog('cached scale=' + baton.scale + '; wanted scale=' + thisB.scale);
//      if ((baton.scale < (thisB.scale/2) && cachedFeatures.length > 200) || (wantedTypes && wantedTypes.length == 1 && wantedTypes.indexOf('density') >= 0) ) {
//          cachedFeatures = downsample(cachedFeatures, thisB.scale);
//      }
//      console.log('Provisioning ' + tier.toString() + ' with ' + cachedFeatures.length + ' features from cache (' + baton.min + ', ' + baton.max + ')');
//      tier.viewFeatures(baton.chr, Math.max(baton.min, this.min), Math.min(baton.max, this.max), baton.scale, cachedFeatures);   // FIXME change scale if downsampling

        thisB.provision(tier, baton.chr, intersection(baton.coverage, new Range(min, max)), baton.scale, wantedTypes, cachedFeatures, baton.status, needsSeq ? awaitedSeq : null);

        var availableScales = source.getScales();
        if (baton.scale <= this.scale || !availableScales) {
//          dlog('used cached features');
            return needsSeq;
        } else {
//          dlog('used cached features (temporarily)');
        }
    }

    if (source.instrument)
        console.log('Starting  fetch ' + viewID + ' (' + min + ', ' + max + ')');
    source.fetch(chr, min, max, this.scale, wantedTypes, this.pool, function(status, features, scale, coverage) {
    	if (source.instrument)
    	    console.log('Finishing fetch ' + viewID);

    	var latestViewID = thisB.latestViews[tier] || -1;
    	if (latestViewID > viewID) {
    	    return;
    	}

        if (!coverage) {
            coverage = new Range(min, max);
        }

        if (!baton || (min < baton.min) || (max > baton.max)) {         // FIXME should be merging in some cases?
            thisB.featureCache[tier] = new KSCacheBaton(chr, min, max, scale, features, status, coverage);
        }

	    thisB.latestViews[tier] = viewID;
        thisB.provision(tier, chr, coverage, scale, wantedTypes, features, status, needsSeq ? awaitedSeq : null);
    });
    return needsSeq;
}

KnownSpace.prototype.provision = function(tier, chr, coverage, actualScale, wantedTypes, features, status, awaitedSeq) {
    tier.updateStatus(status);
   
    if (!status) {
        var mayDownsample = false;
        var src = tier.getSource();
        while (MappedFeatureSource.prototype.isPrototypeOf(src) || CachingFeatureSource.prototype.isPrototypeOf(src) || OverlayFeatureSource.prototype.isPrototypeOf(src)) {
	        if (OverlayFeatureSource.prototype.isPrototypeOf(src)) {
		        src = src.sources[0];
	        } else {
		        src = src.source;
	        }
        }
        if (BWGFeatureSource.prototype.isPrototypeOf(src) || RemoteBWGFeatureSource.prototype.isPrototypeOf(src) || BAMFeatureSource.prototype.isPrototypeOf(src) || RemoteBAMFeatureSource.prototype.isPrototypeOf(src)) {
            mayDownsample = true;
        }

    	if (!src.opts || (!src.opts.forceReduction && !src.opts.noDownsample)) {
            if ((actualScale < (this.scale/2) && features.length > 200)  ||
		        (mayDownsample && wantedTypes && wantedTypes.length == 1 && wantedTypes.indexOf('density') >= 0))
            {
		        features = downsample(features, this.scale);
            }
    	}

        if (awaitedSeq) {
            awaitedSeq.await(function(seq) {
                tier.viewFeatures(chr, coverage, actualScale, features, seq);
            });
        } else {
            tier.viewFeatures(chr, coverage, actualScale, features);
        }
    }
}
/* -*- mode: javascript; c-basic-offset: 4; indent-tabs-mode: nil -*- */

// 
// Dalliance Genome Explorer
// (c) Thomas Down 2006-2010
//
// sample.js: downsampling of quantitative features
//

"use strict";

var __DS_SCALES = [1, 2, 5];

function ds_scale(n) {
    return __DS_SCALES[n % __DS_SCALES.length] * Math.pow(10, (n / __DS_SCALES.length)|0);
}


function DSBin(scale, min, max) {
    this.scale = scale;
    this.tot = 0;
    this.cnt = 0;
    this.hasScore = false;
    this.min = min; this.max = max;
    this.features = [];
}

function _featureOrder(a, b) {
    if (a.min < b.min) {
        return -1;
    } else if (a.min > b.min) {
        return 1;
    } else if (a.max < b.max) {
        return -1;
    } else if (b.max > a.max) {
        return 1;
    } else {
        return 0;
    }
}

DSBin.prototype.score = function() {
    if (this.cnt == 0) {
        return 0;
    } else if (this.hasScore) {
        return this.tot / this.cnt;
    } else {
        var features = this.features;
        features.sort(_featureOrder);

        var maxSeen = -10000000000;
        var cov=0, lap=0;

        for (var fi = 1; fi < features.length; ++fi) {
            var f = features[fi];
            var lMin = Math.max(f.min, this.min);
            var lMax = Math.min(f.max, this.max);
            lap += (lMax - lMin + 1);

            if (lMin > maxSeen) {
                cov += lMax - lMin + 1;
                maxSeen = lMax;
            } else {
                if (lMax > maxSeen) {
                    cov += (lMax - maxSeen);
                    maxSeen = lMax;
                }
            }
        }

        if (cov > 0)
            return (1.0 * lap) / cov;
        else
            return 0;
    }
}

DSBin.prototype.feature = function(f) {
    if (f.score !== undefined) {
        this.tot += f.score;
        this.hasScore = true
    }

    ++this.cnt;
    this.features.push(f);
}

function downsample(features, targetRez) {
    var sn = 0;
    while (ds_scale(sn + 1) < targetRez) {
        ++sn;
    }
    var scale = ds_scale(sn);

    var binTots = [];
    var maxBin = -10000000000;
    var minBin = 10000000000;
    for (var fi = 0; fi < features.length; ++fi) {
        var f = features[fi];
        if (f.groups && f.groups.length > 0) {
            // Don't downsample complex features (?)
            return features;
        }

        var minLap = (f.min / scale)|0;
        var maxLap = (f.max / scale)|0;
        maxBin = Math.max(maxBin, maxLap);
        minBin = Math.min(minBin, minLap);
        for (var b = minLap; b <= maxLap; ++b) {
            var bm = binTots[b];
            if (!bm) {
                bm = new DSBin(scale, b * scale, (b + 1) * scale - 1);
                binTots[b] = bm;
            }
            bm.feature(f);
        }
    }

    var sampledFeatures = [];
    for (var b = minBin; b <= maxBin; ++b) {
        var bm = binTots[b];
        if (bm) {
            var f = new DASFeature();
            f.segment = features[0].segment;
            f.min = (b * scale) + 1;
            f.max = (b + 1) * scale;
            f.score = bm.score();
            f.type = 'density';
            sampledFeatures.push(f);
        }
    }

    var afterDS = Date.now();
    return sampledFeatures;
}
/*
 * A JavaScript implementation of the Secure Hash Algorithm, SHA-1, as defined
 * in FIPS 180-1
 * Version 2.2 Copyright Paul Johnston 2000 - 2009.
 * Other contributors: Greg Holt, Andrew Kepert, Ydnar, Lostinet
 * Distributed under the BSD License
 * See http://pajhome.org.uk/crypt/md5 for details.
 */

 "use strict";

/*
 * Configurable variables. You may need to tweak these to be compatible with
 * the server-side, but the defaults work in most cases.
 */
var hexcase = 0;  /* hex output format. 0 - lowercase; 1 - uppercase        */
var b64pad  = ""; /* base-64 pad character. "=" for strict RFC compliance   */

/*
 * These are the functions you'll usually want to call
 * They take string arguments and return either hex or base-64 encoded strings
 */
function hex_sha1(s)    { return rstr2hex(rstr_sha1(str2rstr_utf8(s))); }
function b64_sha1(s)    { return rstr2b64(rstr_sha1(str2rstr_utf8(s))); }
function any_sha1(s, e) { return rstr2any(rstr_sha1(str2rstr_utf8(s)), e); }
function hex_hmac_sha1(k, d)
  { return rstr2hex(rstr_hmac_sha1(str2rstr_utf8(k), str2rstr_utf8(d))); }
function b64_hmac_sha1(k, d)
  { return rstr2b64(rstr_hmac_sha1(str2rstr_utf8(k), str2rstr_utf8(d))); }
function any_hmac_sha1(k, d, e)
  { return rstr2any(rstr_hmac_sha1(str2rstr_utf8(k), str2rstr_utf8(d)), e); }

/*
 * Perform a simple self-test to see if the VM is working
 */
function sha1_vm_test()
{
  return hex_sha1("abc").toLowerCase() == "a9993e364706816aba3e25717850c26c9cd0d89d";
}

/*
 * Calculate the SHA1 of a raw string
 */
function rstr_sha1(s)
{
  return binb2rstr(binb_sha1(rstr2binb(s), s.length * 8));
}

/*
 * Calculate the HMAC-SHA1 of a key and some data (raw strings)
 */
function rstr_hmac_sha1(key, data)
{
  var bkey = rstr2binb(key);
  if(bkey.length > 16) bkey = binb_sha1(bkey, key.length * 8);

  var ipad = Array(16), opad = Array(16);
  for(var i = 0; i < 16; i++)
  {
    ipad[i] = bkey[i] ^ 0x36363636;
    opad[i] = bkey[i] ^ 0x5C5C5C5C;
  }

  var hash = binb_sha1(ipad.concat(rstr2binb(data)), 512 + data.length * 8);
  return binb2rstr(binb_sha1(opad.concat(hash), 512 + 160));
}

/*
 * Convert a raw string to a hex string
 */
function rstr2hex(input)
{
  try { hexcase } catch(e) { hexcase=0; }
  var hex_tab = hexcase ? "0123456789ABCDEF" : "0123456789abcdef";
  var output = "";
  var x;
  for(var i = 0; i < input.length; i++)
  {
    x = input.charCodeAt(i);
    output += hex_tab.charAt((x >>> 4) & 0x0F)
           +  hex_tab.charAt( x        & 0x0F);
  }
  return output;
}

/*
 * Convert a raw string to a base-64 string
 */
function rstr2b64(input)
{
  try { b64pad } catch(e) { b64pad=''; }
  var tab = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
  var output = "";
  var len = input.length;
  for(var i = 0; i < len; i += 3)
  {
    var triplet = (input.charCodeAt(i) << 16)
                | (i + 1 < len ? input.charCodeAt(i+1) << 8 : 0)
                | (i + 2 < len ? input.charCodeAt(i+2)      : 0);
    for(var j = 0; j < 4; j++)
    {
      if(i * 8 + j * 6 > input.length * 8) output += b64pad;
      else output += tab.charAt((triplet >>> 6*(3-j)) & 0x3F);
    }
  }
  return output;
}

/*
 * Convert a raw string to an arbitrary string encoding
 */
function rstr2any(input, encoding)
{
  var divisor = encoding.length;
  var remainders = Array();
  var i, q, x, quotient;

  /* Convert to an array of 16-bit big-endian values, forming the dividend */
  var dividend = Array(Math.ceil(input.length / 2));
  for(i = 0; i < dividend.length; i++)
  {
    dividend[i] = (input.charCodeAt(i * 2) << 8) | input.charCodeAt(i * 2 + 1);
  }

  /*
   * Repeatedly perform a long division. The binary array forms the dividend,
   * the length of the encoding is the divisor. Once computed, the quotient
   * forms the dividend for the next step. We stop when the dividend is zero.
   * All remainders are stored for later use.
   */
  while(dividend.length > 0)
  {
    quotient = Array();
    x = 0;
    for(i = 0; i < dividend.length; i++)
    {
      x = (x << 16) + dividend[i];
      q = Math.floor(x / divisor);
      x -= q * divisor;
      if(quotient.length > 0 || q > 0)
        quotient[quotient.length] = q;
    }
    remainders[remainders.length] = x;
    dividend = quotient;
  }

  /* Convert the remainders to the output string */
  var output = "";
  for(i = remainders.length - 1; i >= 0; i--)
    output += encoding.charAt(remainders[i]);

  /* Append leading zero equivalents */
  var full_length = Math.ceil(input.length * 8 /
                                    (Math.log(encoding.length) / Math.log(2)))
  for(i = output.length; i < full_length; i++)
    output = encoding[0] + output;

  return output;
}

/*
 * Encode a string as utf-8.
 * For efficiency, this assumes the input is valid utf-16.
 */
function str2rstr_utf8(input)
{
  var output = "";
  var i = -1;
  var x, y;

  while(++i < input.length)
  {
    /* Decode utf-16 surrogate pairs */
    x = input.charCodeAt(i);
    y = i + 1 < input.length ? input.charCodeAt(i + 1) : 0;
    if(0xD800 <= x && x <= 0xDBFF && 0xDC00 <= y && y <= 0xDFFF)
    {
      x = 0x10000 + ((x & 0x03FF) << 10) + (y & 0x03FF);
      i++;
    }

    /* Encode output as utf-8 */
    if(x <= 0x7F)
      output += String.fromCharCode(x);
    else if(x <= 0x7FF)
      output += String.fromCharCode(0xC0 | ((x >>> 6 ) & 0x1F),
                                    0x80 | ( x         & 0x3F));
    else if(x <= 0xFFFF)
      output += String.fromCharCode(0xE0 | ((x >>> 12) & 0x0F),
                                    0x80 | ((x >>> 6 ) & 0x3F),
                                    0x80 | ( x         & 0x3F));
    else if(x <= 0x1FFFFF)
      output += String.fromCharCode(0xF0 | ((x >>> 18) & 0x07),
                                    0x80 | ((x >>> 12) & 0x3F),
                                    0x80 | ((x >>> 6 ) & 0x3F),
                                    0x80 | ( x         & 0x3F));
  }
  return output;
}

/*
 * Encode a string as utf-16
 */
function str2rstr_utf16le(input)
{
  var output = "";
  for(var i = 0; i < input.length; i++)
    output += String.fromCharCode( input.charCodeAt(i)        & 0xFF,
                                  (input.charCodeAt(i) >>> 8) & 0xFF);
  return output;
}

function str2rstr_utf16be(input)
{
  var output = "";
  for(var i = 0; i < input.length; i++)
    output += String.fromCharCode((input.charCodeAt(i) >>> 8) & 0xFF,
                                   input.charCodeAt(i)        & 0xFF);
  return output;
}

/*
 * Convert a raw string to an array of big-endian words
 * Characters >255 have their high-byte silently ignored.
 */
function rstr2binb(input)
{
  var output = Array(input.length >> 2);
  for(var i = 0; i < output.length; i++)
    output[i] = 0;
  for(var i = 0; i < input.length * 8; i += 8)
    output[i>>5] |= (input.charCodeAt(i / 8) & 0xFF) << (24 - i % 32);
  return output;
}

/*
 * Convert an array of big-endian words to a string
 */
function binb2rstr(input)
{
  var output = "";
  for(var i = 0; i < input.length * 32; i += 8)
    output += String.fromCharCode((input[i>>5] >>> (24 - i % 32)) & 0xFF);
  return output;
}

/*
 * Calculate the SHA-1 of an array of big-endian words, and a bit length
 */
function binb_sha1(x, len)
{
  /* append padding */
  x[len >> 5] |= 0x80 << (24 - len % 32);
  x[((len + 64 >> 9) << 4) + 15] = len;

  var w = Array(80);
  var a =  1732584193;
  var b = -271733879;
  var c = -1732584194;
  var d =  271733878;
  var e = -1009589776;

  for(var i = 0; i < x.length; i += 16)
  {
    var olda = a;
    var oldb = b;
    var oldc = c;
    var oldd = d;
    var olde = e;

    for(var j = 0; j < 80; j++)
    {
      if(j < 16) w[j] = x[i + j];
      else w[j] = bit_rol(w[j-3] ^ w[j-8] ^ w[j-14] ^ w[j-16], 1);
      var t = safe_add(safe_add(bit_rol(a, 5), sha1_ft(j, b, c, d)),
                       safe_add(safe_add(e, w[j]), sha1_kt(j)));
      e = d;
      d = c;
      c = bit_rol(b, 30);
      b = a;
      a = t;
    }

    a = safe_add(a, olda);
    b = safe_add(b, oldb);
    c = safe_add(c, oldc);
    d = safe_add(d, oldd);
    e = safe_add(e, olde);
  }
  return Array(a, b, c, d, e);

}

/*
 * Perform the appropriate triplet combination function for the current
 * iteration
 */
function sha1_ft(t, b, c, d)
{
  if(t < 20) return (b & c) | ((~b) & d);
  if(t < 40) return b ^ c ^ d;
  if(t < 60) return (b & c) | (b & d) | (c & d);
  return b ^ c ^ d;
}

/*
 * Determine the appropriate additive constant for the current iteration
 */
function sha1_kt(t)
{
  return (t < 20) ?  1518500249 : (t < 40) ?  1859775393 :
         (t < 60) ? -1894007588 : -899497514;
}

/*
 * Add integers, wrapping at 2^32. This uses 16-bit operations internally
 * to work around bugs in some JS interpreters.
 */
function safe_add(x, y)
{
  var lsw = (x & 0xFFFF) + (y & 0xFFFF);
  var msw = (x >> 16) + (y >> 16) + (lsw >> 16);
  return (msw << 16) | (lsw & 0xFFFF);
}

/*
 * Bitwise rotate a 32-bit number to the left.
 */
function bit_rol(num, cnt)
{
  return (num << cnt) | (num >>> (32 - cnt));
}
/* -*- mode: javascript; c-basic-offset: 4; indent-tabs-mode: nil -*- */

// 
// Dalliance Genome Explorer
// (c) Thomas Down 2006-2013
//
// svg-export.js
//

Browser.prototype.makeSVG = function(opts) {
    opts = opts || {};
    var minTierHeight = opts.minTierHeight || 20;
    var padding = 3;

    var b = this;
    var saveDoc = document.implementation.createDocument(NS_SVG, 'svg', null);

    var saveRoot = makeElementNS(NS_SVG, 'g', null, {
        fontFamily: 'helvetica',
        fontSize: '8pt'
    });
    saveDoc.documentElement.appendChild(saveRoot);

    var margin = 200;

    var dallianceAnchor = makeElementNS(NS_SVG, 'a',
       makeElementNS(NS_SVG, 'text', 'Graphics from Dalliance ' + VERSION, {
           x: (b.featurePanelWidth + margin + 20)/2,
           y: 30,
           strokeWidth: 0,
           fontSize: '12pt',
	       textAnchor: 'middle',
	       fill: 'blue'
       }));
    dallianceAnchor.setAttribute('xmlns:xlink', NS_XLINK);
    dallianceAnchor.setAttribute('xlink:href', 'http://www.biodalliance.org/');
  
    saveRoot.appendChild(dallianceAnchor);
    
    var clipRect = makeElementNS(NS_SVG, 'rect', null, {
    	x: margin,
    	y: 50,
    	width: b.featurePanelWidth,
    	height: 100000
    });
    var clip = makeElementNS(NS_SVG, 'clipPath', clipRect, {id: 'featureClip'});
    saveRoot.appendChild(clip);

    var pos = 70;
    var tierHolder = makeElementNS(NS_SVG, 'g', null, {});

    for (var ti = 0; ti < b.tiers.length; ++ti) {
        var tier = b.tiers[ti];
    	var tierSVG = makeElementNS(NS_SVG, 'g', null, {clipPath: 'url(#featureClip)', clipRule: 'nonzero'});
    	var tierLabels = makeElementNS(NS_SVG, 'g');
    	var tierTopPos = pos;

    	var tierBackground = makeElementNS(NS_SVG, 'rect', null, {x: 0, y: tierTopPos, width: '10000', height: 50, fill: tier.background});
    	tierSVG.appendChild(tierBackground);

    	if (tier.dasSource.tier_type === 'sequence') {
    	    var seqTrack = svgSeqTier(tier, tier.currentSequence);
    	    
    	    tierSVG.appendChild(makeElementNS(NS_SVG, 'g', seqTrack, {transform: 'translate(' + (margin) + ', ' + pos + ')'}));
    	    pos += 80;
    	} else {
            if (!tier.subtiers) {
    		   continue;
            }
    	
    	    var offset = ((tier.glyphCacheOrigin - b.viewStart) * b.scale);
            var hasQuant = false;
            for (var sti = 0; sti < tier.subtiers.length; ++sti) {
                pos += padding;
        		var subtier = tier.subtiers[sti];
                    
        		var glyphElements = [];
        		for (var gi = 0; gi < subtier.glyphs.length; ++gi) {
                    var glyph = subtier.glyphs[gi];
                    glyphElements.push(glyph.toSVG());
        		}

    		    tierSVG.appendChild(makeElementNS(NS_SVG, 'g', glyphElements, {transform: 'translate(' + (margin+offset) + ', ' + pos + ')'}));

        		if (subtier.quant) {
                    hasQuant = true;
        		    var q = subtier.quant;
                    var h = subtier.height;

                    var numTics = 2;
                    if (h > 40) {
                        numTics = 1 + ((h/20) | 0);
                    }
                    var ticSpacing = h / (numTics - 1);
                    var ticInterval = (q.max - q.min) / (numTics - 1);

        		    var path = new SVGPath();
        		    path.moveTo(margin + 5, pos);
        		    path.lineTo(margin, pos);
        		    path.lineTo(margin, pos + subtier.height);
        		    path.lineTo(margin + 5, pos + subtier.height);
                    for (var t = 1; t < numTics-1; ++t) {
                        var ty = t*ticSpacing;
                        path.moveTo(margin, pos + ty);
                        path.lineTo(margin+3, pos + ty);
                    }

        		    tierLabels.appendChild(makeElementNS(NS_SVG, 'path', null, {d: path.toPathData(), fill: 'none', stroke: 'black', strokeWidth: '2px'}));
        		    tierLabels.appendChild(makeElementNS(NS_SVG, 'text', formatQuantLabel(q.max), {x: margin - 3, y: pos + 7, textAnchor: 'end'}));
        		    tierLabels.appendChild(makeElementNS(NS_SVG, 'text', formatQuantLabel(q.min), {x: margin - 3, y: pos +  subtier.height, textAnchor: 'end'}));
                    for (var t = 1; t < numTics-1; ++t) {
                        var ty = t*ticSpacing;
                        tierLabels.appendChild(makeElementNS(NS_SVG, 'text', formatQuantLabel((1.0*q.max) - (t*ticInterval)), 
                            {x: margin - 3, y: pos +  ty + 3, textAnchor: 'end'}));
                    }
        		}

    		    pos += subtier.height + padding;
            }

            if (pos - tierTopPos < minTierHeight) {
                pos = tierTopPos + minTierHeight;
            }
    	}

        var labelName;
        if (typeof tier.config.name === 'string')
            labelName = tier.config.name;
        else
            labelName = tier.dasSource.name;
    	tierLabels.appendChild(
    	    makeElementNS(
    		NS_SVG, 'text',
    		labelName,
    		{x: margin - (hasQuant ? 20 : 12), y: (pos+tierTopPos+8)/2, fontSize: '10pt', textAnchor: 'end'}));

    	
    	tierBackground.setAttribute('height', pos - tierTopPos);
    	tierHolder.appendChild(makeElementNS(NS_SVG, 'g', [tierSVG, tierLabels]));
    }

    if (opts.highlights) {
        var highlights = this.highlights || [];
        for (var hi = 0; hi < highlights.length; ++hi) {
            var h = highlights[hi];
            if ((h.chr == this.chr || h.chr == ('chr' + this.chr)) && h.min < this.viewEnd && h.max > this.viewStart) {
                var tmin = (Math.max(h.min, this.viewStart) - this.viewStart) * this.scale;
                var tmax = (Math.min(h.max, this.viewEnd) - this.viewStart) * this.scale;

                tierHolder.appendChild(makeElementNS(NS_SVG, 'rect', null, {x: margin + tmin, y: 70, width: (tmax-tmin), height: pos-70,
                                                                      stroke: 'none', fill: this.defaultHighlightFill, fillOpacity: this.defaultHighlightAlpha}));
            }
        }
    }

    var rulerPos = -1; 
    if (opts.ruler == 'center') {
        rulerPos = margin + ((this.viewEnd - this.viewStart + 1)*this.scale) / 2;
    } else if (opts.ruler == 'left') {
        rulerPos = margin;
    } else if (opts.ruler == 'right') {
        rulerPos = margin + ((this.viewEnd - this.viewStart + 1)*this.scale);
    }
    if (rulerPos >= 0) {
        tierHolder.appendChild(makeElementNS(NS_SVG, 'line', null, {x1: rulerPos, y1: 70, x2: rulerPos, y2: pos,
                                                              stroke: 'blue'}));
    }

    saveRoot.appendChild(tierHolder);
    saveDoc.documentElement.setAttribute('width', b.featurePanelWidth + 20 + margin);
    saveDoc.documentElement.setAttribute('height', pos + 50);

    var svgBlob = new Blob([new XMLSerializer().serializeToString(saveDoc)], {type: 'image/svg+xml'});
    return svgBlob;
}
/* -*- mode: javascript; c-basic-offset: 4; indent-tabs-mode: nil -*- */

// 
// Dalliance Genome Explorer
// (c) Thomas Down 2006-2010
//
// spans.js: JavaScript Intset/Location port.
//

"use strict";


function Range(min, max)
{
    if (typeof(min) != 'number' || typeof(max) != 'number')
        throw 'Bad range ' + min + ',' + max;
    this._min = min;
    this._max = max;
}

Range.prototype.min = function() {
    return this._min;
}

Range.prototype.max = function() {
    return this._max;
}

Range.prototype.contains = function(pos) {
    return pos >= this._min && pos <= this._max;
}

Range.prototype.isContiguous = function() {
    return true;
}

Range.prototype.ranges = function() {
    return [this];
}

Range.prototype._pushRanges = function(ranges) {
    ranges.push(this);
}

Range.prototype.toString = function() {
    return '[' + this._min + '-' + this._max + ']';
}

function _Compound(ranges) {
    this._ranges = ranges;
    // assert sorted?
}

_Compound.prototype.min = function() {
    return this._ranges[0].min();
}

_Compound.prototype.max = function() {
    return this._ranges[this._ranges.length - 1].max();
}

_Compound.prototype.contains = function(pos) {
    // FIXME implement bsearch if we use this much.
    for (var s = 0; s < this._ranges.length; ++s) {
        if (this._ranges[s].contains(pos)) {
            return true;
        }
    }
    return false;
}

_Compound.prototype.isContiguous = function() {
    return this._ranges.length > 1;
}

_Compound.prototype.ranges = function() {
    return this._ranges;
}

_Compound.prototype._pushRanges = function(ranges) {
    for (var ri = 0; ri < this._ranges.length; ++ri)
        ranges.push(this._ranges[ri]);
}

_Compound.prototype.toString = function() {
    var s = '';
    for (var r = 0; r < this._ranges.length; ++r) {
        if (r>0) {
            s = s + ',';
        }
        s = s + this._ranges[r].toString();
    }
    return s;
}

function union(s0, s1) {
    if (! (s0 instanceof Array)) {
        s0 = [s0];
        if (s1)
            s0.push(s1);
    }

    if (s0.length == 0)
        return null;
    else if (s0.length == 1)
        return s0[0];

    var ranges = [];
    for (var si = 0; si < s0.length; ++si)
        s0[si]._pushRanges(ranges);
    ranges = ranges.sort(_rangeOrder);

    var oranges = [];
    var current = ranges[0];
    current = new Range(current._min, current._max);  // Copy now so we don't have to later.

    for (var i = 1; i < ranges.length; ++i) {
        var nxt = ranges[i];
        if (nxt._min > (current._max + 1)) {
            oranges.push(current);
            current = new Range(nxt._min, nxt._max);
        } else {
            if (nxt._max > current._max) {
                current._max = nxt._max;
            }
        }
    }
    oranges.push(current);

    if (oranges.length == 1) {
        return oranges[0];
    } else {
        return new _Compound(oranges);
    }
}

function intersection(s0, s1) {
    var r0 = s0.ranges();
    var r1 = s1.ranges();
    var l0 = r0.length, l1 = r1.length;
    var i0 = 0, i1 = 0;
    var or = [];

    while (i0 < l0 && i1 < l1) {
        var s0 = r0[i0], s1 = r1[i1];
        var lapMin = Math.max(s0.min(), s1.min());
        var lapMax = Math.min(s0.max(), s1.max());
        if (lapMax >= lapMin) {
            or.push(new Range(lapMin, lapMax));
        }
        if (s0.max() > s1.max()) {
            ++i1;
        } else {
            ++i0;
        }
    }
    
    if (or.length == 0) {
        return null; // FIXME
    } else if (or.length == 1) {
        return or[0];
    } else {
        return new _Compound(or);
    }
}

function coverage(s) {
    var tot = 0;
    var rl = s.ranges();
    for (var ri = 0; ri < rl.length; ++ri) {
        var r = rl[ri];
        tot += (r.max() - r.min() + 1);
    }
    return tot;
}



function rangeOrder(a, b)
{
    if (a.min() < b.min()) {
        return -1;
    } else if (a.min() > b.min()) {
        return 1;
    } else if (a.max() < b.max()) {
        return -1;
    } else if (b.max() > a.max()) {
        return 1;
    } else {
        return 0;
    }
}

function _rangeOrder(a, b)
{
    if (a._min < b._min) {
        return -1;
    } else if (a._min > b._min) {
        return 1;
    } else if (a._max < b._max) {
        return -1;
    } else if (b._max > a._max) {
        return 1;
    } else {
        return 0;
    }
}
/* -*- mode: javascript; c-basic-offset: 4; indent-tabs-mode: nil -*- */

// 
// Dalliance Genome Explorer
// (c) Thomas Down 2006-2013
//
// thub.js: support for track-hub style registries
//

"use strict";

var THUB_STANZA_REGEXP = /\n\s*\n/;
var THUB_PARSE_REGEXP  = /(\w+) +(.+)\n?/;
var THUB_SUBGROUP_REGEXP = /subGroup[1-9]/;

var THUB_PENNANT_PREFIX = 'http://genome.ucsc.edu/images/';

function TrackHub(url) {
    this.genomes = {};
    this.url = url;
}

function TrackHubTrack() {
}

TrackHubTrack.prototype.get = function(k) {
    if (this[k])
        return this[k];
    else if (this._parent) 
        return this._parent.get(k);
}

function TrackHubDB(hub) {
    this.hub = hub;
}

TrackHubDB.prototype.getTracks = function(callback) {
    var thisB = this;
    if (this._tracks) {
        return callback(this._tracks);
    } 
    
    textXHR(this.absURL, function(trackFile, err) {
        if (err) {
            return callback(null, err);
        }
        
        // trackFile = trackFile.replace(/\#.*/g, '');
        trackFile = trackFile.replace('\\\n', ' ');

        var tracks = [];
        var tracksById = {};
        var stanzas = trackFile.split(THUB_STANZA_REGEXP);
        for (var s = 0; s < stanzas.length; ++s) {
            var toks = stanzas[s].replace(/\#.*/g, '').split(THUB_PARSE_REGEXP);
            var track = new TrackHubTrack();
            track._db = thisB;
            for (var l = 0; l < toks.length - 2; l += 3) {
                var k = toks[l+1], v = toks[l+2];
                if (k.match(THUB_SUBGROUP_REGEXP)) {
                    if (!track.subgroups)
                        track.subgroups = {};
                    var sgtoks = v.split(/\s/);
                    var sgtag = sgtoks[0];
                    var sgrecord = {name: sgtoks[1], tags: [], titles: []};
                    for (var sgti = 2; sgti < sgtoks.length; ++sgti) {
                        var grp = sgtoks[sgti].split(/=/);
                        sgrecord.tags.push(grp[0]);
                        sgrecord.titles.push(grp[1]);
                    }
                    track.subgroups[sgtag] = sgrecord;
                } else if (k === 'subGroups') {
                    var sgtoks = v.split(/(\w+)=(\w+)/);
                    track.sgm = {};
                    for (var sgti = 0; sgti < sgtoks.length - 2; sgti += 3) {
                        track.sgm[sgtoks[sgti+1]] = sgtoks[sgti + 2];
                    }
                } else {
                    track[toks[l+1]] = toks[l+2];
                }
            }

            if (track.track && (track.type || track.container || track.view)) {
                tracks.push(track);
                tracksById[track.track] = track;
            }
        }
        
        var toplevels = [];
        var composites = [];
        for (var ti = 0; ti < tracks.length; ++ti) {
            var track = tracks[ti];
            var top = true;
            if (track.parent) {
                var ptoks = track.parent.split(/\s+/);
                var parent = tracksById[ptoks[0]];
                if (parent) {
                    track._parent = parent;

                    if (!parent.children)
                        parent.children = [];
                    parent.children.push(track);

                    if (parent)
                        top = false;
                } else {
                    // console.log("Couldn't find parent " + ptoks[0] + '(' + track.parent + ')');
                }
               
            }
            if (track.compositeTrack) {
                composites.push(track);
            } else if (top) {
                toplevels.push(track);
            }
        }

        for (var ci = 0; ci < composites.length; ++ci) {
            var comp = composites[ci];
            if (!comp.children)
                continue;

            var parentOfViews = false;
            for (var ki = 0; ki < comp.children.length; ++ki) {
                var k = comp.children[ki];
                if (k.view) {
                    k.shortLabel = comp.shortLabel + ": " + k.shortLabel;
                    toplevels.push(k);
                    parentOfViews = true;
                }
            }
            if (!parentOfViews)
                toplevels.push(comp);
        }
            
        thisB._tracks = toplevels;
        return callback(thisB._tracks, null);
    }, {credentials: this.credentials, salt: true});
}

function connectTrackHub(hubURL, callback, opts) {
    opts = opts || {};
    opts.salt = true;

    textXHR(hubURL, function(hubFile, err) {
        if (err) {
            return callback(null, err);
        }

        var toks = hubFile.split(THUB_PARSE_REGEXP);
        var hub = new TrackHub(hubURL);
        if (opts.credentials) {
            hub.credentials = opts.credentials;
        }
        for (var l = 0; l < toks.length - 2; l += 3) {
            hub[toks[l+1]] = toks[l+2];
        }
        
        
        if (hub.genomesFile) {
            var genURL = relativeURL(hubURL, hub.genomesFile);
            textXHR(genURL, function(genFile, err) {
                if (err) {
                    return callback(null, err);
                }

                var stanzas = genFile.split(THUB_STANZA_REGEXP);
                for (var s = 0; s < stanzas.length; ++s) {
                    var toks = stanzas[s].split(THUB_PARSE_REGEXP);
                    var gprops = new TrackHubDB(hub);
                    if (opts.credentials) {
                        gprops.credentials = opts.credentials;
                    }

                    for (var l = 0; l < toks.length - 2; l += 3) {
                        gprops[toks[l+1]] = toks[l+2];
                    }

                    if (gprops.twoBitPath) {
                        gprops.twoBitPath = relativeURL(genURL, gprops.twoBitPath);
                    }

                    if (gprops.genome && gprops.trackDb) {
                        gprops.absURL = relativeURL(genURL, gprops.trackDb);
                        hub.genomes[gprops.genome] = gprops;
                    }
                }

                callback(hub);
                        
            }, opts);
        } else {
            callback(null, 'No genomesFile');
        }
    }, opts);
}


TrackHubTrack.prototype.toDallianceSource = function() {
    var source = {
        name: this.shortLabel,
        desc: this.longLabel
    };
    if (this._db.mapping) {
        source.mapping = this._db.mapping;
    }

    var pennantIcon = this.get('pennantIcon');
    if (pennantIcon) {
        var ptoks = pennantIcon.split(/\s+/);
        source.pennant = THUB_PENNANT_PREFIX + ptoks[0];
    }

    var searchTrix = this.get('searchTrix');
    if (searchTrix) {
        source.trixURI = relativeURL(this._db.absURL, searchTrix);
    }

    if (this.container == 'multiWig') {
        source.merge = 'concat';
        source.overlay = [];
        var children = this.children || [];
        source.style = [];
        source.noDownsample = true;

        for (var ci = 0; ci < children.length; ++ci) {
            var ch = children[ci];
            var cs = ch.toDallianceSource()
            source.overlay.push(cs);

            if (cs.style) {
                for (var si = 0; si < cs.style.length; ++si) {
                    var style = cs.style[si];
                    style.method = ch.shortLabel;  // FIXME
                    if (this.aggregate == 'transparentOverlay')
                        style.style.ALPHA = 0.5;
                    source.style.push(style);
                }
            }
        }
        return source;

        
    } else {
        var typeToks = this.type.split(/\s+/);
        if (typeToks[0] == 'bigBed') {
            var bedTokens = typeToks[1]|0
            var bedPlus = typeToks[2] == '+';

            source.bwgURI = relativeURL(this._db.absURL, this.bigDataUrl);
            source.style = this.bigbedStyles();
            if (this._db.credentials) {
                source.credentials = true;
            }
            if (bedTokens >= 12 && bedPlus)
                source.collapseSuperGroups = true;
            return source;
        } else if (typeToks[0] == 'bigWig') {
            source.bwgURI = relativeURL(this._db.absURL, this.bigDataUrl);
            source.style = this.bigwigStyles();
            source.noDownsample = true;     // FIXME seems like a blunt instrument...
            
            if (this.yLineOnOff && this.yLineOnOff == 'on') {
                source.quantLeapThreshold = this.yLineMark !== undefined ? (1.0 * this.yLineMark) : 0.0;
            }

            if (this._db.credentials) {
                source.credentials = true;
            }

            return source;
        } else if (typeToks[0] == 'bam') {
            source.bamURI = relativeURL(this._db.absURL, this.bigDataUrl);
            if (this._db.credentials) {
                source.credentials = true;
            }
            return source;
        } else {
            console.log('Unsupported ' + this.type);
        }
    }
}

TrackHubTrack.prototype.bigwigStyles = function() {
    var typeToks = this.type.split(/\s+/);
    var min, max;
    if (typeToks.length >= 3) {
        min = 1.0 * typeToks[1];
        max = 1.0 * typeToks[2];
    }

    var height;
    if (this.maxHeightPixels) {
        var mhpToks = this.maxHeightPixels.split(/:/);
        if (mhpToks.length == 3) {
            height = mhpToks[1] | 0;
        } else {
            console.log('maxHeightPixels should be of the form max:default:min');
        }
    }
    
    var gtype = 'bars';
    if (this.graphTypeDefault) {
        gtype = this.graphTypeDefault;
    }
    
    var color = 'black';
    var altColor = null;
    if (this.color) {
        color = 'rgb(' + this.color + ')';
    }
    if (this.altColor) {
        altColor = 'rgb(' + this.altColor + ')';
    }
    
    var stylesheet = new DASStylesheet();
    var wigStyle = new DASStyle();
    if (gtype == 'points') {
        wigStyle.glyph = 'POINT';
    } else {
        wigStyle.glyph = 'HISTOGRAM';
    }

    if (altColor) {
        wigStyle.COLOR1 = color;
        wigStyle.COLOR2 = altColor;
    } else {
        wigStyle.BGCOLOR = color;
    }
    wigStyle.HEIGHT = height || 30;
    if (min || max) {
        wigStyle.MIN = min;
        wigStyle.MAX = max;
    }
    stylesheet.pushStyle({type: 'default'}, null, wigStyle);
    return stylesheet.styles;
}

TrackHubTrack.prototype.bigbedStyles = function() {
    var itemRgb = (''+this.get('itemRgb')).toLowerCase() == 'on';
    var visibility = this.get('visibility') || 'full';
    var color = this.get('color');
    if (color)
        color = 'rgb(' + color + ')';
    else 
        color = 'blue';
    
    var stylesheet = new DASStylesheet();
    var wigStyle = new DASStyle();
    wigStyle.glyph = 'BOX';
    wigStyle.FGCOLOR = 'black';
    wigStyle.BGCOLOR = color;
    wigStyle.HEIGHT = (visibility == 'full' || visibility == 'pack') ? 12 : 8;
    wigStyle.BUMP = (visibility == 'full' || visibility == 'pack');
    wigStyle.LABEL = (visibility == 'full' || visibility == 'pack');
    wigStyle.ZINDEX = 20;
    if (itemRgb)
        wigStyle.BGITEM = true;

    var cbs = this.get('colorByStrand');
    if (cbs) {
        var cbsToks = cbs.split(/\s+/);
        
        var plus = shallowCopy(wigStyle);
        plus.BGCOLOR = 'rgb(' + cbsToks[0] + ')';
        stylesheet.pushStyle({type: 'bigwig', orientation: '+'}, null, plus);

        var minus = shallowCopy(wigStyle);
        minus.BGCOLOR = 'rgb(' + cbsToks[1] + ')';
        stylesheet.pushStyle({type: 'bigwig', orientation: '-'}, null, minus);
    } else {
        stylesheet.pushStyle({type: 'bigwig'}, null, wigStyle);
    }   
    
    var tlStyle = new DASStyle();
    tlStyle.glyph = 'BOX';
    tlStyle.FGCOLOR = 'black';
    if (itemRgb)
        tlStyle.BGITEM = true;
    tlStyle.BGCOLOR = 'red'
    tlStyle.HEIGHT = 10;
    tlStyle.BUMP = true;
    tlStyle.ZINDEX = 20;
    stylesheet.pushStyle({type: 'bb-translation'}, null, tlStyle);
    
    var tsStyle = new DASStyle();
    tsStyle.glyph = 'BOX';
    tsStyle.FGCOLOR = 'black';
    tsStyle.BGCOLOR = 'white';
    tsStyle.HEIGHT = 10;
    tsStyle.ZINDEX = 10;
    tsStyle.BUMP = true;
    tsStyle.LABEL = true;
    stylesheet.pushStyle({type: 'bb-transcript'}, null, tsStyle);

/*
    var densStyle = new DASStyle();
    densStyle.glyph = 'HISTOGRAM';
    densStyle.COLOR1 = 'white';
    densStyle.COLOR2 = 'black';
    densStyle.HEIGHT=30;
    stylesheet.pushStyle({type: 'density'}, null, densStyle); */

    return stylesheet.styles;
}

function THUB_COMPARE(g, h) {
    if (g.priority && h.priority) {
        return (1.0 * g.priority) - (1.0 * h.priority)
    } else if (g.priority) {
        return 1;
    } else if (h.priority) {
        return -1;
    } else {
        return g.shortLabel.localeCompare(h.shortLabel);
    }
}
/* -*- mode: javascript; c-basic-offset: 4; indent-tabs-mode: nil -*- */

// 
// Dalliance Genome Explorer
// (c) Thomas Down 2006-2010
//
// tier.js: (try) to encapsulate the functionality of a browser tier.
//

"use strict";

var __tier_idSeed = 0;

function DasTier(browser, source, config, background)
{
    this.config = config || {};
    this.id = 'tier' + (++__tier_idSeed);
    this.browser = browser;
    this.dasSource = shallowCopy(source);
    this.background = background;

    this.viewport = makeElement('canvas', null, 
                                {width: '' + ((this.browser.featurePanelWidth|0) + 2000), 
                                 height: "30",
                                 className: 'viewport'});
    this.overlay = makeElement('canvas', null,
         {width: + ((this.browser.featurePanelWidth|0) + 2000), 
          height: "30",
          className: 'viewport-overlay'});

    this.notifier = makeElement('div', '', {className: 'notifier'});
    this.notifierHolder = makeElement('div', this.notifier, {className: 'notifier-holder'});
    this.quantOverlay = makeElement(
        'canvas', null, 
        {width: '50', height: "56",
         className: 'quant-overlay'});

    this.removeButton = makeElement('i', null, {className: 'fa fa-times'});
    this.bumpButton = makeElement('i', null, {className: 'fa fa-plus-circle'});
    this.loaderButton = browser.makeLoader(16);
    this.loaderButton.style.display = 'none';
    this.infoElement = makeElement('div', this.dasSource.desc, {className: 'track-label-info'});
    this.nameButton = makeElement('div', [], {className: 'tier-tab'});
    this.nameButton.appendChild(this.removeButton);
    if (source.pennant) {
        this.nameButton.appendChild(makeElement('img', null, {src: source.pennant, width: '16', height: '16'}))
    }
    this.nameElement = makeElement('span', source.name);
    this.nameButton.appendChild(makeElement('span', [this.nameElement, this.infoElement], {className: 'track-name-holder'}));
    this.nameButton.appendChild(this.bumpButton);
    this.nameButton.appendChild(this.loaderButton);
    
    this.label = makeElement('span',
       [this.nameButton],
       {className: 'btn-group track-label'});


    this.row = makeElement('div', [this.viewport,
                                   this.overlay, 
                                   this.quantOverlay],
                            {className: 'tier'});

    if (!background) {
        this.row.style.background = 'none';
    }

    if (!browser.noDefaultLabels)
        this.row.appendChild(this.label);
    this.row.appendChild(this.notifierHolder);
    
    this.layoutHeight = 25;
    this.bumped = true;
    this.styleIdSeed = 0;
    if (source.quantLeapThreshold) {
        this.quantLeapThreshold = source.quantLeapThreshold;
    }
    if (this.dasSource.collapseSuperGroups) {
        this.bumped = false;
    }
    this.layoutWasDone = false;

    if (source.featureInfoPlugin) {
        this.addFeatureInfoPlugin(source.featureInfoPlugin);
    }

    this.initSources();

    var thisB = this;
    if (this.featureSource && this.featureSource.getDefaultFIPs) {
        this.featureSource.getDefaultFIPs(function(fip) {
            if (fip)
                thisB.addFeatureInfoPlugin(fip);
        });
    }

    if (this.featureSource && this.featureSource.addReadinessListener) {
        this.featureSource.addReadinessListener(function(ready) {
            thisB.notify(ready, -1);
        });
    }

    this.listeners = [];
    this.featuresLoadedListeners = [];
}

DasTier.prototype.toString = function() {
    return this.id;
}

DasTier.prototype.addFeatureInfoPlugin = function(p) {
    if (!this.featureInfoPlugins) 
        this.featureInfoPlugins = [];
    this.featureInfoPlugins.push(p);
}

DasTier.prototype.init = function() {
    var tier = this;

    if (tier.dasSource.style) {
        this.setStylesheet({styles: tier.dasSource.style});
        this.browser.refreshTier(this);
    } else {
        var ssSource;
        if (tier.dasSource.stylesheet_uri) {
            ssSource = new DASFeatureSource(tier.dasSource);
        } else {
            ssSource = tier.getSource();
        }
        tier.status = 'Fetching stylesheet';
        
        ssSource.getStyleSheet(function(ss, err) {
            if (err || !ss) {
                tier.error = 'No stylesheet';
                var ss = new DASStylesheet();
                var defStyle = new DASStyle();
                defStyle.glyph = 'BOX';
                defStyle.BGCOLOR = 'blue';
                defStyle.FGCOLOR = 'black';
                ss.pushStyle({type: 'default'}, null, defStyle);
                tier.setStylesheet(ss);
                tier.browser.refreshTier(tier);
            } else {
                tier.setStylesheet(ss);
                if (ss.geneHint) {
                    tier.dasSource.collapseSuperGroups = true;
                    tier.bumped = false;
                    tier.updateLabel();
                }
                tier.browser.refreshTier(tier);
            }
        });
    }
}

DasTier.prototype.setStylesheet = function(ss) {
    this.baseStylesheet = shallowCopy(ss);
    for (var si = 0; si < this.baseStylesheet.styles.length; ++si) {
        var sh = this.baseStylesheet.styles[si] = shallowCopy(this.baseStylesheet.styles[si]);
        sh._methodRE = sh._labelRE = sh._typeRE = null;
        sh.style = shallowCopy(sh.style);
        sh.style.id = 'style' + (++this.styleIdSeed);
    }
    this._updateFromConfig();
}

DasTier.prototype.getSource = function() {
    return this.featureSource;
}

DasTier.prototype.getDesiredTypes = function(scale) {
    var fetchTypes = [];
    var inclusive = false;
    var ssScale = this.browser.zoomForCurrentScale();

    if (this.stylesheet) {
        // dlog('ss = ' + miniJSONify(this.stylesheet));
        var ss = this.stylesheet.styles;
        for (var si = 0; si < ss.length; ++si) {
            var sh = ss[si];
            if (!sh.zoom || sh.zoom == ssScale) {
                if (!sh.type || sh.type == 'default') {
                    inclusive = true;
                    break;
                } else {
                    pushnew(fetchTypes, sh.type);
                }
            }
        }
    } else {
        // inclusive = true;
        return undefined;
    }

    if (inclusive) {
        return null;
    } else {
        return fetchTypes;
    }
}

DasTier.prototype.needsSequence = function(scale ) {
    if (this.dasSource.tier_type === 'sequence' && scale < 5) {
        return true;
    } else if ((this.dasSource.bamURI || this.dasSource.bamBlob) && scale < 20) {
        return true
    }
    return false;
}

DasTier.prototype.viewFeatures = function(chr, coverage, scale, features, sequence) {
    this.currentFeatures = features;
    this.currentSequence = sequence;
    this.notifyFeaturesLoaded();
    
    this.knownChr = chr;
    this.knownCoverage = coverage;

    if (this.status) {
        this.status = null;
        this._notifierToStatus();
    }

    this.draw();
}

DasTier.prototype.draw = function() {
    var features = this.currentFeatures;
    var seq = this.currentSequence;
    if (this.dasSource.tier_type === 'sequence') {
        drawSeqTier(this, seq); 
    } else {
        drawFeatureTier(this);
    }
    this.paint();
    this.originHaxx = 0;
    this.browser.arrangeTiers();
}

DasTier.prototype.findNextFeature = function(chr, pos, dir, fedge, callback) {
    if (this.quantLeapThreshold) {
        var width = this.browser.viewEnd - this.browser.viewStart + 1;
        pos = (pos +  ((width * dir) / 2))|0
        this.featureSource.quantFindNextFeature(chr, pos, dir, this.quantLeapThreshold, callback);
    } else {
        if (this.knownCoverage && pos >= this.knownCoverage.min() && pos <= this.knownCoverage.max()) {
            if (this.currentFeatures) {
                var bestFeature = null;
                for (var fi = 0; fi < this.currentFeatures.length; ++fi) {
                    var f = this.currentFeatures[fi];
                    if (!f.min || !f.max) {
                        continue;
                    }
                    if (f.parents && f.parents.length > 0) {
                        continue;
                    }
                    if (dir < 0) {
                        if (fedge == 1 && f.max >= pos && f.min < pos) {
                            if (!bestFeature || f.min > bestFeature.min ||
                                (f.min == bestFeature.min && f.max < bestFeature.max)) {
                                bestFeature = f;
                            }
                        } else if (f.max < pos) {
                            if (!bestFeature || f.max > bestFeature.max || 
                                (f.max == bestFeature.max && f.min < bestFeature.min) ||
                                (f.min == bestFeature.mmin && bestFeature.max >= pos)) {
                                bestFeature = f;
                            } 
                        }
                    } else {
                        if (fedge == 1 && f.min <= pos && f.max > pos) {
                            if (!bestFeature || f.max < bestFeature.max ||
                                (f.max == bestFeature.max && f.min > bestFeature.min)) {
                                bestFeature = f;
                            }
                        } else if (f.min > pos) {
                            if (!bestFeature || f.min < bestFeature.min ||
                                (f.min == bestFeature.min && f.max > bestFeature.max) ||
                                (f.max == bestFeature.max && bestFeature.min <= pos)) {
                                bestFeature = f;
                            }
                        }
                    }
                }
                if (bestFeature) {
                    //                dlog('bestFeature = ' + miniJSONify(bestFeature));
                    return callback(bestFeature);
                }
                if (dir < 0) {
                    pos = this.knownStart;
                } else {
                    pos = this.knownEnd;
                }
            }
        }

        this.featureSource.findNextFeature(chr, pos, dir, callback);
    }
}


DasTier.prototype.updateLabel = function() {
   this.bumpButton.className = this.bumped ? 'fa fa-minus-circle' : 'fa fa-plus-circle';
   if (this.dasSource.collapseSuperGroups) {
        this.bumpButton.style.display = 'inline-block';
    } else {
        this.bumpButton.style.display = 'none';
    }
}

DasTier.prototype.updateHeight = function() {
    this.currentHeight = Math.max(Math.max(this.layoutHeight, this.label.clientHeight + 2), this.browser.minTierHeight);
    this.row.style.height = '' + this.currentHeight + 'px';
    this.browser.updateHeight();
 }

DasTier.prototype.drawOverlay = function() {
    var t = this;
    var b = this.browser;
    var retina = b.retina && window.devicePixelRatio > 1;
    var g = t.overlay.getContext('2d');
    
    t.overlay.height = t.viewport.height;
    t.overlay.width = t.viewport.width;
    if (retina) {
        g.scale(2, 2);
    }
    
    var origin = b.viewStart - (1000/b.scale);
    var visStart = b.viewStart - (1000/b.scale);
    var visEnd = b.viewEnd + (1000/b.scale);


    for (var hi = 0; hi < b.highlights.length; ++hi) {
        var h = b.highlights[hi];
        if (((h.chr === b.chr) || (h.chr === ('chr' + b.chr))) && h.min < visEnd && h.max > visStart) {
            g.globalAlpha = b.defaultHighlightAlpha;
            g.fillStyle = b.defaultHighlightFill;
            g.fillRect((h.min - origin) * b.scale,
                       0,
                       (h.max - h.min) * b.scale,
                       t.overlay.height);
        }
    }

    t.oorigin = b.viewStart;
    t.overlay.style.width = t.viewport.style.width;
    t.overlay.style.height = t.viewport.style.height;
    t.overlay.style.left = '-1000px'
}

DasTier.prototype.updateStatus = function(status) {
    if (status) {
        this.status = status;
        this.currentFeatures = [];
        this.currentSequence = null;
        this.draw();
        this.updateHeight();
        this._notifierToStatus();
    } else {
        if (this.status) {
            this.status = null
            this._notifierToStatus();
        }
    }
}

DasTier.prototype.notify = function(message, timeout) {
    if (typeof(timeout) !== 'number')
        timeout = 2000;

    if (this.notifierFadeTimeout) {
        clearTimeout(this.notifierFadeTimeout);
        this.notifierFadeTimeout = null;
    }

    if (message) {
        this._notifierOn(message);
        if (timeout > 0) {
            var thisB = this;
            this.notifierFadeTimeout = setTimeout(function() {
                thisB._notifierToStatus();
            }, timeout);
        }
    } else {
        this._notifierToStatus();
    }
}

DasTier.prototype._notifierOn = function(message) {
    this.notifier.textContent = message;
    this.notifier.style.opacity = 0.8;
}

DasTier.prototype._notifierOff = function() {
    this.notifier.style.opacity = 0;
}

DasTier.prototype._notifierToStatus = function() {
    if (this.status) {
        this._notifierOn(this.status)
    } else {
        this._notifierOff();
    }
}

DasTier.prototype.setConfig = function(config) {
    this.config = config || {};
    this._updateFromConfig();
    this.notifyTierListeners();
}

DasTier.prototype.mergeConfig = function(newConfig) {
    for (var k in newConfig) {
        this.config[k] = newConfig[k];
    }
    this._updateFromConfig();
    this.notifyTierListeners();
}

DasTier.prototype._updateFromConfig = function() {
    var needsRefresh = false;
    var needsReorder = false;

    if (typeof this.config.name === 'string')
        this.nameElement.textContent = this.config.name;
    else
        this.nameElement.textContent = this.dasSource.name;

    var wantedHeight = this.config.height || this.dasSource.forceHeight;
    if (wantedHeight != this.forceHeight) {
        this.forceHeight = wantedHeight;
        needsRefresh = true;
    }

    if (this.forceMinDynamic != this.config.forceMinDynamic) {
        this.forceMinDynamic = this.config.forceMinDynamic;
        needsRefresh = true;
    }

    var forceMin = this.config.forceMin != undefined ? this.config.forceMin : this.dasSource.forceMin;
    if (this.forceMin != forceMin) {
        this.forceMin = forceMin;
        needsRefresh = true;
    }

    if (this.forceMaxDynamic != this.config.forceMaxDynamic) {
        this.forceMaxDynamic = this.config.forceMaxDynamic;
        needsRefresh = true;
    }
    
    var forceMax = this.config.forceMax != undefined ? this.config.forceMax : this.dasSource.forceMax;
    if (this.forceMax != forceMax) {
        this.forceMax = forceMax;
        needsRefresh = true;
    }

    var quantLeapThreshold = null;
    if (this.config.quantLeapThreshold !== undefined)
        quantLeapThreshold = this.config.quantLeapThreshold;
    else if (this.dasSource.quantLeapThreshold !== undefined)
        quantLeapThreshold = this.dasSource.quantLeapThreshold;
    if (quantLeapThreshold != this.quantLeapThreshold) {
        this.quantLeapThreshold = quantLeapThreshold;
        needsRefresh = true;
    }
    
    // Possible FIXME -- are there cases where style IDs need to be reassigned?
    var stylesheet = this.config.stylesheet || this.baseStylesheet;
    if (this.stylesheet !== stylesheet) {
        this.stylesheet = stylesheet;
        needsRefresh = true;
    }

    var wantedPinned = this.config.pinned !== undefined ? this.config.pinned : this.dasSource.pinned;
    if (wantedPinned !== this.pinned) {
        this.pinned = wantedPinned;
        needsReorder = true;
    }

    if (needsRefresh)
        this.scheduleRedraw();

    if (needsReorder)
        this.browser.reorderTiers();
}

DasTier.prototype.scheduleRedraw = function() {
    if (!this.currentFeatures)
        return;
    
    var tier = this;

    if (!this.redrawTimeout) {
        this.redrawTimeout = setTimeout(function() {
            tier.draw();
            tier.redrawTimeout = null;
        }, 10);
    }
}


DasTier.prototype.addTierListener = function(l) {
    this.listeners.push(l);
}

DasTier.prototype.notifyTierListeners = function(change) {
    for (var li = 0; li < this.listeners.length; ++li) {
        try {
            this.listeners[li](change);
        } catch (e) {
            console.log(e);
        }
    }
    this.browser.notifyTier();
}

DasTier.prototype.addFeaturesLoadedListener = function(handler) {
    this.featuresLoadedListeners.push(handler);
}

DasTier.prototype.notifyFeaturesLoaded = function() {
    for (var li = 0; li < this.featuresLoadedListeners.length; ++li) {
        try {
            this.featuresLoadedListeners[li].call(this);
        } catch (e) {
            console.log(e);
        }
    }
}
/* -*- mode: javascript; c-basic-offset: 4; indent-tabs-mode: nil -*- */

// 
// Dalliance Genome Explorer
// (c) Thomas Down 2006-2010
//
// track-adder.js
//

"use strict";

Browser.prototype.currentlyActive = function(source) {
    for (var ti = 0; ti < this.tiers.length; ++ti) {
        if (sourcesAreEqual(this.tiers[ti].dasSource, source))
            return this.tiers[ti];
    }
    return false;
}

Browser.prototype.makeButton = function(name, tooltip) {
    var regButton = makeElement('a', name, {href: '#'});
    if (tooltip) {
        this.makeTooltip(regButton, tooltip);
    }
    return makeElement('li', regButton);
}

function activateButton(addModeButtons, which) {
    for (var i = 0; i < addModeButtons.length; ++i) {
        var b = addModeButtons[i];
        if (b === which) {
            b.classList.add('active');
        } else {
            b.classList.remove('active');
        }
    }
}

Browser.prototype.showTrackAdder = function(ev) {
    if (this.uiMode === 'add') {
        this.hideToolPanel();
        this.setUiMode('none');
        return;
    }

    var thisB = this;

    var popup = makeElement('div', null, {className: 'dalliance'} , {width: '100%', display: 'inline-block', boxSizing: 'border-box', MozBoxSizing: 'border-box', verticalAlign: 'top'});

    var addModeButtons = [];
    var makeStab, makeStabObserver;


    if (!this.noRegistryTabs) {
        var regButton = this.makeButton('Registry', 'Browse compatible datasources from the DAS registry');
        addModeButtons.push(regButton);
        
        for (var m in this.mappableSources) {
            var mf  = function(mm) {
                var mapButton = thisB.makeButton(thisB.chains[mm].srcTag, 'Browse datasources mapped from ' + thisB.chains[mm].srcTag);
                addModeButtons.push(mapButton);
                mapButton.addEventListener('click', function(ev) {
                    ev.preventDefault(); ev.stopPropagation();
                    activateButton(addModeButtons, mapButton);
                    makeStab(thisB.mappableSources[mm], mm);
                }, false);
            }; mf(m);
        }
    }

    var groupedDefaults = {};
    for (var si = 0; si < this.defaultSources.length; ++si) {
        var s = this.defaultSources[si];
        var g = s.group || 'Defaults';
        if (groupedDefaults[g]) {
            groupedDefaults[g].push(s);
        } else {
            groupedDefaults[g] = [s];
        }
    }
    

    var makeHubButton = function(tdb) {
        var hub = tdb.hub;
        var hubMenuButton = makeElement('i', null, {className: 'fa fa-list-alt'});
        var label = hub.shortLabel || 'Unknown';
        if (tdb.mapping)
            label = label + ' (' + tdb.genome + ')';
        var hbContent = makeElement('span', [label, ' ', hubMenuButton]);
        var hubButton = thisB.makeButton(hbContent, hub.longLabel);
        addModeButtons.push(hubButton);
        
        hubButton.addEventListener('click', function(ev) {
            ev.preventDefault(); ev.stopPropagation();
            activateButton(addModeButtons, hubButton);
            removeChildren(stabHolder);
            var loader = thisB.makeLoader(24);
            loader.style.marginLeft = 'auto';
            loader.style.marginRight = 'auto';
            loader.style.marginTop = '100px';
            stabHolder.appendChild(makeElement('div', loader, null, {textAlign: 'center'}));

            refreshButton.style.display = 'none';
            addButton.style.display = 'none';
            canButton.style.display = 'none';

            tdb.getTracks(function(tracks, err) {
                if (err) {
                    console.log(err);
                }
                
                makeHubStab(tracks);
            });
        }, false);

        hubMenuButton.addEventListener('click', function(ev) {
            ev.preventDefault(); ev.stopPropagation();
            
            var removeHubItem = makeElement('li', makeElement('a', 'Remove hub'));
            var allOnItem = makeElement('li',  makeElement('a', 'Enable all'));
            var allOffItem = makeElement('li',  makeElement('a', 'Disable all'));
            var hubMenu = makeElement('ul', [removeHubItem, allOnItem, allOffItem], {className: 'dropdown-menu'}, {display: 'block'});

            var mx =  ev.clientX, my = ev.clientY;
            mx +=  document.documentElement.scrollLeft || document.body.scrollLeft;
            my +=  document.documentElement.scrollTop || document.body.scrollTop;

            hubMenu.style.position = 'absolute';
            hubMenu.style.top = '' + (my+10) + 'px';
            hubMenu.style.left = '' + (mx-30) + 'px';
            thisB.hPopupHolder.appendChild(hubMenu);

            var clickCatcher = function(ev) {
                console.log('cc');
                document.body.removeEventListener('click', clickCatcher, true);
                thisB.hPopupHolder.removeChild(hubMenu);
            };
            document.body.addEventListener('click', clickCatcher, true);

            removeHubItem.addEventListener('click', function(ev) {
                for (var hi = 0; hi < thisB.hubObjects.length; ++hi) {
                    if (thisB.hubObjects[hi].absURL == tdb.absURL) {
                        thisB.hubObjects.splice(hi, 1);
                        break;
                    }
                }
                for (var hi = 0; hi < thisB.hubs.length; ++hi) {
                    var hc = thisB.hubs[hi];
                    if (typeof hc === 'string')
                        hc = {url: hc};
                    if (hc.url == tdb.hub.url && !hc.genome || hc.genome == tdb.genome) {
                        thisB.hubs.splice(hi, 1);
                        break;
                    }

                }

                thisB.notifyTier();

                modeButtonHolder.removeChild(hubButton);
                activateButton(addModeButtons, addHubButton);
                switchToHubConnectMode();
            }, false);


            allOnItem.addEventListener('click', function(ev) {
                tdb.getTracks(function(tracks, err) {
                    if (err) {
                        console.log(err);
                    }
                    
                    for (var ti = 0; ti < tracks.length; ++ti) {
                        var ds = tracks[ti].toDallianceSource();
                        if (!thisB.currentlyActive(ds)) {
                            thisB.addTier(ds);
                        }
                    }
                });
            }, false);

            allOffItem.addEventListener('click', function(ev) {
                tdb.getTracks(function(tracks, err) {
                    if (err) {
                        console.log(err);
                    }
                    
                    for (var ti = 0; ti < tracks.length; ++ti) {
                        var ds = tracks[ti].toDallianceSource();
                        if (thisB.currentlyActive(ds)) {
                            thisB.removeTier(ds);
                        }
                    }
                });
            }, false);
        }, false);

        return hubButton;
    }

    var firstDefButton = null;
    var firstDefSources = null;
    for (var g in groupedDefaults) {
        (function(g, ds) {
            var defButton = thisB.makeButton(g, 'Browse the default set of data for this browser');
            defButton.addEventListener('click', function(ev) {
                ev.preventDefault(); ev.stopPropagation();
                activateButton(addModeButtons, defButton);
                makeStab(new Observed(ds));
            }, false);
            addModeButtons.push(defButton);

            if (!firstDefButton) {
                firstDefButton = defButton;
                firstDefSources = ds;
            }
        })(g, groupedDefaults[g]);
    }   
    var custButton = this.makeButton('DAS', 'Add arbitrary DAS data');
    addModeButtons.push(custButton);
    var binButton = this.makeButton('Binary', 'Add data in bigwig or bigbed format');
    addModeButtons.push(binButton);


    for (var hi = 0; hi < this.hubObjects.length; ++hi) {
        var hub = this.hubObjects[hi];
        makeHubButton(hub);
    }

    var addHubButton = this.makeButton('+', 'Connect to a new track-hub');
    addModeButtons.push(addHubButton);


    var modeButtonHolder = makeElement('ul', addModeButtons, {className: 'nav nav-tabs'}, {marginBottom: '0px'});
    popup.appendChild(modeButtonHolder);
    
    var custURL, custName, custCS, custQuant, custFile, custUser, custPass;
    var customMode = false;
    var dataToFinalize = null;

    var asform = makeElement('form', null, {}, {display: 'inline-block', width: '100%'});
    asform.addEventListener('submit', function(ev) {
            ev.stopPropagation(); ev.preventDefault();
            doAdd();
            return false;
    }, true); 
    var stabHolder = makeElement('div');
    stabHolder.style.position = 'relative';
    stabHolder.style.overflow = 'scroll';
    // stabHolder.style.height = '500px';
    asform.appendChild(stabHolder);

    var __mapping;
    var __sourceHolder;


    makeStab = function(msources, mapping) {
        refreshButton.style.display = 'none';
        addButton.style.display = 'none';
        canButton.style.display = 'none';
        if (__sourceHolder) {
            __sourceHolder.removeListener(makeStabObserver);
        }
        __mapping = mapping;
        __sourceHolder = msources;
        __sourceHolder.addListenerAndFire(makeStabObserver);
    }

    makeStabObserver = function(msources) {
        customMode = false;
        var buttons = [];
        removeChildren(stabHolder);
        if (!msources) {
            stabHolder.appendChild(makeElement('p', 'Dalliance was unable to retrieve data source information from the DAS registry, please try again later'));
            return;
        }
        
        var stabBody = makeElement('tbody', null, {className: 'table table-striped table-condensed'}, {width: '100%'});
        var stab = makeElement('table', stabBody, {className: 'table table-striped table-condensed'}, {width: '100%', tableLayout: 'fixed'}); 
        var idx = 0;

        var sources = [];
        for (var i = 0; i < msources.length; ++i) {
            sources.push(msources[i]);
        }
        
        sources.sort(function(a, b) {
            return a.name.toLowerCase().trim().localeCompare(b.name.toLowerCase().trim());
        });

        for (var i = 0; i < sources.length; ++i) {
            var source = sources[i];
            var r = makeElement('tr');

            var bd = makeElement('td', null, {}, {width: '30px'});
            bd.style.textAlign = 'center';
            if (!source.props || source.props.cors) {
                var b = makeElement('input');
                b.type = 'checkbox';
                b.dalliance_source = source;
                if (__mapping) {
                    b.dalliance_mapping = __mapping;
                }
                // b.checked = thisB.currentlyActive(source);
                bd.appendChild(b);
                buttons.push(b);
                b.addEventListener('change', function(ev) {
                    if (ev.target.checked) {
                        thisB.addTier(ev.target.dalliance_source);
                    } else {
                        thisB.removeTier(ev.target.dalliance_source);
                    }
                });
            } else {
                bd.appendChild(document.createTextNode('!'));
                thisB.makeTooltip(bd, makeElement('span', ["This data source isn't accessible because it doesn't support ", makeElement('a', "CORS", {href: 'http://www.w3.org/TR/cors/'}), "."]));
            }
            r.appendChild(bd);
            var ld = makeElement('td');
            ld.appendChild(document.createTextNode(source.name));
            if (source.desc && source.desc.length > 0) {
                thisB.makeTooltip(ld, source.desc);
            }
            r.appendChild(ld);
            stabBody.appendChild(r);
            ++idx;
        }

        var setChecks = function() {
            for (var bi = 0; bi < buttons.length; ++bi) {
                var b = buttons[bi];
                var t = thisB.currentlyActive(b.dalliance_source);
                if (t) {
                    b.checked = true;
                } else {
                    b.checked = false;
                }
            }
        }
        setChecks();
        thisB.addTierListener(function(l) {
            setChecks();
        });

        stabHolder.appendChild(stab);
    };

    function makeHubStab(tracks) {
        refreshButton.style.display = 'none';
        addButton.style.display = 'none';
        canButton.style.display = 'none';

        customMode = false;
        removeChildren(stabHolder);
        
        var ttab = makeElement('div', null, {}, {width: '100%'});
        var sources = [];
        for (var i = 0; i < tracks.length; ++i) {
            sources.push(tracks[i]);
        }
        
        sources.sort(function(a, b) {
            return a.shortLabel.toLowerCase().trim().localeCompare(b.shortLabel.toLowerCase().trim());
        });

        var groups = [];
        var tops = [];
        
        for (var ti = 0; ti < sources.length; ++ti) {
            var track = sources[ti];
            if (track.children && track.children.length > 0 && track.container != 'multiWig') {
                groups.push(track);
            } else {
                tops.push(track);
            }
        }
        if (tops.length > 0) {
            groups.push({
                shortLabel: 'Others',
                priority: -100000000,
                children: tops});
        }

        groups.sort(THUB_COMPARE);
        
        var buttons = [];
        for (var gi = 0; gi < groups.length; ++gi) {
            var group = groups[gi];
            var dg = group;
            if (!dg.dimensions && dg._parent && dg._parent.dimensions)
                dg = dg._parent;

            var dprops = {}
            if (dg.dimensions) {
                var dtoks = dg.dimensions.split(/(\w+)=(\w+)/);
                for (var dti = 0; dti < dtoks.length - 2; dti += 3) {
                    dprops[dtoks[dti + 1]] = dtoks[dti + 2];
                }
            }

            if (dprops.dimX && dprops.dimY) {
                var dimX = dprops.dimX, dimY = dprops.dimY;
                var sgX = dg.subgroups[dimX];
                var sgY = dg.subgroups[dimY];
                
                var trks = {};
                for (var ci = 0; ci < group.children.length; ++ci) {
                    var child = group.children[ci];
                    var vX = child.sgm[dimX], vY = child.sgm[dimY];
                    if (!trks[vX])
                        trks[vX] = {};
                    trks[vX][vY] = child;
                }

                var matrix = makeElement('table', null, {className: 'table table-striped table-condensed'}, {tableLayout: 'fixed'});
                {
                    var header = makeElement('tr');
                    header.appendChild(makeElement('th', null, {}, {width: '150px', height: '100px'}));   // blank corner element
                    for (var si = 0; si < sgX.titles.length; ++si) {
                        var h = makeElement('th', makeElement('div', sgX.titles[si], {}, {transform: 'rotate(-60deg)', 
                                                                       transformOrigin: '0% 100%', 
                                                                       webkitTransform: 'rotate(-60deg) translate(20px,10px)', 
                                                                       webkitTransformOrigin: '0% 100%',
                                                                       textAlign: 'left'}), {}, {width: '35px',
                                                                                                 height: '100px',
                                                                                                 verticalAlign: 'bottom'})
                        header.appendChild(h);
                    }
                    matrix.appendChild(header);
                }

                var mbody = makeElement('tbody', null, {className: 'table table-striped table-condensed'})
                for (var yi = 0; yi < sgY.titles.length; ++yi) {
                    var vY = sgY.tags[yi];
                    var row = makeElement('tr');
                    row.appendChild(makeElement('th', sgY.titles[yi]), {});
                    
                    for (var xi = 0; xi < sgX.titles.length; ++xi) {
                        var vX = sgX.tags[xi];
                        var cell = makeElement('td');
                        if (trks[vX] && trks[vX][vY]) {
                            var track = trks[vX][vY];
                            var ds = track.toDallianceSource();
                            if (!ds)
                                continue;
                            
                            var r = makeElement('tr');
                            var bd = makeElement('td');
                            bd.style.textAlign = 'center';
                            
                            var b = makeElement('input');
                            b.type = 'checkbox';
                            b.dalliance_source = ds;
                            if (__mapping) {
                                b.dalliance_mapping = __mapping;
                            }
                            buttons.push(b);
                            cell.appendChild(b);
                            b.addEventListener('change', function(ev) {
                                if (ev.target.checked) {
                                    thisB.addTier(ev.target.dalliance_source);
                                } else {
                                    thisB.removeTier(ev.target.dalliance_source);
                                }
                            });

                        }
                        row.appendChild(cell);
                    } 
                    mbody.appendChild(row);
                }
                matrix.appendChild(mbody);
                ttab.appendChild(makeTreeTableSection(group.shortLabel, matrix, gi==0));                
            } else {
                var stabBody = makeElement('tbody', null, {className: 'table table-striped table-condensed'});
                var stab = makeElement('table', stabBody, {className: 'table table-striped table-condensed'}, {width: '100%', tableLayout: 'fixed'}); 
                var idx = 0;
            
                group.children.sort(THUB_COMPARE);
                for (var i = 0; i < group.children.length; ++i) {
                    var track = group.children[i];
                    var ds = track.toDallianceSource();
                    if (!ds)
                        continue;

                    var r = makeElement('tr');
                    var bd = makeElement('td', null, {}, {width: '30px'});
                    bd.style.textAlign = 'center';
                    
                    var b = makeElement('input');
                    b.type = 'checkbox';
                    b.dalliance_source = ds;
                    if (__mapping) {
                        b.dalliance_mapping = __mapping;
                    }
                    buttons.push(b);
                    bd.appendChild(b);
                    b.addEventListener('change', function(ev) {
                        if (ev.target.checked) {
                            thisB.addTier(ev.target.dalliance_source);
                        } else {
                            thisB.removeTier(ev.target.dalliance_source);
                        }
                    });

                    r.appendChild(bd);
                    var ld = makeElement('td');
                    ld.appendChild(document.createTextNode(track.shortLabel));
                    if (track.longLabel && track.longLabel.length > 0) {
                        thisB.makeTooltip(ld, track.longLabel);
                    }
                    r.appendChild(ld);
                    stabBody.appendChild(r);
                    ++idx;
                }

                if (groups.length > 1 || group.shortLabel !== 'Others') {
                    ttab.appendChild(makeTreeTableSection(group.shortLabel, stab, gi==0));
                } else {
                    ttab.appendChild(stab);
                }
                
            }
        }

        var setChecks = function() {
            for (var bi = 0; bi < buttons.length; ++bi) {
                var b = buttons[bi];
                var t = thisB.currentlyActive(b.dalliance_source);
                if (t) {
                    b.checked = true;
                    b.disabled = t.sequenceSource != null;
                } else {
                    b.checked = false;
                }
            }
        }
        setChecks();
        thisB.addTierListener(function(l) {
            setChecks();
        });
        
        stabHolder.appendChild(ttab);
    }

    if (regButton) {
        regButton.addEventListener('click', function(ev) {
            ev.preventDefault(); ev.stopPropagation();
            activateButton(addModeButtons, regButton);
            makeStab(thisB.availableSources);
        }, false);
    }
 
    binButton.addEventListener('click', function(ev) {
        ev.preventDefault(); ev.stopPropagation();
        switchToBinMode();
    }, false);
    addHubButton.addEventListener('click', function(ev) {
        ev.preventDefault(); ev.stopPropagation();
        switchToHubConnectMode();
    }, false);


    function switchToBinMode() {
        activateButton(addModeButtons, binButton);
        customMode = 'bin';

        refreshButton.style.display = 'none';
        addButton.style.display = 'inline';
        canButton.style.display = 'none';

        removeChildren(stabHolder);
        var pageHolder = makeElement('div', null, {}, {paddingLeft: '10px', paddingRight: '10px'});
        pageHolder.appendChild(makeElement('h3', 'Add custom URL-based data'));
        pageHolder.appendChild(makeElement('p', ['You can add indexed binary data hosted on an web server that supports CORS (', makeElement('a', 'full details', {href: 'http://www.biodalliance.org/bin.html'}), ').  Currently supported formats are bigwig, bigbed, and indexed BAM.']));

        pageHolder.appendChild(makeElement('br'));
        pageHolder.appendChild(document.createTextNode('URL: '));
        custURL = makeElement('input', '', {size: 80, value: 'http://www.biodalliance.org/datasets/ensGene.bb'}, {width: '100%'});
        pageHolder.appendChild(custURL);
        
        pageHolder.appendChild(makeElement('br'));
        pageHolder.appendChild(makeElement('b', '- or -'));
        pageHolder.appendChild(makeElement('br'));
        pageHolder.appendChild(document.createTextNode('File: '));
        custFile = makeElement('input', null, {type: 'file', multiple: 'multiple'});
        pageHolder.appendChild(custFile);
        
        pageHolder.appendChild(makeElement('p', 'Clicking the "Add" button below will initiate a series of test queries.'));

        stabHolder.appendChild(pageHolder);
        custURL.focus();
    }

    function switchToHubConnectMode() {
        activateButton(addModeButtons, addHubButton);
        refreshButton.style.display = 'none';
        addButton.style.display = 'inline';
        canButton.style.display = 'none';

        customMode = 'hub-connect';
        refreshButton.style.visibility = 'hidden';

        removeChildren(stabHolder);

        var pageHolder = makeElement('div', null, {}, {paddingLeft: '10px', paddingRight: '10px'});
        pageHolder.appendChild(makeElement('h3', 'Connect to a track hub.'));
        pageHolder.appendChild(makeElement('p', ['Enter the top-level URL (usually points to a file called "hub.txt") of a UCSC-style track hub']));
        
        custURL = makeElement('input', '', {size: 120, value: 'http://www.biodalliance.org/datasets/testhub/hub.txt'}, {width: '100%'});
        pageHolder.appendChild(custURL);
        
        stabHolder.appendChild(pageHolder);
        
        custURL.focus();
    }

    custButton.addEventListener('click', function(ev) {
        ev.preventDefault(); ev.stopPropagation();
        switchToCustomMode();
    }, false);

    function switchToCustomMode() {
        activateButton(addModeButtons, custButton);
        refreshButton.style.display = 'none';
        addButton.style.display = 'inline';
        canButton.style.display = 'none';

        customMode = 'das';

        removeChildren(stabHolder);

        var customForm = makeElement('div', null, {},  {paddingLeft: '10px', paddingRight: '10px'});
        customForm.appendChild(makeElement('h3', 'Add custom DAS data'));
        customForm.appendChild(makeElement('p', 'This interface is intended for adding custom or lab-specific data.  Public data can be added more easily via the registry interface.'));
                
        customForm.appendChild(document.createTextNode('URL: '));
        customForm.appendChild(makeElement('br'));
        custURL = makeElement('input', '', {size: 80, value: 'http://www.derkholm.net:8080/das/medipseq_reads/'}, {width: '100%'});
        customForm.appendChild(custURL);

        customForm.appendChild(makeElement('p', 'Clicking the "Add" button below will initiate a series of test queries.  If the source is password-protected, you may be prompted to enter credentials.'));
        stabHolder.appendChild(customForm);

        custURL.focus();
    }



    var addButton = makeElement('button', 'Add', {className: 'btn btn-primary'});
    addButton.addEventListener('click', function(ev) {
        ev.stopPropagation(); ev.preventDefault();
        doAdd();
    }, false);

    function doAdd() {
        if (customMode) {
            if (customMode === 'das') {
                var curi = custURL.value.trim();
                if (!/^.+:\/\//.exec(curi)) {
                    curi = 'http://' + curi;
                }
                var nds = new DASSource({name: 'temporary', uri: curi});
                tryAddDAS(nds);
            } else if (customMode === 'bin') {
                var fileList = custFile.files;

                if (fileList && fileList.length > 0) {
                    tryAddMultiple(fileList);
                } else {
                    var curi = custURL.value.trim();
                    if (!/^.+:\/\//.exec(curi)) {
                        curi = 'http://' + curi;
                    }
                    tryAddBin({uri: curi});
                }
            } else if (customMode === 'reset') {
                switchToCustomMode();
            } else if (customMode === 'reset-bin') {
                switchToBinMode(); 
            } else if (customMode === 'reset-hub') {
                switchToHubConnectMode();
            } else if (customMode === 'prompt-bai') {
                var fileList = custFile.files;
                if (fileList && fileList.length > 0 && fileList[0]) {
                    dataToFinalize.baiBlob = fileList[0];
                    completeBAM(dataToFinalize);
                } else {
                    promptForBAI(dataToFinalize);
                }
            } else if (customMode === 'prompt-tbi') {
                var fileList = custFile.files;
                if (fileList && fileList.length > 0 && fileList[0]) {
                    dataToFinalize.indexBlob = fileList[0];
                    completeTabixVCF(dataToFinalize);
                } else {
                    promptForTabix(dataToFinalize);
                }
            } else if (customMode === 'finalize' || customMode === 'finalize-bin') {
                dataToFinalize.name = custName.value;
                var m = custCS.value;
                if (m != '__default__') {
                    dataToFinalize.mapping = m;
                } else {
                    dataToFinalize.mapping = undefined;
                }
                if (custQuant) {
                    dataToFinalize.maxbins = custQuant.checked;
                }

                if (custUser.value.length > 1 && custPass.value.length > 1) {
                    dlog('password');
                    dataToFinalize.xUser = custUser.value;
                    dataToFinalize.xPass = custPass.value;
                }

                thisB.addTier(dataToFinalize);

                if (customMode == 'finalize-bin')
                    switchToBinMode();
                else
                    switchToCustomMode();
            } else if (customMode === 'hub-connect') {
                var curi = custURL.value.trim();
                if (!/^.+:\/\//.exec(curi)) {
                    curi = 'http://' + curi;
                }
                
                tryAddHub(curi);
            } else if (customMode === 'multiple') {
                for (var mi = 0; mi < multipleSet.length; ++mi) {
                    var s = multipleSet[mi];
                    if (s.hidden)
                        continue;

                    if (s.tier_type == 'bam' && !s.indexBlob && !s.indexUri)
                        continue;
                    if (s.tier_type == 'tabix' && !s.indexBlob && !s.indexUri)
                        continue;

                    var nds = makeSourceConfig(s);
                    if (nds) {
                        nds.noPersist = true;
                        thisB.addTier(nds);
                    }
                }

                switchToBinMode();
            }
        } else {
            thisB.removeAllPopups();
        }
    };

    function tryAddHub(curi, opts, retry) {
        opts = opts || {};
        
        connectTrackHub(curi, function(hub, err) {
            if (err) {
                if (!retry) {
                    return tryAddHub(curi, {credentials: true}, true);
                }
                removeChildren(stabHolder);
                stabHolder.appendChild(makeElement('h2', 'Error connecting to track hub'))
                stabHolder.appendChild(makeElement('p', err));
                customMode = 'reset-hub';
                return;
            } else {
                var bestHub = null;
                var bestHubButton = null;
                for (var genome in hub.genomes) {
                    var mapping = null;
                    var okay = false;

                    if (genome == thisB.coordSystem.ucscName) {
                        okay = true;
                    } else {
                         for (var mid in thisB.chains) {
                            var m = thisB.chains[mid];
                            if (genome == m.coords.ucscName) {
                                mapping = mid;
                                okay = true;
                            }
                         }
                    }

                    if (okay) {
                        var hc = {url: curi, genome: genome};
                        if (opts.credentials)
                            hc.credentials = true;
                        if (mapping) {
                            hc.mapping = mapping;
                            hub.genomes[genome].mapping = mapping;
                        }
                        thisB.hubs.push(hc);
                        thisB.hubObjects.push(hub.genomes[genome]);
                        
                        var hubButton = makeHubButton(hub.genomes[genome]);
                        modeButtonHolder.appendChild(hubButton);

                        if (!mapping || !bestHub) {
                            bestHub = hub.genomes[genome];
                            bestHubButton = hubButton;
                        }
                    }
                }

                if (bestHub) {
                    thisB.notifyTier();
                    activateButton(addModeButtons, bestHubButton);
                    bestHub.getTracks(function(tracks, err) {
                        makeHubStab(tracks);
                    });
                } else {
                    removeChildren(stabHolder);
                    stabHolder.appendChild(makeElement('h2', 'No data for this genome'))
                    stabHolder.appendChild(makeElement('p', 'This URL appears to be a valid track-hub, but it doesn\'t contain any data for the coordinate system of this browser'));
                    stabHolder.appendChild(makeElement('p', 'coordSystem.ucscName = ' + thisB.coordSystem.ucscName));
                    customMode = 'reset-hub';
                    return;
                }
            }
        }, opts);
    }

    var tryAddDAS = function(nds, retry) {
        var knownSpace = thisB.knownSpace;
        if (!knownSpace) {
            alert("Can't confirm track-addition to an uninit browser.");
            return;
        }
        var tsm = Math.max(knownSpace.min, (knownSpace.min + knownSpace.max - 100) / 2)|0;
        var testSegment = new DASSegment(knownSpace.chr, tsm, Math.min(tsm + 99, knownSpace.max));
        nds.features(testSegment, {}, function(features, status) {
            if (status) {
                if (!retry) {
                    dlog('retrying with credentials');
                    nds.credentials = true;
                    tryAddDAS(nds, true);
                } else {
                    removeChildren(stabHolder);
                    stabHolder.appendChild(makeElement('h2', 'Custom data not found'));
                    stabHolder.appendChild(makeElement('p', 'DAS uri: ' + nds.uri + ' is not answering features requests'));
                    customMode = 'reset';
                    return;
                }
            } else {
                var nameExtractPattern = new RegExp('/([^/]+)/?$');
                var match = nameExtractPattern.exec(nds.uri);
                if (match) {
                    nds.name = match[1];
                }

                tryAddDASxSources(nds);
                return;
            }
        });
    }

    function tryAddDASxSources(nds, retry) {
        var uri = nds.uri;
        if (retry) {
            var match = /(.+)\/[^\/]+\/?/.exec(uri);
            if (match) {
                uri = match[1] + '/sources';
            }
        }
        function sqfail() {
            if (!retry) {
                return tryAddDASxSources(nds, true);
            } else {
                return addDasCompletionPage(nds);
            }
        }
        new DASRegistry(uri, {credentials: nds.credentials}).sources(
            function(sources) {
                if (!sources || sources.length == 0) {
                    return sqfail();
                } 

                var fs = null;
                if (sources.length == 1) {
                    fs = sources[0];
                } else {
                    for (var i = 0; i < sources.length; ++i) {
                        if (sources[i].uri === nds.uri) {
                            fs = sources[i];
                            break;
                        }
                    }
                }

                var coordsDetermined = false, quantDetermined = false;
                if (fs) {
                    nds.name = fs.name;
                    nds.desc = fs.desc;
                    if (fs.maxbins) {
                        nds.maxbins = true;
                    } else {
                        nds.maxbins = false;
                    }
                    if (fs.capabilities) {
                        nds.capabilities = fs.capabilities;
                    }
                    quantDetermined = true
                    
                    if (fs.coords && fs.coords.length == 1) {
                        var coords = fs.coords[0];
                        if (coordsMatch(coords, thisB.coordSystem)) {
                            coordsDetermined = true;
                        } else if (thisB.chains) {
                            for (var k in thisB.chains) {
                                if (coordsMatch(coords, thisB.chains[k].coords)) {
                                    nds.mapping = k;
                                    coordsDetermined = true;
                                }
                            }
                        }
                    }
                    
                }
                return addDasCompletionPage(nds, coordsDetermined, quantDetermined);
            },
            function() {
                return sqfail();
            }
        );
    }

    var makeSourceConfig = function(s) {
        var nds = {name: s.name};
        if (s.mapping && s.mapping != '__default__')
            nds.mapping = s.mapping;

        if (s.tier_type == 'bwg') {
            if (s.blob)
                nds.bwgBlob = s.blob;
            else if (s.uri)
                nds.bwgURI = s.uri;
            return nds;
        } else if (s.tier_type == 'bam') {
            if (s.blob) {
                nds.bamBlob = s.blob;
                nds.baiBlob = s.indexBlob;
            } else {
                nds.bamURI = s.uri;
                nds.baiURI = s.indexUri;
            }
            return nds;
        } else if (s.tier_type == 'tabix') {
            nds.tier_type = 'tabix';
            nds.payload = s.payload;
            if (s.blob) {
                nds.blob = s.blob;
                nds.indexBlob = s.indexBlob;
            } else {
                nds.uri = s.uri;
                nds.indexUri = s.indexUri;
            }
            return nds;
        } else if (s.tier_type == 'memstore') {
            nds.tier_type = 'memstore';
            nds.payload = s.payload;
            if (s.blob)
                nds.blob = s.blob;
            else
                nds.uri = s.uri;
            return nds;
        }
    }

    var tryAddBin = function(source) {
        probeResource(source, function(source, err) {
            if (err) {
                removeChildren(stabHolder);
                stabHolder.appendChild(makeElement('h2', "Couldn't access custom data"));
                stabHolder.appendChild(makeElement('p', '' + err));
                customMode = 'reset-bin';
            } else {
                var nds = makeSourceConfig(source);
                if (source.tier_type == 'bam') {
                    return completeBAM(nds);
                } else if (source.tier_type == 'tabix') {
                    return completeTabixVCF(nds);
                } else {
                    return addDasCompletionPage(nds, false, false, true);
                }
            }
        });
    }

    function promptForBAI(nds) {
        refreshButton.style.display = 'none';
        addButton.style.display = 'inline';
        canButton.style.display = 'inline';

        removeChildren(stabHolder);
        customMode = 'prompt-bai'
        stabHolder.appendChild(makeElement('h2', 'Select an index file'));
        stabHolder.appendChild(makeElement('p', 'Dalliance requires a BAM index (.bai) file when displaying BAM data.  These normally accompany BAM files.  For security reasons, web applications like Dalliance can only access local files which you have explicity selected.  Please use the file chooser below to select the appropriate BAI file'));

        stabHolder.appendChild(document.createTextNode('Index file: '));
        custFile = makeElement('input', null, {type: 'file'});
        stabHolder.appendChild(custFile);
        dataToFinalize = nds;
    }

    function promptForTabix(nds) {
        refreshButton.style.display = 'none';
        addButton.style.display = 'inline';
        canButton.style.display = 'inline';

        removeChildren(stabHolder);
        customMode = 'prompt-tbi'
        stabHolder.appendChild(makeElement('h2', 'Select an index file'));
        stabHolder.appendChild(makeElement('p', 'Dalliance requires a Tabix index (.tbi) file when displaying VCF data.  For security reasons, web applications like Dalliance can only access local files which you have explicity selected.  Please use the file chooser below to select the appropriate BAI file'));

        stabHolder.appendChild(document.createTextNode('Index file: '));
        custFile = makeElement('input', null, {type: 'file'});
        stabHolder.appendChild(custFile);
        dataToFinalize = nds;
    }

    function completeBAM(nds) {
        var indexF;
        if (nds.baiBlob) {
            indexF = new BlobFetchable(nds.baiBlob);
        } else {
            indexF = new URLFetchable(nds.bamURI + '.bai');
        }
        indexF.slice(0, 256).fetch(function(r) {
                var hasBAI = false;
                if (r) {
                    var ba = new Uint8Array(r);
                    var magic2 = readInt(ba, 0);
                    hasBAI = (magic2 == BAI_MAGIC);
                }
                if (hasBAI) {
                    return addDasCompletionPage(nds, false, false, true);
                } else {
                    return binFormatErrorPage('You have selected a valid BAM file, but a corresponding index (.bai) file was not found.  Please index your BAM (samtools index) and place the BAI file in the same directory');
                }
        });
    }

    function completeTabixVCF(nds) {
        var indexF;
        if (nds.indexBlob) {
            indexF = new BlobFetchable(nds.indexBlob);
        } else {
            indexF = new URLFetchable(nds.uri + '.tbi');
        }
        indexF.slice(0, 1<<16).fetch(function(r) {
            var hasTabix = false;
            if (r) {
                var ba = new Uint8Array(r);
                if (ba[0] == 31 || ba[1] == 139) {
                    var unc = unbgzf(r);
                    ba = new Uint8Array(unc);
                    var m2 = readInt(ba, 0);
                    hasTabix = (m2 == TABIX_MAGIC);
                }
            }
            if (hasTabix) {
                return addDasCompletionPage(nds, false, false, true);
            } else {
                return binFormatErrorPage('You have selected a valid VCF file, but a corresponding index (.tbi) file was not found.  Please index your VCF ("tabix -p vcf -f myfile.vcf.gz") and place the .tbi file in the same directory');
            }
        });
    }

    function binFormatErrorPage(message) {
        refreshButton.style.display = 'none';
        addButton.style.display = 'inline';
        canButton.style.display = 'inline';

        removeChildren(stabHolder);
        message = message || 'Custom data format not recognized';
        stabHolder.appendChild(makeElement('h2', 'Error adding custom data'));
        stabHolder.appendChild(makeElement('p', message));
        stabHolder.appendChild(makeElement('p', 'Currently supported formats are bigBed, bigWig, and BAM.'));
        customMode = 'reset-bin';
        return;
    }
                     
    var addDasCompletionPage = function(nds, coordsDetermined, quantDetermined, quantIrrelevant) {
        refreshButton.style.display = 'none';
        addButton.style.display = 'inline';
        canButton.style.display = 'inline';

        removeChildren(stabHolder);
        stabHolder.appendChild(makeElement('h2', 'Add custom data: step 2'));
        stabHolder.appendChild(document.createTextNode('Label: '));
        custName = makeElement('input', '', {value: nds.name});
        stabHolder.appendChild(custName);


        // stabHolder.appendChild(document.createTextNode('User: '));
        custUser = makeElement('input', '');
        // stabHolder.appendChild(custUser);
        //stabHolder.appendChild(document.createTextNode('Pass: '));
        custPass = makeElement('input', '');
        // stabHolder.appendChild(custPass);
        

        stabHolder.appendChild(makeElement('br'));
        stabHolder.appendChild(makeElement('br'));
        stabHolder.appendChild(makeElement('h4', 'Coordinate system: '));
        custCS = makeElement('select', null);
        custCS.appendChild(makeElement('option', thisB.nameForCoordSystem(thisB.coordSystem), {value: '__default__'}));
        if (thisB.chains) {
            for (var csk in thisB.chains) {
                var cs = thisB.chains[csk].coords;
                custCS.appendChild(makeElement('option', thisB.nameForCoordSystem(cs), {value: csk}));
            }
        }
        custCS.value = nds.mapping || '__default__';
        stabHolder.appendChild(custCS);

        if (coordsDetermined) {
            stabHolder.appendChild(makeElement('p', "(Based on server response, probably doesn't need changing.)"));
        } else {
            stabHolder.appendChild(makeElement('p', [makeElement('b', 'Warning: '), "unable to determine the correct value from server responses.  Please check carefully."]));
            stabHolder.appendChild(makeElement('p', "If you don't see the mapping you're looking for, please contact thomas@biodalliance.org"));
        }

        if (!quantIrrelevant) {
            stabHolder.appendChild(document.createTextNode('Quantitative: '));
            custQuant = makeElement('input', null, {type: 'checkbox', checked: true});
            if (typeof nds.maxbins !== 'undefined') {
                custQuant.checked = nds.maxbins;
            }
            stabHolder.appendChild(custQuant);
            if (quantDetermined) {
                stabHolder.appendChild(makeElement('p', "(Based on server response, probably doesn't need changing.)"));
            } else {
                stabHolder.appendChild(makeElement('p', [makeElement('b', "Warning: "), "unable to determine correct value.  If in doubt, leave checked."]));
            }
        }

        if (nds.bwgBlob) {
            stabHolder.appendChild(makeElement('p', [makeElement('b', 'Warning: '), 'data added from local file.  Due to the browser security model, the track will disappear if you reload Dalliance.']));
        }

        custName.focus();

        if (customMode === 'bin' || customMode === 'prompt-bai' || customMode === 'prompt-tbi')
            customMode = 'finalize-bin';
        else
            customMode = 'finalize';
        dataToFinalize = nds;
    }

    var multipleSet = null;
    var tryAddMultiple = function(fileList) {
        var newSources = multipleSet = [];
        customMode = 'multiple';
        for (var fi = 0; fi < fileList.length; ++fi) {
            var f = fileList[fi];
            if (f) {
                newSources.push({blob: f});
            }
        }

        for (var fi = 0; fi < newSources.length; ++fi) {
            probeMultiple(newSources[fi]);
        }
        updateMultipleStatus();
    }

    var probeMultiple = function(ns) {
        probeResource(ns, function(source, err) {
            if (err) {
                source.error = err;
            }

            var usedIndices = [];
            var bams = {}, tabixes = {};
            for (var si = 0; si < multipleSet.length; ++si) {
                var s = multipleSet[si];
                if (s.tier_type == 'bam' && !s.indexBlob) {
                    bams[s.blob.name] = s;
                }
                if (s.tier_type == 'tabix' && !s.indexBlob) {
                    tabixes[s.blob.name] = s;
                }
            }

            for (var si = 0; si < multipleSet.length; ++si) {
                var s = multipleSet[si];
                if (s.tier_type === 'bai') {
                    var baiPattern = new RegExp('(.+)\\.bai$');
                    var match = baiPattern.exec(s.blob.name);
                    if (match && bams[match[1]]) {
                        bams[match[1]].indexBlob = s.blob;
                        usedIndices.push(si);
                    }
                } else if (s.tier_type === 'tabix-index') {
                    var tbiPattern = new RegExp('(.+)\\.tbi$');
                    var match = tbiPattern.exec(s.blob.name);
                    if (match && tabixes[match[1]]) {
                        tabixes[match[1]].indexBlob = s.blob;
                        usedIndices.push(si);
                    }
                }
            }

            for (var bi = usedIndices.length - 1; bi >= 0; --bi) {
                multipleSet.splice(usedIndices[bi], 1);
            }

            updateMultipleStatus();
        });
    }

    var updateMultipleStatus = function() {
        removeChildren(stabHolder);
        var needsIndex = false;
        var multTable = makeElement('table', multipleSet
          .filter(function(s) {return !s.hidden})
          .map(function(s) {
            var row = makeElement('tr');
            row.appendChild(makeElement('td', s.name || s.blob.name));
            var typeContent;
            if (s.error) {
                typeContent = makeElement('span', 'Error', null, {color: 'red'});
            } else if (s.tier_type) {
                typeContent = s.payload || s.tier_type;
            } else {
                typeContent = thisB.makeLoader(16);
            }

            var ccs;
            var state = 'unknown';
            if (s.tier_type == 'bwg' || s.tier_type == 'memstore') {
                state = 'okay';
            } else if (s.tier_type == 'bam') {
                state = s.indexBlob ? 'okay' : 'needs-index';
            } else if (s.tier_type == 'tabix') {
                state = s.indexBlob ? 'okay' : 'needs-index';
            }

            if (state == 'okay') {
                ccs = makeElement('select', null, null, {width: '150px'});
                ccs.appendChild(makeElement('option', thisB.nameForCoordSystem(thisB.coordSystem), {value: '__default__'}));
                if (thisB.chains) {
                    for (var csk in thisB.chains) {
                        var cs = thisB.chains[csk].coords;
                        ccs.appendChild(makeElement('option', thisB.nameForCoordSystem(cs), {value: csk}));
                    }
                }
                ccs.value = s.mapping || '__default__';

                ccs.addEventListener('change', function(ev) {
                    s.mapping = ccs.value;
                    console.log(s);
                }, false);
            } else if (state == 'needs-index') {
                ccs = makeElement('span', 'Needs index', {}, {color: 'red'});
                needsIndex = true;
            }

            return makeElement('tr', [makeElement('td', s.name || s.blob.name),
                                      makeElement('td', typeContent),
                                      makeElement('td', ccs)]);

        }), {className: 'table table-striped table-condensed'});
        stabHolder.appendChild(multTable);

        if (needsIndex) {
            stabHolder.appendChild(makeElement('p', 'Some of these files are missing required index (.bai or .tbi) files.  For security reasons, web applications like Dalliance can only access local files which you have explicity selected.  Please use the file chooser below to select the appropriate index file'));
            stabHolder.appendChild(document.createTextNode('Index file(s): '));
            var indexFile = makeElement('input', null, {type: 'file', multiple: 'multiple'});
            stabHolder.appendChild(indexFile);
            indexFile.addEventListener('change', function(ev) {
                console.log('fileset changed');
                var fileList = indexFile.files || [];
                for (var fi = 0; fi < fileList.length; ++fi) {
                    var f = fileList[fi];
                    if (f) {
                        var ns = {blob: f, hidden: true};
                        multipleSet.push(ns);
                        probeMultiple(ns);
                    }
                }
            }, false);
        }
    }

    var canButton = makeElement('button', 'Cancel', {className: 'btn'});
    canButton.addEventListener('click', function(ev) {
        ev.stopPropagation(); ev.preventDefault();
        if (customMode === 'finalize-bin')
            switchToBinMode();
        else
            switchToCustomMode();
    }, false);

    var refreshButton = makeElement('button', 'Refresh', {className: 'btn'});
    refreshButton.addEventListener('click', function(ev) {
        ev.stopPropagation(); ev.preventDefault();
        thisB.queryRegistry(__mapping);
    }, false);
    this.makeTooltip(refreshButton, 'Click to re-fetch data from the DAS registry');

    var buttonHolder = makeElement('div', [addButton, ' ', canButton, ' ', refreshButton]);
    buttonHolder.style.margin = '10px';
    asform.appendChild(buttonHolder);

    popup.appendChild(asform);
    makeStab(thisB.availableSources);

    this.showToolPanel(popup);
    this.setUiMode('add');

    if (firstDefButton) {
        activateButton(addModeButtons, firstDefButton);
        makeStab(new Observed(firstDefSources));
    }
}
/* -*- mode: javascript; c-basic-offset: 4; indent-tabs-mode: nil -*- */

// 
// Dalliance Genome Explorer
// (c) Thomas Down 2006-2010
//
// twoBit.js: packed-binary reference sequences
//

"use strict";

var TWOBIT_MAGIC = 0x1a412743;
var TWOBIT_MAGIC_BE = 0x4327411a;

function TwoBitFile() {
}

function makeTwoBit(fetchable, cnt) {
    var tb = new TwoBitFile();
    tb.data = fetchable;
    var headerBlockSize = 12500;

    tb.data.slice(0, headerBlockSize).fetch(function(r) {
        if (!r) {
            return cnt(null, "Couldn't access data");
        }
        var ba = new Uint8Array(r);
        var magic = readInt(ba, 0);
        if (magic == TWOBIT_MAGIC) {
            tb.readInt = readInt;
        } else if (magic == TWOBIT_MAGIC_BE) {
            tb.readInt = readIntBE;
        } else {
            return cnt(null, "Not a .2bit file, magic=0x" + magic.toString(16));
        }

        var version = tb.readInt(ba, 4);
        if (version != 0) {
            return cnt(null, 'Unsupported version ' + version);
        }

        tb.seqCount = tb.readInt(ba, 8);
        tb.seqDict = {};

        var p = 16, i=0;
        var o = 0;  // Offset of the current block if we need to fetch multiple header blocks.

        var parseSeqInfo = function() {
            while (i < tb.seqCount) {
                var ns = ba[p];
                if (p + ns + 6 >= ba.length) {
                    return tb.data.slice(o + p, headerBlockSize).fetch(function (r) {
                        o += p;
                        p = 0;
                        ba = new Uint8Array(r);
                        parseSeqInfo();
                    });
                } else {
                    ++p;
                    var name = '';
                    for (var j = 1; j <= ns; ++j) {
                        name += String.fromCharCode(ba[p++]);
                    }
                    var offset = tb.readInt(ba, p);
                    p += 4;
                    tb.seqDict[name] = new TwoBitSeq(tb, offset);
                    ++i;
                }
            }
            return cnt(tb);
        }

        parseSeqInfo();
        
    });
}

TwoBitFile.prototype.getSeq = function(chr) {
    var seq = this.seqDict[chr];
    if (!seq) {
        seq = this.seqDict['chr' + chr];
    }
    return seq;
}

TwoBitFile.prototype.fetch = function(chr, min, max, cnt) {
    var seq = this.getSeq(chr);
    if (!seq) {
        return cnt(null, "Couldn't find " + chr);
    } else if (max <= min) {
        return cnt('');
    } else {
        seq.fetch(min, max, cnt);
    }
}

function TwoBitSeq(tbf, offset) {
    this.tbf = tbf;
    this.offset = offset;
}

TwoBitSeq.prototype.init = function(cnt) {
    if (this.seqOffset) {
        return cnt();
    }

    var thisB = this;
    thisB.tbf.data.slice(thisB.offset, 8).fetch(function(r1) {
        if (!r1) {
            return cnt('Fetch failed');
        }
        var ba = new Uint8Array(r1);
        thisB._length = thisB.tbf.readInt(ba, 0);
        thisB.nBlockCnt = thisB.tbf.readInt(ba, 4);
        thisB.tbf.data.slice(thisB.offset + 8, thisB.nBlockCnt*8 + 4).fetch(function(r2) {
            if (!r2) {
                return cnt('Fetch failed');
            }
            var ba = new Uint8Array(r2);
            var nbs = null;
            for (var b = 0; b < thisB.nBlockCnt; ++b) {
                var nbMin = thisB.tbf.readInt(ba, b * 4);
                var nbLen = thisB.tbf.readInt(ba, (b + thisB.nBlockCnt) * 4);
                var nb = new Range(nbMin, nbMin + nbLen - 1);
                if (!nbs) {
                    nbs = nb;
                } else {
                    nbs = union(nbs, nb);
                }
            }
            thisB.nBlocks = nbs;
            thisB.mBlockCnt = thisB.tbf.readInt(ba, thisB.nBlockCnt*8);
            thisB.seqLength = ((thisB._length + 3)/4)|0;
            thisB.seqOffset = thisB.offset + 16 + ((thisB.nBlockCnt + thisB.mBlockCnt) * 8);
            return cnt();
        });
    });
}

var TWOBIT_TABLE = ['T', 'C', 'A', 'G'];

TwoBitSeq.prototype.fetch = function(min, max, cnt) {
    --min; --max;       // Switch to zero-based.
    var thisB = this;
    this.init(function(error) {
        if (error) {
            return cnt(null, error);
        }

        var fetchMin = min >> 2;
        var fetchMax = max + 3 >> 2;
        if (fetchMin < 0 || fetchMax > thisB.seqLength) {
            return cnt('Coordinates out of bounds: ' + min + ':' + max);
        }

        thisB.tbf.data.slice(thisB.seqOffset + fetchMin, fetchMax - fetchMin).fetch(function(r) {
            if (r == null) {
                return cnt('SeqFetch failed');
            }
            var seqData = new Uint8Array(r);

            var nSpans = [];
            if (thisB.nBlocks) {
                var intr = intersection(new Range(min, max), thisB.nBlocks);
                if (intr) {
                    nSpans = intr.ranges();
                }
            }
            
            var seqstr = '';
            var ptr = min;
            function fillSeq(fsm) {
                while (ptr <= fsm) {
                    var bb = (ptr >> 2) - fetchMin;
                    var ni = ptr & 0x3;
                    var bv = seqData[bb];
                    var n;
                    if (ni == 0) {
                        n = (bv >> 6) & 0x3;
                    } else if (ni == 1) {
                        n = (bv >> 4) & 0x3;
                    } else if (ni == 2) {
                        n = (bv >> 2) & 0x3;
                    } else {
                        n = (bv) & 0x3;
                    }
                    seqstr += TWOBIT_TABLE[n];
                    ++ptr;
                }
            }
            
            for (var b = 0; b < nSpans.length; ++b) {
                var nb = nSpans[b];
                if (ptr > nb.min()) {
                    throw 'N mismatch...';
                }
                if (ptr < nb.min()) {
                    fillSeq(nb.min() - 1);
                }
                while (ptr <= nb.max()) {
                    seqstr += 'N';
                    ++ptr;
                }
            }
            if (ptr <= max) {
                fillSeq(max);
            }

            return cnt(seqstr);
        });
    });
}

TwoBitSeq.prototype.length = function(cnt) {
    var thisB = this;
    this.init(function(error) {
        if (error) {
            return cnt(null, error);
        } else {
            return cnt(thisB._length);
        }
    });
}
/* -*- mode: javascript; c-basic-offset: 4; indent-tabs-mode: nil -*- */

// 
// Dalliance Genome Explorer
// (c) Thomas Down 2006-2010
//
// utils.js: odds, sods, and ends.
//

"use strict";

var NUM_REGEXP = new RegExp('[0-9]+');

function stringToNumbersArray(str) {
    var nums = new Array();
    var m;
    while (m = NUM_REGEXP.exec(str)) {
        nums.push(m[0]);
        str=str.substring(m.index + (m[0].length));
    }
    return nums;
}

var STRICT_NUM_REGEXP = new RegExp('^[0-9]+$');

function stringToInt(str) {
    str = str.replace(new RegExp(',', 'g'), '');
    if (!STRICT_NUM_REGEXP.test(str)) {
        return null;
    }
    return str|0;
}

function pushnew(a, v) {
    for (var i = 0; i < a.length; ++i) {
        if (a[i] == v) {
            return;
        }
    }
    a.push(v);
}

function pusho(obj, k, v) {
    if (obj[k]) {
        obj[k].push(v);
    } else {
        obj[k] = [v];
    }
}

function pushnewo(obj, k, v) {
    var a = obj[k];
    if (a) {
        for (var i = 0; i < a.length; ++i) {    // indexOf requires JS16 :-(.
            if (a[i] == v) {
                return;
            }
        }
        a.push(v);
    } else {
        obj[k] = [v];
    }
}


function pick(a, b, c, d)
{
    if (a) {
        return a;
    } else if (b) {
        return b;
    } else if (c) {
        return c;
    } else if (d) {
        return d;
    }
}

function pushnew(l, o)
{
    for (var i = 0; i < l.length; ++i) {
        if (l[i] == o) {
            return;
        }
    }
    l.push(o);
}

function maybeConcat(a, b) {
    var l = [];
    if (a) {
        for (var i = 0; i < a.length; ++i) {
            pushnew(l, a[i]);
        }
    }
    if (b) {
        for (var i = 0; i < b.length; ++i) {
            pushnew(l, b[i]);
        }
    }
    return l;
}

function arrayIndexOf(a, x) {
    if (!a) {
        return -1;
    }

    for (var i = 0; i < a.length; ++i) {
        if (a[i] === x) {
            return i;
        }
    }
    return -1;
}

function arrayRemove(a, x) {
    var i = arrayIndexOf(a, x);
    if (i >= 0) {
        a.splice(i, 1);
        return true;
    }
    return false;
}

//
// DOM utilities
//


function makeElement(tag, children, attribs, styles)
{
    var ele = document.createElement(tag);
    if (children) {
        if (! (children instanceof Array)) {
            children = [children];
        }
        for (var i = 0; i < children.length; ++i) {
            var c = children[i];
            if (typeof c == 'string') {
                c = document.createTextNode(c);
            }
            ele.appendChild(c);
        }
    }
    
    if (attribs) {
        for (var l in attribs) {
            try {
                ele[l] = attribs[l];
            } catch (e) {
                console.log('error setting ' + l);
                throw(e);
            }
        }
    }
    if (styles) {
        for (var l in styles) {
            ele.style[l] = styles[l];
        }
    }
    return ele;
}

function makeElementNS(namespace, tag, children, attribs)
{
    var ele = document.createElementNS(namespace, tag);
    if (children) {
        if (! (children instanceof Array)) {
            children = [children];
        }
        for (var i = 0; i < children.length; ++i) {
            var c = children[i];
            if (typeof c == 'string') {
                c = document.createTextNode(c);
            }
            ele.appendChild(c);
        }
    }
    
    setAttrs(ele, attribs);
    return ele;
}

var attr_name_cache = {};

function setAttr(node, key, value)
{
    var attr = attr_name_cache[key];
    if (!attr) {
        var _attr = '';
        for (var c = 0; c < key.length; ++c) {
            var cc = key.substring(c, c+1);
            var lcc = cc.toLowerCase();
            if (lcc != cc) {
                _attr = _attr + '-' + lcc;
            } else {
                _attr = _attr + cc;
            }
        }
        attr_name_cache[key] = _attr;
        attr = _attr;
    }
    node.setAttribute(attr, value);
}

function setAttrs(node, attribs)
{
    if (attribs) {
        for (var l in attribs) {
            setAttr(node, l, attribs[l]);
        }
    }
}



function removeChildren(node)
{
    if (!node || !node.childNodes) {
        return;
    }

    while (node.childNodes.length > 0) {
        node.removeChild(node.firstChild);
    }
}



//
// WARNING: not for general use!
//

function miniJSONify(o, exc) {
    if (typeof o === 'undefined') {
        return 'undefined';
    } else if (o == null) {
        return 'null';
    } else if (typeof o == 'string') {
        return "'" + o + "'";
    } else if (typeof o == 'number') {
        return "" + o;
    } else if (typeof o == 'boolean') {
        return "" + o;
    } else if (typeof o == 'object') {
        if (o instanceof Array) {
            var s = null;
            for (var i = 0; i < o.length; ++i) {
                s = (s == null ? '' : (s + ', ')) + miniJSONify(o[i], exc);
            }
            return '[' + (s?s:'') + ']';
        } else {
            exc = exc || {};
            var s = null;
            for (var k in o) {
                if (exc[k])
                    continue;
                if (k != undefined && typeof(o[k]) != 'function') {
                    s = (s == null ? '' : (s + ', ')) + k + ': ' + miniJSONify(o[k], exc);
                }
            }
            return '{' + (s?s:'') + '}';
        }
    } else {
        return (typeof o);
    }
}

function shallowCopy(o) {
    var n = {};
    for (var k in o) {
        n[k] = o[k];
    }
    return n;
}

function Observed(x) {
    this.value = x;
    this.listeners = [];
}

Observed.prototype.addListener = function(f) {
    this.listeners.push(f);
}

Observed.prototype.addListenerAndFire = function(f) {
    this.listeners.push(f);
    f(this.value);
}

Observed.prototype.removeListener = function(f) {
    arrayRemove(this.listeners, f);
}

Observed.prototype.get = function() {
    return this.value;
}

Observed.prototype.set = function(x) {
    this.value = x;
    for (var i = 0; i < this.listeners.length; ++i) {
        this.listeners[i](x);
    }
}

function Awaited() {
    this.queue = [];
}

Awaited.prototype.provide = function(x) {
    if (this.res !== undefined) {
        throw "Resource has already been provided.";
    }

    this.res = x;
    for (var i = 0; i < this.queue.length; ++i) {
        this.queue[i](x);
    }
    this.queue = null;   // avoid leaking closures.
}

Awaited.prototype.await = function(f) {
    if (this.res !== undefined) {
        f(this.res);
        return this.res;
    } else {
        this.queue.push(f);
    }
}

var __dalliance_saltSeed = 0;

function saltURL(url) {
    return url + '?salt=' + b64_sha1('' + Date.now() + ',' + (++__dalliance_saltSeed));
}

function textXHR(url, callback, opts) {
    if (opts.salt) 
        url = saltURL(url);

    var req = new XMLHttpRequest();
    req.onreadystatechange = function() {
    	if (req.readyState == 4) {
    	    if (req.status >= 300) {
    		    callback(null, 'Error code ' + req.status);
    	    } else {
    		    callback(req.responseText);
    	    }
    	}
    };
    
    req.open('GET', url, true);
    req.responseType = 'text';

    if (opts && opts.credentials) {
        req.withCredentials = true;
    }
    req.send('');
}

function relativeURL(base, rel) {
    // FIXME quite naive -- good enough for trackhubs?

    if (rel.indexOf('http:') == 0 || rel.indexOf('https:') == 0) {
        return rel;
    }

    var li = base.lastIndexOf('/');
    if (li >= 0) {
        return base.substr(0, li + 1) + rel;
    } else {
        return rel;
    }
}


//
// Missing APIs
// 

if (!('trim' in String.prototype)) {
    String.prototype.trim = function() {
        return this.replace(/^\s+/, '').replace(/\s+$/, '');
    };
}
/* -*- mode: javascript; c-basic-offset: 4; indent-tabs-mode: nil -*- */

// 
// Dalliance Genome Explorer
// (c) Thomas Down 2006-2010
//
// version.js
//

"use strict";

var VERSION = {
    CONFIG: 5,
    MAJOR:  0,
    MINOR:  11,
    MICRO:  10,
    PATCH:  '',
    BRANCH: 'dev'
}

VERSION.toString = function() {
    var vs = '' + this.MAJOR + '.' + this.MINOR + '.' + this.MICRO;
    if (this.PATCH) {
        vs = vs + this.PATCH;
    }
    if (this.BRANCH && this.BRANCH != '') {
        vs = vs + '-' + this.BRANCH;
    }
    return vs;
}
/* -*- mode: javascript; c-basic-offset: 4; indent-tabs-mode: nil -*- */

// 
// Dalliance Genome Explorer
// (c) Thomas Down 2006-2013
//
// browser-us.js: standard UI wiring
//

"use strict";

function formatLongInt(n) {
    return (n|0).toString().replace(/\B(?=(\d{3})+(?!\d))/g, ',')
}

/*
 * Quite a bit of this ought to be done using a templating system, but
 * since web-components isn't quite ready for prime time yet we'll stick
 * with constructing it all in Javascript for now...
 */

Browser.prototype.initUI = function(holder, genomePanel) {
    document.head.appendChild(makeElement('link', '', {rel: 'stylesheet', href: this.uiPrefix + 'css/bootstrap-scoped.css'}));
    document.head.appendChild(makeElement('link', '', {rel: 'stylesheet', href: this.uiPrefix + 'css/dalliance-scoped.css'}));
    document.head.appendChild(makeElement('link', '', {rel: 'stylesheet', href: this.uiPrefix + 'css/font-awesome.min.css'}));

    var b = this;

    if (!b.disableDefaultFeaturePopup) {
        this.addFeatureListener(function(ev, feature, hit, tier) {
            b.featurePopup(ev, feature, hit, tier);
        });
    }

    holder.classList.add('dalliance');
    var toolbar = b.toolbar = makeElement('div', null, {className: 'btn-toolbar toolbar'});

    var title = b.coordSystem.speciesName + ' ' + b.nameForCoordSystem(b.coordSystem);
    if (this.setDocumentTitle) {
        document.title = title + ' :: dalliance';
    }
    
    var locField = makeElement('input', '', {className: 'loc-field'});
    b.makeTooltip(locField, 'Enter a genomic location or gene name');
    var locStatusField = makeElement('p', '', {className: 'loc-status'});


    var zoomInBtn = makeElement('a', [makeElement('i', null, {className: 'fa fa-search-plus'})], {className: 'btn'});
    var zoomSlider = makeElement('input', '', {type: 'range', min: 100, max: 250}, {className: 'zoom-slider'}, {width: '150px'});  // NB min and max get overwritten.
    var zoomOutBtn = makeElement('a', [makeElement('i', null, {className: 'fa fa-search-minus'})], {className: 'btn'});

    var clearHighlightsButton = makeElement('a', [makeElement('i', null, {className: 'fa fa-eraser'})], {className: 'btn'});

    var addTrackBtn = makeElement('a', [makeElement('i', null, {className: 'fa fa-plus'})], {className: 'btn'});
    var favBtn = makeElement('a', [makeElement('i', null, {className: 'fa fa-bookmark'})], {className: 'btn'});
    var svgBtn = makeElement('a', [makeElement('i', null, {className: 'fa fa-print'})], {className: 'btn'});
    var resetBtn = makeElement('a', [makeElement('i', null, {className: 'fa fa-refresh'})], {className: 'btn'});
    var optsButton = makeElement('a', [makeElement('i', null, {className: 'fa fa-cogs'})], {className: 'btn'});
    var helpButton = makeElement('a', [makeElement('i', null, {className: 'fa fa-question'})], {className: 'btn'});

    var tierEditButton = makeElement('a', [makeElement('i', null, {className: 'fa fa-road'})], {className: 'btn'});
    b.makeTooltip(tierEditButton, 'Configure currently selected track(s) (E)')

    var leapLeftButton = makeElement('a', [makeElement('i', null, {className: 'fa fa-chevron-left'})], {className: 'btn'});
    var leapRightButton = makeElement('a', [makeElement('i', null, {className: 'fa fa-chevron-right'})], {className: 'btn pull-right'});

    var modeButtons = makeElement('div', null, {className: 'btn-group pull-right'});
    if (!this.noTrackAdder)
        modeButtons.appendChild(addTrackBtn);
    if (!this.noTrackEditor)
        modeButtons.appendChild(tierEditButton);
    if (!this.noExport)
        modeButtons.appendChild(svgBtn);
    if (!this.noOptions)
        modeButtons.appendChild(optsButton);
    if (!this.noHelp)
        modeButtons.appendChild(helpButton);

    this.setUiMode = function(m) {
        this.uiMode = m;
        var mb = {help: helpButton, add: addTrackBtn, opts: optsButton, 'export': svgBtn, tier: tierEditButton};
        for (var x in mb) {
            if (x == m)
                mb[x].classList.add('active');
            else
                mb[x].classList.remove('active');
        }
    }


    if (!this.noLeapButtons)
        toolbar.appendChild(leapRightButton);

    if (modeButtons.firstChild)
        toolbar.appendChild(modeButtons);
    
    if (!this.noLeapButtons)
        toolbar.appendChild(leapLeftButton);
    if (!this.noTitle) {
        toolbar.appendChild(makeElement('div', makeElement('h4', title, {}, {margin: '0px'}), {className: 'btn-group title'}));
    }
    if (!this.noLocationField)
        toolbar.appendChild(makeElement('div', [locField, locStatusField], {className: 'btn-group loc-group'}));
    if (!this.noClearHighlightsButton)
        toolbar.appendChild(clearHighlightsButton);

    if (!this.noZoomSlider) {
        toolbar.appendChild(makeElement('div', [zoomInBtn,
                                                makeElement('span', zoomSlider, {className: 'btn'}),
                                                zoomOutBtn], {className: 'btn-group'}));
    }
    

    holder.appendChild(toolbar);
    holder.appendChild(genomePanel);

    this.addViewListener(function(chr, min, max, _oldZoom, zoom) {
        locField.value = (chr + ':' + formatLongInt(min) + '..' + formatLongInt(max));
        zoomSlider.min = zoom.min|0;
        zoomSlider.max = zoom.max|0;
        zoomSlider.value = zoom.current|0;
        if (b.storeStatus) {
            b.storeViewStatus();
        }

        if (b.highlights.length > 0) {
            clearHighlightsButton.style.display = 'inline-block';
        } else {
            clearHighlightsButton.style.display = 'none';
        }
    });

    this.addTierListener(function() {
        if (b.storeStatus) {
            b.storeTierStatus();
        }
    });

    locField.addEventListener('keydown', function(ev) {
        if (ev.keyCode == 40) {
            ev.preventDefault(); ev.stopPropagation();
            b.setSelectedTier(0);
        } if (ev.keyCode == 10 || ev.keyCode == 13) {
            ev.preventDefault();


            var g = locField.value;
            b.search(g, function(err) {
                if (err) {
                    locStatusField.textContent = '' + err;
                } else {
                    locStatusField.textContent = '';
                }
            });
        }
    }, false);
    
    var trackAddPopup;
    addTrackBtn.addEventListener('click', function(ev) {
        if (trackAddPopup && trackAddPopup.displayed) {
            b.removeAllPopups();
        } else {
            trackAddPopup = b.showTrackAdder(ev);
        }
    }, false);
    b.makeTooltip(addTrackBtn, 'Add a new track from the registry or an indexed file. (A)');

    zoomInBtn.addEventListener('click', function(ev) {
      ev.stopPropagation(); ev.preventDefault();

      b.zoomStep(-10);
    }, false);
    b.makeTooltip(zoomInBtn, 'Zoom in (+)');

    zoomOutBtn.addEventListener('click', function(ev) {
      ev.stopPropagation(); ev.preventDefault();

      b.zoomStep(10);
    }, false);
    b.makeTooltip(zoomOutBtn, 'Zoom out (-)');

    zoomSlider.addEventListener('change', function(ev) {
    	b.zoomSliderValue = (1.0 * zoomSlider.value);
    	b.zoom(Math.exp((1.0 * zoomSlider.value) / b.zoomExpt));
    }, false);

    favBtn.addEventListener('click', function(ev) {
       ev.stopPropagation(); ev.preventDefault();
    }, false);
    b.makeTooltip(favBtn, 'Favourite regions');

    svgBtn.addEventListener('click', function(ev) {
       ev.stopPropagation(); ev.preventDefault();
        b.openExportPanel();
    }, false);
    b.makeTooltip(svgBtn, 'Export publication-quality SVG. (X)');

    var optsPopup;
    optsButton.addEventListener('click', function(ev) {
        ev.stopPropagation(); ev.preventDefault();

        b.toggleOptsPopup(ev);
    }, false);
    b.makeTooltip(optsButton, 'Configure options.');

    helpButton.addEventListener('click', function(ev) {
        ev.stopPropagation(); ev.preventDefault();
        b.toggleHelpPopup(ev);
    });
    b.makeTooltip(helpButton, 'Help; Keyboard shortcuts. (H)');

    tierEditButton.addEventListener('click', function(ev) {
        ev.stopPropagation(); ev.preventDefault();
        if (b.selectedTiers.length == 1) {
            b.openTierPanel(b.tiers[b.selectedTiers[0]]);
        }
    }, false);

    leapLeftButton.addEventListener('click', function(ev) {
        b.leap(b.reverseKeyScrolling ? -1 : 1, false);
    }, false);
    b.makeTooltip(leapLeftButton, function(ev) {
        var st = b.getSelectedTier();
        var tier;
        if (st >= 0)
            tier = b.tiers[st];

        if (tier && tier.featureSource && sourceAdapterIsCapable(tier.featureSource, 'quantLeap') && typeof(tier.quantLeapThreshold) == 'number') {
            return 'Jump to the next region with a score above the threshold in the selected track "' + (tier.config.name || tier.dasSource.name) + '"" (ctrl+LEFT)';
        } else if (tier && tier.featureSource && sourceAdapterIsCapable(tier.featureSource, 'leap')) {
            return 'Jump to the next feature in the selected track "' + (tier.config.name || tier.dasSource.name) + '" (ctrl+LEFT)';
        } else {
            return 'Jump left (shift+LEFT)';
        }
    });

    leapRightButton.addEventListener('click', function(ev) {
        b.leap(b.reverseKeyScrolling ? 1 : -1, false);
    }, false);
    b.makeTooltip(leapRightButton, function(ev) {
        var st = b.getSelectedTier();
        var tier;
        if (st >= 0)
            tier = b.tiers[st];

        if (tier && tier.featureSource && sourceAdapterIsCapable(tier.featureSource, 'quantLeap') && typeof(tier.quantLeapThreshold) == 'number') {
            return 'Jump to the next region with a score above the threshold in the selected track "' + (tier.config.name || tier.dasSource.name) + '"" (ctrl+RIGHT)';
        } else if (tier && tier.featureSource && sourceAdapterIsCapable(tier.featureSource, 'leap')) {
            return 'Jump to the next feature in the selected track "' + (tier.config.name || tier.dasSource.name) + '" (ctrl+RIGHT)';
        } else {
            return 'Jump right (shift+RIGHT)';
        }
    });

    clearHighlightsButton.addEventListener('click', function(ev) {
        b.clearHighlights();
    }, false);
    b.makeTooltip(clearHighlightsButton, 'Clear highlights (C)');

    b.addTierSelectionWrapListener(function(dir) {
        if (dir < 0) {
            b.setSelectedTier(null);
            locField.focus();
        }
    });

    b.addTierSelectionListener(function(sel) {
        if (b.uiMode === 'tier') {
            if (sel.length == 0) {
                b.hideToolPanel();
                b.manipulatingTier = null;
                b.uiMode = 'none';
            } else {
                var ft = b.tiers[sel[0]];
                if (ft != b.manipulatingTier) {
                    b.openTierPanel(ft);
                }
            }
        }
    });

    var uiKeyHandler = function(ev) {
        // console.log('bukh: ' + ev.keyCode);
        if (ev.keyCode == 65 || ev.keyCode == 97) {  // a
            ev.preventDefault(); ev.stopPropagation();
            b.showTrackAdder();
        } else if (ev.keyCode == 72 || ev.keyCode == 104) { // h
            ev.stopPropagation(); ev.preventDefault();
            b.toggleHelpPopup(ev);
        } else if (ev.keyCode == 69 || ev.keyCode == 101) { //e
            ev.stopPropagation(); ev.preventDefault();
            if (b.selectedTiers.length == 1) {
                b.openTierPanel(b.tiers[b.selectedTiers[0]]);
            }
        } else if (ev.keyCode == 88 || ev.keyCode == 120) { // x
            ev.stopPropagation(); ev.preventDefault();
            b.openExportPanel();
        } else if (ev.keyCode == 67 || ev.keyCode == 99) { // c
            ev.stopPropagation(); ev.preventDefault();
            b.clearHighlights();
        }
    };

    holder.addEventListener('focus', function(ev) {
        holder.addEventListener('keydown', uiKeyHandler, false);
    }, false);
    holder.addEventListener('blur', function(ev) {
        holder.removeEventListener('keydown', uiKeyHandler, false);
    }, false);

    holder.addEventListener('keydown', function(ev) {
        if (ev.keyCode === 27) {
            ev.preventDefault();
            ev.stopPropagation();
            b.uiMode = 'none';
            b.hideToolPanel();

            if (b.selectedTiers && b.selectedTiers.length > 0) {
                b.browserHolder.focus();
            }
        }
    }, false);
}

Browser.prototype.showToolPanel = function(panel, nowrap) {
    var thisB = this;

    if (this.activeToolPanel) {
        this.activeToolPanel.parentElement.removeChild(this.activeToolPanel);
    }

    var content;
    if (nowrap)
        content = panel;
    else
        content = makeElement('div', panel, {}, {overflowY: 'auto', width: '100%'});


    var divider = makeElement('div', makeElement('i', null, {className: 'fa fa-caret-right'}), {className: 'tool-divider'});
    divider.addEventListener('click', function(ev) {
        thisB.hideToolPanel();
        thisB.setUiMode('none');
    }, false);
    this.makeTooltip(divider, 'Close tool panel (ESC)');
    this.activeToolPanel = makeElement('div', [divider, content], {className: 'tool-holder'});
    this.svgHolder.appendChild(this.activeToolPanel);
    this.resizeViewer();

    var thisB = this;
}

Browser.prototype.hideToolPanel = function() {
    if (this.activeToolPanel) {
        this.activeToolPanel.parentElement.removeChild(this.activeToolPanel);
    }
    this.svgHolder.style.width = '100%';
    this.activeToolPanel = null;
    this.resizeViewer();
}

Browser.prototype.toggleHelpPopup = function(ev) {
    if (this.uiMode === 'help') {
        this.hideToolPanel();
        this.setUiMode('none');
    } else {
        var helpFrame = makeElement('iframe', null, {scrolling: 'yes', seamless: 'seamless', src: this.uiPrefix + 'help/index.html', className: 'help-panel'});
        this.showToolPanel(helpFrame, false);
        this.setUiMode('help');
    }
}

Browser.prototype.toggleOptsPopup = function(ev) {
    var b = this;

    if (this.uiMode === 'opts') {
        this.hideToolPanel();
        this.setUiMode('none');
    } else {
        var optsForm = makeElement('div', null, {className: 'form-horizontal'}, {boxSizing: 'border-box', MozBoxSizing: 'border-box', display: 'inline-block', verticalAlign: 'top'});
        var optsTable = makeElement('table');
        optsTable.cellPadding = 5;

        var scrollModeButton = makeElement('input', '', {type: 'checkbox', checked: b.reverseScrolling});
        scrollModeButton.addEventListener('change', function(ev) {
            b.reverseScrolling = scrollModeButton.checked;
            b.storeStatus();
        }, false);
        optsTable.appendChild(makeElement('tr', [makeElement('td', 'Reverse trackpad scrolling', {align: 'right'}), makeElement('td', scrollModeButton)]));

        var scrollKeyButton = makeElement('input', '', {type: 'checkbox', checked: b.reverseKeyScrolling});
        scrollKeyButton.addEventListener('change', function(ev) {
            b.reverseKeyScrolling = scrollKeyButton.checked;
            b.storeStatus();
        }, false);
        optsTable.appendChild(makeElement('tr', [makeElement('td', 'Reverse scrolling buttons and keys', {align: 'right'}), makeElement('td', scrollKeyButton)]));


        var rulerSelect = makeElement('select');
        rulerSelect.appendChild(makeElement('option', 'Left', {value: 'left'}));
        rulerSelect.appendChild(makeElement('option', 'Center', {value: 'center'}));
        rulerSelect.appendChild(makeElement('option', 'Right', {value: 'right'}));
        rulerSelect.appendChild(makeElement('option', 'None', {value: 'none'}));
        rulerSelect.value = b.rulerLocation;
        rulerSelect.addEventListener('change', function(ev) {
            b.rulerLocation = rulerSelect.value;
            b.positionRuler();
            for (var ti = 0; ti < b.tiers.length; ++ti) {
                b.tiers[ti].paintQuant();
            }
            b.storeStatus();
        }, false);
        optsTable.appendChild(makeElement('tr', [makeElement('td', 'Vertical guideline', {align: 'right'}), makeElement('td', rulerSelect)]));

        optsForm.appendChild(optsTable);

        var resetButton = makeElement('button', 'Reset browser', {className: 'btn'}, {marginLeft: 'auto', marginRight: 'auto', display: 'block'});
        resetButton.addEventListener('click', function(ev) {
            b.reset();
        }, false);
        optsForm.appendChild(resetButton);

        this.showToolPanel(optsForm);
        this.setUiMode('opts');
    }
}// 
// Dalliance Genome Explorer
// (c) Thomas Down 2006-2010
//
// glyphs.js: components which know how to draw themselves
//

"use strict";

function SVGPath() {
    this.ops = [];
}

SVGPath.prototype.moveTo = function(x, y) {
    this.ops.push('M ' + x + ' ' + y);
}

SVGPath.prototype.lineTo = function(x, y) {
    this.ops.push('L ' + x + ' ' + y);
}

SVGPath.prototype.closePath = function() {
    this.ops.push('Z');
}

SVGPath.prototype.toPathData = function() {
    return this.ops.join(' ');
}

function PathGlyphBase(stroke, fill) {
    this._stroke = stroke;
    this._fill = fill;
}

PathGlyphBase.prototype.draw = function(g) {
    g.beginPath();
    this.drawPath(g);

    if (this._fill) {
        g.fillStyle = this._fill;
        g.fill();
    }
    if (this._stroke) {
        g.strokeStyle = this._stroke;
        g.stroke();
    }
}

PathGlyphBase.prototype.toSVG = function() {
    var g = new SVGPath();
    this.drawPath(g);
    
    return makeElementNS(
        NS_SVG, 'path',
        null,
        {d: g.toPathData(),
         fill: this._fill || 'none',
         stroke: this._stroke || 'none'});
}

PathGlyphBase.prototype.drawPath = function(g) {
    throw 'drawPath method on PathGlyphBase must be overridden';
}

function BoxGlyph(x, y, width, height, fill, stroke, alpha, radius) {
    this.x = x;
    this.y = y;
    this._width = width;
    this._height = height;
    this.fill = fill;
    this.stroke = stroke;
    this._alpha = alpha;
    this._radius = radius || 0;
}

BoxGlyph.prototype.draw = function(g) {
    var r = this._radius;

    if (r > 0) {
        g.beginPath();
        g.moveTo(this.x + r, this.y);
        g.lineTo(this.x + this._width - r, this.y);
        g.arcTo(this.x + this._width, this.y, this.x + this._width, this.y + r, r);
        g.lineTo(this.x + this._width, this.y + this._height - r);
        g.arcTo(this.x + this._width, this.y + this._height, this.x + this._width - r, this.y + this._height, r);
        g.lineTo(this.x + r, this.y + this._height);
        g.arcTo(this.x, this.y + this._height, this.x, this.y + this._height - r, r);
        g.lineTo(this.x, this.y + r);
        g.arcTo(this.x, this.y, this.x + r, this.y, r);
        g.closePath();

        if (this._alpha != null) {
            g.save();
            g.globalAlpha = this._alpha;
        }
        
        if (this.fill) {
            g.fillStyle = this.fill;
            g.fill();
        }
        if (this.stroke) {
            g.strokeStyle = this.stroke;
            g.lineWidth = 0.5;
            g.stroke();
        }

        if (this._alpha != null) {
            g.restore();
        }
    } else {
        if (this._alpha != null) {
            g.save();
            g.globalAlpha = this._alpha;
        }

        if (this.fill) {
            g.fillStyle = this.fill;
            g.fillRect(this.x, this.y, this._width, this._height);
        }

        if (this.stroke) {
            g.strokeStyle = this.stroke;
            g.lineWidth = 0.5;
            g.strokeRect(this.x, this.y, this._width, this._height)
        }

        if (this._alpha != null) {
            g.restore();
        }
    }
}

BoxGlyph.prototype.toSVG = function() {
    var s = makeElementNS(NS_SVG, 'rect', null,
                         {x: this.x, 
                          y: this.y, 
                          width: this._width, 
                          height: this._height,
                          stroke: this.stroke || 'none',
                          strokeWidth: 0.5,
                          fill: this.fill || 'none'});
    if (this._alpha != null) {
        s.setAttribute('opacity', this._alpha);
    }

    return s;
}

BoxGlyph.prototype.min = function() {
    return this.x;
}

BoxGlyph.prototype.max = function() {
    return this.x + this._width;
}

BoxGlyph.prototype.height = function() {
    return this.y + this._height;
}


function GroupGlyph(glyphs, connector) {
    this.glyphs = glyphs;
    this.connector = connector;
    this.h = glyphs[0].height();

    var covList = [];
    for (var g = 0; g < glyphs.length; ++g) {
        var gg = glyphs[g];
        covList.push(new Range(gg.min(), gg.max()));
        this.h = Math.max(this.h, gg.height());
    }
    this.coverage = union(covList);
}

GroupGlyph.prototype.drawConnectors = function(g) {
    var ranges = this.coverage.ranges();
    for (var r = 1; r < ranges.length; ++r) {
        var gl = ranges[r];
        var last = ranges[r - 1];
        if (last && gl.min() > last.max()) {
            var start = last.max();
            var end = gl.min();
            var mid = (start+end)/2
            
            if (this.connector === 'hat+') {
                g.moveTo(start, this.h/2);
                g.lineTo(mid, 0);
                g.lineTo(end, this.h/2);
            } else if (this.connector === 'hat-') {
                g.moveTo(start, this.h/2);
                g.lineTo(mid, this.h);
                g.lineTo(end, this.h/2);
            } else if (this.connector === 'collapsed+') {
                g.moveTo(start, this.h/2);
                g.lineTo(end, this.h/2);
                if (end - start > 4) {
                    g.moveTo(mid - 2, (this.h/2) - 3);
                    g.lineTo(mid + 2, this.h/2);
                    g.lineTo(mid - 2, (this.h/2) + 3);
                }
            } else if (this.connector === 'collapsed-') {
                g.moveTo(start, this.h/2);
                g.lineTo(end, this.h/2);
                if (end - start > 4) {
                    g.moveTo(mid + 2, (this.h/2) - 3);
                    g.lineTo(mid - 2, this.h/2);
                    g.lineTo(mid + 2, (this.h/2) + 3);
                }
            } else {
                g.moveTo(start, this.h/2);
                g.lineTo(end, this.h/2);
            }
        }
        last = gl;
    }
}

GroupGlyph.prototype.draw = function(g) {
    for (var i = 0; i < this.glyphs.length; ++i) {
        var gl = this.glyphs[i];
        gl.draw(g);
    }

    g.strokeStyle = 'black';
    g.beginPath();
    this.drawConnectors(g);
    g.stroke();
}

GroupGlyph.prototype.toSVG = function() {
    var g = makeElementNS(NS_SVG, 'g');
    for (var i = 0; i < this.glyphs.length; ++i) {
        g.appendChild(this.glyphs[i].toSVG());
    }

    var p = new SVGPath();
    this.drawConnectors(p);

    var pathData = p.toPathData();
    if (pathData.length > 0) {
        var path = makeElementNS(
            NS_SVG, 'path',
            null,
            {d: p.toPathData(),
             fill: 'none',
             stroke: 'black',
             strokeWidth: 0.5});
        g.appendChild(path);
    }

    return g;
}

GroupGlyph.prototype.min = function() {
    return this.coverage.min();
}

GroupGlyph.prototype.max = function() {
    return this.coverage.max();
}

GroupGlyph.prototype.height = function() {
    return this.h;
}


function LineGraphGlyph(points, color, height) {
    this.points = points;
    this.color = color;
    this._height = height || 50;
}

LineGraphGlyph.prototype.min = function() {
    return this.points[0];
};

LineGraphGlyph.prototype.max = function() {
    return this.points[this.points.length - 2];
};

LineGraphGlyph.prototype.height = function() {
    return this._height;
}

LineGraphGlyph.prototype.draw = function(g) {
    g.save();
    g.strokeStyle = this.color;
    g.lineWidth = 2;
    g.beginPath();
    for (var i = 0; i < this.points.length; i += 2) {
        var x = this.points[i];
        var y = this.points[i + 1];
        if (i == 0) {
            g.moveTo(x, y);
        } else {
            g.lineTo(x, y);
        }
    }
    g.stroke();
    g.restore();
}

LineGraphGlyph.prototype.toSVG = function() {
    var p = new SVGPath();
    for (var i = 0; i < this.points.length; i += 2) {
        var x = this.points[i];
        var y = this.points[i + 1];
        if (i == 0) {
            p.moveTo(x, y);
        } else {
            p.lineTo(x, y);
        }
    }
    
    return makeElementNS(
        NS_SVG, 'path',
        null,
        {d: p.toPathData(),
         fill: 'none',
         stroke: this.color,
         strokeWidth: '2px'});
}

function LabelledGlyph(glyph, text, unmeasured, anchor, align, font) {
    this.glyph = glyph;
    this.text = text;
    this.anchor = anchor || 'left';
    this.align = align || 'below';
    if (font) {
        this.font = font;
    }
    if (this.font) {
        GLOBAL_GC.save();
        GLOBAL_GC.font = this.font;
    }
    var metrics = GLOBAL_GC.measureText(text);
    if (this.font) {
        GLOBAL_GC.restore();
    }
    this.textLen = metrics.width;
    this.textHeight = 5;
    this.bump = glyph.bump;
    this.measured = !unmeasured;
}

LabelledGlyph.prototype.toSVG = function() {
    var child = this.glyph.toSVG();
    var opts = {};
    
    if (this.align == 'above') {
        child = makeElementNS(NS_SVG, 'g', child, {transform: "translate(0, " + (this.textHeight|0 + 2) + ")"});
        opts.y = this.textHeight;
    } else {
        opts.y = this.glyph.height() + 15;
    }

    if (this.font) {
        opts.fontSize  = 7;
    }

    if ('center' == this.anchor) {
        opts.x = (this.glyph.min() + this.glyph.max() - this.textLen) / 2;
    } else {
        opts.x = this.glyph.min();
    }

    return makeElementNS(NS_SVG, 'g',
        [child,
         makeElementNS(NS_SVG, 'text', this.text, opts)]);
}

LabelledGlyph.prototype.min = function() {
    return this.glyph.min();
}

LabelledGlyph.prototype.max = function() {
    if (this.measured)
        return Math.max(this.glyph.max(), (1.0*this.glyph.min()) + this.textLen + 10);
    else
        return this.glyph.max();
}

LabelledGlyph.prototype.height = function() {
    var h = this.glyph.height();
    if (this.measured) {
        if (this.align == 'above') {
            h += this.textHeight + 2;
        } else {
            h += 20;
        }
    }
    return h;
}

LabelledGlyph.prototype.draw = function(g) {
    if (this.align == 'above') {
        g.save();
        g.translate(0, this.textHeight + 2);
    }
    this.glyph.draw(g);
    if (this.align == 'above') {
        g.restore();
    }

    g.fillStyle = 'black';
    if (this.font) {
        g.save();
        g.font = this.font;
    }
    var p;
    if ('center' == this.anchor) {
        p = (this.glyph.min() + this.glyph.max() - this.textLen) / 2;
    } else {
        p = this.glyph.min();
    }
    g.fillText(this.text, p, this.align == 'above' ? this.textHeight : this.glyph.height() + 15);
    if (this.font) {
        g.restore();
    }
}



function CrossGlyph(x, height, stroke) {
    this._x = x;
    this._height = height;
    this._stroke = stroke;
}

CrossGlyph.prototype.draw = function(g) {
    var hh = this._height/2;
    
    g.beginPath();
    g.moveTo(this._x, 0);
    g.lineTo(this._x, this._height);
    g.moveTo(this._x - hh, hh);
    g.lineTo(this._x + hh, hh);

    g.strokeStyle = this._stroke;
    g.lineWidth = 1;

    g.stroke();
}

CrossGlyph.prototype.toSVG = function() {
    var hh = this._height/2;

    var g = new SVGPath();
    g.moveTo(this._x, 0);
    g.lineTo(this._x, this._height);
    g.moveTo(this._x - hh, hh);
    g.lineTo(this._x + hh, hh);
    
    return makeElementNS(
        NS_SVG, 'path',
        null,
        {d: g.toPathData(),
         fill: 'none',
         stroke: this._stroke,
         strokeWidth: '1px'});
}

CrossGlyph.prototype.min = function() {
    return this._x - this._height/2;
}

CrossGlyph.prototype.max = function() {
    return this._x + this._height/2;
}

CrossGlyph.prototype.height = function() {
    return this._height;
}

function ExGlyph(x, height, stroke) {
    this._x = x;
    this._height = height;
    this._stroke = stroke;
}

ExGlyph.prototype.draw = function(g) {
    var hh = this._height/2;
    
    g.beginPath();
    g.moveTo(this._x - hh, 0);
    g.lineTo(this._x + hh, this._height);
    g.moveTo(this._x - hh, this._height);
    g.lineTo(this._x + hh, 0);

    g.strokeStyle = this._stroke;
    g.lineWidth = 1;

    g.stroke();
}

ExGlyph.prototype.toSVG = function() {
    var hh = this._height/2;

    var g = new SVGPath();
    g.moveTo(this._x - hh, 0);
    g.lineTo(this._x + hh, this._height);
    g.moveTo(this._x - hh, this._height);
    g.lineTo(this._x + hh, 0);
    
    return makeElementNS(
        NS_SVG, 'path',
        null,
        {d: g.toPathData(),
         fill: 'none',
         stroke: this._stroke,
         strokeWidth: '1px'});
}

ExGlyph.prototype.min = function() {
    return this._x - this._height/2;
}

ExGlyph.prototype.max = function() {
    return this._x + this._height/2;
}

ExGlyph.prototype.height = function() {
    return this._height;
}



function TriangleGlyph(x, height, dir, width, fill, stroke) {
    PathGlyphBase.call(this, stroke, fill);

    this._x = x;
    this._height = height;
    this._dir = dir;
    this._width = width;
}

TriangleGlyph.prototype = Object.create(PathGlyphBase.prototype);

TriangleGlyph.prototype.drawPath = function(g) {
    var hh = this._height/2;
    var hw = this._width/2;

    if (this._dir === 'S') {
        g.moveTo(this._x, this._height);
        g.lineTo(this._x - hw, 0);
        g.lineTo(this._x + hw, 0);
    } else if (this._dir === 'W') {
        g.moveTo(this._x + hw, hh);
        g.lineTo(this._x - hw, 0);
        g.lineTo(this._x - hw, this._height);
    } else if (this._dir === 'E') {
        g.moveTo(this._x - hw, hh);
        g.lineTo(this._x + hw, 0);
        g.lineTo(this._x + hw, this._height);
    } else {
        g.moveTo(this._x , 0);
        g.lineTo(this._x + hw, this._height);
        g.lineTo(this._x - hw, this._height);
    }

    g.closePath();
}

TriangleGlyph.prototype.min = function() {
    return this._x - this._height/2;
}

TriangleGlyph.prototype.max = function() {
    return this._x + this._height/2;
}

TriangleGlyph.prototype.height = function() {
    return this._height;
}




function DotGlyph(x, height, fill, stroke) {
    this._x = x;
    this._height = height;
    this._fill = fill;
    this._stroke = stroke;
}

DotGlyph.prototype.draw = function(g) {
    var hh = this._height/2;
    g.fillStyle = this._stroke;
    g.beginPath();
    g.arc(this._x, hh, hh, 0, 6.29);

    if (this._fill) {
        g.fillStyle = this._fill;
        g.fill();
    }

    if (this._stroke) {
        g.strokeStyle = this._stroke;
        g.stroke();
    }
}

DotGlyph.prototype.toSVG = function() {
    var hh = this._height/2;
    return makeElementNS(
        NS_SVG, 'circle',
        null,
        {cx: this._x, cy: hh, r: hh,
         fill: this._fill || 'none',
         stroke: this._stroke || 'none',
         strokeWidth: '1px'});
}

DotGlyph.prototype.min = function() {
    return this._x - this._height/2;
}

DotGlyph.prototype.max = function() {
    return this._x + this._height/2;
}

DotGlyph.prototype.height = function() {
    return this._height;
}


function PaddedGlyph(glyph, minp, maxp) {
    this.glyph = glyph;
    this._min = minp;
    this._max = maxp;
    if (glyph) {
        this.bump = glyph.bump;
    }
}

PaddedGlyph.prototype.draw = function(g) {
    if (this.glyph) 
        this.glyph.draw(g);
}

PaddedGlyph.prototype.toSVG = function() {
    if (this.glyph) {
        return this.glyph.toSVG();
    } else {
        return makeElementNS(NS_SVG, 'g');
    }
}

PaddedGlyph.prototype.min = function() {
    return this._min;
}

PaddedGlyph.prototype.max = function() {
    return this._max;
}

PaddedGlyph.prototype.height = function() {
    if (this.glyph) {
        return this.glyph.height();
    } else {
        return 1;
    }
}


function AArrowGlyph(min, max, height, fill, stroke, ori) {
    PathGlyphBase.call(this, stroke, fill);
    this._min = min;
    this._max = max;
    this._height = height;
    this._ori = ori;
}

AArrowGlyph.prototype = Object.create(PathGlyphBase.prototype);

AArrowGlyph.prototype.min = function() {
    return this._min;
}

AArrowGlyph.prototype.max = function() {
    return this._max;
}

AArrowGlyph.prototype.height = function() {
    return this._height;
}

AArrowGlyph.prototype.drawPath = function(g) {
    var maxPos = this._max;
    var minPos = this._min;
    var height = this._height;
    var lInset = 0;
    var rInset = 0;
    var minLength = this._height + 2;
    var instep = 0.333333 * this._height;
    var y = 0;

    if (this._ori) {
        if (this._ori === '+') {
            rInset = 0.5 * this._height;
        } else if (this._ori === '-') {
            lInset = 0.5 * this._height;
        }
    }

    if (maxPos - minPos < minLength) {
        minPos = (maxPos + minPos - minLength) / 2;
        maxPos = minPos + minLength;
    }

    g.moveTo(minPos + lInset, y+instep);
    g.lineTo(maxPos - rInset, y+instep);
    g.lineTo(maxPos - rInset, y);
    g.lineTo(maxPos, y + this._height/2);
    g.lineTo(maxPos - rInset, y+height);
    g.lineTo(maxPos - rInset, y+instep+instep);
    g.lineTo(minPos + lInset, y+instep+instep);
    g.lineTo(minPos + lInset, y+height);
    g.lineTo(minPos, y+height/2);
    g.lineTo(minPos + lInset, y);
    g.lineTo(minPos + lInset, y+instep);
}

function SpanGlyph(min, max, height, stroke) {
    PathGlyphBase.call(this, stroke, null);
    this._min = min;
    this._max = max;
    this._height = height;
}

SpanGlyph.prototype = Object.create(PathGlyphBase.prototype);

SpanGlyph.prototype.min = function() {return this._min};
SpanGlyph.prototype.max = function() {return this._max};
SpanGlyph.prototype.height = function() {return this._height};

SpanGlyph.prototype.drawPath = function(g) {
    var minPos = this._min, maxPos = this._max;
    var height = this._height, hh = height/2;
    g.moveTo(minPos, hh);
    g.lineTo(maxPos, hh);
    g.moveTo(minPos, 0);
    g.lineTo(minPos, height);
    g.moveTo(maxPos, 0);
    g.lineTo(maxPos, height);
}

function LineGlyph(min, max, height, style, strand, stroke) {
    this._min = min;
    this._max = max;
    this._height = height;
    this._style = style;
    this._strand = strand;
    this._stroke = stroke;
}

LineGlyph.prototype.min = function() {return this._min};
LineGlyph.prototype.max = function() {return this._max};
LineGlyph.prototype.height = function() {return this._height};

LineGlyph.prototype.drawPath = function(g) {
    var minPos = this._min, maxPos = this._max;
    var height = this._height, hh = height/2;

    if (this._style === 'hat') {
        g.moveTo(minPos, hh);
        g.lineTo((minPos + maxPos)/2, this._strand === '-' ? height : 0);
        g.lineTo(maxPos, hh);
    } else {
        g.moveTo(minPos, hh);
        g.lineTo(maxPos, hh);
    }
}


LineGlyph.prototype.draw = function(g) {
    g.beginPath();
    this.drawPath(g);
    g.strokeStyle = this._stroke;
    if (this._style === 'dashed' && g.setLineDash) {
        g.save();
        g.setLineDash([3]);
        g.stroke();
        g.restore();
    } else {
        g.stroke();
    }
}

LineGlyph.prototype.toSVG = function() {
    var g = new SVGPath();
    this.drawPath(g);
    
    var opts = {d: g.toPathData(),
            stroke: this._stroke || 'none'};
    if (this._style === 'dashed') {
        opts['strokeDasharray'] = '3';
    }

    return makeElementNS(
        NS_SVG, 'path',
        null, opts
    );
}





function PrimersGlyph(min, max, height, fill, stroke) {
    this._min = min;
    this._max = max;
    this._height = height;
    this._fill = fill;
    this._stroke = stroke;
}

PrimersGlyph.prototype.min = function() {return this._min};
PrimersGlyph.prototype.max = function() {return this._max};
PrimersGlyph.prototype.height = function() {return this._height};


PrimersGlyph.prototype.drawStemPath = function(g) {
    var minPos = this._min, maxPos = this._max;
    var height = this._height, hh = height/2;
    g.moveTo(minPos, hh);
    g.lineTo(maxPos, hh);
}

PrimersGlyph.prototype.drawTrigsPath = function(g) {
    var minPos = this._min, maxPos = this._max;
    var height = this._height, hh = height/2;
    g.moveTo(minPos, 0);
    g.lineTo(minPos + height, hh);
    g.lineTo(minPos, height);
    g.lineTo(minPos, 0);
    g.moveTo(maxPos, 0);
    g.lineTo(maxPos - height, hh);
    g.lineTo(maxPos, height);
    g.lineTo(maxPos, 0);
}


PrimersGlyph.prototype.draw = function(g) {
    g.beginPath();
    this.drawStemPath(g);
    g.strokeStyle = this._stroke;
    g.stroke();
    g.beginPath();
    this.drawTrigsPath(g);
    g.fillStyle = this._fill;
    g.fill();
}

PrimersGlyph.prototype.toSVG = function() {
    var s = new SVGPath();
    this.drawStemPath(s);
    var t = new SVGPath();
    this.drawTrigsPath(t);
    
    return makeElementNS(
        NS_SVG, 'g',
        [makeElementNS(
            NS_SVG, 'path',
            null,
            {d: s.toPathData(),
             stroke: this._stroke || 'none'}),
         makeElementNS(
             NS_SVG, 'path',
             null,
             {d: t.toPathData(),
              fill: this._fill || 'none'})]);
}

function ArrowGlyph(min, max, height, color, parallel, sw, ne) {
    PathGlyphBase.call(this, null, color);
    this._min = min;
    this._max = max;
    this._height = height;
    this._color = color;
    this._parallel = parallel;
    this._sw = sw;
    this._ne = ne;
}

ArrowGlyph.prototype = Object.create(PathGlyphBase.prototype);

ArrowGlyph.prototype.min = function() {return this._min};
ArrowGlyph.prototype.max = function() {return this._max};
ArrowGlyph.prototype.height = function() {return this._height};

ArrowGlyph.prototype.drawPath = function(g) {
    var min = this._min, max = this._max, height = this._height;
    
    if (this._parallel) {
        var hh = height/2;
        var instep = 0.4 * height;
        if (this._sw) {
            g.moveTo(min + hh, height-instep);
            g.lineTo(min + hh, height);
            g.lineTo(min, hh);
            g.lineTo(min + hh, 0);
            g.lineTo(min + hh, instep);
        } else {
            g.moveTo(min, height-instep);
            g.lineTo(min, instep);
        }
        if (this._ne) {
            g.lineTo(max - hh, instep);
            g.lineTo(max - hh, 0);
            g.lineTo(max, hh);
            g.lineTo(max - hh, height);
            g.lineTo(max - hh, height - instep);
        } else {
            g.lineTo(max, instep);
            g.lineTo(max, height-instep);
        }
        g.closePath();
    } else {
        var mid = (min+max)/2;
        var instep = 0.4*(max-min);
        var th = height/3;

        if (this._ne) {
            g.moveTo(min + instep, th);
            g.lineTo(min, th);
            g.lineTo(mid, 0);
            g.lineTo(max, th);
            g.lineTo(max - instep, th);
        } else {
            g.moveTo(min+instep, 0);
            g.lineTo(max-instep, 0);
        }
        if (this._sw) {
            g.lineTo(max - instep, height-th);
            g.lineTo(max, height-th);
            g.lineTo(mid, height);
            g.lineTo(min, height-th)
            g.lineTo(min + instep, height-th);
        } else {
            g.lineTo(max - instep, height);
            g.lineTo(min + instep, height);
        }
        g.closePath();
    }
}


function TooManyGlyph(min, max, height, fill, stroke) {
    this._min = min;
    this._max = max;
    this._height = height;
    this._fill = fill;
    this._stroke = stroke;
}

TooManyGlyph.prototype.min = function() {return this._min};
TooManyGlyph.prototype.max = function() {return this._max};
TooManyGlyph.prototype.height = function() {return this._height};

TooManyGlyph.prototype.toSVG = function() {
    return makeElementNS(NS_SVG, 'rect', null,
                         {x: this._min, 
                          y: 0, 
                          width: this._max - this._min, 
                          height: this._height,
                          stroke: this._stroke || 'none',
                          fill: this._fill || 'none'});
}

TooManyGlyph.prototype.draw = function(g) {
    if (this._fill) {
        g.fillStyle = this._fill;
        g.fillRect(this._min, 0, this._max - this._min, this._height);
    }
    if (this._stroke) {
        g.strokeStyle = this._stroke;
        g.strokeRect(this._min, 0, this._max - this._min, this._height);
        g.beginPath();
        for (var n = 2; n < this._height; n += 3) {
            g.moveTo(this._min, n);
            g.lineTo(this._max, n);
        }
        g.stroke();
    }
}

function TextGlyph(min, max, height, fill, string) {
    this._min = min;
    this._max = max;
    this._height = height;
    this._fill = fill;
    this._string = string;
    this._textLen = GLOBAL_GC.measureText(string).width;
}

TextGlyph.prototype.min = function() {return this._min};
TextGlyph.prototype.max = function() {return Math.max(this._max, this._min + this._textLen)};
TextGlyph.prototype.height = function() {return this._height};

TextGlyph.prototype.draw = function(g) {
    g.fillStyle = this._fill;
    g.fillText(this._string, this._min, this._height - 4);
}

TextGlyph.prototype.toSVG = function() {
    return makeElementNS(NS_SVG, 'text', this._string, {x: this._min, y: this._height - 4});
};

(function(scope) {

var isRetina = window.devicePixelRatio > 1;
var __dalliance_SequenceGlyphCache = {};

function SequenceGlyph(min, max, height, seq, ref, scheme, quals) {
    this._min = min;
    this._max = max;
    this._height = height;
    this._seq = seq;
    this._ref = ref;
    this._scheme = scheme;
    this._quals = quals;
}

SequenceGlyph.prototype.min = function() {return this._min};
SequenceGlyph.prototype.max = function() {return this._max};
SequenceGlyph.prototype.height = function() {return this._height};


SequenceGlyph.prototype.alphaForQual = function(qual) {
    return 0.1 + 0.9*Math.max(0.0, Math.min((1.0 * qual) / 30.0, 1.0));
}

SequenceGlyph.prototype.draw = function(gc) {
    var seq = this._seq;
    if (!seq)
        seq = 'NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN';    // FIXME       

    var scale = (this._max - this._min + 1) / seq.length;

    for (var p = 0; p < seq.length; ++p) {
        var base = seq.substr(p, 1).toUpperCase();
        var color = baseColors[base];
        if (!color) {
            color = 'gray';
        }

        if (this._scheme === 'mismatch' && this._ref) {
            var refbase = this._ref.substr(p, 1).toUpperCase();
            if (refbase === 'N') {
                color = 'gray';
            } else if (refbase === base) {
                color = 'black';
            } else {
                color = 'red';
            }
        }
        
        if (this._quals) {
            var qc = this._quals.charCodeAt(p) - 33;
            var oldAlpha = gc.globalAlpha;            // NB hoisted!
            gc.globalAlpha = this.alphaForQual(qc);
        }

        if (scale >= 8) {
            var key = color + '_' + base
            var img = __dalliance_SequenceGlyphCache[key];
            if (!img) {
                img = document.createElement('canvas');
                if (isRetina) {
                    img.width = 16;
                    img.height = 20;
                } else {
                    img.width = 8;
                    img.height = 10;
                }
                var imgGc = img.getContext('2d');
                if (isRetina) {
                    imgGc.scale(2, 2);
                }
                imgGc.fillStyle = color;
                imgGc.fillText(base, 0, 8);
                __dalliance_SequenceGlyphCache[key] = img;
            }
            if (isRetina)
                gc.drawImage(img, this._min + p*scale, 0, 8, 10);
            else
                gc.drawImage(img, this._min + p*scale, 0);
        } else {
            gc.fillStyle = color;
            gc.fillRect(this._min + p*scale, 0, scale, this._height);
        }

        if (this._quals) {
            gc.globalAlpha = oldAlpha;
        }
    }
}

SequenceGlyph.prototype.toSVG = function() {
    var seq = this._seq;
    var scale = (this._max - this._min + 1) / this._seq.length;
    var  g = makeElementNS(NS_SVG, 'g'); 

    for (var p = 0; p < seq.length; ++p) {
        var base = seq.substr(p, 1).toUpperCase();
        var color = baseColors[base];
        if (!color) {
            color = 'gray';
        }

        if (this._scheme === 'mismatch' && this._ref) {
            var refbase = this._ref.substr(p, 1).toUpperCase();
            if (refbase === 'N') {
                color = 'gray';
            } else if (refbase === base) {
                color = 'black';
            } else {
                color = 'red';
            }
        }

        var alpha = 1.0;
        if (this._quals) {
            var qc = this._quals.charCodeAt(p) - 33;
            alpha = this.alphaForQual(qc);
        }

        if (scale >= 8) {
            g.appendChild(
                    makeElementNS(NS_SVG, 'text', base, {
                        x: this._min + p*scale,
                        y: 8,
                        fill: color,
                        fillOpacity: alpha}));
        } else {
            g.appendChild(
                    makeElementNS(NS_SVG, 'rect', null, {
                        x:this._min + p*scale,
                        y: 0,
                        width: scale,
                        height: this._height,
                        fill: color,
                        fillOpacity: alpha}));

        }
    }

    return g;
}

scope.SequenceGlyph = SequenceGlyph;

}(this));

function TranslatedGlyph(glyph, x, y, height) {
    this.glyph = glyph;
    this._height = height;
    this._x = x;
    this._y = y;
}

TranslatedGlyph.prototype.height = function() {
    if (this._height) {
        return this._height;
    } else {
        return this.glyph.height() + this._y;
    }
}

TranslatedGlyph.prototype.min = function() {
    return this.glyph.min() + this._x;
}

TranslatedGlyph.prototype.max = function() {
    return this.glyph.max() + this._x;
}

TranslatedGlyph.prototype.minY = function() {
    return this._y;
}

TranslatedGlyph.prototype.maxY = function() {
    return this._y + this.glyph.height();
}

TranslatedGlyph.prototype.draw = function(g) {
    g.save();
    g.translate(this._x, this._y);
    this.glyph.draw(g);
    g.restore();
}

TranslatedGlyph.prototype.toSVG = function() {
    var s =  this.glyph.toSVG();
    s.setAttribute('transform', 'translate(' + this._x + ',' + this._y + ')');
    return s;
}

function PointGlyph(x, y, height, fill) {
    this._x = x;
    this._y = y;
    this._height = height;
    this._fill = fill;
}

PointGlyph.prototype.min = function() {
    return this._x - 2;
}

PointGlyph.prototype.max = function() {
    return this._x + 2;
}

PointGlyph.prototype.height = function() {
    return this._height;
}

PointGlyph.prototype.draw = function(g) {
    g.save();
    g.globalAlpha = 0.3;
    g.fillStyle = this._fill;
    g.beginPath();
    g.arc(this._x, this._y, 1.5, 0, 6.29);
    g.fill();
    g.restore();
}

PointGlyph.prototype.toSVG = function() {
    return makeElementNS(
        NS_SVG, 'circle',
        null,
        {cx: this._x, cy: this._y, r: 2,
         fill: this._fill,
         stroke: 'none'});
}


function GridGlyph(height) {
    this._height = height || 50;
}

GridGlyph.prototype.notSelectable = true;

GridGlyph.prototype.min = function() {
    return -100000;
};

GridGlyph.prototype.max = function() {
    return 100000;
};

GridGlyph.prototype.height = function() {
    return this._height;
}

GridGlyph.prototype.draw = function(g) {
    g.save();
    g.strokeStyle = 'black'
    g.lineWidth = 0.1;

    g.beginPath();
    for (var y = 0; y <= this._height; y += 10) {
        g.moveTo(-5000, y);
        g.lineTo(5000, y);
    }
    g.stroke();
    g.restore();
}

GridGlyph.prototype.toSVG = function() {
    var p = new SVGPath();
    for (var y = 0; y <= this._height; y += 10) {
        p.moveTo(-5000, y);
        p.lineTo(5000, y);
    }
    
    return makeElementNS(
        NS_SVG, 'path',
        null,
        {d: p.toPathData(),
         fill: 'none',
         stroke: 'black',
         strokeWidth: '0.1px'});
}

function StarGlyph(x, r, points, fill, stroke) {
    PathGlyphBase.call(this, stroke, fill);
    this._x = x;
    this._r = r;
    this._points = points;
}

StarGlyph.prototype = Object.create(PathGlyphBase.prototype);

StarGlyph.prototype.min = function() {
    return this._x - this._r;
}

StarGlyph.prototype.max = function() {
    return this._x + this._r;
}

StarGlyph.prototype.height = function() {
    return 2 * this._r;
}

StarGlyph.prototype.drawPath = function(g) {
    var midX = this._x, midY = this._r, r = this._r;
    for (var p = 0; p < this._points; ++p) {
        var theta = (p * 6.28) / this._points;
        var px = midX + r*Math.sin(theta);
        var py = midY - r*Math.cos(theta);
        if (p == 0) {
            g.moveTo(px, py);
        } else {
            g.lineTo(px, py);
        }
        theta = ((p+0.5) * 6.28) / this._points;
        px = midX + 0.4*r*Math.sin(theta);
        py = midY - 0.4*r*Math.cos(theta);
        g.lineTo(px, py);
    }
    g.closePath();
}

function PlimsollGlyph(x, height, overhang, fill, stroke) {
    this._x = x;
    this._height = height;
    this._overhang = overhang;
    this._fill = fill;
    this._stroke = stroke;
    this._hh = height / 2;
}

PlimsollGlyph.prototype.draw = function(g) {
    var hh = this._height/2;
    g.fillStyle = this._stroke;
    g.beginPath();
    g.arc(this._x, hh, hh - this._overhang, 0, 6.29);
    g.moveTo(this._x, 0);
    g.lineTo(this._x, this._height);

    if (this._fill) {
        g.fillStyle = this._fill;
        g.fill();
    }

    if (this._stroke) {
        g.strokeStyle = this._stroke;
        g.stroke();
    }
}

PlimsollGlyph.prototype.toSVG = function() {
    var hh = this._hh;
    return makeElementNS(NS_SVG, 'g', 
        [makeElementNS(NS_SVG, 'circle', null, {cx: this._x, cy: hh, r: hh - this._overhang}),
         makeElementNS(NS_SVG, 'line', null, {x1: this._x, y1: 0, x2: this._x, y2: this._height})],
        {fill: this._fill || 'none',
         stroke: this._stroke || 'none',
         strokeWidth: '1px'});
}

PlimsollGlyph.prototype.min = function() {
    return this._x - this._hh;
}

PlimsollGlyph.prototype.max = function() {
    return this._x + this._hh;
}

PlimsollGlyph.prototype.height = function() {
    return this._height;
}

/* -*- mode: javascript; c-basic-offset: 4; indent-tabs-mode: nil -*- */

// 
// Dalliance Genome Explorer
// (c) Thomas Down 2006-2013
//
// session.js
//

"use strict";

Browser.prototype.nukeStatus = function() {
    delete localStorage['dalliance.' + this.cookieKey + '.view-chr'];
    delete localStorage['dalliance.' + this.cookieKey + '.view-start'];
    delete localStorage['dalliance.' + this.cookieKey + '.view-end'];
    delete localStorage['dalliance.' + this.cookieKey + '.current-seq-length'];
    delete localStorage['dalliance.' + this.cookieKey + '.showing-alt-zoom'];
    delete localStorage['dalliance.' + this.cookieKey + '.saved-zoom'];

    delete localStorage['dalliance.' + this.cookieKey + '.sources'];
    delete localStorage['dalliance.' + this.cookieKey + '.hubs'];
    delete localStorage['dalliance.' + this.cookieKey + '.version'];

    delete localStorage['dalliance.' + this.cookieKey + '.reverse-scrolling'];
    delete localStorage['dalliance.' + this.cookieKey + '.reverse-key-scrolling'];
    delete localStorage['dalliance.' + this.cookieKey + '.ruler-location'];
}

Browser.prototype.storeStatus = function() {
    this.storeViewStatus();
    this.storeTierStatus();
}

Browser.prototype.storeViewStatus = function() {
    if (!this.cookieKey || this.noPersist || this.noPersistView) {
        return;
    }

    localStorage['dalliance.' + this.cookieKey + '.view-chr'] = this.chr;
    localStorage['dalliance.' + this.cookieKey + '.view-start'] = this.viewStart|0;
    localStorage['dalliance.' + this.cookieKey + '.view-end'] = this.viewEnd|0
    localStorage['dalliance.' + this.cookieKey + '.showing-alt-zoom'] = '' + this.isSnapZooming;
    localStorage['dalliance.' + this.cookieKey + '.saved-zoom'] = this.savedZoom;
    if (this.currentSeqMax) {
	   localStorage['dalliance.' + this.cookieKey + '.current-seq-length'] = this.currentSeqMax;
    }
}


Browser.prototype.storeTierStatus = function() {
    if (!this.cookieKey || this.noPersist) {
        return;
    }

    var currentSourceList = [];
    for (var t = 0; t < this.tiers.length; ++t) {
        var tt = this.tiers[t];
        var ts = tt.dasSource;
        if (!ts.noPersist) {
            currentSourceList.push({source: tt.dasSource, config: tt.config || {}});
        }
    }
    localStorage['dalliance.' + this.cookieKey + '.sources'] = JSON.stringify(currentSourceList);


    var coveredHubURLs = {};
    var currentHubList = [];
    for (var hi = 0; hi < this.hubObjects.length; ++hi) {
        var tdb = this.hubObjects[hi];
        var hc = {url: tdb.hub.url, genome: tdb.genome};
        if (tdb.credentials)
            hc.credentials = tdb.credentials;
        if (tdb.mapping)
            hc.mapping = tdb.mapping;
        coveredHubURLs[hc.url] = true;
        currentHubList.push(hc);
    }

    // Needed to handle hubs that failed to connect, or hubs that haven't
    // connected yet when we're called soon after startup.
    for (var hi = 0; hi < this.hubs.length; ++hi) {
        var hc = this.hubs[hi];
        if (typeof hc === 'string')
            hc = {url: hc};
        if (!coveredHubURLs[hc.url])
            currentHubList.push(hc);
    }

    localStorage['dalliance.' + this.cookieKey + '.hubs'] = JSON.stringify(currentHubList);

    localStorage['dalliance.' + this.cookieKey + '.reverse-scrolling'] = this.reverseScrolling;
    localStorage['dalliance.' + this.cookieKey + '.reverse-key-scrolling'] = this.reverseKeyScrolling;
    localStorage['dalliance.' + this.cookieKey + '.ruler-location'] = this.rulerLocation;

    localStorage['dalliance.' + this.cookieKey + '.export-ruler'] = this.exportRuler;
    localStorage['dalliance.' + this.cookieKey + '.export-highlights'] = this.exportHighlights;
    
    localStorage['dalliance.' + this.cookieKey + '.version'] = VERSION.CONFIG;
}

Browser.prototype.restoreStatus = function() {
    if (this.noPersist)
        return;
    
    var storedConfigVersion = localStorage['dalliance.' + this.cookieKey + '.version'];
    if (storedConfigVersion) {
        storedConfigVersion = storedConfigVersion|0;
    } else {
        storedConfigVersion = -100;
    }
    if (VERSION.CONFIG != storedConfigVersion) {
        return;
    }

    var storedConfigHash = localStorage['dalliance.' + this.cookieKey + '.configHash'] || '';
    var pageConfigHash = hex_sha1(miniJSONify({sources: this.sources, hubs: this.hubs}));
    if (pageConfigHash != storedConfigHash) {
        localStorage['dalliance.' + this.cookieKey + '.configHash'] = pageConfigHash;
        return;
    }

    var defaultSourcesByURI = {};
    for (var si = 0; si < this.sources.length; ++si) {
        var source = this.sources[si];
        var uri = sourceDataURI(source);
        var ul = defaultSourcesByURI[uri];
        if (!ul)
            defaultSourcesByURI[uri] = ul = [];
        ul.push(source);
        
    }

    if (!this.noPersistView) {
        var qChr = localStorage['dalliance.' + this.cookieKey + '.view-chr'];
        var qMin = localStorage['dalliance.' + this.cookieKey + '.view-start']|0;
        var qMax = localStorage['dalliance.' + this.cookieKey + '.view-end']|0;
        if (qChr && qMin && qMax) {
        	this.chr = qChr;
        	this.viewStart = qMin;
        	this.viewEnd = qMax;
        	
        	var csm = localStorage['dalliance.' + this.cookieKey + '.current-seq-length'];
        	if (csm) {
        	    this.currentSeqMax = csm|0;
        	}

            this.isSnapZooming = (localStorage['dalliance.' + this.cookieKey + '.showing-alt-zoom']) == 'true';

            var sz = parseFloat(localStorage['dalliance.' + this.cookieKey + '.saved-zoom']);
            if (typeof sz === 'number' && !isNaN(sz)) {
                this.savedZoom = sz;
            }
        }
    }

    var rs = localStorage['dalliance.' + this.cookieKey + '.reverse-scrolling'];
    this.reverseScrolling = (rs && rs == 'true');
    var rks = localStorage['dalliance.' + this.cookieKey + '.reverse-key-scrolling'];
    this.reverseKeyScrolling = (rks && rks == 'true');

    var rl = localStorage['dalliance.' + this.cookieKey + '.ruler-location'];
    if (rl)
        this.rulerLocation = rl;

    var x = localStorage['dalliance.' + this.cookieKey + '.export-ruler'];
    if (x)
        this.exportRuler = (x === 'true');
    var x = localStorage['dalliance.' + this.cookieKey + '.export-highlights'];
    if (x)
        this.exportHighlights = (x === 'true');

    var sourceStr = localStorage['dalliance.' + this.cookieKey + '.sources'];
    if (sourceStr) {
	    var storedSources = JSON.parse(sourceStr);
        this.sources = [];
        this.restoredConfigs = [];
        for (var si = 0; si < storedSources.length; ++si) {
            var source = this.sources[si] = storedSources[si].source;
            this.restoredConfigs[si] = storedSources[si].config;
            var uri = sourceDataURI(source);
            var ul = defaultSourcesByURI[uri] || [];
            for (var osi = 0; osi < ul.length; ++osi) {    
                var oldSource = ul[osi];
                if (sourcesAreEqual(source, oldSource)) {
                    if (oldSource.featureInfoPlugin) {
                        source.featureInfoPlugin = oldSource.featureInfoPlugin;
                    }
                }
            }
        }
    }

    var hubStr = localStorage['dalliance.' + this.cookieKey + '.hubs'];
    if (hubStr) {
        this.hubs = JSON.parse(hubStr);
    }

    return true;
}

Browser.prototype.reset = function() {
    for (var i = this.tiers.length - 1; i >= 0; --i) {
       this.removeTier({index: i}, true);
    }
    for (var i = 0; i < this.defaultSources.length; ++i) {
        var s = this.defaultSources[i];
        if (!s.disabled) 
            this.addTier(this.defaultSources[i]);
    }

    this.highlights.splice(0, this.highlights.length);

    this.setLocation(this.defaultChr, this.defaultStart, this.defaultEnd);
}/* -*- mode: javascript; c-basic-offset: 4; indent-tabs-mode: nil -*- */

// 
// Dalliance Genome Explorer
// (c) Thomas Down 2006-2013
//
// sourceadapters.js
//

"use strict";

var __dalliance_sourceAdapterFactories = {};

function dalliance_registerSourceAdapterFactory(type, factory) {
    __dalliance_sourceAdapterFactories[type] = factory;
};


(function(global) {
    var __dalliance_parserFactories = {};

    global.dalliance_registerParserFactory = function(type, factory) {
        __dalliance_parserFactories[type] = factory;
    };

    global.dalliance_makeParser = function(type) {
        if (__dalliance_parserFactories[type]) {
            return __dalliance_parserFactories[type](type);
        }
    };
}(this));



DasTier.prototype.initSources = function() {
    var thisTier = this;
    var fs = new DummyFeatureSource(), ss;

    if (this.dasSource.tier_type == 'sequence') {
        if (this.dasSource.twoBitURI) {
            ss = new TwoBitSequenceSource(this.dasSource);
        } else {
            ss = new DASSequenceSource(this.dasSource);
        }
    } else {
        fs = this.browser.createFeatureSource(this.dasSource);
    }

    this.featureSource = fs;
    this.sequenceSource = ss;

    if (this.featureSource && this.featureSource.addChangeListener) {
        this.featureSource.addChangeListener(function() {
            thisTier.browser.refreshTier(thisTier);
        });
    }

}

Browser.prototype.createFeatureSource = function(config) {
    var fs = this.sourceCache.get(config);
    if (fs) {
        return fs;
    }

    if (config.tier_type && __dalliance_sourceAdapterFactories[config.tier_type]) {
        var saf = __dalliance_sourceAdapterFactories[config.tier_type];
        fs = saf(config).features;
    } else if (config.bwgURI || config.bwgBlob) {
        var worker = this.getWorker();
        if (worker)
            fs = new RemoteBWGFeatureSource(config, worker);
        else
            fs = new BWGFeatureSource(config);
    } else if (config.bamURI || config.bamBlob) {
        var worker = this.getWorker();
        if (worker)
            fs = new RemoteBAMFeatureSource(config, worker);
        else
            fs = new BAMFeatureSource(config);
    } else if (config.bamblrURI) {
        fs = new BamblrFeatureSource(config);
    } else if (config.jbURI) {
        fs = new JBrowseFeatureSource(config);
    } else if (config.uri || config.features_uri) {
        fs = new DASFeatureSource(config);
    }

    if (config.overlay) {
        var sources = [];
        if (fs)
            sources.push(new CachingFeatureSource(fs));

        for (var oi = 0; oi < config.overlay.length; ++oi) {
            sources.push(this.createFeatureSource(config.overlay[oi]));
        }
        fs = new OverlayFeatureSource(sources, config);
    }

    if (config.mapping) {
        fs = new MappedFeatureSource(fs, this.chains[config.mapping]);
    }

    if (config.name && !fs.name) {
        fs.name = config.name;
    }

    if (fs != null) {
        fs = new CachingFeatureSource(fs);
        this.sourceCache.put(config, fs);
    }
    return fs;
}

function SourceCache() {
    this.sourcesByURI = {}
}

SourceCache.prototype.get = function(conf) {
    var scb = this.sourcesByURI[sourceDataURI(conf)];
    if (scb) {
        for (var si = 0; si < scb.configs.length; ++si) {
            if (sourcesAreEqual(scb.configs[si], conf)) {
                return scb.sources[si];
            }
        }
    }
}

SourceCache.prototype.put = function(conf, source) {
    var uri = sourceDataURI(conf);
    var scb = this.sourcesByURI[uri];
    if (!scb) {
        scb = {configs: [], sources: []};
        this.sourcesByURI[uri] = scb;
    }
    scb.configs.push(conf);
    scb.sources.push(source);
}

var __cfs_id_seed = 0;

function CachingFeatureSource(source) {
    var thisB = this;

    this.source = source;
    this.cfsid = 'cfs' + (++__cfs_id_seed);
    if (source.name) {
        this.name = source.name;
    }
    if (source.addChangeListener) {
        source.addChangeListener(function() {
            thisB.cfsid = 'cfs' + (++__cfs_id_seed);
        });
    }
}

CachingFeatureSource.prototype.addReadinessListener = function(listener) {
    if (this.source.addReadinessListener)
        return this.source.addReadinessListener(listener);
    else
        listener(null);
}

CachingFeatureSource.prototype.search = function(query, callback) {
    if (this.source.search)
        return this.source.search(query, callback);
}

CachingFeatureSource.prototype.getDefaultFIPs = function(callback) {
    if (this.source.getDefaultFIPs)
        return this.source.getDefaultFIPs(callback); 
}

CachingFeatureSource.prototype.getStyleSheet = function(callback) {
    this.source.getStyleSheet(callback);
}

CachingFeatureSource.prototype.getScales = function() {
    return this.source.getScales();
}

CachingFeatureSource.prototype.addActivityListener = function(l) {
    if (this.source.addActivityListener) {
        this.source.addActivityListener(l);
    }
}

CachingFeatureSource.prototype.addChangeListener = function(l) {
    if (this.source.addChangeListener)
        this.source.addChangeListener(l);
}

CachingFeatureSource.prototype.findNextFeature = function(chr, pos, dir, callback) {
    this.source.findNextFeature(chr, pos, dir, callback);
}

CachingFeatureSource.prototype.quantFindNextFeature = function(chr, pos, dir, threshold, callback) {
    this.source.quantFindNextFeature(chr, pos, dir, threshold, callback);
}

CachingFeatureSource.prototype.capabilities = function() {
    if (this.source.capabilities) {
        return this.source.capabilities();
    } else {
        return {};
    }
}

CachingFeatureSource.prototype.fetch = function(chr, min, max, scale, types, pool, callback) {
    if (pool == null) {
        throw 'pool is null...';
    }

    var awaitedFeatures = pool.awaitedFeatures[this.cfsid];
    if (!awaitedFeatures) {
        var awaitedFeatures = new Awaited();
        pool.awaitedFeatures[this.cfsid] = awaitedFeatures;
        this.source.fetch(chr, min, max, scale, types, pool, function(status, features, scale, coverage) {
            if (!awaitedFeatures.res)
                awaitedFeatures.provide({status: status, features: features, scale: scale, coverage: coverage});
        });
    } 

    awaitedFeatures.await(function(af) {
        callback(af.status, af.features, af.scale, af.coverage);
    });
}
    
function FeatureSourceBase() {
    this.busy = 0;
    this.activityListeners = [];
    this.readinessListeners = [];
    this.readiness = null;
}

FeatureSourceBase.prototype.addReadinessListener = function(listener) {
    this.readinessListeners.push(listener);
    listener(this.readiness);
}

FeatureSourceBase.prototype.notifyReadiness = function() {
    for (var li = 0; li < this.readinessListeners.length; ++li) {
        try {
            this.readinessListeners[li](this.readiness);
        } catch (e) {
            console.log(e);
        }
    }
}

FeatureSourceBase.prototype.addActivityListener = function(listener) {
    this.activityListeners.push(listener);
}

FeatureSourceBase.prototype.notifyActivity = function() {
    for (var li = 0; li < this.activityListeners.length; ++li) {
        try {
            this.activityListeners[li](this.busy);
        } catch (e) {
            console.log(e);
        }
    }
}

FeatureSourceBase.prototype.getScales = function() {
    return null;
}

FeatureSourceBase.prototype.fetch = function(chr, min, max, scale, types, pool, cnt) {
    return cnt(null, [], 1000000000);
}

FeatureSourceBase.prototype.getStyleSheet = function(callback) {
    var stylesheet = new DASStylesheet();
    var defStyle = new DASStyle();
    defStyle.glyph = 'BOX';
    defStyle.BGCOLOR = 'blue';
    defStyle.FGCOLOR = 'black';
    stylesheet.pushStyle({type: 'default'}, null, defStyle);
    return callback(stylesheet);
}



function DASFeatureSource(dasSource) {
    this.dasSource = new DASSource(dasSource);
    this.busy = 0;
    this.activityListeners = [];
}

DASFeatureSource.prototype.addActivityListener = function(listener) {
    this.activityListeners.push(listener);
}

DASFeatureSource.prototype.notifyActivity = function() {
    for (var li = 0; li < this.activityListeners.length; ++li) {
        try {
            this.activityListeners[li](this.busy);
        } catch (e) {
            console.log(e);
        }
    }
}

DASFeatureSource.prototype.getStyleSheet = function(callback) {
    this.dasSource.stylesheet(function(stylesheet) {
	callback(stylesheet);
    }, function() {
	callback(null, "Couldn't fetch DAS stylesheet");
    });
}

DASFeatureSource.prototype.fetch = function(chr, min, max, scale, types, pool, callback) {
    if (types && types.length == 0) {
        callback(null, [], scale);
        return;
    }

    if (!this.dasSource.uri && !this.dasSource.features_uri) {
        // FIXME should this be making an error callback???
        return;
    }

    if (this.dasSource.dasStaticFeatures && this.cachedStaticFeatures) {
        return callback(null, this.cachedStaticFeatures, this.cachedStaticScale);
    }

    var tryMaxBins = (this.dasSource.maxbins !== false);
    var fops = {
        type: types
    };
    if (tryMaxBins) {
        fops.maxbins = 1 + (((max - min) / scale) | 0);
    }
    
    var thisB = this;
    thisB.busy++;
    thisB.notifyActivity();

    this.dasSource.features(
        new DASSegment(chr, min, max),
        fops,
        function(features, status) {
            
            thisB.busy--;
            thisB.notifyActivity();

            var retScale = scale;
            if (!tryMaxBins) {
                retScale = 0.1;
            }
            if (!status && thisB.dasSource.dasStaticFeatures) {
                thisB.cachedStaticFeatures = features;
                thisB.cachedStaticScale = retScale;
            }
            callback(status, features, retScale);
        }
    );
}

DASFeatureSource.prototype.findNextFeature = this.sourceFindNextFeature = function(chr, pos, dir, callback) {
    if (this.dasSource.capabilities && arrayIndexOf(this.dasSource.capabilities, 'das1:adjacent-feature') >= 0) {
        var thisB = this;
        if (this.dasAdjLock) {
            return dlog('Already looking for a next feature, be patient!');
        }
        this.dasAdjLock = true;
        var fops = {
            adjacent: chr + ':' + (pos|0) + ':' + (dir > 0 ? 'F' : 'B')
        }
        var types = thisTier.getDesiredTypes(thisTier.browser.scale);
        if (types) {
            fops.types = types;
        }
        thisTier.dasSource.features(null, fops, function(res) {
            thisB.dasAdjLock = false;
            if (res.length > 0 && res[0] != null) {
                dlog('DAS adjacent seems to be working...');
                callback(res[0]);
            }
        });
    }
};

function DASSequenceSource(dasSource) {
    this.dasSource = new DASSource(dasSource);
    this.awaitedEntryPoints = new Awaited();

    var thisB = this;
    this.dasSource.entryPoints(
        function(ep) {
            thisB.awaitedEntryPoints.provide(ep);
        });
}


DASSequenceSource.prototype.fetch = function(chr, min, max, pool, callback) {
    this.dasSource.sequence(
        new DASSegment(chr, min, max),
        function(seqs) {
            if (seqs.length == 1) {
                return callback(null, seqs[0]);
            } else {
                return callback("Didn't get sequence");
            }
        }
    );
}

DASSequenceSource.prototype.getSeqInfo = function(chr, cnt) {
    this.awaitedEntryPoints.await(function(ep) {
        for (var epi = 0; epi < ep.length; ++epi) {
            if (ep[epi].name == chr) {
                return cnt({length: ep[epi].end});
            }
        }
        return cnt();
    });
}
    

function TwoBitSequenceSource(source) {
    var thisB = this;
    this.source = source;
    this.twoBit = new Awaited();
    makeTwoBit(new URLFetchable(source.twoBitURI), function(tb, error) {
        if (error) {
            dlog(error);
        } else {
            thisB.twoBit.provide(tb);
        }
    });
}

TwoBitSequenceSource.prototype.fetch = function(chr, min, max, pool, callback) {
        this.twoBit.await(function(tb) {
            tb.fetch(chr, min, max,
                     function(seq, err) {
                         if (err) {
                             return callback(err, null);
                         } else {
                             var sequence = new DASSequence(chr, min, max, 'DNA', seq);
                             return callback(null, sequence);
                         }
                     })
        });
}

TwoBitSequenceSource.prototype.getSeqInfo = function(chr, cnt) {
    this.twoBit.await(function(tb) {
        var seq = tb.getSeq(chr);
        if (seq) {
            tb.getSeq(chr).length(function(l) {
                cnt({length: l});
            });
        } else {
            cnt();
        }
    });
}

DASFeatureSource.prototype.getScales = function() {
    return [];
}

var bwg_preflights = {};

function BWGFeatureSource(bwgSource) {
    FeatureSourceBase.call(this);

    var thisB = this;
    this.readiness = 'Connecting';
    this.bwgSource = this.opts = bwgSource;    
    thisB.bwgHolder = new Awaited();

    if (this.opts.preflight) {
        var pfs = bwg_preflights[this.opts.preflight];
        if (!pfs) {
            pfs = new Awaited();
            bwg_preflights[this.opts.preflight] = pfs;

            var req = new XMLHttpRequest();
            req.onreadystatechange = function() {
                if (req.readyState == 4) {
                    if (req.status == 200) {
                        pfs.provide('success');
                    } else {
                        pfs.provide('failure');
                    }
                }
            };
            req.open('get', this.opts.preflight + '?' + hex_sha1('salt' + Date.now()), true);    // Instead, ensure we always preflight a unique URI.
            if (this.opts.credentials) {
                req.withCredentials = true;
            }
            req.send('');
        }
        pfs.await(function(status) {
            if (status === 'success') {
                thisB.init();
            }
        });
    } else {
        thisB.init();
    }
}

BWGFeatureSource.prototype = Object.create(FeatureSourceBase.prototype);

BWGFeatureSource.prototype.init = function() {
    var thisB = this;
    var make, arg;
    if (this.bwgSource.bwgURI) {
        make = makeBwgFromURL;
        arg = this.bwgSource.bwgURI;
    } else {
        make = makeBwgFromFile;
        arg = this.bwgSource.bwgBlob;
    }

    make(arg, function(bwg, err) {
        if (err) {
            thisB.error = err;
            thisB.readiness = null;
            thisB.notifyReadiness();
            thisB.bwgHolder.provide(null);
        } else {
            thisB.bwgHolder.provide(bwg);
            thisB.readiness = null;
            thisB.notifyReadiness();
            if (bwg.type == 'bigbed') {
                bwg.getExtraIndices(function(ei) {
                    thisB.extraIndices = ei;
                });
            }
        }
    }, this.opts.credentials);
}

BWGFeatureSource.prototype.capabilities = function() {
    var caps = {leap: true};
    if (this.bwgHolder.res && this.bwgHolder.res.type == 'bigwig')
        caps.quantLeap = true;
    if (this.extraIndices && this.extraIndices.length > 0) {
        caps.search = [];
        for (var eii = 0; eii < this.extraIndices.length; ++eii) {
            caps.search.push(this.extraIndices[eii].field);
        }
    }
    return caps;
}

BWGFeatureSource.prototype.fetch = function(chr, min, max, scale, types, pool, callback) {
    var thisB = this;
    this.bwgHolder.await(function(bwg) {
        if (bwg == null) {
            return callback(thisB.error || "Can't access binary file", null, null);
        }

        // dlog('bwg: ' + bwg.name + '; want scale: ' + scale);
        var data;
        // dlog(miniJSONify(types));
        var wantDensity = !types || types.length == 0 || arrayIndexOf(types, 'density') >= 0;
/*        if (wantDensity) {
            dlog('want density; scale=' + scale);
        } */
        if (thisB.opts.clientBin) {
            wantDensity = false;
        }
        if (bwg.type == 'bigwig' || wantDensity || (typeof thisB.opts.forceReduction !== 'undefined')) {
            var zoom = -1;
            for (var z = 0; z < bwg.zoomLevels.length; ++z) {
                if (bwg.zoomLevels[z].reduction <= scale) {
                    zoom = z;
                } else {
                    break;
                }
            }
            if (typeof thisB.opts.forceReduction !== 'undefined') {
                zoom = thisB.opts.forceReduction;
            }
           // dlog('selected zoom: ' + zoom);
            if (zoom < 0) {
                data = bwg.getUnzoomedView();
            } else {
                data = bwg.getZoomedView(zoom);
            }
        } else {
            data = bwg.getUnzoomedView();
        }
        
        thisB.busy++;
        thisB.notifyActivity();
        data.readWigData(chr, min, max, function(features) {
            thisB.busy--;
            thisB.notifyActivity();

            var fs = 1000000000;
            if (bwg.type === 'bigwig') {
                var is = (max - min) / features.length / 2;
                if (is < fs) {
                    fs = is;
                }
            }
            if (thisB.opts.link) {
                for (var fi = 0; fi < features.length; ++fi) {
                    var f = features[fi];
                    if (f.label) {
                        f.links = [new DASLink('Link', thisB.opts.link.replace(/\$\$/, f.label))];
                    }
                }
            }
            callback(null, features, fs);
        });
    });
}

BWGFeatureSource.prototype.quantFindNextFeature = function(chr, pos, dir, threshold, callback) {
    // var beforeQFNF = Date.now()|0;
    var thisB = this;
    thisB.busy++;
    thisB.notifyActivity();
    this.bwgHolder.res.thresholdSearch(chr, pos, dir, threshold, function(a, b) {
        thisB.busy--;
        thisB.notifyActivity();
        // var afterQFNF = Date.now()|0;
        // console.log('QFNF took ' + (afterQFNF - beforeQFNF) + 'ms');
        return callback(a, b);
    });
}

BWGFeatureSource.prototype.findNextFeature = function(chr, pos, dir, callback) {
    var thisB = this;
    thisB.busy++;
    thisB.notifyActivity();
    this.bwgHolder.res.getUnzoomedView().getFirstAdjacent(chr, pos, dir, function(res) {
        thisB.busy--;
        thisB.notifyActivity();
        if (res.length > 0 && res[0] != null) {
            callback(res[0]);
        }
    });
}

BWGFeatureSource.prototype.getScales = function() {
    var bwg = this.bwgHolder.res;
    if (bwg /* && bwg.type == 'bigwig' */) {
        var scales = [1];  // Can we be smarter about inferring baseline scale?
        for (var z = 0; z < bwg.zoomLevels.length; ++z) {
            scales.push(bwg.zoomLevels[z].reduction);
        }
        return scales;
    } else {
        return null;
    }
}

BWGFeatureSource.prototype.search = function(query, callback) {
    if (!this.extraIndices || this.extraIndices.length == 0) {
        return callback(null, 'No indices available');
    }

    var index = this.extraIndices[0];
    return index.lookup(query, callback);
}

BWGFeatureSource.prototype.getDefaultFIPs = function(callback) {
    if (this.opts.noExtraFeatureInfo)
        return true;

    this.bwgHolder.await(function(bwg) {
        if (!bwg) return;

        if (bwg.schema && bwg.definedFieldCount < bwg.schema.fields.length) {
            var fip = function(feature, featureInfo) {
                for (var fi = bwg.definedFieldCount; fi < bwg.schema.fields.length; ++fi) {
                    var f = bwg.schema.fields[fi];
                    featureInfo.add(f.comment, feature[f.name]);
                }
            };

            callback(fip);
        } else {
            // No need to do anything.
        }
    });
}

BWGFeatureSource.prototype.getStyleSheet = function(callback) {
    var thisB = this;

    this.bwgHolder.await(function(bwg) {
        if (!bwg) {
            return callback(null, 'bbi error');
        }

    	var stylesheet = new DASStylesheet();
        if (bwg.type == 'bigbed') {
            var wigStyle = new DASStyle();
            wigStyle.glyph = 'BOX';
            wigStyle.FGCOLOR = 'black';
            wigStyle.BGCOLOR = 'blue'
            wigStyle.HEIGHT = 8;
            wigStyle.BUMP = true;
            wigStyle.LABEL = true;
            wigStyle.ZINDEX = 20;
            stylesheet.pushStyle({type: 'bigwig'}, null, wigStyle);
	    
            wigStyle.glyph = 'BOX';
            wigStyle.FGCOLOR = 'black';
            wigStyle.BGCOLOR = 'red'
            wigStyle.HEIGHT = 10;
            wigStyle.BUMP = true;
            wigStyle.ZINDEX = 20;
            stylesheet.pushStyle({type: 'bb-translation'}, null, wigStyle);
                    
            var tsStyle = new DASStyle();
            tsStyle.glyph = 'BOX';
            tsStyle.FGCOLOR = 'black';
            tsStyle.BGCOLOR = 'white';
            tsStyle.HEIGHT = 10;
            tsStyle.ZINDEX = 10;
            tsStyle.BUMP = true;
            tsStyle.LABEL = true;
            stylesheet.pushStyle({type: 'bb-transcript'}, null, tsStyle);

            var densStyle = new DASStyle();
            densStyle.glyph = 'HISTOGRAM';
            densStyle.COLOR1 = 'white';
            densStyle.COLOR2 = 'black';
            densStyle.HEIGHT=30;
            stylesheet.pushStyle({type: 'density'}, null, densStyle);
        } else {
            var wigStyle = new DASStyle();
            wigStyle.glyph = 'HISTOGRAM';
            wigStyle.COLOR1 = 'white';
            wigStyle.COLOR2 = 'black';
            wigStyle.HEIGHT=30;
            stylesheet.pushStyle({type: 'default'}, null, wigStyle);
        }

        if (bwg.definedFieldCount == 12 && bwg.fieldCount >= 14) {
            stylesheet.geneHint = true;
        }

    	return callback(stylesheet);
    });
}

function RemoteBWGFeatureSource(bwgSource, worker) {
    FeatureSourceBase.call(this);

    this.worker = worker;
    this.readiness = 'Connecting';
    this.bwgSource = this.opts = bwgSource;
    this.keyHolder = new Awaited();
    this.init();
}

RemoteBWGFeatureSource.prototype = Object.create(FeatureSourceBase.prototype);

RemoteBWGFeatureSource.prototype.init = function() {
    var thisB = this;
    var uri = this.bwgSource.uri || this.bwgSource.bwgURI;
    var blob = this.bwgSource.blob || this.bwgSource.bwgBlob;

    var cnt = function(key, err) {
        thisB.readiness = null;
        thisB.notifyReadiness();

        if (key) {
            thisB.worker.postCommand({command: 'meta', connection: key}, function(meta, err) {
                if (err) {
                    thisB.error = err;
                    thisB.keyHolder.provide(null);
                } else {
                    thisB.meta = meta;
                    thisB.keyHolder.provide(key);
                }
            });
        } else {
            thisB.error = err;
            thisB.keyHolder.provide(null);
        }
    };

    if (blob) {
        this.worker.postCommand({command: 'connectBBI', blob: blob}, cnt);
    } else {
        this.worker.postCommand({command: 'connectBBI', uri: uri}, cnt); 
    }
}

RemoteBWGFeatureSource.prototype.capabilities = function() {
    var caps = {leap: true};

    if (this.meta && this.meta.type == 'bigwig')
        caps.quantLeap = true;
    if (this.meta && this.meta.extraIndices && this.meta.extraIndices.length > 0) {
        caps.search = [];
        for (var eii = 0; eii < this.meta.extraIndices.length; ++eii) {
            caps.search.push(this.meta.extraIndices[eii].field);
        }
    }
    return caps;
}

RemoteBWGFeatureSource.prototype.fetch = function(chr, min, max, scale, types, pool, callback) {
    var thisB = this;

    thisB.busy++;
    thisB.notifyActivity();

    this.keyHolder.await(function(key) {
        if (!key) {
            thisB.busy--;
            thisB.notifyActivity();
            return callback(thisB.error || "Can't access binary file", null, null);
        }

        var zoom = -1;
        var wantDensity = !types || types.length == 0 || arrayIndexOf(types, 'density') >= 0;
        if (thisB.opts.clientBin) {
            wantDensity = false;
        }
        if (thisB.meta.type == 'bigwig' || wantDensity || (typeof thisB.opts.forceReduction !== 'undefined')) {
            for (var z = 1; z < thisB.meta.zoomLevels.length; ++z) {
                if (thisB.meta.zoomLevels[z] <= scale) {
                    zoom = z - 1; // Scales returned in metadata start at 1, unlike "real" zoom levels.
                } else {
                    break;
                }
            }
            if (typeof thisB.opts.forceReduction !== 'undefined') {
                zoom = thisB.opts.forceReduction;
            }
        }
        
        thisB.worker.postCommand({command: 'fetch', connection: key, chr: chr, min: min, max: max, zoom: zoom}, function(features, error) {
            thisB.busy--;
            thisB.notifyActivity();

            var fs = 1000000000;
            if (thisB.meta.type === 'bigwig') {
                var is = (max - min) / features.length / 2;
                if (is < fs) {
                    fs = is;
                }
            } 
            if (thisB.opts.link) {
                for (var fi = 0; fi < features.length; ++fi) {
                    var f = features[fi];
                    if (f.label) {
                        f.links = [new DASLink('Link', thisB.opts.link.replace(/\$\$/, f.label))];
                    }
                }
            } 
            callback(error, features, fs);
        });
    });
}


RemoteBWGFeatureSource.prototype.quantFindNextFeature = function(chr, pos, dir, threshold, callback) {
    var thisB = this;
    this.busy++;
    this.notifyActivity();
    this.worker.postCommand({command: 'quantLeap', connection: this.keyHolder.res, chr: chr, pos: pos, dir: dir, threshold: threshold, under: false}, function(result, err) {
        console.log(result, err);
        thisB.busy--;
        thisB.notifyActivity();
        return callback(result, err);
    });
}

RemoteBWGFeatureSource.prototype.findNextFeature = function(chr, pos, dir, callback) {
    var thisB = this;
    this.busy++;
    this.notifyActivity();
    this.worker.postCommand({command: 'leap', connection: this.keyHolder.res, chr: chr, pos: pos, dir: dir}, function(result, err) {
        thisB.busy--;
        thisB.notifyActivity();
        if (result.length > 0 && result[0] != null) {
            callback(result[0]);
        }
    });
}

RemoteBWGFeatureSource.prototype.getScales = function() {
    var meta = this.meta;
    if (meta) {
        return meta.zoomLevels;
    } else {
        return null;
    }
}

RemoteBWGFeatureSource.prototype.search = function(query, callback) {
    if (!this.meta.extraIndices || this.meta.extraIndices.length == 0) {
        return callback(null, 'No indices available');
    }

    var thisB = this;
    this.busy++;
    this.notifyActivity();
    var index = this.meta.extraIndices[0];
    this.worker.postCommand({command: 'search', connection: this.keyHolder.res, query: query, index: index}, function(result, err) {
        thisB.busy--;
        thisB.notifyActivity();

        callback(result, err);
    });
}

RemoteBWGFeatureSource.prototype.getDefaultFIPs = function(callback) {
    if (this.opts.noExtraFeatureInfo)
        return true;

    var thisB = this;
    this.keyHolder.await(function(key) {
        var bwg = thisB.meta;
        if (!bwg) return;

        if (bwg.schema && bwg.definedFieldCount < bwg.schema.fields.length) {
            var fip = function(feature, featureInfo) {
                for (var fi = bwg.definedFieldCount; fi < bwg.schema.fields.length; ++fi) {
                    var f = bwg.schema.fields[fi];
                    featureInfo.add(f.comment, feature[f.name]);
                }
            };

            callback(fip);
        } else {
            // No need to do anything.
        }
    });
} 

RemoteBWGFeatureSource.prototype.getStyleSheet = function(callback) {
    var thisB = this;

    this.keyHolder.await(function(key) {
        var bwg = thisB.meta;
        if (!bwg) {
            return callback(null, 'bbi error');
        } 

        var stylesheet = new DASStylesheet();
        if (bwg.type == 'bigbed') {
            var wigStyle = new DASStyle();
            wigStyle.glyph = 'BOX';
            wigStyle.FGCOLOR = 'black';
            wigStyle.BGCOLOR = 'blue'
            wigStyle.HEIGHT = 8;
            wigStyle.BUMP = true;
            wigStyle.LABEL = true;
            wigStyle.ZINDEX = 20;
            stylesheet.pushStyle({type: 'bigwig'}, null, wigStyle);
        
            wigStyle.glyph = 'BOX';
            wigStyle.FGCOLOR = 'black';
            wigStyle.BGCOLOR = 'red'
            wigStyle.HEIGHT = 10;
            wigStyle.BUMP = true;
            wigStyle.ZINDEX = 20;
            stylesheet.pushStyle({type: 'bb-translation'}, null, wigStyle);
                    
            var tsStyle = new DASStyle();
            tsStyle.glyph = 'BOX';
            tsStyle.FGCOLOR = 'black';
            tsStyle.BGCOLOR = 'white';
            tsStyle.HEIGHT = 10;
            tsStyle.ZINDEX = 10;
            tsStyle.BUMP = true;
            tsStyle.LABEL = true;
            stylesheet.pushStyle({type: 'bb-transcript'}, null, tsStyle);

            var densStyle = new DASStyle();
            densStyle.glyph = 'HISTOGRAM';
            densStyle.COLOR1 = 'white';
            densStyle.COLOR2 = 'black';
            densStyle.HEIGHT=30;
            stylesheet.pushStyle({type: 'density'}, null, densStyle);
        } else {
            var wigStyle = new DASStyle();
            wigStyle.glyph = 'HISTOGRAM';
            wigStyle.COLOR1 = 'white';
            wigStyle.COLOR2 = 'black';
            wigStyle.HEIGHT=30;
            stylesheet.pushStyle({type: 'default'}, null, wigStyle);
        }


        if (bwg.definedFieldCount == 12 && bwg.fieldCount >= 14) {
            stylesheet.geneHint = true;
        } 

        return callback(stylesheet);
    });
}


function BamblrFeatureSource(bamblrSource) {
    this.bamblr = bamblrSource.bamblrURI;
}

BamblrFeatureSource.prototype.getScales = function() {
    return [];
}

BamblrFeatureSource.prototype.getStyleSheet = function(callback) {
    var stylesheet = new DASStylesheet();

    var densStyle = new DASStyle();
    densStyle.glyph = 'HISTOGRAM';
    densStyle.COLOR1 = 'black';
    densStyle.COLOR2 = 'red';
    densStyle.HEIGHT=30;
    stylesheet.pushStyle({type: 'default'}, null, densStyle);

    return callback(stylesheet);
}

BamblrFeatureSource.prototype.fetch = function(chr, min, max, scale, types, pool, callback) {
    var rez = scale|0;
    if (rez < 1) {
        rez = 1;
    }
    var url = this.bamblr + '?seq=' + chr + '&min=' + min + '&max=' + max + '&rez=' + rez;
    new URLFetchable(url).fetch(function(data) {
        if (data == null) {
            dlog('failing bamblr');
            return;
        } else {
            var id = new Int32Array(data);
            var features = [];
            for (var ri = 0; ri < id.length; ++ri) {
                var f = new DASFeature();
                f.min = min + (ri * rez)
                f.max = f.min + rez - 1;
                f.segment = chr;
                f.type = 'bamblr';
                f.score = id[ri];
                features.push(f);
            }
            callback(null, features, rez);
            return;
        }
    });
}

function bamRecordToFeature(r) {
    if (r.flag & 0x4)
        return; 
    
    var len;
    if (r.seq)
        len = r.seq.length;
    else 
        len = r.seqLength;
    
    if (r.cigar) {
        len = 0;
        var ops = parseCigar(r.cigar);
        for (var ci = 0; ci < ops.length; ++ci) {
            var co = ops[ci];
            if (co.op == 'M' || co.op == 'D')
                len += co.cnt;
        }
    }

    var f = new DASFeature();
    f.min = r.pos + 1;
    f.max = r.pos + len;
    f.segment = r.segment;
    f.type = 'bam';
    f.id = r.readName;
    f.notes = [/* 'Sequence=' + r.seq, 'CIGAR=' + r.cigar, */ 'MQ=' + r.mq];
    f.cigar = r.cigar;
    f.seq = r.seq;
    f.quals = r.quals;

    return f;
}

function BAMFeatureSource(bamSource) {
    FeatureSourceBase.call(this);

    var thisB = this;
    this.bamSource = bamSource;
    this.opts = {credentials: bamSource.credentials, preflight: bamSource.preflight};
    this.bamHolder = new Awaited();
    
    if (this.opts.preflight) {
        var pfs = bwg_preflights[this.opts.preflight];
        if (!pfs) {
            pfs = new Awaited();
            bwg_preflights[this.opts.preflight] = pfs;

            var req = new XMLHttpRequest();
            req.onreadystatechange = function() {
                if (req.readyState == 4) {
                    if (req.status == 200) {
                        pfs.provide('success');
                    } else {
                        pfs.provide('failure');
                    }
                }
            };
            // req.setRequestHeader('cache-control', 'no-cache');    /* Doesn't work, not an allowed request header in CORS */
            req.open('get', this.opts.preflight + '?' + hex_sha1('salt' + Date.now()), true);    // Instead, ensure we always preflight a unique URI.
            if (this.opts.credentials) {
                req.withCredentials = 'true';
            }
            req.send('');
        }
        pfs.await(function(status) {
            if (status === 'success') {
                thisB.init();
            }
        });
    } else {
        thisB.init();
    }
}

BAMFeatureSource.prototype = Object.create(FeatureSourceBase.prototype);

BAMFeatureSource.prototype.init = function() {
    var thisB = this;
    var bamF, baiF;
    if (this.bamSource.bamBlob) {
        bamF = new BlobFetchable(this.bamSource.bamBlob);
        baiF = new BlobFetchable(this.bamSource.baiBlob);
    } else {
        bamF = new URLFetchable(this.bamSource.bamURI, {credentials: this.opts.credentials});
        baiF = new URLFetchable(this.bamSource.baiURI || (this.bamSource.bamURI + '.bai'), {credentials: this.opts.credentials});
    }
    makeBam(bamF, baiF, function(bam, err) {
        thisB.readiness = null;
        thisB.notifyReadiness();

        if (bam) {
            thisB.bamHolder.provide(bam);
        } else {
            thisB.error = err;
            thisB.bamHolder.provide(null);
        }
    });
}

BAMFeatureSource.prototype.fetch = function(chr, min, max, scale, types, pool, callback) {
    var light = types && (types.length == 1) && (types[0] == 'density');

    var thisB = this;
    
    thisB.busy++;
    thisB.notifyActivity();
    
    this.bamHolder.await(function(bam) {
        if (!bam) {
            thisB.busy--;
            thisB.notifyActivity();
            return callback(thisB.error || "Couldn't fetch BAM");
        }

        bam.fetch(chr, min, max, function(bamRecords, error) {
            thisB.busy--;
            thisB.notifyActivity();

            if (error) {
                callback(error, null, null);
            } else {
                var features = [];
                for (var ri = 0; ri < bamRecords.length; ++ri) {
                    var r = bamRecords[ri];

                    var f = bamRecordToFeature(r);
                    if (f)
                        features.push(f);
                }
                callback(null, features, 1000000000);
            }
        }, {light: light});
    });
}

BAMFeatureSource.prototype.getScales = function() {
    return 1000000000;
}

BAMFeatureSource.prototype.getStyleSheet = function(callback) {
    this.bamHolder.await(function(bam) {
	    var stylesheet = new DASStylesheet();
                
        var densStyle = new DASStyle();
        densStyle.glyph = 'HISTOGRAM';
        densStyle.COLOR1 = 'black';
        densStyle.COLOR2 = 'red';
        densStyle.HEIGHT=30;
        stylesheet.pushStyle({type: 'density'}, 'low', densStyle);
        stylesheet.pushStyle({type: 'density'}, 'medium', densStyle);

        var wigStyle = new DASStyle();
        wigStyle.glyph = '__SEQUENCE';
        wigStyle.FGCOLOR = 'black';
        wigStyle.BGCOLOR = 'blue'
        wigStyle.HEIGHT = 8;
        wigStyle.BUMP = true;
        wigStyle.LABEL = false;
        wigStyle.ZINDEX = 20;
        stylesheet.pushStyle({type: 'bam'}, 'high', wigStyle);

	    return callback(stylesheet);
    });
}


function RemoteBAMFeatureSource(bamSource, worker) {
    FeatureSourceBase.call(this);

    var thisB = this;
    this.bamSource = bamSource;
    this.worker = worker;
    this.opts = {credentials: bamSource.credentials, preflight: bamSource.preflight};
    this.keyHolder = new Awaited();
    
    this.init();
}

RemoteBAMFeatureSource.prototype = Object.create(FeatureSourceBase.prototype);

RemoteBAMFeatureSource.prototype.init = function() {    var thisB = this;
    var uri = this.bamSource.uri || this.bamSource.bamURI;
    var indexUri = this.bamSource.indexUri || this.bamSource.baiURI || uri + '.bai';

    var blob = this.bamSource.bamBlob || this.bamSource.blob;
    var indexBlob = this.bamSource.baiBlob || this.bamSource.indexBlob;

    var cnt = function(result, err) {
        thisB.readiness = null;
        thisB.notifyReadiness();

        if (result) {
            thisB.keyHolder.provide(result);
        } else {
            thisB.error = err;
            thisB.keyHolder.provide(null);
        }
    };

    if (blob) {
        this.worker.postCommand({command: 'connectBAM', blob: blob, indexBlob: indexBlob}, cnt /* , [blob, indexBlob] */);
    } else {
        this.worker.postCommand({command: 'connectBAM', uri: uri, indexUri: indexUri}, cnt); 
    }
}

RemoteBAMFeatureSource.prototype.fetch = function(chr, min, max, scale, types, pool, callback) {
    var light = types && (types.length == 1) && (types[0] == 'density');
    var thisB = this;
    
    thisB.busy++;
    thisB.notifyActivity();
    
    this.keyHolder.await(function(key) {
        if (!key) {
            thisB.busy--;
            thisB.notifyActivity();
            return callback(thisB.error || "Couldn't fetch BAM");
        }

        thisB.worker.postCommand({command: 'fetch', connection: key, chr: chr, min: min, max: max, opts: {light: light}}, function(bamRecords, error) {
            // console.log('retrieved ' + bamRecords.length + ' via worker.');

            thisB.busy--;
            thisB.notifyActivity();

            if (error) {
                callback(error, null, null);
            } else {
                var features = [];
                for (var ri = 0; ri < bamRecords.length; ++ri) {
                    var r = bamRecords[ri];
                    var f = bamRecordToFeature(r);
                    if (f)
                        features.push(f);
                }
                callback(null, features, 1000000000);
            }
        });
    });
}

RemoteBAMFeatureSource.prototype.getScales = function() {
    return 1000000000;
}

RemoteBAMFeatureSource.prototype.getStyleSheet = function(callback) {
    this.keyHolder.await(function(bam) {
        var stylesheet = new DASStylesheet();
                
        var densStyle = new DASStyle();
        densStyle.glyph = 'HISTOGRAM';
        densStyle.COLOR1 = 'black';
        densStyle.COLOR2 = 'red';
        densStyle.HEIGHT=30;
        stylesheet.pushStyle({type: 'density'}, 'low', densStyle);
        stylesheet.pushStyle({type: 'density'}, 'medium', densStyle);

        var wigStyle = new DASStyle();
        wigStyle.glyph = '__SEQUENCE';
        wigStyle.FGCOLOR = 'black';
        wigStyle.BGCOLOR = 'blue'
        wigStyle.HEIGHT = 8;
        wigStyle.BUMP = true;
        wigStyle.LABEL = false;
        wigStyle.ZINDEX = 20;
        stylesheet.pushStyle({type: 'bam'}, 'high', wigStyle);
        return callback(stylesheet);
    });
}


function MappedFeatureSource(source, mapping) {
    this.source = source;
    this.mapping = mapping;
    
    this.activityListeners = [];
    this.busy = 0;
}

MappedFeatureSource.prototype.addActivityListener = function(listener) {
    this.activityListeners.push(listener);
}

MappedFeatureSource.prototype.notifyActivity = function() {
    for (var li = 0; li < this.activityListeners.length; ++li) {
        try {
            this.activityListeners[li](this.busy);
        } catch (e) {
            console.log(e);
        }
    }
}

MappedFeatureSource.prototype.getStyleSheet = function(callback) {
    return this.source.getStyleSheet(callback);
}

MappedFeatureSource.prototype.getScales = function() {
    return this.source.getScales();
}

MappedFeatureSource.prototype.simplifySegments = function(segs, minGap) {
    if (segs.length == 0) return segs;

    segs.sort(function(s1, s2) {
        var d = s1.name - s2.name;
        if (d)
            return d;
        d = s1.start - s2.start;
        if (d)
            return d;
        return s1.end - s2.end;   // Should never come to this...?
    });

    var ssegs = [];
    var currentSeg = segs[0];
    for (var si = 0; si < segs.length; ++si) {
        var ns = segs[si];

        // console.log(ns.name + ' ' + ns.start + ' ' + ns.end);
        if (ns.name != currentSeg.name || ns.start > (currentSeg.end + minGap)) {
            ssegs.push(currentSeg);
            currentSeg = ns;
        } else {
            currentSeg = new DASSegment(currentSeg.name, Math.min(currentSeg.start, ns.start), Math.max(currentSeg.end, ns.end));
        }
    }
    ssegs.push(currentSeg);
    return ssegs;
}

MappedFeatureSource.prototype.fetch = function(chr, min, max, scale, types, pool, callback) {
    var thisB = this;
    var fetchLength = max - min + 1;

    thisB.busy++;
    thisB.notifyActivity();

    this.mapping.sourceBlocksForRange(chr, min, max, function(mseg) {
        if (mseg.length == 0) {
            thisB.busy--;
            thisB.notifyActivity();

            callback("No mapping available for this regions", [], scale);
        } else {
            mseg = thisB.simplifySegments(mseg, Math.max(100, 0.05 * fetchLength));

            var mappedFeatures = [];
            var mappedLoc = null;
            var count = mseg.length;
            var finalStatus;

            mseg.map(function(seg) {
                thisB.source.fetch(seg.name, seg.start, seg.end, scale, types, pool, function(status, features, fscale) {
                    if (status && !finalStatus)
                        finalStatus = status;

                    if (features) {
                        for (var fi = 0; fi < features.length; ++fi) {
                            var f = features[fi];
                            var sn = f.segment;
                            if (sn.indexOf('chr') == 0) {
                                sn = sn.substr(3);
                            }
                            var mmin = thisB.mapping.mapPoint(sn, f.min);
                            var mmax = thisB.mapping.mapPoint(sn, f.max);
                            if (!mmin || !mmax || mmin.seq != mmax.seq || mmin.seq != chr) {
                                // Discard feature.
                                // dlog('discarding ' + miniJSONify(f));
                                if (f.parts && f.parts.length > 0) {    // FIXME: Ugly hack to make ASTD source map properly.
                                     mappedFeatures.push(f);
                                }
                            } else {
                                f.segment = mmin.seq;
                                f.min = mmin.pos;
                                f.max = mmax.pos;
                                if (f.min > f.max) {
                                    var tmp = f.max;
                                    f.max = f.min;
                                    f.min = tmp;
                                }
                                if (mmin.flipped) {
                                    if (f.orientation == '-') {
                                        f.orientation = '+';
                                    } else if (f.orientation == '+') {
                                        f.orientation = '-';
                                    }
                                }
                                mappedFeatures.push(f);
                            }
                        }
                    }

                    var segDestCoverage = new Range(thisB.mapping.mapPoint(seg.name, seg.start).pos, thisB.mapping.mapPoint(seg.name, seg.end).pos);
                    if (mappedLoc)
                        mappedLoc = union(mappedLoc, segDestCoverage);
                    else
                        mappedLoc = segDestCoverage;
                    --count;

                    if (count == 0) {
                        thisB.busy--;
                        thisB.notifyActivity();
                        callback(finalStatus, mappedFeatures, fscale, mappedLoc);
                    }
                });
            });
        }
    });
}

function DummyFeatureSource() {
}

DummyFeatureSource.prototype.getScales = function() {
    return null;
}

DummyFeatureSource.prototype.fetch = function(chr, min, max, scale, types, pool, cnt) {
    return cnt(null, [], 1000000000);
}

DummyFeatureSource.prototype.getStyleSheet = function(callback) {
    var stylesheet = new DASStylesheet();
    var defStyle = new DASStyle();
    defStyle.glyph = 'BOX';
    defStyle.BGCOLOR = 'blue';
    defStyle.FGCOLOR = 'black';
    stylesheet.pushStyle({type: 'default'}, null, defStyle);
    return callback(stylesheet);
}

function DummySequenceSource() {
}

DummySequenceSource.prototype.fetch = function(chr, min, max, pool, cnt) {
    return cnt(null, null);
}

function JBrowseFeatureSource(source) {
    this.store = new JBrowseStore(source.jbURI, source.jbQuery);
}

JBrowseFeatureSource.prototype.getScales = function() {
    return null;
}

JBrowseFeatureSource.prototype.getStyleSheet = function(callback) {
    var stylesheet = new DASStylesheet();
    var wigStyle = new DASStyle();
    wigStyle.glyph = 'BOX';
    wigStyle.FGCOLOR = 'black';
    wigStyle.BGCOLOR = 'green'
    wigStyle.HEIGHT = 8;
    wigStyle.BUMP = true;
    wigStyle.LABEL = true;
    wigStyle.ZINDEX = 20;
    stylesheet.pushStyle({type: 'default'}, null, wigStyle);

    return callback(stylesheet);
}

JBrowseFeatureSource.prototype.fetch = function(chr, min, max, scale, types, pool, callback) {
    if (types && types.length == 0) {
        callback(null, [], scale);
        return;
    }
    
    var fops = {};

    this.store.features(
        new DASSegment(chr, min, max),
        fops,
        function(features, status) {
            callback(status, features, 100000);
        }
    );
}

function sourceAdapterIsCapable(s, cap) {
    if (!s.capabilities)
        return false;
    else return s.capabilities()[cap];
}
/* -*- mode: javascript; c-basic-offset: 4; indent-tabs-mode: nil -*- */

// 
// Dalliance Genome Explorer
// (c) Thomas Down 2006-2013
//
// jbjson.js -- query JBrowse-style REST data stores
//

function JBrowseStore(base, query) {
    this.base = base;
    this.query = query;
}

var topLevelResp;

JBrowseStore.prototype.features = function(segment, opts, callback) {
    opts = opts || {};

    url = this.base + '/features/' + segment.name;

    var filters = [];
    if (this.query) {
	filters.push(this.query);
    }
    if (segment.isBounded) {
	filters.push('start=' + segment.start);
	filters.push('end=' + segment.end);
    }
    if (filters.length > 0) {
	url = url + '?' + filters.join('&');
    }

    var req = new XMLHttpRequest();
    req.onreadystatechange = function() {
	if (req.readyState == 4) {
	    if (req.status >= 300) {
		callback(null, 'Error code ' + req.status);
	    } else {
		var jf = JSON.parse(req.response)['features'];
		var features = [];
		for (fi = 0; fi < jf.length; ++fi) {
		    var j = jf[fi];
		    
		    var f = new DASFeature();
		    f.segment = segment;
		    f.min = (j['start'] | 0) + 1;
		    f.max = j['end'] | 0;
		    if (j.name) {
			f.label = j.name;
		    }
		    f.type = j.type || 'unknown';
		    
		    features.push(f);
		}
		callback(features);
	    }
	}
	
    };
    
    req.open('GET', url, true);
    req.responseType = 'text';
    req.send('');
}
/* -*- mode: javascript; c-basic-offset: 4; indent-tabs-mode: nil -*- */

// 
// Dalliance Genome Explorer
// (c) Thomas Down 2006-2013
//
// ensembljson.js -- query the Ensembl REST API.
//

"use strict";

function EnsemblFeatureSource(source) {
    this.source = source;
    this.base = source.uri || 'http://beta.rest.ensembl.org';
    this.species = source.species || 'human';

    this.activityListeners = [];
    this.busy = 0;

    if (typeof source.type === 'string') {
        this.type = [source.type];
    } else {
        this.type = source.type || ['regulatory'];
    }
}

EnsemblFeatureSource.prototype.addActivityListener = function(listener) {
    this.activityListeners.push(listener);
}

EnsemblFeatureSource.prototype.notifyActivity = function() {
    for (var li = 0; li < this.activityListeners.length; ++li) {
        try {
            this.activityListeners[li](this.busy);
        } catch (e) {
            console.log(e);
        }
    }
}


EnsemblFeatureSource.prototype.getStyleSheet = function(callback) {
    var stylesheet = new DASStylesheet();

    var tsStyle = new DASStyle();
    tsStyle.glyph = '__NONE';
    if (this.type.indexOf('exon') >= 0)
        stylesheet.pushStyle({type: 'transcript'}, null, tsStyle);
    if (this.type.indexOf('exon') >= 0 || this.type.indexOf('transcript') >= 0)
        stylesheet.pushStyle({type: 'gene'}, null, tsStyle);

    var cdsStyle = new DASStyle();
    cdsStyle.glyph = 'BOX';
    cdsStyle.FGCOLOR = 'black';
    cdsStyle.BGCOLOR = 'red'
    cdsStyle.HEIGHT = 8;
    cdsStyle.BUMP = true;
    cdsStyle.LABEL = true;
    cdsStyle.ZINDEX = 10;
    stylesheet.pushStyle({type: 'cds'}, null, cdsStyle);

    {
        var varStyle = new DASStyle();
        varStyle.glyph = 'SQUARE';
        varStyle.BUMP = 'yes';
        varStyle.LABEL = 'no';
        // varStyle.BGCOLOR = '#888888';
        varStyle.FGCOLOR = 'blue';
        stylesheet.pushStyle({type: 'variation', method: '.+_UTR_variant'}, null, varStyle);
    }
    {
        var varStyle = new DASStyle();
        varStyle.glyph = 'TRIANGLE';
        varStyle.DIRECTION = 'S';
        varStyle.BUMP = 'yes';
        varStyle.LABEL = 'no';
        // varStyle.BGCOLOR = '#888888';
        varStyle.FGCOLOR = 'blue';
        stylesheet.pushStyle({type: 'variation', method: 'missense_variant'}, null, varStyle);
    }
    {
        var varStyle = new DASStyle();
        varStyle.glyph = 'TRIANGLE';
        varStyle.DIRECTION = 'N';
        varStyle.BUMP = 'yes';
        varStyle.LABEL = 'no';
        // varStyle.BGCOLOR = '#888888';
        varStyle.FGCOLOR = 'blue';
        stylesheet.pushStyle({type: 'variation', method: 'splice_.+_variant'}, null, varStyle);
    }
    {
        var varStyle = new DASStyle();
        varStyle.glyph = 'STAR';
        varStyle.POINTS = 6;
        varStyle.BUMP = 'yes';
        varStyle.LABEL = 'no';
        // varStyle.BGCOLOR = '#888888';
        varStyle.FGCOLOR = 'blue';
        stylesheet.pushStyle({type: 'variation', method: 'regulatory_region_variant'}, null, varStyle);
    }
    {
        var varStyle = new DASStyle();
        varStyle.glyph = 'PLIMSOLL';
        varStyle.BUMP = 'yes';
        varStyle.LABEL = 'no';
        // varStyle.BGCOLOR = '#888888';
        varStyle.FGCOLOR = 'rgb(50,80,255)';
        varStyle.STROKECOLOR = 'black';
        stylesheet.pushStyle({type: 'variation'}, null, varStyle);
    }
        {
        var varStyle = new DASStyle();
        varStyle.glyph = 'SQUARE';
        varStyle.BUMP = 'yes';
        varStyle.LABEL = 'no';
        varStyle.BGCOLOR = '#888888';
        varStyle.FGCOLOR = 'red';
        stylesheet.pushStyle({type: 'indel', method: '.+_UTR_variant'}, null, varStyle);
    }
    {
        var varStyle = new DASStyle();
        varStyle.glyph = 'TRIANGLE';
        varStyle.DIRECTION = 'S';
        varStyle.BUMP = 'yes';
        varStyle.LABEL = 'no';
        varStyle.BGCOLOR = '#888888';
        varStyle.FGCOLOR = 'red';
        stylesheet.pushStyle({type: 'indel', method: 'missense_variant'}, null, varStyle);
    }
    {
        var varStyle = new DASStyle();
        varStyle.glyph = 'TRIANGLE';
        varStyle.DIRECTION = 'N';
        varStyle.BUMP = 'yes';
        varStyle.LABEL = 'no';
        varStyle.BGCOLOR = '#888888';
        varStyle.FGCOLOR = 'red';
        stylesheet.pushStyle({type: 'indel', method: 'splice_.+_variant'}, null, varStyle);
    }
    {
        var varStyle = new DASStyle();
        varStyle.glyph = 'STAR';
        varStyle.POINTS = 6;
        varStyle.BUMP = 'yes';
        varStyle.LABEL = 'no';
        varStyle.BGCOLOR = '#888888';
        varStyle.FGCOLOR = 'red';
        stylesheet.pushStyle({type: 'indel', method: 'regulatory_region_variant'}, null, varStyle);
    }
    {
        var varStyle = new DASStyle();
        varStyle.glyph = 'PLIMSOLL';
        varStyle.BUMP = 'yes';
        varStyle.LABEL = 'no';
        varStyle.BGCOLOR = '#888888';
        varStyle.FGCOLOR = 'red';
        varStyle.STROKECOLOR = 'black';
        stylesheet.pushStyle({type: 'indel'}, null, varStyle);
    }

    var wigStyle = new DASStyle();
    wigStyle.glyph = 'BOX';
    wigStyle.FGCOLOR = 'black';
    wigStyle.BGCOLOR = 'orange'
    wigStyle.HEIGHT = 8;
    wigStyle.BUMP = true;
    wigStyle.LABEL = true;
    wigStyle.ZINDEX = 20;
    stylesheet.pushStyle({type: 'default'}, null, wigStyle);



    return callback(stylesheet);
}


EnsemblFeatureSource.prototype.getScales = function() {
    return [];
}

EnsemblFeatureSource.prototype.fetch = function(chr, min, max, scale, types, pool, callback) {
    var thisB = this;
    url = this.base + '/feature/region/' + this.species + '/' + chr + ':' + min + '-' + max;

    var filters = [];
    for (var ti = 0; ti < this.type.length; ++ti) {
        filters.push('feature=' + this.type[ti]);
    }
    filters.push('content-type=application/json');
    url = url + '?' + filters.join(';');

    var req = new XMLHttpRequest();
    req.onreadystatechange = function() {
    	if (req.readyState == 4) {
            thisB.busy--;
            thisB.notifyActivity();

    	    if (req.status >= 300) {
                var err = 'Error code ' + req.status;
                try {
                    var jr = JSON.parse(req.response);
                    if (jr.error) {
                        err = jr.error;
                    }
                } catch (ex) {};

    		    callback(err, null);
    	    } else {
        		var jf = JSON.parse(req.response);
        		var features = [];
        		for (fi = 0; fi < jf.length; ++fi) {
        		    var j = jf[fi];

        		    var notes = [];
        		    var f = new DASFeature();
        		    f.segment = chr;
        		    f.min = j['start'] | 0;
        		    f.max = j['end'] | 0;
        		    f.type = j.feature_type || 'unknown';
        		    f.id = j.ID;

                    if (j.Parent) {
                        var grp = new DASGroup();
                        grp.id = j.Parent;
                        f.groups = [grp];
                    }

                    if (j.strand) {
                        if (j.strand < 0) 
                            f.orientation = '-';
                        else if (j.strand > 0) 
                            f.orientation = '+';
                    }

                    if (j.consequence_type)
                        f.method = j.consequence_type;

                    if (j.alt_alleles) {
                        notes.push('Alleles=' + j.alt_alleles.join('/'));
                        if (j.alt_alleles.length > 1) {
                            if (j.alt_alleles[1].length != j.alt_alleles[0].length || j.alt_alleles[1] == '-') {
                                f.type = 'indel';
                            }
                        }
                    }
        		    
                    if (notes.length > 0) {
                        f.notes = notes;
                    }
        		    features.push(f);
        		}
        		callback(null, features);
    	    }
    	}
	
    };
    
    thisB.busy++;
    thisB.notifyActivity();

    req.open('GET', url, true);
    req.responseType = 'text';
    req.send('');
}

dalliance_registerSourceAdapterFactory('ensembl', function(source) {
    return {features: new EnsemblFeatureSource(source)};
});
/* -*- mode: javascript; c-basic-offset: 4; indent-tabs-mode: nil -*- */

// 
// Dalliance Genome Explorer
// (c) Thomas Down 2006-2013
//
// overlay.js: featuresources composed from multiple underlying sources
//

"use strict";

function OverlayFeatureSource(sources, opts) {
    this.sources = sources;
    this.opts = opts || {};
    this.activityListeners = [];
    this.readinessListeners = [];
    this.changeListeners = [];
    this.business = [];
    this.readiness = [];

    for (var i = 0; i < this.sources.length; ++i) {
        this.initN(i);
    }

    if (opts.merge == 'concat') {
        this.merge = OverlayFeatureSource_merge_concat;
    } else {
        this.merge = OverlayFeatureSource_merge_byKey;
    }
}

OverlayFeatureSource.prototype.initN = function(n) {
    var s = this.sources[n];
    var thisB = this;
    this.business[n] = 0;

    if (s.addActivityListener) {
        s.addActivityListener(function(b) {
            thisB.business[n] = b;
            thisB.notifyActivity();
        });
    }
    if (s.addChangeListener) {
        s.addChangeListener(function() {
            thisB.notifyChange();
        });
    }
    if (s.addReadinessListener) {
        s.addReadinessListener(function(r) {
            thisB.readiness[n] = r;
            thisB.notifyReadiness();
        });
    }
}

OverlayFeatureSource.prototype.addReadinessListener = function(l) {
    this.readinessListeners.push(l);
    this.notifyReadinessListener(l);
}

OverlayFeatureSource.prototype.notifyReadiness = function() {
    for (var i = 0; i < this.readinessListeners.length; ++i) {
        this.notifyReadinessListener(this.readinessListeners[i]);
    }
}

OverlayFeatureSource.prototype.notifyReadinessListener = function(l) {
    var r = null;
    for (var i = 0; i < this.readiness.length; ++i) {
        if (this.readiness[i] != null) {
            r = this.readiness[i]; break;
        }
    }
    try {
        l(r);
    } catch (e) {
        console.log(e);
    }
}

OverlayFeatureSource.prototype.addActivityListener = function(l) {
    this.activityListeners.push(l);
}

OverlayFeatureSource.prototype.notifyActivity = function() {
    var busy = 0;
    for (var i = 0; i < this.business.length; ++i) {
        busy += this.business[i];
    }

    for (var li = 0; li < this.activityListeners.length; ++li) {
        try {
            this.activityListeners[li](busy);
        } catch (e) {
            console.log(e);
        }
    }
}

OverlayFeatureSource.prototype.addChangeListener = function(listener) {
    this.changeListeners.push(listener);
}

OverlayFeatureSource.prototype.notifyChange = function() {
    for (var li = 0; li < this.changeListeners.length; ++li) {
        try {
            this.changeListeners[li](this.busy);
        } catch (e) {
            console.log(e);
        }
    }
}

OverlayFeatureSource.prototype.getScales = function() {
    return this.sources[0].getScales();
}

OverlayFeatureSource.prototype.getStyleSheet = function(callback) {
    return this.sources[0].getStyleSheet(callback);
}

OverlayFeatureSource.prototype.capabilities = function() {
    var caps = {};
    var s0 = this.sources[0];
    if (s0.capabilities) 
        caps = shallowCopy(s0.capabilities());

    for (var i = 1; i < this.sources.length; ++i) {
        var si = this.sources[i];
        if (si.capabilities) {
            var co = si.capabilities();
            if (co.search) {
                caps.search = co.search;
            }
        }
    }

    return caps;
}

OverlayFeatureSource.prototype.search = function(query, callback) {
    for (var i = 0; i < this.sources.length; ++i) {
        if (sourceAdapterIsCapable(this.sources[i], 'search')) {
            return this.sources[i].search(query, callback);
        }
    }
}

OverlayFeatureSource.prototype.fetch = function(chr, min, max, scale, types, pool, callback) {
    var baton = new OverlayBaton(this, callback, this.sources.length);
    for (var si = 0; si < this.sources.length; ++si) {
	this.fetchN(baton, si, chr, min, max, scale, types, pool);
    }
}

OverlayFeatureSource.prototype.fetchN = function(baton, si, chr, min, max, scale, types, pool) {
    this.sources[si].fetch(chr, min, max, scale, types, pool, function(status, features, scale) {
	return baton.completed(si, status, features, scale);
    });
}

OverlayFeatureSource.prototype.quantFindNextFeature = function(chr, pos, dir, threshold, callback) {
    return this.sources[0].quantFindNextFeature(chr, pos, dir, threshold, callback);
}

OverlayFeatureSource.prototype.findNextFeature = function(chr, pos, dir, callback) {
    return this.sources[0].findNextFeature(chr, pos, dir, callback);
}

function OverlayBaton(source, callback, count) {
    this.source = source;
    this.callback = callback;
    this.count = count;

    this.returnCount = 0;
    this.statusCount = 0;
    this.returns = [];
    this.features = []
    this.statuses = [];
    this.scale = null;
}

OverlayBaton.prototype.completed = function(index, status, features, scale) {
    if (this.scale == null || index == 0) 
	   this.scale = scale;

    if (this.returns[index])
	   throw 'Multiple returns for source ' + index;

    this.returns[index] = true;
    this.returnCount++;

    this.features[index] = features;

    if (status) {
    	this.statuses[index] = status;
    	this.statusCount++;
    }


    if (this.returnCount == this.count) {
    	if (this.statusCount > 0) {
    	    var message = '';
    	    for (var si = 0; si < this.count; ++si) {
        		var s = this.statuses[si];
        		if (s) {
        		    if (message.length > 0) 
        			message += ', ';
        		    message += s;
        		}
    	    }
    	    return this.callback(message, null, this.scale);
    	} else {
    	    this.callback(null, this.source.merge(this.features), this.scale);
    	}
    }
}

OverlayFeatureSource.prototype.getDefaultFIPs = function(callback) {
    for (var si = 0; si < this.sources.length; ++si) {
        var s = this.sources[si];
        if (s.getDefaultFIPs)
            s.getDefaultFIPs(callback);
    }
}

OverlayFeatureSource.prototype.keyForFeature = function(feature) {
    return '' + feature.min + '..' + feature.max;
}

function OverlayFeatureSource_merge_byKey(featureSets) {
    var omaps = [];

    for (var fsi = 1; fsi < featureSets.length; ++fsi) {
        var om = {};
        var of = featureSets[fsi];
        for (var fi = 0; fi < of.length; ++fi) {
    	   om[this.keyForFeature(of[fi])] = of[fi];
        }
        omaps.push(om);
    }


    var mf = [];
    var fl = featureSets[0];
    for (var fi = 0; fi < fl.length; ++fi) {
    	var f = fl[fi];

        for (var oi = 0; oi < omaps.length; ++oi) {
            var om = omaps[oi];
        	of = om[this.keyForFeature(f)]
        	if (of) {
                for (var k in of) {
                    if (k === 'score') {
                        f.score2 = of.score;
                    } else if (k === 'min' || k === 'max' || k === 'segment' || k === '_cachedStyle') {
                        // do nothing
                    } else {
                        f[k] = of[k];
                    }
                }
        	}
        }
    	mf.push(f);
    }
    return mf;
}

function OverlayFeatureSource_merge_concat(featureSets) {
    var features = [];
    for (var fsi = 0; fsi < featureSets.length; ++fsi) {
        var fs = featureSets[fsi];
        var name = this.sources[fsi].name;
        for (var fi = 0; fi < fs.length; ++fi) {
            var f = fs[fi];
            f.method = name;
            features.push(f);
        }
    }
    return features;
}
/* -*- mode: javascript; c-basic-offset: 4; indent-tabs-mode: nil -*- */

// 
// Dalliance Genome Explorer
// (c) Thomas Down 2006-2013
//
// tier-actions.js
//

"use strict";

Browser.prototype.mergeSelectedTiers = function() {
    var sources = [];
    var styles = [];

    for (var sti = 0; sti < this.selectedTiers.length; ++sti) {
        var tier = this.tiers[this.selectedTiers[sti]];
	    sources.push(shallowCopy(tier.dasSource));
        var ss = tier.stylesheet.styles;
        for (var si = 0; si < ss.length; ++si) {
            var sh = ss[si];
            var nsh = shallowCopy(sh);
            nsh.method = tier.dasSource.name.replace(/[()+*?]/g, '\\$&');
            nsh._methodRE = null;
            nsh.style = shallowCopy(sh.style);
            if (nsh.style.ZINDEX === undefined)
                nsh.style.ZINDEX = sti;

            if (tier.forceMin) {
                nsh.style.MIN = tier.forceMin;
            }
            if (tier.forceMax) {
                nsh.style.MAX = tier.forceMax;
            }

            styles.push(nsh);
        }
    }
    
    this.addTier(
	{name: 'Merged',
	 merge: 'concat',
	 overlay: sources,
	 noDownsample: true,
     style: styles});

    this.setSelectedTier(this.tiers.length - 1);
}
/* -*- mode: javascript; c-basic-offset: 4; indent-tabs-mode: nil -*- */

// 
// Dalliance Genome Explorer
// (c) Thomas Down 2006-2011
//
// bin.js general binary data support
//

"use strict";

var REGION_PATTERN = /([\d+,\w,\.,\_,\-]+):([0-9,]+)([KkMmGg])?([\-,\,.]+([0-9,]+)([KkMmGg])?)?/;

function parseLocCardinal(n, m) {
    var i = n.replace(/,/g, '');
    if (m === 'k' || m === 'K') {
        return i * 1000;
    } else if (m == 'm' || m === 'M') {
        return i * 1000000;
    } else {
        return i;
    }
}

Browser.prototype.search = function(g, statusCallback) {
    var thisB = this;
    var m = REGION_PATTERN.exec(g);

    if (m) {
        var chr = m[1], start, end;
        if (m[5]) {
            start = parseLocCardinal(m[2],  m[3]);
            end = parseLocCardinal(m[5], m[6]);
        } else {
            var width = this.viewEnd - this.viewStart + 1;
            start = (parseLocCardinal(m[2], m[3]) - (width/2))|0;
            end = start + width - 1;
        }
        this.setLocation(chr, start, end, statusCallback);
    } else {
        if (!g || g.length == 0) {
            return false;
        }

        var searchCount = 0;
        var foundLatch = false;

        var searchCallback = function(found, err) {
            --searchCount;
            if (err) {
                return statusCallback(err);
            }

            if (!found) found = [];
            var min = 500000000, max = -100000000;
            var nchr = null;
            for (var fi = 0; fi < found.length; ++fi) {
                var f = found[fi];
            
                if (nchr == null) {
                    nchr = f.segment;
                }
                min = Math.min(min, f.min);
                max = Math.max(max, f.max);
            }

            if (!nchr) {
                if (searchCount == 0 && !foundLatch)
                    return statusCallback("no match for '" + g + "'");
            } else {
                foundLatch = true;
                thisB.highlightRegion(nchr, min, max);
            
                var padding = Math.max(2500, (0.3 * (max - min + 1))|0);
                thisB.setLocation(nchr, min - padding, max + padding, statusCallback);
            }
        }

        var doTrixSearch = function(tier, trix) {
            trix.lookup(g, function(result, status) {
                if (result == null || result.length < 2) {
                    return tier.featureSource.search(g, searchCallback);
                } else {
                    var hit = result[1].split(',')[0];
                    return tier.featureSource.search(hit, searchCallback);
                }
            });
        }

        if (this.searchEndpoint) {
            searchCount = 1;
            return this.doDasSearch(thisB.searchEndpoint, g, searchCallback);
        }

        for (var ti = 0; ti < this.tiers.length; ++ti) {
            (function(tier) {
                if (sourceAdapterIsCapable(tier.featureSource, 'search')) {
                    if (tier.dasSource.trixURI) {
                        ++searchCount;
                        if (tier.trix) {
                            doTrixSearch(tier, tier.trix);
                        } else {
                            connectTrix(new URLFetchable(tier.dasSource.trixURI), new URLFetchable(tier.dasSource.trixURI + 'x'), function(trix) {
                                tier.trix = trix;
                                doTrixSearch(tier, trix);
                            });
                        }
                    } else {
                        ++searchCount;
                        tier.featureSource.search(g, searchCallback);
                    }
                } else if (tier.dasSource.provides_search) {
                    ++searchCount;
                    thisB.doDasSearch(tier.dasSource, g, searchCallback);
                }
            })(this.tiers[ti]);
        }
    }
}

Browser.prototype.doDasSearch = function(source, g, searchCallback) {
    var thisB = this;
    source.features(null, {group: g, type: 'transcript'}, function(found) {
        if (!found) found = [];
        var min = 500000000, max = -100000000;
        var nchr = null;

        var found2 = [];
        for (var fi = 0; fi < found.length; ++fi) {
            var f = found[fi];
            
            if (f.label.toLowerCase() != g.toLowerCase()) {
                // ...because Dazzle can return spurious overlapping features.
                continue;
            }
            found2.push(f);
        }

        return searchCallback(found2);
    }, false);
}/* -*- mode: javascript; c-basic-offset: 4; indent-tabs-mode: nil -*- */

// 
// Dalliance Genome Explorer
// (c) Thomas Down 2006-2013
//
// tier-edit.js
//

"use strict";

var __dalliance_smallGlyphs = {DOT: true, EX: true, STAR: true, SQUARE: true, CROSS: true, TRIANGLE: true, PLIMSOLL: true}

Browser.prototype.openTierPanel = function(tier) {
    var b = this;

    if (this.uiMode === 'tier' && this.manipulatingTier === tier) {
        this.hideToolPanel();
        this.setUiMode('none');
    } else {
        var setStyleColors = function(style) {
            if (style.BGGRAD)
                return;

            if (numColors == 1) {
                if (style.glyph == 'LINEPLOT' || __dalliance_smallGlyphs[style.glyph]) {
                    style.FGCOLOR = tierColorField.value;
                } else {
                    style.BGCOLOR = tierColorField.value;
                }
                style.COLOR1 = style.COLOR2 = style.COLOR3 = null;
            } else {
                style.COLOR1 = tierColorField.value;
                style.COLOR2 = tierColorField2.value;
                if (numColors > 2) {
                    style.COLOR3 = tierColorField3.value;
                } else {
                    style.COLOR3 = null;
                }
            }
            style._gradient = null;
        }

        var mutateStylesheet = function(visitor) {
            var nss = copyStylesheet(tier.stylesheet);
            var ssScale = tier.browser.zoomForCurrentScale();

            for (var i = 0; i < nss.styles.length; ++i) {
                var sh = nss.styles[i];
                if (sh.zoom && sh.zoom != ssScale) {
                    continue;
                }

                visitor(sh.style);
            }

            return nss;
        }

        var changeColor = function(ev) {
            tier.mergeConfig({stylesheet: mutateStylesheet(setStyleColors)});
        }
        
        this.manipulatingTier = tier;

        var tierForm = makeElement('div', null, {className: 'tier-edit'});
        if (tier.dasSource.mapping) {
            var coords = this.chains[tier.dasSource.mapping].coords;
            tierForm.appendChild(makeElement('div', 'Mapped from ' + coords.auth + coords.version , null, 
                {background: 'gray', paddingBottom: '5px', marginBottom: '5px', textAlign: 'center'}));
        }
        var semanticBanner = makeElement('div', 'Editing styles for current zoom level', null,
                {background: 'gray', paddingBottom: '5px', marginBottom: '5px', textAlign: 'center', display: 'none'});
        tierForm.appendChild(semanticBanner);

        var tierNameField = makeElement('input', null, {type: 'text'});
        var tierPinnedToggle = makeElement('input', null, {type: 'checkbox'});

        var glyphField = makeElement('select');
        glyphField.appendChild(makeElement('option', 'Histogram', {value: 'HISTOGRAM'}));
        glyphField.appendChild(makeElement('option', 'Line Plot', {value: 'LINEPLOT'}));
        glyphField.appendChild(makeElement('option', 'Ribbon', {value: 'GRADIENT'}));
        glyphField.appendChild(makeElement('option', 'Scatter', {value: 'SCATTER'}));

        var tierColorField = makeElement('input', null, {type: 'text', value: '#dd00dd'});
        var tierColorField2 = makeElement('input', null, {type: 'text', value: '#dd00dd'});
        var tierColorField3 = makeElement('input', null, {type: 'text', value: '#dd00dd'});
        try {
            tierColorField.type = tierColorField2.type = tierColorField3.type = 'color';
        } catch (e) {
            // IE throws if attempt to set type to 'color'.
        }

        var tierColorFields = [tierColorField, tierColorField2, tierColorField3];
        var colorListPlus = makeElement('i', null, {className: 'fa fa-plus-circle'});
        var colorListMinus = makeElement('i', null, {className: 'fa fa-minus-circle'});
        var numColors = 1;
        var colorListElement = makeElement('td', tierColorFields);
        var setNumColors = function(n) {
            numColors = n;
            for (var i = 0; i < n; ++i) 
                tierColorFields[i].style.display = 'block';
            for (var i = n; i < tierColorFields.length; ++i)
                tierColorFields[i].style.display = 'none';
        }
        colorListPlus.addEventListener('click', function(ev) {
            if (numColors < 3) {
                setNumColors(numColors + 1);
                changeColor(null);
            }
        }, false);
        colorListMinus.addEventListener('click', function(ev) {
            if (numColors > 1) {
                setNumColors(numColors - 1);
                changeColor(null);
            }
        }, false);

        var tierMinField = makeElement('input', null, {type: 'text', value: '0.0'});
        var tierMaxField = makeElement('input', null, {type: 'text', value: '10.0'});
        var tierMinToggle = makeElement('input', null, {type: 'checkbox'});
        var tierMaxToggle = makeElement('input', null, {type: 'checkbox'});

        var quantLeapToggle = makeElement('input', null, {type: 'checkbox', checked: tier.quantLeapThreshold !== undefined});
        var quantLeapThreshField = makeElement('input', null, {type: 'text', value: tier.quantLeapThreshold, disabled: !quantLeapToggle.checked});

        var tierHeightField = makeElement('input', null, {type: 'text', value: '50'});

        var bumpToggle = makeElement('input', null, {type: 'checkbox'});
        var labelToggle = makeElement('input', null, {type: 'checkbox'});

        var mainStyle = null;
        if (tier.stylesheet.styles.length > 0) {
            var s = mainStyle = tier.stylesheet.styles[0].style;
        }

        var refresh = function() {
            if (typeof tier.config.name === 'string')
                tierNameField.value = tier.config.name;
            else 
                tierNameField.value = tier.dasSource.name;

            tierPinnedToggle.checked = tier.pinned;

            if (tier.forceHeight) {
                tierHeightField.value = '' + tier.forceHeight;
            } else if (mainStyle && mainStyle.HEIGHT) {
                tierHeightField.value = '' + mainStyle.HEIGHT;
            }

            if (typeof tier.quantLeapThreshold == 'number') {
                quantLeapToggle.checked = true;
                quantLeapThreshField.disabled = false;
                if (parseFloat(quantLeapThreshField.value) != tier.quantLeapThreshold)
                    quantLeapThreshField.value = tier.quantLeapThreshold;
            } else {
                quantLeapToggle.checked = false;
                quantLeapThreshField.disabled = true;
            }

            if (tier.stylesheet.styles.length > 0) {
                var s = mainStyle = tier.stylesheet.styles[0].style;
                var isQuantitative=false, isSimpleQuantitative = false;
                var ssScale = tier.browser.zoomForCurrentScale();
                var activeStyleCount = 0;

                for (var si = 0; si < tier.stylesheet.styles.length; ++si) {
                    var sh = tier.stylesheet.styles[si];  
                    if (sh.zoom && sh.zoom != ssScale) {
                        continue;
                    }
                    ++activeStyleCount;

                    var ss = tier.stylesheet.styles[si].style;
                    if (ss.glyph == 'LINEPLOT' || ss.glyph == 'HISTOGRAM' || ss.glyph == 'GRADIENT' || isDasBooleanTrue(ss.SCATTER)) {
                        if (!isQuantitative)
                            s = mainStyle = ss;
                        isQuantitative = true;
                    }
                }

                semanticBanner.style.display = (activeStyleCount == tier.stylesheet.styles.length) ? 'none' : 'block';

                isSimpleQuantitative = isQuantitative && activeStyleCount == 1;
                var isGradient = s.COLOR2 || s.BGGRAD;

                if (isQuantitative) {
                    minRow.style.display = 'table-row';
                    maxRow.style.display = 'table-row';
                    bumpRow.style.display = 'none';
                    labelRow.style.display = 'none';
                } else {
                    minRow.style.display = 'none';
                    maxRow.style.display = 'none';
                    bumpRow.style.display = 'table-row';
                    bumpToggle.checked = isDasBooleanTrue(mainStyle.BUMP);
                    labelRow.style.display = 'table-row';
                    labelToggle.checked = isDasBooleanTrue(mainStyle.LABEL);
                }

                if (isSimpleQuantitative) {
                    styleRow.style.display = 'table-row';
                    colorRow.style.display = 'table-row';
                } else {
                    styleRow.style.display = 'none';
                    colorRow.style.display = 'none';

                }

                var numColors = 1;
                if (s.COLOR1) {
                    tierColorField.value = dasColourForName(s.COLOR1).toHexString();
                    if (s.COLOR2) {
                        tierColorField2.value = dasColourForName(s.COLOR2).toHexString();
                        if (s.COLOR3) {
                            tierColorField3.value = dasColourForName(s.COLOR3).toHexString();
                            numColors = 3;
                        } else {
                            numColors = 2;
                        }
                    }
                } else {
                    if (s.glyph == 'LINEPLOT' || s.glyph == 'DOT' && s.FGCOLOR) {
                        tierColorField.value = dasColourForName(s.FGCOLOR).toHexString();
                    } else if (s.BGCOLOR) {
                        tierColorField.value = dasColourForName(s.BGCOLOR).toHexString();
                    }
                } 
                setNumColors(numColors);

                if (isDasBooleanTrue(s.SCATTER)) {
                    glyphField.value = 'SCATTER';
                } else {
                    glyphField.value = s.glyph;
                } 

                var setMinValue, setMaxValue;
                if (s.MIN !== undefined) {
                    var x = parseFloat(s.MIN);
                    if (!isNaN(x))
                        setMinValue = x;
                }
                if (!tier.forceMinDynamic && (s.MIN !== undefined || tier.forceMin !== undefined)) {
                    tierMinToggle.checked = true;
                    tierMinField.disabled = false;
                } else {
                    tierMinToggle.checked = false;
                    tierMinField.disabled = true;
                }

                if (s.MAX !== undefined) {
                    var x = parseFloat(s.MAX)
                    if (!isNaN(x))
                        setMaxValue = x;
                }
                if (!tier.forceMaxDynamic && (s.MAX !== undefined || tier.forceMax !== undefined)) {
                    tierMaxToggle.checked = true;
                    tierMaxField.disabled = false;
                } else {
                    tierMaxToggle.checked = false;
                    tierMaxField.disabled = true;
                }

                if (tier.forceMin != undefined) {
                    setMinValue = tier.forceMin;
                }
                if (tier.forceMax != undefined) {
                    setMaxValue = tier.forceMax;
                }
                if (typeof(setMinValue) == 'number' && setMinValue != parseFloat(tierMinField.value)) {
                    tierMinField.value = setMinValue;
                }
                if (typeof(setMaxValue) == 'number' && setMaxValue != parseFloat(tierMaxField.value)) {
                    tierMaxField.value = setMaxValue;
                }

                var seqStyle = getSeqStyle(tier.stylesheet);
                if (seqStyle) {
                    seqMismatchRow.style.display = 'table-row';
                    seqMismatchToggle.checked = (seqStyle.__SEQCOLOR === 'mismatch');
                    seqInsertRow.style.display = 'table-row';
                    seqInsertToggle.checked =  isDasBooleanTrue(seqStyle.__INSERTIONS);
                } else {
                    seqMismatchRow.style.display = 'none';
                    seqInsertRow.style.display = 'none';
                }
            }

            if (isQuantitative && sourceAdapterIsCapable(tier.featureSource, 'quantLeap'))
                quantLeapRow.style.display = 'table-row';
            else 
                quantLeapRow.style.display = 'none';
        }

        var seqMismatchToggle = makeElement('input', null, {type: 'checkbox'});
        var seqMismatchRow = makeElement('tr',
            [makeElement('th', 'Color mismatches'),
             makeElement('td', seqMismatchToggle)]);
        seqMismatchToggle.addEventListener('change', function(ev) {
            var nss = copyStylesheet(tier.stylesheet);
            var seqStyle = getSeqStyle(nss);
            seqStyle.__SEQCOLOR = seqMismatchToggle.checked ? 'mismatch' : 'base';
            tier.mergeConfig({stylesheet: nss});
        });

        var seqInsertToggle = makeElement('input', null, {type: 'checkbox'});
        var seqInsertRow = makeElement('tr',
            [makeElement('th', 'Show insertions'),
             makeElement('td', seqInsertToggle)]);
        seqInsertToggle.addEventListener('change', function(ev) {
            var nss = copyStylesheet(tier.stylesheet);
            var seqStyle = getSeqStyle(nss);
            seqStyle.__INSERTIONS = seqInsertToggle.checked ? 'yes' : 'no';
            tier.mergeConfig({stylesheet: nss});
        });

        var styleRow = makeElement('tr',
                [makeElement('th', 'Style'),
                 makeElement('td', glyphField)]);
        var colorRow = makeElement('tr',
                [makeElement('th', ['Colour(s)', colorListPlus, colorListMinus]),
                 colorListElement]);
        var minRow = makeElement('tr',
                [makeElement('th', 'Min value'),
                 makeElement('td', [tierMinToggle, ' ', tierMinField])]);
        var maxRow = makeElement('tr',
                [makeElement('th', 'Max value'),
                 makeElement('td', [tierMaxToggle, ' ', tierMaxField])]);
        var quantLeapRow = 
             makeElement('tr',
                [makeElement('th', 'Threshold leap:'),
                 makeElement('td', [quantLeapToggle, ' ', quantLeapThreshField])]);
        var bumpRow = makeElement('tr',
                [makeElement('th', 'Bump overlaps'),
                 makeElement('td', bumpToggle)]);
        var labelRow = makeElement('tr',
                [makeElement('th', 'Label features'),
                 makeElement('td', labelToggle)]);


        var tierTable = makeElement('table',
            [makeElement('tr',
                [makeElement('th', 'Name', {}, {width: '150px', textAlign: 'right'}),
                 tierNameField]),

             makeElement('tr',
                [makeElement('th', 'Pin to top'),
                 tierPinnedToggle]),

             makeElement('tr',
                [makeElement('th', 'Height'),
                 makeElement('td', tierHeightField)]),

            styleRow,
            colorRow,
            minRow,
            maxRow,
            quantLeapRow,
            bumpRow,
            labelRow,
            seqMismatchRow,
            seqInsertRow
             ]);


        refresh();

        tierForm.appendChild(tierTable);

        var resetButton = makeElement('button', 'Reset track', {className: 'btn'}, {marginLeft: 'auto', marginRight: 'auto', display: 'block'});
        resetButton.addEventListener('click', function(ev) {
            tier.setConfig({});
        }, false);
        tierForm.appendChild(resetButton);

        tierNameField.addEventListener('input', function(ev) {
            tier.mergeConfig({name: tierNameField.value});
        }, false);

        tierPinnedToggle.addEventListener('change', function(ev) {
            tier.mergeConfig({pinned: tierPinnedToggle.checked});
        }, false);

        for (var ci = 0; ci < tierColorFields.length; ++ci) {
            tierColorFields[ci].addEventListener('change', changeColor, false);
        }

        glyphField.addEventListener('change', function(ev) {
            var nss = mutateStylesheet(function(ts) {
                if (glyphField.value === 'SCATTER') {
                    ts.SCATTER = true;
                    ts.glyph = 'DOT';
                    ts.SIZE = '3';
                } else {
                    ts.glyph = glyphField.value;
                    ts.SCATTER = undefined;
                }
                setStyleColors(ts);
            });
            tier.mergeConfig({stylesheet: nss});
        }, false);

        tierMinToggle.addEventListener('change', function(ev) {
            var conf = {forceMinDynamic: !tierMinToggle.checked};
            tierMinField.disabled = !tierMinToggle.checked;
            var x = parseFloat(tierMinField.value);
            if (tierMinToggle.checked && typeof(x) == 'number' && !isNaN(x))
                conf.forceMin = parseFloat(x);
            tier.mergeConfig(conf);
        });
        tierMinField.addEventListener('input', function(ev) {
            var x = parseFloat(tierMinField.value);
            if (typeof(x) == 'number' && !isNaN(x))
                tier.mergeConfig({forceMin: x});
        }, false);

        tierMaxToggle.addEventListener('change', function(ev) {
            var conf = {forceMaxDynamic: !tierMaxToggle.checked};
            tierMaxField.disabled = !tierMaxToggle.checked;
            var x = parseFloat(tierMaxField.value);
            if (tierMaxToggle.checked && typeof(x) == 'number' && !isNaN(x))
                conf.forceMax = parseFloat(x);
            tier.mergeConfig(conf);
        });
        tierMaxField.addEventListener('input', function(ev) {
            var x = parseFloat(tierMaxField.value);
            if (typeof(x) == 'number' && !isNaN(x))
                tier.mergeConfig({forceMax: x});
        }, false);

        tierHeightField.addEventListener('input', function(ev) {
            var x = parseFloat(tierHeightField.value);
            if (typeof(x) == 'number' && !isNaN(x))
                tier.mergeConfig({height: Math.min(500, x|0)});
        }, false);

        var updateQuant = function() {
            quantLeapThreshField.disabled = !quantLeapToggle.checked;
            if (quantLeapToggle.checked) {
                var x = parseFloat(quantLeapThreshField.value);
                if (typeof(x) == 'number' && !isNaN(x)) {
                    tier.mergeConfig({quantLeapThreshold: parseFloat(quantLeapThreshField.value)});
                }
            } else {
                tier.mergeConfig({quantLeapThreshold: null});
            }
        }
        quantLeapToggle.addEventListener('change', function(ev) {
            updateQuant();
        }, false);
        quantLeapThreshField.addEventListener('input', function(ev) {
            updateQuant();
        }, false);

        labelToggle.addEventListener('change', function(ev) {
            var nss = mutateStylesheet(function(style) {
                style.LABEL = labelToggle.checked ? 'yes' : 'no';
            });
            tier.mergeConfig({stylesheet: nss});
        }, false);
        bumpToggle.addEventListener('change', function(ev) {
            var nss = mutateStylesheet(function(style) {
                style.BUMP = bumpToggle.checked ? 'yes' : 'no';
            });
            tier.mergeConfig({stylesheet: nss});
        }, false);


        this.showToolPanel(tierForm);
        this.setUiMode('tier');

        tier.addTierListener(refresh);

        var currentScale = tier.browser.scale;
        tier.browser.addViewListener(function() {
            if (tier.browser.scale != currentScale) {
                currentScale = tier.browser.scale;
                refresh();
            }
        });
    }
}

function getSeqStyle(stylesheet) {
    for (var si = 0; si < stylesheet.styles.length; ++si) {
        var ss = stylesheet.styles[si].style;
        if (ss.glyph === '__SEQUENCE') {
            return ss;
        }
    }
}

function copyStylesheet(ss) {
    var nss = shallowCopy(ss);
    nss.styles = [];
    for (var si = 0; si < ss.styles.length; ++si) {
        var sh = nss.styles[si] = shallowCopy(ss.styles[si]);
        sh._methodRE = sh._labelRE = sh._typeRE = undefined;
        sh.style = shallowCopy(sh.style);
        sh.style.id = undefined;
    }
    return nss;
}/* -*- mode: javascript; c-basic-offset: 4; indent-tabs-mode: nil -*- */

// 
// Dalliance Genome Explorer
// (c) Thomas Down 2006-2010
//
// trix.js: UCSC-style free text indices
//

"use strict";

function connectTrix(ix, ixx, callback) {
    ixx.fetchAsText(function(ixxData) {
        if (!ixxData) 
            return callback(null, "Couldn't fetch index-index");

        var toks = ixxData.split(/(.+)([0-9A-F]{10})\n/);

        var keys = [];
        var offsets = [];
        for (var ti = 1; ti < toks.length; ti += 3) {
            keys.push(toks[ti]);
            offsets.push(parseInt(toks[ti+1], 16));
        }

        return callback(new TrixIndex(keys, offsets, ix));
    });
}

function TrixIndex(keys, offsets, ix) {
    this.keys = keys;
    this.offsets = offsets;
    this.ix = ix;
}

TrixIndex.prototype.lookup = function(query, callback) {
    var ixslice;

    for (var i = 0; i < this.keys.length; ++i) {
        if (query.localeCompare(this.keys[i]) < 0) {
            ixslice = this.ix.slice(this.offsets[i - 1], this.offsets[i] - this.offsets[i - 1]);
            break;
        }
    }

    if (!ixslice) {
        ixslice = this.ix.slice(this.offsets[this.offsets.length - 1]);
    }

    ixslice.fetchAsText(function(ist) {
        var lines = ist.split('\n');
        for (var li = 0; li < lines.length; ++li) {
            if (lines[li].indexOf(query.toLowerCase() + ' ') == 0) {
                return callback(lines[li].split(' '));
            }
        }
        return callback(null);
    });
}/* -*- mode: javascript; c-basic-offset: 4; indent-tabs-mode: nil -*- */

// 
// Dalliance Genome Explorer
// (c) Thomas Down 2006-2011
//
// tabix.js: basic support for tabix-indexed flatfiles
//

"use strict";

var TABIX_MAGIC = 0x01494254;

function TabixFile() {
}

function connectTabix(data, tbi, callback) {
    var tabix = new TabixFile();
    tabix.data = data;
    tabix.tbi = tbi;

    tabix.tbi.fetch(function(header) {   // Do we really need to fetch the whole thing? :-(
        if (!header) {
            return dlog("Couldn't access Tabix");
        }

        var unchead = unbgzf(header, header.byteLength);
        var uncba = new Uint8Array(unchead);
        var magic = readInt(uncba, 0);
        if (magic != TABIX_MAGIC) {
            return callback(null, 'Not a tabix index');
        }

        var nref = readInt(uncba, 4);
        tabix.format = readInt(uncba, 8);
        tabix.colSeq = readInt(uncba, 12);
        tabix.colStart = readInt(uncba, 16);
        tabix.colEnd = readInt(uncba, 20);
        tabix.meta = readInt(uncba, 24);
        tabix.skip = readInt(uncba, 28);
        var nameLength = readInt(uncba, 32);

        tabix.indices = [];

        var p = 36;
        tabix.chrToIndex = {};
        tabix.indexToChr = [];
        for (var i = 0; i < nref; ++i) {
            var name = ''

            while (true) {
                var ch = uncba[p++];
                if (ch == 0)
                    break;

                name += String.fromCharCode(ch);
            }

            tabix.chrToIndex[name] = i;
            if (name.indexOf('chr') == 0) {
                tabix.chrToIndex[name.substring(3)] = i;
            } else {
                tabix.chrToIndex['chr' + name] = i;
            }
            tabix.indexToChr.push(name);
        }

        var minBlockIndex = 1000000000;
        for (var ref = 0; ref < nref; ++ref) {
            var blockStart = p;
            var nbin = readInt(uncba, p); p += 4;
            for (var b = 0; b < nbin; ++b) {
                var bin = readInt(uncba, p);
                var nchnk = readInt(uncba, p+4);
                p += 8 + (nchnk * 16);
            }
            var nintv = readInt(uncba, p); p += 4;
            
            var q = p;
            for (var i = 0; i < nintv; ++i) {
                var v = readVob(uncba, q); q += 8;
                if (v) {
                    var bi = v.block;
                    if (v.offset > 0)
                        bi += 65536;

                    if (bi < minBlockIndex)
                        minBlockIndex = bi;
                    break;
                }
            }
            p += (nintv * 8);


            var ub = uncba;
            if (nbin > 0) {
                tabix.indices[ref] = new Uint8Array(unchead, blockStart, p - blockStart);
            }                     
        }

        tabix.headerMax = minBlockIndex;

        callback(tabix);
    });
}

// Copy-paste from BamFile

TabixFile.prototype.blocksForRange = function(refId, min, max) {
    var index = this.indices[refId];
    if (!index) {
        return [];
    }

    var intBinsL = reg2bins(min, max);
    var intBins = [];
    for (var i = 0; i < intBinsL.length; ++i) {
        intBins[intBinsL[i]] = true;
    }
    var leafChunks = [], otherChunks = [];

    var nbin = readInt(index, 0);
    var p = 4;
    for (var b = 0; b < nbin; ++b) {
        var bin = readInt(index, p);
        var nchnk = readInt(index, p+4);
        p += 8;
        if (intBins[bin]) {
            for (var c = 0; c < nchnk; ++c) {
                var cs = readVob(index, p);
                var ce = readVob(index, p + 8);
                (bin < 4681 ? otherChunks : leafChunks).push(new Chunk(cs, ce));
                p += 16;
            }
        } else {
            p +=  (nchnk * 16);
        }
    }

    var nintv = readInt(index, p);
    var lowest = null;
    var minLin = Math.min(min>>14, nintv - 1), maxLin = Math.min(max>>14, nintv - 1);
    for (var i = minLin; i <= maxLin; ++i) {
        var lb =  readVob(index, p + 4 + (i * 8));
        if (!lb) {
            continue;
        }
        if (!lowest || lb.block < lowest.block || lb.offset < lowest.offset) {
            lowest = lb;
        }
    }
    
    var prunedOtherChunks = [];
    if (lowest != null) {
        for (var i = 0; i < otherChunks.length; ++i) {
            var chnk = otherChunks[i];
            if (chnk.maxv.block >= lowest.block && chnk.maxv.offset >= lowest.offset) {
                prunedOtherChunks.push(chnk);
            }
        }
    } 
    otherChunks = prunedOtherChunks;

    var intChunks = [];
    for (var i = 0; i < otherChunks.length; ++i) {
        intChunks.push(otherChunks[i]);
    }
    for (var i = 0; i < leafChunks.length; ++i) {
        intChunks.push(leafChunks[i]);
    }

    intChunks.sort(function(c0, c1) {
        var dif = c0.minv.block - c1.minv.block;
        if (dif != 0) {
            return dif;
        } else {
            return c0.minv.offset - c1.minv.offset;
        }
    });
    var mergedChunks = [];
    if (intChunks.length > 0) {
        var cur = intChunks[0];
        for (var i = 1; i < intChunks.length; ++i) {
            var nc = intChunks[i];
            if (nc.minv.block == cur.maxv.block /* && nc.minv.offset == cur.maxv.offset */) { // no point splitting mid-block
                cur = new Chunk(cur.minv, nc.maxv);
            } else {
                mergedChunks.push(cur);
                cur = nc;
            }
        }
        mergedChunks.push(cur);
    }

    return mergedChunks;
}

TabixFile.prototype.fetch = function(chr, min, max, callback) {
    var thisB = this;

    var chrId = this.chrToIndex[chr];
    if (chrId == undefined)
        return callback([]);

    var canonicalChr = this.indexToChr[chrId];

    var chunks;
    if (chrId === undefined) {
        chunks = [];
    } else {
        chunks = this.blocksForRange(chrId, min, max);
        if (!chunks) {
            callback(null, 'Error in index fetch');
        }
    }

    var records = [];
    var index = 0;
    var data;

    function tramp() {
        if (index >= chunks.length) {
            return callback(records);
        } else if (!data) {
            var c = chunks[index];
            var fetchMin = c.minv.block;
            var fetchMax = c.maxv.block + (1<<16); // *sigh*
            thisB.data.slice(fetchMin, fetchMax - fetchMin).fetch(function(r) {
                data = unbgzf(r, c.maxv.block - c.minv.block + 1);
                return tramp();
            });
        } else {
            var ba = new Uint8Array(data);
            thisB.readRecords(ba, chunks[index].minv.offset, records, min, max, canonicalChr);
            data = null;
            ++index;
            return tramp();
        }
    }
    tramp();
}

TabixFile.prototype.readRecords = function(ba, offset, sink, min, max, chr) {
   LINE_LOOP:
    while (true) {
        var line = '';
        while (offset < ba.length) {
            var ch = ba[offset++];
            if (ch == 10) {
                var toks = line.split('\t');

                if (toks[this.colSeq - 1] == chr) {
                    var fmin = parseInt(toks[this.colStart - 1]);
                    var fmax = fmin;
                    if (this.colEnd > 0)
                        fmax = parseInt(toks[this.colEnd - 1]);
                    if (this.format & 0x10000) ++fmin;

                    if (fmin <= max && fmax >= min)
                        sink.push(line);
                }
                continue LINE_LOOP;
            } else {
                line += String.fromCharCode(ch);
            }
        }
        return;
    }
}

TabixFile.prototype.fetchHeader = function(callback) {
    var self = this;
    var fetchPtr = 0, ptr = 0, line='';
    var lines = [];

    self.data.slice(0, self.headerMax).fetch(function(chnk) {
        if (!chnk) {
            return callback(null, "Fetch failed");
        }
        var ba = new Uint8Array(unbgzf(chnk, chnk.byteLength));
        var ptr = 0, line = '', lines = [];
        while (ptr < ba.length) {
            var ch = ba[ptr++]
            if (ch == 10) {
                if (line.charCodeAt(0) == self.meta) {
                    lines.push(line);
                    line = '';
                } else {
                    return callback(lines);
                }
            } else {
                line += String.fromCharCode(ch);
            }
        }
        callback(lines);
    });
}/* -*- mode: javascript; c-basic-offset: 4; indent-tabs-mode: nil -*- */

// 
// Dalliance Genome Explorer
// (c) Thomas Down 2006-2013
//
// tabix-source.js
//

"use strict";

function TabixFeatureSource(source) {
    FeatureSourceBase.call(this);
    this.readiness = 'Connecting';
    this.source = source;

    this.tabixHolder = new Awaited();
    var thisB = this;


    var parser = dalliance_makeParser(source.payload);
    if (!parser) {
        throw 'Unsuported tabix payload ' + source.payload;
    } else {
        this.parser = parser;
    }

    var data, index;
    if (this.source.blob) {
        data = new BlobFetchable(this.source.blob);
        index = new BlobFetchable(this.source.indexBlob);
    } else {
        data = new URLFetchable(this.source.uri, {credentials: this.source.credentials});
        index = new URLFetchable(this.source.indexURI || (this.source.uri + '.tbi'), {credentials: this.source.credentials});
    }
    connectTabix(data, index, function(tabix, err) {
        thisB.tabixHolder.provide(tabix);
        tabix.fetchHeader(function(lines, err) {
            if (lines) {
                var session = parser.createSession(function() { /* Null sink because we shouldn't get records */ });
                for (var li = 0; li < lines.length; ++li) {
                    session.parse(lines[li]);
                }
                session.flush();
            }
        });
        thisB.readiness = null
        thisB.notifyReadiness();
    });
}

TabixFeatureSource.prototype = Object.create(FeatureSourceBase.prototype);

TabixFeatureSource.prototype.fetch = function(chr, min, max, scale, types, pool, callback) {
    var thisB = this;
    
    thisB.busy++;
    thisB.notifyActivity();
    
    this.tabixHolder.await(function(tabix) {
        tabix.fetch(chr, min, max, function(records, error) {
            thisB.busy--;
            thisB.notifyActivity();

            var features = [];
            var session = thisB.parser.createSession(function(f) {features.push(f)});
            for (var ri = 0; ri < records.length; ++ri) {
                var f = session.parse(records[ri]);
            }
            session.flush();
            callback(null, features, 1000000000);
        });
    });
}


TabixFeatureSource.prototype.getStyleSheet = function(callback) {
    if (this.parser && this.parser.getStyleSheet)
        this.parser.getStyleSheet(callback)
}

TabixFeatureSource.prototype.getDefaultFIPs = function(callback) {
    if (this.parser && this.parser.getDefaultFIPs)
        this.parser.getDefaultFIPs(callback);
}


dalliance_registerSourceAdapterFactory('tabix', function(source) {
    return {features: new TabixFeatureSource(source)};
});/* -*- mode: javascript; c-basic-offset: 4; indent-tabs-mode: nil -*- */

// 
// Dalliance Genome Explorer
// (c) Thomas Down 2006-2014
//
// memstore.js
//

"use strict";

function MemStore() {
    this.featuresByChr = {};
    this.maxLength = 1;
    this.chrRing = null;
}

MemStore.prototype.addFeatures = function(features) {
    var dirty = {};
    for (var fi = 0; fi < features.length; ++fi) {
        var f = features[fi];
        var chr = f.segment || f.chr;
        var fa = this.featuresByChr[chr];
        if (!fa) {
            fa = [];
            this.featuresByChr[chr] = fa;
        }
        fa.push(f);
        dirty[chr] = true;

        var len = f.max - f.min + 1;
        if (len > this.maxLength)
            this.maxLength = len;
    }

    for (chr in dirty) {
        var fa = this.featuresByChr[chr];
        fa.sort(function(f1, f2) {
            var d = f1.min - f2.min;
            if (d != 0)
                return d;
            return f1.max - f2.max;
        });
    }
    this.chrRing = null;
}

MemStore.prototype._indexFor = function(fa, p) {
    var lb = 0, ub = fa.length;
    while (ub > lb) {
        var mid = ((lb + ub)/2)|0;
        if (mid >= fa.length)
            return fa.length;
        var mg = fa[mid];
        if (p < mg.min) {
            ub = mid;
        } else {
            lb = mid + 1;
        }
    }
    return ub;
}

MemStore.prototype.fetch = function(chr, min, max) {
    var fa = this.featuresByChr[chr];
    if (!fa) {
        if (chr.indexOf('chr') == 0)
            fa = this.featuresByChr[chr.substring(3)];
        else
            fa = this.featuresByChr['chr' + chr];
    }
    if (!fa)
        return [];

    var mini = Math.max(0, this._indexFor(fa, min - this.maxLength - 1));
    var maxi = Math.min(fa.length - 1, this._indexFor(fa, max));

    var res = [];
    for (var fi = mini; fi <= maxi; ++fi) {
        var f = fa[fi];
        if (f.min <= max && f.max >= min)
            res.push(f);
    }
    return res;
}

MemStore.prototype.findNextFeature = function(chr, pos, dir) {
    if (this.chrRing == null) {
        this.chrRing = [];
        for (var chr in this.featuresByChr) {
            this.chrRing.push(chr);
        }
        this.chrRing.sort();
    }

    var fa = this.featuresByChr[chr];
    if (!fa) {
        if (chr.indexOf('chr') == 0) {
            chr = chr.substring(3);
            fa = this.featuresByChr[chr];
        } else {
            chr = 'chr' + chr;
            fa = this.featuresByChr[chr];
        }
    }
    if (!fa)
        return null;

    var i = Math.max(0, Math.min(this._indexFor(fa, pos), fa.length - 1));
    if (dir > 0) {
        while (i < fa.length) {
            var f = fa[i++];
            if (f.min > pos)
                return f;
        }
        var chrInd = this.chrRing.indexOf(chr) + 1;
        if (chrInd >= this.chrRing.length)
            chrInd = 0;
        return this.findNextFeature(this.chrRing[chrInd], 0, dir);
    } else {
        while (i >= 0) {
            var f = fa[i--];
            if (f.max < pos)
                return f;
        }
        var chrInd = this.chrRing.indexOf(chr) - 1;
        if (chrInd < 0)
            chrInd = this.chrRing.length - 1;
        return this.findNextFeature(this.chrRing[chrInd], 10000000000, dir);
    }
}

function MemStoreFeatureSource(source) {
    this.source = source;
    FeatureSourceBase.call(this);
    this.storeHolder = new Awaited();
    this.parser = dalliance_makeParser(source.payload);
    if (!this.parser) {
        throw "Unsupported memstore payload: " + source.payload;
    }

    var thisB = this;
    this._load(function(resp, err) {
        if (!resp) {
            thisB.error = err || "No data"
            thisB.storeHolder.provide(null);
        } else {
            var store = new MemStore();
            var features = [];
            var lines = resp.split('\n');

            var session = thisB.parser.createSession(function(f) {features.push(f)});
            for (var li = 0; li < lines.length; ++li) {
                var line = lines[li];
                if (line.length > 0) {
                    session.parse(line);
                }
            }
            session.flush();

            store.addFeatures(features);

            thisB.storeHolder.provide(store);
        }
    });
}

MemStoreFeatureSource.prototype = Object.create(FeatureSourceBase.prototype);

MemStoreFeatureSource.prototype._load = function(callback) {
    if (this.source.blob) {
        var r = new FileReader();
        r.onloadend = function() {
            return callback(r.result, r.error);
        }
        r.readAsText(this.source.blob);
    } else {
        textXHR(this.source.uri, callback, {});
    }
}

MemStoreFeatureSource.prototype.fetch = function(chr, min, max, scale, types, pool, cnt) {
    var thisB = this;
    this.storeHolder.await(function(store) {
        if (store) {
            var f = store.fetch(chr, min, max);
            return cnt(null, f, 100000000);
        } else {
            return cnt(thisB.error)
        }
    });
}

MemStoreFeatureSource.prototype.getStyleSheet = function(callback) {
    if (this.parser && this.parser.getStyleSheet)
        this.parser.getStyleSheet(callback)
}

MemStoreFeatureSource.prototype.getDefaultFIPs = function(callback) {
    if (this.parser && this.parser.getDefaultFIPs)
        this.parser.getDefaultFIPs(callback);
}

MemStoreFeatureSource.prototype.getScales = function() {
    return 100000000;
}

MemStoreFeatureSource.prototype.findNextFeature = function(chr, pos, dir, callback) {
    var thisB = this;
    this.storeHolder.await(function(store) {
        if (store) {
            return callback(store.findNextFeature(chr, pos, dir));
        } else {
            return callback(null, thisB.error);
        }
    });
}


MemStoreFeatureSource.prototype.capabilities = function() {
    var caps = {leap: true};
    return caps;
}

dalliance_registerSourceAdapterFactory('memstore', function(source) {
    return {features: new MemStoreFeatureSource(source)};
});
/* -*- mode: javascript; c-basic-offset: 4; indent-tabs-mode: nil -*- */

// 
// Dalliance Genome Explorer
// (c) Thomas Down 2006-2014
//
// vcf.js
//

"use strict";

function VCFParser() {
    this.info = [];
}

var VCF_INFO_RE = /([^;=]+)(=([^;]+))?;?/;
var VCF_INFO_HEADER = /##INFO=<([^>]+)>/;
var VCF_INFO_HEADER_TOK = /([^,=]+)=([^,]+|"[^"]+"),?/

VCFParser.prototype.createSession = function(sink) {
    return new VCFParseSession(this, sink);
}

function VCFParseSession(parser, sink) {
    this.parser = parser;
    this.sink  = sink;
}

VCFParseSession.prototype.parse = function(line) {
    if (line.length == 0)
        return;
    if (line[0] == '#') {
        if (line.length > 1 && line[1] == '#') {
            var m = VCF_INFO_HEADER.exec(line);
            if (m) {
                var toks = m[1].split(VCF_INFO_HEADER_TOK);
                var id = null, desc = null;
                for (var ti = 0; ti < toks.length - 1; ti += 3) {
                    var key = toks[ti + 1];
                    var value = toks[ti + 2].replace(/"/g, '');
                    if (key == 'ID') {
                        id = value;
                    } else if (key == 'Description') {
                        desc = value;
                    }
                }
                if (id && desc) {
                    this.parser.info.push(
                        {id: id,
                         desc: desc}
                    );
                }
            }
            return;
        } else {
            return;
        }
    }

    var toks = line.split('\t');
    var f = new DASFeature();
    f.segment = toks[0];
    f.id = toks[2]
    f.refAllele = toks[3];
    f.altAlleles = toks[4].split(',');
    f.min = parseInt(toks[1]);
    f.max = f.min + f.refAllele.length - 1;

    var infoToks = toks[7].split(VCF_INFO_RE);
    f.info = {};
    for (var ti = 0; ti < infoToks.length; ti += 4) {
        f.info[infoToks[ti + 1]] = infoToks[ti + 3];
    }


    var alt = f.altAlleles[0];
    var ref = f.refAllele;
    if (alt.length > ref.length) {
        f.type = "insertion";
        if (alt.indexOf(ref) == 0) {
            f.insertion = alt.substr(ref.length);
            f.min += ref.length;
            f.max = f.min - 1; // Effectively "between" bases.
        } else {
            f.insertion = alt;
        }
    } else if (alt.length < ref.length) {
        f.type = "deletion";
    } else {
        f.type = 'substitution';
    }

    this.sink(f);
}

VCFParseSession.prototype.flush = function() {};

VCFParser.prototype.getStyleSheet = function(callback) {
    var stylesheet = new DASStylesheet();

    {
        var varStyle = new DASStyle();
        varStyle.glyph = '__INSERTION';
        varStyle.BUMP = 'yes';
        varStyle.LABEL = 'no';
        varStyle.FGCOLOR = 'rgb(50,80,255)';
        varStyle.BGCOLOR = '#888888';
        varStyle.STROKECOLOR = 'black';
        stylesheet.pushStyle({type: 'insertion'}, null, varStyle);
    }
    {
        var varStyle = new DASStyle();
        varStyle.glyph = 'PLIMSOLL';
        varStyle.BUMP = 'yes';
        varStyle.LABEL = 'no';
        varStyle.FGCOLOR = 'rgb(255, 60, 60)';
        varStyle.BGCOLOR = '#888888';
        varStyle.STROKECOLOR = 'black';
        stylesheet.pushStyle({type: 'deletion'}, null, varStyle);
    }
    {
        var varStyle = new DASStyle();
        varStyle.glyph = 'PLIMSOLL';
        varStyle.BUMP = 'yes';
        varStyle.LABEL = 'no';
        varStyle.FGCOLOR = 'rgb(50,80,255)';
        varStyle.BGCOLOR = '#888888';
        varStyle.STROKECOLOR = 'black';
        stylesheet.pushStyle({type: 'default'}, null, varStyle);
    }

    return callback(stylesheet);
}

VCFParser.prototype.getDefaultFIPs = function(callback) {
    var self = this;
    var fip = function(feature, featureInfo) {
        featureInfo.add("Ref. allele", feature.refAllele);
        featureInfo.add("Alt. alleles", feature.altAlleles.join(','));

        if (feature.info) {
            for (var ii = 0; ii < self.info.length; ++ii) {
                var info = self.info[ii];
                var val = feature.info[info.id];
                if (val !== undefined) {
                    featureInfo.add(info.desc, val);
                }
            }
        }
    };
    callback(fip);
}

dalliance_registerParserFactory('vcf', function() {return new VCFParser()});/* -*- mode: javascript; c-basic-offset: 4; indent-tabs-mode: nil -*- */

// 
// Dalliance Genome Explorer
// (c) Thomas Down 2006-2014
//
// bedwig.js
//

"use strict";

function BedWigParser(type) {
    this.type = type;
}

BedWigParser.prototype.createSession = function(sink) {
    if (this.type == 'wig')
        return new WigParseSession(this, sink);
    else
        return new BedParseSession(this, sink);
}

var __KV_REGEXP=/([^=]+)=(.+)/;
var __SPACE_REGEXP=/\s/;


function BedParseSession(parser, sink) {
    this.parser = parser;
    this.sink = sink;
}

BedParseSession.prototype.parse = function(line) {
    var toks = line.split(__SPACE_REGEXP);
    if (toks.length < 3)
        return;

    var start = parseInt(toks[1]) + 1;
    var end = parseInt(toks[2]);

    var f = {segment: toks[0], 
             min: start,
             max: end};

    if (toks.length > 3 && toks[3] !== '.') {
        f.label = toks[3];
    }

    if (toks.length > 4) {
        f.score = parseFloat(toks[4])
    }

    if (toks.length > 5) {
        f.orientation = toks[5];
    }

    if (toks.length > 8) {
        var color = toks[8];
        if (BED_COLOR_REGEXP.test(color)) {
            f.itemRgb = 'rgb(' + color + ')';
        }
    }

    if (toks.length >= 12) {
        var thickStart = parseInt(toks[6]);
        var thickEnd   = parseInt(toks[7]);
        var blockCount = parseInt(toks[9]);
        var blockSizes = toks[10].split(',').map(function(x) {return parseInt(x)});
        var blockStarts = toks[11].split(',').map(function(x) {return parseInt(x)});

        f.type = 'transcript'
        var grp = new DASGroup();
        grp.id = toks[3];
        grp.type = 'transcript'
        grp.notes = [];
        f.groups = [grp];

        if (toks.length > 12) {
            var geneId = toks[12];
            var geneName = geneId;
            if (toks.length > 13) {
                geneName = toks[13];
            }
            var gg = new DASGroup();
            gg.id = geneId;
            gg.label = geneName;
            gg.type = 'gene';
            f.groups.push(gg);
        }  

        var spans = null;
        for (var b = 0; b < blockCount; ++b) {
            var bmin = blockStarts[b] + start;
            var bmax = bmin + blockSizes[b];
            var span = new Range(bmin, bmax);
            if (spans) {
                spans = union(spans, span);
            } else {
                spans = span;
            }
        }
                    
        var tsList = spans.ranges();
        for (var s = 0; s < tsList.length; ++s) {
            var ts = tsList[s];
            var bf = shallowCopy(f);
            bf.min = ts.min();
            bf.max = ts.max();
            this.sink(bf);
        }

        if (thickEnd > thickStart) {
            var tl = intersection(spans, new Range(thickStart, thickEnd));
            if (tl) {
                f.type = 'translation';
                var tlList = tl.ranges();
                for (var s = 0; s < tlList.length; ++s) {
                    var ts = tlList[s];
                    var bf = shallowCopy(f);
                    bf.min = ts.min();
                    bf.max = ts.max();
                    this.sink(bf);
                }
            }
        }
    } else {
        this.sink(f);
    }
}

BedParseSession.prototype.flush = function() {};

function WigParseSession(parser, sink) {
    this.parser = parser;
    this.sink = sink;
    this.wigState = null;
}

WigParseSession.prototype.parse = function(line) {
    var toks = line.split(__SPACE_REGEXP);

    if (toks[0] == 'fixedStep') {
        this.wigState = 'fixedStep';
        this.chr = this.pos = this.step = null;
        this.span = 1;

        for (var ti = 1; ti < toks.length; ++ti) {
            var m = __KV_REGEXP.exec(toks[ti]);
            if (m) {
                if (m[1] == 'chrom') {
                    this.chr = m[2];
                } else if (m[1] == 'start') {
                    this.pos = parseInt(m[2]);
                } else if (m[1] == 'step') {
                    this.step = parseInt(m[2]);
                } else if (m[1] == 'span') {
                    this.span = parseInt(m[2]);
                }
            }
        }
    } else if (toks[0] == 'variableStep') {
        this.wigState = 'variableStep';
        this.chr = null;
        this.span = 1;

        for (var ti = 1; ti < toks.length; ++ti) {
            var m = __KV_REGEXP.exec(toks[ti]);
            if (m[1] == 'chrom') {
                this.chr = m[2];
            } else if (m[1] == 'span') {
                this.span = parseInt(m[2]);
            }
        }
    } else {
        if (!this.wigState) {
            if (toks.length < 4)
                return;

            var f = {segment: toks[0], 
                     min: parseInt(toks[1]) + 1, 
                     max: parseInt(toks[2]),
                     score: parseFloat(toks[3])};

            this.sink(f);
        } else if (this.wigState == 'fixedStep') {
            if (toks.length != 1)
                return;
            var score = parseFloat(toks[0]);
            var f = {segment: this.chr, min: this.pos, max: this.pos + this.span - 1, score: score};
            this.pos += this.step;
            this.sink(f);
        } else if (this.wigState == 'variableStep') {
            if (toks.length != 2)
                return;
            var pos = parseInt(toks[0]);
            var score = parseFloat(toks[1]);
            var f = {segment: this.chr, min: pos, max: pos + this.span - 1, score: score};
            this.sink(f);
        }
    }
}

WigParseSession.prototype.flush = function() {};

BedWigParser.prototype.getStyleSheet = function(callback) {
    var thisB = this;
    var stylesheet = new DASStylesheet();

    if (this.type == 'wig') {
        var wigStyle = new DASStyle();
        wigStyle.glyph = 'HISTOGRAM';
        wigStyle.BGCOLOR = 'blue';
        wigStyle.HEIGHT=30;
        stylesheet.pushStyle({type: 'default'}, null, wigStyle);
    } else {
        var wigStyle = new DASStyle();
        wigStyle.glyph = 'BOX';
        wigStyle.FGCOLOR = 'black';
        wigStyle.BGCOLOR = 'blue'
        wigStyle.HEIGHT = 8;
        wigStyle.BUMP = true;
        wigStyle.LABEL = true;
        wigStyle.ZINDEX = 20;
        stylesheet.pushStyle({type: 'default'}, null, wigStyle);

        var wigStyle = new DASStyle();
        wigStyle.glyph = 'BOX';
        wigStyle.FGCOLOR = 'black';
        wigStyle.BGCOLOR = 'red'
        wigStyle.HEIGHT = 10;
        wigStyle.BUMP = true;
        wigStyle.ZINDEX = 20;
        stylesheet.pushStyle({type: 'translation'}, null, wigStyle);
                
        var tsStyle = new DASStyle();
        tsStyle.glyph = 'BOX';
        tsStyle.FGCOLOR = 'black';
        tsStyle.BGCOLOR = 'white';
        tsStyle.HEIGHT = 10;
        tsStyle.ZINDEX = 10;
        tsStyle.BUMP = true;
        tsStyle.LABEL = true;
        stylesheet.pushStyle({type: 'transcript'}, null, tsStyle);

        var densStyle = new DASStyle();
        densStyle.glyph = 'HISTOGRAM';
        densStyle.COLOR1 = 'white';
        densStyle.COLOR2 = 'black';
        densStyle.HEIGHT=30;
        stylesheet.pushStyle({type: 'density'}, null, densStyle);
    }

    return callback(stylesheet);
}

dalliance_registerParserFactory('bed', function(t) {return new BedWigParser(t)});
dalliance_registerParserFactory('wig', function(t) {return new BedWigParser(t)});/* -*- mode: javascript; c-basic-offset: 4; indent-tabs-mode: nil -*- */

// 
// Dalliance Genome Explorer
// (c) Thomas Down 2006-2014
//
// bedwig.js
//

"use strict";

function probeResource(source, listener, retry) {
    var BED_REGEXP = new RegExp('^.+\t[0-9]+\t[0-9]+.*$');
    var KV_REGEXP=/([^=]+)=\"?([^\"]+)\"?/;
    var VCFHEAD_RE = /^##\s*fileformat=VCFv4\..+/;

    var fetchable;
    if (source.blob)
        fetchable = new BlobFetchable(source.blob);
    else 
        fetchable = new URLFetchable(source.uri);

    fetchable.slice(0, 1<<16).salted().fetch(function(result, error) {
        if (!result) {
            if (!retry) {
                source.credentials = true;
                probeResource(source, listener, true)
            }

            return listener(source, "Couldn't fetch data");
        }

        var ba = new Uint8Array(result);
        var la = new Uint32Array(result, 0, 1);
        var magic = la[0];
        if (magic == BIG_WIG_MAGIC || magic == BIG_BED_MAGIC) {
            source.tier_type = 'bwg';
            var nameExtractPattern = new RegExp('/?([^/]+?)(.bw|.bb|.bigWig|.bigBed)?$');
            var match = nameExtractPattern.exec(source.uri || source.blob.name);
            if (match) {
                source.name = match[1];
            }

            return listener(source, null);
        } else if (magic == BAI_MAGIC) {
            source.tier_type = 'bai';
            return listener(source, null);
        } else if (ba[0] == 31 || ba[1] == 139) {
            var unc = unbgzf(result);
            var uncba = new Uint8Array(unc);
            magic = readInt(uncba, 0);
            if (magic == BAM_MAGIC) {
                source.tier_type = 'bam';
                var nameExtractPattern = new RegExp('/?([^/]+?)(.bam)?$');
                var match = nameExtractPattern.exec(source.uri || source.blob.name);
                if (match) {
                    source.name = match[1];
                }

                return listener(source, null);
            } else if (magic == TABIX_MAGIC) {
                source.tier_type = 'tabix-index';
                return listener(source, null);
            } else if (magic == 0x69662323) {
                source.tier_type = 'tabix';
                source.payload = 'vcf';
                var nameExtractPattern = new RegExp('/?([^/]+?)(.vcf)?(.gz)?$');
                var match = nameExtractPattern.exec(source.uri || source.blob.name);
                if (match) {
                    source.name = match[1];
                }

                return listener(source, null);
            } else {
                console.log('magic = ' + magic.toString(16));
               return listener(source, "Unsupported format");
            }
        } else {
            var text = String.fromCharCode.apply(null, ba);
            var lines = text.split("\n");

            if (lines.length > 0 && VCFHEAD_RE.test(lines[0])) {
                source.tier_type = 'memstore';
                source.payload = 'vcf';
                var nameExtractPattern = new RegExp('/?([^/]+?)(\.vcf)?$');
                var match = nameExtractPattern.exec(source.uri || source.blob.name);
                if (match && !source.name) {
                    source.name = match[1];
                }
                return listener(source, null);
            }

            for (var li = 0; li < lines.length; ++li) {
                var line = lines[li].replace('\r', '');
                if (line.length == 0) continue;

                if (line.indexOf('browser') == 0) continue;

                if (line.indexOf('track') == 0) {
                    var maybeType = 'bed';
                    var toks = line.split(/\s/);
                    for (var ti = 1; ti < toks.length; ++ti) {
                        var m = KV_REGEXP.exec(toks[ti]);
                        if (m) {
                            if (m[1] == 'type' && m[2] == 'wiggle_0') {
                                maybeType = 'wig'
                            } else if (m[0] == 'name') {
                                source.name = m[2];
                            }
                        }
                    }

                    finishProbeBedWig(source, maybeType);
                    return listener(source, null);
                }

                if (line.indexOf('fixedStep') == 0) {
                    finishProbeBedWig(source, 'wig');
                    return listener(source, null);
                }

                if (line.indexOf('variableStep') == 0) {
                    finishProbeBedWig(source, 'wig');
                    return listener(source, null);
                }

                if (BED_REGEXP.test(line)) {
                    finishProbeBedWig(source, null);
                    return listener(source, null);
                }

                break;
            }

            return listener(source, "Unsupported format");
        }
    });
}

function finishProbeBedWig(source, maybeType) {
    source.tier_type = 'memstore';
    var nameExtractPattern = new RegExp('/?([^/]+?)(.(bed|wig))?$');
    var match = nameExtractPattern.exec(source.uri || source.blob.name);
    if (match) {
        if (!source.name)
            source.name = match[1];
        if (!maybeType && match[3]) {
            maybeType = match[3];
        }
    }
    source.payload = maybeType || 'bed';
}/* -*- mode: javascript; c-basic-offset: 4; indent-tabs-mode: nil -*- */

// 
// Dalliance Genome Explorer
// (c) Thomas Down 2006-2014
//
// export-ui.js
//

Browser.prototype.openExportPanel = function() {
    var b = this;
    if (this.uiMode === 'export') {
        this.hideToolPanel();
        this.setUiMode('none');
    } else {
        var exportForm = makeElement('div', null, {className: 'export-form'});
        var exportSelect = makeElement('select');
        exportSelect.appendChild(makeElement('option', 'SVG', {value: 'svg'}));
        exportSelect.appendChild(makeElement('option', 'Dalliance config', {value: 'config'}));
        exportSelect.appendChild(makeElement('option', 'Dalliance sources', {value: 'sources'}));
        exportSelect.appendChild(makeElement('option', 'Dalliance page', {value: 'page'}));
        exportSelect.value = 'svg';
        exportForm.appendChild(makeElement('p', ['Export as: ', exportSelect]));

        var exportHighlightsToggle = makeElement('input', null, {type: 'checkbox', checked: this.exportHighlights});
        exportHighlightsToggle.addEventListener('change', function(ev) {
            b.exportHighlights = exportHighlightsToggle.checked;
            b.storeStatus();
        }, false);
        var exportRulerToggle = makeElement('input', null, {type: 'checkbox', checked: this.exportRuler});
        exportRulerToggle.addEventListener('change', function(ev) {
            b.exportRuler = exportRulerToggle.checked;
            b.storeStatus();
        }, false);

        var exportButton = makeElement('button', 'Export', {className: 'btn btn-primary'});
        exportButton.addEventListener('click', function(ev) {
            removeChildren(exportContent);

            var blobURL;
            var note, type, name;
            if (exportSelect.value === 'svg') {
                blobURL = URL.createObjectURL(b.makeSVG({highlights: exportHighlightsToggle.checked,
                                                         ruler: exportRulerToggle.checked ? b.rulerLocation : 'none'}));
                note = 'SVG';
                type = 'image/svg';
                name = 'dalliance-view.svg';
            } else if (exportSelect.value === 'config') {
                var config = JSON.stringify(b.exportFullConfig(), null, 2);
                var blob = new Blob([config], {type: 'text/plain'});
                blobURL = URL.createObjectURL(blob);
                note = 'Configuration';
                type = 'text/plain';
                name = 'dalliance-config.json';
            } else if (exportSelect.value === 'sources') {
                var config = JSON.stringify(b.exportSourceConfig(), null, 2);
                var blob = new Blob([config], {type: 'text/plain'});
                blobURL = URL.createObjectURL(blob);
                note = 'Source array';
                type = 'text/plain';
                name = 'dalliance-sources.json';
            } else if (exportSelect.value === 'page') {
                var page = b.exportPageTemplate();
                var type = 'text/html';
                var blob = new Blob([page], {type: type});
                blobURL = URL.createObjectURL(blob);
                note = 'Page template';
                name = 'dalliance-view.html';
            }

            if (blobURL) {
                var downloadLink = makeElement('a', '[Download]', {
                    href: blobURL,
                    download: 'dalliance-view.svg',
                    type: type
                });

                var previewLink = makeElement('a', '[Preview in browser]', {
                    href: blobURL,
                    type: type,
                    target: '_new'
                });

                exportContent.appendChild(makeElement('p', ['' + note + ' created: ', downloadLink, previewLink]));
            }
        }, false);

        b.addViewListener(function() {
            removeChildren(exportContent);
        });
        b.addTierListener(function() {
            removeChildren(exportContent);
        });

        var exportContent = makeElement('p', '');

        var exportOptsTable = makeElement('table',
            [makeElement('tr',
                [makeElement('th', 'Include highlights', {}, {width: '200px', textAlign: 'right'}),
                 makeElement('td', exportHighlightsToggle)]),
             makeElement('tr',
                [makeElement('th', 'Include vertical guideline'),
                 makeElement('td', exportRulerToggle)])
            ]);

        exportForm.appendChild(exportOptsTable);
        exportForm.appendChild(exportButton);
        exportForm.appendChild(exportContent);

        if (this.uiMode !== 'none')
            this.hideToolPanel();
        this.browserHolder.insertBefore(exportForm, this.svgHolder);
        this.activeToolPanel = exportForm;

        this.setUiMode('export');
    }
}
/* -*- mode: javascript; c-basic-offset: 4; indent-tabs-mode: nil -*- */

// 
// Dalliance Genome Explorer
// (c) Thomas Down 2006-2014
//
// export-config.js
//

Browser.prototype.exportFullConfig = function(opts) {
    opts = opts || {};

    var config = {
        chr: this.chr,
        viewStart: this.viewStart|0,
        viewEnd: this.viewEnd|0,
        cookieKey: 'dalliance_' + hex_sha1(Date.now()),

        coordSystem: this.coordSystem,

        sources: this.exportSourceConfig()
    };

    if (this.uiPrefix)
        config.uiPrefix = this.uiPrefix;

    return config;
}

Browser.prototype.exportSourceConfig = function(opts) {
    opts = opts || {};

    var sourceConfig = [];
    for (var ti = 0; ti < this.tiers.length; ++ti) {
        var tier = this.tiers[ti];
        var source = shallowCopy(tier.dasSource);
        source.coords = undefined;
        source.props = undefined;
        if (!source.disabled)
            source.disabled = undefined;

        if (tier.config.stylesheet) {
            source.style = tier.config.stylesheet.styles;
            source.stylesheet_uri = undefined;
        } else if (source.style) {
            source.style = copyStylesheet({styles: source.style}).styles;
        }

        if (typeof(tier.config.name) === 'string') {
            source.name = tier.config.name;
        }

        if (tier.config.height !== undefined) {
            source.forceHeight = tier.config.height;
        }
        if (tier.config.forceMin !== undefined) {
            source.forceMin = tier.config.forceMin;
        }
        if (tier.config.forceMinDynamic)
            source.forceMinDynamic = tier.config.forceMinDynamic;
        if (tier.config.forceMax !== undefined) {
            source.forceMax = tier.config.forceMax;
        }
        if (tier.config.forceMaxDynamic)
            source.forceMaxDynamic = tier.config.forceMaxDynamic;

        sourceConfig.push(source);
    }

    return sourceConfig;
}

Browser.prototype.exportPageTemplate = function(opts) {
    opts = opts || {};
    var template = '<html>\n' +
                   '  <head>\n' +
                   '    <script language="javascript" src="' + this.uiPrefix + 'dalliance-compiled.js"></script>\n' +
                   '    <script language="javascript">\n' +
                   '      var dalliance_browser = new Browser(' + JSON.stringify(this.exportFullConfig(opts), null, 2) + ');\n' +
                   '    </script>\n' +  
                   '  </head>\n' +
                   '  <body>\n' +
                   '    <div id="svgHolder">Dalliance goes here</div>\n' +
                   '  </body>\n' +
                   '<html>\n';

    return template;
}/* -*- mode: javascript; c-basic-offset: 4; indent-tabs-mode: nil -*- */

// 
// Javascript ZLib
// By Thomas Down 2010-2011
//
// Based very heavily on portions of jzlib (by ymnk@jcraft.com), who in
// turn credits Jean-loup Gailly and Mark Adler for the original zlib code.
//
// inflate.js: ZLib inflate code
//

//
// Shared constants
//

var MAX_WBITS=15; // 32K LZ77 window
var DEF_WBITS=MAX_WBITS;
var MAX_MEM_LEVEL=9;
var MANY=1440;
var BMAX = 15;

// preset dictionary flag in zlib header
var PRESET_DICT=0x20;

var Z_NO_FLUSH=0;
var Z_PARTIAL_FLUSH=1;
var Z_SYNC_FLUSH=2;
var Z_FULL_FLUSH=3;
var Z_FINISH=4;

var Z_DEFLATED=8;

var Z_OK=0;
var Z_STREAM_END=1;
var Z_NEED_DICT=2;
var Z_ERRNO=-1;
var Z_STREAM_ERROR=-2;
var Z_DATA_ERROR=-3;
var Z_MEM_ERROR=-4;
var Z_BUF_ERROR=-5;
var Z_VERSION_ERROR=-6;

var METHOD=0;   // waiting for method byte
var FLAG=1;     // waiting for flag byte
var DICT4=2;    // four dictionary check bytes to go
var DICT3=3;    // three dictionary check bytes to go
var DICT2=4;    // two dictionary check bytes to go
var DICT1=5;    // one dictionary check byte to go
var DICT0=6;    // waiting for inflateSetDictionary
var BLOCKS=7;   // decompressing blocks
var CHECK4=8;   // four check bytes to go
var CHECK3=9;   // three check bytes to go
var CHECK2=10;  // two check bytes to go
var CHECK1=11;  // one check byte to go
var DONE=12;    // finished check, done
var BAD=13;     // got an error--stay here

var inflate_mask = [0x00000000, 0x00000001, 0x00000003, 0x00000007, 0x0000000f, 0x0000001f, 0x0000003f, 0x0000007f, 0x000000ff, 0x000001ff, 0x000003ff, 0x000007ff, 0x00000fff, 0x00001fff, 0x00003fff, 0x00007fff, 0x0000ffff];

var IB_TYPE=0;  // get type bits (3, including end bit)
var IB_LENS=1;  // get lengths for stored
var IB_STORED=2;// processing stored block
var IB_TABLE=3; // get table lengths
var IB_BTREE=4; // get bit lengths tree for a dynamic block
var IB_DTREE=5; // get length, distance trees for a dynamic block
var IB_CODES=6; // processing fixed or dynamic block
var IB_DRY=7;   // output remaining window bytes
var IB_DONE=8;  // finished last block, done
var IB_BAD=9;   // ot a data error--stuck here

var fixed_bl = 9;
var fixed_bd = 5;

var fixed_tl = [
    96,7,256, 0,8,80, 0,8,16, 84,8,115,
    82,7,31, 0,8,112, 0,8,48, 0,9,192,
    80,7,10, 0,8,96, 0,8,32, 0,9,160,
    0,8,0, 0,8,128, 0,8,64, 0,9,224,
    80,7,6, 0,8,88, 0,8,24, 0,9,144,
    83,7,59, 0,8,120, 0,8,56, 0,9,208,
    81,7,17, 0,8,104, 0,8,40, 0,9,176,
    0,8,8, 0,8,136, 0,8,72, 0,9,240,
    80,7,4, 0,8,84, 0,8,20, 85,8,227,
    83,7,43, 0,8,116, 0,8,52, 0,9,200,
    81,7,13, 0,8,100, 0,8,36, 0,9,168,
    0,8,4, 0,8,132, 0,8,68, 0,9,232,
    80,7,8, 0,8,92, 0,8,28, 0,9,152,
    84,7,83, 0,8,124, 0,8,60, 0,9,216,
    82,7,23, 0,8,108, 0,8,44, 0,9,184,
    0,8,12, 0,8,140, 0,8,76, 0,9,248,
    80,7,3, 0,8,82, 0,8,18, 85,8,163,
    83,7,35, 0,8,114, 0,8,50, 0,9,196,
    81,7,11, 0,8,98, 0,8,34, 0,9,164,
    0,8,2, 0,8,130, 0,8,66, 0,9,228,
    80,7,7, 0,8,90, 0,8,26, 0,9,148,
    84,7,67, 0,8,122, 0,8,58, 0,9,212,
    82,7,19, 0,8,106, 0,8,42, 0,9,180,
    0,8,10, 0,8,138, 0,8,74, 0,9,244,
    80,7,5, 0,8,86, 0,8,22, 192,8,0,
    83,7,51, 0,8,118, 0,8,54, 0,9,204,
    81,7,15, 0,8,102, 0,8,38, 0,9,172,
    0,8,6, 0,8,134, 0,8,70, 0,9,236,
    80,7,9, 0,8,94, 0,8,30, 0,9,156,
    84,7,99, 0,8,126, 0,8,62, 0,9,220,
    82,7,27, 0,8,110, 0,8,46, 0,9,188,
    0,8,14, 0,8,142, 0,8,78, 0,9,252,
    96,7,256, 0,8,81, 0,8,17, 85,8,131,
    82,7,31, 0,8,113, 0,8,49, 0,9,194,
    80,7,10, 0,8,97, 0,8,33, 0,9,162,
    0,8,1, 0,8,129, 0,8,65, 0,9,226,
    80,7,6, 0,8,89, 0,8,25, 0,9,146,
    83,7,59, 0,8,121, 0,8,57, 0,9,210,
    81,7,17, 0,8,105, 0,8,41, 0,9,178,
    0,8,9, 0,8,137, 0,8,73, 0,9,242,
    80,7,4, 0,8,85, 0,8,21, 80,8,258,
    83,7,43, 0,8,117, 0,8,53, 0,9,202,
    81,7,13, 0,8,101, 0,8,37, 0,9,170,
    0,8,5, 0,8,133, 0,8,69, 0,9,234,
    80,7,8, 0,8,93, 0,8,29, 0,9,154,
    84,7,83, 0,8,125, 0,8,61, 0,9,218,
    82,7,23, 0,8,109, 0,8,45, 0,9,186,
    0,8,13, 0,8,141, 0,8,77, 0,9,250,
    80,7,3, 0,8,83, 0,8,19, 85,8,195,
    83,7,35, 0,8,115, 0,8,51, 0,9,198,
    81,7,11, 0,8,99, 0,8,35, 0,9,166,
    0,8,3, 0,8,131, 0,8,67, 0,9,230,
    80,7,7, 0,8,91, 0,8,27, 0,9,150,
    84,7,67, 0,8,123, 0,8,59, 0,9,214,
    82,7,19, 0,8,107, 0,8,43, 0,9,182,
    0,8,11, 0,8,139, 0,8,75, 0,9,246,
    80,7,5, 0,8,87, 0,8,23, 192,8,0,
    83,7,51, 0,8,119, 0,8,55, 0,9,206,
    81,7,15, 0,8,103, 0,8,39, 0,9,174,
    0,8,7, 0,8,135, 0,8,71, 0,9,238,
    80,7,9, 0,8,95, 0,8,31, 0,9,158,
    84,7,99, 0,8,127, 0,8,63, 0,9,222,
    82,7,27, 0,8,111, 0,8,47, 0,9,190,
    0,8,15, 0,8,143, 0,8,79, 0,9,254,
    96,7,256, 0,8,80, 0,8,16, 84,8,115,
    82,7,31, 0,8,112, 0,8,48, 0,9,193,

    80,7,10, 0,8,96, 0,8,32, 0,9,161,
    0,8,0, 0,8,128, 0,8,64, 0,9,225,
    80,7,6, 0,8,88, 0,8,24, 0,9,145,
    83,7,59, 0,8,120, 0,8,56, 0,9,209,
    81,7,17, 0,8,104, 0,8,40, 0,9,177,
    0,8,8, 0,8,136, 0,8,72, 0,9,241,
    80,7,4, 0,8,84, 0,8,20, 85,8,227,
    83,7,43, 0,8,116, 0,8,52, 0,9,201,
    81,7,13, 0,8,100, 0,8,36, 0,9,169,
    0,8,4, 0,8,132, 0,8,68, 0,9,233,
    80,7,8, 0,8,92, 0,8,28, 0,9,153,
    84,7,83, 0,8,124, 0,8,60, 0,9,217,
    82,7,23, 0,8,108, 0,8,44, 0,9,185,
    0,8,12, 0,8,140, 0,8,76, 0,9,249,
    80,7,3, 0,8,82, 0,8,18, 85,8,163,
    83,7,35, 0,8,114, 0,8,50, 0,9,197,
    81,7,11, 0,8,98, 0,8,34, 0,9,165,
    0,8,2, 0,8,130, 0,8,66, 0,9,229,
    80,7,7, 0,8,90, 0,8,26, 0,9,149,
    84,7,67, 0,8,122, 0,8,58, 0,9,213,
    82,7,19, 0,8,106, 0,8,42, 0,9,181,
    0,8,10, 0,8,138, 0,8,74, 0,9,245,
    80,7,5, 0,8,86, 0,8,22, 192,8,0,
    83,7,51, 0,8,118, 0,8,54, 0,9,205,
    81,7,15, 0,8,102, 0,8,38, 0,9,173,
    0,8,6, 0,8,134, 0,8,70, 0,9,237,
    80,7,9, 0,8,94, 0,8,30, 0,9,157,
    84,7,99, 0,8,126, 0,8,62, 0,9,221,
    82,7,27, 0,8,110, 0,8,46, 0,9,189,
    0,8,14, 0,8,142, 0,8,78, 0,9,253,
    96,7,256, 0,8,81, 0,8,17, 85,8,131,
    82,7,31, 0,8,113, 0,8,49, 0,9,195,
    80,7,10, 0,8,97, 0,8,33, 0,9,163,
    0,8,1, 0,8,129, 0,8,65, 0,9,227,
    80,7,6, 0,8,89, 0,8,25, 0,9,147,
    83,7,59, 0,8,121, 0,8,57, 0,9,211,
    81,7,17, 0,8,105, 0,8,41, 0,9,179,
    0,8,9, 0,8,137, 0,8,73, 0,9,243,
    80,7,4, 0,8,85, 0,8,21, 80,8,258,
    83,7,43, 0,8,117, 0,8,53, 0,9,203,
    81,7,13, 0,8,101, 0,8,37, 0,9,171,
    0,8,5, 0,8,133, 0,8,69, 0,9,235,
    80,7,8, 0,8,93, 0,8,29, 0,9,155,
    84,7,83, 0,8,125, 0,8,61, 0,9,219,
    82,7,23, 0,8,109, 0,8,45, 0,9,187,
    0,8,13, 0,8,141, 0,8,77, 0,9,251,
    80,7,3, 0,8,83, 0,8,19, 85,8,195,
    83,7,35, 0,8,115, 0,8,51, 0,9,199,
    81,7,11, 0,8,99, 0,8,35, 0,9,167,
    0,8,3, 0,8,131, 0,8,67, 0,9,231,
    80,7,7, 0,8,91, 0,8,27, 0,9,151,
    84,7,67, 0,8,123, 0,8,59, 0,9,215,
    82,7,19, 0,8,107, 0,8,43, 0,9,183,
    0,8,11, 0,8,139, 0,8,75, 0,9,247,
    80,7,5, 0,8,87, 0,8,23, 192,8,0,
    83,7,51, 0,8,119, 0,8,55, 0,9,207,
    81,7,15, 0,8,103, 0,8,39, 0,9,175,
    0,8,7, 0,8,135, 0,8,71, 0,9,239,
    80,7,9, 0,8,95, 0,8,31, 0,9,159,
    84,7,99, 0,8,127, 0,8,63, 0,9,223,
    82,7,27, 0,8,111, 0,8,47, 0,9,191,
    0,8,15, 0,8,143, 0,8,79, 0,9,255
];
var fixed_td = [
    80,5,1, 87,5,257, 83,5,17, 91,5,4097,
    81,5,5, 89,5,1025, 85,5,65, 93,5,16385,
    80,5,3, 88,5,513, 84,5,33, 92,5,8193,
    82,5,9, 90,5,2049, 86,5,129, 192,5,24577,
    80,5,2, 87,5,385, 83,5,25, 91,5,6145,
    81,5,7, 89,5,1537, 85,5,97, 93,5,24577,
    80,5,4, 88,5,769, 84,5,49, 92,5,12289,
    82,5,13, 90,5,3073, 86,5,193, 192,5,24577
];

  // Tables for deflate from PKZIP's appnote.txt.
  var cplens = [ // Copy lengths for literal codes 257..285
        3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 15, 17, 19, 23, 27, 31,
        35, 43, 51, 59, 67, 83, 99, 115, 131, 163, 195, 227, 258, 0, 0
  ];

  // see note #13 above about 258
  var cplext = [ // Extra bits for literal codes 257..285
        0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2,
        3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 0, 112, 112  // 112==invalid
  ];

 var cpdist = [ // Copy offsets for distance codes 0..29
        1, 2, 3, 4, 5, 7, 9, 13, 17, 25, 33, 49, 65, 97, 129, 193,
        257, 385, 513, 769, 1025, 1537, 2049, 3073, 4097, 6145,
        8193, 12289, 16385, 24577
  ];

  var cpdext = [ // Extra bits for distance codes
        0, 0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6,
        7, 7, 8, 8, 9, 9, 10, 10, 11, 11,
        12, 12, 13, 13];

//
// ZStream.java
//

function ZStream() {
}


ZStream.prototype.inflateInit = function(w, nowrap) {
    if (!w) {
	w = DEF_WBITS;
    }
    if (nowrap) {
	nowrap = false;
    }
    this.istate = new Inflate();
    return this.istate.inflateInit(this, nowrap?-w:w);
}

ZStream.prototype.inflate = function(f) {
    if(this.istate==null) return Z_STREAM_ERROR;
    return this.istate.inflate(this, f);
}

ZStream.prototype.inflateEnd = function(){
    if(this.istate==null) return Z_STREAM_ERROR;
    var ret=istate.inflateEnd(this);
    this.istate = null;
    return ret;
}
ZStream.prototype.inflateSync = function(){
    // if(istate == null) return Z_STREAM_ERROR;
    return istate.inflateSync(this);
}
ZStream.prototype.inflateSetDictionary = function(dictionary, dictLength){
    // if(istate == null) return Z_STREAM_ERROR;
    return istate.inflateSetDictionary(this, dictionary, dictLength);
}

/*

  public int deflateInit(int level){
    return deflateInit(level, MAX_WBITS);
  }
  public int deflateInit(int level, boolean nowrap){
    return deflateInit(level, MAX_WBITS, nowrap);
  }
  public int deflateInit(int level, int bits){
    return deflateInit(level, bits, false);
  }
  public int deflateInit(int level, int bits, boolean nowrap){
    dstate=new Deflate();
    return dstate.deflateInit(this, level, nowrap?-bits:bits);
  }
  public int deflate(int flush){
    if(dstate==null){
      return Z_STREAM_ERROR;
    }
    return dstate.deflate(this, flush);
  }
  public int deflateEnd(){
    if(dstate==null) return Z_STREAM_ERROR;
    int ret=dstate.deflateEnd();
    dstate=null;
    return ret;
  }
  public int deflateParams(int level, int strategy){
    if(dstate==null) return Z_STREAM_ERROR;
    return dstate.deflateParams(this, level, strategy);
  }
  public int deflateSetDictionary (byte[] dictionary, int dictLength){
    if(dstate == null)
      return Z_STREAM_ERROR;
    return dstate.deflateSetDictionary(this, dictionary, dictLength);
  }

*/

/*
  // Flush as much pending output as possible. All deflate() output goes
  // through this function so some applications may wish to modify it
  // to avoid allocating a large strm->next_out buffer and copying into it.
  // (See also read_buf()).
  void flush_pending(){
    int len=dstate.pending;

    if(len>avail_out) len=avail_out;
    if(len==0) return;

    if(dstate.pending_buf.length<=dstate.pending_out ||
       next_out.length<=next_out_index ||
       dstate.pending_buf.length<(dstate.pending_out+len) ||
       next_out.length<(next_out_index+len)){
      System.out.println(dstate.pending_buf.length+", "+dstate.pending_out+
			 ", "+next_out.length+", "+next_out_index+", "+len);
      System.out.println("avail_out="+avail_out);
    }

    System.arraycopy(dstate.pending_buf, dstate.pending_out,
		     next_out, next_out_index, len);

    next_out_index+=len;
    dstate.pending_out+=len;
    total_out+=len;
    avail_out-=len;
    dstate.pending-=len;
    if(dstate.pending==0){
      dstate.pending_out=0;
    }
  }

  // Read a new buffer from the current input stream, update the adler32
  // and total number of bytes read.  All deflate() input goes through
  // this function so some applications may wish to modify it to avoid
  // allocating a large strm->next_in buffer and copying from it.
  // (See also flush_pending()).
  int read_buf(byte[] buf, int start, int size) {
    int len=avail_in;

    if(len>size) len=size;
    if(len==0) return 0;

    avail_in-=len;

    if(dstate.noheader==0) {
      adler=_adler.adler32(adler, next_in, next_in_index, len);
    }
    System.arraycopy(next_in, next_in_index, buf, start, len);
    next_in_index  += len;
    total_in += len;
    return len;
  }

  public void free(){
    next_in=null;
    next_out=null;
    msg=null;
    _adler=null;
  }
}
*/


//
// Inflate.java
//

function Inflate() {
    this.was = [0];
}

Inflate.prototype.inflateReset = function(z) {
    if(z == null || z.istate == null) return Z_STREAM_ERROR;
    
    z.total_in = z.total_out = 0;
    z.msg = null;
    z.istate.mode = z.istate.nowrap!=0 ? BLOCKS : METHOD;
    z.istate.blocks.reset(z, null);
    return Z_OK;
}

Inflate.prototype.inflateEnd = function(z){
    if(this.blocks != null)
      this.blocks.free(z);
    this.blocks=null;
    return Z_OK;
}

Inflate.prototype.inflateInit = function(z, w){
    z.msg = null;
    this.blocks = null;

    // handle undocumented nowrap option (no zlib header or check)
    var nowrap = 0;
    if(w < 0){
      w = - w;
      nowrap = 1;
    }

    // set window size
    if(w<8 ||w>15){
      this.inflateEnd(z);
      return Z_STREAM_ERROR;
    }
    this.wbits=w;

    z.istate.blocks=new InfBlocks(z, 
				  z.istate.nowrap!=0 ? null : this,
				  1<<w);

    // reset state
    this.inflateReset(z);
    return Z_OK;
  }

Inflate.prototype.inflate = function(z, f){
    var r, b;

    if(z == null || z.istate == null || z.next_in == null)
      return Z_STREAM_ERROR;
    f = f == Z_FINISH ? Z_BUF_ERROR : Z_OK;
    r = Z_BUF_ERROR;
    while (true){
      switch (z.istate.mode){
      case METHOD:

        if(z.avail_in==0)return r;r=f;

        z.avail_in--; z.total_in++;
        if(((z.istate.method = z.next_in[z.next_in_index++])&0xf)!=Z_DEFLATED){
          z.istate.mode = BAD;
          z.msg="unknown compression method";
          z.istate.marker = 5;       // can't try inflateSync
          break;
        }
        if((z.istate.method>>4)+8>z.istate.wbits){
          z.istate.mode = BAD;
          z.msg="invalid window size";
          z.istate.marker = 5;       // can't try inflateSync
          break;
        }
        z.istate.mode=FLAG;
      case FLAG:

        if(z.avail_in==0)return r;r=f;

        z.avail_in--; z.total_in++;
        b = (z.next_in[z.next_in_index++])&0xff;

        if((((z.istate.method << 8)+b) % 31)!=0){
          z.istate.mode = BAD;
          z.msg = "incorrect header check";
          z.istate.marker = 5;       // can't try inflateSync
          break;
        }

        if((b&PRESET_DICT)==0){
          z.istate.mode = BLOCKS;
          break;
        }
        z.istate.mode = DICT4;
      case DICT4:

        if(z.avail_in==0)return r;r=f;

        z.avail_in--; z.total_in++;
        z.istate.need=((z.next_in[z.next_in_index++]&0xff)<<24)&0xff000000;
        z.istate.mode=DICT3;
      case DICT3:

        if(z.avail_in==0)return r;r=f;

        z.avail_in--; z.total_in++;
        z.istate.need+=((z.next_in[z.next_in_index++]&0xff)<<16)&0xff0000;
        z.istate.mode=DICT2;
      case DICT2:

        if(z.avail_in==0)return r;r=f;

        z.avail_in--; z.total_in++;
        z.istate.need+=((z.next_in[z.next_in_index++]&0xff)<<8)&0xff00;
        z.istate.mode=DICT1;
      case DICT1:

        if(z.avail_in==0)return r;r=f;

        z.avail_in--; z.total_in++;
        z.istate.need += (z.next_in[z.next_in_index++]&0xff);
        z.adler = z.istate.need;
        z.istate.mode = DICT0;
        return Z_NEED_DICT;
      case DICT0:
        z.istate.mode = BAD;
        z.msg = "need dictionary";
        z.istate.marker = 0;       // can try inflateSync
        return Z_STREAM_ERROR;
      case BLOCKS:

        r = z.istate.blocks.proc(z, r);
        if(r == Z_DATA_ERROR){
          z.istate.mode = BAD;
          z.istate.marker = 0;     // can try inflateSync
          break;
        }
        if(r == Z_OK){
          r = f;
        }
        if(r != Z_STREAM_END){
          return r;
        }
        r = f;
        z.istate.blocks.reset(z, z.istate.was);
        if(z.istate.nowrap!=0){
          z.istate.mode=DONE;
          break;
        }
        z.istate.mode=CHECK4;
      case CHECK4:

        if(z.avail_in==0)return r;r=f;

        z.avail_in--; z.total_in++;
        z.istate.need=((z.next_in[z.next_in_index++]&0xff)<<24)&0xff000000;
        z.istate.mode=CHECK3;
      case CHECK3:

        if(z.avail_in==0)return r;r=f;

        z.avail_in--; z.total_in++;
        z.istate.need+=((z.next_in[z.next_in_index++]&0xff)<<16)&0xff0000;
        z.istate.mode = CHECK2;
      case CHECK2:

        if(z.avail_in==0)return r;r=f;

        z.avail_in--; z.total_in++;
        z.istate.need+=((z.next_in[z.next_in_index++]&0xff)<<8)&0xff00;
        z.istate.mode = CHECK1;
      case CHECK1:

        if(z.avail_in==0)return r;r=f;

        z.avail_in--; z.total_in++;
        z.istate.need+=(z.next_in[z.next_in_index++]&0xff);

        if(((z.istate.was[0])) != ((z.istate.need))){
          z.istate.mode = BAD;
          z.msg = "incorrect data check";
          z.istate.marker = 5;       // can't try inflateSync
          break;
        }

        z.istate.mode = DONE;
      case DONE:
        return Z_STREAM_END;
      case BAD:
        return Z_DATA_ERROR;
      default:
        return Z_STREAM_ERROR;
      }
    }
  }


Inflate.prototype.inflateSetDictionary = function(z,  dictionary, dictLength) {
    var index=0;
    var length = dictLength;
    if(z==null || z.istate == null|| z.istate.mode != DICT0)
      return Z_STREAM_ERROR;

    if(z._adler.adler32(1, dictionary, 0, dictLength)!=z.adler){
      return Z_DATA_ERROR;
    }

    z.adler = z._adler.adler32(0, null, 0, 0);

    if(length >= (1<<z.istate.wbits)){
      length = (1<<z.istate.wbits)-1;
      index=dictLength - length;
    }
    z.istate.blocks.set_dictionary(dictionary, index, length);
    z.istate.mode = BLOCKS;
    return Z_OK;
  }

//  static private byte[] mark = {(byte)0, (byte)0, (byte)0xff, (byte)0xff};
var mark = [0, 0, 255, 255]

Inflate.prototype.inflateSync = function(z){
    var n;       // number of bytes to look at
    var p;       // pointer to bytes
    var m;       // number of marker bytes found in a row
    var r, w;   // temporaries to save total_in and total_out

    // set up
    if(z == null || z.istate == null)
      return Z_STREAM_ERROR;
    if(z.istate.mode != BAD){
      z.istate.mode = BAD;
      z.istate.marker = 0;
    }
    if((n=z.avail_in)==0)
      return Z_BUF_ERROR;
    p=z.next_in_index;
    m=z.istate.marker;

    // search
    while (n!=0 && m < 4){
      if(z.next_in[p] == mark[m]){
        m++;
      }
      else if(z.next_in[p]!=0){
        m = 0;
      }
      else{
        m = 4 - m;
      }
      p++; n--;
    }

    // restore
    z.total_in += p-z.next_in_index;
    z.next_in_index = p;
    z.avail_in = n;
    z.istate.marker = m;

    // return no joy or set up to restart on a new block
    if(m != 4){
      return Z_DATA_ERROR;
    }
    r=z.total_in;  w=z.total_out;
    this.inflateReset(z);
    z.total_in=r;  z.total_out = w;
    z.istate.mode = BLOCKS;
    return Z_OK;
}

  // Returns true if inflate is currently at the end of a block generated
  // by Z_SYNC_FLUSH or Z_FULL_FLUSH. This function is used by one PPP
  // implementation to provide an additional safety check. PPP uses Z_SYNC_FLUSH
  // but removes the length bytes of the resulting empty stored block. When
  // decompressing, PPP checks that at the end of input packet, inflate is
  // waiting for these length bytes.
Inflate.prototype.inflateSyncPoint = function(z){
    if(z == null || z.istate == null || z.istate.blocks == null)
      return Z_STREAM_ERROR;
    return z.istate.blocks.sync_point();
}


//
// InfBlocks.java
//

var INFBLOCKS_BORDER = [16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15];

function InfBlocks(z, checkfn, w) {
    this.hufts=new Int32Array(MANY*3);
    this.window=new Uint8Array(w);
    this.end=w;
    this.checkfn = checkfn;
    this.mode = IB_TYPE;
    this.reset(z, null);

    this.left = 0;            // if STORED, bytes left to copy 

    this.table = 0;           // table lengths (14 bits) 
    this.index = 0;           // index into blens (or border) 
    this.blens = null;         // bit lengths of codes 
    this.bb=new Int32Array(1); // bit length tree depth 
    this.tb=new Int32Array(1); // bit length decoding tree 

    this.codes = new InfCodes();

    this.last = 0;            // true if this block is the last block 

  // mode independent information 
    this.bitk = 0;            // bits in bit buffer 
    this.bitb = 0;            // bit buffer 
    this.read = 0;            // window read pointer 
    this.write = 0;           // window write pointer 
    this.check = 0;          // check on output 

    this.inftree=new InfTree();
}




InfBlocks.prototype.reset = function(z, c){
    if(c) c[0]=this.check;
    if(this.mode==IB_CODES){
      this.codes.free(z);
    }
    this.mode=IB_TYPE;
    this.bitk=0;
    this.bitb=0;
    this.read=this.write=0;

    if(this.checkfn)
      z.adler=this.check=z._adler.adler32(0, null, 0, 0);
  }

 InfBlocks.prototype.proc = function(z, r){
    var t;              // temporary storage
    var b;              // bit buffer
    var k;              // bits in bit buffer
    var p;              // input data pointer
    var n;              // bytes available there
    var q;              // output window write pointer
    var m;              // bytes to end of window or read pointer

    // copy input/output information to locals (UPDATE macro restores)
    {p=z.next_in_index;n=z.avail_in;b=this.bitb;k=this.bitk;}
    {q=this.write;m=(q<this.read ? this.read-q-1 : this.end-q);}

    // process input based on current state
    while(true){
      switch (this.mode){
      case IB_TYPE:

	while(k<(3)){
	  if(n!=0){
	    r=Z_OK;
	  }
	  else{
	    this.bitb=b; this.bitk=k; 
	    z.avail_in=n;
	    z.total_in+=p-z.next_in_index;z.next_in_index=p;
	    this.write=q;
	    return this.inflate_flush(z,r);
	  };
	  n--;
	  b|=(z.next_in[p++]&0xff)<<k;
	  k+=8;
	}
	t = (b & 7);
	this.last = t & 1;

	switch (t >>> 1){
        case 0:                         // stored 
          {b>>>=(3);k-=(3);}
          t = k & 7;                    // go to byte boundary

          {b>>>=(t);k-=(t);}
          this.mode = IB_LENS;                  // get length of stored block
          break;
        case 1:                         // fixed
          {
              var bl=new Int32Array(1);
	      var bd=new Int32Array(1);
              var tl=[];
	      var td=[];

	      inflate_trees_fixed(bl, bd, tl, td, z);
              this.codes.init(bl[0], bd[0], tl[0], 0, td[0], 0, z);
          }

          {b>>>=(3);k-=(3);}

          this.mode = IB_CODES;
          break;
        case 2:                         // dynamic

          {b>>>=(3);k-=(3);}

          this.mode = IB_TABLE;
          break;
        case 3:                         // illegal

          {b>>>=(3);k-=(3);}
          this.mode = BAD;
          z.msg = "invalid block type";
          r = Z_DATA_ERROR;

	  this.bitb=b; this.bitk=k; 
	  z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
	  this.write=q;
	  return this.inflate_flush(z,r);
	}
	break;
      case IB_LENS:
	while(k<(32)){
	  if(n!=0){
	    r=Z_OK;
	  }
	  else{
	    this.bitb=b; this.bitk=k; 
	    z.avail_in=n;
	    z.total_in+=p-z.next_in_index;z.next_in_index=p;
	    this.write=q;
	    return this.inflate_flush(z,r);
	  };
	  n--;
	  b|=(z.next_in[p++]&0xff)<<k;
	  k+=8;
	}

	if ((((~b) >>> 16) & 0xffff) != (b & 0xffff)){
	  this.mode = BAD;
	  z.msg = "invalid stored block lengths";
	  r = Z_DATA_ERROR;

	  this.bitb=b; this.bitk=k; 
	  z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
	  this.write=q;
	  return this.inflate_flush(z,r);
	}
	this.left = (b & 0xffff);
	b = k = 0;                       // dump bits
	this.mode = this.left!=0 ? IB_STORED : (this.last!=0 ? IB_DRY : IB_TYPE);
	break;
      case IB_STORED:
	if (n == 0){
	  this.bitb=b; this.bitk=k; 
	  z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
	  write=q;
	  return this.inflate_flush(z,r);
	}

	if(m==0){
	  if(q==end&&read!=0){
	    q=0; m=(q<this.read ? this.read-q-1 : this.end-q);
	  }
	  if(m==0){
	    this.write=q; 
	    r=this.inflate_flush(z,r);
	    q=this.write; m = (q < this.read ? this.read-q-1 : this.end-q);
	    if(q==this.end && this.read != 0){
	      q=0; m = (q < this.read ? this.read-q-1 : this.end-q);
	    }
	    if(m==0){
	      this.bitb=b; this.bitk=k; 
	      z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
	      this.write=q;
	      return this.inflate_flush(z,r);
	    }
	  }
	}
	r=Z_OK;

	t = this.left;
	if(t>n) t = n;
	if(t>m) t = m;
	arrayCopy(z.next_in, p, window, q, t);
	p += t;  n -= t;
	q += t;  m -= t;
	if ((this.left -= t) != 0)
	  break;
	this.mode = (this.last != 0 ? IB_DRY : IB_TYPE);
	break;
      case IB_TABLE:

	while(k<(14)){
	  if(n!=0){
	    r=Z_OK;
	  }
	  else{
	    this.bitb=b; this.bitk=k; 
	    z.avail_in=n;
	    z.total_in+=p-z.next_in_index;z.next_in_index=p;
	    this.write=q;
	    return this.inflate_flush(z,r);
	  };
	  n--;
	  b|=(z.next_in[p++]&0xff)<<k;
	  k+=8;
	}

	this.table = t = (b & 0x3fff);
	if ((t & 0x1f) > 29 || ((t >> 5) & 0x1f) > 29)
	  {
	    this.mode = IB_BAD;
	    z.msg = "too many length or distance symbols";
	    r = Z_DATA_ERROR;

	    this.bitb=b; this.bitk=k; 
	    z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
	    this.write=q;
	    return this.inflate_flush(z,r);
	  }
	t = 258 + (t & 0x1f) + ((t >> 5) & 0x1f);
	if(this.blens==null || this.blens.length<t){
	    this.blens=new Int32Array(t);
	}
	else{
	  for(var i=0; i<t; i++){
              this.blens[i]=0;
          }
	}

	{b>>>=(14);k-=(14);}

	this.index = 0;
	mode = IB_BTREE;
      case IB_BTREE:
	while (this.index < 4 + (this.table >>> 10)){
	  while(k<(3)){
	    if(n!=0){
	      r=Z_OK;
	    }
	    else{
	      this.bitb=b; this.bitk=k; 
	      z.avail_in=n;
	      z.total_in+=p-z.next_in_index;z.next_in_index=p;
	      this.write=q;
	      return this.inflate_flush(z,r);
	    };
	    n--;
	    b|=(z.next_in[p++]&0xff)<<k;
	    k+=8;
	  }

	  this.blens[INFBLOCKS_BORDER[this.index++]] = b&7;

	  {b>>>=(3);k-=(3);}
	}

	while(this.index < 19){
	  this.blens[INFBLOCKS_BORDER[this.index++]] = 0;
	}

	this.bb[0] = 7;
	t = this.inftree.inflate_trees_bits(this.blens, this.bb, this.tb, this.hufts, z);
	if (t != Z_OK){
	  r = t;
	  if (r == Z_DATA_ERROR){
	    this.blens=null;
	    this.mode = IB_BAD;
	  }

	  this.bitb=b; this.bitk=k; 
	  z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
	  write=q;
	  return this.inflate_flush(z,r);
	}

	this.index = 0;
	this.mode = IB_DTREE;
      case IB_DTREE:
	while (true){
	  t = this.table;
	  if(!(this.index < 258 + (t & 0x1f) + ((t >> 5) & 0x1f))){
	    break;
	  }

	  var h; //int[]
	  var i, j, c;

	  t = this.bb[0];

	  while(k<(t)){
	    if(n!=0){
	      r=Z_OK;
	    }
	    else{
	      this.bitb=b; this.bitk=k; 
	      z.avail_in=n;
	      z.total_in+=p-z.next_in_index;z.next_in_index=p;
	      this.write=q;
	      return this.inflate_flush(z,r);
	    };
	    n--;
	    b|=(z.next_in[p++]&0xff)<<k;
	    k+=8;
	  }

//	  if (this.tb[0]==-1){
//            dlog("null...");
//	  }

	  t=this.hufts[(this.tb[0]+(b & inflate_mask[t]))*3+1];
	  c=this.hufts[(this.tb[0]+(b & inflate_mask[t]))*3+2];

	  if (c < 16){
	    b>>>=(t);k-=(t);
	    this.blens[this.index++] = c;
	  }
	  else { // c == 16..18
	    i = c == 18 ? 7 : c - 14;
	    j = c == 18 ? 11 : 3;

	    while(k<(t+i)){
	      if(n!=0){
		r=Z_OK;
	      }
	      else{
		this.bitb=b; this.bitk=k; 
		z.avail_in=n;
		z.total_in+=p-z.next_in_index;z.next_in_index=p;
		this.write=q;
		return this.inflate_flush(z,r);
	      };
	      n--;
	      b|=(z.next_in[p++]&0xff)<<k;
	      k+=8;
	    }

	    b>>>=(t);k-=(t);

	    j += (b & inflate_mask[i]);

	    b>>>=(i);k-=(i);

	    i = this.index;
	    t = this.table;
	    if (i + j > 258 + (t & 0x1f) + ((t >> 5) & 0x1f) ||
		(c == 16 && i < 1)){
	      this.blens=null;
	      this.mode = IB_BAD;
	      z.msg = "invalid bit length repeat";
	      r = Z_DATA_ERROR;

	      this.bitb=b; this.bitk=k; 
	      z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
	      this.write=q;
	      return this.inflate_flush(z,r);
	    }

	    c = c == 16 ? this.blens[i-1] : 0;
	    do{
	      this.blens[i++] = c;
	    }
	    while (--j!=0);
	    this.index = i;
	  }
	}

	this.tb[0]=-1;
	{
	    var bl=new Int32Array(1);
	    var bd=new Int32Array(1);
	    var tl=new Int32Array(1);
	    var td=new Int32Array(1);
	    bl[0] = 9;         // must be <= 9 for lookahead assumptions
	    bd[0] = 6;         // must be <= 9 for lookahead assumptions

	    t = this.table;
	    t = this.inftree.inflate_trees_dynamic(257 + (t & 0x1f), 
					      1 + ((t >> 5) & 0x1f),
					      this.blens, bl, bd, tl, td, this.hufts, z);

	    if (t != Z_OK){
	        if (t == Z_DATA_ERROR){
	            this.blens=null;
	            this.mode = BAD;
	        }
	        r = t;

	        this.bitb=b; this.bitk=k; 
	        z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
	        this.write=q;
	        return this.inflate_flush(z,r);
	    }
	    this.codes.init(bl[0], bd[0], this.hufts, tl[0], this.hufts, td[0], z);
	}
	this.mode = IB_CODES;
      case IB_CODES:
	this.bitb=b; this.bitk=k;
	z.avail_in=n; z.total_in+=p-z.next_in_index;z.next_in_index=p;
	this.write=q;

	if ((r = this.codes.proc(this, z, r)) != Z_STREAM_END){
	  return this.inflate_flush(z, r);
	}
	r = Z_OK;
	this.codes.free(z);

	p=z.next_in_index; n=z.avail_in;b=this.bitb;k=this.bitk;
	q=this.write;m = (q < this.read ? this.read-q-1 : this.end-q);

	if (this.last==0){
	  this.mode = IB_TYPE;
	  break;
	}
	this.mode = IB_DRY;
      case IB_DRY:
	this.write=q; 
	r = this.inflate_flush(z, r); 
	q=this.write; m = (q < this.read ? this.read-q-1 : this.end-q);
	if (this.read != this.write){
	  this.bitb=b; this.bitk=k; 
	  z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
	  this.write=q;
	  return this.inflate_flush(z, r);
	}
	mode = DONE;
      case IB_DONE:
	r = Z_STREAM_END;

	this.bitb=b; this.bitk=k; 
	z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
	this.write=q;
	return this.inflate_flush(z, r);
      case IB_BAD:
	r = Z_DATA_ERROR;

	this.bitb=b; this.bitk=k; 
	z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
	this.write=q;
	return this.inflate_flush(z, r);

      default:
	r = Z_STREAM_ERROR;

	this.bitb=b; this.bitk=k; 
	z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
	this.write=q;
	return this.inflate_flush(z, r);
      }
    }
  }

InfBlocks.prototype.free = function(z){
    this.reset(z, null);
    this.window=null;
    this.hufts=null;
}

InfBlocks.prototype.set_dictionary = function(d, start, n){
    arrayCopy(d, start, window, 0, n);
    this.read = this.write = n;
}

  // Returns true if inflate is currently at the end of a block generated
  // by Z_SYNC_FLUSH or Z_FULL_FLUSH. 
InfBlocks.prototype.sync_point = function(){
    return this.mode == IB_LENS;
}

  // copy as much as possible from the sliding window to the output area
InfBlocks.prototype.inflate_flush = function(z, r){
    var n;
    var p;
    var q;

    // local copies of source and destination pointers
    p = z.next_out_index;
    q = this.read;

    // compute number of bytes to copy as far as end of window
    n = ((q <= this.write ? this.write : this.end) - q);
    if (n > z.avail_out) n = z.avail_out;
    if (n!=0 && r == Z_BUF_ERROR) r = Z_OK;

    // update counters
    z.avail_out -= n;
    z.total_out += n;

    // update check information
    if(this.checkfn != null)
      z.adler=this.check=z._adler.adler32(this.check, this.window, q, n);

    // copy as far as end of window
    arrayCopy(this.window, q, z.next_out, p, n);
    p += n;
    q += n;

    // see if more to copy at beginning of window
    if (q == this.end){
      // wrap pointers
      q = 0;
      if (this.write == this.end)
        this.write = 0;

      // compute bytes to copy
      n = this.write - q;
      if (n > z.avail_out) n = z.avail_out;
      if (n!=0 && r == Z_BUF_ERROR) r = Z_OK;

      // update counters
      z.avail_out -= n;
      z.total_out += n;

      // update check information
      if(this.checkfn != null)
	z.adler=this.check=z._adler.adler32(this.check, this.window, q, n);

      // copy
      arrayCopy(this.window, q, z.next_out, p, n);
      p += n;
      q += n;
    }

    // update pointers
    z.next_out_index = p;
    this.read = q;

    // done
    return r;
  }

//
// InfCodes.java
//

var IC_START=0;  // x: set up for LEN
var IC_LEN=1;    // i: get length/literal/eob next
var IC_LENEXT=2; // i: getting length extra (have base)
var IC_DIST=3;   // i: get distance next
var IC_DISTEXT=4;// i: getting distance extra
var IC_COPY=5;   // o: copying bytes in window, waiting for space
var IC_LIT=6;    // o: got literal, waiting for output space
var IC_WASH=7;   // o: got eob, possibly still output waiting
var IC_END=8;    // x: got eob and all data flushed
var IC_BADCODE=9;// x: got error

function InfCodes() {
}

InfCodes.prototype.init = function(bl, bd, tl, tl_index, td, td_index, z) {
    this.mode=IC_START;
    this.lbits=bl;
    this.dbits=bd;
    this.ltree=tl;
    this.ltree_index=tl_index;
    this.dtree = td;
    this.dtree_index=td_index;
    this.tree=null;
}

InfCodes.prototype.proc = function(s, z, r){ 
    var j;              // temporary storage
    var t;              // temporary pointer (int[])
    var tindex;         // temporary pointer
    var e;              // extra bits or operation
    var b=0;            // bit buffer
    var k=0;            // bits in bit buffer
    var p=0;            // input data pointer
    var n;              // bytes available there
    var q;              // output window write pointer
    var m;              // bytes to end of window or read pointer
    var f;              // pointer to copy strings from

    // copy input/output information to locals (UPDATE macro restores)
    p=z.next_in_index;n=z.avail_in;b=s.bitb;k=s.bitk;
    q=s.write;m=q<s.read?s.read-q-1:s.end-q;

    // process input and output based on current state
    while (true){
      switch (this.mode){
	// waiting for "i:"=input, "o:"=output, "x:"=nothing
      case IC_START:         // x: set up for LEN
	if (m >= 258 && n >= 10){

	  s.bitb=b;s.bitk=k;
	  z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
	  s.write=q;
	  r = this.inflate_fast(this.lbits, this.dbits, 
			   this.ltree, this.ltree_index, 
			   this.dtree, this.dtree_index,
			   s, z);

	  p=z.next_in_index;n=z.avail_in;b=s.bitb;k=s.bitk;
	  q=s.write;m=q<s.read?s.read-q-1:s.end-q;

	  if (r != Z_OK){
	    this.mode = r == Z_STREAM_END ? IC_WASH : IC_BADCODE;
	    break;
	  }
	}
	this.need = this.lbits;
	this.tree = this.ltree;
	this.tree_index=this.ltree_index;

	this.mode = IC_LEN;
      case IC_LEN:           // i: get length/literal/eob next
	j = this.need;

	while(k<(j)){
	  if(n!=0)r=Z_OK;
	  else{

	    s.bitb=b;s.bitk=k;
	    z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
	    s.write=q;
	    return s.inflate_flush(z,r);
	  }
	  n--;
	  b|=(z.next_in[p++]&0xff)<<k;
	  k+=8;
	}

	tindex=(this.tree_index+(b&inflate_mask[j]))*3;

	b>>>=(this.tree[tindex+1]);
	k-=(this.tree[tindex+1]);

	e=this.tree[tindex];

	if(e == 0){               // literal
	  this.lit = this.tree[tindex+2];
	  this.mode = IC_LIT;
	  break;
	}
	if((e & 16)!=0 ){          // length
	  this.get = e & 15;
	  this.len = this.tree[tindex+2];
	  this.mode = IC_LENEXT;
	  break;
	}
	if ((e & 64) == 0){        // next table
	  this.need = e;
	  this.tree_index = tindex/3 + this.tree[tindex+2];
	  break;
	}
	if ((e & 32)!=0){               // end of block
	  this.mode = IC_WASH;
	  break;
	}
	this.mode = IC_BADCODE;        // invalid code
	z.msg = "invalid literal/length code";
	r = Z_DATA_ERROR;

	s.bitb=b;s.bitk=k;
	z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
	s.write=q;
	return s.inflate_flush(z,r);

      case IC_LENEXT:        // i: getting length extra (have base)
	j = this.get;

	while(k<(j)){
	  if(n!=0)r=Z_OK;
	  else{

	    s.bitb=b;s.bitk=k;
	    z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
	    s.write=q;
	    return s.inflate_flush(z,r);
	  }
	  n--; b|=(z.next_in[p++]&0xff)<<k;
	  k+=8;
	}

	this.len += (b & inflate_mask[j]);

	b>>=j;
	k-=j;

	this.need = this.dbits;
	this.tree = this.dtree;
	this.tree_index = this.dtree_index;
	this.mode = IC_DIST;
      case IC_DIST:          // i: get distance next
	j = this.need;

	while(k<(j)){
	  if(n!=0)r=Z_OK;
	  else{

	    s.bitb=b;s.bitk=k;
	    z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
	    s.write=q;
	    return s.inflate_flush(z,r);
	  }
	  n--; b|=(z.next_in[p++]&0xff)<<k;
	  k+=8;
	}

	tindex=(this.tree_index+(b & inflate_mask[j]))*3;

	b>>=this.tree[tindex+1];
	k-=this.tree[tindex+1];

	e = (this.tree[tindex]);
	if((e & 16)!=0){               // distance
	  this.get = e & 15;
	  this.dist = this.tree[tindex+2];
	  this.mode = IC_DISTEXT;
	  break;
	}
	if ((e & 64) == 0){        // next table
	  this.need = e;
	  this.tree_index = tindex/3 + this.tree[tindex+2];
	  break;
	}
	this.mode = IC_BADCODE;        // invalid code
	z.msg = "invalid distance code";
	r = Z_DATA_ERROR;

	s.bitb=b;s.bitk=k;
	z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
	s.write=q;
	return s.inflate_flush(z,r);

      case IC_DISTEXT:       // i: getting distance extra
	j = this.get;

	while(k<(j)){
	  if(n!=0)r=Z_OK;
	  else{

	    s.bitb=b;s.bitk=k;
	    z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
	    s.write=q;
	    return s.inflate_flush(z,r);
	  }
	  n--; b|=(z.next_in[p++]&0xff)<<k;
	  k+=8;
	}

	this.dist += (b & inflate_mask[j]);

	b>>=j;
	k-=j;

	this.mode = IC_COPY;
      case IC_COPY:          // o: copying bytes in window, waiting for space
        f = q - this.dist;
        while(f < 0){     // modulo window size-"while" instead
          f += s.end;     // of "if" handles invalid distances
	}
	while (this.len!=0){

	  if(m==0){
	    if(q==s.end&&s.read!=0){q=0;m=q<s.read?s.read-q-1:s.end-q;}
	    if(m==0){
	      s.write=q; r=s.inflate_flush(z,r);
	      q=s.write;m=q<s.read?s.read-q-1:s.end-q;

	      if(q==s.end&&s.read!=0){q=0;m=q<s.read?s.read-q-1:s.end-q;}

	      if(m==0){
		s.bitb=b;s.bitk=k;
		z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
		s.write=q;
		return s.inflate_flush(z,r);
	      }  
	    }
	  }

	  s.window[q++]=s.window[f++]; m--;

	  if (f == s.end)
            f = 0;
	  this.len--;
	}
	this.mode = IC_START;
	break;
      case IC_LIT:           // o: got literal, waiting for output space
	if(m==0){
	  if(q==s.end&&s.read!=0){q=0;m=q<s.read?s.read-q-1:s.end-q;}
	  if(m==0){
	    s.write=q; r=s.inflate_flush(z,r);
	    q=s.write;m=q<s.read?s.read-q-1:s.end-q;

	    if(q==s.end&&s.read!=0){q=0;m=q<s.read?s.read-q-1:s.end-q;}
	    if(m==0){
	      s.bitb=b;s.bitk=k;
	      z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
	      s.write=q;
	      return s.inflate_flush(z,r);
	    }
	  }
	}
	r=Z_OK;

	s.window[q++]=this.lit; m--;

	this.mode = IC_START;
	break;
      case IC_WASH:           // o: got eob, possibly more output
	if (k > 7){        // return unused byte, if any
	  k -= 8;
	  n++;
	  p--;             // can always return one
	}

	s.write=q; r=s.inflate_flush(z,r);
	q=s.write;m=q<s.read?s.read-q-1:s.end-q;

	if (s.read != s.write){
	  s.bitb=b;s.bitk=k;
	  z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
	  s.write=q;
	  return s.inflate_flush(z,r);
	}
	this.mode = IC_END;
      case IC_END:
	r = Z_STREAM_END;
	s.bitb=b;s.bitk=k;
	z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
	s.write=q;
	return s.inflate_flush(z,r);

      case IC_BADCODE:       // x: got error

	r = Z_DATA_ERROR;

	s.bitb=b;s.bitk=k;
	z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
	s.write=q;
	return s.inflate_flush(z,r);

      default:
	r = Z_STREAM_ERROR;

	s.bitb=b;s.bitk=k;
	z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
	s.write=q;
	return s.inflate_flush(z,r);
      }
    }
  }

InfCodes.prototype.free = function(z){
    //  ZFREE(z, c);
}

  // Called with number of bytes left to write in window at least 258
  // (the maximum string length) and number of input bytes available
  // at least ten.  The ten bytes are six bytes for the longest length/
  // distance pair plus four bytes for overloading the bit buffer.

InfCodes.prototype.inflate_fast = function(bl, bd, tl, tl_index, td, td_index, s, z) {
    var t;                // temporary pointer
    var   tp;             // temporary pointer (int[])
    var tp_index;         // temporary pointer
    var e;                // extra bits or operation
    var b;                // bit buffer
    var k;                // bits in bit buffer
    var p;                // input data pointer
    var n;                // bytes available there
    var q;                // output window write pointer
    var m;                // bytes to end of window or read pointer
    var ml;               // mask for literal/length tree
    var md;               // mask for distance tree
    var c;                // bytes to copy
    var d;                // distance back to copy from
    var r;                // copy source pointer

    var tp_index_t_3;     // (tp_index+t)*3

    // load input, output, bit values
    p=z.next_in_index;n=z.avail_in;b=s.bitb;k=s.bitk;
    q=s.write;m=q<s.read?s.read-q-1:s.end-q;

    // initialize masks
    ml = inflate_mask[bl];
    md = inflate_mask[bd];

    // do until not enough input or output space for fast loop
    do {                          // assume called with m >= 258 && n >= 10
      // get literal/length code
      while(k<(20)){              // max bits for literal/length code
	n--;
	b|=(z.next_in[p++]&0xff)<<k;k+=8;
      }

      t= b&ml;
      tp=tl; 
      tp_index=tl_index;
      tp_index_t_3=(tp_index+t)*3;
      if ((e = tp[tp_index_t_3]) == 0){
	b>>=(tp[tp_index_t_3+1]); k-=(tp[tp_index_t_3+1]);

	s.window[q++] = tp[tp_index_t_3+2];
	m--;
	continue;
      }
      do {

	b>>=(tp[tp_index_t_3+1]); k-=(tp[tp_index_t_3+1]);

	if((e&16)!=0){
	  e &= 15;
	  c = tp[tp_index_t_3+2] + (b & inflate_mask[e]);

	  b>>=e; k-=e;

	  // decode distance base of block to copy
	  while(k<(15)){           // max bits for distance code
	    n--;
	    b|=(z.next_in[p++]&0xff)<<k;k+=8;
	  }

	  t= b&md;
	  tp=td;
	  tp_index=td_index;
          tp_index_t_3=(tp_index+t)*3;
	  e = tp[tp_index_t_3];

	  do {

	    b>>=(tp[tp_index_t_3+1]); k-=(tp[tp_index_t_3+1]);

	    if((e&16)!=0){
	      // get extra bits to add to distance base
	      e &= 15;
	      while(k<(e)){         // get extra bits (up to 13)
		n--;
		b|=(z.next_in[p++]&0xff)<<k;k+=8;
	      }

	      d = tp[tp_index_t_3+2] + (b&inflate_mask[e]);

	      b>>=(e); k-=(e);

	      // do the copy
	      m -= c;
	      if (q >= d){                // offset before dest
		//  just copy
		r=q-d;
		if(q-r>0 && 2>(q-r)){           
		  s.window[q++]=s.window[r++]; // minimum count is three,
		  s.window[q++]=s.window[r++]; // so unroll loop a little
		  c-=2;
		}
		else{
		  s.window[q++]=s.window[r++]; // minimum count is three,
		  s.window[q++]=s.window[r++]; // so unroll loop a little
		  c-=2;
		}
	      }
	      else{                  // else offset after destination
                r=q-d;
                do{
                  r+=s.end;          // force pointer in window
                }while(r<0);         // covers invalid distances
		e=s.end-r;
		if(c>e){             // if source crosses,
		  c-=e;              // wrapped copy
		  if(q-r>0 && e>(q-r)){           
		    do{s.window[q++] = s.window[r++];}
		    while(--e!=0);
		  }
		  else{
		    arrayCopy(s.window, r, s.window, q, e);
		    q+=e; r+=e; e=0;
		  }
		  r = 0;                  // copy rest from start of window
		}

	      }

	      // copy all or what's left
              do{s.window[q++] = s.window[r++];}
		while(--c!=0);
	      break;
	    }
	    else if((e&64)==0){
	      t+=tp[tp_index_t_3+2];
	      t+=(b&inflate_mask[e]);
	      tp_index_t_3=(tp_index+t)*3;
	      e=tp[tp_index_t_3];
	    }
	    else{
	      z.msg = "invalid distance code";

	      c=z.avail_in-n;c=(k>>3)<c?k>>3:c;n+=c;p-=c;k-=c<<3;

	      s.bitb=b;s.bitk=k;
	      z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
	      s.write=q;

	      return Z_DATA_ERROR;
	    }
	  }
	  while(true);
	  break;
	}

	if((e&64)==0){
	  t+=tp[tp_index_t_3+2];
	  t+=(b&inflate_mask[e]);
	  tp_index_t_3=(tp_index+t)*3;
	  if((e=tp[tp_index_t_3])==0){

	    b>>=(tp[tp_index_t_3+1]); k-=(tp[tp_index_t_3+1]);

	    s.window[q++]=tp[tp_index_t_3+2];
	    m--;
	    break;
	  }
	}
	else if((e&32)!=0){

	  c=z.avail_in-n;c=(k>>3)<c?k>>3:c;n+=c;p-=c;k-=c<<3;
 
	  s.bitb=b;s.bitk=k;
	  z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
	  s.write=q;

	  return Z_STREAM_END;
	}
	else{
	  z.msg="invalid literal/length code";

	  c=z.avail_in-n;c=(k>>3)<c?k>>3:c;n+=c;p-=c;k-=c<<3;

	  s.bitb=b;s.bitk=k;
	  z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
	  s.write=q;

	  return Z_DATA_ERROR;
	}
      } 
      while(true);
    } 
    while(m>=258 && n>= 10);

    // not enough input or output--restore pointers and return
    c=z.avail_in-n;c=(k>>3)<c?k>>3:c;n+=c;p-=c;k-=c<<3;

    s.bitb=b;s.bitk=k;
    z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
    s.write=q;

    return Z_OK;
}

//
// InfTree.java
//

function InfTree() {
}

InfTree.prototype.huft_build = function(b, bindex, n, s, d, e, t, m, hp, hn, v) {

    // Given a list of code lengths and a maximum table size, make a set of
    // tables to decode that set of codes.  Return Z_OK on success, Z_BUF_ERROR
    // if the given code set is incomplete (the tables are still built in this
    // case), Z_DATA_ERROR if the input is invalid (an over-subscribed set of
    // lengths), or Z_MEM_ERROR if not enough memory.

    var a;                       // counter for codes of length k
    var f;                       // i repeats in table every f entries
    var g;                       // maximum code length
    var h;                       // table level
    var i;                       // counter, current code
    var j;                       // counter
    var k;                       // number of bits in current code
    var l;                       // bits per table (returned in m)
    var mask;                    // (1 << w) - 1, to avoid cc -O bug on HP
    var p;                       // pointer into c[], b[], or v[]
    var q;                       // points to current table
    var w;                       // bits before this table == (l * h)
    var xp;                      // pointer into x
    var y;                       // number of dummy codes added
    var z;                       // number of entries in current table

    // Generate counts for each bit length

    p = 0; i = n;
    do {
      this.c[b[bindex+p]]++; p++; i--;   // assume all entries <= BMAX
    }while(i!=0);

    if(this.c[0] == n){                // null input--all zero length codes
      t[0] = -1;
      m[0] = 0;
      return Z_OK;
    }

    // Find minimum and maximum length, bound *m by those
    l = m[0];
    for (j = 1; j <= BMAX; j++)
      if(this.c[j]!=0) break;
    k = j;                        // minimum code length
    if(l < j){
      l = j;
    }
    for (i = BMAX; i!=0; i--){
      if(this.c[i]!=0) break;
    }
    g = i;                        // maximum code length
    if(l > i){
      l = i;
    }
    m[0] = l;

    // Adjust last length count to fill out codes, if needed
    for (y = 1 << j; j < i; j++, y <<= 1){
      if ((y -= this.c[j]) < 0){
        return Z_DATA_ERROR;
      }
    }
    if ((y -= this.c[i]) < 0){
      return Z_DATA_ERROR;
    }
    this.c[i] += y;

    // Generate starting offsets into the value table for each length
    this.x[1] = j = 0;
    p = 1;  xp = 2;
    while (--i!=0) {                 // note that i == g from above
      this.x[xp] = (j += this.c[p]);
      xp++;
      p++;
    }

    // Make a table of values in order of bit lengths
    i = 0; p = 0;
    do {
      if ((j = b[bindex+p]) != 0){
        this.v[this.x[j]++] = i;
      }
      p++;
    }
    while (++i < n);
    n = this.x[g];                     // set n to length of v

    // Generate the Huffman codes and for each, make the table entries
    this.x[0] = i = 0;                 // first Huffman code is zero
    p = 0;                        // grab values in bit order
    h = -1;                       // no tables yet--level -1
    w = -l;                       // bits decoded == (l * h)
    this.u[0] = 0;                     // just to keep compilers happy
    q = 0;                        // ditto
    z = 0;                        // ditto

    // go through the bit lengths (k already is bits in shortest code)
    for (; k <= g; k++){
      a = this.c[k];
      while (a--!=0){
	// here i is the Huffman code of length k bits for value *p
	// make tables up to required level
        while (k > w + l){
          h++;
          w += l;                 // previous table always l bits
	  // compute minimum size table less than or equal to l bits
          z = g - w;
          z = (z > l) ? l : z;        // table size upper limit
          if((f=1<<(j=k-w))>a+1){     // try a k-w bit table
                                      // too few codes for k-w bit table
            f -= a + 1;               // deduct codes from patterns left
            xp = k;
            if(j < z){
              while (++j < z){        // try smaller tables up to z bits
                if((f <<= 1) <= this.c[++xp])
                  break;              // enough codes to use up j bits
                f -= this.c[xp];           // else deduct codes from patterns
              }
	    }
          }
          z = 1 << j;                 // table entries for j-bit table

	  // allocate new table
          if (this.hn[0] + z > MANY){       // (note: doesn't matter for fixed)
            return Z_DATA_ERROR;       // overflow of MANY
          }
          this.u[h] = q = /*hp+*/ this.hn[0];   // DEBUG
          this.hn[0] += z;
 
	  // connect to last table, if there is one
	  if(h!=0){
            this.x[h]=i;           // save pattern for backing up
            this.r[0]=j;     // bits in this table
            this.r[1]=l;     // bits to dump before this table
            j=i>>>(w - l);
            this.r[2] = (q - this.u[h-1] - j);               // offset to this table
            arrayCopy(this.r, 0, hp, (this.u[h-1]+j)*3, 3); // connect to last table
          }
          else{
            t[0] = q;               // first table is returned result
	  }
        }

	// set up table entry in r
        this.r[1] = (k - w);
        if (p >= n){
          this.r[0] = 128 + 64;      // out of values--invalid code
	}
        else if (v[p] < s){
          this.r[0] = (this.v[p] < 256 ? 0 : 32 + 64);  // 256 is end-of-block
          this.r[2] = this.v[p++];          // simple code is just the value
        }
        else{
          this.r[0]=(e[this.v[p]-s]+16+64); // non-simple--look up in lists
          this.r[2]=d[this.v[p++] - s];
        }

        // fill code-like entries with r
        f=1<<(k-w);
        for (j=i>>>w;j<z;j+=f){
          arrayCopy(this.r, 0, hp, (q+j)*3, 3);
	}

	// backwards increment the k-bit code i
        for (j = 1 << (k - 1); (i & j)!=0; j >>>= 1){
          i ^= j;
	}
        i ^= j;

	// backup over finished tables
        mask = (1 << w) - 1;      // needed on HP, cc -O bug
        while ((i & mask) != this.x[h]){
          h--;                    // don't need to update q
          w -= l;
          mask = (1 << w) - 1;
        }
      }
    }
    // Return Z_BUF_ERROR if we were given an incomplete table
    return y != 0 && g != 1 ? Z_BUF_ERROR : Z_OK;
}

InfTree.prototype.inflate_trees_bits = function(c, bb, tb, hp, z) {
    var result;
    this.initWorkArea(19);
    this.hn[0]=0;
    result = this.huft_build(c, 0, 19, 19, null, null, tb, bb, hp, this.hn, this.v);

    if(result == Z_DATA_ERROR){
      z.msg = "oversubscribed dynamic bit lengths tree";
    }
    else if(result == Z_BUF_ERROR || bb[0] == 0){
      z.msg = "incomplete dynamic bit lengths tree";
      result = Z_DATA_ERROR;
    }
    return result;
}

InfTree.prototype.inflate_trees_dynamic = function(nl, nd, c, bl, bd, tl, td, hp, z) {
    var result;

    // build literal/length tree
    this.initWorkArea(288);
    this.hn[0]=0;
    result = this.huft_build(c, 0, nl, 257, cplens, cplext, tl, bl, hp, this.hn, this.v);
    if (result != Z_OK || bl[0] == 0){
      if(result == Z_DATA_ERROR){
        z.msg = "oversubscribed literal/length tree";
      }
      else if (result != Z_MEM_ERROR){
        z.msg = "incomplete literal/length tree";
        result = Z_DATA_ERROR;
      }
      return result;
    }

    // build distance tree
    this.initWorkArea(288);
    result = this.huft_build(c, nl, nd, 0, cpdist, cpdext, td, bd, hp, this.hn, this.v);

    if (result != Z_OK || (bd[0] == 0 && nl > 257)){
      if (result == Z_DATA_ERROR){
        z.msg = "oversubscribed distance tree";
      }
      else if (result == Z_BUF_ERROR) {
        z.msg = "incomplete distance tree";
        result = Z_DATA_ERROR;
      }
      else if (result != Z_MEM_ERROR){
        z.msg = "empty distance tree with lengths";
        result = Z_DATA_ERROR;
      }
      return result;
    }

    return Z_OK;
}
/*
  static int inflate_trees_fixed(int[] bl,  //literal desired/actual bit depth
                                 int[] bd,  //distance desired/actual bit depth
                                 int[][] tl,//literal/length tree result
                                 int[][] td,//distance tree result 
                                 ZStream z  //for memory allocation
				 ){

*/

function inflate_trees_fixed(bl, bd, tl, td, z) {
    bl[0]=fixed_bl;
    bd[0]=fixed_bd;
    tl[0]=fixed_tl;
    td[0]=fixed_td;
    return Z_OK;
}

InfTree.prototype.initWorkArea = function(vsize){
    if(this.hn==null){
        this.hn=new Int32Array(1);
        this.v=new Int32Array(vsize);
        this.c=new Int32Array(BMAX+1);
        this.r=new Int32Array(3);
        this.u=new Int32Array(BMAX);
        this.x=new Int32Array(BMAX+1);
    }
    if(this.v.length<vsize){ 
        this.v=new Int32Array(vsize); 
    }
    for(var i=0; i<vsize; i++){this.v[i]=0;}
    for(var i=0; i<BMAX+1; i++){this.c[i]=0;}
    for(var i=0; i<3; i++){this.r[i]=0;}
//  for(int i=0; i<BMAX; i++){u[i]=0;}
    arrayCopy(this.c, 0, this.u, 0, BMAX);
//  for(int i=0; i<BMAX+1; i++){x[i]=0;}
    arrayCopy(this.c, 0, this.x, 0, BMAX+1);
}

var testArray = new Uint8Array(1);
var hasSubarray = (typeof testArray.subarray === 'function');
var hasSlice = false; /* (typeof testArray.slice === 'function'); */ // Chrome slice performance is so dire that we're currently not using it...

function arrayCopy(src, srcOffset, dest, destOffset, count) {
    if (count == 0) {
        return;
    } 
    if (!src) {
        throw "Undef src";
    } else if (!dest) {
        throw "Undef dest";
    }

    if (srcOffset == 0 && count == src.length) {
        arrayCopy_fast(src, dest, destOffset);
    } else if (hasSubarray) {
        arrayCopy_fast(src.subarray(srcOffset, srcOffset + count), dest, destOffset); 
    } else if (src.BYTES_PER_ELEMENT == 1 && count > 100) {
        arrayCopy_fast(new Uint8Array(src.buffer, src.byteOffset + srcOffset, count), dest, destOffset);
    } else { 
        arrayCopy_slow(src, srcOffset, dest, destOffset, count);
    }

}

function arrayCopy_slow(src, srcOffset, dest, destOffset, count) {

    // dlog('_slow call: srcOffset=' + srcOffset + '; destOffset=' + destOffset + '; count=' + count);

     for (var i = 0; i < count; ++i) {
        dest[destOffset + i] = src[srcOffset + i];
    }
}

function arrayCopy_fast(src, dest, destOffset) {
    dest.set(src, destOffset);
}


  // largest prime smaller than 65536
var ADLER_BASE=65521; 
  // NMAX is the largest n such that 255n(n+1)/2 + (n+1)(BASE-1) <= 2^32-1
var ADLER_NMAX=5552;

function adler32(adler, /* byte[] */ buf,  index, len){
    if(buf == null){ return 1; }

    var s1=adler&0xffff;
    var s2=(adler>>16)&0xffff;
    var k;

    while(len > 0) {
      k=len<ADLER_NMAX?len:ADLER_NMAX;
      len-=k;
      while(k>=16){
        s1+=buf[index++]&0xff; s2+=s1;
        s1+=buf[index++]&0xff; s2+=s1;
        s1+=buf[index++]&0xff; s2+=s1;
        s1+=buf[index++]&0xff; s2+=s1;
        s1+=buf[index++]&0xff; s2+=s1;
        s1+=buf[index++]&0xff; s2+=s1;
        s1+=buf[index++]&0xff; s2+=s1;
        s1+=buf[index++]&0xff; s2+=s1;
        s1+=buf[index++]&0xff; s2+=s1;
        s1+=buf[index++]&0xff; s2+=s1;
        s1+=buf[index++]&0xff; s2+=s1;
        s1+=buf[index++]&0xff; s2+=s1;
        s1+=buf[index++]&0xff; s2+=s1;
        s1+=buf[index++]&0xff; s2+=s1;
        s1+=buf[index++]&0xff; s2+=s1;
        s1+=buf[index++]&0xff; s2+=s1;
        k-=16;
      }
      if(k!=0){
        do{
          s1+=buf[index++]&0xff; s2+=s1;
        }
        while(--k!=0);
      }
      s1%=ADLER_BASE;
      s2%=ADLER_BASE;
    }
    return (s2<<16)|s1;
}



function jszlib_inflate_buffer(buffer, start, length, afterUncOffset) {
    if (!start) {
        buffer = new Uint8Array(buffer);
    } else {
        buffer = new Uint8Array(buffer, start, length);
    }

    var z = new ZStream();
    z.inflateInit(DEF_WBITS, true);
    z.next_in = buffer;
    z.next_in_index = 0;
    z.avail_in = buffer.length;

    var oBlockList = [];
    var totalSize = 0;
    while (true) {
        var obuf = new Uint8Array(32000);
        z.next_out = obuf;
        z.next_out_index = 0;
        z.avail_out = obuf.length;
        var status = z.inflate(Z_NO_FLUSH);
        if (status != Z_OK && status != Z_STREAM_END && status != Z_BUF_ERROR) {
            throw z.msg;
        }
        if (z.avail_out != 0) {
            var newob = new Uint8Array(obuf.length - z.avail_out);
            arrayCopy(obuf, 0, newob, 0, (obuf.length - z.avail_out));
            obuf = newob;
        }
        oBlockList.push(obuf);
        totalSize += obuf.length;
        if (status == Z_STREAM_END || status == Z_BUF_ERROR) {
            break;
        }
    }

    if (afterUncOffset) {
        afterUncOffset[0] = (start || 0) + z.next_in_index;
    }

    if (oBlockList.length == 1) {
        return oBlockList[0].buffer;
    } else {
        var out = new Uint8Array(totalSize);
        var cursor = 0;
        for (var i = 0; i < oBlockList.length; ++i) {
            var b = oBlockList[i];
            arrayCopy(b, 0, out, cursor, b.length);
            cursor += b.length;
        }
        return out.buffer;
    }
}/*
html5slider - a JS implementation of <input type=range> for Firefox 16 and up
https://github.com/fryn/html5slider

Copyright (c) 2010-2012 Frank Yan, <http://frankyan.com>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/

(function() {

// test for native support
var test = document.createElement('input');
try {
  test.type = 'range';
  if (test.type == 'range')
    return;
} catch (e) {
  return;
}

// test for required property support
test.style.background = 'linear-gradient(red, red)';
if (!test.style.backgroundImage || !('MozAppearance' in test.style) ||
    !document.mozSetImageElement || !this.MutationObserver)
  return;

var scale;
var isMac = navigator.platform == 'MacIntel';
var thumb = {
  radius: isMac ? 9 : 6,
  width: isMac ? 22 : 12,
  height: isMac ? 16 : 20
};
var track = 'linear-gradient(transparent ' + (isMac ?
  '6px, #999 6px, #999 7px, #ccc 8px, #bbb 9px, #bbb 10px, transparent 10px' :
  '9px, #999 9px, #bbb 10px, #fff 11px, transparent 11px') +
  ', transparent)';
var styles = {
  'min-width': thumb.width + 'px',
  'min-height': thumb.height + 'px',
  'max-height': thumb.height + 'px',
  padding: '0 0 ' + (isMac ? '2px' : '1px'),
  border: 0,
  'border-radius': 0,
  cursor: 'default',
  'text-indent': '-999999px' // -moz-user-select: none; breaks mouse capture
};
var options = {
  attributes: true,
  attributeFilter: ['min', 'max', 'step', 'value']
};
var forEach = Array.prototype.forEach;
var onChange = document.createEvent('HTMLEvents');
onChange.initEvent('change', true, false);

if (document.readyState == 'loading')
  document.addEventListener('DOMContentLoaded', initialize, true);
else
  initialize();

function initialize() {
  // create initial sliders
  forEach.call(document.querySelectorAll('input[type=range]'), transform);
  // create sliders on-the-fly
  new MutationObserver(function(mutations) {
    mutations.forEach(function(mutation) {
      if (mutation.addedNodes)
        forEach.call(mutation.addedNodes, function(node) {
          check(node);
          if (node.childElementCount)
            forEach.call(node.querySelectorAll('input'), check);
        });
    });
  }).observe(document, { childList: true, subtree: true });
}

function check(input) {
  if (input.localName == 'input' && input.type != 'range' &&
      input.getAttribute('type') == 'range')
    transform(input);
}

function transform(slider) {

  var isValueSet, areAttrsSet, isChanged, isClick, prevValue, rawValue, prevX;
  var min, max, step, range, value = slider.value;

  // lazily create shared slider affordance
  if (!scale) {
    scale = document.body.appendChild(document.createElement('hr'));
    style(scale, {
      '-moz-appearance': isMac ? 'scale-horizontal' : 'scalethumb-horizontal',
      display: 'block',
      visibility: 'visible',
      opacity: 1,
      position: 'fixed',
      top: '-999999px'
    });
    document.mozSetImageElement('__sliderthumb__', scale);
  }

  // reimplement value and type properties
  var getValue = function() { return '' + value; };
  var setValue = function setValue(val) {
    value = '' + val;
    isValueSet = true;
    draw();
    delete slider.value;
    slider.value = value;
    slider.__defineGetter__('value', getValue);
    slider.__defineSetter__('value', setValue);
  };
  slider.__defineGetter__('value', getValue);
  slider.__defineSetter__('value', setValue);
  slider.__defineGetter__('type', function() { return 'range'; });

  // sync properties with attributes
  ['min', 'max', 'step'].forEach(function(prop) {
    if (slider.hasAttribute(prop))
      areAttrsSet = true;
    slider.__defineGetter__(prop, function() {
      return this.hasAttribute(prop) ? this.getAttribute(prop) : '';
    });
    slider.__defineSetter__(prop, function(val) {
      val === null ? this.removeAttribute(prop) : this.setAttribute(prop, val);
    });
  });

  // initialize slider
  slider.readOnly = true;
  style(slider, styles);
  update();

  new MutationObserver(function(mutations) {
    mutations.forEach(function(mutation) {
      if (mutation.attributeName != 'value') {
        update();
        areAttrsSet = true;
      }
      // note that value attribute only sets initial value
      else if (!isValueSet) {
        value = slider.getAttribute('value');
        draw();
      }
    });
  }).observe(slider, options);

  slider.addEventListener('mousedown', onDragStart, true);
  slider.addEventListener('keydown', onKeyDown, true);
  slider.addEventListener('focus', onFocus, true);
  slider.addEventListener('blur', onBlur, true);

  function onDragStart(e) {
    isClick = true;
    setTimeout(function() { isClick = false; }, 0);
    if (e.button || !range)
      return;
    var width = parseFloat(getComputedStyle(this, 0).width);
    var multiplier = (width - thumb.width) / range;
    if (!multiplier)
      return;
    // distance between click and center of thumb
    var dev = e.clientX - this.getBoundingClientRect().left - thumb.width / 2 -
              (value - min) * multiplier;
    // if click was not on thumb, move thumb to click location
    if (Math.abs(dev) > thumb.radius) {
      isChanged = true;
      this.value -= -dev / multiplier;
    }
    rawValue = value;
    prevX = e.clientX;
    this.addEventListener('mousemove', onDrag, true);
    this.addEventListener('mouseup', onDragEnd, true);
  }

  function onDrag(e) {
    var width = parseFloat(getComputedStyle(this, 0).width);
    var multiplier = (width - thumb.width) / range;
    if (!multiplier)
      return;
    rawValue += (e.clientX - prevX) / multiplier;
    prevX = e.clientX;
    isChanged = true;
    this.value = rawValue;
  }

  function onDragEnd() {
    this.removeEventListener('mousemove', onDrag, true);
    this.removeEventListener('mouseup', onDragEnd, true);
  }

  function onKeyDown(e) {
    if (e.keyCode > 36 && e.keyCode < 41) { // 37-40: left, up, right, down
      onFocus.call(this);
      isChanged = true;
      this.value = value + (e.keyCode == 38 || e.keyCode == 39 ? step : -step);
    }
  }

  function onFocus() {
    if (!isClick)
      this.style.boxShadow = !isMac ? '0 0 0 2px #fb0' :
        'inset 0 0 20px rgba(0,127,255,.1), 0 0 1px rgba(0,127,255,.4)';
  }

  function onBlur() {
    this.style.boxShadow = '';
  }

  // determines whether value is valid number in attribute form
  function isAttrNum(value) {
    return !isNaN(value) && +value == parseFloat(value);
  }

  // validates min, max, and step attributes and redraws
  function update() {
    min = isAttrNum(slider.min) ? +slider.min : 0;
    max = isAttrNum(slider.max) ? +slider.max : 100;
    if (max < min)
      max = min > 100 ? min : 100;
    step = isAttrNum(slider.step) && slider.step > 0 ? +slider.step : 1;
    range = max - min;
    draw(true);
  }

  // recalculates value property
  function calc() {
    if (!isValueSet && !areAttrsSet)
      value = slider.getAttribute('value');
    if (!isAttrNum(value))
      value = (min + max) / 2;;
    // snap to step intervals (WebKit sometimes does not - bug?)
    value = Math.round((value - min) / step) * step + min;
    if (value < min)
      value = min;
    else if (value > max)
      value = min + ~~(range / step) * step;
  }

  // renders slider using CSS background ;)
  function draw(attrsModified) {
    calc();
    if (isChanged && value != prevValue)
      slider.dispatchEvent(onChange);
    isChanged = false;
    if (!attrsModified && value == prevValue)
      return;
    prevValue = value;
    var position = range ? (value - min) / range * 100 : 0;
    var bg = '-moz-element(#__sliderthumb__) ' + position + '% no-repeat, ';
    style(slider, { background: bg + track });
  }

}

function style(element, styles) {
  for (var prop in styles)
    element.style.setProperty(prop, styles[prop], 'important');
}

})();
