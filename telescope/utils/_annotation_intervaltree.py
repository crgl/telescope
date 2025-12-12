# -*- coding: utf-8 -*-
from __future__ import print_function
from __future__ import absolute_import
from builtins import object

import re
from collections import defaultdict, namedtuple, Counter, OrderedDict
import logging as lg
import pickle

import pandas as pd


from intervaltree import Interval, IntervalTree


__author__ = 'Matthew L. Bendall'
__copyright__ = "Copyright (C) 2019 Matthew L. Bendall"


GTFRow = namedtuple('GTFRow', ['chrom','source','feature','start','end','score','strand','frame','attribute'])

def overlap_length(a,b):
    return max(0, min(a.end,b.end) - max(a.begin,b.begin))

def merge_intervals(a, b, d=None):
    return Interval(min(a.begin,b.begin), max(a.end,b.end), d)

class _AnnotationIntervalTree(object):

    def __init__(self, gtf_file, attribute_name, stranded_mode, feature_type='exon'):
        lg.debug('Using intervaltree for annotation.')
        self.loci = set()
        self.key = attribute_name
        self.itree = defaultdict(IntervalTree)
        self.run_stranded = True if stranded_mode != 'None' else False

        # GTF filehandle
        attribute_df = pd.read_csv(gtf_file, sep='\t', comment='#', header=None, names=['chrom','source','feature','start','end','score','strand','frame','attribute'])
        feature_df = attribute_df[(attribute_df['feature'] == feature_type)]
        annotation_df = feature_df[(feature_df['attribute'].str.contains(self.key))].copy()
        skipped_annotations = feature_df[~(feature_df['attribute'].str.contains(self.key))].index
        annotation_df['key'] = annotation_df['attribute'].apply(lambda x: re.search(r'%s\s+"(.+?)";' % self.key, x).group(1))
        annotation_df = annotation_df.drop_duplicates(subset=['chrom','start','end','strand','key']).copy()
        annotation_df['chrom'] = annotation_df['chrom'].astype(str)
        annotation_df['start'] = annotation_df['start'].astype(int)
        annotation_df['end'] = annotation_df['end'].astype(int)
        annotation_df['key'] = annotation_df['key'].astype(str)
        annotation_df['strand'] = annotation_df['strand'].astype(str)
        self.loci = set(annotation_df['key'].unique())
        for rownum in skipped_annotations:
            lg.warning('Skipping row %d: missing attribute "%s"' % (rownum, self.key))
        for f in annotation_df.itertuples(index=False, name='GTFRow'):
            ''' Add to interval tree '''
            new_iv = Interval(f.start, f.end+1, {self.key: f.key, 'strand': f.strand})
            # Merge overlapping intervals from same locus
            overlap = self.itree[str(f.chrom)].overlap(new_iv)
            if len(overlap) > 0:
                mergeable = [iv for iv in overlap if iv.data[self.key]==f.key]
                if mergeable:
                    assert len(mergeable) == 1, "Error"
                    new_iv = merge_intervals(mergeable[0], new_iv, {self.key: f.key, 'strand': f.strand})
                    self.itree[f.chrom].remove(mergeable[0])
            self.itree[f.chrom].add(new_iv)

    def feature_length(self):
        """ Get feature lengths

        Returns:
            (dict of str: int): Feature names to feature lengths

        """
        ret = Counter()
        for chrom in list(self.itree.keys()):
            for iv in list(self.itree[chrom].items()):
                ret[iv.data[self.key]] += iv.length()
        return ret

    def subregion(self, ref, start_pos=None, end_pos=None):
        _subannot = type(self).__new__(type(self))
        _subannot.key = self.key
        _subannot.itree = defaultdict(IntervalTree)
        _subannot.run_stranded = self.run_stranded

        if ref in self.itree:
            _subtree = self.itree[ref].copy()
            if start_pos is not None:
                _subtree.chop(_subtree.begin(), start_pos)
            if end_pos is not None:
                _subtree.chop(end_pos, _subtree.end() + 1)
            _subannot.itree[ref] = _subtree
        return _subannot

    def intersect_blocks(self, ref, blocks, frag_strand):
        _result = Counter()
        for b_start, b_end in blocks:
            query = Interval(b_start, (b_end + 1))
            for iv in self.itree[ref].overlap(query):
                if self.run_stranded == True:
                    if iv.data['strand'] == frag_strand:
                        _result[iv.data[self.key]] += overlap_length(iv, query)
                else:
                    _result[iv.data[self.key]] += overlap_length(iv, query)
        return _result

    def save(self, filename):
        with open(filename, 'wb') as outh:
            pickle.dump({
                'key': self.key,
                'loci': self.loci,
                'itree': self.itree,
                'run_stranded': self.run_stranded,
            }, outh)

    @classmethod
    def load(cls, filename):
        with open(filename, 'rb') as fh:
            loader = pickle.load(fh)
        obj = cls.__new__(cls)
        obj.key = loader['key']
        obj.loci = loader['loci']
        obj.itree = loader['itree']
        obj.run_stranded = loader['run_stranded']

        return obj
