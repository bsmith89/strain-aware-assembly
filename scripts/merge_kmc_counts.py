#!/usr/bin/env python3

import sys


def parse_kv_file(inpath):
    with open(inpath) as f:
        for line in f:
            out = line.split()
            yield out

class Stream():
    def __init__(self, iterator):
        self.iterator = iterator
        self.next_key, self.next_value = next(self.iterator)
        self.exhausted = False

    def get_if_next(self, key, default=None):
        if self.exhausted:
            return default
        elif self.next_key == key:
            this_key, this_value = self.next_key, self.next_value
            try:
                self.next_key, self.next_value = next(self.iterator)
            except StopIteration:
                self.exhausted = True
            else:
                assert self.next_key > this_key  # Assert strictly increasing.
            return this_value
        else:
            return default

    def __lt__(self, other):
        if self.exhausted:
            return False
        elif other.exhausted:
            return True
        else:
            return self.next_key < other.next_key


def iter_merged_streams(streams):
    while not all([s.exhausted for s in streams]):
        next_key = min(streams).next_key
        out = [s.get_if_next(next_key, default=0) for s in streams]
        yield next_key, out


if __name__ == "__main__":
    streams = [Stream(parse_kv_file(path)) for path in sys.argv[1:]]
    for key, values in iter_merged_streams(streams):
        print(key, *values, sep='\t')
