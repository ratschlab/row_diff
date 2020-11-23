#!/usr/bin/env python3

import struct
import re
from collections import deque
import argparse


class Node:
    def __init__(self, val, left, right):
        self.val = val
        self.left = left
        self.right = right


class HuffmanTree:
    def frequencies(self, log_file):
        log = open(log_file, 'r')
        line = ''
        magic_str = '<Row multiplicity: num unique rows>:';
        for line in log.readlines():
            if magic_str in line:
                break
        if line == '':
            print('No row multiplicity lines found')
            return

        line = line[line.index(magic_str) + len(magic_str):]
        mult_count = [int(s) for s in re.split(',|:|\n|\s', line) if s.isdigit()]
        self.total_size = 0
        probs = []
        for i in range(0, len(mult_count) - 1, 2):
            self.total_size += mult_count[i] * mult_count[i + 1]
        print(f'No of symbols: {self.total_size}')
        # number are pairs of (mult_count : row_count)
        for v in range(len(mult_count) - 1, -1, -2):
            for i in range(mult_count[v]):
                probs.append(mult_count[v - 1] / self.total_size)

        return probs

    def get_min(self, q1, q2):
        if not q1:
            return q2.popleft()
        if not q2:
            return q1.popleft()

        return q1.popleft() if q1[0].val < q2[0].val else q2.popleft()

    def traverse(self, node, depth):
        if not node.left or not node.right:
            return depth * round(self.total_size * node.val)
        else:
            return self.traverse(node.left, depth + 1) + self.traverse(node.right, depth + 1)

    def build_tree(self, filename):
        probs = self.frequencies(filename)
        q1 = deque()
        for p in probs:
            q1.append(Node(p, None, None))

        if len(probs) == 1:
            return 1

        q2 = deque()
        while q1 or len(q2) != 1:
            left = self.get_min(q1, q2)
            right = self.get_min(q1, q2)
            top = Node(left.val + right.val, left, right)
            # print(f'Adding parent with weight {top.val} for {left.val} and {right.val}')
            q2.append(top)

        return self.traverse(q2[0], 0)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('filename', help='Metagraph log to parse')

    args = parser.parse_args()

    tree = HuffmanTree()
    print(f'Encoded size in bits: {tree.build_tree(args.filename)}')
