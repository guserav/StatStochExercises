#! /usr/bin/env python3
import random
from statistics import mean

def trie(k, n):
    s = set()
    i = 0
    while len(s) < k:
        i += 1
        s.add(random.randint(1,n))
    return i

def main():
    n = 10
    l = []
    tries = 100000
    for i in range(n):
        l.append([])
        for t in range(tries):
            l[i].append(trie(i+1,n))
    l = list(map(lambda x: mean(x), l))
    print(l)
    l = []
    for t in range(100000):
        l.append(trie(63,63))
    print(mean(l))

if __name__ == "__main__":
    main()
