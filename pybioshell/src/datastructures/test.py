import unittest
from pybioshell.src.datastructures import euclidean_distance_squared, depth_first_inorder, BinaryTreeNode, collect_values
from pybioshell.src.datastructures.test import create_kd_tree, find_nearest, find_within, KdTreeData, mark_partitions
from pybioshell.src.datastructures.vec3 import Vec3
import random

class TestKdTree(unittest.TestCase):

    def test_create_kd_tree(self):
        data = [[0.0, 0.5], [0.7, 0.5], [0.8, 1.1]]
        create_kd_tree(data, 2)

    def test_find_nearest(self):
        N = 128
        random.seed(0)
        data = [[random.random()] for _ in range(N)]
        query = [0.5]
        min_d, min_e = euclidean_distance_squared(query, data[0], 1), data[0]
        for e in data:
            d = euclidean_distance_squared(query, e, 1)
            if d < min_d:
                min_d, min_e = d, e
        root = create_kd_tree(data.copy(), 1)
        node_cnt = 0
        depth_first_inorder(root, lambda _n: (node_cnt := node_cnt + 1))
        self.assertEqual(N, node_cnt)
        d, e = find_nearest(root, query, 1, euclidean_distance_squared)
        self.assertAlmostEqual(d, min_d, delta=0.000001)
        self.assertAlmostEqual(e[0], min_e[0], delta=0.000001)

    def test_find_within(self):
        N = 7
        data = [Vec3(0.0, 0.0, 0.0) for _ in range(N * N)]
        for i in range(N * N):
            data[i].x = (i % N) * 0.1 + 0.1
            data[i].y = (i // N) * 0.1 + 0.1
        root = create_kd_tree(data.copy(), 2)
        query = Vec3(0.3, 0.3, 0.0)
        neighbors = find_within(root, query, 2, 0.021, euclidean_distance_squared)
        self.assertEqual(len(neighbors), 9)

    def test_find_within_3d(self):
        N = 10
        data = [Vec3(0.0, 0.0, 0.0) for _ in range(N * N * N)]
        for k in range(N):
            for j in range(N):
                for i in range(N):
                    idx = N * N * k + N * j + i
                    data[idx].x = i * 0.1 + 0.1
                    data[idx].y = j * 0.1 + 0.1
                    data[idx].z = k * 0.1 + 0.1
        root = create_kd_tree(data.copy(), 3)
        query = Vec3(0.3, 0.3, 0.3)
        neighbors = find_within(root, query, 3, 0.031, euclidean_distance_squared)
        self.assertEqual(len(neighbors), 27)

    def test_mark_partitions(self):
        N = 31
        random.seed(0)
        data = [[2.0 * random.random() - 1.0, 2.0 * random.random() - 1.0] for _ in range(N)]
        root = create_kd_tree(data.copy(), 2)
        mark_partitions(root, 3, 0)
        min_v = [[2.0, 2.0] for _ in range(8)]
        max_v = [[-2.0, -2.0] for _ in range(8)]
        for n in collect_values(root):
            min_v[n.level][0] = min(min_v[n.level][0], n.value[0])
            min_v[n.level][1] = min(min_v[n.level][1], n.value[1])
            max_v[n.level][0] = max(max_v[n.level][0], n.value[0])
            max_v[n.level][1] = max(max_v[n.level][1], n.value[1])
        self.assertGreater(min_v[7][1], max_v[6][1])
        self.assertGreater(min_v[5][1], max_v[4][1])
        self.assertGreater(min_v[7][0], max_v[5][0])
        self.assertGreater(min_v[6][0], max_v[4][0])

if __name__ == '__main__':
    unittest.main()