# calc.py
from pybioshell.src.datastructures.vec3 import Vec3

class TestVec3:
    def test_addition(self):
        vec1 = Vec3(1, 2, 3)
        vec2 = Vec3(4, 5, 6)
        result = vec1 + vec2
        assert result.x == 5
        assert result.y == 7
        assert result.z == 9

    def test_subtraction(self):
        vec1 = Vec3(1, 2, 3)
        vec2 = Vec3(4, 5, 6)
        result = vec1 - vec2
        assert result.x == -3
        assert result.y == -3
        assert result.z == -3

    def test_multiplication(self):
        vec = Vec3(1, 2, 3)
        scalar = 2
        result = vec * scalar
        assert result.x == 2
        assert result.y == 4
        assert result.z == 6

    def test_division(self):
        vec = Vec3(4, 6, 8)
        scalar = 2
        result = vec / scalar
        assert result.x == 2
        assert result.y == 3
        assert result.z == 4

    def test_dot_product(self):
        vec1 = Vec3(1, 2, 3)
        vec2 = Vec3(4, 5, 6)
        result = vec1.dot(vec2)
        assert result == 32

    def test_cross_product(self):
        vec1 = Vec3(1, 2, 3)
        vec2 = Vec3(4, 5, 6)
        result = vec1.cross(vec2)
        assert result.x == -3
        assert result.y == 6
        assert result.z == -3