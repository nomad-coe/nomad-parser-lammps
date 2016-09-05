import numpy as np

class A:
    def __init__(self):
        self.data = np.array([])

    def update(self, row):
        self.data = np.append(self.data, row)

    def finalize(self):
        return np.reshape(self.data, newshape=(self.data.shape[0]/5, 5))



class B:

    def __init__(self):
        self.data = []

    def update(self, row):
        for r in row:
            self.data.append(r)

    def finalize(self):
        return np.reshape(self.data, newshape=(len(self.data)/5, 5))


class C:

    def __init__(self):
        self.data = np.zeros((100,))
        self.capacity = 100
        self.size = 0

    def update(self, row):
        for r in row:
            self.add(r)

    def add(self, x):
        if self.size == self.capacity:
            self.capacity *= 4
            newdata = np.zeros((self.capacity,))
            newdata[:self.size] = self.data
            self.data = newdata

        self.data[self.size] = x
        self.size += 1

    def finalize(self):
        data = self.data[:self.size]
        return np.reshape(data, newshape=(len(data)/5, 5))



if __name__ == '__main__':
    from timeit import timeit
    def test1():
        x = A()
        for i in range(10000):
            x.update([i])


    print(timeit(test1, number=1))

    def test2():
        x = B()
        for i in range(10000):
            x.update([i])

    print(timeit(test2, number=1))

    def test3():
        x = C()
        for i in range(10000):
            x.update([i])

    print(timeit(test3, number=1))

