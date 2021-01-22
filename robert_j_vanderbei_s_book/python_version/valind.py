class valind:
    d = 0
    i = 0

    def __init__(self):
        d = 0
        i = 0

    def toString(self):
        str(self.d) + ' ' + str(self.i)

    def println(self):
        print('d={:d}, i={:d}'.format(self.d, self.i), end=' ')
