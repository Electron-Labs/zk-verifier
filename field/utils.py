# Extended Euclid Algorithm to implement multiplicative modular inverse
def egcd(a, b):
    if a == 0:
        return (b, 0, 1)

    else:
        g, y, x = egcd(b % a, a)
        return (g, x - (b // a) * y, y)


def mod_inverse(a, m):
    g, x, y = egcd(a, m)
    if g != 1:
        raise Exception('modular inverse does not exist')
    else:
        return x % m


def sub_noborrow(num, z):
    other = [0] * len(num)
    other[0] = z
    borrow = 0
    idx = 0

    for a, b in zip(num, other):
        tmp = (1 << 64) + a - b - borrow
        if (tmp >> 64) == 0:
            borrow = 1
        else:
            borrow = 0
        num[idx] = tmp % 2**64
        idx += 1

    return num


def add_nocarry(num, z):
    other = [0] * len(num)
    other[0] = z
    carry = 0
    idx = 0

    for a, b in zip(num, other):
        tmp = a + b + carry
        carry = (tmp >> 64) % 2**64
        num[idx] = tmp % 2**64
        idx += 1

    return num


def div2(num):
    t = 0
    idx = len(num) - 1

    for i in reversed(num):
        t2 = i << 63
        i >>= 1
        i |= t
        t = t2
        num[idx] = i
        idx -= 1

    return num


def is_zero(num):
    for val in num:
        if val != 0:
            return False
    return True


def is_odd(num):
    if num[0] & 1 == 1:
        return True
    return False


# A helper function to find Non-Adjacent form of a number
def find_naf(num):
    bin_vec = [num]
    result = []

    while not is_zero(bin_vec):
        z = None
        if is_odd(bin_vec):
            z = 2 - (bin_vec[0] % 4)
            if z >= 0:
                bin_vec = sub_noborrow(bin_vec, z)
            else:
                bin_vec = add_nocarry(bin_vec, -z)
        else:
            z = 0
        result.append(z)
        bin_vec = div2(bin_vec)

    return result