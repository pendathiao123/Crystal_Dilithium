from collections import deque

def reduce_mod_pm(n, a):
    """
    Prend un entier n et le représente 
    comme un entier dans la plage

    r = n % a

    pour a impair :
        -(a-1)/2 < r <= (a-1)/2
    pour a pair :
        - a / 2  < r <= a / 2
    """
    r = n % a
    if r > (a >> 1):
        r -= a

    # assert r > -(a >> 1)
    # assert r <= (a >> 1)
    # assert (n % a) == (r % a)
    return r

def decompose(r, a, q):
    """
    Prend un élément r et le représente 
    comme :

    r = r1*a + r0

    Avec r0 dans la plage

    -(a << 1) < r0 <= (a << 1)
    """
    r  = r % q
    r0 = reduce_mod_pm(r, a)
    r1 = r - r0
    if r1 == q - 1:
        return 0, r0 - 1
    r1 = r1 // a
    assert r == r1 * a + r0
    return r1, r0

def high_bits(r, a, q):
    r1, _ = decompose(r, a, q)
    return r1
    
def low_bits(r, a, q):
    _, r0 = decompose(r, a, q)
    return r0

# def __broken_make_hint(z, r, a, q):
#     r1 = high_bits(r, a, q)
#     v1 = high_bits(r + z, a, q)        
#     return int(r1 != v1)

def make_hint(z0, r1, a, q):
    """
    La fonction ci-dessus de la documentation
    échoue parfois, mais cela semble fonctionner...

    Cela suppose que 

    TODO : découvrir quel est le cas limite pour la fonction ci-dessus
    """
    gamma2 = (a >> 1)
    if z0 <= gamma2 or z0 > (q - gamma2) or (z0 == (q - gamma2) and r1 == 0):
        return 0
    return 1

def use_hint(h, r, a, q):
    m = (q - 1) // a
    r1, r0 = decompose(r, a, q)
    if h == 1:
        if r0 > 0:
            return (r1 + 1) % m
        return (r1 - 1) % m
    return r1
    
def check_norm_bound(n, b, q):
    x = n % q                # x ∈ {0,        ...,                    ...,     q-1}
    x = ((q - 1) >> 1) - x   # x ∈ {-(q-1)/2, ...,       -1,       0, ..., (q-1)/2}
    x = x ^ (x >> 31)        # x ∈ { (q-3)/2, ...,        0,       0, ..., (q-1)/2}
    x = ((q - 1) >> 1) - x   # x ∈ {0, 1,     ...,  (q-1)/2, (q-1)/2, ...,       1}
    return x >= b

def get_n_blocks(xof, n, blocks_read):
    blocks_read += n
    # extrait les derniers n blocs
    total_bytes = 136 * blocks_read
    xof_bytes   = xof.digest(total_bytes)[-136 * n:]
    # On utilise `deque` car il a une méthode .popleft() rapide
    return deque(xof_bytes), blocks_read
    
def get_mask_integers(bit_count, xof, n, blocks_read):
    blocks_read += n * bit_count
    # extrait les derniers n * bit_count blocs
    total_bytes = 136 * blocks_read
    xof_bytes   = xof.digest(total_bytes)[-136 * n * bit_count:]
    
    r = int.from_bytes(xof_bytes, 'little')
    mask = (1 << bit_count) - 1
    mask_integers = []
    for _ in range(256):
        mask_integers.append(r & mask)
        r >>= bit_count
    return deque(mask_integers), blocks_read

def xor_bytes(a, b):
    """
    XOR deux tableaux de bytes, suppose qu'ils 
    ont la même longueur
    """
    return bytes(a ^ b for a, b in zip(a, b))
