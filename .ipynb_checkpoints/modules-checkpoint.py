class Module:
    def __init__(self, ring):
        self.ring = ring
        
    def __bit_unpack(self, input_bytes, m, n, alg, packed_len, *args, is_ntt=False):
        # Divise les octets d'entrée en morceaux de longueur fixed (packed_len)
        poly_bytes = [input_bytes[i:i+packed_len] for i in range(0, len(input_bytes), packed_len)]
        # Applique l'algorithme de décompression à chaque morceau pour créer une matrice de polynômes
        matrix = [[alg(poly_bytes[n*i+j], *args, is_ntt=is_ntt)
                   for j in range(n)]
                   for i in range(m)]
        return self(matrix)
        
    def bit_unpack_t0(self, input_bytes, m, n, is_ntt=False):
        # Décompresse les octets d'entrée en utilisant l'algorithme bit_unpack_t0
        packed_len = 416
        algorithm = self.ring.bit_unpack_t0
        return self.__bit_unpack(input_bytes, m, n, algorithm, packed_len, is_ntt=False)
    
    def bit_unpack_t1(self, input_bytes, m, n, is_ntt=False):
        # Décompresse les octets d'entrée en utilisant l'algorithme bit_unpack_t1
        packed_len = 320
        algorithm = self.ring.bit_unpack_t1
        return self.__bit_unpack(input_bytes, m, n, algorithm, packed_len, is_ntt=False)
        
    def bit_unpack_s(self, input_bytes, m, n, eta, is_ntt=False):
        # Décompresse les octets d'entrée en utilisant l'algorithme bit_unpack_s avec un paramètre eta
        if eta == 2:
            packed_len = 96
        elif eta == 4:
            packed_len = 128
        else:
            raise ValueError("eta doit être soit 2 soit 4")
        algorithm = self.ring.bit_unpack_s
        return self.__bit_unpack(input_bytes, m, n, algorithm, packed_len, eta, is_ntt=False)
        
    def bit_unpack_w(self, input_bytes, m, n, gamma_2, is_ntt=False):
        # Décompresse les octets d'entrée en utilisant l'algorithme bit_unpack_w avec un paramètre gamma_2
        if gamma_2 == 95232:
            packed_len = 192
        elif gamma_2 == 261888:
            packed_len = 128
        else:
            raise ValueError("gamma_2 doit être soit (q-1)/88 soit (q-1)/32")
        algorithm = self.ring.bit_unpack_w
        return self.__bit_unpack(input_bytes, m, n, algorithm, packed_len, gamma_2, is_ntt=False)
    
    def bit_unpack_z(self, input_bytes, m, n, gamma_1, is_ntt=False):
        # Décompresse les octets d'entrée en utilisant l'algorithme bit_unpack_z avec un paramètre gamma_1
        if gamma_1 == (1 << 17):
            packed_len = 576
        elif gamma_1 == (1 << 19):
            packed_len = 640
        else:
            raise ValueError("gamma_1 doit être soit 2^17 soit 2^19")
        algorithm = self.ring.bit_unpack_z
        return self.__bit_unpack(input_bytes, m, n, algorithm, packed_len, gamma_1, is_ntt=False)
        
    def __repr__(self):
        return f"Module sur l'anneau commutatif : {self.ring}"

    def __str__(self):
        return f"Module sur l'anneau commutatif : {self.ring}"
        
    def __eq__(self, other):
        return self.ring == other.ring

    def __call__(self, matrix_elements):
        # Vérifie que les éléments fournis sont une liste de listes de polynômes
        if not isinstance(matrix_elements, list):
            raise TypeError("Les éléments d'un module sont des matrices.")
        if isinstance(matrix_elements[0], list):
            for element_list in matrix_elements:
                if not all(isinstance(aij, self.ring.element) for aij in element_list):
                    raise TypeError(f"Tous les éléments de la matrice doivent être des éléments de l'anneau : {self.ring}")
            return Module.Matrix(self, matrix_elements)
        elif isinstance(matrix_elements[0], self.ring.element):
            if not all(isinstance(aij, self.ring.element) for aij in matrix_elements):
                raise TypeError(f"Tous les éléments de la matrice doivent être des éléments de l'anneau : {self.ring}")
            return Module.Matrix(self, [matrix_elements])
        else:
            raise TypeError("Les éléments d'un module sont des matrices, construites à partir des éléments de l'anneau de base.")

    class Matrix:
        def __init__(self, parent, matrix_elements):
            self.parent = parent
            self.rows = matrix_elements
            self.m = len(matrix_elements)
            self.n = len(matrix_elements[0])
            if not self.check_dimensions():
                raise ValueError("Longueurs des lignes incohérentes dans la matrice")

        def get_dim(self):
            # Retourne les dimensions de la matrice
            return self.m, self.n

        def check_dimensions(self):
            # Vérifie que toutes les lignes de la matrice ont la même longueur
            return all(len(row) == self.n for row in self.rows)

        def transpose(self):
            # Retourne la transposée de la matrice
            new_rows = [list(item) for item in zip(*self.rows)]
            return self.parent(new_rows)

        def transpose_self(self):
            # Transpose la matrice en place
            self.m, self.n = self.n, self.m
            self.rows = [list(item) for item in zip(*self.rows)]
            return self
            
        def reduce_coefficents(self):
            # Réduit tous les coefficients des polynômes de la matrice modulo q
            for row in self.rows:
                for ele in row:
                    ele.reduce_coefficents()
            return self
            
        def to_montgomery(self):
            # Convertit tous les éléments de la matrice en forme de Montgomery
            for row in self.rows:
                for ele in row:
                    ele.to_montgomery()
            return self
            
        def from_montgomery(self):
            # Convertit tous les éléments de la matrice de la forme de Montgomery à la forme standard
            for row in self.rows:
                for ele in row:
                    ele.from_montgomery()
            return self
            
        def scale(self, other):
            # Multiplie chaque élément de la matrice par un polynôme ou un entier
            if not (isinstance(other, self.parent.ring.Polynomial) 
                    or isinstance(other, int)):
                raise TypeError("La multiplication des éléments peut seulement se faire avec des polynômes ou des entiers")
            matrix = [[other * ele for ele in row] for row in self.rows]
            return self.parent(matrix)
    
        def check_norm_bound(self, bound):
            # Vérifie si la norme de chaque polynôme de la matrice respecte une limite donnée
            for row in self.rows:
                if any(p.check_norm_bound(bound) for p in row):
                    return True
            return False

        def power_2_round(self, d):
            # Applique l'opération `power_2_round` à chaque élément de la matrice pour créer deux matrices
            m1_elements = [[0 for _ in range(self.n)] for _ in range(self.m)]
            m0_elements = [[0 for _ in range(self.n)] for _ in range(self.m)]
            for i in range(self.m):
                for j in range(self.n):
                    m1_ele, m0_ele = self[i][j].power_2_round(d)
                    m1_elements[i][j] = m1_ele
                    m0_elements[i][j] = m0_ele
            return self.parent(m1_elements), self.parent(m0_elements)
            
        def decompose(self, alpha):
            # Applique l'opération `decompose` à chaque élément de la matrice pour créer deux matrices
            m1_elements = [[0 for _ in range(self.n)] for _ in range(self.m)]
            m0_elements = [[0 for _ in range(self.n)] for _ in range(self.m)]
            for i in range(self.m):
                for j in range(self.n):
                    m1_ele, m0_ele = self[i][j].decompose(alpha)
                    m1_elements[i][j] = m1_ele
                    m0_elements[i][j] = m0_ele
            return self.parent(m1_elements), self.parent(m0_elements)

        def __bit_pack(self, algorithm, *args):
            # Combine les éléments compressés en une séquence d'octets
            return b"".join(algorithm(poly, *args) for row in self.rows for poly in row)
            
        def bit_pack_t1(self):
            # Compresse les éléments de la matrice en utilisant l'algorithme bit_pack_t1
            algorithm = self.parent.ring.Polynomial.bit_pack_t1
            return self.__bit_pack(algorithm)
            
        def bit_pack_t0(self):
            # Compresse les éléments de la matrice en utilisant l'algorithme bit_pack_t0
            algorithm = self.parent.ring.Polynomial.bit_pack_t0
            return self.__bit_pack(algorithm)
            
        def bit_pack_s(self, eta):
            # Compresse les éléments de la matrice en utilisant l'algorithme bit_pack_s avec un paramètre eta
            algorithm = self.parent.ring.Polynomial.bit_pack_s
            return self.__bit_pack(algorithm, eta)
        
        def bit_pack_w(self, gamma_2):
            # Compresse les éléments de la matrice en utilisant l'algorithme bit_pack_w avec un paramètre gamma_2
            algorithm = self.parent.ring.Polynomial.bit_pack_w
            return self.__bit_pack(algorithm, gamma_2)
        
        def bit_pack_z(self, gamma_1):
            # Compresse les éléments de la matrice en utilisant l'algorithme bit_pack_z avec un paramètre gamma_1
            algorithm = self.parent.ring.Polynomial.bit_pack_z
            return self.__bit_pack(algorithm, gamma_1)
            
        def to_ntt(self):
            # Convertit tous les éléments de la matrice en forme NTT
            for row in self.rows:
                for ele in row:
                    ele.to_ntt()
            return self
    
        def from_ntt(self):
            # Convertit tous les éléments de la matrice de la forme NTT à la forme standard
            for row in self.rows:
                for ele in row:
                    ele.from_ntt()
            return self
            
        def copy_to_ntt(self):
            # Crée une copie de la matrice avec tous les éléments en forme NTT
            matrix = [[ele.copy_to_ntt() for ele in row] for row in self.rows]
            return self.parent(matrix)
            
        def copy_from_ntt(self):
            # Crée une copie de la matrice avec tous les éléments en forme standard
            matrix = [[ele.copy_from_ntt() for ele in row] for row in self.rows]
            return self.parent(matrix)
        
        def high_bits(self, alpha, is_ntt=False):
            # Obtient les bits de poids fort des éléments de la matrice
            matrix = [[ele.high_bits(alpha, is_ntt=is_ntt) for ele in row] for row in self.rows]
            return self.parent(matrix)
            
        def low_bits(self, alpha, is_ntt=False):
            # Obtient les bits de poids faible des éléments de la matrice
            matrix = [[ele.low_bits(alpha, is_ntt=is_ntt) for ele in row] for row in self.rows]
            return self.parent(matrix)            
            
        def __getitem__(self, i):
            return self.rows[i]

        def __eq__(self, other):
            return other.rows == self.rows

        def __add__(self, other):
            # Addition de deux matrices
            if not isinstance(other, Module.Matrix):
                raise TypeError("Les matrices peuvent seulement être ajoutées à d'autres matrices")
            if self.parent != other.parent:
                raise TypeError("Les matrices doivent avoir le même anneau de base")
            if self.get_dim() != other.get_dim():
                raise ValueError("Les matrices n'ont pas les mêmes dimensions")
            new_elements = []
            for i in range(self.m):
                new_elements.append([a+b for a,b in zip(self.rows[i], other.rows[i])])
            return self.parent(new_elements)

        def __radd__(self, other):
            return self.__add__(other)

        def __iadd__(self, other):
            self = self + other
            return self

        def __sub__(self, other):
            # Soustraction de deux matrices
            if not isinstance(other, Module.Matrix):
                raise TypeError("Les matrices peuvent seulement être soustraites à d'autres matrices")
            if self.parent != other.parent:
                raise TypeError("Les matrices doivent avoir le même anneau de base")
            if self.get_dim() != other.get_dim():
                raise ValueError("Les matrices n'ont pas les mêmes dimensions")
            new_elements = []
            for i in range(self.m):
                new_elements.append([a-b for a,b in zip(self.rows[i], other.rows[i])])
            return self.parent(new_elements)

        def __rsub__(self, other):
            return self.__sub__(other)

        def __isub__(self, other):
            self = self - other
            return self

        def __matmul__(self, other):
            # Produit matriciel (A @ B)
            if not isinstance(other, Module.Matrix):
                raise TypeError("Les matrices peuvent seulement être multipliées par d'autres matrices")
            if self.parent != other.parent:
                raise TypeError("Les matrices doivent avoir le même anneau de base")
            if self.n != other.m:
                raise ValueError("Les matrices ont des dimensions incompatibles")
            new_elements = [[sum(a*b for a,b in zip(A_row, B_col)) for B_col in other.transpose().rows] for A_row in self.rows]
            return self.parent(new_elements)

        def __repr__(self):
            if len(self.rows) == 1:
                return str(self.rows[0])
            max_col_width = []
            for n_col in range(self.n):
                max_col_width.append(max(len(str(row[n_col])) for row in self.rows))
            info = ']\n['.join([', '.join([f'{str(x):>{max_col_width[i]}}' for i,x in enumerate(r)]) for r in self.rows])
            return f"[{info}]"
