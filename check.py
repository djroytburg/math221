from sympy import Matrix

# state_vector: returns a state vector after x iterations through a transition matrix
# Arguments:
    # tm: a transition matrix P.
    # sv: a kernel state vector, s(n).
    # x: the quantity of times to iterate the state vector through the transition matrix.
#   Output:
    # returns an array with the # of rows as sv and one column. This is the state vector s(n + x).

def state_vector (tm, sv, x):

    for _ in range(x):
        nsv = []
        for i in range(len(tm)):
            nsv.append(sv[i])
        for a in range(len(tm)):
            row = 0
            for b in range(len(tm)):
                row = row + (tm[a][b] * sv[b])
            nsv[a] = row
        sv = nsv
    return sv

# steady_state: returns the steady state vector.
# Arguments:
    # tsm: a transition matrix P.

#   Output:
    # returns a list, which is the steady state vector.

def steady_state(tsm):

    ntsm = []
    for a in range(len(tsm) - 1):
        ntsm.append(tsm[a])

    for x in range(len(tsm) - 1):
        for y in range(len(tsm)):
            if x == y:
                ntsm[x][y] = ntsm[x][y] - 1
            else:
                ntsm[x][y] = ntsm[x][y]
    a = Matrix(ntsm)
    gen_sol_unformatted = a.nullspace()

    gen_sol = []
    for i in range(len(gen_sol_unformatted[0])):
        x = float(gen_sol_unformatted[0][i])
        gen_sol.append(x)

    sum = 0.0
    for j in range(len(gen_sol)):
        sum = sum + gen_sol[j]
    t = 1.0 / sum

    output = []
    for k in range(len(gen_sol)):
        output.append(gen_sol[k] * t)
    return output

# compare_vector: compares two state vectors under the same transition matrix
# Arguments:
    # tm: a transition matrix P.
    # sv1: a kernel state vector, s(n).
    # sv2: an optional second kernel state vector, s(n). Default produces just the state_vector.
    # x: the n in question.
#   Output:
    # returns the difference between the two compared state vectors

def compare_vector (tm, sv1, x, sv2 = [-1]):
    gem1 = state_vector(tm, sv1, x)
    if sv2[0] != -1:
        gem2 = state_vector(tm, sv2, x)
    else:
        return gem1
    output = []
    for i in range(x):
        output.append(gem1[i] - gem2[i])
    return output



