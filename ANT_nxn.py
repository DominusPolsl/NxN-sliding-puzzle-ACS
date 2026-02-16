import random
from time import perf_counter
import math
import os

def h(heur_cache, state_t):
    v = heur_cache.get(state_t)
    if v is None:
        v = manhattan_LC(state_t)
        heur_cache[state_t] = v
    return v

# Count inversions: comparing elements 'a' and 'b' at indices 'i' and 'j' respectively; when a > b & i > j, it's called an inversion
def countInversions(table, n):
        counter = 0
        elements = n*n
        for k in range(elements-1):
            t1 = table[k] 
            for i in range(k+1, elements):
                t2 = table[i]
                if t1 > t2 and t2 != 0 and t1 != 0:
                    counter += 1
        return counter

# Determine the row containing zero. Row indices here are 1..4 starting from the bottom row
def detectRank(table, n):
    pos = table.index(0)    
    emptyRank = n - pos//n
    return emptyRank

def isSolvable(table, n):
    inv = countInversions(table, n)
    if n % 2 == 1:
        return inv % 2 == 0
    else:
        row_from_bottom = detectRank(table, n)
        return (inv + row_from_bottom) % 2 == 1

# Randomly shuffle the table
def shuffle(table, n):
    random.shuffle(table)
    while not isSolvable(table, n):
        random.shuffle(table)
    return table

# Heuristic
def inicializeCriteriumFunc(n):
    def manhattan_LC(table):
            h = 0
            l = n * n
            for i in range(l):
                tile = table[i]
                if tile == i + 1 or tile == 0:
                    continue
                else:
                    raw_target = (tile - 1) // n
                    column_target = (tile - 1) % n
                    raw = i // n
                    column = i % n
                    h += abs(raw - raw_target) + abs(column - column_target)
            
            # Linear Conflict
            for row in range(n):
                row_tiles = []
                for col in range(n):
                    tile = table[row * n + col]
                    if tile != 0 and (tile - 1) // n == row:
                        row_tiles.append(tile)
                
                for j in range(len(row_tiles)):
                    for k in range(j + 1, len(row_tiles)):
                        if row_tiles[j] > row_tiles[k]:
                            h += 2

            for col in range(n):
                col_tiles = []
                for row in range(n):
                    tile = table[row * n + col]
                    if tile != 0 and (tile - 1) % n == col:
                        col_tiles.append(tile)
                
                for j in range(len(col_tiles)):
                    for k in range(j + 1, len(col_tiles)):
                        if col_tiles[j] > col_tiles[k]:
                            h += 2

            return h
    return manhattan_LC

manhattan_LC = inicializeCriteriumFunc(4)

def build_moves_for_ant(n: int) -> tuple[tuple[int, ...], ...]:
    if n < 2:
        raise ValueError("n musi być >= 2")

    moves = []
    for pos in range(n * n):
        r, c = divmod(pos, n)
        neigh = []

        # góra
        if r > 0:
            neigh.append(pos - n)
        # dół
        if r < n - 1:
            neigh.append(pos + n)
        # lewo
        if c > 0:
            neigh.append(pos - 1)
        # prawo
        if c < n - 1:
            neigh.append(pos + 1)

        moves.append(tuple(neigh))

    return tuple(moves)
# ======== Ant Possible Moves ========
movesForAnt = build_moves_for_ant(4)


class Node:
    def __init__(self, state, zeroPos, parent, h):
        self.heuristic = h
        self.state = state
        self.zeroPos = zeroPos
        self.parent = parent

class Ant:
    def __init__(self, node, vs, moves):
        self.visitedStates = vs
        self.currentNode = node
        self.found = False 
        self.bestNode = None
        self.moves = moves


def detectMove(ant, dim, heur_cache, pher, ksi, alpha, T, beta, w):
    node = ant.currentNode
    best_state = None
    tau_0 = 0.1
    h_curr = node.heuristic
    visited = ant.visitedStates
    moves = ant.moves + 1
    prevZero = node.zeroPos
    parent = node
    weights = []
    candidates = []

    for dir in movesForAnt[prevZero]:
        new_state = node.state.copy()
        new_state[dir], new_state[prevZero] = new_state[prevZero], new_state[dir]
        state_t = tuple(new_state)

        if state_t in visited:
            continue
        if node.parent is not None and tuple(new_state) == tuple(node.parent.state):
            continue

        hh = w * h(heur_cache, state_t) + (1 - w) * moves
        tau = tau_0
        tile1 = new_state[prevZero]
        raw1 = prevZero // dim
        col1 = prevZero % dim
        for neigh in movesForAnt[prevZero]:
            tile2 = new_state[neigh]
            raw2 = neigh // dim
            col2 = neigh % dim
            if col1 > col2:
                tau += pher.get((tile1,tile2,3), 0)
            elif col1 < col2:
                tau += pher.get((tile1,tile2,1), 0)
            elif raw1 > raw2:
                tau += pher.get((tile1,tile2,0), 0)
            else:
                tau += pher.get((tile1,tile2,2), 0)
        
        delta = hh - h_curr
        eta = pow(2.718281828, -delta / T)

        weight = (tau ** alpha) * (eta ** beta)
        weights.append(weight)
        candidates.append((new_state, dir, hh))

    # roulette wheel
    W = sum(weights)
    r = random.uniform(0, W)
    c = 0
    for i in range(len(weights)):
        c += weights[i]
        if r <= c:
            best_state = candidates[i]
            break

    if best_state is None:
        if parent.parent is not None:
            ant.currentNode = parent.parent
            ant.moves -= 1
        return ant.currentNode

    tile1 = best_state[0][prevZero]
    raw1 = prevZero // dim
    col1 = prevZero % dim
    for neigh in movesForAnt[prevZero]:
        tile2 = best_state[0][neigh]
        raw2 = neigh // dim
        col2 = neigh % dim
        if col1 > col2:
            key = (tile1,tile2,3)
            if key in pher:
                pher[key] = pher[key] * (1-ksi) + ksi * tau_0
        elif col1 < col2:
            key = (tile1,tile2,1)
            if key in pher:
                pher[key] = pher[key] * (1-ksi) + ksi * tau_0
        elif raw1 > raw2:
            key = (tile1,tile2,0)
            if key in pher:
                pher[key] = pher[key] * (1-ksi) + ksi * tau_0
        else:
            key = (tile1,tile2,2)
            if key in pher:
                pher[key] = pher[key] * (1-ksi) + ksi * tau_0

    new_node = Node(best_state[0], best_state[1], parent, best_state[2])
    ant.moves += 1
    visited.add(tuple(best_state[0]))
    return new_node

def AntSearch(InitialState, dim, N, top, s, R, tau_0, ro, ksi, alpha, delta0, p_min, p_max, beta, w):
    solution = [i for i in range(1,dim*dim)]
    solution.append(0)
    pheromonesDict = {} # 0 - top, 1 - rigth, 2 - down, 3 - left
    # {(1,2,1): 0.1, (1,4,2): 0.1}
    for i in range(dim*dim-1):
        raw = i//dim
        col = i%dim
        if raw != 0:
            t = solution[(raw-1) * dim + col]
            pheromonesDict[(i+1, t, 0)] = tau_0
        if col != dim-1:
            r = solution[raw * dim + col + 1]
            if r != 0:
                pheromonesDict[(i+1, r, 1)] = tau_0
        if raw != dim-1:
            d = solution[(raw + 1) * dim + col]
            if d != 0:
                pheromonesDict[(i+1, d, 2)] = tau_0
        if col != 0:
            l = solution[raw * dim + col - 1]
            pheromonesDict[(i+1, l, 3)] = tau_0

    heur_cache = {}
    solved = False
    global_best_h = math.inf
    stagnation = 0

    # Parametry załeżny od danych wejściowych
    D = -delta0/math.log((p_min+p_max/2)) * beta # 2.88
    D_min = -delta0/math.log(p_min) * beta  # Tmin1 = -Δ0 / ln(p_min) * beta
    D_max = -delta0/math.log(p_max) * beta  # Tmax1 = -Δ0 / ln(p_max) * beta
    tau_max = 10 * tau_0
    tau_min = tau_max/50
    
    initial_t = tuple(InitialState)
    InitialNode = Node(InitialState, InitialState.index(0), None, h(heur_cache, tuple(InitialState)))
    print(manhattan_LC(InitialState))
    ants = [Ant(InitialNode, {initial_t}, 0) for _ in range(N)]
    while not solved:
        bestNodesList = []
        visit_count = {}
        for _ in range(s):
            for ant in ants:
                ant.currentNode = detectMove(ant, dim, heur_cache, pheromonesDict, ksi, alpha, D, beta, w)
                st = tuple(ant.currentNode.state)
                if ant.bestNode == None or ant.currentNode.heuristic < ant.bestNode.heuristic:
                    ant.bestNode = ant.currentNode
                if ant.currentNode.state == solution:
                    return ant.currentNode
                visit_count[st] = visit_count.get(st, 0) + 1
            
        bestAnt = ants[0]
        minheuristic = bestAnt.bestNode.heuristic
        for ant in ants:
            bestNodesList.append((ant.bestNode, ant.bestNode.heuristic))
                
            if ant.bestNode.heuristic < minheuristic:
                minheuristic = ant.bestNode.heuristic
                bestAnt = ant

        bestState = bestAnt.bestNode.state
        distance = 1 / (bestAnt.bestNode.heuristic + 1)
        for i in range(dim*dim):
            tile = bestState[i]
            raw = i//dim
            col = i%dim
            if raw != 0:
                t = bestState[(raw-1) * dim + col]
                key = (tile, t, 0)
                if key in pheromonesDict:
                    pheromonesDict[key] += ro*distance/(1-ro)
            if col != dim-1:
                r = bestState[raw * dim + col + 1]
                key = (tile, r, 1)
                if key in pheromonesDict:
                    pheromonesDict[key] += ro*distance/(1-ro)
            if raw != dim-1:
                d = bestState[(raw + 1) * dim + col]
                key = (tile, d, 2)
                if key in pheromonesDict:
                    pheromonesDict[key] += ro*distance/(1-ro)
            if col != 0:
                l = bestState[raw * dim + col - 1]
                key = (tile, l, 3)
                if key in pheromonesDict:
                    pheromonesDict[key] += ro*distance/(1-ro)

        for key in pheromonesDict:
            pheromonesDict[key] = max(tau_min, min(pheromonesDict[key]* (1-ro), tau_max))
        print(bestAnt.bestNode.heuristic, bestAnt.bestNode.state)
        if minheuristic < global_best_h:
            global_best_h = minheuristic
            D *= 0.95
            stagnation = 0
        stagnation += 1
        if stagnation > R:
            stagnation = stagnation//2
            D *= 1.2
        D = min(D_max, max(D_min, D))
        ants = []
        bestNodesList = sorted(bestNodesList, key=lambda x: x[1])[:top]
        for node in bestNodesList:
            ants.extend([Ant(node[0], {tuple(node[0].state)}, 0) for _ in range(N//top + 1)])

def printTrace(path):
    for st in path:
        print(st)

def PathTrace(res):
    path = []
    while res != None:
        path.insert(0, res)
        res = res.parent

    states = set()
    out = []
    for p in path:
        t = tuple(p.state)
        if t in states:
            out = out[:out.index(t) + 1] # 9 12 23 41 0 
            states = set(out)
        else:
            states.add(t)
            out.append(t)
    
    return out

start6x6 = [i for i in range(1, 36)]
start6x6.append(0)
start7x7 = [i for i in range(1, 49)]
start7x7.append(0)
start8x8 = [i for i in range(1, 64)]
start8x8.append(0)
start9x9 = [i for i in range(1, 81)]
start9x9.append(0)
test80 = [0,12,9,13,15,11,10,14,7,8,5,6,4,3,2,1]
test80_2 = [0,12,9,13,15,11,10,14,3,7,2,5,4,8,6,1]
bad_conf = [1,5,9,13,2,6,10,14,3,7,11,15,4,8,12,0]
testState = [1,2,3,4,5,6,7,8,10,11,12,15,9,13,0,14]
testState56 = [12,2,6,13,1,8,15,11,0,9,14,4,5,3,10,7]
testState51 = [11,2,5,6,14,10,3,1,13,0,9,15,7,8,4,12]
testhz = [1, 2, 3, 4, 5, 10, 7, 8, 9, 13, 11, 12, 0, 14, 15, 6]
test5x5 = [2,17,1,5,23,15,10,7,8,4,21,20,19,0,24,3,11,22,9,12,18,13,16,6,14]

test6x6 = shuffle(start6x6, 6)
test7x7 = shuffle(start7x7, 7)
test8x8 = shuffle(start8x8, 8)
test9x9 = shuffle(start9x9, 9)
n = 5
movesForAnt = build_moves_for_ant(n)
manhattan_LC = inicializeCriteriumFunc(n)

#           InitialState, dim, N, top, s, R, tau_0, ro, ksi, alpha, delta0, p_min, p_max, beta, w
res1 = AntSearch(test5x5, n, 2000, 200, 50, 10, 0.2, 0.05, 0.2, 2, 3, 0.05, 0.7, 1, 0.8) 
# movesForAnt = build_moves_for_ant(9)
# manhattan_LC = inicializeCriteriumFunc(9)
# res2 = AntSearch(test9x9, 9, 2000, 200, 50, 10, 0.2, 0.2, 0.3, 1, 4, 0.05, 0.7, 1, 1)
# delta0 - jaką złą zmianę heurystyki chciałbym dopuścić
# p_min - jak dużo chcę dopuszczać zmianę heurystyki delta0 z małym niepokojem kolonii
# p_max - jak dużo chcę dopuszczać zmianę heurystyki delta0 z dużym niepokojem kolonii

res1 = PathTrace(res1)
print(len(res1))
# res2 = PathTrace(res2)
# print(len(res2))
