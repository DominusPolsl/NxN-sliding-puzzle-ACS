import random
from time import perf_counter
import math
import os

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

# Function that checks the solvability of sliding puzzle. When the side length of sliding puzzle is odd then
# only the parity of inversions are taken into an account. When the side length is even then the parity of number of inversions
# and the zero raw rank are taken into account
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

        # up
        if r > 0:
            neigh.append(pos - n)
        # down
        if r < n - 1:
            neigh.append(pos + n)
        # left
        if c > 0:
            neigh.append(pos - 1)
        # right
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
    def __init__(self, node, vs):
        # Everything resets after ant return to the colony
        self.visitedStates = vs # states that ant visited
        self.currentNode = node # node where ant is currently located
        self.bestNode = None # the best state which ant

def detectMove(ant, n, pher, tau_0, xi, alpha, D, beta):
    node = ant.currentNode
    best_state = None
    h_curr = node.heuristic
    visited = ant.visitedStates
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
        # Keep the ant form returning to the state it came from
        if node.parent is not None and tuple(new_state) == tuple(node.parent.state):
            continue

        hh = manhattan_LC(state_t) # the criteria function value

        # ==== Getting pheromones from every correct relation and adding to tau, which is a pheromone coefficient
        tile1 = new_state[prevZero]
        raw1 = prevZero // n
        col1 = prevZero % n
        pher_sum = 0
        for neigh in movesForAnt[prevZero]:
            tile2 = new_state[neigh]
            raw2 = neigh // n
            col2 = neigh % n
            if col1 > col2:
                pher_sum += pher.get((tile1,tile2,3), 0)
            elif col1 < col2:
                pher_sum += pher.get((tile1,tile2,1), 0)
            elif raw1 > raw2:
                pher_sum += pher.get((tile1,tile2,0), 0)
            else:
                pher_sum += pher.get((tile1,tile2,2), 0)
        # tau = tau_0 + pher_sum/len(movesForAnt[prevZero])
        tau = tau_0 + pher_sum
        
        delta = hh - h_curr # criteria function difference between new value and current
        if delta <= 0:
            eta = math.exp(-delta) # the coefficient of criteria function, depends on criteria function improvement and colony disturbance
        else:
            eta = math.exp(-math.log(delta) / math.log(D))

        weight = (tau ** alpha) * (eta ** beta)
        weights.append(weight)
        candidates.append((new_state, dir, hh))

    if len(candidates) == 0:
        if parent.parent is not None:
            ant.currentNode = parent.parent
        return ant.currentNode

    # roulette wheel, to choose the best among candidates but without being deterministic
    W = sum(weights)
    r = random.uniform(0, W)
    c = 0
    for i in range(len(weights)):
        c += weights[i]
        if r <= c:
            best_state = candidates[i]
            break

    tile1 = best_state[0][prevZero]
    raw1 = prevZero // n
    col1 = prevZero % n
    # Local evaporation mechanism, prevents ant from doing same moves for long time
    for neigh in movesForAnt[prevZero]:
        tile2 = best_state[0][neigh]
        raw2 = neigh // n
        col2 = neigh % n
        if col1 > col2:
            key = (tile1,tile2,3)
            if key in pher:
                pher[key] = pher[key] * (1-xi) + xi * tau_0 # + xi * tau_0 - Relaxation to the initial value
        elif col1 < col2:
            key = (tile1,tile2,1)
            if key in pher:
                pher[key] = pher[key] * (1-xi) + xi * tau_0
        elif raw1 > raw2:
            key = (tile1,tile2,0)
            if key in pher:
                pher[key] = pher[key] * (1-xi) + xi * tau_0
        else:
            key = (tile1,tile2,2)
            if key in pher:
                pher[key] = pher[key] * (1-xi) + xi * tau_0

    new_node = Node(best_state[0], best_state[1], parent, best_state[2])
    visited.add(tuple(best_state[0]))
    return new_node

def AntSearch(initialState, n, N, top, s_min, s_max, tau_0, rho, xi, alpha, R, a_min, a_max, Dc, Dp, beta):
    
    N = N - (N%top)
    s = s_min
    scale = n * n
    solution = [i for i in range(1,n*n)]
    solution.append(0)
    pheromonesDict = {} # 0 - up, 1 - rigth, 2 - down, 3 - left
    # Setting pheromones dictionary which keys are correct relations in sliding puzzle
    # {(1,2,1): 0.1, (1,4,2): 0.1}
    for i in range(n*n-1):
        raw = i//n
        col = i%n
        if raw != 0:
            t = solution[(raw-1) * n + col]
            pheromonesDict[(i+1, t, 0)] = tau_0
        if col != n-1:
            r = solution[raw * n + col + 1]
            if r != 0:
                pheromonesDict[(i+1, r, 1)] = tau_0
        if raw != n-1:
            d = solution[(raw + 1) * n + col]
            if d != 0:
                pheromonesDict[(i+1, d, 2)] = tau_0
        if col != 0:
            l = solution[raw * n + col - 1]
            pheromonesDict[(i+1, l, 3)] = tau_0

    global_best_h = math.inf # The best criterium function value 
    stagnation = 0 # stagnation is counter that keeps track how close colony is to the increment of disturbance

    # parametrs dependent on input data
    delta0 = 2*(n - 1) - 1
    D_min = (math.exp(-math.log(1+delta0, a_min))) * beta
    D_max = (math.exp(-math.log(1+delta0, a_max))) * beta 
    D = (math.sqrt(D_min*D_max)) 

    
    tau_min = tau_0# minimal pheromone value
    
    initial_t = tuple(initialState)
    InitialNode = Node(initialState, initialState.index(0), None, manhattan_LC(tuple(initialState)))
    print(manhattan_LC(initialState))
    ants = [Ant(InitialNode, {initial_t}) for _ in range(N)]
    while True:
        bestNodesList = []
        for _ in range(s): # s steps before ants return to the base
            for ant in ants:
                ant.currentNode = detectMove(ant, n, pheromonesDict, tau_0, xi, alpha, D, beta)
                if ant.bestNode == None or ant.currentNode.heuristic < ant.bestNode.heuristic:
                    ant.bestNode = ant.currentNode
                if ant.currentNode.state == solution:
                    return ant.currentNode

        # ==== Choosing the best node among all ants ==== 
        bestAnt = ants[0]
        minheuristic = bestAnt.bestNode.heuristic
        for ant in ants:
            bestNodesList.append((ant.bestNode, ant.bestNode.heuristic))
                
            if ant.bestNode.heuristic < minheuristic:
                minheuristic = ant.bestNode.heuristic
                bestAnt = ant

        bestState = bestAnt.bestNode.state

        # ==== Pheromone increment for correct relations in the node with the best criteria function value ====
        distance = 1 / (1 + minheuristic/scale)
        for i in range(n*n):
            tile = bestState[i]
            if tile == 0:
                continue
            raw = i//n
            col = i%n
            if raw != 0:
                t = bestState[(raw-1) * n + col]
                key = (tile, t, 0)
                if key in pheromonesDict:
                    pheromonesDict[key] += distance/(1-rho)
            if col != n-1:
                r = bestState[raw * n + col + 1]
                key = (tile, r, 1)
                if key in pheromonesDict:
                    pheromonesDict[key] += distance/(1-rho)
            if raw != n-1:
                d = bestState[(raw + 1) * n + col]
                key = (tile, d, 2)
                if key in pheromonesDict:
                    pheromonesDict[key] += distance/(1-rho)
            if col != 0:
                l = bestState[raw * n + col - 1]
                key = (tile, l, 3)
                if key in pheromonesDict:
                    pheromonesDict[key] += distance/(1-rho)

        # ==== maximum pheromone value level adaptaion to the solving process progress ====
        tau_max = 1 / (1 + global_best_h / scale)

        # ==== Evaporating all pheromones respectfully to evaporation coefficient ====
        for key in pheromonesDict:
            pheromonesDict[key] = max(tau_min, min(pheromonesDict[key] * (1-rho), tau_max))

        # ==== Stagnation and disturbance actualization ====
        if minheuristic < global_best_h:
            global_best_h = minheuristic
            stagnation = 0
            s = s_min
            s = int(round(s))
            D = D_min + (1-Dc) * (D - D_min) 
        else:
            stagnation += 1
            q = min(1.0, stagnation / R)
            s = int(math.floor(s_min + (s_max - s_min) * (q ** 2)))

            D_target = D_min + (D_max - D_min) * (q ** 2)
            D = (1-Dp) * D + Dp * D_target   # wygładzanie
        D = min(D_max, max(D_min, D))
        
        print(bestAnt.bestNode.heuristic, bestAnt.bestNode.state)

        # ==== Choosing top best nodes and then making them the starting points for ants ====
        ants = []
        bestNodesList = sorted(bestNodesList, key=lambda x: x[1])[:top]
        for node in bestNodesList:
            ants.extend([Ant(node[0], {tuple(node[0].state)}) for _ in range(N//top)])

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

start5x5 = [i for i in range(1, 25)]
start5x5.append(0)
start6x6 = [i for i in range(1, 36)]
start6x6.append(0)
start7x7 = [i for i in range(1, 49)]
start7x7.append(0)
start8x8 = [i for i in range(1, 64)]
start8x8.append(0)
start9x9 = [i for i in range(1, 81)]
start9x9.append(0)
start10x10 = [i for i in range(1, 100)]
start10x10.append(0)
start11x11 = [i for i in range(1, 121)]
start11x11.append(0)
test80 = [0,12,9,13,15,11,10,14,7,8,5,6,4,3,2,1]
test80_2 = [0,12,9,13,15,11,10,14,3,7,2,5,4,8,6,1]
bad_conf = [1,5,9,13,2,6,10,14,3,7,11,15,4,8,12,0]
testState = [1,2,3,4,5,6,7,8,10,11,12,15,9,13,0,14]
testState56 = [12,2,6,13,1,8,15,11,0,9,14,4,5,3,10,7]
testState51 = [11,2,5,6,14,10,3,1,13,0,9,15,7,8,4,12]
testhz = [1, 2, 3, 4, 5, 10, 7, 8, 9, 13, 11, 12, 0, 14, 15, 6]
test5x5 = shuffle(start5x5, 5)
# test5x5 = [2,17,1,5,23,15,10,7,8,4,21,20,19,0,24,3,11,22,9,12,18,13,16,6,14]

test6x6 = [10, 31, 27, 2, 19, 16, 21, 15, 28, 22, 20, 9, 30, 29, 1, 5, 34, 26, 33, 14, 25, 24, 17, 4, 13, 32, 35, 0, 7, 23, 3, 18, 6, 8, 12, 11]
test7x7 = shuffle(start7x7, 7)
test7x7 = [13, 25, 20, 29, 26, 23, 43, 10, 8, 5, 44, 32, 15, 22, 42, 34, 28, 30, 3, 27, 6, 45, 19, 11, 46, 47, 14, 17, 18, 9, 24, 37, 4, 33, 21, 31, 1, 0, 2, 40, 39, 12, 48, 36, 16, 38, 41, 7, 35]
test8x8 = shuffle(start8x8, 8)
# test8x8 = [14, 6, 0, 26, 11, 37, 34, 19, 36, 47, 60, 2, 55, 63, 56, 25, 24, 16, 29, 15, 39, 42, 49, 48, 22, 8, 45, 5, 38, 33, 28, 40, 52, 57, 44, 9, 18, 1, 10, 31, 13, 43, 41, 51, 59, 32, 46, 21, 54, 35, 12, 7, 62, 3, 27, 61, 58, 53, 23, 20, 50, 4, 17, 30]
test9x9 = shuffle(start9x9, 9)
test10x10 = shuffle(start10x10, 10)
test11x11 = shuffle(start11x11, 11)


n = 8
movesForAnt = build_moves_for_ant(n)
manhattan_LC = inicializeCriteriumFunc(n)

#           initialState, n, N, top, s_min, s_max, tau_0, rho, xi, alpha, R, a_min, a_max, Dc, Dp, beta
start = perf_counter()
res1 = AntSearch(test8x8, n, 100, 50, 50, 90, 0.1, 0.2, 0.4, 1, 10, 0.1, 0.85, 0.6, 0.05, 3) 
end = perf_counter()
# initialState - The beginning of ants journery
# n - sliding puzzle dimension(more of a size e.g. 3x3, 4x4, 5x5)
# N - number of ants in colony
# top - number of best nodes left after each global iteration for beam search
# s_min - minimum number of steps that ant make exploring states in one global iteration
# s_max - maximum number of steps that wil make exploring states in one global iteration
# tau_0 - the initial pheromone impact coefficient
# rho - the global pheramone evaporation coefficient. Applied to every relation in PheramoneDict
# xi - the local pheramone evaporation coefficient. Applied localy to one-four correct relations per move of an ant
# R - A stagnation scaling factor that regulates how heavily non-improving iterations impact the increase in ant exploration steps and the overall disturbance of the colony.
# alpha - the power of pheromone impact
# a_min - how much of allowance for ant during low disturbance
# a_max - allowance, parametr which determine how much of exploration we want to allow during high distrubance
# Dc - how fast ants will become less disturbed in case of criteria function improvment (0-1)
# Dp - how fast ants will become more disturbed in case of criteria function stagnation (0-1)
# beta - the power of critteria function impact on ant move decision

res1 = PathTrace(res1)
print(f"Number of steps: {len(res1)}")
print("Time (s): ", f"{end-start:.3f}")
