import random
from time import perf_counter
import math
import os

# Caching the specific state and criteria function pair, so that it can be reused in future without necessecity of calculating it again
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
    def __init__(self, node, vs, moves):
        # Everything resets after ant return to the colony
        self.visitedStates = vs # states that ant visited
        self.currentNode = node # node where ant is currently located
        self.bestNode = None # the best state which ant
        self.moves = moves # number of moves that ant made


def detectMove(ant, dim, heur_cache, pher, tau_0, ksi, alpha, T, beta, w):
    node = ant.currentNode
    best_state = None
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
        # Keep the ant form returning to the state it came from
        if node.parent is not None and tuple(new_state) == tuple(node.parent.state):
            continue

        hh = w * h(heur_cache, state_t) + (1 - w) * moves # the criteria function value

        # ==== Getting pheromones from every correct relation and adding to tau, which is a pheromone coefficient
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
        
        delta = hh - h_curr # criteria function difference between current value and new one
        eta = pow(2.718281828, -delta / T) # the coefficient of criteria function, depends on criteria function improvement and colony disturbance

        weight = (tau ** alpha) * (eta ** beta)
        weights.append(weight)
        candidates.append((new_state, dir, hh))

    # roulette wheel, to choose the best among candidates but without being deterministic
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
    # Local evaporation mechanism, prevents ant from doing same moves for long time
    for neigh in movesForAnt[prevZero]:
        tile2 = best_state[0][neigh]
        raw2 = neigh // dim
        col2 = neigh % dim
        if col1 > col2:
            key = (tile1,tile2,3)
            if key in pher:
                pher[key] = pher[key] * (1-ksi) + ksi * tau_0 # + ksi * tau_0 - Relaxation to the initial value
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
    N = N - (N%top)
    solution = [i for i in range(1,dim*dim)]
    solution.append(0)
    pheromonesDict = {} # 0 - up, 1 - rigth, 2 - down, 3 - left
    # Setting pheromones dictionary which keys are correct relations in sliding puzzle
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
    global_best_h = math.inf # The best criterium function value 
    stagnation = 0 # stagnation is counter that keeps track how close colony is to the increment of disturbance

    # parametrs dependent on input data
    D = -delta0/math.log(math.sqrt(p_min*p_max)) * beta
    D_min = -delta0/math.log(p_min) * beta  # Tmin1 = -Δ0 / ln(p_min) * beta
    D_max = -delta0/math.log(p_max) * beta  # Tmax1 = -Δ0 / ln(p_max) * beta
    # -delta0/math.log((p_min+p_max/2)) - in other words it is how much of allowance we want to give to bad criteria function change

    tau_max = 10 * tau_0 # maximal pheromone value
    tau_min = tau_max/50 # minimal pheromone value
    
    initial_t = tuple(InitialState)
    InitialNode = Node(InitialState, InitialState.index(0), None, w * h(heur_cache, tuple(InitialState)))
    print(manhattan_LC(InitialState))
    ants = [Ant(InitialNode, {initial_t}, 0) for _ in range(N)]
    while not solved:
        bestNodesList = []
        for _ in range(s): # s steps before ants return to the base
            for ant in ants:
                ant.currentNode = detectMove(ant, dim, heur_cache, pheromonesDict, ksi, alpha, D, beta, w)
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
        distance = 1 / (bestAnt.bestNode.heuristic + 1)
        for i in range(dim*dim):
            tile = bestState[i]
            if tile == 0:
                continue
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

        # ==== Evaporating all pheromones respectfully to evaporation coefficient ====
        for key in pheromonesDict:
            pheromonesDict[key] = max(tau_min, min(pheromonesDict[key]* (1-ro), tau_max))
        print(bestAnt.bestNode.heuristic, bestAnt.bestNode.state)

        # ==== Decreasment of disturbance coefficient ====
        if minheuristic < global_best_h:
            global_best_h = minheuristic
            D *= 0.95
            stagnation = 0
        stagnation += 1

        # ==== Increasment of disturbance coefficient ====
        if stagnation > R:
            stagnation = stagnation//2
            D *= 1.2
        D = min(D_max, max(D_min, D))

        # ==== Choosing top best nodes and then making them the starting points for ants ====
        ants = []
        bestNodesList = sorted(bestNodesList, key=lambda x: x[1])[:top]
        for node in bestNodesList:
            ants.extend([Ant(node[0], {tuple(node[0].state)}, 0) for _ in range(N//top)])

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
start10x10 = [i for i in range(1, 100)]
start10x10.append(0)
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
test10x10 = shuffle(start10x10, 10)
n = 6
movesForAnt = build_moves_for_ant(n)
manhattan_LC = inicializeCriteriumFunc(n)

#           InitialState, dim, N, top, s, R, tau_0, ro, ksi, alpha, delta0, p_min, p_max, beta, w
start = perf_counter()
res1 = AntSearch(test6x6, n, 200, 40, 100, 5, 0.2, 0.05, 0.2, 2, 10, 0.001, 0.95, 1, 1) 
end = perf_counter()
# InitialState - The beginning of ants journery
# dim - sliding puzzle dimension(more of a size e.g. 3x3, 4x4, 5x5)
# N - number of ants in colony
# top - number of best nodes left after each global iteration for beam search
# s - number of steps ants would made exploring states
# R - the threshold representing how many global iterations without heuristic improvment we would tolerate
# tau_0 - the initial pheromone impact coefficient
# ro - the global pheramone evaporation coefficient. Applied to every relation in PheramoneDict
# ksi - the local pheramone evaporation coefficient. Applied localy to one-four correct relations per move of an ant
# alpha - the power of pheromone impact
# delta0 - the worst critteria function change we want to tolerate
# p_min - the toleration of critteria function bad difference with low colony disturbance coefficient
# p_max - the toleration of critteria function bad difference with high colony disturbance coefficient
# beta - the power of critteria function impact on ant move decision
# w - the weight which regulates importance of Manhattan Distance + Linear Conflict critterium and moves at the same time


res1 = PathTrace(res1)
print(f"Number of steps: {len(res1)}")
print("Time (s): ", f"{end-start:.2f}")


# def ZeroTrace(path):
#     nodesZeroList = []
#     nodesCounter = len(path)
#     for p in path:
#         nodesZeroList.append(p.index(0))

#     return nodesZeroList, nodesCounter


# with os.scandir("In") as entries:
#     for entry in entries:
#         if entry.is_file():
#             with open(entry.path, 'r') as f:
#                 t = [int(i) for i in f.readline().split(', ')]
#                 start = perf_counter()
#                 solution = AntSearch(t, n, 200, 40, 30, 5, 0.2, 0.05, 0.2, 2, 3, 0.001, 0.95, 1, 1) 
#                 end = perf_counter()
#                 path = PathTrace(solution)
#                 nodesStats = ZeroTrace(path)
#                 nodesZeroList = nodesStats[0]
#                 nodesCounter = nodesStats[1]
#                 raw1 = f"Blank element moves: {'->'.join([str(i) for i in nodesZeroList])}\n"
#                 raw2 = f"Number of moves: {nodesCounter - 1}\n"
#                 raw3 = f"Solution search time: {end-start:.5f}\n"
#                 with open(f"{'Out'}/{entry.name[:-4]}_out.txt", "w", encoding="utf-8") as f:
#                     f.writelines([raw1, raw2, raw3])