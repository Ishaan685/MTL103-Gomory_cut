import numpy as np
np.set_printoptions(linewidth=200)

def solve_LP_simplex(D):

    """
    input:
    D is the full tableau
    output:
    uodated D
    """

    red_costs = D[0, 2:]
    min_arg = np.argmin(red_costs)
    min_value = red_costs[min_arg]
    chosen_arg = np.argmax(red_costs<0)
    chosen_value = red_costs[chosen_arg]

    while (min_value < 0):
        first_col = D[:, 0].copy() #var_id
        # print(first_col)
        u = D[1:, chosen_arg+2] 
        if (np.max(u) <= 0):
            print("optimal solution is unbounded")
            break

        x_B = D[1:, 1] #x_b
        indices = np.array(range(1, len(x_B)+1))

        pos_u = u[u > 0]
        pos_x_B = x_B[u > 0]
        pos_indices = indices[u > 0]

        ratios = pos_x_B / pos_u
        pos_min_ind = np.argmin(ratios)
        pivot_ind = pos_indices[pos_min_ind]
        
        pivot_col = D[:, chosen_arg+2]
        pivot_row = D[pivot_ind, :] / pivot_col[pivot_ind]
        D[pivot_ind, :] = pivot_row

        for i in range(len(D)):
            if (i == pivot_ind):
                continue
            new_row = D[i, :] - (pivot_col[i] * pivot_row)
            D[i, :] = new_row

        print(first_col)
        first_col[pivot_ind] = chosen_arg+1
        D[:, 0] = first_col
        # print(first_col)

        # print("=============================================")
        # print(D)
        red_costs = D[0, 2:]
        min_arg = np.argmin(red_costs)
        min_value = red_costs[min_arg]
        chosen_arg = np.argmax(red_costs<0)
        chosen_value = red_costs[chosen_arg]

    print("func LP solvers")
    print(D)
    return D # need to look at what to return

def init_dual_feasible_tableau(A, B, C):

    """
    input:
    A is the constrain matrix m x n
    B is the constant matrix
    C is the cost matrix associated with the objective function
    returns:
    D is the initial dual feasible tableau
    """
    C = -1*C.copy()
    mod_A = np.concatenate((A, np.identity(len(A))), axis = 1) # for slack vars
    mod_B = B
    indices = np.array(range(0, len(mod_B)))
    neg_ind = indices[B < 0]

    for i in neg_ind: # make b > 0
        mod_A[i, :] = -1 * mod_A[i, :]
        mod_B[i] = -1 * mod_B[i]

    aux_A = np.concatenate((mod_A, np.identity(len(A))), axis = 1) # for initial basis finding
    aux_B = B.reshape((len(A), 1))
    init_costs = -1*(np.sum(mod_A, axis = 0)).reshape(1,-1)
    # print(mod_A)
    # print(init_costs)
    aux_C = np.concatenate((init_costs, np.zeros((1, len(mod_A)))), axis = 1)
    # aux_C = np.concatenate((np.zeros((1, len(mod_A[0]))), np.zeros((1, len(mod_A)))), axis = 1)

    first_cell = np.array([[-1 * np.sum(B)]]) #initial cost
    top_row = np.concatenate((first_cell, aux_C), axis = 1) #top row, (n+1)
    aux_D = np.concatenate((aux_B, aux_A), axis = 1) #initial soln
    aux_D = np.concatenate((top_row, aux_D), axis = 0)

    var_ind = np.array(range(len(mod_A[0])+1, len(mod_A[0])+1+len(mod_A))) #range of artificial vars
    var_ind = var_ind.reshape((len(var_ind), 1))
    top_cell = np.array([[0]])
    first_col = np.concatenate((top_cell, var_ind), axis = 0)
    aux_D = np.concatenate((first_col, aux_D), axis = 1)

    # print(aux_D) #correct
    D = solve_LP_simplex(aux_D)

    print(D[0,1])
    if (D[0,1] > 1e-10):
        print("infeasible original solution")
        return
    elif (D[0,1] < -1*1e10):
        print("something went wrong")
        return
    # else:
    
    print("tp")

    m = len(D)-1
    n = len(D[0])-1
    # print(len(D))
    # print(D)
    for i in range(1, len(D)):
        print("hello :", i)
        ind = D[i,0]
        if (ind <= n-m):
            continue
        vec = D[2:, ind]
        valid_rows = D[i, 2:n+2-m]
        if (not np.any(valid_rows)):
            # remove that row?
            continue
        #elimination
        index_vals = np.array(range(1, len(valid_rows)+1))
        valid_index_vals = index_vals[valid_rows != 0]
        change_index = valid_index_vals[0]
        first_col = D[:, 0].copy()
        first_col[i] = change_index

        pivot_ind = np.argmax(vec) + 1
        pivot_col = D[:, change_index+1]
        pivot_row = D[pivot_ind, :] / pivot_col[pivot_ind]
        D[pivot_ind] = pivot_row

        for i in range(len(D)):
            if (i == pivot_ind):
                continue
            new_row = D[i, :] - (pivot_col[i] * pivot_row)
            D[i, :] = new_row
        D[0, :] = first_col

    print("artificial variables removed")
    # print(D)

    new_D = D[:, :len(D[0])-len(D)+1]
    print(new_D)
    c_B = np.zeros((1, m))
    C = np.concatenate((C, np.zeros(len(A))))
    cost = 0
    for i in range(1, len(D)):
        # print((D[i, 0]-1))
        c_val = C[int(D[i, 0]-1)]
        c_B[0, i-1] = c_val
        cost += D[i, 1] * c_val
    # new_D[0, 1] = -1 * cost

    mult = new_D[1:, 2:]
    top_c = [[0, -1*cost]]
    # print(C.shape)
    red_costs = C - np.dot(c_B, mult)
    # print(red_costs.shape)
    top_row = np.concatenate((top_c, red_costs), axis = 1)
    new_D[0, :] = top_row
    # print("printing new D")
    # print(new_D)
    final_D = solve_LP_simplex(new_D)
    print("func initial basis")
    print(final_D)
    return final_D

def solve_LP_dual_simplex(D):

    """
    input:
    D is the initial dual feasible tableau
    """
    """
    returns : 
    D is the final dual and primal feasible tableau
    X is the final optimal solution
    """

    D = D.copy()
    X = np.zeros(D.shape[1])

    D[0,1] = float('inf')
    min_index = np.argmin(D[:,1] >= 0) #this row is changed and the variable corresponding to it exits the basis
    min_value = D[min_index,1]
    print(min_value)
    done = True
    optimal_finite = True

    while (min_value < 0):
        done = False
        min_ind_v = -1
        min_val_v = -1e10
        for j in range(2, D.shape[1]):
            if(D[min_index][j] < 0):
                if(D[0,j]/D[min_index,j] > min_val_v):
                    min_val_v = D[0,j]/D[min_index,j]
                    min_ind_v = j #index of column with largest negative ratio, this column enters the basis

        if(min_ind_v == -1e10):
            done = True
            optimal_finite = False
            print("The optimal cost is infinite")
            break
        else:
            D[min_index] = D[min_index]/D[min_index,min_ind_v]
            for row in range(0, D.shape[0]):
                if(row != min_index):
                    D[row,1:] -= ((D[min_index]*D[row,min_ind_v])/(D[min_index,min_ind_v]))[1:] #converting into unit vector column
            D[min_index,0] = min_ind_v - 1 #changing the variable id
            # print(min_ind_v-1)

            D[0,1] = float('inf')
            min_index = np.argmin(D[:,1] >= 0) #this row is changed and the variable corresponding to it exits the basis
            min_value = D[min_index,1]
    print("solved dual simplex =========================================================================")
    print(D)

    done = True
    X = np.zeros(D.shape[1]-2)
    for row in range(1,D.shape[0]):
        var_id = D[row, 0]
        var_value = D[row, 1]
        my_ind = var_id - 1
        if(abs(var_value <= 1e-10)):
            var_value = 0
        # print(my_ind, var_value)
        X[int(my_ind)] = var_value

    print("values of variables are :",X)
    return D, X


def gomory(filename):
    #Take input from file "filename"
    with open(filename, "r") as file:
        
        n_str, m_str = file.readline().strip().split()
        n = int(n_str) #columns, no of vars
        m = int(m_str) #rows, no of constraints

        b_str = file.readline().strip().split()
        B = np.array([int(b_s) for b_s in b_str])

        c_str = file.readline().strip().split()
        C = np.array([int(c_s) for c_s in c_str])

        A = np.zeros((m, n))
        for i in range(m):
            a_temp = file.readline().strip().split()
            for j in range(n):
                A[i,j] = a_temp[j]

        print("Input read successfully")

        print(A)        
        print(B)        
        print(C)        
        D = init_dual_feasible_tableau(A, B, C)
        print(D)
        #size of D must be (m+1) x (n+2)
        #(costs, main_tableau) x (vars_index, vars_value, main_tableau)

        D, X = solve_LP_dual_simplex(D)
        not_int = True

        while(not_int):
        # for i in range(5):

            for x_ind in range(1,D.shape[0]):
                # print("heyyyyyyyyyyyyyyy",D[x_ind,1])
                if(D[x_ind,1] - np.floor(D[x_ind,1]) > 1e-10 and np.ceil(D[x_ind,1]) - D[x_ind,1] > 1e-10):
                    print("the var not int is : ",D[x_ind][0])

                    # C.append(0) #the cost vector
                    a_temp = np.zeros(D.shape[1])
                    for temp_ind in range(2, D.shape[1]):
                        a_temp[temp_ind] = -1*(D[x_ind,temp_ind] - np.floor(D[x_ind,temp_ind]))
                        if(abs(a_temp[temp_ind]) <= 1e-10):
                            a_temp[temp_ind] = 0
                        # a_temp[temp_ind] = (-1 * np.modf(D[x_ind,temp_ind])[0])

                    a_temp[0] = D.shape[1]-1 #var_id in the tableau
                    a_temp[1] = -1*(D[x_ind,1] - np.floor(D[x_ind,1])) #the value of the new var in the tableau 
                    
                    D = np.append(D, [a_temp], axis = 0) #the new constraint has been added

                    col_temp = np.zeros(D.shape[0]) #(m+1) size
                    col_temp[0] = 0 #cost associated with the new variable
                    col_temp[-1] = 1 #corresponding to the last constraint

                    D = np.append(D, np.reshape(col_temp, (-1, 1)), axis=1) #the new variable has been added

                    break
            
            print("updated d with gomory cut")
            print(D)
            D, X = solve_LP_dual_simplex(D)
            not_int = False
            for x in X:
                if(x - np.floor(x) > 1e-10 and np.ceil(x) - x > 1e-10):
                    not_int = True
        
        ans = []
        for x in X:
            ans.append(round(x,2))
        print(ans[0:n])
        return ans[0:n]
    
# gomory("input1.txt")
# print(float("inf"))

            
# Rounding till 9 digits
# Degeneracy, cycling rule in both simplex, and dual simplex
# Remove artificial unused row
# Testing       






    

