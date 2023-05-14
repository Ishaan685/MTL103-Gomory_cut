# Authors

2020CS10494 Goonjan Saha  
2020CS50497 Ishaan Govil  

---

# Introduction 

In this assignment, we implement the `Gomory-cut` algorithm for finding the solution to an integer linear programming problem which in general is known to be *NP-hard*. The `Gomory-cut` method makes use of the `Cutting-plane` technique to find a cutting plane such that the inequality introduced by the cutting plane is satisfied by all the integer points in teh feasible region of the IP but not by the currente optimal solution of the relaxed LP.

---

# Input   

The *gomory-cut* function takes as input a filename of a `.txt` file which must contain the input in the following format -   
The *first line* of input contains two space separated integers `n` and `m` - The dimensions of `c` and `b` respectively.  
The *second line* has `m` space separated integers - `b1`, `b2`,..., `bm`.  
The *third line* has `n` space separated integers - `c1`, `c2`,..., `cn`.  
The *ith line* of the next `m` lines contains `n` space separated integers `a_i1`, `a_i2`,..., `a_in`.

---

# Output

The function returns the answer as an array of `n` integers `x∗1`, `x∗2`,..., `x∗n` - where `x∗` is the solution to the given *ILP*.  

---



