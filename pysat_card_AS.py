'''
@Descripttion: CNF Encoding of Boolean Cardinality Constraints(Xoodoo_AS)
@Author: Huina Li
@Date: 2021-04-27
@LastEditTime:
'''
# from pysat.card import *

# cnf = CardEnc.atmost(lits=[1, 2, 3], bound=2, encoding=EncType.seqcounter)
# print (cnf.clauses)

from pysat.card import *
import argparse

parser = argparse.ArgumentParser(description='CNF Encoding of Boolean Cardinality Constraints(Xoodoo_AS)')
parser.add_argument('--bound_num', '-b', default=25, type=int, help='This argument is max AS node number. default is 25')
parser.add_argument('--var_num', '-v', default=384, type=int, nargs='+', help='This argument is variable number list. default is 384,\
                                                                            meaning [1,2,...,384]')
parser.add_argument('--mode', '-m', default=0, type=int, choices=[0,1,2], help='This argument is card sum mode, 0 for atmost, 1 for atleast,\
                                                                            2 for equals. default is 0')
args = parser.parse_args()
bound_num = args.bound_num
var_num = args.var_num

lits = []
total = ''
for var in var_num:
    lits += [i for i in range(1,var+1)]
    total += str(var) + '_'

cnf = []
filepath = ['<=', '>=', '=']
filepath = 'CNF_'+total+'AS'+filepath[args.mode]+str(bound_num)+'.txt'
if args.mode == 0:
    cnf = CardEnc.atmost(lits, bound=bound_num, encoding=EncType.seqcounter)
elif args.mode == 1:
    cnf = CardEnc.atleast(lits, bound=bound_num, encoding=EncType.seqcounter)
elif args.mode == 2:
    cnf = CardEnc.equals(lits, bound=bound_num, encoding=EncType.seqcounter)

print(cnf.nv)
with open(filepath,'w+') as f:
    for i in cnf.clauses:
        temp = []
        for ii in i:
            ii = str(ii)
            temp.append(ii)
        f.write(' '.join(temp))
        f.write('\n')