#!/usr/bin/env python3
"""Rebuild the 15 saved multistart records; --run computes any missing fits."""
from paired_benchmark import *

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run',action='store_true')
    args=parser.parse_args();rows=[]
    for rep in range(1,4):
        for start in range(1,6):
            cid,data=create_case(rep,50,start)
            if args.run:run_python(cid,data,rep,50,'pot_exact',0)
            row=evaluate(cid,data,'pot_exact');row['start']=start;rows.append(row)
    keys=sorted(set().union(*(row.keys() for row in rows)))
    write_rows(BASE/'multistart_results.csv',[{k:row.get(k,'') for k in keys} for row in rows])
    print('Rebuilt 15 multistart records.')
