import os
import sys

isinfo=sys.argv[1]
out=sys.argv[2]
with open(isinfo) as fin_is, open(out,"w") as fout_is:
    b_mean=False
    b_st=False
    rg_id=""
    mean_is=0.0
    std_is=0.0
    for line in fin_is:
        fields=line.split(":")
        if fields[0]=="Mean insert size":
            b_mean=True
            ffields=fields[1].split()
            rg_id=ffields[-1]
            continue
        elif fields[0]=="Standard deviation of insert size":
            b_st=True
            continue

        if b_mean==True:
            mean_is=float(fields[0])
            b_mean=False
        if b_st==True:
            std_is=float(fields[0])
            b_st=False

            fout_is.write(rg_id+" "+str(mean_is)+" "+str(std_is)+"\n")
