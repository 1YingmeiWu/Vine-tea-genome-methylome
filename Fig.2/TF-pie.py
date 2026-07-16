# -*- coding: utf-8 -*-
"""
v1.0.0.202205    jjl@j-jl.com
"""
#%%---------------------------------------------------------------
import pandas as pd
import os,sys,argparse,re,glob
#from loguru import logger
#%%-----------------------------------------------------------------
Bin,filename = os.path.split(os.path.realpath(__file__))
parser = argparse.ArgumentParser(description='test\n')
parser.add_argument('-i', '--indata', dest='indata', required=True, help='配置文件') # required=True , 
parser.add_argument('-l', '--label',  dest='label',  required=True, help='标签列列名') # required=True , 
parser.add_argument('-n', '--number', dest='number', default="-",   help='数字列列名，可以为-，此时会按照标签列频数计算') # required=True , 
parser.add_argument('-t', '--title',  dest='title',  default="-", help='标题,默认为输出文件名') # required=True , 
parser.add_argument('-o', '--prefix', dest='prefix', default=f"{os.getcwd()}/pie",  help='输出文件前缀')
args   = parser.parse_args()
if args.title == "-": args.title = os.path.basename(os.path.abspath(args.prefix))
#%%---------------------------------------------------------------------------
def main():

    if not os.path.exists(os.path.abspath(os.path.dirname(args.prefix))): os.makedirs(os.path.abspath(os.path.dirname(args.prefix)))
    #---------------------------------------------------------------------------
    df_data = pd.read_csv(args.indata,sep="\t")
    info_dict = {}
    if args.number == "-":
        data = pd.DataFrame(df_data[args.label].value_counts())
        if data.shape[0] > 20 : data = data.head(20)
        for i,row in data.iterrows(): 
            name = i
            info_dict[name] = {"name": name, "value": row[args.label]}
    else:
        if df_data.shape[0] > 20 : df_data = df_data.head(20)
        for i,row in df_data.iterrows(): 
            name =  row[args.label]
            info_dict[name] = {"name": name, "value": row[args.number]}

    #---------------------------------------------------------------------------
    w = open(f"{args.prefix}.html","w")
    with open(f"{Bin}/pie.echart.html","r") as f:
        for line in f.readlines():
            if re.search("_replace",line):
                line = line.replace("data_replace",str([value for key,value in info_dict.items()]))
                line = line.replace("legend_replace",str(list(info_dict.keys())))
                line = line.replace("title_replace",args.title)
            w.write(line)
    w.close()
    #logger.info(f"outfile: {args.prefix}.html")
    #---------------------------------------------------------------------------
    #---------------------------------------------------------------------------
    #---------------------------------------------------------------------------
    #---------------------------------------------------------------------------
#%%-----------------------------------------------------------------------------
if __name__ == "__main__":
    #logger.success("开始运行")
    main()
    #logger.success('运行完成' )