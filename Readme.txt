This code implements the algorithm NST+HT+f-FB for recovering sparse signals. 
 
The key of this theory is the concept that the number of indices selected at each iteration should be considered in order to speed up the convergence. 

The particular index selection 

$T_{k}=s$, $T_{k}=6k$, $T_{k}=8k$, $T_{k}=10k$, $T_{k}=12k$ ，$T_{k}=14k$ and $T_{k}=k^{2}$ are contained in the folder 'NST+HT+f-FB'. You can try other 
functional functions. 

The main code 'Demo' is the performance comparison within the class of NST+HT+f-FB algorithms, where m=500,n =1000, and sparsities vary from 1 to 280.
You can use this code for other settings and real data.  

The folder ''fig1_08-02-2020-02-31'' contains experimental results produced by 'Demo'. You can try by your self.

The code is provided for research. If you encounter any problem using this code, please feel free to email：
ningninghan@szu.edu.cn. Thanks! 