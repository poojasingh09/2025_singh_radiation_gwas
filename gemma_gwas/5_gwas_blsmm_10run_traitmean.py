#  gwasblsmm_trait_mean.py  
#  blsmm_trait_mean.loop.sh will loop over all traits and calculate mean gamma
# calculate mean blsmm parameters per SNP from 10 blsmm runs for a trait
#modified from joe mcgirr's script at https://github.com/joemcgirr/fishfASE/blob/master/gwas/gemma_mean_parameters_dentition.py

import pandas as pd
import numpy as np
import sys, argparse


working_dir = './'
trait_name = sys.argv[1]
print(trait_name)

# start loop with first param file
j=str(1)
param = pd.read_table(working_dir+'all.SNPs.filt.dentition.imp.bslmm.dentition.'+trait_name+'.'+j+'.param.txt')
#print(param)
param['beta'+j]= param['beta']
param['gamma'+j]= param['gamma']
param['alpha'+j]= param['alpha']
param = param[['chr','ps','alpha'+j,'beta'+j,'gamma'+j]]
#print(param)

# loop through all param files

iters = np.arange(2,11)

for i in iters: 
    
    j=str(i)
    parami = pd.read_table(working_dir+'all.SNPs.filt.dentition.imp.bslmm.dentition.'+trait_name+'.'+j+'.param.txt')
    parami['beta'+j]= parami['beta']
    parami['gamma'+j]= parami['gamma']
    parami['alpha'+j]= parami['alpha']
    parami = parami[['chr','ps','alpha'+j,'beta'+j,'gamma'+j]]
    param = param.merge(parami, on = ['chr', 'ps'])
#print(param)


alphas = param.filter(regex='alpha')
alphas = alphas.assign(alpha_mean=alphas.mean(axis=1))
alphas['chr']=param['chr']
alphas['ps']=param['ps']
#print(alphas)

betas = param.filter(regex='beta')
betas = betas.assign(beta_mean=betas.mean(axis=1))
betas['chr']=param['chr']
betas['ps']=param['ps']
betas = betas[['chr','ps','beta_mean']]

gammas = param.filter(regex='gamma')
gammas = alphas.assign(gamma_mean=gammas.mean(axis=1))
gammas['chr']=param['chr']
gammas['ps']=param['ps']
gammas = gammas[['chr', 'ps','gamma_mean']]

mean_params = alphas.merge(betas, on = ['chr','ps'])
mean_params = mean_params.merge(gammas, on = ['chr','ps'])
#print(mean_params)


mean_params = mean_params[['chr', 'ps', 'alpha_mean','beta_mean','gamma_mean']]
#print(mean_params)
# output files with means

mean_params.to_csv(working_dir + 'all.SNPs.filt.dentition.imp.bslmm.dentition.'+trait_name+'.mean.params.txt',index=False, sep = "\t", header = True)

