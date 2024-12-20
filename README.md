# SLtest

This is a repo for testing the simplified likelihood using third moments, as outlined in the paper https://arxiv.org/abs/1809.05548

The test is based on the datacards of the [EFT tutorial](https://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/latest/model_building_tutorial2024/model_building_exercise/?h=eft#eft-model)

The idea is to have the fit done with combine, and then compare the likelihood scan with the simplified likelihood approach. The logical flow is as follows:

1. Run the fit with combine
2. Taking N replicas of the dataset using the bootstrapping technique
3. Computing the three static moments of the dataset
4. Computing the likelihood scan with the simplified likelihood approach

## Theory
The starting point is to write the test statistic $-2\Delta\text{NLL}$ as a N-dimensional normal distribution:

$$
-2\Delta\text{NLL} = \mathbf{x}^T \rho^{-1} \mathbf{x}
$$

where $\mathbf{x}$ is a gaussian variable $\mathbf{x}=\mathcal{N}(\mathbf{0}, \rho)$, and $\rho$ is the covariance matrix of the elements $\rho_{ij}=\text{Cov}[x_i,x_j]$.

Suppose we have our measurement $\mathbf{\mu}$. In general, the measurements $\mathbf{\mu}$ do not have a Normal distribution. However, we can use the CTL:

$$
\mu = a + b\cdot x + c\cdot x^2
$$

where x are our variables with a Normal distribution. One fact to keep in mind is that $\text{Cov}[\mu_i,\mu_j]\neq\text{Cov}[x_i,x_j]=\rho_{ij}$. The above formula represents the CTL going further than the first order.

Given $\mu$, we can express it as a function of the Normal variable $x$ and the parameters $a,b,c$.

$$
x(\mu) = \frac{\sqrt{b^2-4(a-\mu)\cdot c}-b}{2c}
$$

Using the replicas from the bootstrapping technique, I can compute a series of $\mu_i$. From this collection of $\mu_i$, I can compute the three static moments of the dataset $m_1,m_2,m_3$. Using Eq. 2.9, 2.10, 2.11, and 2.12 of the paper linked above, I can conbute $a,b,c$ and the matrix $\rho$. Then, using the Normal variables $x_i$, I can reconstruct the simplified likelihood function $-2\Delta\text{NLL} = \mathbf{x}^T \rho^{-1} \mathbf{x}$.

In our specif case, the simplified likelihood can be written as:

$$
-2\Delta\text{NLL} = (\mathbf{\hat{x}}-\mathbf{x})^T \rho^{-1} (\mathbf{\hat{x}}-\mathbf{x})
$$
where $\mathbf{\hat{x}}$ represents the results of the Combine fit.

## Environment

This repo should be run inside Combinev10
```
cmsrel CMSSW_14_1_0_pre4
cd CMSSW_14_1_0_pre4/src
cmsenv
git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
cd HiggsAnalysis/CombinedLimit
git fetch origin
git checkout v10.0.2
scramv1 b clean; scramv1 b
```
After cloning the repo:
```
mkdir ws
```
and you should change the paths in the ```ws_combine/condor_fit.sh```, line 6 and line 7.

## Inputs
The inputs are stored in the folder ```inputs```:

- htt_tt.input.root
- htt_tt_125_8TeV.txt

The only difference in the datacard, compared to the one used in the EFT tutorial, is that the numbers in the ```observation``` line are replaced by ```-1```. This is done as in the bootstap replicas, the number of events is always different, and this make the following steps easier.

## Running the fit with Combine
```
text2workspace.py inputs/htt_tt_125_8TeV.txt -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel --PO verbose --PO 'map=htt_tt_0_8TeV/ggH:r_high[1,-5,10]' --PO 'map=htt_tt_1_8TeV/ggH:r_low[1,-5,10]' --PO 'map=htt_tt_2_8TeV/ggH:r_low[1,-5,10]' -m 125 -o ws_mu_ggH.root

combine -M MultiDimFit --algo grid -d ws_mu_ggH.root -P r_high -n .high -m 125
combine -M MultiDimFit --algo grid -d ws_mu_ggH.root -P r_low -n .low -m 125

plot1DScan.py higgsCombine.high.MultiDimFit.mH125.root --POI r_high -o scan_ggH_high
plot1DScan.py higgsCombine.low.MultiDimFit.mH125.root --POI r_low -o scan_ggH_low
```
I scale the most populated category by ```r_high``` and the other two by ```r_low```.

This represents our reference. The simplified likelihood should be able to reproduce this result.

## Extracting bootstrap replicas
The bootstrap replicas are extracted with the following ROOT macro:
```
root -l 'bootstrapping.C(1000)'
```
I have tried to use pyroot, but there are issues when opening many files. The old C++ macro works much much better, without any doubt. 

The argument represents the number of bootstrap samples to be extracted.

The bootstrap is done by sampling the weights to associate to each event. For each event, a set of $N$ ```int``` numbers are extracted from a Poisson distribution with mean 1, where $N$ is the number of replicas to generate.

Tricky part: the dataset is stored as a histogram. Thus, I decided to split the content of each bin into a set of events, and then I sample the weights to associate to each event. 

[TO THINK ABOUT IT] The number of events in this dataset is not very large, especially in the category ```2```. Is this an issue? Are my replicas really representative of the original dataset?

## Running the fits on the bootstrap replicas
The fits are submitted on condor.
```
cd ws_combine
condor_submit condor_fit.sub
```
The number of fits to perform is to be specified in the ```queue``` line of the ```condor_fit.sub``` file. Of course, this number should be less or equal than the bootstrap replicas generated before.

[ISSUE] Often, one of the fits fails and the job goes to ```HOLD```. At the moment, I am simply ignoring (implemented in the script used in the following step) the cases where this happens. 

## Simplified likelihood
```
python3 extract_values.py 
```
This script does the magic:
- It extracts the values of ```r_high``` and ```r_low``` from the fits on the bootstrap replicas
- It computes the three moments of ```r_high``` and ```r_low```
- It computes the parameters $a,b,c$ and the covariance matrix $\rho$
- It computes the simplified likelihood function $-2\Delta\text{NLL} = \mathbf{x}^T \rho^{-1} \mathbf{x}$
- It compared the simplified likelihood with the combine fit

