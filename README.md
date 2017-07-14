# Predicting Arrhenius parameters from experiments

Chemical Percolation Model for Coal Devolatization (CP..D) uses Arrhenius parameters to predict statistical reaction kinetics. This module does the opposite - given time evolved reaction rate data what are the properties of the material?

---
##### **USAGE**
###### Neural Network
To train the neural net:
* Run mycpd(*sample_size>1000,experimental R vs. t*)

To predict:
* Using mytest(*experimental R vs. t*) should return the three parameters - prexponent, activation energy, and standar deviation.

**In CPD.m, include the full path to CPD binary. In the same file, changes to 'cpdoutput' filename may be necessary**
###### Gradient Descent
 
  ```matlab
  myvars=[];vars_cpd=[];output=[];
  global output,vars_cpd;
  INP=[3.1E15,6.1E4,0.81E4]'; %Starting guess
  exp1=csvread('your_exp_results');
  exp1(:,2)=exp1(:,2)*1000; %MSD will be too small on unscaled values
  fun=@Optimize;
  options=optimoptions('fminunc','Algorithm','trust-region','SpecifyObjectiveGradient',true,'OptimalityTolerance', 0.5e-5);
  [a,b]=fminunc(fun,INP,options);
  ```
###### Genetic/Particle swarm:

###### Genetic
 ```matlab
   output=[];myvars=[];
   ga(@mygenetic,3,[],[],[],[],[1E14 5000 6000],[1E15 65000 10000]); %lower and upped bound search space values
   global output,myvars;
   [n,args]=sort(-1*output); %Sort by fitness
  %myvars(args,:) should contain optimized params
```
###### Particle swarm
```matlab
  output=[];myvars=[];
   particleswarm(@mygenetic,3,[1E14 40000 6000]',[1E15 65000 11000]); %Search space
  %On-completion:
  global output,myvars;
  [n,args]=sort(-1*output); %Sort by fitness
  %myvars(args,:) should contain optimized params
```
---
