# Predicting Arrhenius parameters from experiments

Chemical Percolation Model for Coal Devolatization (CP..D) uses Arrhenius parameters to predict statistical reaction kinetics. This module does the opposite - given time evolved reaction rate data what would are the properties of the material?

---
**USAGE**
To train the neural net:
* Run mycpd(*sample_size>1000,experimental R vs. t*)

To predict:
* Using mytest(*experimental R vs. t*) should return the three parameters - prexponent, activation energy, and standar deviation.
---
