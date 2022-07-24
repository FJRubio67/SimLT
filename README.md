# SimLT R package

In many scenarios, we may be interested in simulating survival times from a population, based on the corresponding life tables. For instance, in relative survival analysis simulation studies require the simulation of times to event based on life tables. These life tables are stratified on certain characteristics such as age, year, and deprivation level. The life table defines a piecewise-constant rate hazard function, for each combination of characteristics, which can be used to simulate the times to event.

The `SimLT` R package allows for simulatinb of times to event, based on the information in a life table. It requires constructing the piecewise-constant hazard function first as shown in:

[Simulating times to event from a Life Table](https://rpubs.com/FJRubio/LTSim)

See also: [GHSurv](https://github.com/FJRubio67/GHSurv), [LBANS](https://github.com/FJRubio67/LBANS)

```
library(devtools)
install_github("FJRubio67/SimLT")

library(SimLT)
?sim_pophaz
?simDesMatrix
?hpwexp
?chpwexp
?rsim.pwexp
```
