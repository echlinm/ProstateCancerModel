module Networkanalysis

using SparseArrays
include("BNrules.jl")

include("BN_toolbox\\bnAsparse.jl")
include("BN_toolbox\\bnNextState.jl")
include("BN_toolbox\\bnAttractor.jl")
#Generate Network Model
model_path="model.json"
GeneNames,nv,varF,F = BNrules.getNetwork(model_path)


x = fill(0,11);
#Fog1=5,Gata1=6,Gata2=7
x[[5;6;7]] = [1;0; 0]

## CAN'T USE ANY OF THESE BECAUSE THEY AREN'T ASYNCHRONOUS
bnNextState(x,F,varF,nv)[[5;6;7]]

(),Avec=bnAsparse(F,varF,nv)
bnAttractor(Avec)



end#module
