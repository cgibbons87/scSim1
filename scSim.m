%Spacecraft Simulator(scSim)

clear cf A store
global cf;
global A;
global store;

%name of config m-file to use(leave off ".m")
cf = 'config01';
run(cf);

%get planet motion data
store = scPlanets(cf);
A = cell2mat(store(1,:));

[sStore] = scSpacecraft(scV, scP);
[stageStore] = scStages(sStore);

scVideo2D(sStore,stageStore);