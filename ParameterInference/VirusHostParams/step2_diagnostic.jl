# Code by Marian Dominguez-Mirazo, 2024
using CSV, MCMCChains, DataFrames, DelimitedFiles
dir = "step2_MCMC/round2/"
dataids = ["data1","data2","data3"]
for dataid in dataids
	this_diagno = zeros(7,9)
	for i in collect(1:7)
		id = string(i);
		this_file = dir*dataid*"/viruschain_"*id*"_round2.csv"
		chain = DataFrame(CSV.File(this_file,header=false));
		lechain = Chains(Array(chain));
		diagno = ess_rhat(lechain);
		Nstep = Nstep = size(lechain.value)[1]
		ESSratio = diagno[1:4,2]/Nstep
		Rhat = diagno[1:4,3]
		this_diagno[i,1] = i;
		this_diagno[i,2:5] = Rhat;
		this_diagno[i,6:9] = ESSratio;
	end
	this_out_file = dir*dataid*"/diagnostic.csv"
	print(this_out_file)
	writedlm(this_out_file,this_diagno,",");
end