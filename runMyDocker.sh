docker run -it 	-v $(pwd)/costResults:/Cost_Estimation/costResults \
				-v $(pwd)/latticeExperiments:/fplll_experiments/plots \
	 			--name qenum-container qenum-image

