# ------------------------
# BUILD
# ------------------------
# Download base image ubuntu 22.04
FROM ubuntu:22.04

# LABEL about the custom image
LABEL maintainer=""
LABEL version="0.1"
LABEL description="This is a custom Docker Image for the Quantum Enumeration Toolkit."

# Disable Prompt During Packages Installation
ARG DEBIAN_FRONTEND=noninteractive

# Update Ubuntu Software repository
RUN apt update

# Install from ubuntu repository
RUN apt install -y python3 python3-tqdm 
RUN apt install -y sagemath
RUN apt install -y git autoconf libtool parallel
RUN rm -rf /var/lib/apt/lists/*
RUN apt clean
RUN pip3 install cffi

# Add Code root path
COPY ./Cost_Estimation/ /home/Cost_Estimation/
COPY ./LatticeHeuristics/ /home/LatticeHeuristics/
COPY ./fplll_experiments/ /home/fplll_experiments/


# Build library
RUN cd /home/LatticeHeuristics/alfk/ && python3 build_alfk.py

# ------------------------
# PREPARE COST ESTIMATION
# ------------------------
VOLUME /home/Cost_Estimation/costResults

# ------------------------
# PREPARE LATTICE EXPERIMENTS
# ------------------------
RUN cd /home/fplll_experiments/ && make deep_clean fetch_fplll fetch_fplll first_compilation prepare_experiments

VOLUME /home/fplll_experiments/plots
VOLUME /home/fplll_experiments/subtree-plots
VOLUME /home/fplll_experiments/pruned/jensen-plots

# ------------------------
# RUN
# ------------------------
CMD \
  # cd /home/fplll_experiments/ ; \
	# echo "Run subtree experiments" ; make run_subtree_experiments ; \
	# echo "Plot subtree experiments" ; make plot_subtree_experiments ; \
	# echo "Run Jensen's Gap experiments" ; make run_jensen_experiments ; \
	# make plot_jensen_experiments ; \
	cd /home/Cost_Estimation/ ; \
	echo "Run Cost estimation" ; \
	python3 main_interface.py -bound lower -cache; \
	python3 main_interface.py -bound upper -cache; \
	python3 main_interface.py -bound lowerlower -cache; \
	sh 

# EXCHANGE THE ABOVE TO REDUCE TIME
# python3 main_interface.py -bound lower -cache -num-cores 24; \
# python3 main_interface.py -bound upper -cache -num-cores 24; \
# python3 main_interface.py -bound lowerlower -cache -num-cores 24; \


# # DUMMY FOR DEBUGGING
# CMD cd /home/fplll_experiments/ ; \
# 	#echo "Run subtree experiments" ; make run_subtree_experiments ; \
# 	#echo "Plot subtree experiments" ; make plot_subtree_experiments ; \ 
# 	#echo "Run Jensen's Gap experiments" ; make run_jensen_experiments ; \
# 	#make plot_jensen_experiments ; \
# 	cd /home/Cost_Estimation/ ; \
# 	#echo "Run Cost estimation" ; \
# 	#echo "Note: The 'Deprecation Warning' comes from matplotlib's LaTeX parsing and can be ignored." ; \
# 	#python3 main_interface.py -parameter-set kyber512 -circuits Query -max-depth 40 -z 1 -m 1 -y 1 ; \
# 	sh 
