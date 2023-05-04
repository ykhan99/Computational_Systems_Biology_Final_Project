# Flux Balance Analysis of Protein Lactate Dehydrogenase

# imports
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import cobra as cb


Loading the Model

# We will use the E. coli genome-scale metabolic model. The E. coli model has been curated for over 20 years through several publications and can be downloaded from the #BiGG Database. We will use the iML1515.xml model which is one of the latest versions, published in 2017. We will load the model below using COBRApy.

# load model
model = cb.io.read_sbml_model('iML1515.xml')

model

# Define environment

# Set a blank slate
for ex in model.exchanges:
    ex.lower_bound = 0
    ex.upper_bound = 1000
    
# glucose minimal medium
# glc_min_med = ['EX_pi_e','EX_co2_e','EX_fe3_e','EX_h_e','EX_mn2_e','EX_fe2_e','EX_glc__D_e','EX_zn2_e',
#                'EX_mg2_e','EX_ca2_e','EX_ni2_e','EX_cu2_e','EX_sel_e','EX_cobalt2_e','EX_h2o_e','EX_mobd_e',
#                'EX_so4_e','EX_nh4_e','EX_k_e','EX_na1_e','EX_cl_e','EX_o2_e','EX_tungs_e','EX_slnt_e']
glc_min_med = ['EX_pi_e','EX_co2_e','EX_fe3_e','EX_h_e','EX_mn2_e','EX_fe2_e','EX_glc__D_e','EX_zn2_e',
               'EX_mg2_e','EX_ca2_e','EX_ni2_e','EX_cu2_e','EX_sel_e','EX_cobalt2_e','EX_h2o_e','EX_mobd_e',
               'EX_so4_e','EX_nh4_e','EX_k_e','EX_na1_e','EX_cl_e','EX_tungs_e','EX_slnt_e']



for ex_id in glc_min_med:
    model.exchanges.get_by_id(ex_id).lower_bound = -1000
    print(model.exchanges.get_by_id(ex_id).id, model.exchanges.get_by_id(ex_id).name)

# set glucose exchange lower bound to -18.5
model.exchanges.get_by_id('EX_glc__D_e').lower_bound = -18.5


# Now that we have defined the environment for our simulation we can simulate the flux through the metabolic network. To do this, flux balance analysis will typically optimize a reaction called the biomass reaction. The biomass reaction consumes all of the essential building blocks for the cell (amino acids, nucleotides, lipids, etc...). Optimizing this reaction makes the assumption that the metabolic network flux is optimally producing the biomass precursors. We will implement the simulation by optimizing biomass flux below:


# This code prints out the metabolites that are consumed by the biomass reaction
# It is not necessary to run but gives some insight into the model assumptions
model.objective = 'BIOMASS_Ec_iML1515_core_75p37M'
bio_mets = list(model.reactions.get_by_id(str(model.objective.expression).split(' ')[0].split('*')[1]).metabolites.keys())
bio_stoich = list(model.reactions.get_by_id(str(model.objective.expression).split(' ')[0].split('*')[1]).metabolites.values())
bio_mets_c = []
for i in range(len(bio_mets)):
    if bio_stoich[i] < 0:
        bio_mets_c.append(bio_mets[i].name)
n_M = len(bio_mets_c)
bio_mets_c



# Simulate Flux
model.objective = 'BIOMASS_Ec_iML1515_core_75p37M'
solution = model.optimize()
print(solution.objective_value)



#Knocking Out Genes

#Genome-scale metabolic models contain a boolean mapping from genes to metabolic reactions. With this information we can simulate gene knockouts and their effect on metabolism. Let's try knocking out some genes and see what happens.
#We will knock out the gene HisD which is used for the production of the amino acid histidine in E. coli. We will then simulate growth with FBA to see if this gene was essential.

# Set a blank slate
for ex in model.exchanges:
    ex.lower_bound = 0
    ex.upper_bound = 1000
    
# glucose minimal medium
glc_min_med = ['EX_pi_e','EX_co2_e','EX_fe3_e','EX_h_e','EX_mn2_e','EX_fe2_e','EX_glc__D_e','EX_zn2_e',
               'EX_mg2_e','EX_ca2_e','EX_ni2_e','EX_cu2_e','EX_sel_e','EX_cobalt2_e','EX_h2o_e','EX_mobd_e',
               'EX_so4_e','EX_nh4_e','EX_k_e','EX_na1_e','EX_cl_e','EX_o2_e','EX_tungs_e','EX_slnt_e']

for ex_id in glc_min_med:
    model.exchanges.get_by_id(ex_id).lower_bound = -1000

# set glucose exchange lower bound to -10
model.exchanges.get_by_id('EX_glc__D_e').lower_bound = -18.5

# Simulate Gene Deletion
# Delete Gene ldhA which is essential for lactate metabolism in E. coli
with model:
    model.genes.get_by_id('b1380').knock_out()
    solution2 = model.optimize()
    print('Biomass Flux with LDH KO:', solution2.objective_value)
    
    
    
    
# Double Gene Knockouts

#In general, metabolic models allow us to simulate metabolic phenotypes at a much higher scale than what can be feasible accomplished through experiments. For example, we can simulate any combination of gene knockouts. For an organism with 1516 metabolic genes that is  2^1516 possible experiments (actually exploring this whole space is even beyond what we can simulate in a reasonable amount of compute time). Experimentally, researchers have measured the effects of double gene knockouts in various organisms. Occasionally, a double gene knockout will have a phenotype that is unexpected based on the single gene knock outs. For example, neither gene is essential individually but when they are both knocked out they become essential. This unexpected effect is generally termed "epistasis" and the specific example given here is known as "synthetic lethality". Let's try simulating some random double gene knockouts to see what the model predicts.



    
# Simulation 1 Gene Deletions
with model:
    model.genes.get_by_id('b1380').knock_out()
    solution1 = model.optimize()
    print('Biomass Flux Knock out ldhA', solution1.objective_value)
    
with model:
    model.genes.get_by_id('b3187').knock_out()
    solution1 = model.optimize()
    print('Biomass Flux Knock out ispB', solution1.objective_value)
    
with model:
    model.genes.get_by_id('b1380').knock_out()
    model.genes.get_by_id('b3187').knock_out()
    solution1 = model.optimize()
    print('Biomass Flux Knock out ldhA and ispB', solution1.objective_value)
    
    
# Simulation 2 Gene Deletions
with model:
    model.genes.get_by_id('b1380').knock_out()
    solution1 = model.optimize()
    print('Biomass Flux Knock out ldhA', solution1.objective_value)
    
with model:
    model.genes.get_by_id('b1208').knock_out()
    solution1 = model.optimize()
    print('Biomass Flux Knock out ispE', solution1.objective_value)
    
with model:
    model.genes.get_by_id('b1380').knock_out()
    model.genes.get_by_id('b1208').knock_out()
    solution1 = model.optimize()
    print('Biomass Flux Knock out ldhA and ispE', solution1.objective_value)
    
 

# Simulation 3 Gene Deletions
with model:
    model.genes.get_by_id('b1380').knock_out()
    solution1 = model.optimize()
    print('Biomass Flux Knock out ldhA', solution1.objective_value)
    
with model:
    model.genes.get_by_id('b3972').knock_out()
    solution1 = model.optimize()
    print('Biomass Flux Knock out murB', solution1.objective_value)
    
with model:
    model.genes.get_by_id('b1380').knock_out()
    model.genes.get_by_id('b3972').knock_out()
    solution1 = model.optimize()
    print('Biomass Flux Knock out ldhA and murB', solution1.objective_value)

  

# Simulation 4 Gene Deletions
with model:
    model.genes.get_by_id('b1380').knock_out()
    solution1 = model.optimize()
    print('Biomass Flux Knock out ldhA', solution1.objective_value)
    
with model:
    model.genes.get_by_id('b0085').knock_out()
    solution1 = model.optimize()
    print('Biomass Flux Knock out murE', solution1.objective_value)
    
with model:
    model.genes.get_by_id('b1380').knock_out()
    model.genes.get_by_id('b0085').knock_out()
    solution1 = model.optimize()
    print('Biomass Flux Knock out ldhA and murE', solution1.objective_value)    
    
    
    
# Simulation 4 Gene Deletions
with model:
    model.genes.get_by_id('b1380').knock_out()
    solution1 = model.optimize()
    print('Biomass Flux Knock out ldhA', solution1.objective_value)
    
with model:
    model.genes.get_by_id('b0142').knock_out()
    solution1 = model.optimize()
    print('Biomass Flux Knock out folK', solution1.objective_value)
    
with model:
    model.genes.get_by_id('b1380').knock_out()
    model.genes.get_by_id('b0142').knock_out()
    solution1 = model.optimize()
    print('Biomass Flux Knock out ldhA and folK', solution1.objective_value)    
    
    
    
    
# Simulation 5 Gene Deletions
with model:
    model.genes.get_by_id('b1380').knock_out()
    solution1 = model.optimize()
    print('Biomass Flux Knock out ldhA', solution1.objective_value)
    
with model:
    model.genes.get_by_id('b2153').knock_out()
    solution1 = model.optimize()
    print('Biomass Flux Knock out folE', solution1.objective_value)
    
with model:
    model.genes.get_by_id('b1380').knock_out()
    model.genes.get_by_id('b2153').knock_out()
    solution1 = model.optimize()
    print('Biomass Flux Knock out ldhA and folE', solution1.objective_value)    
    
    
# The models above simiulates growth when we Knock them both genes out individually but simiulates no growth when both are knocked out together, this is because neither gene is essential individually but when they are both knocked out they become essential. This unexpected effect is generally termed "epistasis" and is known as "synthetic lethality".    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
