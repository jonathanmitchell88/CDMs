####################
# Article: Convergence-divergence models: Generalizations of phylogenetic trees modeling gene flow over time
# Journal: For peer review in the Bulletin of Mathematical Biology
# Authors: Jonathan D. Mitchell (1,2*) and Barbara R. Holland (1,2)
# Affiliations: 1) School of Natural Sciences (Mathematics), University of Tasmania, Hobart, TAS, Australia, 2) ARC Centre of Excellence for Plant Success in Nature and Agriculture, University of Tasmania, Hobart, TAS, Australia.
# Corresponding author e-mail: jonathanmitchell88@gmail.com
####################

-- CDM 5

-- Set phylogenetic tensors for two leaf labelings of CDM 5 to be equal. Determine which constraints this requires on the phylogenetic tensors. If these constraints correspond to points that are not in the parameter spaces (e.g. boundaries), then the two CDMs are distinguishable.

S=QQ;
R=S[y_1..y_9,z_1..z_9,MonomialOrder=>{GRevLex=>9,GRevLex=>9}];
I1=ideal(y_4*y_5*y_6*y_8-z_4*z_5*z_6*z_8,y_9*(1-y_8)+y_2*y_3*y_4*y_6*y_8-z_8*(z_7*(1-z_6)+z_2*z_3*z_5*z_6),y_8*(y_7*(1-y_6)+y_2*y_3*y_5*y_6)-(z_9*(1-z_8)+z_2*z_3*z_4*z_6*z_8),y_2*y_3*y_4*y_5*y_6*y_8-z_2*z_3*z_4*z_5*z_6*z_8,y_1*y_2*y_4*y_8-z_1*z_2*z_5*z_6,y_1*y_2*y_5*y_6-z_1*z_2*z_4*z_8,y_1*y_2*y_4*y_5*y_6*y_8-z_1*z_2*z_4*z_5*z_6*z_8,y_1*y_3*y_6*y_8-z_1*z_3*z_6*z_8,y_1*y_2*y_3*y_4*y_6*y_8-z_1*z_2*z_3*z_5*z_6*z_8,y_1*y_2*y_3*y_5*y_6*y_8-z_1*z_2*z_3*z_4*z_6*z_8,y_1*(y_4*y_8*(y_2*y_7*(1-y_6)+y_3*y_5*y_6)+y_2*y_5*y_6*y_9*(1-y_8))-z_1*(z_4*z_8*(z_2*z_7*(1-z_6)+z_3*z_5*z_6)+z_2*z_5*z_6*z_9*(1-z_8)),y_1*y_2*y_3*y_4*y_5*y_6*y_8-z_1*z_2*z_3*z_4*z_5*z_6*z_8);
I2=ideal(gens gb I1);
I3=eliminate(I2,{y_1,y_2,y_3,y_4,y_5,y_6,y_7,y_8,y_9});
transpose gens I3

-- CDM 4

-- Set phylogenetic tensors for two leaf labelings of CDM 4 to be equal. Determine which constraints this requires on the phylogenetic tensors. If these constraints correspond to points that are not in the parameter spaces (e.g. boundaries), then the two CDMs are distinguishable.

S=QQ;
R=S[y_1..y_9,z_1..z_9,MonomialOrder=>{GRevLex=>9,GRevLex=>9}];
I1=ideal(y_4*y_5*y_6*y_8-z_4*z_5*z_6*z_8,y_9*(1-y_8)+y_2*y_3*y_4*y_6*y_8-z_8*(z_7*(1-z_6)+z_2*z_3*z_5*z_6),y_8*(y_7*(1-y_6)+y_2*y_3*y_5*y_6)-(z_9*(1-z_8)+z_2*z_3*z_4*z_6*z_8),y_2*y_3*y_4*y_5*y_6*y_8-z_2*z_3*z_4*z_5*z_6*z_8,y_1*y_2*y_4*y_8-z_1*z_2*z_5*z_6,y_1*y_2*y_5*y_6-z_1*z_2*z_4*z_8,y_1*y_2*y_4*y_5*y_6*y_8-z_1*z_2*z_4*z_5*z_6*z_8,y_1*y_3*y_6*y_8-z_1*z_3*z_6*z_8,y_1*y_2*y_3*y_4*y_6*y_8-z_1*z_2*z_3*z_5*z_6*z_8,y_1*y_2*y_3*y_5*y_6*y_8-z_1*z_2*z_3*z_4*z_6*z_8,y_1*(y_4*y_8*(y_2*y_7*(1-y_6)+y_3*y_5*y_6)+y_2*y_5*y_6*y_9*(1-y_8))-z_1*(z_4*z_8*(z_2*z_7*(1-z_6)+z_3*z_5*z_6)+z_2*z_5*z_6*z_9*(1-z_8)),y_1*y_2*y_3*y_4*y_5*y_6*y_8-z_1*z_2*z_3*z_4*z_5*z_6*z_8,y_9-1,z_9-1);
I2=ideal(gens gb I1);
I3=eliminate(I2,{y_1,y_2,y_3,y_4,y_5,y_6,y_7,y_8,y_9});
transpose gens I3
