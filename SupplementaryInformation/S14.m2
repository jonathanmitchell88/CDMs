####################
# Article: Convergence-divergence models: Generalizations of phylogenetic trees modeling gene flow over time
# Journal: For peer review in the Bulletin of Mathematical Biology
# Authors: Jonathan D. Mitchell (1,2*) and Barbara R. Holland (1,2)
# Affiliations: 1) School of Natural Sciences (Mathematics), University of Tasmania, Hobart, TAS, Australia, 2) ARC Centre of Excellence for Plant Success in Nature and Agriculture, University of Tasmania, Hobart, TAS, Australia.
# Corresponding author e-mail: jonathanmitchell88@gmail.com
####################

-- Groebner basis for CDM 3 under the assumptions of Theorem 13.

S=QQ;
R=S[y_1..y_7,r_1..r_7,MonomialOrder=>{Lex=>7,Lex=>7}];
I=ideal(r_1-y_4*y_5*y_6,r_2-y_2*y_3*y_4*y_6,r_3-(y_7*(1-y_6)+y_2*y_3*y_5*y_6),r_4-y_1*y_2*y_4,r_5-y_1*y_2*y_5*y_6,r_6-y_1*y_3*y_6,r_7-y_1*(y_2*y_4*y_7*(1-y_6)+y_3*y_4*y_5*y_6));
transpose gens gb I
