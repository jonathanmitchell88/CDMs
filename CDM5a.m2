-- CDM 5
-- We want to find the set of equations that can be used to solve for the parameters of CDM 5.

-- r_1=r0011
-- r_2=r0101
-- r_3=r0110
-- r_4=r0111
-- r_5=r1001
-- r_6=r1010
-- r_7=r1011
-- r_8=r1100
-- r_9=r1101
-- r_10=r1110
-- r_11=r1111
-- r_12=delta

-- To solve for the y parameters we want a set of rs with as large a cardinality as possible such that no rs are functions of other rs in the set. q1111 depends on both r_11 and r_12, with no other qs depending on r_11 or r_12. Thus, at most only one of r_11 and r_12 can be included in the set. We find the largest cardinality set that includes at most one of r_11 and r_12.

S=QQ
R=S[y_1..y_9,r_1..r_12,MonomialOrder=>{GRevLex=>9,GRevLex=>12}]
I1=ideal(r_1-y_4*y_5*y_6*y_8,r_2-(y_9*(1-y_8)+y_2*y_3*y_4*y_6*y_8),r_3-(y_7*y_8*(1-y_6)+y_2*y_3*y_5*y_6*y_8),r_4-y_2*y_3*y_4*y_5*y_6*y_8,r_5-y_1*y_2*y_4*y_8,r_6-y_1*y_2*y_5*y_6,r_7-y_1*y_2*y_4*y_5*y_6*y_8,r_8-y_1*y_3*y_6*y_8,r_9-y_1*y_2*y_3*y_4*y_6*y_8,r_10-y_1*y_2*y_3*y_5*y_6*y_8,r_11-y_1*(y_4*y_8*(y_2*y_7*(1-y_6)+y_3*y_5*y_6)+y_2*y_5*y_6*y_9*(1-y_8)),r_12-y_1*y_2*y_3*y_4*y_5*y_6*y_8)
I2=ideal(gens gb I1)

-- Eliminate y_1..y_9. Since we can include at most one of r_11 and r_12, first try eliminating r_11.

I3=eliminate(I2,{y_1,y_2,y_3,y_4,y_5,y_6,y_7,y_8,y_9,r_11})

-- I3 has no generators that are functions of r_12. Thus, r_12 is in the set. (We have not checked whether r_11 could be included in the set without r_12, but this is not necessary.) Furthermore, we can eliminate r_7 and r_10.

I4=eliminate(I3,{r_7,r_10})

-- There are no more variables that can be eliminated. We have eliminated 12 variables. Thus, the cardinality of the set is:

dim I4-12

-- gamma is another variable to include in the set. Thus, the cardinality of the set including gamma is:

dim I4-12+1

-- This is the set {r_1,r_2,r_3,r_4,r_5,r_6,r_8,r_9,r_12,gamma}.

-- Solve for y_1..y_9 as functions of r_1,r_2,r_3,r_4,r_5,r_6,r_8,r_9,r_12.

S=QQ
R=S[y_1..y_9,r_1,r_2,r_3,r_4,r_5,r_6,r_8,r_9,r_12,MonomialOrder=>{GRevLex=>9,GRevLex=>9}]
I=ideal(r_1-y_4*y_5*y_6*y_8,r_2-(y_9*(1-y_8)+y_2*y_3*y_4*y_6*y_8),r_3-(y_7*y_8*(1-y_6)+y_2*y_3*y_5*y_6*y_8),r_4-y_2*y_3*y_4*y_5*y_6*y_8,r_5-y_1*y_2*y_4*y_8,r_6-y_1*y_2*y_5*y_6,r_8-y_1*y_3*y_6*y_8,r_9-y_1*y_2*y_3*y_4*y_6*y_8,r_12-y_1*y_2*y_3*y_4*y_5*y_6*y_8)
transpose gens gb I

