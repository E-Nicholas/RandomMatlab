function grams = solutions(mw,finalconcMM,finalvolML,stockX)
% A little function to help with those damned solution calculations.
% Function will return the required weight of a molecule given it's MW, the
% desired final concentration in MM, the final volume in mL, and a desired
% stock X value - if you want a concentrated stock solution(say you want to 
% make a 10X solution for dilution to a final [desired] at finalconcMM)
% -Eric Nicholas
finalconcMM = finalconcMM/1000;
sf = 1000/finalvolML;
prior = (mw/sf)*finalconcMM;
grams = prior * stockX;

end