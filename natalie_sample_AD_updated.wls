#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Section::Closed:: *)
(*Introduction*)


(* ::Text:: *)
(*Hi Grace! This script has the skeletons of functions to find ecological equilibria, invasion fitness, and the derivative of invasion fitness.*)
(*It's based on a model of organisms evolving to transmit a symbiont horizontally (contagiously) or vertically (from parent to offspring).*)
(*Like our model, the organisms can be in different "states," in this case locations P and Q (vs. heterotrophic and phototrophic in our model)*)


(* ::Text:: *)
(*State variables: (equivalent to our W, H, P)*)


(* ::Item:: *)
(*iP: the fraction of individuals in location P who are infected*)


(* ::Item:: *)
(*iQ: the fraction of individuals in location Q who are infected *)


(* ::Text:: *)
(*Traits under selection: (equivalent to our \[Eta], \[Phi])*)


(* ::Item:: *)
(*h: the horizontal transmission rate*)


(* ::Item:: *)
(*v: the vertical transmission rate*)


(* ::Text:: *)
(*Parameters: (equivalent to our d, c, etc.)*)


(* ::Item:: *)
(*d: the probability of dispersal between patches*)


(* ::Section::Closed:: *)
(*Function descriptions*)


(* ::Text:: *)
(*This is totally optional. It just makes ?Function say something about your function. NC: So cool.*)
(*(You don't have to put them all together at the top like this. That's just a trick for making package-writing easier.)*)


(* Helper function *)
myJacobian::usage = "myJacobian finds the Jacobian of dRdt, dHdt and dMdt, evaluated with r,h and m.";
negativeEigQ::usage = "negativeEigQ returns true if all eigenvalues of a matrix are negative.";

(* Ecological equilibria *)
drdt::usage = "Rate of change in Resource Abundance";
dhdt::usage = "Rate of change in Host Abundance";
dmdt::usage = "Rate of change in Resident Abundance";
stableEquil::usage = "stableEquil finds the stable ecological equilbria of Resources R, Host H, and Resident M."; 
stableEquilList::usage = "stableEquilList find stable ecological equilibria for a list of resident trait values."; 

(* Mutant growth rate *)
mutX::usage = "mutX is our invasion fitness equation. It is the change in mutant population per unit time, per individual mutant.";
 (*mutGrowth::usage == "mutGrowth returns the growth rate of a mutant when rare";*)

(* Derivative of mutant growth rate *)
(*dGrowth::usage = "dGrowth estimates the partial derivative of the mutant growth rate, evaluated at the resident's transmission rates. Returns {\[PartialD]growth/\[PartialD]h, \[PartialD]growth/\[PartialD]v}.";
dGrowthList::usage = "dGrowthTable produces a list of derivatives of mutant growth rates when rare.";*)

(*d \[Rule] 1/2, r0 \[Rule] 25, c \[Rule] 7/10, v \[Rule] 1/10, s \[Rule] 1/20, ca \[Rule] 3/10, cc \[Rule] 2/10, lx \[Rule] 2/10, g \[Rule] 4/5, l \[Rule] 1/20, z \[Rule] 3/10, b \[Rule] ?*)

(*Adjustment to host traits.*)
(*d \[Rule] 1/2, r0 = 25, c = 74/100, v \[Rule] 1/10, s \[Rule] 3/100, ca \[Rule] 32/100, cc \[Rule] 18/100, lx \[Rule] 15/100, g \[Rule] 84/100, l \[Rule] 1/10, z \[Rule] 3/10*)
(*!!!!highlight + cmd + / to comment things!!!!*)


(* ::Text:: *)
(*Set Directory to where our notebook is! :) wahoo*)


SetDirectory[NotebookDirectory[]]


(* ::Section::Closed:: *)
(*Helper functions*)


(* ::Text:: *)
(*myJacobian finds the Jacobian of drdt, dhdt and dmdt, the functions describing the ecological dynamics, for use in stability analysis.*)


myJacobian[r_, h_, m_, d_, r0_, c_, v_, s_, ca_, cc_, b_, lx_, g_, l_, z_, um_(* parameters & um *)]:=
 {{D[drdt[o, h, m, d, r0, c, v, s, ca, cc, b, lx, g, l, z, um (* parameters & um *)], o]/. o -> r,  (* partial derivative of R with respect to R *)
 D[drdt[r, o, m, d, r0, c, v, s, ca, cc, b, lx, g, l, z, um (* parameters & um *)], o]/. o -> h,
 D[drdt[r, h, o, d, r0, c, v, s, ca, cc, b, lx, g, l, z, um (* parameters & um *)], o]/. o -> m},
 {D[dhdt[o, h, m, d, r0, c, v, s, ca, cc, b, lx, g, l, z, um (* parameters & um *)], o]/. o -> r,
 D[dhdt[r, o, m, d, r0, c, v, s, ca, cc, b, lx, g, l, z, um (* parameters & um *)], o]/. o -> h, 
 D[dhdt[r, h, o, d, r0, c, v, s, ca, cc, b, lx, g, l, z, um (* parameters & um *)], o]/. o -> m},
 {D[dmdt[o, h, m, d, r0, c, v, s, ca, cc, b, lx, g, l, z, um (* parameters & um *)], o]/. o -> r, (*After taking the derivative, change all o's to m's*)
 D[dmdt[r, o, m, d, r0, c, v, s, ca, cc, b, lx, g, l, z, um (* parameters & um *)], o]/. o -> h,
D[dmdt[r, h, o, d, r0, c, v, s, ca, cc, b, lx, g, l, z, um (* parameters & um *)], o]/. o -> m}}


(* ::Text:: *)
(*negativeEigQ takes a Jacobian and returns True if all eigenvalues are negative.*)
(*Note: I know there's sometimes cycling in the mixotroph model, so potentially might also want to check for that.*)


negativeEigQ[jacobian_] :=
AllTrue[Eigenvalues[jacobian], Negative]


(* ::Section:: *)
(*Ecological Equilibria (pre-invasion!)*)


(* ::Text:: *)
(*The functions drdt, dhdt and dmdt give the change in the abundance of each variable(r,h and m) per time step (time in units of host births).*)


drdt[r_, h_, m_, d_, r0_, c_, v_, s_, ca_, cc_, b_, lx_, g_, l_, z_, um_] := d (r0 - r) + s h 1/c - um r m + (lx um m r)/(lx + (g h))(* whole equation for dR/dt goes here *)
dhdt[r_, h_, m_, d_, r0_, c_, v_, s_, ca_, cc_, b_, lx_, g_, l_, z_, um_] := c g ((um m r)/(lx + (g h))) h - s h - ((b (um r)/(lx + (g h)) ca) + cc) m h 
dmdt[r_, h_, m_, d_, r0_, c_, v_, s_, ca_, cc_, b_, lx_, g_, l_, z_, um_] := v ((1 - (z um)) ((b (um r)/(lx + (g h)) ca) + cc)) m h - l m


(* ::Text:: *)
(*drdt*)


d (r0 - r) + s h 1/c - um r m + (lx um m r)/(lx + (g h))


(* ::Text:: *)
(*dhdt*)


c g ((um m r)/(lx + (g h))) h - s h - ((b (um r)/(lx + (g h)) ca) + cc) m h 


(* ::Text:: *)
(*dmdt*)


v ((1 - (z um)) ((b (um r)/(lx + (g h)) ca) + cc)) m h - l m


(* ::Text:: *)
(*Here lies our original stable equil code! *)


(* ::Text:: *)
(*stableEquil finds the stable equilibria of our system for each variable; Resources, r ; Residents, m; and Hosts, h.*)
(*"precision" is for the input "WorkingPrecision" in NSolve, which can be helpful if NSolve is being frustrating.*)
(*I'm using a module because I want to store the equilibria temporarily so I can do something with them (namely, discard the unstable equilibria).*)
(**)
(*The last bit of this function has an anonymous function. Anonymous functions use # as the variable (or #1, #2, etc for multiple variables) and end with &*)


stableEquil[um_, b_, d_, r0_, c_, v_, s_, ca_, cc_, lx_, g_, l_, z_,  precision_: MachinePrecision] := 
Module[{equilibria},
  (* Find all equilibria *)
  equilibria = NSolve[{drdt[r, h, m, d, r0, c, v, s, ca, cc, b, lx, g, l, z, um] == 0, dhdt[r, h, m, d, r0, c, v, s, ca, cc, b, lx, g, l, z, um] == 0, dmdt[r, h, m, d, r0, c, v, s, ca, cc, b, lx, g, l, z, um] == 0, r >= 0, h >= 0, m >= 0} , {r, h, m}, Reals, WorkingPrecision -> precision];
  (*Reals throws out imaginary roots*)
  
  (* Select the stable equilibria *)
  equilibria = Select[equilibria, negativeEigQ[myJacobian[r/.#, h/.#, m/.#, d, r0, c, v, s, ca, cc, b, lx, g, l, z, um]]&];
  
  (* Have equilibria be returned as numbers instead of iP \[Rule] number *)
  equilibria = Map[{r/.#, h/.#, m/.#}&, equilibria];
  
  Return[equilibria]
]


Clear[stableEquil]


(* ::Subsubsection::Closed:: *)
(*** *THIS WAS A TEST, DO NOT HAVE TO RUN ME IN SUBSEQUENT RUNS *** *)


(* ::Text:: *)
(*Getting equilibrias without checking stability:*)


equil[d_, r0_, c_, v_, s_, ca_, cc_, lx_, g_, l_, z_, um_, b_, precision_: MachinePrecision] := (*When defining a function make sure first alphabet is lower-cased*)
Module[{equilibria},
  (* Find all equilibria *)
  equilibria = NSolve[{drdt[r, h, m, d, r0, c, v, s, ca, cc, b, lx, g, l, z, um] == 0, dhdt[r, h, m, d, r0, c, v, s, ca, cc, b, lx, g, l, z, um] == 0, dmdt[r, h, m, d, r0, c, v, s, ca, cc, b, lx, g, l, z, um] == 0, r >= 0, h >= 0, m >= 0} , {r, h, m}, Reals, WorkingPrecision -> precision];
  
  (* Have equilibria be returned as numbers instead of iP \[Rule] number *)
  equilibria = Map[{r/.#, h/.#, m/.#}&, equilibria];
  
  Return[equilibria]
]


(* ::Subsubsection:: *)
(*stableEquilList function! takes  a list of resident trait values as input at which to find stable equilibria . *)


(* ::Text:: *)
(*-> Does not apply to this model but good to know! Each entry of the list should be of the form {h, v}, where h is the horizontal transmission rate and v is the vertical transmission rate; e.g. {{0, 0}, {0, 0.1}, {0, 0.2}}*)
(**)
(*Output is returned in the form of a list with entries of the form {um,  equilibria}; e.g. *)
(*{{0, 0, {{0, 0}}}, (* only one stable equilibrium here *)*)
(*{0.5, 0.1, {{0.1,  0.5}, {0.24, 0.73}}}} (* Multiple stable equilibria here *)*)
(**)
(*I think the error message came from accessing parts of traitList that didn't exist. The original traitList was a list of pairs of traits (because we were dealing with two traits evolving) but ours is a list of just one. So accessing trait[[1]] and trait[[2]] was accessing things that didn't exist.*)
(**)
(*I also changed some other things while I was trying to figure out what was going on, just to make sure they weren't causing any issues.*)
(*ps. precision_:MachinePrecision means that's an optional input (the ":" starts the input it will give if you don't put anything)*)


stableEquilList[traitList_, d_, r0_, c_, v_, s_, ca_, cc_, lx_, g_, l_, z_, precision_:MachinePrecision] := 
Module[{tmpEquil}, 

(* Make a table of stable equilibria *)
tmpEquil=Table[stableEquil[trait[[1]], trait[[2]], d, r0, c, v, s, ca, cc, lx, g, l, z, precision], {trait, traitList}]; (*um of traitList, previous*) (*now its trait 1 and 2 of traitList!*)
(* Table[ stuff, {um, traitList}] means make a list where each entries does a thing with one entry (called um) of traitList. It's a bit like lapply in R. I removed um, because we're starting to look at two traits. *)

(* Sort the equilibria for each pair of uptake rates *)
tmpEquil = Map[Sort, tmpEquil];

(* Add the trait values that produced the equilibrium/equilibria *)
tmpEquil = Table[{traitList[[i, 1]], traitList[[i, 2]],tmpEquil[[i]]}, {i, Length[traitList]}];

Return[tmpEquil]]


?Sort
?Part


(* ::Text:: *)
(*Naming traitList to be two values. *)
(*Note that previous analyses, used to be 1/10 of a fidelity value at every run, I'm assuming that to now they need to be of equal length? So I'm slicing the values up to an even finer degree.*)


(* ::Text:: *)
(*Here, we are using Mathematica to solve our jacobian, and testing whether our eigenvalue turned out negative (negative eigenvalue indicates that our equilibrium is stable.)*)
(*Our Jacobian is a square matrix, and there will exist 3 eigenvalues(i.e. three roots maybe?) that satisfies the stability equation.*)
(*myJacobian[r_, h_, m_, d_, r0_, c_, v_, s_, ca_, cc_, b_, lx_, g_, l_, z_, um_]*)


myJacobian[25, 0, 0, 1/2, 25, 7/10, 1/10, 1/20, 5/10, 2/10, 2/10, 2/10, 4/5, 1/20, 3/10, 1/20]
negativeEigQ[%] 


myJacobian[10, 0, 0, 6/10, 10, 9/10, 3/100, 3/10, 3/10, 1/5, 0.75, 1/20, 17/20, 1/10, 2/10, 1/20]
negativeEigQ[%] 


(* ::Text:: *)
(*negativeEigQ outputs true, meaning that the eigenvalue is negative, and that this is a stable equilibrium.*)


(* ::Subsection:: *)
(*traitListuno initialization, traitLister table and baseline parameter values here:*)


traitListuno = {0,1/40, 2/40, 3/40, 4/40, 5/40, 6/40, 7/40, 8/40, 9/40, 10/40, 11/40, 12/40, 13/40, 14/40, 15/40, 16/40, 17/40, 18/40, 19/40, 20/40, 21/40, 22/40, 23/40, 24/40, 25/40, 26/40, 27/40, 28/40, 29/40, 30/40, 31/40, 32/40, 33/40, 34/40, 35/40, 36/40, 37/40, 38/40, 39/40, 1}


(* ::Text:: *)
(*From original model code in Rstudio, Parameter values are as follows:*)


(*d \[Rule] 0.5, r0 \[Rule] 25, c \[Rule] 7/10, v \[Rule] 1/10, s \[Rule] 1/20, ca \[Rule] 3/10, cc \[Rule] 2/10, lx \[Rule] 2/10, g \[Rule] 4/5, l \[Rule] 1/20, z \[Rule] 3/10, b \[Rule] ?*)


traitLister = 
Flatten[Table[{traitListuno[[i]], traitListuno[[j]]}, 

{i, 1, 41, 1}, {j, 1, 41, 1}],1]


(* ::Chapter:: *)
(*ListDensityPlots of Resources, Host, and Symbiont.*)


(* ::Section:: *)
(*Resources :Making ListDensityPlots with newfound dual trait function*)


(* ::Subsubsection::Closed:: *)
(*Hopeless attempts of editing our lists, keeping it for the laughs..*)


r1eq[p_,q_] := ldpr[traitLister, 1/2, 25, 7/10, 1/10, 1/20, 3/10, 2/10, 2/10, 4/5, 1/20, 3/10][[p,q]]


r2eq[o_,u_,t_] := ldpr[traitLister, 1/2, 25, 7/10, 1/10, 1/20, 3/10, 2/10, 2/10, 4/5, 1/20, 3/10][[o, u, {t}]]


(* ::Text:: *)
(*- SELECT Instructions here! - *)


Clear[rae]


(*Resolves having to use r1eq and r2eq by utilizing "Select"!*)
rae = {{1, 2}, {1, 3}, {2, 1}, {1, 4}}
Select[rae, #[[2]] > 2&] (*check the second entry within the list, and select for the entries that has numbers that are greateer than 2*)
(*# variable that you're not naming, takes each item from the list. When creating an anonymous function with a # it has to end with &(!!!)*)
(*#1 + #2&[3, 4]*) 


sillyFun[x_] := x[[2]] > 2
Select[{{1, 2}, {1, 3}, {2, 1}, {1, 4}}, sillyFun]
(*They are equilvalent*)


(* ::Text:: *)
(*Select instructions end .*)


?Flatten


(* ::Text:: *)
(*?Couldn't figure out how to define levels. *)


?Select 
(*Double check again what you can do with Select!*)


(* ::Text:: *)
(*Testing parenthesis's effect on Select, select feels very strongly about how your list is structured. Which bothers me :/ ykw, probably have to do select within our function huh. *sobs**)
(*Think about how to tell it to do this: if there are two equilibria select for the one that's not 25/ Select for the first one......... U w U*)


test1 = {{{1}},{{2},{3}},{{4}},{{5}}}
Select[test1, #[[1]] < 5&]


(* ::Text:: *)
(*Okay. It isn't working. This might work: Flattening to level 3, and then doing "takelist" to allow specifying the lengths of the sublists; and then try doing select. is the best i can come up with lolll.*)


(*Extract[%76[[1,2]],{1,1}]*)


(* ::Text:: *)
(*Let' s try dropping, drop every other element .*)


Drop[fldpr, {2, 30, 2}] (*element 2 to 30, drop every other element.*)


(* ::Text:: *)
(*Okay, but now how do I combine this back to the original list? Natalie is again taking a longwinded route. Natalie is better off doing mindless copying and pasting. *)
(*But it's okay let's try collecting the actual list first. *)


?Part


?Flatten


(* ::Input:: *)
(*(* eekhmmm ... test.com *)*)


r2eq[3,3,1]
Clear[r3];ldpr[traitLister, 1/2, 25, 7/10, 1/10, 1/20, 3/10, 2/10, 2/10, 4/5, 1/20, 3/10][[3,3, {1}]]


(*r1 := ldpr[traitLister, 1/2, 25, 7/10, 1/10, 1/20, 3/10, 2/10, 2/10, 4/5, 1/20, 3/10][[1,3]]
r2 := ldpr[traitLister, 1/2, 25, 7/10, 1/10, 1/20, 3/10, 2/10, 2/10, 4/5, 1/20, 3/10][[2,3, {1}]]
r3 := ldpr[traitLister, 1/2, 25, 7/10, 1/10, 1/20, 3/10, 2/10, 2/10, 4/5, 1/20, 3/10][[3,3, {1}]]
r4 := ldpr[traitLister, 1/2, 25, 7/10, 1/10, 1/20, 3/10, 2/10, 2/10, 4/5, 1/20, 3/10][[4,3]]
r5 := ldpr[traitLister, 1/2, 25, 7/10, 1/10, 1/20, 3/10, 2/10, 2/10, 4/5, 1/20, 3/10][[5,3]]
r6 := ldpr[traitLister, 1/2, 25, 7/10, 1/10, 1/20, 3/10, 2/10, 2/10, 4/5, 1/20, 3/10][[6,3]]
r7 := ldpr[traitLister, 1/2, 25, 7/10, 1/10, 1/20, 3/10, 2/10, 2/10, 4/5, 1/20, 3/10][[7,3]]
r8 := ldpr[traitLister, 1/2, 25, 7/10, 1/10, 1/20, 3/10, 2/10, 2/10, 4/5, 1/20, 3/10][[8,3]]
r9 := ldpr[traitLister, 1/2, 25, 7/10, 1/10, 1/20, 3/10, 2/10, 2/10, 4/5, 1/20, 3/10][[9,3]]*)


(*?Select or Gather or Map*)


req (*:*)= {r1eq[1,3],r2eq[2,3,1],r2eq[3,3,1],r1eq[4,3],r1eq[5,3],r1eq[6,3],r1eq[7,3],r1eq[8,3],r1eq[9,3],r1eq[10,3],
r1eq[11,3],r1eq[12,3],r1eq[13,3],r1eq[14,3],r1eq[15,3],r1eq[16,3],r1eq[17,3],r1eq[18,3],r1eq[19,3],r1eq[20,3],
r1eq[21,3],r1eq[22,3],r1eq[23,3],r1eq[24,3],r1eq[25,3],r1eq[26,3],r1eq[27,3],
r2eq[28, 3, 1],r2eq[29, 3, 1],r2eq[30, 3, 1],r2eq[31, 3, 1],r2eq[32, 3, 1],r2eq[33, 3, 1],r2eq[34, 3, 1],r2eq[35, 3, 1],r2eq[36, 3, 1],r2eq[37, 3, 1],
r2eq[38, 3, 1],r2eq[39, 3, 1],r2eq[40, 3, 1]}


Flatten[req] (*How to remove excess curly braces*)


(*sreq := req /.List\[Rule]Sequence; good way to make a list into a sequence.*)
Clear[q] 


(* ::Text:: *)
(*List Density Plot creation whoopers .. isn't working......, U w U, so weird why aren't the data points showing up :/*)


(*Plot example using x, y, f format: ListDensityPlot[{{0,0,1},{1,0,0},{0,1,0}},Mesh\[Rule]All] Maybe try data? I have no idea...*)
(*Table, calculations; {items}, {pace}*)			                                     
(*data=Table[With[{r=RandomReal[{0,5}],t=RandomReal[{0,2Pi}]},{r Cos[t],r Sin[t],Sin[r^2]/r^2}],{10^4}];*)
(*With r as a randomreal value between 0 and 5; t the same; put that value in the equations ^, creating a list, and repeat for 10 to the power of 4 times.*)


(* ::Text:: *)
(*Axis labels aren't working either. Sol: use FrameLabel instead when doing LDP.*)


ListDensityPlot[{{0, 0, 25.}, {1/40, 1/40, 8.73537461967574}, {1/20, 1/20, 4.4285422517036315}, {3/40,3/40, 24.999999999999993}, {1/10,1/10, 25.000000000000004}, {1/8,1/8, 25.}, {3/20,3/20, 25.000000000000004}, {7/40,7/40, 25.000000000000004}, {1/5,1/5, 24.99999999999999},{9/40,9/40, 25.},{1/4,1/4, 25.},{11/40,11/40, 25.},{3/10,3/10, 25.},{13/40,13/40, 25.},{7/20,7/20, 25.},{3/8,3/8, 25.},{2/5,2/5, 25.},{17/40,17/40, 25.},{9/20,9/20,25.},{1/2,1/2, 25.},{21/40,21/40, 25.},{11/20,11/20, 25.},{23/40,23/40, 25.},{3/5,3/5, 25.},{5/8,5/8, 25.},{13/20,13/20, 25.},{27/40,27/40, 25.},{7/10,7/10, 0.515941072142575},{29/40,29/40, 0.5105976815333032},{3/4,3/4, 0.5062322035237973},{31/40,31/40, 0.5027982246647518},{4/5,4/5, 0.5002603020489491},{33/40,33/40, 0.4985930189507557},{17/20,17/20, 0.4977803312471313},{7/8,7/8, 0.49781517245114737},{9/10,9/10, 0.4986992987974486},{37/40,37/40, 0.5004433679945054},{19/20,19/20, 0.5030672570090124},{39/40,39/40, 0.5066006365392539},{1,1, 0.5110838337102814}}, ColorFunction->ColorData[{"MintColors", {0.498,25}}], DataRange -> Automatic, ColorFunctionScaling-> True, PlotLegends-> Automatic, PlotRange-> Automatic, Mesh -> All, FrameLabel -> {"fidelity, b", "Uptake Rate, um"}]


(* ::Text:: *)
(*Maybe need multiple data points? I have no idea...*)


(* ::Text:: *)
(*When um == b! Run code below to find values for variables at equilibrium.*)


(*stableEquilList[traitList_, d_, r0_, c_, v_, s_, ca_, cc_, lx_, g_, l_, z_,*)


(*stableEquilList[traitList_, d_, r0_, c_, v_, s_, ca_, cc_, lx_, g_, l_, z_,*)


(*d \[Rule] 1/2, r0 = 25, c = 74/100, v \[Rule] 1/10, s \[Rule] 3/100, ca \[Rule] 32/100, cc \[Rule] 18/100, lx \[Rule] 15/100, g \[Rule] 84/100, l \[Rule] 1/10, z \[Rule] 3/10*)


stableEquilList[traitList, 1/2, 25, 74/100, 1/10, 3/100, 32/10, 18/100, 15/100, 84/100, 1/10, 3/10]


?Nest


?Part


(* ::Subsubsection:: *)
(*Anonymous  function SAMPLEe!*)


#1 + #2&[3, 4]

{{1, 2}, {3, 4}} // MatrixForm (*two slashes lets you apply a function that you forgot to write or soemthing to the stuff before the function.*)
3 // #^2 & (*when you forget to write a function*)


(* ::Subsubsection:: *)
(*1. Baseline R* value generator here.*)


stableEquilr[um_, b_, d_, r0_, c_, v_, s_, ca_, cc_, lx_, g_, l_, z_,  precision_: MachinePrecision] := 
Module[{equilibriar},
  (* Find all equilibria *)
  equilibriar = NSolve[{drdt[r, h, m, d, r0, c, v, s, ca, cc, b, lx, g, l, z, um] == 0, dhdt[r, h, m, d, r0, c, v, s, ca, cc, b, lx, g, l, z, um] == 0, dmdt[r, h, m, d, r0, c, v, s, ca, cc, b, lx, g, l, z, um] == 0, r >= 0, h >= 0, m >= 0} , {r, h, m}, Reals, WorkingPrecision -> precision];
  (*Reals throws out imaginary roots*)
  
  (* Select the stable equilibria *)
  equilibriar =  Select[equilibriar, negativeEigQ[myJacobian[r/.#, h/.#, m/.#, d, r0, c, v, s, ca, cc, b, lx, g, l, z, um]]&];
 
 (* While our the first term of equilibriar is < 25, Drop the second term.
  equilbriar = While[equilibriar #[[1]]& < 25 &; Drop #[[2]]&
  but this will likely not work if parantheses issue is not resolved*)
  
  (* Have equilibria be returned as numbers instead of r \[Rule] number *)
  equilibriar = Map[{r/.#}&, equilibriar];
  
  Return[equilibriar]
]


(*Noppers.*)
?FlattenAt
(*After that flatten them at those positions?*)
?Position 
(*Use position to figure out where 25 occurs. *)


(*Module[{temp1, temp2},temp1 = dervmut[stuff];temp2 = temp1[[1]]; Nsolve[temp2\[Equal]0];*)

(*ldp.r[traitList_, d_, r0_, c_, v_, s_, ca_, cc_, lx_, g_, l_, z_, precision_:MachinePrecision] := Module[{tempo, tempo2}, tempo = *)
ldpr[traitLister, d_, r0_, c_, v_, s_, ca_, cc_, lx_, g_, l_, z_, precision_:MachinePrecision] := 
Module[{tempo}, 
(* Make a table of stable equilibria *)
tempo=Table[stableEquilr[traiter[[1]], traiter[[2]], d, r0, c, v, s, ca, cc, lx, g, l, z, precision], {traiter, traitLister}]; (*goes through each entry of traitLister and define it as traiter*)
(* Change length of traitlister inputs [1], [2] to test which one is which axis.*)

(* Sort the equilibria for each pair of uptake rates *)
tempo = Map[Sort, tempo];

(* Add the trait values that produced the equilibrium/equilibria *)
tempo = Table[{(*traitLister[[i, 1]], traitLister[[i, 2]],*)tempo[[i]]}, {i, Length[traitLister]}];

Return[tempo]]


(* ::Text:: *)
(*original parameter values*)


(*d \[Rule] 1/2, r0 \[Rule] 25, c \[Rule] 7/10, v \[Rule] 1/10, s \[Rule] 1/20, ca \[Rule] 3/10, cc \[Rule] 2/10, lx \[Rule] 2/10, g \[Rule] 4/5, l \[Rule] 1/20, z \[Rule] 3/10, b \[Rule] ?*)


(* ::Text:: *)
(*Adjusted parameter values.*)


(*d \[Rule] 1/2, r0 = 25, c = 74/100, v \[Rule] 1/10, s \[Rule] 3/100, ca \[Rule] 32/100, cc \[Rule] 18/100, lx \[Rule] 15/100, g \[Rule] 84/100, l \[Rule] 1/10, z \[Rule] 3/10*)


ldprr = ldpr[traitLister, 1/2, 25, 7/10, 1/10, 1/20, 3/10, 2/10, 2/10, 4/5, 1/20, 3/10]


(* ::Subsubsection:: *)
(*2. UPDATED: A function that takes the first equilibria of all stuff within the list.*)


(* ::Text:: *)
(*Okay now you have all equilibrias, but with potentially useless ones embedded (situation where there's 2 equilibrias and one is 25.). Figure out how to remove them is the next step . *)


ldprr


(*Flatten[ldprr, {1,2,4}] : means flatten ldprr, but only at levels 1,2 and 4.*) 


(* ::Text:: *)
(*final ListDensityPlotR data for all R*s. Here, specifying 1,2,4 means to only flatten those levels.*)


finldprr = Flatten[ldprr, {1,2,4}][[All, 1]] (*Take everything of the 1st level but only take the first entry.*)


(* ::Text:: *)
(*Generating a list of finalist equilibrias that would go into our LDP_r.*)


listthatgoesintolpdr[traitLister, d_, r0_, c_, v_, s_, ca_, cc_, lx_, g_, l_, z_, precision_:MachinePrecision] := 
Module[{tempo}, (*it doesn't matter if you'd used tempo before because you're in a module, how cool!*)

(* Make a table of final stable equilibrias *)
tempo=finldprr;

(* Sort the equilibria for each pair of uptake rates *)
(*tempo = Map[Sort, tempo]; \[Rule] don't think we'll need to do Sort.*)

(* Add the trait values that produced the equilibrium/equilibria *)
tempo = Table[{traitLister[[i, 1]], traitLister[[i, 2]], tempo[[i]]}, {i, Length[traitLister]}];

Return[tempo]]


(* ::Text:: *)
(*d -> 0.5, r0 -> 25, c -> 7/10, v -> 1/10, s -> 1/20, ca -> 3/10, cc -> 2/10, lx -> 2/10, g -> 4/5, l -> 1/20, z -> 3/10, b -> ?*)


listthatgoesintolpdr[traitLister, 1/2, 25, 7/10, 1/10, 1/20, 3/10, 2/10, 2/10, 4/5, 1/20, 3/10]


(* ::Subsubsection:: *)
(*3. Moment of truth! Running code to create LDP. IT WORKED.*)


(* ::Text:: *)
(*Narrowing our plot range for R**)


(*Only for 0*)
?Chop 


Select[finldprr, # != 25. &]


{minRstar, maxRstar} = {Min[%], Max[%]}


(* ::Text:: *)
(*Yay LDP! *)


?ListDensityPlot


ListDensityPlot[{listthatgoesintolpdr[traitLister, 1/2, 25, 7/10, 1/10, 1/20, 3/10, 2/10, 2/10, 4/5, 1/20, 3/10]}, PlotRange -> {Full, Full, {minRstar, maxRstar}}, ClippingStyle -> {Black, White}, ColorFunction->ColorData["MintColors"], DataRange -> Automatic, ColorFunctionScaling-> True, PlotLegends-> Automatic, PlotRange-> Automatic, FrameLabel -> {"Uptake Rate, um", "fidelity, b" }]


Export["rstartrue.png", %]


(* ::Section:: *)
(*Host : Making ListDensity Plots for H* changing parameter values b, and u *)


(* ::Subsubsection:: *)
(*1. finding equilibrium values by using stableEquil   *)


stableEquilh[um_, b_, d_, r0_, c_, v_, s_, ca_, cc_, lx_, g_, l_, z_,  precision_: MachinePrecision] := 
Module[{equilibriah},
  (* Find all equilibria *)
  equilibriah = NSolve[{drdt[r, h, m, d, r0, c, v, s, ca, cc, b, lx, g, l, z, um] == 0, dhdt[r, h, m, d, r0, c, v, s, ca, cc, b, lx, g, l, z, um] == 0, dmdt[r, h, m, d, r0, c, v, s, ca, cc, b, lx, g, l, z, um] == 0, r >= 0, h >= 0, m >= 0} , {r, h, m}, Reals, WorkingPrecision -> precision];
  (*Reals throws out imaginary roots*)
  
  (* Select the stable equilibria *)
  equilibriah =  Select[equilibriah, negativeEigQ[myJacobian[r/.#, h/.#, m/.#, d, r0, c, v, s, ca, cc, b, lx, g, l, z, um]]&];

  (* Have equilibria be returned as numbers instead of r \[Rule] number *)
  equilibriah = Map[{h/.#}&, equilibriah];
  
  Return[equilibriah]
]


(* ::Subsubsection:: *)
(*2. Generate a list of All possible Hstars that Jacobian produces negative eigenvalues (real roots)*)


(* ::Text:: *)
(*its cool how you never even define traiter but matha knows what traiter[[1]] means .*)


ldph[traitLister, d_, r0_, c_, v_, s_, ca_, cc_, lx_, g_, l_, z_, precision_:MachinePrecision] := 
Module[{tempo}, 
(* Make a table of stable equilibria *)
tempo=Table[stableEquilh[traiter[[1]], traiter[[2]], d, r0, c, v, s, ca, cc, lx, g, l, z, precision], {traiter, traitLister}];

(* Sort the equilibria for each pair of uptake rates *)
tempo = Map[Sort, tempo];

(* Add the trait values that produced the equilibrium/equilibria *)
tempo = Table[{(*traitLister[[i, 1]], traitLister[[i, 2]], *)tempo[[i]]}, {i, Length[traitLister]}];

Return[tempo]]


(* ::Text:: *)
(*Yay! Now we have all of our equilibrias! :)*)


templdph = ldph[traitLister, 1/2, 25, 7/10, 1/10, 1/20, 3/10, 2/10, 2/10, 4/5, 1/20, 3/10]


?Drop


?Cases


(* ::Text:: *)
(*Have to select for the second or third equil (i.e. just the non-zero ones hehe).*)


(* ::Subsubsection:: *)
(*3. Select for "useful equilibrias" *)


rae = Flatten[templdph, {1,2,4}]


(* ::Text:: *)
(*Trying to filter through the list and find the minimum and maximum H abundance at equilibrium. *)


Flatten[Select[%, (Length[#] == 2 && !(# == {0, 0})) || Length[#] ==3&]]


Select[%, Chop[#] > 0 &] (*!! Alexandra Knowledge !! Here > doesn't mean greater than, but everything 0 or close to 0 will be chopped. :)*)


{minHstar, maxHstar} = {Min[%], Max[%]} (*!! Alexandra Knowledge Finding the minimum and maximum of H*!!*)


?Partition


rae2 = Map[Last, rae]  (*!! Alexandra Knowledge !!*)


?Last


(* ::Subsubsection::Closed:: *)
(*4. Create a list that includes trait values! *)


listthatgoesintolpdh[traitLister, d_, r0_, c_, v_, s_, ca_, cc_, lx_, g_, l_, z_, precision_:MachinePrecision] := 
Module[{tempo}, (*it doesn't matter if you'd used tempo before because you're in a module, how cool!*)

(* Make a table of final stable equilibrias *)
tempo=rae2;

(* Sort the equilibria for each pair of uptake rates *)
(*tempo = Map[Sort, tempo]; \[Rule] don't think we'll need to do Sort.*)

(* Add the trait values that produced the equilibrium/equilibria *)
tempo = Table[{traitLister[[i, 1]], traitLister[[i, 2]], tempo[[i]]}, {i, Length[traitLister]}];

Return[tempo]]


(* ::Text:: *)
(*d -> 0.5, r0 -> 25, c -> 7/10, v -> 1/10, s -> 1/20, ca -> 3/10, cc -> 2/10, lx -> 2/10, g -> 4/5, l -> 1/20, z -> 3/10, b -> ?*)


listthatgoesintolpdh[traitLister, 1/2, 25, 7/10, 1/10, 1/20, 3/10, 2/10, 2/10, 4/5, 1/20, 3/10]


(*Alexandra Knowledge : PlotRange -> {Full(x), Full(y), {minHstar, maxHstar}(f)}, ClippingStyle -> {"Green", "Pink"}(*Given that value is NOT within plot range do green, for min; and do pink for max*)*)


ListDensityPlot[{listthatgoesintolpdh[traitLister, 1/2, 25, 7/10, 1/10, 1/20, 3/10, 2/10, 2/10, 4/5, 1/20, 3/10]}, ColorFunction->ColorData["MintColors"], DataRange -> Automatic, ColorFunctionScaling-> True, PlotLegends-> Automatic, PlotRange -> {Full, Full, {minHstar, maxHstar}}, ClippingStyle -> {White, Pink}, FrameLabel -> {"Uptake Rate, um","fidelity, b"}]


Export["hstartrue.png", %]


(* ::Subsubsection::Closed:: *)
(*Whatt went wrong : ( Nope! Nothing to see here.*)


(*for every element within rae(already flattened equilibrias), if there are 2 || 3 equils pick the one that isn't zero, otherwise just take the one 0 equil
Finally match them up with their respective um and b values, and create a list out of that.*)
finldph = Module[{rae},

rae = Flatten[templdph, {1,2,4}];

rae = ;

Return[rae]
]


yrae = Cases[rae, {i_ , j_, k_} /; k > 0  :> k ] (*Nah*)


(*Pick[#, #[[1]] > 8 && #[[2]] > 8 & /@ #] &@list[[All, {2, 3}]] or Select[#[[1]] > 8 && #[[2]] > 8 &]@list[[All, {2, 3}]]*)


?DeleteCases


?If


Cases[rae, {_, _}][[All, 2]] (*Except[0]*)


If[Cases[rae, {_, _}], Select[#[[2]] > 0 &] , Select[ #2[[1]] = 0 &]] (*Wants to do the if for every single element within rae*)


(* ::Text:: *)
(*!REC: Think about how to do if there's two equil select for the second one. Ans: LAST! :)*)


?Position


?Select


?Sort


Sort[templdph, {1, 2, 3}]


(* ::Section:: *)
(*Resident : Making ListDensity Plots for M* changing parameter values b, and u *)


(* ::Subsubsection:: *)
(*1. finding equilibrium values by using stableEquil   *)


stableEquilm[um_, b_, d_, r0_, c_, v_, s_, ca_, cc_, lx_, g_, l_, z_,  precision_: MachinePrecision] := 
Module[{equilibriam},
  (* Find all equilibria *)
  equilibriam = NSolve[{drdt[r, h, m, d, r0, c, v, s, ca, cc, b, lx, g, l, z, um] == 0, dhdt[r, h, m, d, r0, c, v, s, ca, cc, b, lx, g, l, z, um] == 0, dmdt[r, h, m, d, r0, c, v, s, ca, cc, b, lx, g, l, z, um] == 0, r >= 0, h >= 0, m >= 0} , {r, h, m}, Reals, WorkingPrecision -> precision];
  (*Reals throws out imaginary roots*)
  
  (* Select the stable equilibria *)
  equilibriam =  Select[equilibriam, negativeEigQ[myJacobian[r/.#, h/.#, m/.#, d, r0, c, v, s, ca, cc, b, lx, g, l, z, um]]&];
 
 (* While our the first term of equilibriar is < 25, Drop the second term.
  equilbriar = While[equilibriar #[[1]]& < 25 &; Drop #[[2]]&
  but this will likely not work if parantheses issue is not resolved*)
  
  (* Have equilibria be returned as numbers instead of r \[Rule] number *)
  equilibriam = Map[{m/.#}&, equilibriam];
  
  Return[equilibriam]
]


(* ::Subsubsection:: *)
(*2. Generate a list of All possible Hstars that Jacobian produces negative eigenvalues (real roots)*)


ldpm[traitLister, d_, r0_, c_, v_, s_, ca_, cc_, lx_, g_, l_, z_, precision_:MachinePrecision] := 
Module[{tempo}, 
(* Make a table of stable equilibria *)
tempo=Table[stableEquilm[traiter[[1]], traiter[[2]], d, r0, c, v, s, ca, cc, lx, g, l, z, precision], {traiter, traitLister}];

(* Sort the equilibria for each pair of uptake rates *)
tempo = Map[Sort, tempo];

(* Add the trait values that produced the equilibrium/equilibria *)
tempo = Table[{tempo[[i]]}, {i, Length[traitLister]}];

Return[tempo]]


templdpm = ldpm[traitLister, 1/2, 25, 7/10, 1/10, 1/20, 3/10, 2/10, 2/10, 4/5, 1/20, 3/10]


raeb = Flatten[templdpm, {1,2,4}]


(* ::Text:: *)
(*Got our finalists :*)


raeb2 = Map[Last, raeb] 


(* ::Text:: *)
(*finding max and min Mstar . That are of value to us .*)


Flatten[Select[raeb, (Length[#] == 2 && !(# == {0, 0})) || Length[#] ==3&]] (*! is NOT*) 


traeb = Select[%, Chop[#] > 0 &] (*This is tr\[EAcute]s cool. Not even kidding*) 


{minMstar, maxMstar} = {Min[traeb], Max[traeb]} 


(* ::Subsubsection:: *)
(*3. Create a list that includes trait values! *)


finalistm[traitLister, d_, r0_, c_, v_, s_, ca_, cc_, lx_, g_, l_, z_, precision_:MachinePrecision] := 
Module[{tempo}, (*it doesn't matter if you'd used tempo before because you're in a module, how cool!*)

(* Make a table of final stable equilibrias *)
tempo=raeb2;

(* Add the trait values that produced the equilibrium/equilibria *)
tempo = Table[{traitLister[[i, 1]], traitLister[[i, 2]], tempo[[i]]}, {i, Length[traitLister]}];

Return[tempo]]


finalistm[traitLister, 1/2, 25, 7/10, 1/10, 1/20, 3/10, 2/10, 2/10, 4/5, 1/20, 3/10]


?Select


(* ::Subsubsection::Closed:: *)
(*4. Finding umMin, for evol function. ignore this subsection.*)


(*Rather than doing this, it'll be more efficient to just stare at finalistm! after umMin, the string of consequtive elements would also have to produce positive Mstars.*)


(* ::Text:: *)
(*{0,1/40,1/20,3/40,1/10,1/8,3/20,7/40,1/5,9/40,1/4,11/40,3/10,13/40,7/20,3/8,2/5,17/40,9/20,19/40,1/2,21/40,11/20,23/40,3/5,5/8,        13/20,27/40,7/10,29/40,3/4,31/40,4/5,33/40,17/20,7/8,9/10,37/40, 19/20,       39/40,1}*)


(*Select the value of um that is the lowest.*)


(* ::Text:: *)
(*umMinList*)


(*prosp = Select[%42, #[[3]] !=  0&]*)


Extract[umMinList, 27]


Clear[way]


Select[%42, #[[2]] == 7/10 && #[[3]] != 0 &][[1]] (*%42 was finalistm*)


(* ::Text:: *)
(*This isn' t umMin!!!!*)


Clear[umMin]


umMin = Select[%, #[[3]] > 0 & , 1][[All, 1]]
(*umMin = 1/40*)


(*Nevermind, she found out that stableequil takes um_, b_ as parameter inputs! :) Wheee*)


(* ::Subsubsection:: *)
(*5. Create LDP_M!*)


ListDensityPlot[{finalistm[traitLister, 1/2, 25, 7/10, 1/10, 1/20, 3/10, 2/10, 2/10, 4/5, 1/20, 3/10]}, ColorFunction->ColorData["MintColors"], DataRange -> Automatic, ColorFunctionScaling-> True, PlotLegends-> Automatic, PlotRange -> {Full, Full, {minMstar, maxMstar}}, ClippingStyle -> {White, Pink}, FrameLabel -> {"Uptake Rate, um", "fidelity, b"}]


Export["mstarreal.png", %]


(* ::Subsection::Closed:: *)
(*Don't need this table anymore, but just keeping it cuz it's cool.*)


(*data=Flatten[Table[{i,j,Sin[i+ j]},{j,0,2Pi,0.5},{i,-3,3}],1];*)

axisdpoints4r = Flatten[Table[ (*p2 = updateval[0];*)
{traitListuno[[i]], traitListuno[[j]], (*sreq2[[i]]*)}, 

{i, 1, 40, 1}, {j, 1, 40, 1}

],1]


(* ::Text:: *)
(*List Density Plot attempt with data form table "axisdpoinst4r"*)


ListDensityPlot[axisdpoints4r]


ldpr[traitLister, 1/2, 25, 7/10, 1/10, 1/20, 3/10, 2/10, 2/10, 4/5, 1/20, 3/10][[2,3, {1}]]


(* ::Subsection::Closed:: *)
(*Jacobian below is written in a fancier way.*)


myJacobian[25, 0, 0, 1/2, 25, 7/10, 1/10, 1/20, 5/10, 2/10, 2/10, 2/10, 4/5, 1/20, 3/10, 11/20]//MatrixForm
negativeEigQ[%] 


(* ::Text:: *)
(*ListDensityPlot*)
(*ColorFunctions!*)


?ListDensityPlot 
?With
(*x and y are datapoints, f is the color that you want your datapoints to be.*)
(*One Plot for each variable, only plot equilibrium with zero when there isn't any other equilibrias.*)
(*List density plot, for each variable r, h and m, as a function of um? That'll be cool. Goood luckkk.*)


(* ::Subsection::Closed:: *)
(*Tried making listdensityplots using extract, and it didn't work well*)


(* ::Text:: *)
(*Doesn't work... An option must be a rule or a list of rules. -__-"*)


ListDensityPlot[{{%9[[1, 1]], Extract[%9[[1, 2]], {1, 1}], ColorData["MintColors", 1]}, {%9[[2, 1]], Extract[%9[[2, 2]], {1, 1}], ColorData["MintColors", 1]}, {%9[[3, 1]], Extract[%9[[3, 2]], {1, 1}], ColorData["MintColors", 0]}, {%9[[4, 1]], Extract[%9[[4, 2]], {1, 1}], 1}, {%9[[5, 1]], Extract[%9[[5, 2]], {1, 1}]}, ColorFunction->ColorData["MintColors"]}]


(* ::Text:: *)
(*stableEquil[d_, r0_, c_, v_, s_, ca_, cc_, lx_, g_, l_, z_, b_]*)


%76[[1,1]]
Extract[%76[[1,2]],{1,1}]
Extract[%76[[2,2]],{1,1}]


ColorData[Gradients]


(* ::Text:: *)
(**)


(* ::Subsection::Closed:: *)
(*Can ignore this section - I was trying other parameter values for stabliEquilList. *)


(*d \[Rule] 0.6, r0 \[Rule] 10, c \[Rule] 0.9, v \[Rule] 0.03, s \[Rule] 0.3, ca \[Rule] 3/10, cc \[Rule] 2/10, lx \[Rule] 1/20, g \[Rule] 17/20, l \[Rule] 1/10, z \[Rule] 2/10, b \[Rule] 0.9*)


stableEquilList[traitList, 0.6, 10, 0.9, 0.03, 0.3, 0.3, 0.2, 0.05, 0.85, 0.1, 0.2, 0.9]


(* ::Text:: *)
(*Running stableEquilList to see what the equilibrium looks like! :- When b = 0.1*)


stableEquilList[traitList, 0.6, 10, 0.9, 0.03, 0.3, 0.3, 0.2, 0.05, 0.85, 0.1, 0.2, 0.45]


(* ::Text:: *)
(*Running stableEquilList to see what the equilibrium looks like! :- When b = 0.2*)


stableEquilList[traitList, 0.5, 25, 0.7, 0.1, 0.05, 0.3, 0.2, 0.2, 0.8, 0.05, 0.3, 1.]


(* ::Text:: *)
(*b doesn't seem to be changing our equilibrium points.*)


(*BIG LARGE SPACE BEFORE EVOLUTION HAPPENS. WHEEE!*)








(* ::Chapter:: *)
(*Moving onto Evolution Dynamics.*)


(* ::Section:: *)
(*Mutant growth rate (invasion fitness)*)


(* ::Text:: *)
(*NOTE: We don't have to do this so painfully! We don't need mutX! Getting growth rate derivative is way simpler. This is still adaptive dynamics, just the painful version for when things can exist in multiple states. *)


(* ::Text:: *)
(*To do adaptive dynamics, we need a formula for a mutant's growth rate when rare. The mutant's growth rate will probably depend on  its own traits, the resident's traits, and the environment.*)
(*Finding this growth rate falls in the gray area between coming up with a model and doing math on something you already have.*)


(* ::Text:: *)
(*In both the mixotroph model and the transmission model, one difficulty is that the mutant' s growth rate probably depends on what state it is in.*)
(*Is it a heterotroph or phototroph? Or, in the transmission model, is the mutant in location P or location Q?*)


(* ::Text:: *)
(*One way to find mutant fitness in this case to write down something like a Leslie matrix, that includes mutant birth, death, and transition between different states.*)
(*The matrix should look something like this:*)
(*X = ({*)
(* {H -> H, P -> H},*)
(* {H -> P, P -> P}*)
(*})*)
(*Where H -> H represents the expected number of heterotrophs in the next time step produced by 1 heterotroph in the current time step. That includes a heterotroph "producing" itself in the next time step by remaining alive and a heterotroph.*)
(*In the transmission model, I had to do some extra math at this point, because X depended on the population size, and I was interested in the infinite population case.*)
(**)
(*One you have X, its leading eigenvalue gives you the mutant growth rate.*)


(* ::Text:: *)
(*I have stored the function below into dervmutX.*)
(*make dervmutx take r and h  at equil.*)


(*v ((1 - (z um)) ((b (um r)/(lx + (g h)) ca) + cc)) m h - l m*)


v ((1 - (z un)) ((b (un r)/(lx + (g h)) ca) + cc)) h - l


Simplify[D[mutX[r, h, d, r0, c, v, s, ca, cc, lx, g, l, z, b], un]]


(* ::Text:: *)
(*The derivative of mutX will come out to be : constant - um, causing a negative slope with um at the x - axis .*)
(*Currently system is evolving towards lower um== lower host fitness, but there's still fidelity that might changes things.*)
(*-> the host wants a good partner it'll have to be faithful even though the partner would be the worse it could possibly be.*)


(* ::Subsubsection:: *)
(*DervmutX here,*)


mutX[r_, h_, d_, r0_, c_, v_, s_, ca_, cc_, lx_, g_, l_, z_, b_] := v ((1 - (z un)) ((b (un r)/(lx + (g h)) ca) + cc)) h - l (*mutX = Invasion fitness. *)
dervmutX[r_, h_, d_, r0_, c_, v_, s_, ca_, cc_, lx_, g_, l_, z_, b_]:= D[mutX[r, h, d, r0, c, v, s, ca, cc, lx, g, l, z, b], un] (*This function takes the derivative of the Invasion Fitness*)


(*take derivative of invasion fitness, inv fitness when resident trait value equals mutant*)
(*function below should give me a number that changes depending on um. Numerically solve for where that is 0.*)
(*D[mutX[r, h, d, r0, c, v, s, ca, cc, lx, g, l, z, b], un]/. un -> um*)


(*d \[Rule] 0.5, r0 \[Rule] 25, c \[Rule] 7/10, v \[Rule] 1/10, s \[Rule] 1/20, ca \[Rule] 3/10, cc \[Rule] 2/10, lx \[Rule] 2/10, g \[Rule] 4/5, l \[Rule] 1/20, z \[Rule] 3/10, b \[Rule] ?*)


(*dervmutX[r_, h_, d_, r0_, c_, v_, s_, ca_, cc_, lx_, g_, l_, z_, b_]*)


(* ::Subsection::Closed:: *)
(*Nothing to see here, keep going! lets test b=1.0, and um = (1/40)! Two distinct equilibria exists here.*)


dervmutX[Extract[%9[[2,2]],{1,1}],Extract[%9[[2,2]],{1,2}], 1/2, 25, 7/10, 1/10, 1/20, 3/10, 2/10, 1/10, 4/5, 1/20, 3/10, 1]


(*julia = stableEquil[blah, blah];
Extract[julia[[2, 2]], {1, 1}]*)


?Interpolation 
(*Guesses what's in the middle when you have a bunch of data of data points*)


NSolve[%14 == 0] 


(*Woohoo! Our first successful output wheeeeee*) (*{um \[Rule] 1.66272}*)


(* ::Subsubsection::Closed:: *)
(*b = 1.; um = (27/40)*)


dervmutX[Extract[%9[[28,2]],{1,1}],Extract[%9[[28,2]],{1,2}], 1/2, 25, 7/10, 1/10, 1/20, 3/10, 2/10, 1/10, 4/5, 1/20, 3/10, 1]


NSolve[%17 == 0]


 (*{um \[Rule] 1.56608}*)


(* ::Subsubsection::Closed:: *)
(*b = 1.; um = (7/10)*)


dervmutX[Extract[%9[[29,2]],{1,1}],Extract[%9[[29,2]],{1,2}], 1/2, 25, 7/10, 1/10, 1/20, 3/10, 2/10, 1/10, 4/5, 1/20, 3/10, 1]


NSolve[%19 == 0]


 (*{um \[Rule] 1.56268}*)


(* ::Subsubsection:: *)
(*b = 1.; um = (29/40)*)


dervmutX[Extract[%9[[30,2]],{1,1}],Extract[%9[[30,2]],{1,2}], 1/2, 25, 7/10, 1/10, 1/20, 3/10, 2/10, 1/10, 4/5, 1/20, 3/10, 1]


NSolve[%21 == 0]


(*um \[Rule] 1.5593*)


(* ::Subsubsection:: *)
(*b =1.; um =(3/4)*)


dervmutX[Extract[%9[[31,2]],{1,1}],Extract[%9[[31,2]],{1,2}], 1/2, 25, 7/10, 1/10, 1/20, 3/10, 2/10, 1/10, 4/5, 1/20, 3/10, 1]


NSolve[%23 == 0]


(*um \[Rule] 1.55595*)


(* ::Subsubsection:: *)
(*b =1.; um =(31/40)*)


dervmutX[Extract[%9[[32,2]],{1,1}],Extract[%9[[32,2]],{1,2}], 1/2, 25, 7/10, 1/10, 1/20, 3/10, 2/10, 1/10, 4/5, 1/20, 3/10, 1]


NSolve[%30 == 0]


(*um \[Rule] 1.55262*)


(* ::Subsubsection:: *)
(*b = 1., um = (4/5)*)


dervmutX[Extract[%9[[33,2]],{1,1}],Extract[%9[[33,2]],{1,2}], 1/2, 25, 7/10, 1/10, 1/20, 3/10, 2/10, 1/10, 4/5, 1/20, 3/10, 1]
NSolve[%32 == 0]


(*um \[Rule] 1.54932*)


(* ::Subsubsection:: *)
(*b =1.; um = (33/40)*)


dervmutX[Extract[%9[[34,2]],{1,1}],Extract[%9[[34,2]],{1,2}], 1/2, 25, 7/10, 1/10, 1/20, 3/10, 2/10, 1/10, 4/5, 1/20, 3/10, 1]
NSolve[%34 == 0]


(*um \[Rule] 1.54604*)


(* ::Subsubsection:: *)
(*b= 1., um =(17/20)*)


dervmutX[Extract[%9[[35,2]],{1,1}],Extract[%9[[35,2]],{1,2}], 1/2, 25, 7/10, 1/10, 1/20, 3/10, 2/10, 1/10, 4/5, 1/20, 3/10, 1]
NSolve[%36 == 0]


(*um \[Rule] 1.54278*)


(* ::Subsubsection:: *)
(*b = 1., um =  (7/8)*)


dervmutX[Extract[%9[[36,2]],{1,1}],Extract[%9[[36,2]],{1,2}], 1/2, 25, 7/10, 1/10, 1/20, 3/10, 2/10, 1/10, 4/5, 1/20, 3/10, 1]
NSolve[%38 == 0]


(*um \[Rule] 1.53954*)


(* ::Subsubsection:: *)
(*b = .1, um = (9/10)*)


dervmutX[Extract[%9[[37,2]],{1,1}],Extract[%9[[37,2]],{1,2}], 1/2, 25, 7/10, 1/10, 1/20, 3/10, 2/10, 1/10, 4/5, 1/20, 3/10, 1]
NSolve[%40 == 0]


(*um \[Rule] 1.53633*)


(* ::Subsubsection:: *)
(*b = 1., um = (37/40)*)


dervmutX[Extract[%9[[38,2]],{1,1}],Extract[%9[[38,2]],{1,2}], 1/2, 25, 7/10, 1/10, 1/20, 3/10, 2/10, 1/10, 4/5, 1/20, 3/10, 1]
NSolve[%50 == 0]


(*um \[Rule] 1.53315*)


(* ::Subsubsection:: *)
(*b = 1., um = (19/20)*)


dervmutX[Extract[%9[[39,2]],{1,1}],Extract[%9[[39,2]],{1,2}], 1/2, 25, 7/10, 1/10, 1/20, 3/10, 2/10, 1/10, 4/5, 1/20, 3/10, 1]
NSolve[%52 == 0]


(*um \[Rule] 1.52998*)


(* ::Subsubsection:: *)
(*b = 1., um = (39/40)!!!*)


dervmutX[Extract[%9[[40,2]],{1,1}],Extract[%9[[40,2]],{1,2}], 1/2, 25, 7/10, 1/10, 1/20, 3/10, 2/10, 1/10, 4/5, 1/20, 3/10, 1]
NSolve[%54 == 0]


(*um \[Rule] 1.52684*)


(* ::Subsubsection:: *)
(*b = 1.; um = (1)*)


dervmutX[Extract[%9[[40,2]],{1,1}],Extract[%9[[40,2]],{1,2}], 1/2, 25, 7/10, 1/10, 1/20, 3/10, 2/10, 1/10, 4/5, 1/20, 3/10, 1]
NSolve[% == 0] (** %%,%%%,%%%% goes back 2-4 previous output.)


(*um \[Rule] 1.52684*)













































































(* ::Subsubsection::Closed:: *)
(*b =.9; um = 1/40*)


(*Trying our new function w/Extract! failed.*) 
dervmutXup[17, 2, 1/2, 25, 7/10, 1/10, 1/20, 3/10, 2/10, 1/10, 4/5, 1/20, 3/10, .9]


(* ::Text:: *)
(*p.s. details of doing AD through mutant increase from rare, and a second derivative test, is under section 18 of our overleaf notebook!*)


(*Tests*)
Extract[%17[[2,2]],{1,1}]
Extract[%9[[39,2]],{1,2}]
Extract[%10[[2,2]],{1,1}]



?ListDensityPlot


(*Test.1*)


ListDensityPlot[{{0,25,1},{1/40,9.109329294527447,1},{1,0.5110838337102811,0}}, ColorFunction->ColorData["MintColors"]]


(*Trial ListDensityPlot*)
ListDensityPlot[{1,0,1},{0,1,0}, Mesh-> All]


(*MORE BIG SPACE BEFORE HARDER EVOLUTION HAPPENS.*)


(* ::Text:: *)
(*weird way of doing it by hand. *)


(* ::Section:: *)
(*Finding um_ESS.*)


(* ::Subsubsection:: *)
(*Ecological Equil 2.0.*)


stableEquilr[um_, b_, d_, r0_, c_, v_, s_, ca_, cc_, lx_, g_, l_, z_,  precision_: MachinePrecision] := 
Module[{equilibriar},
  (* Find all equilibria *)
  equilibriar = NSolve[{drdt[r, h, m, d, r0, c, v, s, ca, cc, b, lx, g, l, z, um] == 0, dhdt[r, h, m, d, r0, c, v, s, ca, cc, b, lx, g, l, z, um] == 0, dmdt[r, h, m, d, r0, c, v, s, ca, cc, b, lx, g, l, z, um] == 0, r >= 0, h >= 0, m >= 0} , {r, h, m}, Reals, WorkingPrecision -> precision];
  (*Reals throws out imaginary roots*)
  
  (* Select the stable equilibria *)
  equilibriar =  Select[equilibriar, negativeEigQ[myJacobian[r/.#, h/.#, m/.#, d, r0, c, v, s, ca, cc, b, lx, g, l, z, um]]&];
 
 (* While our the first term of equilibriar is < 25, Drop the second term.
  equilbriar = While[equilibriar #[[1]]& < 25 &; Drop #[[2]]&
  but this will likely not work if parantheses issue is not resolved*)
  
  (* Have equilibria be returned as numbers instead of r \[Rule] number *)
  equilibriar = Map[{r/.#}&, equilibriar];
  
  Return[equilibriar]
]


stableEquilh[um_, b_, d_, r0_, c_, v_, s_, ca_, cc_, lx_, g_, l_, z_,  precision_: MachinePrecision] := 
Module[{equilibriah},
  (* Find all equilibria, Reals throws out imaginary roots *)
  equilibriah = NSolve[{drdt[r, h, m, d, r0, c, v, s, ca, cc, b, lx, g, l, z, um] == 0, dhdt[r, h, m, d, r0, c, v, s, ca, cc, b, lx, g, l, z, um] == 0, dmdt[r, h, m, d, r0, c, v, s, ca, cc, b, lx, g, l, z, um] == 0, r >= 0, h >= 0, m >= 0} , {r, h, m}, Reals, WorkingPrecision -> precision];
  
  
  (* Select the stable equilibria *)
  equilibriah = Select[equilibriah, negativeEigQ[myJacobian[r/.#, h/.#, m/.#, d, r0, c, v, s, ca, cc, b, lx, g, l, z, um]]&]; 
  
  (* Have equilibria be returned as numbers instead of h \[Rule] number *)
  equilibriah = Map[{h/.#}&, equilibriah];
  
  Return[equilibriah]
]


(*um = 0.5; d = 0.5; r0 = 25; c = 7/10; v = 1/10; s = 1/20; ca = 3/10; cc = 2/10; lx = 2/10; g = 4/5; l = 1/20; z = 3/10; b = .5 *)


(*1/2, 1/2, 1/2, 25, 7/10, 1/10, 1/20, 3/10, 2/10, 2/10, 4/5, 1/20, 3/10*)


(* ::Text:: *)
(*Testing out what Sort does*)


Sort[{{25, 0, 0}, {14, 1, 1}}, #1[[1]] < #2[[1]]&] 
(*sorts list so that the first entry of the element is always smaller than the second one *)


?Module
(*Module[{temp1, temp2},temp1 = dervmut[stuff];temp2 = temp1[[1]]; Nsolve[temp2\[Equal]0];*)


(* ::Text:: *)
(*Continue writing function for evaluating evo.*)
(*Traitvalue here is just traitListuno. *)


?Table


(*traitListuno = Table[x , {x, 0, 1, 1/40}]*)


(*d \[Rule] 0.5, r0 \[Rule] 25, c \[Rule] 7/10, v \[Rule] 1/10, s \[Rule] 1/20, ca \[Rule] 3/10, cc \[Rule] 2/10, lx \[Rule] 2/10, g \[Rule] 4/5, l \[Rule] 1/20, z \[Rule] 3/10, b \[Rule] ?*)


(* ::Text:: *)
(*This function should print the um value (that equals to the input[um], ) that makes dervmutx == 0*)
(*	? How to do select  one single "um" from a list of um-s, at each run?? *)


(* Add the trait values that produced the equilibrium/equilibria *)
(*tmpEquil = Table[{traitList[[i, 1]], traitList[[i, 2]],tmpEquil[[i]]} (*{expr}*), {i, Length[traitList]}*) (*generates a list of the values of expr when i runs from 1 to i_max*)];
(* Select the stable equilibria *)
 (* equilibriah = Select[equilibriah, negativeEigQ[myJacobian[r/.#, h/.#, m/.#, d, r0, c, v, s, ca, cc, b, lx, g, l, z, um]]&]; *)
 (*mutX[r_, h_, d_, r0_, c_, v_, s_, ca_, cc_, lx_, g_, l_, z_, b_] := v (1- (un z)) (cc + ((b un r)/(lx + g h)) h - l)*)
 (*dervmutX[r_, h_, d_, r0_, c_, v_, s_, ca_, cc_, lx_, g_, l_, z_, b_]:= D[mutX[r, h, d, r0, c, v, s, ca, cc, lx, g, l, z, b], un]/.un -> um *)


(*/@ Map*)


(*um = 0.5; d = 0.5; r0 = 25; c = 7/10; v = 1/10; s = 1/20; ca = 3/10; cc = 2/10; lx = 2/10; g = 4/5; l = 1/20; z = 3/10; b = .5 *)


dervmutX1


(* ::Subsubsection:: *)
(*While loop notes.*)


(* ::Text:: *)
(*stuff to add:*)


(*This will make um increase if larger um has a higher invasion fitness (and decrease if smaller um has a higher invasion fitness)*)


?While


?Print


?Break 


(* ::Text:: *)
(*Examples of break and while.*)


a = -1;
ohNoes = 0;
While[a < 100, a = 1.5 * a; ohNoes = ohNoes + 1;
If[ohNoes > 1000, Print["a is <= 0, isn't it?"];Break[]]]
ohNoes


a = 1;
While[a < 100, a = 1.5*a(*; 
Return[a]*)
]
a


(* ::Text:: *)
(*Example of what you can define um as.*)


(*um = um + 0.00001 * (dervmutX1/. um -> um)*) (*test the 0.00001 with different degrees of large/smallness*)


?Sort


(* ::Text:: *)
(*how to loop back into finding stabr, stabh and dervmut....*)


(* ::Text:: *)
(*Interpolating function. *)


?Module


(*Sample code :) to show how /. corrects for x being \[Rule] to something*)
stuff = x -> 1
(x /. stuff) == 1
x


(* ::Subsection:: *)
(*um_ESS module~!*)


(* ::Text:: *)
(*um_ess module writing notes :*)
(*need instructions to assign a new um and recall steps 1 and 2, after NSolve is done.*)
(*Define um within module and again within while.*)


(* ::Text:: *)
(*|| or! make an interpolating function outside of esstraitval.*)


(* ::Text:: *)
(*check to see whether function below would fully stop we come across dervmutX = 0*)


1/20*10


Length[2]


stableEquilr[1/20, 7/20, 1/2, 25, 7/10, 1/10, 1/20, 3/10, 2/10, 2/10, 4/5, 1/20, 3/10]


stableEquilh[1/20, 7/20, 1/2, 25, 7/10, 1/10, 1/20, 3/10, 2/10, 2/10, 4/5, 1/20, 3/10]


um_, b_, d_, r0_, c_, v_, s_, ca_, cc_, lx_, g_, l_, z_


Sort[Flatten[stableEquilr[1/20, 7/20, 1/2, 25, 7/10, 1/10, 1/20, 3/10, 2/10, 2/10, 4/5, 1/20, 3/10]]][[1]]


Sort[Flatten[stableEquilh[1/20, 7/20, 1/2, 25, 7/10, 1/10, 1/20, 3/10, 2/10, 2/10, 4/5, 1/20, 3/10]],Greater][[1]]


Sort[{{1}, {3}, {2}}, Greater]


dervmutX[16.860396736734263, 0, 1/2, 25, 7/10, 1/10, 1/20, 3/10, 2/10, 2/10, 4/5, 1/20, 3/10, 7/20]


(* ::Text:: *)
(*We removed this part, and decided that all esstraitval needs is the while loop! *)


(*Sort takes elements in order of their magnitude/value, [[1]] tells it to take n + 1 input of the list with each run*) 
(*stableEquilr1 = Sort[stableEquilr[um, b, d, r0, c, v, s, ca, cc, lx, g, l, z]][[1,1]]; 
stableEquilh1 = Sort[stableEquilh[um, b, d, r0, c, v, s, ca, cc, lx, g, l, z], Greater][[1,1]];

(*Call dervmutX*)
dervmutX1 = dervmutX[stableEquilr1,stableEquilh1, d, r0, c, v, s, ca, cc, lx, g, l, z, b]; 

(*Solving the invasion fitness for when it is zero; and naming it exp*)
exp = NSolve[dervmutX1 == 0, un]; (*outputs un \[Rule] value*)
Print["got to NSolve outside of while."];

(*If if statements holds true, the code would stop here and not run the while loop.*)
(*This line of code says that if the un that you solved from exp equals to the um that you inputed, return um.*) (*Maybe this is why it always outputs umMin?*)
If[(un /. exp) == um, Print["found un == um! "];Return[um]]; (*/. \[Equal] um Corrects for \[Rule]*)  (*Natalie added print but it didn't print. maybe it run past here, and shouldn't be returning um?*)

(*Adjust um in the direction where the fitness is increasing*)
um = um + 0.01 * (dervmutX1/. un -> um);*) (*Take all of the un in dervmutX1 and put in the value of um in it*)

(*Record our ultimate um_ess-es (ultium).Where should I be calling um?*)


(* ::Text:: *)
(*Note about esstraitval : - All the "Print"s that might or might not' ve been commented out, are good ways to check how far your function ran till. *)
(*Think about what changing all of these things changes. *)
(*	i. initial umMin.*)
(*	ii. counter.*)
(*	iii. Chop tolerance. *)
(*	xiv. how far you move away from um!*)


esstraitval[umMinList_, b_, d_, r0_, c_, v_, s_, ca_, cc_, lx_, g_, l_, z_] := 
(*perks of using a module, what happens in the module stays. In the module; i.e. if there ever was something named b, forget it and do the following to it instead.*)
(*Perk above only applies to the stuff that you put into curly braces below! *)
Module[{stableEquilr1, stableEquilh1, dervmutX1, exp, um, ultium},

umIn = 7/20; (*Lowest um value on list, take that as our initial um.*) (*first iteration was done with 3/4!*)
um = umIn;

counter = 0;
While[counter < 3000, (*le test, counter: don't do it more than 5000 times*)
(*Everything below is the body.*)
counter = counter +1; 

(*Calling stable equilibrias*)
stableEquilr1 =stableEquilr[um, b, d, r0, c, v, s, ca, cc, lx, g, l, z];
stableEquilh1 = stableEquilh[um, b, d, r0, c, v, s, ca, cc, lx, g, l, z];

(*Check that there are equilibrias present!*)
If[Length[stableEquilr1] == 0, Print["we tested if there are resource equilibriums, and theres not!"];Print[um]; Break[]];
If[Length[stableEquilh1] == 0, Print["we tested if there are host equilibriums, and theres not!"];Print[um]; Break[]];
(*Print["length is zero?"];*)
(*Chop[Sort[Flatten[stableEquilr1], Greater][[1]], 10^-15] == r0*)
If[Length[stableEquilr1] == 1, Print["Warning: stableequilr is very low!"];Print[um]]; (*Note*)
If[Length[stableEquilh1] == 1, Print["Warning: stableequilh is very low!"];Print[um]];

(*Then, sort and select the right ones.*)
stableEquilr1 =Sort[Flatten[stableEquilr1]][[1]];(*Sort takes elements in order of their magnitude/value, [[1]] tells it to take n + 1 input of the list with each counter run and this only happens when combined with while (because we're doing counter!)*)
stableEquilh1 = Sort[Flatten[stableEquilh1], Greater][[1]];

(*Call dervmutX!*)
dervmutX1 =(*Flatten[*)dervmutX[stableEquilr1,stableEquilh1, d, r0, c, v, s, ca, cc, lx, g, l, z, b](*]*); 

(*Noppers we no longer need this! :X*)
(*/. \[Equal] um Corrects for \[Rule]*)
(* CHECK for dervmutx = 0 works yay! (to check, we did this : dervmutX1 = 0); and the message printed when that happened.*)
(*Solving the invasion fitness for when it is zero. Name it exp*)
(*exp = NSolve[dervmutX1 == 0, un];*) (*outputs um \[Rule] value*)
(*Print["Got to NSolve within while."];*)


(*Message: we're here, means that we got to equilibria! *)
If[Chop[(dervmutX1/. un -> um), 10^-10] == 0, ultium = um ; Print["we're here"]; Break[]];  
(*Note: to change Chop tolerance: Chop[dervmutX1, 10^-15]*)

(*if the un that you solved for doesn't match, continue on with a different um *)
um = um + 0.1 * (dervmutX1/. un -> um);
(*Print["unknown"];Print[dervmutX1]; Print[stableEquilr1];*)

(*When the function should stop*)
If[um >= 1, ultium = 1; Break[]]];
(*Code below says that, without exceeding the maximum counter count, search for the equilibria. If equilibria isn't found when you reach counter, print message and return the um that you got!*)
If[counter == 3000, Print["Searched for equilibrium for 3000 steps. No equil found yet :("]; Return[um], Return[ultium]];

(*exp[um_] \[Equal]  um (*test*), (*body*) Print[exp]; If[exp \[Equal]um, Break[], um = um + 0.00001 * (dervmutX1/. um -> um)]](*otherwise assign new um and find dervmutx and equilbriar&h again*);*)
(*decrease step size and increase counter to check reliabiliity of um value.*)

(*traitvalue = dervmutX[[1]];*)
(*whileloop goes through steps 1,2, until find dervmutx equals zero, looks kinda like modules/ifelse on steriods. probably need a while loop here as well.*)
(*]*)
(*b = traitListuno[[i]];  i = {1, 40};*)
 (*returns derivative when um is equal to the input um value.find um that makes dervmutx \[Equal] 0 *)
Return[um]
]


(* ::Text:: *)
(*Flatten doesn't work for the stableEquils, but doing [[1,1]] would get you what you want! ;)*)


(*Flatten[*)Sort[stableEquilh[1/40, 1/8, 1/2, 25, 7/10, 1/10, 1/20, 3/10, 2/10, 2/10, 4/5, 1/20, 3/10], Greater][[1,1]](*]*)


?While


?Table


(* ::Text:: *)
(*todo 3/15:- make a table for b:- Table[everything you do for a particular value of b, {b, {a list of the b values you want to check}] .*)
(*you can adjust the movement of um, or the total time steps you're going to run through.*)


(*create umMinList*)  (*for each b value take the lowest um that produces a positive m*, skip the stuff before white regions*)


(* ::Text:: *)
(*The five different initial um_s that we're trying .*)


1/20, 7/20, 3/8, 2/5, 17/40


umMinList = {1/20, 1/20, 1/20, 1/20, 1/20, 1/20, 1/20, 1/20, 1/20, 1/20, 1/20, 1/20, 1/20, 1/20, 1/20, 1/20, 1/20, 1/20, 1/20, 1/20, 1/20, 1/20, 1/20, 1/20, 1/20, 1/20, 7/20, 7/20, 3/8, 3/8, 3/8, 3/8, 3/8, 3/8, 3/8, 2/5, 2/5, 2/5, 2/5, 2/5, 17/40}


abe = Table[b, {b, 0, 1, 1/40}]


Extract[abe, 1]


(* ::Text:: *)
(*Doesn' t work, trying again tomorrow :) - problem with esstraitval.*)


(*Sample*)
(*Table[ess[traitval[abe[[i]], umMinList[[i]], ..,], {i, 1, 41}]*)


esstraitval[umMinList[[i]], abe[[i]], 1/2, 25, 7/10, 1/10, 1/20, 3/10, 2/10, 2/10, 4/5, 1/20, 3/10]


?Table


esstraitval[1/20, 7/20, 1/2, 25, 7/10, 1/10, 1/20, 3/10, 2/10, 2/10, 4/5, 1/20, 3/10]


(*Testing flatten's capabilities.*)
Flatten[.50]


(* ::Text:: *)
(*adding .25 to be further away from edges*)


(* ::Text:: *)
(*3/29 :- check for a length of zero, break and output that um. GOOGLE Try-Catch Mathematica.*)


(*Next steps: *)
1. Run it for longer (increasing counter...) doesn't work.. 
2. Checking thoroughly by making up our fitness function, so that the function is able to stop when dervmutX is zero . Because if there is mulitple strategies that could evolve we wanna know!
Specifically making sure that this -> If[dervmutX == 0, then stop] is working well! using print, and making a fake dervmutx that is simple making a mock trial .
3. Pondering biological meaning of these results . (changing equations? (**maybe, if they don't make sense .**) symbiont is also looking out for itself .
	x . other desires! roam and be free!


(*d \[Rule] 1/2, r0 \[Rule] 25, c \[Rule] 7/10, v \[Rule] 1/10, s \[Rule] 1/20, ca \[Rule] 3/10, cc \[Rule] 2/10, lx \[Rule] 2/10, g \[Rule] 4/5, l \[Rule] 1/20, z \[Rule] 3/10, b \[Rule] ?*)


umMinList; abe










(* ::Subsubsection:: *)
(*4/5 :- Todo list, and esstraitval outputs! *)


(* ::Text:: *)
(*Tables are so useful.*)
(*Check values higher from final 2 umMins how many stable equilibrias there are .*)


Check
Table[Length[stableEquilr[um, 39/40, 1/2, 25, 7/10, 1/10, 1/20, 3/10, 2/10, 2/10, 4/5, 1/20, 3/10]], {um, 2/5, 1, 1/1000}]


(* ::Text:: *)
(*Step 1: For now Natalie is fixing everything else while tryna increase the umMin that we've initiated the function in.*)
(*1st: Chop = 10^-15, counter at 100, umMinList at 1/20. *)
(*2nd: This run is for when Chop = 10^-10, and counter is at , umMinList is at 1/20! output has 41 we're heres and it equilibrates at 1/20. *)
(*-- Adjusting Chop doesn't affect output, counter 3000 seems to work well for it.*)
(* *initial at 3/8 shows interesting results! should try increasing counter for that one! increased till 5000 & chop at 10^-5. *)
(* 6th iteration: uses 17/40*)
(*7th iteration: uses 5000 counter, chops at 10^-5, worked out for 3/8, and 2/5. *)
(*8th: 5000 counter but chops at 10^10, gunna run after this but with increased chop to 10^-5. :)*)
(* After alll iterations, should figure out how to create a Listlineplot outta this.*)


ott = Table[esstraitval[umMinList[[i]](*+.5*), abe[[i]], 1/2, 25, 7/10, 1/10, 1/20, 3/10, 2/10, 2/10, 4/5, 1/20, 3/10],{i, 1, 41}]


(*Next Steps*)
1. Mess around with iterations and chop for it to not output "No equil found"! [\done]
	i . Increasing counter max count; [\done, 5000]
	ii . Increasing the size of chop (10^-15(guesstimated (alexandra word)) < chop < 10^-5(fixed)), [\done, chop ]
5. Annotate and upload everything onto GitHub . 
6. Organize Overleaf! 


(* ::Subsubsection:: *)
(*Side Quest! Function for the invasion fitness.*)


dervmutEx[umMinList_, b_, d_, r0_, c_, v_, s_, ca_, cc_, lx_, g_, l_, z_] := 
(*perks of using a module, what happens in the module stays. In the module; i.e. if there ever was something named b, forget it and do the following to it instead.*)
(*Perk above only applies to the stuff that you put into curly braces below! *)
Module[{stableEquilr1, stableEquilh1, dervmutX1, exp, um, un1, ultium},

umIn = 7/20; (*Lowest um value on list, take that as our initial um.*) (*first iteration was done with 3/4!*)
um = umIn;

counter = 0;
While[counter < 3000, (*le test, counter: don't do it more than 5000 times*)
(*Everything below is the body.*)
counter = counter +1; 

(*Calling stable equilibrias*)
stableEquilr1 =stableEquilr[um, b, d, r0, c, v, s, ca, cc, lx, g, l, z];
stableEquilh1 = stableEquilh[um, b, d, r0, c, v, s, ca, cc, lx, g, l, z];

(*Check that there are equilibrias present!*)
If[Length[stableEquilr1] == 0, Print["we tested if there are resource equilibriums, and theres not!"];Print[um]; Break[]];
If[Length[stableEquilh1] == 0, Print["we tested if there are host equilibriums, and theres not!"];Print[um]; Break[]];

(*Then, sort and select the right ones.*)
stableEquilr1 =Sort[Flatten[stableEquilr1]][[1]];(*Sort takes elements in order of their magnitude/value, [[1]] tells it to take n + 1 input of the list with each counter run and this only happens when combined with while (because we're doing counter!)*)
stableEquilh1 = Sort[Flatten[stableEquilh1], Greater][[1]];

(*Call dervmutX!*)
dervmutX1 =(*Flatten[*)dervmutX[stableEquilr1,stableEquilh1, d, r0, c, v, s, ca, cc, lx, g, l, z, b](*]*); 
(*Learn using the cool slant dot thingie. *)
un1 = Table[dervmutX1/.un -> p, {p, 0, 1, 1/40}]; 

Return[un1]
]]


25, .7, .1, .05, .3, .2, .2, .8, .05, .3, .3


jaco = dervmutX[2, 1, 1/2, 25, 7/10, 1/10, 1/20, 3/10, 2/10, 2/10, 4/5, 1/20, 3/10, 9/10]
un1 = Table[jaco/.un -> .3]


?CloudEvaluate


$CloudCreditsAvailable


?Table


umMinList
abe


?Repeat


Table[dervmutEx[umMinList[[i]], abe[[i]], 1/2, 25, 7/10, 1/10, 1/20, 3/10, 2/10, 2/10, 4/5, 1/20, 3/10],{i, 1, 41}]


(* ::Text:: *)
(*Useful run below . Wooo!*)


dervmutExTry = Flatten[Table[dervmutEx[traitListuno[[i]], abe[[i]], 1/2, 25, 7/10, 1/10, 1/20, 3/10, 2/10, 2/10, 4/5, 1/20, 3/10],{i, 1, 41}]]


dervmutExList = Flatten[Table[dervmutEx[umMinList[[i]], abe[[i]], 1/2, 25, 7/10, 1/10, 1/20, 3/10, 2/10, 2/10, 4/5, 1/20, 3/10],{i, 1, 41}]]


?Extract


Extract[%81, 1681]


?ListDensityPlot


(*Modify this list to make it dervmutEx; umMin; and abe*)
finalistm[traitLister, d_, r0_, c_, v_, s_, ca_, cc_, lx_, g_, l_, z_, precision_:MachinePrecision] := 
Module[{tempo}, (*it doesn't matter if you'd used tempo before because you're in a module, how cool!*)

(* Make a table of final stable equilibrias *)
tempo=raeb2;

(* Add the trait values that produced the equilibrium/equilibria *)
tempo = Table[{traitLister[[i, 1]], traitLister[[i, 2]], tempo[[i]]}, {i, Length[traitLister]}];

Return[tempo]]


(* ::Text:: *)
(*4/18: Okay, current issue is that umMinList is too short! Ask Alexandra on Monday on how to repeat an element! :D*)


(*to repeat an element:*)
(*Flatten[Table[coffee, {i, 6}], 1]*)


coffee = {1, 2, 3, 4}
Flatten[Table[coffee, {i, 6}], 1]
Flatten[Table[coffee[[j]], {j, Length[coffee]},{i, 6}](*, 2*)]


umMinListLong = Flatten[Table[umMinList, {i, 41}], 1]


abeLong = Flatten[Table[abe[[j]], {j, Length[abe]},{i, 41}](*, 2*)]


traitlistunoLong = Flatten[Table[traitListuno, {i, 41}], 1]


(*tempo = Table[{tempo[[i]]}, {i, Length[traitLister]}]*)


finaldervmutList[traitlistunoLong, abeLong, d_, r0_, c_, v_, s_, ca_, cc_, lx_, g_, l_, z_] := 
Module[{dervmutEx}, (*it doesn't matter if you'd used the stuff inside curley braces before because you're in a module, how cool!*)

(* Make a table of final stable equilibrias *)
dervmutEx = dervmutExTry;

(* Add the trait values that produced the equilibrium/equilibria *)
dervmutEx = Table[{traitlistunoLong[[i]], abeLong[[i]], dervmutExTry[[i]]}, {i, 1, 1681}];

Return[dervmutEx]]


Length[abeLong]


finaldervmutList[traitlistunoLong, abeLong, 1/2, 25, 7/10, 1/10, 1/20, 3/10, 2/10, 2/10, 4/5, 1/20, 3/10]


0


?ListDensityPlot


Histogram[Flatten[finaldervmutList[traitlistunoLong, abeLong, 1/2, 25, 7/10, 1/10, 1/20, 3/10, 2/10, 2/10, 4/5, 1/20, 3/10][[All, 2]]]] (*3 means just give the third entry.*)


ListDensityPlot[finaldervmutList[traitlistunoLong, abeLong, 1/2, 25, 7/10, 1/10, 1/20, 3/10, 2/10, 2/10, 4/5, 1/20, 3/10]]


(* ::Text:: *)
(*#CandyColors*)


(*multiplication stretches 0 \[Rule] 1 so that you can have axes values at higher values, +shifts center to .5*)


mycolorfun[col_] := If[col != 0, ColorData["GreenPinkTones"][col*10+0.5] , "Black"]


(*At certain fidelities, no candidate ESS at aprox .05, and high fidelity .*)

ListDensityPlot[finaldervmutList[traitlistunoLong, abeLong, 1/2, 25, 7/10, 1/10, 1/20, 3/10, 2/10, 2/10, 4/5, 1/20, 3/10], ColorFunction-> mycolorfun, DataRange -> Automatic, ColorFunctionScaling-> False, PlotLegends-> Automatic, PlotRange -> {Full, Full, {-0.05, 0.05}},ClippingStyle -> {White, Pink}, FrameLabel -> {"Uptake Rate, um", "fidelity, b"}]


finaldervmutList[traitlistunoLong, abeLong, 1/2, 25, 7/10, 1/10, 1/20, 3/10, 2/10, 2/10, 4/5, 1/20, 3/10]//MatrixForm


(* ::Text:: *)
(*dervmutX is zero by gray region where b = .2*)


?If


(* ::Subsection:: *)
(*stuff! imp Stuff~!*)


(* ::Text:: *)
(*4/16:- *)
(*! Print out resident(*resource or host would work too) biomass at um_ess to see what's happening at the region between color and blankness at low um.*)
(*- For esstraitval, find counter that produces the first warning message for our stableequils, printing out the um for slightly before and slightly after that happens, *)
(*- !! Did we take the the derivative of the mutant growth rate correctly? DOUBLE- CHECK. that the derv-mutX is in face the derivative of the mutant growth rate, taken in terms of the mutant trait value.  *)
(*- Stray minus signs, double checking equations..*)
(*? If h is zero then dervmutX should be zero, but for some reason it isn't?!*)


(* ::Text:: *)
(*4/30:-*)
(*-Ponder what that gray region from side quest means. If something is capping the mutant growth rate (the derivative), that might be problematic.*)
(*Ways: *)
(*1. Stare at equations and see if anything occurs.... *)
(**)
(*2. Check code that calculated dervmutX, if its unclear recode it and see if it does the same thing. *)
(*- try numbers where b is high/ numbers that you know the solution of dervmutX to, and compare to see if its the same solution as what mathematica calculated out (from running dervmutX with those numbers)*)
(*might be that true dervmutX is that value we just didn't think through it when we wrote the model. *)
(**)
(*3. Plots of H,M,R, see if they change a lot at high b (*especially M*), -> this might be a problem. OTHERWISE MAKE NEW LDPs!*)
(*doesn't change a lot at high b? According to LDP made before..double check how code was written for LDPs on overleaf, they seemed flipped. *)
(**)
(**)
(*-Archived.-*)
(*A line that shows how the zero values change, with b would be helpful. Trace it in code!*)
(*With this! *)
(*If*)


mycolorfun[col_] := If[col != 0, ColorData["GreenPinkTones"][#*10+0.5]& , "Black"]


(* ::Text:: *)
(*Possible issues : (1) typo, (2) should be Ifelse, (3) "Black" might have to be a number/RGB value ?How to implement within..*)


(* ::Text:: *)
(*helps mark specific points on the plot . *)


(* ::Text:: *)
(*Needs to make interpolating functions for them.*)


ListPlot[{{x1, y1}, {x2, y2}}]


Show[LDP, ListPlot]


(* ::Text:: *)
(*5/21 : - FIGURE OUT WHY b BEING HIGH MAKES EVERYTHING ZERO?! *)
(*- check out 4/30!*)
(**)
(*Future: try different parameter values! *)


(* ::Section::Closed:: *)
(*Plotting *)


(* ::Subsubsection:: *)
(*MutGrowth. we'd already done this. *)


(* ::Text:: *)
(*Because finding eigenvalues exactly does not always go smoothly, I put a time constraint on it. If the exact eigenvalues can't be found quickly, it gives up and tries to find the eigenvalues of a numerical approximation.*)


mutGrowth[iP_, iQ_, hmut_, vmut_, hres_, vres_, d_] := 
TimeConstrained[
	Eigenvalues[mutX[iP, iQ, hmut, vmut, hres, vres, d], 1], (* Eigenvalues[matrix, 1] gives the leading eigenvalue *)
	timeOut, (* try to find eigenvalues for timeOut number of seconds *)
	Eigenvalues[N[mutX[iP, iQ, hmut, vmut, hres, vres, d]], 1] (* Do this if the previous calculation timed out *)
][[1]] (* This [[1]] is because eigenvalues are always returned as a list, even if there's only one *)


(* ::Text:: *)
(*To plot stable equilibria, I use ListDensityPlot. The code below will make 2 separate plots:*)
(*One of the iP values for the first equilibrium for each entry in myStableEquil*)
(*One of the iQ values for the first equilibrium for each entry in myStableEquil*)


myStableEquil = {}; (* results of running stableEquilList *)

(* Get the values of iP and iQ at the first stable equilibrium listed for each trait combination *)
toPlotiP = Map[{#[[1]], #[[2]], #[[3, 1, 1]]}&, myStableEquil];
toPlotiQ = Map[{#[[1]], #[[2]], #[[3, 1, 2]]}&, myStableEquil];

ListDensityPlot[toPlotiP, PlotLabel -> "Ecological equilibrium values of iP", AxesLabel -> {"h", "v"}]
ListDensityPlot[toPlotiQ, PlotLabel -> "Ecological equilibrium values of iQ", AxesLabel -> {"h", "v"}]


(* Get any points with multiple stable equilibria *)
Select[myStableEquil, Length[#[[3]]] > 1&]

(* You can plot the second (or 3rd, 4th, etc.) equilibria of a point with this code. *)
myStableEquil2 = Select[myStableEquil, Length[#[[3]]] > 1&];
toPlotiP2 = Map[{#[[1]], #[[2]], #[[3, 2, 1]]}&, myStableEquil2];
toPlotiQ2 = Map[{#[[1]], #[[2]], #[[3, 2, 2]]}&, myStableEquil2];

ListDensityPlot[toPlotiP2, PlotLabel -> "2nd stable ecological equilibrium values of iP", AxesLabel -> {"h", "v"}]
ListDensityPlot[toPlotiQ2, PlotLabel -> "2nd stable ecological equilibrium values of iQ", AxesLabel -> {"h", "v"}]


(* ::Text:: *)
(*To plot evolutionary trajectories, I use ListStreamPlot*)


mydGrowth = {}; (* results of running dGrowthList *)

(* Get the values of dGrowthdh and dGrowthdV at the first stable equilibrium listed for each trait combination *)
toPlotGrowth = Map[{{#[[1]], #[[2]]}, {#[[3, 1, 1]], #[[3, 1, 2]]}}&, mydGrowth];

ListStreamPlot[toPlotGrowth, PlotLabel -> "Derivative of invasion fitness", AxesLabel -> {"h", "v"}]

(* Plot on top of the ecological equilibria for a fancy look *)
Show[ListDensityPlot[toPlotiP, PlotLabel -> "Ecol. equil. iP, arrows show derivative of invasion fitness", AxesLabel -> {"h", "v"}], 
ListStreamPlot[toPlotGrowth]]






