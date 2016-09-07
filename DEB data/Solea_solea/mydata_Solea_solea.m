function [data, auxData, metaData, txtData, weights] = mydata_Solea_solea

%% set metaData
metaData.phylum     = 'Chordata'; 
metaData.class      = 'Actinopterygii'; 
metaData.order      = 'Pleuronectiformes'; 
metaData.family     = 'Soleidae';
metaData.species    = 'Solea_solea'; 
metaData.species_en = 'common sole'; 
metaData.T_typical  = C2K(10); % K, body temp
metaData.data_0     = {'ab'; 'aj'; 'ap'; 'am'; 'Lb'; 'Lj'; 'Lp'; 'Li'; 'Wwb'; 'Wwp'; 'Wwi'; 'Ri'}; 
metaData.data_1     = {'t-L'; 'L-N'; 'T-ab'; 'T-aj'}; 

metaData.COMPLETE = 2.5; % using criteria of LikaKear2011

metaData.author   = {'Lorna Teal'};    
metaData.date_subm = [2011 08 24];              
metaData.email    = {'lorna.teal@wur.nl'};            
metaData.address  = {'IMARES, IJmuiden'}; 

metaData.date_acc = [2015 08 28];          

%% set data
% zero-variate data

data.ab = 15;      units.ab = 'd';    label.ab = 'age at birth';              bibkey.ab = 'Fond1979';   
  temp.ab = C2K(16);  units.temp.ab = 'K'; label.temp.ab = 'temperature';
data.tj = 52;      units.tj = 'd';    label.tj = 'time since birth at metam'; bibkey.tj = 'Fond1979';
  temp.tj = C2K(10);  units.temp.tj = 'K'; label.temp.tj = 'temperature';
data.ap = 3 * 365; units.ap = 'd';    label.ap = 'age at puberty';            bibkey.ap = 'MollKraa2007';
  temp.ap = C2K(10);  units.temp.ap = 'K'; label.temp.ap = 'temperature';
  comment.ap = '50% of age 3 are mature';
data.am = 26 * 365; units.am = 'd';   label.am = 'life span';                 bibkey.am = 'Deni1990';   
  temp.am = C2K(10);  units.temp.am = 'K'; label.temp.am = 'temperature'; 
  comment.am = 'temp is Northsea estimate'; 

data.Lb  = 0.425;  units.Lb  = 'cm';  label.Lb  = 'total length at birth';    bibkey.Lb  = 'Fond1979';  
data.Lj  = 0.95;   units.Lj  = 'cm';  label.Lj  = 'total length at metam';    bibkey.Lj  = 'Fond1979';  
data.Lp  = 25;     units.Lp  = 'cm';  label.Lp  = 'total length at puberty';  bibkey.Lp  = 'MollKraa2007'; 
data.Li  = 70;     units.Li  = 'cm';  label.Li  = 'ultimate total length';    bibkey.Li  = 'Whee1978';
  comment.Li = 'size rarely goes above but specimen of 70 has been reported; North Sea data rarely above 45cm!';

data.Wwb = 3e3*(0.425/70)^3; units.Wwb = 'g'; label.Wwb = 'wet weight at birth'; bibkey.Wwb = 'guessed';
  comment.Wwb = 'computed as 3e3*(0.425/70)^3';
data.Wwp = 200;   units.Wwp = 'g';    label.Wwp = 'wet weight at puberty';    bibkey.Wwp = 'MollKraa2007';
  comment.Wwp = 'f ~0.6, T = 10 (North sea guestimates)';
data.Wwi = 3000;  units.Wwi = 'g';    label.Wwi = 'ultimate wet weight';      bibkey.Wwi = 'fishbase';

data.R45  = 1.1e6/365; units.R45  = '#/d'; label.R45  = 'reprod rate at 45 cm'; bibkey.R45  = 'WittGree1995';   
temp.R45 = C2K(10); units.temp.R45 = 'K'; label.temp.R45 = 'temperature';
 
% uni-variate data
% t-L data from otolith back-calculation data for females!(IMARES data)
data.tL = [ ... % time since birth (a), length (cm)
1	8.909404598
2	19.91239656
3	26.95918876
4	31.44085313
5	34.39816096
6	36.2481578
7	37.61539273
8	38.32918647
9	39.16698018
10	39.22343839
11	39.87342987
12	40.96987132
13	42.20480305
14	42.43587572
15	44.03172643];
data.tL(:,1) = 365 * data.tL(:,1); % covert a to d
units.tL   = {'d', 'cm'};  label.tL = {'time since birth', 'total length'};  
temp.tL    = C2K(10);  units.temp.tL = 'K'; label.temp.tL = 'temperature';
bibkey.tL = 'Teal2011';
  
% T-ab data  from Fond1979
% laboratory experiment - reliable data
data.Tab = [... temperature (C), incubation time (d; fertilisation to first feeding)
10 27.5
13 19.75
16 15
19 11.5];
units.Tab = {'deg C', 'd'}; label.Tab = {'temperature', 'age at birth'};  
bibkey.Tab = 'Fond1979';
comment.Tab = 'incubation time: from fertilisation to first feeding)';

% data  from Fond1979
% laboratory experiment - reliable data, assume f = 1
data.Tabj = [... % temperature (C), Development time (first feeding to metamorphisis)
10 52 
13 28
16 20
19 19
22 17.5];
units.Tabj = {'deg C', 'd'}; label.Tabj = {'temperature', 'time since birth at metam'};  
bibkey.Tabj = 'Fond1979';
comment.Tabj = 'development time: from first feeding to metamorphisis)';

% data from WittGree1995 for the Southern Bigth (NS)
%  assume f ~0.6, T = 273 + 10 K (North sea guestimates)
data.LN = [... % Length (cm) vs fecundity (annual egg prod) 
25	100000
30	200000
35	400000
40	700000
45	1100000];
units.LN   = {'cm', '#'};  label.LN = {'total length', 'fecundity'};  
temp.LN    = C2K(10);  units.temp.LN = 'K'; label.temp.LN = 'temperature';
bibkey.LN = 'WittGree1995';

%% set weights for all real data
weights = setweights(data, []);

%% set pseudodata and respective weights
[data, units, label, weights] = addpseudodata(data, units, label, weights);

%% pack auxData and txtData for output
auxData.temp = temp;
txtData.units = units;
txtData.label = label;
txtData.bibkey = bibkey;
txtData.comment = comment;

%% Discussion points
D1 = 'Author_mod_1: I found information on the number of eggs per female as a function of length in Anon2013 that was much higher than in Anon2015 but chose to not include it as the temperature was not provided';
D2 = 'Author_mod_1: I was surprised to observe that the weights coefficient for ab changed so much the parameter values';     
metaData.discussion = struct('D1', D1, 'D2', D2);

%% Facts
F1 = 'The larval stage lasts 202 days and no feeding occurs';
metaData.bibkey.F1 = 'Wiki'; % optional bibkey
metaData.facts = struct('F1',F1);

%% References
bibkey = 'Wiki'; type = 'Misc'; bib = ...
'howpublished = {\url{http://en.wikipedia.org/wiki/Solea_solea}}';
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'Kooy2010'; type = 'Book'; bib = [ ...  % used in setting of chemical parameters and pseudodata
'author = {Kooijman, S.A.L.M.}, ' ...
'year = {2010}, ' ...
'title  = {Dynamic Energy Budget theory for metabolic organisation}, ' ...
'publisher = {Cambridge Univ. Press, Cambridge}, ' ...
'pages = {Table 4.2 (page 150), 8.1 (page 300)}, ' ...
'howpublished = {\url{http://www.bio.vu.nl/thb/research/bib/Kooy2010.html}}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'Fond1979'; type = 'Article'; bib = [ ... 
'author = {Fonds, M.}, ' ... 
'year = {1979}, ' ...
'title = {Laboratory observations on the influence of temperature and salinity on the devlopment of the eggs and groth of the larvae of Solea solea (Pisces)}, ' ...
'journal = {Marine Ecology Progress Series}, ' ...
'volume = {1}, ' ...
'pages = {91--99}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'MollKraa2007'; type = 'Article'; bib = [ ... 
'author = {Mollet, F. and Kraak, S. and Rijnsdorp, A. }, ' ... 
'year = {2007}, ' ...
'title = {Fisheries-induced evolutionary changes in maturation reaction norms in North Sea sole Solea solea.}, ' ...
'journal = {Marine Ecology Progress Series}, ' ...
'volume = {351}, ' ...
'pages = {189--199}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'Whee1978'; type = 'Article'; bib = [ ... 
'author = {Wheeler, A.}, ' ... 
'year = {1978}, ' ...
'title = {Key to the fishes of northern Europe.}, ' ...
'publisher = {Frederik Warne}, ' ...
'address = {London}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'Deni1990'; type = 'Article'; bib = [ ... 
'author = {Deniel, C.}, ' ... 
'year = {1990}, ' ...
'title = {Comparative study of growth of flatfishes on the west coast of Brittany.}, ' ...
'journal = {J. Fish Biol.}, ' ...
'volume = {37}, ' ...
'number = {1}, '...
'pages = {149--166}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'WittGree1995'; type = 'Article'; bib = [ ... 
'author = {Witthames, P. R. and Green Walker, M. and Dinis, M. T. and Whiting, C. L.}, ' ... 
'year = {1995}, ' ...
'title = {The geographical variation in the potential annual fecundity of Dover sole Solea solea (L.) from european shelf waters during 1991}, ' ...
'journal = {Netherlands Journal of Sea Research}, ' ...
'volume = {34}, ' ...
'number = {1-3}, '...
'pages = {45--58}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'Teal2011'; type = 'Misc'; bib = ...
'note = {IMARES data base frisbe}';
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];

