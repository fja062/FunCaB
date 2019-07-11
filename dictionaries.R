# dictionaries
dict_Site <- read_delim(delim = " ", file = 
"old v2 v3 new
Arh ARH arh Arhelleren
Ovs OVS ovs Ovstedal
Ves VES ves Veskre
Skj SKJ skj Skjellingahaugen
Lav LAV lav Lavisdalen
Gud GUD gud Gudmedalen
Ulv ULV ulv Ulvhaugen
Vik VIK vik Vikesland
Hog HOG hog Hogsete
Alr ALR alr Alrust
Fau FAU fau Fauske
Ram RAM ram Rambera")


sitePref <- read_delim(delim = " ", file = 
"siteID pref alternative
Arhelleren 1 Ovstedal
Ovstedal 1 Arhelleren
Veskre 1 Arhelleren
Veskre 2 Ovstedal
Skjellingahaugen 1 Rambera
Skjellingahaugen 2 Arhelleren
Rambera 1 Skjellingahaugen
Gudmedalen 1 Vikesland
Gudmedalen 2 Lavisdalen
Lavisdalen 1 Vikesland
Lavisdalen 2 Gudmedalen
Vikesland 1 Hogsete
Vikesland 2 Gudmedalen
Hogsete 1 Vikesland
Hogsete 2 Gudmedalen
Alrust 1 Fauske
Fauske 1 Alrust
Ulvhaugen 1 Alrust
Ulvhaugen 2 Fauske")


dict_TTC_turf <- read_delim(delim = ";", file =
"TTtreat;turfID
51 TTC;Fau1C
57 TTC;Fau2C
68 TTC;Fau4C
73 TTC;Fau5C
29 TTC;Alr1C
31 TTC;Alr2C
37 TTC;Alr3C
134 TTC;Vik2C
140 TTC;Vik3C
141 TTC;Vik4C
146 TTC;Vik5C
101 TTC;Hog1C
110 TTC;Hog2C
115 TTC;Hog3C
286 TTC;Ovs1C
291 TTC;Ovs2C
297 TTC;Ovs3C
211 TTC;Arh1C
222 TTC;Arh3C
226 TTC;Arh4C
263 TTC;Ves1C
281 TTC;Ves4C
194 TTC;Ram4C
198 TTC;Ram5C
6 TTC;Ulv2C
11 TTC;Ulv3C
236 TTC;Skj1C
243 TTC;Skj2C
246 TTC;Skj3C
251 TTC;Skj4C
511 TTC;Gud12C
") #516 TTC;Gud13C | 506 TTC;Gud5C <- temporarily taken out; could be that it's not actually one of ours
